"""
Differential expression analysis for proteomics data.

Provides statistical testing to identify proteins with significant abundance
differences between experimental conditions.
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import digamma, polygamma
from statsmodels.stats.multitest import multipletests

from mokume.core.logger import get_logger

logger = get_logger("mokume.analysis.de")


def _trigamma(x):
    """Compute the trigamma function (second derivative of log-gamma)."""
    return polygamma(1, x)


def _tetragamma(x):
    """Compute the tetragamma function (third derivative of log-gamma)."""
    return polygamma(2, x)


class DifferentialExpression:
    """
    Differential expression analysis for protein intensity data.

    Parameters
    ----------
    method : str
        ``"limrots"``, ``"deqms"``, or ``"proda"``.
    log2fc_threshold : float
        Minimum |log2FC| to call a protein UP or DOWN (default 0.5).
    fdr_threshold : float
        Maximum adjusted p-value (FDR) for significance (default 0.05).
    skip_log2 : bool
        If True, skip the log2 transformation (data is already in log2 space, e.g., from ratio-based quantification).
    fdr_method : str
        Multiple-testing correction method: ``"bh"`` (Benjamini-Hochberg, default) or ``"ihw"`` (Independent Hypothesis Weighting).
    n_boot : int
        Number of bootstrap iterations for LimROTS (default 100).
    n_threads : int or None
        Number of parallel workers for LimROTS bootstrap.
        ``None`` (default) uses all available CPU cores.
    """

    SUPPORTED_METHODS = ("limrots", "deqms", "proda")

    def __init__(
        self,
        method: str = "limrots",
        log2fc_threshold: float = 0.5,
        fdr_threshold: float = 0.05,
        skip_log2: bool = False,
        fdr_method: str = "bh",
        n_boot: int = 100,
        n_threads: int | None = None,
        peptide_counts: pd.Series | None = None,
    ):
        self.method = method.lower()
        if self.method not in self.SUPPORTED_METHODS:
            raise ValueError(
                f"Unknown DE method '{method}'. "
                f"Supported: {self.SUPPORTED_METHODS}"
            )
        self.log2fc_threshold = log2fc_threshold
        self.fdr_threshold = fdr_threshold
        self.skip_log2 = skip_log2
        self.fdr_method = fdr_method.lower()
        self.n_boot = n_boot
        self.n_threads = n_threads
        self.peptide_counts = peptide_counts

    def run(
        self,
        protein_df: pd.DataFrame,
        sample_to_condition: dict[str, str],
        contrast: tuple[str, str],
    ) -> pd.DataFrame:
        """
        Run differential expression analysis for a single contrast.

        Parameters
        ----------
        protein_df : pd.DataFrame
            Wide-format protein matrix (ProteinName | sample1 | sample2 | ...).
        sample_to_condition : dict
            Mapping from sample name to condition.
        contrast : tuple[str, str]
            Pair of conditions to compare, e.g. ``("NASH", "HL")``.
            The fold change is computed as ``contrast[0] / contrast[1]``.

        Returns
        -------
        pd.DataFrame
            Results with columns: ProteinName, log2FC, pvalue, adj_pvalue,
            significance, mean_condition1, mean_condition2.
        """
        cond_a, cond_b = contrast
        log2_matrix, samples_a, samples_b = self._prepare_matrix(
            protein_df, sample_to_condition, cond_a, cond_b,
        )

        dispatch = {
            "limrots": self._run_limrots,
            "deqms": self._run_deqms,
            "proda": self._run_proda,
        }
        runner = dispatch[self.method]
        return runner(log2_matrix, samples_a, samples_b, cond_a, cond_b)

    def _prepare_matrix(
        self,
        protein_df: pd.DataFrame,
        sample_to_condition: dict[str, str],
        cond_a: str,
        cond_b: str,
    ) -> tuple[pd.DataFrame, list[str], list[str]]:
        """Validate inputs, split samples, and log2-transform."""
        protein_col = protein_df.columns[0]
        sample_cols = [c for c in protein_df.columns if c != protein_col]
        samples_a, samples_b = _split_samples(
            sample_cols, sample_to_condition, cond_a, cond_b,
        )
        logger.info(
            "DE analysis: %s (%d samples) vs %s (%d samples)",
            cond_a, len(samples_a), cond_b, len(samples_b),
        )
        intensity = protein_df.set_index(protein_col)[sample_cols]
        log2_mat = _log2_transform(intensity, self.skip_log2)
        return log2_mat, samples_a, samples_b

    def _run_limrots(
        self,
        log2_matrix: pd.DataFrame,
        samples_a: list[str],
        samples_b: list[str],
        cond_a: str,
        cond_b: str,
    ) -> pd.DataFrame:
        """Run LimROTS: limma + ROTS bootstrap-optimized test.

        Delegates to :func:`mokume.analysis.limrots.run_limrots`.
        """
        from mokume.analysis.limrots import run_limrots

        result = run_limrots(
            log2_matrix, samples_a, samples_b,
            cond_a, cond_b, n_boot=self.n_boot, n_threads=self.n_threads,
        )
        if result.empty:
            return result
        return self._finalize_results(result)

    def _run_deqms(
        self,
        log2_matrix: pd.DataFrame,
        samples_a: list[str],
        samples_b: list[str],
        cond_a: str,
        cond_b: str,
    ) -> pd.DataFrame:
        """Run DEqMS analysis with peptide-count-weighted prior variance.

        Delegates to :func:`mokume.analysis.deqms.run_deqms`.
        """
        from mokume.analysis.deqms import run_deqms

        result = run_deqms(
            log2_matrix, samples_a, samples_b,
            cond_a, cond_b,
            peptide_counts=self.peptide_counts,
        )
        if result.empty:
            return result
        return self._finalize_results(result)

    def _run_proda(
        self,
        log2_matrix: pd.DataFrame,
        samples_a: list[str],
        samples_b: list[str],
        cond_a: str,
        cond_b: str,
    ) -> pd.DataFrame:
        """Run proDA probabilistic dropout analysis.

        Delegates to :func:`mokume.analysis.proda.run_proda`.
        """
        from mokume.analysis.proda import run_proda

        result = run_proda(
            log2_matrix, samples_a, samples_b, cond_a, cond_b
        )
        if result.empty:
            return result
        return self._finalize_results(result)

    def _finalize_results(self, de_df: pd.DataFrame) -> pd.DataFrame:
        """Apply FDR correction and classify significance."""
        if self.fdr_method == "ihw" and "adj_pvalue" not in de_df.columns:
            de_df["adj_pvalue"] = _ihw_correction(
                de_df["pvalue"].values,
                de_df,
                alpha=self.fdr_threshold,
            )
        elif "adj_pvalue" not in de_df.columns:
            # Benjamini-Hochberg (default)
            de_df["adj_pvalue"] = multipletests(
                de_df["pvalue"].values, method="fdr_bh"
            )[1]

        # Classify significance
        de_df["significance"] = "Unchanged"
        de_df.loc[
            (de_df["adj_pvalue"] < self.fdr_threshold)
            & (de_df["log2FC"] > self.log2fc_threshold),
            "significance",
        ] = "UP"
        de_df.loc[
            (de_df["adj_pvalue"] < self.fdr_threshold)
            & (de_df["log2FC"] < -self.log2fc_threshold),
            "significance",
        ] = "DOWN"

        # Sort by adjusted p-value
        de_df = de_df.sort_values("adj_pvalue").reset_index(drop=True)

        n_up = (de_df["significance"] == "UP").sum()
        n_down = (de_df["significance"] == "DOWN").sum()
        logger.info(
            "DE results: %d proteins tested, %d UP, %d DOWN "
            "(|log2FC| > %.1f, FDR < %.2f)",
            len(de_df), n_up, n_down,
            self.log2fc_threshold, self.fdr_threshold,
        )

        return de_df

    def run_comparisons(
        self,
        protein_df: pd.DataFrame,
        sample_to_condition: dict[str, str],
        contrasts: list[tuple[str, str]],
    ) -> dict[str, pd.DataFrame]:
        """
        Run DE analysis for multiple contrasts.

        Parameters
        ----------
        protein_df : pd.DataFrame
            Wide-format protein matrix.
        sample_to_condition : dict
            Sample-to-condition mapping.
        contrasts : list[tuple[str, str]]
            List of contrast pairs.

        Returns
        -------
        dict[str, pd.DataFrame]
            Results keyed by contrast name (e.g., "NASH-HL").
        """
        results = {}
        for contrast in contrasts:
            key = f"{contrast[0]}-{contrast[1]}"
            results[key] = self.run(protein_df, sample_to_condition, contrast)
        return results


def _split_samples(
    sample_cols: list[str],
    sample_to_condition: dict[str, str],
    cond_a: str,
    cond_b: str,
) -> tuple[list[str], list[str]]:
    """Split sample columns into two groups by condition label."""
    sa = [s for s in sample_cols if sample_to_condition.get(s) == cond_a]
    sb = [s for s in sample_cols if sample_to_condition.get(s) == cond_b]
    avail = sorted(set(sample_to_condition.values()))
    if not sa:
        raise ValueError(f"No samples for '{cond_a}'. Available: {avail}")
    if not sb:
        raise ValueError(f"No samples for '{cond_b}'. Available: {avail}")
    return sa, sb


def _log2_transform(
    intensity: pd.DataFrame, skip: bool,
) -> pd.DataFrame:
    """Optionally log2-transform an intensity matrix."""
    if skip:
        return intensity
    return np.log2(intensity.replace(0, np.nan))


def _collect_protein_stats(
    log2_matrix: pd.DataFrame,
    samples_a: list[str],
    samples_b: list[str],
    cond_a: str,
    cond_b: str,
) -> pd.DataFrame:
    """Compute per-protein pooled statistics for limma/DEqMS."""
    rows = []
    for protein in log2_matrix.index:
        va = log2_matrix.loc[protein, samples_a].dropna().values
        vb = log2_matrix.loc[protein, samples_b].dropna().values
        if len(va) < 2 or len(vb) < 2:
            continue
        na, nb = len(va), len(vb)
        ma, mb = np.mean(va), np.mean(vb)
        ss = np.sum((va - ma) ** 2) + np.sum((vb - mb) ** 2)
        df_res = (na - 1) + (nb - 1)
        rows.append({
            "ProteinName": protein,
            "log2FC": ma - mb,
            f"mean_{cond_a}": ma,
            f"mean_{cond_b}": mb,
            "n_a": na, "n_b": nb,
            "s2": ss / df_res,
            "df_residual": df_res,
        })
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def _moderated_test(
    stats_df: pd.DataFrame,
    d0: float,
    s0_sq: float,
    cond_a: str,
    cond_b: str,
) -> pd.DataFrame:
    """Compute moderated t-statistics from per-protein stats."""
    results = []
    for _, row in stats_df.iterrows():
        df_res = row["df_residual"]
        s2_mod = (d0 * s0_sq + df_res * row["s2"]) / (d0 + df_res)
        se = np.sqrt(s2_mod * (1.0 / row["n_a"] + 1.0 / row["n_b"]))
        if se == 0:
            continue
        t_stat = row["log2FC"] / se
        pvalue = 2.0 * stats.t.sf(np.abs(t_stat), df=d0 + df_res)
        if np.isnan(pvalue):
            continue
        results.append({
            "ProteinName": row["ProteinName"],
            "log2FC": row["log2FC"],
            "pvalue": pvalue,
            f"mean_{cond_a}": row[f"mean_{cond_a}"],
            f"mean_{cond_b}": row[f"mean_{cond_b}"],
            "n_a": int(row["n_a"]),
            "n_b": int(row["n_b"]),
        })
    return pd.DataFrame(results)


def _ihw_covariate(de_df: pd.DataFrame) -> np.ndarray | None:
    """Extract an informative covariate for IHW from DE results."""
    mean_cols = [c for c in de_df.columns if c.startswith("mean_")]
    if mean_cols:
        return de_df[mean_cols].mean(axis=1).values
    if "n_a" in de_df.columns and "n_b" in de_df.columns:
        return (de_df["n_a"] + de_df["n_b"]).values.astype(float)
    return None


def _ihw_optimize_weights(
    pvalues: np.ndarray,
    bins: np.ndarray,
    unique_bins: np.ndarray,
    alpha: float,
) -> np.ndarray:
    """Iteratively optimise per-bin weights for IHW."""
    n_bins = len(unique_bins)
    weights = np.ones(n_bins)
    for _ in range(10):
        w_exp = np.array([weights[b] for b in bins])
        w_p = pvalues / np.clip(w_exp, 0.01, None)
        rej = np.zeros(n_bins)
        for idx, bval in enumerate(unique_bins):
            _, adj, _, _ = multipletests(w_p[bins == bval], method="fdr_bh", alpha=alpha)
            rej[idx] = (adj < alpha).sum()
        total = rej.sum()
        if total > 0:
            nw = np.clip(rej / total * n_bins, 0.1, 10.0)
            nw = nw / nw.mean()
            if np.allclose(weights, nw, atol=0.01):
                break
            weights = nw
    return weights


def _ihw_correction(
    pvalues: np.ndarray,
    de_df: pd.DataFrame,
    alpha: float = 0.05,
    n_bins: int = 5,
) -> np.ndarray:
    """IHW multiple-testing correction (Ignatiadis et al. 2016)."""
    n = len(pvalues)
    if n == 0:
        return pvalues.copy()

    covariate = _ihw_covariate(de_df)
    if covariate is None:
        logger.warning("IHW: no covariate, falling back to BH")
        return multipletests(pvalues, method="fdr_bh")[1]

    valid = np.isfinite(pvalues) & np.isfinite(covariate)
    if valid.sum() < n_bins * 2:
        return multipletests(pvalues, method="fdr_bh")[1]

    bins = pd.qcut(covariate[valid], q=n_bins, labels=False, duplicates="drop")
    unique_bins = np.unique(bins)
    weights = _ihw_optimize_weights(pvalues[valid], bins, unique_bins, alpha)

    w_exp = np.array([weights[b] for b in bins])
    w_p = np.clip(pvalues[valid] / np.clip(w_exp, 0.01, None), 0, 1)
    _, adj_valid, _, _ = multipletests(w_p, method="fdr_bh")

    adj = np.full(n, np.nan)
    adj[valid] = adj_valid

    n_bh = multipletests(pvalues[valid], method="fdr_bh")[0].sum()
    n_ihw = (adj_valid < alpha).sum()
    logger.info(
        "IHW: %d bins, BH=%d, IHW=%d (gain=%+d)",
        len(unique_bins), n_bh, n_ihw, n_ihw - n_bh,
    )
    return adj


def _fit_f_prior(
    s2: np.ndarray, df: np.ndarray
) -> tuple[float, float]:
    """Estimate limma prior hyperparameters (d0, s0^2) via moment matching."""
    valid = np.isfinite(s2) & (s2 > 0) & np.isfinite(df) & (df > 0)
    s2, df = s2[valid], df[valid]

    if len(s2) < 3:
        logger.warning("Too few proteins for eBayes estimation, using defaults")
        return 3.0, float(np.median(s2))

    unique_df = np.unique(df)
    if len(unique_df) == 1:
        d0, s0_sq = _fit_f_prior_balanced(s2, unique_df[0])
    else:
        d0, s0_sq = _fit_f_prior_unbalanced(s2, df)

    return max(d0, 0.01), max(s0_sq, 1e-15)


def _fit_f_prior_balanced(
    s2: np.ndarray, d: float,
) -> tuple[float, float]:
    """Balanced case: all proteins share the same residual df."""
    z = np.log(s2)
    z_mean, z_var = np.mean(z), np.var(z, ddof=1)
    target = z_var - _trigamma(d / 2.0)
    d0 = 1e6 if target <= 0 else _solve_trigamma(target) * 2.0
    log_s0 = (
        z_mean
        - digamma(d / 2.0) + np.log(d / 2.0)
        + digamma(d0 / 2.0) - np.log(d0 / 2.0)
    )
    return d0, np.exp(log_s0)


def _fit_f_prior_unbalanced(
    s2: np.ndarray, df: np.ndarray,
) -> tuple[float, float]:
    """Unbalanced case: proteins have different residual df."""
    z = np.log(s2)
    e_z = digamma(df / 2.0) - np.log(df / 2.0)
    z_res = z - e_z
    target = np.var(z_res, ddof=1) - np.mean(_trigamma(df / 2.0))
    d0 = 1e6 if target <= 0 else _solve_trigamma(target) * 2.0
    log_s0 = np.mean(z_res) + digamma(d0 / 2.0) - np.log(d0 / 2.0)
    return d0, np.exp(log_s0)


def _solve_trigamma(target: float) -> float:
    """
    Solve trigamma(x) = target for x > 0.

    Uses Newton's method following limma's approach.

    Parameters
    ----------
    target : float
        Target value of trigamma function.

    Returns
    -------
    float
        Solution x such that trigamma(x) ≈ target.
    """
    if target <= 0:
        return 1e10

    # Initial guess from asymptotic expansion: trigamma(x) ~ 1/x + 1/(2x^2)
    # For large x, trigamma(x) ~ 1/x, so x ~ 1/target
    if target > 1e6:
        return 1.0 / target

    # Better starting value
    x = 0.5 + 1.0 / target

    # Newton's method: x_{n+1} = x_n - (trigamma(x_n) - target) / tetragamma(x_n)
    for _ in range(50):
        fx = _trigamma(x) - target
        dfx = _tetragamma(x)
        if abs(dfx) < 1e-20:
            break
        step = fx / dfx
        # Damped Newton to stay positive
        while x - step <= 0:
            step /= 2.0
        x = x - step
        if abs(fx) < 1e-12:
            break

    return max(x, 1e-10)
