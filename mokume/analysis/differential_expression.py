"""Differential expression analysis for proteomics data.

Provides statistical testing to identify proteins with significant abundance
differences between experimental conditions.
"""

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from mokume.analysis import _shared_stats
from mokume.analysis.deqms import run_deqms
from mokume.analysis.limrots import run_limrots
from mokume.analysis.proda import run_proda
from mokume.core.logger import get_logger

logger = get_logger("mokume.analysis.de")


_DE_OPTION_DEFAULTS = {
    "log2fc_threshold": 0.5,
    "fdr_threshold": 0.05,
    "skip_log2": False,
    "fdr_method": "bh",
    "n_boot": 100,
    "n_threads": None,
}


def _trigamma(x):
    """Compute the trigamma function (second derivative of log-gamma)."""
    return _shared_stats._trigamma(x)


def _tetragamma(x):
    """Compute the tetragamma function (third derivative of log-gamma)."""
    return _shared_stats._tetragamma(x)


class DifferentialExpression:
    """Differential expression analysis for protein intensity data."""

    SUPPORTED_METHODS = ("limrots", "deqms", "proda")

    def __init__(
        self,
        method: str = "limrots",
        peptide_counts: pd.Series | None = None,
        **options,
    ):
        """Initialize differential expression analysis settings."""
        unknown = set(options) - set(_DE_OPTION_DEFAULTS)
        if unknown:
            unknown_args = ", ".join(sorted(unknown))
            raise TypeError(f"Unexpected DifferentialExpression options: {unknown_args}")

        config = {**_DE_OPTION_DEFAULTS, **options}
        self.method = method.lower()
        if self.method not in self.SUPPORTED_METHODS:
            raise ValueError(
                f"Unknown DE method '{method}'. "
                f"Supported: {self.SUPPORTED_METHODS}"
            )
        self.log2fc_threshold = config["log2fc_threshold"]
        self.fdr_threshold = config["fdr_threshold"]
        self.skip_log2 = config["skip_log2"]
        self.fdr_method = config["fdr_method"].lower()
        self.n_boot = config["n_boot"]
        self.n_threads = config["n_threads"]
        self.peptide_counts = peptide_counts

    def run(
        self,
        protein_df: pd.DataFrame,
        sample_to_condition: dict[str, str],
        contrast: tuple[str, str],
    ) -> pd.DataFrame:
        """Run differential expression analysis for a single contrast."""
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
        return runner(log2_matrix, (samples_a, samples_b), contrast)

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
        sample_groups: tuple[list[str], list[str]],
        contrast: tuple[str, str],
    ) -> pd.DataFrame:
        """Run LimROTS via :func:`mokume.analysis.limrots.run_limrots`."""
        result = run_limrots(
            log2_matrix,
            sample_groups[0],
            sample_groups[1],
            contrast,
            n_boot=self.n_boot,
            n_threads=self.n_threads,
        )
        if result.empty:
            return result
        return self._finalize_results(result)

    def _run_deqms(
        self,
        log2_matrix: pd.DataFrame,
        sample_groups: tuple[list[str], list[str]],
        contrast: tuple[str, str],
    ) -> pd.DataFrame:
        """Run DEqMS via :func:`mokume.analysis.deqms.run_deqms`."""
        result = run_deqms(
            log2_matrix,
            sample_groups[0],
            sample_groups[1],
            contrast,
            peptide_counts=self.peptide_counts,
        )
        if result.empty:
            return result
        return self._finalize_results(result)

    def _run_proda(
        self,
        log2_matrix: pd.DataFrame,
        sample_groups: tuple[list[str], list[str]],
        contrast: tuple[str, str],
    ) -> pd.DataFrame:
        """Run proDA via :func:`mokume.analysis.proda.run_proda`."""
        cond_a, cond_b = contrast
        result = run_proda(
            log2_matrix, sample_groups[0], sample_groups[1], cond_a, cond_b
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
        """Run differential expression analysis for multiple contrasts."""
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
    """Optimise per-bin weights for IHW."""
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


def _bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    """Apply Benjamini-Hochberg correction."""
    return multipletests(pvalues, method="fdr_bh")[1]


def _ihw_adjust_valid(
    pvalues: np.ndarray,
    covariate: np.ndarray,
    alpha: float,
    n_bins: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """Apply IHW to the subset with finite p-values and covariates."""
    valid = np.isfinite(pvalues) & np.isfinite(covariate)
    if valid.sum() < n_bins * 2:
        return None

    bins = pd.qcut(covariate[valid], q=n_bins, labels=False, duplicates="drop")
    unique_bins = np.unique(bins)
    weights = _ihw_optimize_weights(pvalues[valid], bins, unique_bins, alpha)
    weighted_p = np.clip(pvalues[valid] / np.clip(np.array([weights[b] for b in bins]), 0.01, None), 0, 1)
    return valid, unique_bins, _bh_adjust(weighted_p)


def _ihw_correction(
    pvalues: np.ndarray,
    de_df: pd.DataFrame,
    alpha: float = 0.05,
    n_bins: int = 5,
) -> np.ndarray:
    """IHW multiple-testing correction (Ignatiadis et al. 2016)."""
    if len(pvalues) == 0:
        return pvalues.copy()

    covariate = _ihw_covariate(de_df)
    if covariate is None:
        logger.warning("IHW: no covariate, falling back to BH")
        return _bh_adjust(pvalues)

    ihw_result = _ihw_adjust_valid(pvalues, covariate, alpha, n_bins)
    if ihw_result is None:
        return _bh_adjust(pvalues)

    valid, unique_bins, adj_valid = ihw_result
    adj = np.full(len(pvalues), np.nan)
    adj[valid] = adj_valid

    n_ihw = (adj_valid < alpha).sum()
    logger.info(
        "IHW: %d bins, BH=%d, IHW=%d (gain=%+d)",
        len(unique_bins), multipletests(pvalues[valid], method="fdr_bh")[0].sum(), n_ihw,
        n_ihw - multipletests(pvalues[valid], method="fdr_bh")[0].sum(),
    )
    return adj
