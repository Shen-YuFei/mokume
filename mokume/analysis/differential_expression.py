"""
Differential expression analysis for proteomics data.

Provides statistical testing to identify proteins with significant abundance
differences between experimental conditions.
"""

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from mokume.core.logger import get_logger

logger = get_logger("mokume.analysis.de")


class DifferentialExpression:
    """
    Differential expression analysis for protein intensity data.

    Parameters
    ----------
    method : str
        Statistical test method: ``"ttest"`` (Welch's t-test) or
        ``"limma"`` (moderated t-test with empirical Bayes variance shrinkage).
    log2fc_threshold : float
        Minimum absolute log2 fold change for significance.
    fdr_threshold : float
        Maximum adjusted p-value (FDR) for significance.
    skip_log2 : bool
        If True, skip the log2 transformation (data is already in log2 space,
        e.g., from ratio-based quantification).
    """

    def __init__(
        self,
        method: str = "ttest",
        log2fc_threshold: float = 0.5,
        fdr_threshold: float = 0.05,
        skip_log2: bool = False,
    ):
        self.method = method.lower()
        self.log2fc_threshold = log2fc_threshold
        self.fdr_threshold = fdr_threshold
        self.skip_log2 = skip_log2

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
        protein_col = protein_df.columns[0]
        sample_cols = [c for c in protein_df.columns if c != protein_col]

        # Split samples by condition
        samples_a = [s for s in sample_cols if sample_to_condition.get(s) == cond_a]
        samples_b = [s for s in sample_cols if sample_to_condition.get(s) == cond_b]

        if not samples_a:
            raise ValueError(
                f"No samples found for condition '{cond_a}'. "
                f"Available conditions: {sorted(set(sample_to_condition.values()))}"
            )
        if not samples_b:
            raise ValueError(
                f"No samples found for condition '{cond_b}'. "
                f"Available conditions: {sorted(set(sample_to_condition.values()))}"
            )

        logger.info(
            f"DE analysis: {cond_a} ({len(samples_a)} samples) vs "
            f"{cond_b} ({len(samples_b)} samples)"
        )

        intensity_matrix = protein_df.set_index(protein_col)[sample_cols]

        # Log2 transform (handle zeros/negatives) — skip if data is already log2
        if not self.skip_log2:
            intensity_matrix = intensity_matrix.replace(0, np.nan)
        log2_matrix = intensity_matrix if self.skip_log2 else np.log2(intensity_matrix)

        if self.method == "limma":
            return self._run_limma(
                log2_matrix, samples_a, samples_b, cond_a, cond_b
            )

        # Standard per-protein t-test
        results = []
        for protein in log2_matrix.index:
            vals_a = log2_matrix.loc[protein, samples_a].dropna().values
            vals_b = log2_matrix.loc[protein, samples_b].dropna().values

            if len(vals_a) < 2 or len(vals_b) < 2:
                continue

            log2fc = np.mean(vals_a) - np.mean(vals_b)
            stat, pvalue = stats.ttest_ind(vals_a, vals_b, equal_var=False)

            if np.isnan(pvalue):
                continue

            results.append({
                "ProteinName": protein,
                "log2FC": log2fc,
                "pvalue": pvalue,
                f"mean_{cond_a}": np.mean(vals_a),
                f"mean_{cond_b}": np.mean(vals_b),
                "n_a": len(vals_a),
                "n_b": len(vals_b),
            })

        if not results:
            logger.warning("No proteins passed DE filtering")
            return pd.DataFrame(
                columns=["ProteinName", "log2FC", "pvalue", "adj_pvalue", "significance"]
            )

        de_df = pd.DataFrame(results)
        return self._finalize_results(de_df)

    def _run_limma(
        self,
        log2_matrix: pd.DataFrame,
        samples_a: list[str],
        samples_b: list[str],
        cond_a: str,
        cond_b: str,
    ) -> pd.DataFrame:
        """
        Run limma-style moderated t-test with proper empirical Bayes.

        Implements the eBayes procedure from Smyth (2004):
        1. Fit ordinary linear model per protein (compute residual variances)
        2. Estimate prior hyperparameters (d0, s0^2) from the distribution
           of residual variances across all proteins
        3. Compute moderated variances and t-statistics
        """
        # Step 1: Compute per-protein statistics
        protein_stats = []
        for protein in log2_matrix.index:
            vals_a = log2_matrix.loc[protein, samples_a].dropna().values
            vals_b = log2_matrix.loc[protein, samples_b].dropna().values

            if len(vals_a) < 2 or len(vals_b) < 2:
                continue

            n_a, n_b = len(vals_a), len(vals_b)
            mean_a, mean_b = np.mean(vals_a), np.mean(vals_b)
            log2fc = mean_a - mean_b

            # Residual variance (pooled, equal variance assumption like lmFit)
            ss_a = np.sum((vals_a - mean_a) ** 2)
            ss_b = np.sum((vals_b - mean_b) ** 2)
            df_residual = (n_a - 1) + (n_b - 1)
            s2 = (ss_a + ss_b) / df_residual

            protein_stats.append({
                "ProteinName": protein,
                "log2FC": log2fc,
                f"mean_{cond_a}": mean_a,
                f"mean_{cond_b}": mean_b,
                "n_a": n_a,
                "n_b": n_b,
                "s2": s2,
                "df_residual": df_residual,
            })

        if not protein_stats:
            logger.warning("No proteins passed DE filtering")
            return pd.DataFrame(
                columns=["ProteinName", "log2FC", "pvalue", "adj_pvalue", "significance"]
            )

        stats_df = pd.DataFrame(protein_stats)

        # Step 2: Estimate prior hyperparameters via empirical Bayes
        s2_values = stats_df["s2"].values
        df_values = stats_df["df_residual"].values

        d0, s0_sq = _fit_f_prior(s2_values, df_values)
        logger.info(f"eBayes prior: d0={d0:.3f}, s0^2={s0_sq:.6f}")

        # Step 3: Compute moderated statistics
        results = []
        for _, row in stats_df.iterrows():
            df_res = row["df_residual"]
            s2 = row["s2"]

            # Moderated variance: shrink toward prior
            s2_mod = (d0 * s0_sq + df_res * s2) / (d0 + df_res)

            # Moderated t-statistic
            se = np.sqrt(s2_mod * (1.0 / row["n_a"] + 1.0 / row["n_b"]))
            if se == 0:
                continue

            t_stat = row["log2FC"] / se
            df_mod = d0 + df_res
            pvalue = 2.0 * stats.t.sf(np.abs(t_stat), df=df_mod)

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

        de_df = pd.DataFrame(results)
        return self._finalize_results(de_df)

    def _finalize_results(self, de_df: pd.DataFrame) -> pd.DataFrame:
        """Apply FDR correction and classify significance."""
        # FDR correction (Benjamini-Hochberg)
        reject, adj_pvalues, _, _ = multipletests(
            de_df["pvalue"].values, method="fdr_bh"
        )
        de_df["adj_pvalue"] = adj_pvalues

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
            f"DE results: {len(de_df)} proteins tested, "
            f"{n_up} UP, {n_down} DOWN "
            f"(|log2FC| > {self.log2fc_threshold}, FDR < {self.fdr_threshold})"
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


def _fit_f_prior(
    s2: np.ndarray, df: np.ndarray
) -> tuple[float, float]:
    """
    Estimate prior hyperparameters (d0, s0^2) for the empirical Bayes
    moderated t-test, following the limma approach (Smyth, 2004).

    Under the limma model, the residual variances s^2_g follow a scaled
    inverse chi-squared distribution:

        s^2_g | sigma^2_g ~ (sigma^2_g / d_g) * chi^2(d_g)

    with a prior:

        1/sigma^2_g ~ (1 / (d0 * s0^2)) * chi^2(d0)

    This results in the marginal distribution:

        s^2_g ~ s0^2 * F(d_g, d0)

    The method estimates d0 and s0^2 by matching the first two moments
    of log(s^2) to the expected moments of log(F) distribution (the
    method used in limma::fitFDist).

    Parameters
    ----------
    s2 : np.ndarray
        Per-protein residual variances.
    df : np.ndarray
        Per-protein residual degrees of freedom.

    Returns
    -------
    tuple[float, float]
        (d0, s0_sq): prior degrees of freedom and prior variance.
    """
    # Filter out zeros/inf
    valid = np.isfinite(s2) & (s2 > 0) & np.isfinite(df) & (df > 0)
    s2 = s2[valid]
    df = df[valid]

    if len(s2) < 3:
        logger.warning("Too few proteins for eBayes estimation, using defaults")
        return 3.0, float(np.median(s2))

    # Use the method of moments on log(s2) following limma::fitFDist
    # E[log(s2)] = log(s0^2) + digamma(df/2) - log(df/2) - digamma(d0/2) + log(d0/2)
    # Var[log(s2)] = trigamma(df/2) + trigamma(d0/2)

    from scipy.special import digamma, polygamma

    # trigamma = polygamma(1, x)
    def trigamma(x):
        return polygamma(1, x)

    # If all df are the same (common case for balanced designs):
    unique_df = np.unique(df)

    if len(unique_df) == 1:
        d = unique_df[0]
        # Mean and variance of log(s2)
        z = np.log(s2)
        z_mean = np.mean(z)
        z_var = np.var(z, ddof=1)

        # Var[log(s2)] = trigamma(d/2) + trigamma(d0/2)
        # => trigamma(d0/2) = z_var - trigamma(d/2)
        target_trigamma = z_var - trigamma(d / 2.0)

        if target_trigamma <= 0:
            # Very large d0 (essentially infinite prior df)
            d0 = 1e6
        else:
            # Solve trigamma(d0/2) = target_trigamma for d0
            d0 = _solve_trigamma(target_trigamma) * 2.0

        # E[log(s2)] = log(s0^2) + digamma(d/2) - log(d/2) - digamma(d0/2) + log(d0/2)
        # => log(s0^2) = z_mean - digamma(d/2) + log(d/2) + digamma(d0/2) - log(d0/2)
        log_s0_sq = (
            z_mean
            - digamma(d / 2.0) + np.log(d / 2.0)
            + digamma(d0 / 2.0) - np.log(d0 / 2.0)
        )
        s0_sq = np.exp(log_s0_sq)

    else:
        # Unbalanced case: different df per protein
        # Use weighted approach
        z = np.log(s2)
        e_z = digamma(df / 2.0) - np.log(df / 2.0)  # E[log(chi2/df)] for each protein

        # Estimate d0 by matching variance
        z_residual = z - e_z
        z_var = np.var(z_residual, ddof=1)
        mean_trigamma_d = np.mean(trigamma(df / 2.0))
        target_trigamma = z_var - mean_trigamma_d

        if target_trigamma <= 0:
            d0 = 1e6
        else:
            d0 = _solve_trigamma(target_trigamma) * 2.0

        # Estimate s0^2
        log_s0_sq = np.mean(z_residual) + digamma(d0 / 2.0) - np.log(d0 / 2.0)
        s0_sq = np.exp(log_s0_sq)

    # Ensure reasonable bounds
    d0 = max(d0, 0.01)
    s0_sq = max(s0_sq, 1e-15)

    return d0, s0_sq


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
    from scipy.special import polygamma

    def trigamma(x):
        return polygamma(1, x)

    def tetragamma(x):
        return polygamma(2, x)

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
        fx = trigamma(x) - target
        dfx = tetragamma(x)
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
