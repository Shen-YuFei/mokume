"""
Differential expression analysis for proteomics data.

Provides statistical testing to identify proteins with significant abundance
differences between experimental conditions.
"""

import numpy as np
import pandas as pd

from mokume.analysis.ensemble_config import SUPPORTED_DE_METHODS
from mokume.analysis.deqms import run_deqms
from mokume.analysis.limma import run_limma
from mokume.analysis.limrots import run_limrots
from mokume.analysis.proda import run_proda
from mokume.analysis.rots import run_rots
from mokume.core.logger import get_logger

logger = get_logger("mokume.analysis.de")


_DE_OPTION_DEFAULTS = {
    "log2fc_threshold": 0.5,
    "fdr_threshold": 0.05,
    "skip_log2": False,
    "fdr_method": "bh",
    "n_boot": 100,
    "seed": 42,
}


def _multipletests():
    """Lazily import statsmodels' ``multipletests`` (needs ``mokume[analysis]``)."""
    try:
        from statsmodels.stats.multitest import multipletests
    except ImportError as exc:
        raise ImportError(
            "statsmodels is required for differential-expression FDR correction. "
            "Install it with: pip install mokume[analysis]"
        ) from exc
    return multipletests


class DifferentialExpression:
    """Differential expression analysis for protein intensity data."""

    SUPPORTED_METHODS = SUPPORTED_DE_METHODS

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
            raise TypeError(
                f"Unexpected DifferentialExpression options: {unknown_args}"
            )

        config = {**_DE_OPTION_DEFAULTS, **options}
        self.method = method.lower()
        if self.method not in self.SUPPORTED_METHODS:
            raise ValueError(
                f"Unknown DE method '{method}'. Supported: {self.SUPPORTED_METHODS}"
            )
        self.log2fc_threshold = config["log2fc_threshold"]
        self.fdr_threshold = config["fdr_threshold"]
        self.skip_log2 = config["skip_log2"]
        self.fdr_method = config["fdr_method"].lower()
        self.n_boot = config["n_boot"]
        self.seed = config["seed"]
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
            protein_df,
            sample_to_condition,
            cond_a,
            cond_b,
        )

        dispatch = {
            "limrots": self._run_limrots,
            "deqms": self._run_deqms,
            "proda": self._run_proda,
            "limma": self._run_limma,
            "rots": self._run_rots,
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
            sample_cols,
            sample_to_condition,
            cond_a,
            cond_b,
        )
        logger.info(
            "DE analysis: %s (%d samples) vs %s (%d samples)",
            cond_a,
            len(samples_a),
            cond_b,
            len(samples_b),
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
        """Run LimROTS via :func:`mokume.analysis.limrots.run_limrots`.

        LimROTS already provides a permutation-based FDR in ``adj_pvalue``; keep
        that estimate and only add significance labels, instead of overwriting it
        with BH/IHW like the parametric methods.
        """
        result = run_limrots(
            log2_matrix,
            sample_groups[0],
            sample_groups[1],
            contrast,
            n_boot=self.n_boot,
            seed=self.seed,
        )
        if result.empty:
            return result
        return self._classify_and_sort(result)

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

    def _run_limma(
        self,
        log2_matrix: pd.DataFrame,
        sample_groups: tuple[list[str], list[str]],
        contrast: tuple[str, str],
    ) -> pd.DataFrame:
        """Run limma via :func:`mokume.analysis.limma.run_limma`."""
        cond_a, cond_b = contrast
        result = run_limma(
            log2_matrix, sample_groups[0], sample_groups[1], cond_a, cond_b
        )
        if result.empty:
            return result
        return self._finalize_results(result)

    def _run_rots(
        self,
        log2_matrix: pd.DataFrame,
        sample_groups: tuple[list[str], list[str]],
        contrast: tuple[str, str],
    ) -> pd.DataFrame:
        """Run ROTS via :func:`mokume.analysis.rots.run_rots`.

        Like LimROTS, ROTS provides its own permutation-based FDR in
        ``adj_pvalue``; preserve it rather than re-adjusting with BH/IHW.
        """
        cond_a, cond_b = contrast
        result = run_rots(
            log2_matrix,
            sample_groups[0],
            sample_groups[1],
            cond_a,
            cond_b,
            n_boot=self.n_boot,
            seed=self.seed,
        )
        if result.empty:
            return result
        return self._classify_and_sort(result)

    def _finalize_results(self, de_df: pd.DataFrame) -> pd.DataFrame:
        """Apply BH/IHW FDR correction over ``pvalue``, then classify and sort.

        Methods that already carry their own permutation-based FDR in
        ``adj_pvalue`` (LimROTS, ROTS) call :meth:`_classify_and_sort` directly
        instead, so their method-specific FDR is preserved rather than
        overwritten with BH/IHW.
        """
        pvalues = de_df["pvalue"].values
        valid = np.isfinite(pvalues)
        adj = np.full(len(pvalues), np.nan)
        if valid.any():
            if self.fdr_method == "ihw":
                adj[valid] = _ihw_correction(
                    pvalues[valid],
                    de_df.loc[valid],
                    alpha=self.fdr_threshold,
                )
            else:
                adj[valid] = _multipletests()(pvalues[valid], method="fdr_bh")[1]
        de_df["adj_pvalue"] = adj
        return self._classify_and_sort(de_df)

    def _classify_and_sort(self, de_df: pd.DataFrame) -> pd.DataFrame:
        """Label UP/DOWN/Unchanged from ``adj_pvalue``, sort, and log.

        Operates on the FDR estimate already present in ``adj_pvalue`` (whether
        written by :meth:`_finalize_results` or by a method's own permutation
        FDR), so it never re-adjusts p-values.
        """
        # A protein that is absent (or imputed to a degenerate/extreme value)
        # in one group can yield a non-finite log2 fold change, e.g. when an
        # imputer overflows the intensity to +inf or produces a non-positive
        # value. Such a fold change is undefined, so coerce it to NaN: this
        # keeps it out of the significance classification below (NaN
        # comparisons are False) instead of letting an Inf silently pass the
        # ``|log2FC| > threshold`` filter and be counted as differential.
        if "log2FC" in de_df.columns:
            nonfinite = ~np.isfinite(de_df["log2FC"].to_numpy(dtype=float))
            if nonfinite.any():
                logger.warning(
                    "DE results: %d proteins had a non-finite log2FC "
                    "(coerced to NaN, excluded from significance)",
                    int(nonfinite.sum()),
                )
                de_df.loc[nonfinite, "log2FC"] = np.nan

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
            len(de_df),
            n_up,
            n_down,
            self.log2fc_threshold,
            self.fdr_threshold,
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
    intensity: pd.DataFrame,
    skip: bool,
) -> pd.DataFrame:
    """Optionally log2-transform an intensity matrix."""
    if skip:
        return intensity
    return np.log2(intensity.replace(0, np.nan))


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
            _, adj, _, _ = _multipletests()(
                w_p[bins == bval], method="fdr_bh", alpha=alpha
            )
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
    return _multipletests()(pvalues, method="fdr_bh")[1]


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
    weighted_p = np.clip(
        pvalues[valid] / np.clip(np.array([weights[b] for b in bins]), 0.01, None), 0, 1
    )
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
        len(unique_bins),
        _multipletests()(pvalues[valid], method="fdr_bh")[0].sum(),
        n_ihw,
        n_ihw - _multipletests()(pvalues[valid], method="fdr_bh")[0].sum(),
    )
    return adj
