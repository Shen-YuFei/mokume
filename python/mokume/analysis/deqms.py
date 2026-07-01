"""
DEqMS-style differential expression analysis (pure Python).

Reimplements the DEqMS pipeline: limma (lmFit → contrasts.fit → eBayes) plus
spectraCounteBayes variance moderation, without rpy2 or R dependencies.

Reference
---------
Zhu Y, Orre LM, Zhou Tran Y, et al. DEqMS: a method for accurate variance
estimation in differential protein expression analysis. *Mol Cell Proteomics*.
2020;19(6):1047-1057.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.special import digamma, polygamma
from scipy.stats import t as t_dist

from mokume.analysis._helpers import (
    bh_adjust,
    filter_testable,
    per_group_summary,
)
from mokume.analysis.limma import (
    _contrasts_fit,
    _ebayes,
    _lm_fit_with_na,
    EBayesResult,
)
from mokume.core.logger import get_logger

logger = get_logger("mokume.analysis.deqms")


# ---------------------------------------------------------------------------
# spectraCounteBayes: count-aware variance moderation (Zhu et al. 2020)
# ---------------------------------------------------------------------------


def _spectra_count_ebayes(
    fit: EBayesResult,
    counts: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Apply spectra-count-aware Bayes moderation (R DEqMS::spectraCounteBayes).

    Returns (sca_t, sca_pvalue, sca_var_post, d0).
    """
    sigma2 = fit.sigma**2
    log_var = np.log(sigma2)
    df = fit.df_residual.copy()
    n_genes = len(log_var)

    df_valid = df.copy()
    df_valid[df_valid == 0] = np.nan

    counts_f = counts.astype(float).copy()
    if np.nanmin(counts_f) == 0:
        counts_f = counts_f + 1

    x = np.log2(counts_f)

    valid = np.isfinite(log_var) & np.isfinite(x) & np.isfinite(df_valid)
    if valid.sum() < 10:
        logger.warning(
            "spectraCounteBayes: too few valid points (%d), "
            "falling back to standard eBayes",
            valid.sum(),
        )
        return fit.t_stat[:, 0], fit.p_value[:, 0], fit.s2_post, fit.df_prior

    x_valid = x[valid]
    x_range = float(np.ptp(x_valid))
    if x_range < 1e-10:
        logger.warning(
            "spectraCounteBayes: all counts identical, falling back to standard eBayes",
        )
        return fit.t_stat[:, 0], fit.p_value[:, 0], fit.s2_post, fit.df_prior

    try:
        from statsmodels.nonparametric.smoothers_lowess import lowess
    except ImportError as exc:
        raise ImportError(
            "statsmodels is required for the DEqMS method. "
            "Install it with: pip install mokume[analysis]"
        ) from exc

    try:
        fitted = lowess(
            log_var[valid],
            x_valid,
            frac=0.75,
            return_sorted=False,
        )
    except Exception:
        logger.warning("LOESS failed, falling back to standard eBayes")
        return fit.t_stat[:, 0], fit.p_value[:, 0], fit.s2_post, fit.df_prior

    y_pred = np.full(n_genes, np.nan)
    y_pred[valid] = fitted

    non_valid = ~valid & np.isfinite(x)
    if non_valid.any():
        # statsmodels.lowess does not predict at new points, so linearly
        # interpolate the fitted curve at the unseen x values.
        order = np.argsort(x_valid)
        y_pred[non_valid] = np.interp(x[non_valid], x_valid[order], fitted[order])

    eg = log_var - digamma(df_valid / 2) + np.log(df_valid / 2)
    egpred = y_pred - digamma(df_valid / 2) + np.log(df_valid / 2)

    myfct = (eg - egpred) ** 2 - polygamma(1, df_valid / 2)
    mean_myfct = float(np.nanmean(myfct))

    d0 = _grid_search_d0(mean_myfct, n_genes)

    s02 = np.exp(egpred + digamma(d0 / 2) - np.log(d0 / 2))

    post_var = (d0 * s02 + df * sigma2) / (d0 + df)
    post_df = d0 + df

    sca_t = fit.coefficients[:, 0] / (fit.stdev_unscaled[:, 0] * np.sqrt(post_var))
    sca_pvalue = 2.0 * t_dist.sf(np.abs(sca_t), post_df)

    return sca_t, sca_pvalue, post_var, d0


def _grid_search_d0(mean_myfct: float, n_genes: int) -> float:
    """Grid search for prior df d0 (R DEqMS internal algorithm)."""
    max_iter = n_genes * 10
    best_d0 = 0.1
    best_diff = np.inf
    prev_prev = np.inf
    prev = np.inf
    for i in range(1, max_iter + 1):
        test_d0 = i / 10.0
        diff = abs(mean_myfct - float(polygamma(1, test_d0 / 2)))
        if diff < best_diff:
            best_diff = diff
            best_d0 = test_d0
        if i > 2 and prev_prev < prev:
            break
        prev_prev = prev
        prev = diff
    return best_d0


# ---------------------------------------------------------------------------
# Count vector builder (shared logic with deqms.py)
# ---------------------------------------------------------------------------


def _build_count_vector(
    proteins: list[str],
    peptide_counts: pd.Series | None,
) -> np.ndarray:
    """Align peptide counts with the matrix row order; default to 1."""
    counts = np.ones(len(proteins), dtype=int)
    if peptide_counts is None:
        return counts
    aligned = peptide_counts.reindex(proteins)
    valid = aligned.notna()
    counts[valid.values] = aligned[valid].astype(int).values
    return counts


# ---------------------------------------------------------------------------
# Public API: run_deqms
# ---------------------------------------------------------------------------


def run_deqms(
    log2_matrix: pd.DataFrame,
    samples_a: list[str],
    samples_b: list[str],
    contrast: tuple[str, str],
    **options,
) -> pd.DataFrame:
    """Run DEqMS differential expression (pure Python implementation)."""
    peptide_counts = options.pop("peptide_counts", None)
    if options:
        unknown = ", ".join(sorted(options))
        raise TypeError(f"Unexpected DEqMS options: {unknown}")

    cond_a, cond_b = contrast
    sub_matrix = filter_testable(log2_matrix, samples_a, samples_b)
    if sub_matrix.empty:
        logger.warning("No proteins passed DE filtering for DEqMS")
        return pd.DataFrame()

    sample_order = list(samples_a) + list(samples_b)
    mat = sub_matrix[sample_order].values
    gene_names = list(sub_matrix.index)

    design = np.zeros((len(sample_order), 2))
    design[: len(samples_a), 0] = 1.0
    design[len(samples_a) :, 1] = 1.0
    coef_names = [cond_a, cond_b]

    fit = _lm_fit_with_na(mat, design, gene_names, coef_names)
    contrasts = np.array([[1.0], [-1.0]])
    fit_c = _contrasts_fit(fit, contrasts)
    fit_eb = _ebayes(fit_c)

    counts = _build_count_vector(gene_names, peptide_counts)

    sca_t, sca_pvalue, _, _ = _spectra_count_ebayes(fit_eb, counts)

    raw = pd.DataFrame(
        {
            "ProteinName": gene_names,
            "log2FC": fit_eb.coefficients[:, 0],
            "sca_t": sca_t,
            "sca_pvalue": sca_pvalue,
            "peptide_count": counts,
        }
    )

    summary = per_group_summary(sub_matrix, samples_a, samples_b, cond_a, cond_b)
    merged = raw.merge(summary, on="ProteinName", how="left")
    merged["pvalue"] = merged["sca_pvalue"]
    merged["adj_pvalue"] = bh_adjust(merged["sca_pvalue"].values)
    merged["sca_adj_pvalue"] = merged["adj_pvalue"]

    columns = [
        "ProteinName",
        "log2FC",
        "pvalue",
        "adj_pvalue",
        "sca_t",
        "sca_pvalue",
        "sca_adj_pvalue",
        f"mean_{cond_a}",
        f"mean_{cond_b}",
        "n_a",
        "n_b",
        "peptide_count",
    ]
    available = [c for c in columns if c in merged.columns]
    return merged[available].sort_values("adj_pvalue").reset_index(drop=True)
