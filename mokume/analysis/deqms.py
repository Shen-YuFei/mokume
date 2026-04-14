"""
DEqMS-style differential expression analysis.

Implements the spectraCounteBayes method from Zhu et al. (2020) which extends
limma's empirical Bayes moderated t-test by incorporating peptide/spectra
counts to model the variance-count relationship.

Reference
---------
Zhu Y, Orre LM, Zhou Tran Y, et al. DEqMS: a method for accurate variance
estimation in differential protein expression analysis. *Mol Cell Proteomics*.
2020;19(6):1047-1057.

"""

import warnings
from typing import NamedTuple

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import digamma, polygamma
from statsmodels.nonparametric.smoothers_lowess import lowess
from statsmodels.stats.multitest import multipletests

from mokume.core.logger import get_logger

logger = get_logger("mokume.analysis.deqms")


class _DEqMSInputs(NamedTuple):
    log2fc: np.ndarray
    sigma2: np.ndarray
    df_residual: np.ndarray
    peptide_counts: np.ndarray
    stdev_unscaled: np.ndarray


def _trigamma(x):
    """Compute the trigamma function (second derivative of log-gamma)."""
    return polygamma(1, x)


def _estimate_d0_grid(mean_sq_residual, max_genes):
    """Estimate prior degrees of freedom d0 by grid search."""
    best_d0 = 1.0
    best_diff = np.inf
    prev_diff = np.inf
    prev_prev_diff = np.inf

    for i in range(1, max_genes * 10 + 1):
        test_d0 = i / 10.0
        current_diff = abs(mean_sq_residual - _trigamma(test_d0 / 2.0))

        if i > 2 and prev_prev_diff < prev_diff:
            break

        if current_diff < best_diff:
            best_diff = current_diff
            best_d0 = test_d0

        prev_prev_diff = prev_diff
        prev_diff = current_diff

    return best_d0


def _constant_variance_fit(log_var: np.ndarray, valid: np.ndarray) -> np.ndarray:
    y_pred = np.full_like(log_var, np.nan)
    if valid.any():
        y_pred[valid] = np.nanmean(log_var[valid])
    return y_pred


def _fit_variance_curve(
    log_var: np.ndarray,
    peptide_counts: np.ndarray,
    fit_method: str,
) -> np.ndarray:
    """Fit log-variance ~ log2(count) via LOESS and return predictions."""
    counts = peptide_counts.copy().astype(float)
    if np.isfinite(counts).any() and np.nanmin(counts) == 0:
        logger.info("Minimum peptide count is 0, adding pseudocount 1")
        counts = counts + 1
    if fit_method != "loess":
        raise ValueError(f"Unsupported fit_method: {fit_method}")
    x = np.log2(counts)
    valid = np.isfinite(x) & np.isfinite(log_var)
    if valid.sum() < 3 or np.unique(x[valid]).size < 2:
        logger.info("DEqMS LOESS fallback to constant variance fit")
        return _constant_variance_fit(log_var, valid)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            fitted = lowess(log_var[valid], x[valid], frac=0.75, return_sorted=False)
    except RuntimeWarning as exc:
        logger.warning("DEqMS LOESS unstable, using constant variance fit: %s", exc)
        return _constant_variance_fit(log_var, valid)
    y_pred = np.full_like(log_var, np.nan)
    y_pred[valid] = fitted
    return y_pred


def _build_deqms_inputs(stats_df: pd.DataFrame) -> _DEqMSInputs:
    return _DEqMSInputs(
        log2fc=stats_df["log2FC"].values,
        sigma2=stats_df["s2"].values,
        df_residual=stats_df["df_residual"].values,
        peptide_counts=stats_df["peptide_count"].values,
        stdev_unscaled=stats_df["stdev_unscaled"].values,
    )


def spectra_count_ebayes(
    inputs: _DEqMSInputs,
    fit_method: str = "loess",
) -> dict:
    """Apply spectra-count-adjusted eBayes moderation with DEqMS."""
    log_var = np.log(inputs.sigma2)
    df_safe = inputs.df_residual.copy().astype(float)
    df_safe[df_safe == 0] = np.nan

    eg = log_var - digamma(df_safe / 2.0) + np.log(df_safe / 2.0)
    y_pred = _fit_variance_curve(log_var, inputs.peptide_counts, fit_method)
    eg_pred = y_pred - digamma(df_safe / 2.0) + np.log(df_safe / 2.0)

    myfct = (eg - eg_pred) ** 2 - _trigamma(df_safe / 2.0)
    d0 = _estimate_d0_grid(np.nanmean(myfct), int(np.sum(df_safe > 0)))
    logger.info("DEqMS prior d0 = %.2f", d0)

    s02 = np.exp(eg_pred + digamma(d0 / 2.0) - np.log(d0 / 2.0))
    post_var = (d0 * s02 + df_safe * inputs.sigma2) / (d0 + df_safe)
    post_df = d0 + df_safe
    sca_t = inputs.log2fc / (inputs.stdev_unscaled * np.sqrt(post_var))
    sca_p = 2.0 * stats.t.sf(np.abs(sca_t), df=post_df)

    return {
        "sca_t": sca_t, "sca_p": sca_p,
        "post_var": post_var, "prior_var": s02,
        "d0": d0, "post_df": post_df,
    }


def _collect_deqms_stats(
    log2_matrix: pd.DataFrame,
    sample_groups: tuple[list[str], list[str]],
    contrast: tuple[str, str],
    peptide_counts: pd.Series | None,
) -> pd.DataFrame:
    """Compute per-protein stats with peptide counts for DEqMS."""
    samples_a, samples_b = sample_groups
    cond_a, cond_b = contrast
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
        count = 1
        if peptide_counts is not None and protein in peptide_counts.index:
            count = peptide_counts.loc[protein]
        rows.append({
            "ProteinName": protein, "log2FC": ma - mb,
            f"mean_{cond_a}": ma, f"mean_{cond_b}": mb,
            "n_a": na, "n_b": nb,
            "s2": ss / df_res, "df_residual": df_res,
            "stdev_unscaled": np.sqrt(1.0 / na + 1.0 / nb),
            "peptide_count": count,
        })
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def _finalize_deqms(
    stats_df: pd.DataFrame, result: dict, cond_a: str, cond_b: str,
) -> pd.DataFrame:
    """Attach eBayes results and apply BH correction."""
    for key in ("sca_t", "post_var", "prior_var", "d0", "post_df"):
        stats_df[key] = result[key]
    stats_df["sca_pvalue"] = result["sca_p"]

    valid = np.isfinite(stats_df["sca_pvalue"])
    adj = np.full(len(stats_df), np.nan)
    if valid.any():
        adj[valid] = multipletests(stats_df.loc[valid, "sca_pvalue"].values, method="fdr_bh")[1]
    stats_df["sca_adj_pvalue"] = adj
    stats_df["pvalue"] = stats_df["sca_pvalue"]
    stats_df["adj_pvalue"] = stats_df["sca_adj_pvalue"]

    cols = [
        "ProteinName", "log2FC", "pvalue", "adj_pvalue",
        "sca_t", "sca_pvalue", "sca_adj_pvalue",
        f"mean_{cond_a}", f"mean_{cond_b}",
        "n_a", "n_b", "peptide_count",
        "post_var", "prior_var", "d0", "post_df",
    ]
    return stats_df[[c for c in cols if c in stats_df.columns]].sort_values(
        "adj_pvalue"
    ).reset_index(drop=True)


def run_deqms(
    log2_matrix: pd.DataFrame,
    samples_a: list[str],
    samples_b: list[str],
    contrast: tuple[str, str],
    **options,
) -> pd.DataFrame:
    """Run DEqMS differential expression analysis."""
    peptide_counts = options.pop("peptide_counts", None)
    fit_method = options.pop("fit_method", "loess")
    if options:
        unknown = ", ".join(sorted(options))
        raise TypeError(f"Unexpected DEqMS options: {unknown}")

    cond_a, cond_b = contrast
    stats_df = _collect_deqms_stats(
        log2_matrix,
        (samples_a, samples_b),
        contrast,
        peptide_counts,
    )
    if stats_df.empty:
        logger.warning("No proteins passed DE filtering for DEqMS")
        return pd.DataFrame()

    result = spectra_count_ebayes(_build_deqms_inputs(stats_df), fit_method=fit_method)
    return _finalize_deqms(stats_df, result, cond_a, cond_b)
