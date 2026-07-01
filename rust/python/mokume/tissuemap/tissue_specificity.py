"""AdaTiSS tissue specificity scoring (Jiang et al. 2020, PMID 32916130).

Algorithm summary
-----------------
For each protein *i*:

1. Collect all non-NaN log2 expression values *x*.
2. Estimate population parameters (μ, σ):
   - **Pure MAD mode** (default, ``use_pure_mad=True``): use median / MAD.
     Most robust across TMT and LFQ datasets.
   - **DPD mode** (``use_pure_mad=False``): fit a Gaussian population using
     Density Power Divergence (DPD, γ = 1).  Falls back to median / MAD
     for sparse proteins (< ``min_samples``).
3. Compute robust z-scores:  ``Z_ij = (X_ij − μ_i) / σ_i``
4. Tissue specificity score: ``TS_it = median(Z_ij  for j ∈ tissue t)``

Enrichment categories (per protein):
- **tissue-specific**: TS ≥ 4 in one tissue AND no other tissue has TS ∈ (2.5, 4)
- **tissue-enriched**: at least one TS ≥ 2.5 but not specific
- **house-keeping**: detected in all tissues AND all |TS| < 2
- **other**: everything else
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import norm as _norm
from sklearn.mixture import GaussianMixture

from mokume.tissuemap.config import TissueSpecificityConfig
from mokume.tissuemap.enrichment import (
    FLOOR_ENRICHED,
    FLOOR_HOUSEKEEPING,
    FLOOR_SPECIFIC,
    _classify_enrichment,
)

logger = logging.getLogger(__name__)


@dataclass
class PopulationParams:
    """Fitted population parameters for one protein."""

    mu: float
    sigma: float
    pi: float  # proportion of inliers


def _fit_population_dpd(
    x: np.ndarray,
    gamma: float = 1.0,
    max_iter: int = 100,
    tol: float = 1e-4,
    sigma_floor: float = 0.01,
) -> PopulationParams:
    """Fit Gaussian population via Density Power Divergence (DPD).

    Parameters
    ----------
    x : np.ndarray
        1-D array of non-NaN expression values.
    gamma : float
        DPD tuning parameter (γ = 1 recommended by paper).
    max_iter : int
        Maximum EM-like iterations.
    tol : float
        Convergence tolerance on μ and σ.
    sigma_floor : float
        Minimum allowed σ.

    Returns
    -------
    PopulationParams
    """
    mu = float(np.median(x))
    mad = float(np.median(np.abs(x - mu)))
    sigma = max(mad * 1.4826, sigma_floor)  # MAD → σ estimate

    for _ in range(max_iter):
        # Density-power weights
        z = (x - mu) / sigma
        w = np.exp(-gamma / 2.0 * z * z)
        w_sum = w.sum()
        if w_sum < 1e-12:
            break

        mu_new = float(np.dot(w, x) / w_sum)
        sigma_new = float(np.sqrt(np.dot(w, (x - mu_new) ** 2) / w_sum))
        sigma_new = max(sigma_new, sigma_floor)

        if abs(mu_new - mu) < tol and abs(sigma_new - sigma) < tol:
            mu, sigma = mu_new, sigma_new
            break
        mu, sigma = mu_new, sigma_new

    # Population proportion: fraction within ±2σ
    pi = float(np.mean(np.abs(x - mu) <= 2.0 * sigma))
    return PopulationParams(mu=mu, sigma=sigma, pi=pi)


def _fit_population_fallback(
    x: np.ndarray,
    sigma_floor: float = 0.01,
) -> PopulationParams:
    """Fallback for sparse proteins: median / MAD."""
    mu = float(np.median(x))
    mad = float(np.median(np.abs(x - mu)))
    sigma = max(mad * 1.4826, sigma_floor)
    pi = float(np.mean(np.abs(x - mu) <= 2.0 * sigma))
    return PopulationParams(mu=mu, sigma=sigma, pi=pi)


def _fit_pure_mad(
    x: np.ndarray,
    sigma_floor: float = 0.01,
) -> PopulationParams:
    """Pure MAD estimator — most robust across TMT and LFQ datasets.

    Uses median as μ and raw MAD (without 1.4826 scaling) as σ.
    This avoids DPD's tendency to collapse σ on TMT data while
    remaining consistent with LFQ data where σ is naturally large.
    """
    mu = float(np.median(x))
    mad = float(np.median(np.abs(x - mu)))
    sigma = max(mad, sigma_floor)
    pi = float(np.mean(np.abs(x - mu) <= 2.0 * sigma))
    return PopulationParams(mu=mu, sigma=sigma, pi=pi)


_SIGMA_FLOOR_ABSOLUTE_MIN = 0.01
_SIGMA_FLOOR_PERCENTILE = 5


def _resolve_sigma_floor(
    log2_matrix: np.ndarray,
    configured: float | None,
) -> float:
    """Determine the effective sigma floor.

    When *configured* is ``None`` (auto mode), computes the 5th percentile
    of per-protein MAD values.  This adapts the floor to the data's natural
    scale — higher for LFQ (large variance), lower for TMT (small variance).

    Parameters
    ----------
    log2_matrix : np.ndarray
        Samples (rows) x proteins (columns), log2 scale with NaN.
    configured : float | None
        User-specified floor, or ``None`` for auto-detection.

    Returns
    -------
    float
        Effective sigma floor (always >= ``_SIGMA_FLOOR_ABSOLUTE_MIN``).
    """
    if configured is not None:
        return max(configured, _SIGMA_FLOOR_ABSOLUTE_MIN)

    # Auto: 5th percentile of per-protein MAD
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        mads = np.nanmedian(
            np.abs(log2_matrix - np.nanmedian(log2_matrix, axis=0, keepdims=True)),
            axis=0,
        )
    mads = mads[~np.isnan(mads)]
    mads = mads[mads > 0]

    if len(mads) < 10:
        logger.warning(
            "Too few proteins for auto sigma_floor, using %.4f",
            _SIGMA_FLOOR_ABSOLUTE_MIN,
        )
        return _SIGMA_FLOOR_ABSOLUTE_MIN

    floor = float(np.percentile(mads, _SIGMA_FLOOR_PERCENTILE))
    floor = max(floor, _SIGMA_FLOOR_ABSOLUTE_MIN)
    logger.info(
        "Auto sigma_floor: %.4f (5th percentile of %d protein MADs)", floor, len(mads)
    )
    return floor


@dataclass
class GMMThresholds:
    """GMM-derived adaptive thresholds and fit parameters."""

    ts_enriched: float
    ts_specific: float
    bg_mean: float
    bg_std: float
    bg_weight: float
    sp_mean: float
    sp_std: float
    sp_weight: float


_GMM_MIN_PROTEINS = 50
_GMM_MAX_ITER = 200
_GMM_POSTERIOR_CUTOFF = 0.95


def _fallback_gmm(
    fallback_enriched: float,
    fallback_specific: float,
    means=None,
    stds=None,
    weights=None,
) -> GMMThresholds:
    """Return GMMThresholds with fixed fallback values."""
    return GMMThresholds(
        ts_enriched=fallback_enriched,
        ts_specific=fallback_specific,
        bg_mean=float(means[0]) if means is not None else 0.0,
        bg_std=float(stds[0]) if stds is not None else 1.0,
        bg_weight=float(weights[0]) if weights is not None else 0.5,
        sp_mean=float(means[1]) if means is not None else fallback_specific,
        sp_std=float(stds[1]) if stds is not None else 1.0,
        sp_weight=float(weights[1]) if weights is not None else 0.5,
    )


def _derive_gmm_thresholds(
    valid: np.ndarray,
    means,
    stds,
    weights,
) -> tuple[float, float]:
    """Compute enriched/specific thresholds from GMM components."""
    bg_idx = int(np.argmin(means))
    sp_idx = int(np.argmax(means))

    x_grid = np.linspace(float(valid.min()), float(valid.max()), 10000)
    p_bg = weights[bg_idx] * _norm.pdf(x_grid, means[bg_idx], stds[bg_idx])
    p_sp = weights[sp_idx] * _norm.pdf(x_grid, means[sp_idx], stds[sp_idx])

    # Enriched: intersection between means
    diff = p_bg - p_sp
    sign_changes = np.where(np.diff(np.sign(diff)))[0]
    intersections = x_grid[sign_changes]
    between = intersections[
        (intersections > means[bg_idx])
        & (intersections < means[sp_idx] + 3 * stds[sp_idx])
    ]
    ts_enriched = (
        float(between[0])
        if len(between) > 0
        else float(means[bg_idx] + 2 * stds[bg_idx])
    )

    # Specific: posterior P(specific|x) > cutoff
    posterior_sp = p_sp / (p_bg + p_sp + 1e-12)
    above_cutoff = x_grid[posterior_sp > _GMM_POSTERIOR_CUTOFF]
    ts_specific = float(above_cutoff[0]) if len(above_cutoff) > 0 else ts_enriched * 1.5
    ts_specific = max(ts_specific, ts_enriched + 0.1)

    return ts_enriched, ts_specific


def _auto_thresholds_gmm(
    max_ts: np.ndarray,
    *,
    fallback_enriched: float = FLOOR_ENRICHED,
    fallback_specific: float = FLOOR_SPECIFIC,
) -> GMMThresholds:
    """Fit a 2-component GMM and derive enriched/specific thresholds."""
    valid = max_ts[np.isfinite(max_ts)]
    if len(valid) < _GMM_MIN_PROTEINS:
        logger.warning(
            "Too few proteins (%d) for GMM, using fixed thresholds", len(valid)
        )
        return _fallback_gmm(fallback_enriched, fallback_specific)

    gmm = GaussianMixture(n_components=2, random_state=42, max_iter=_GMM_MAX_ITER)
    gmm.fit(valid.reshape(-1, 1))

    means = gmm.means_.flatten()
    stds = np.sqrt(gmm.covariances_.flatten())
    weights = gmm.weights_.flatten()

    bg_idx = int(np.argmin(means))
    sp_idx = int(np.argmax(means))

    if np.isclose(means[bg_idx], means[sp_idx], atol=1e-6):
        logger.warning("GMM components collapsed, using fixed thresholds")
        return _fallback_gmm(fallback_enriched, fallback_specific, means, stds, weights)

    ts_enriched, ts_specific = _derive_gmm_thresholds(valid, means, stds, weights)

    result = GMMThresholds(
        ts_enriched=ts_enriched,
        ts_specific=ts_specific,
        bg_mean=float(means[bg_idx]),
        bg_std=float(stds[bg_idx]),
        bg_weight=float(weights[bg_idx]),
        sp_mean=float(means[sp_idx]),
        sp_std=float(stds[sp_idx]),
        sp_weight=float(weights[sp_idx]),
    )
    logger.info(
        "GMM auto thresholds: enriched=%.3f, specific=%.3f "
        "(bg: mu=%.2f, sigma=%.2f, w=%.1f%%; sp: mu=%.2f, sigma=%.2f, w=%.1f%%)",
        result.ts_enriched,
        result.ts_specific,
        result.bg_mean,
        result.bg_std,
        result.bg_weight * 100,
        result.sp_mean,
        result.sp_std,
        result.sp_weight * 100,
    )
    return result


def _compute_ts_vectorized_mad(
    log2_matrix: np.ndarray,
    tissue_labels: np.ndarray,
    unique_tissues: list[str],
    sigma_floor: float,
) -> tuple[np.ndarray, list[PopulationParams]]:
    """Fully vectorized TS computation for the pure MAD path.

    Replaces the per-protein Python loop with bulk numpy operations.
    Loops only over tissues (~20-40), not proteins (~9000).
    """
    _n_samples, n_proteins = log2_matrix.shape
    n_tissues = len(unique_tissues)

    # Vectorized μ, MAD, σ for all proteins at once
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        mu = np.nanmedian(log2_matrix, axis=0)  # (n_proteins,)
        mad = np.nanmedian(
            np.abs(log2_matrix - mu[np.newaxis, :]),
            axis=0,
        )  # (n_proteins,)
    sigma = np.maximum(mad, sigma_floor)  # (n_proteins,)

    # Sparse proteins: < 3 non-NaN values → invalid
    n_valid = np.sum(~np.isnan(log2_matrix), axis=0)
    too_sparse = n_valid < 3
    mu[too_sparse] = np.nan
    sigma[too_sparse] = np.nan

    # Population proportion: fraction within ±2σ
    within_2s = np.abs(log2_matrix - mu[np.newaxis, :]) <= 2.0 * sigma[np.newaxis, :]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        pi = np.nanmean(within_2s.astype(np.float32), axis=0)
    pi[too_sparse] = 0.0

    # Full z-score matrix  (NaN propagates naturally)
    z_matrix = (log2_matrix - mu[np.newaxis, :]) / sigma[np.newaxis, :]

    # Per-tissue median z-scores — loop over tissues, not proteins
    ts_matrix = np.full((n_proteins, n_tissues), np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for t_idx, tissue in enumerate(unique_tissues):
            t_mask = tissue_labels == tissue
            z_tissue = z_matrix[t_mask, :]  # (n_tissue_samples, n_proteins)
            ts_matrix[:, t_idx] = np.nanmedian(z_tissue, axis=0)

    # Build PopulationParams list (lightweight Python objects)
    pop_params = [
        PopulationParams(
            mu=float(mu[i]) if not np.isnan(mu[i]) else np.nan,
            sigma=float(sigma[i]) if not np.isnan(sigma[i]) else np.nan,
            pi=float(pi[i]),
        )
        for i in range(n_proteins)
    ]
    return ts_matrix, pop_params


_DPD_GAMMA = 1.0
_DPD_MIN_SAMPLES = 50
_DPD_MAX_ITER = 100
_DPD_TOL = 1e-4


def _compute_ts_loop(
    log2_matrix: np.ndarray,
    tissue_labels: np.ndarray,
    unique_tissues: list[str],
    sigma_floor: float,
) -> tuple[np.ndarray, list[PopulationParams]]:
    """Per-protein loop path for DPD fitting (non-vectorizable)."""
    n_proteins = log2_matrix.shape[1]
    n_tissues = len(unique_tissues)
    ts_matrix = np.full((n_proteins, n_tissues), np.nan)
    pop_params: list[PopulationParams] = []

    for p_idx in range(n_proteins):
        col = log2_matrix[:, p_idx]
        x = col[~np.isnan(col)]

        if len(x) < 3:
            pop_params.append(PopulationParams(mu=np.nan, sigma=np.nan, pi=0.0))
            continue

        if len(x) >= _DPD_MIN_SAMPLES:
            params = _fit_population_dpd(
                x,
                gamma=_DPD_GAMMA,
                max_iter=_DPD_MAX_ITER,
                tol=_DPD_TOL,
                sigma_floor=sigma_floor,
            )
        else:
            params = _fit_population_fallback(x, sigma_floor=sigma_floor)
        pop_params.append(params)

        z_all = (col - params.mu) / params.sigma
        for t_idx, tissue in enumerate(unique_tissues):
            t_mask = tissue_labels == tissue
            z_tissue = z_all[t_mask]
            z_valid = z_tissue[~np.isnan(z_tissue)]
            if len(z_valid) >= 1:
                ts_matrix[p_idx, t_idx] = float(np.median(z_valid))

    return ts_matrix, pop_params


def _resolve_thresholds(
    config: TissueSpecificityConfig,
    max_ts_values: np.ndarray,
) -> tuple[float, float, float, GMMThresholds | None]:
    """Resolve enriched/specific/housekeeping thresholds from config or GMM."""
    need_auto = (
        config.ts_enriched_threshold is None
        or config.ts_specific_threshold is None
        or config.ts_housekeeping_threshold is None
    )

    gmm_result: GMMThresholds | None = None
    if need_auto:
        gmm_result = _auto_thresholds_gmm(
            max_ts_values,
            fallback_enriched=FLOOR_ENRICHED,
            fallback_specific=FLOOR_SPECIFIC,
        )

    if config.ts_enriched_threshold is not None:
        ts_enriched = config.ts_enriched_threshold
    else:
        ts_enriched = max(gmm_result.ts_enriched, FLOOR_ENRICHED)
        logger.info(
            "Auto enriched threshold: %.3f (GMM=%.3f, floor=%.1f)",
            ts_enriched,
            gmm_result.ts_enriched,
            FLOOR_ENRICHED,
        )

    if config.ts_specific_threshold is not None:
        ts_specific = config.ts_specific_threshold
    else:
        ts_specific = max(gmm_result.ts_specific, FLOOR_SPECIFIC)
        logger.info(
            "Auto specific threshold: %.3f (GMM=%.3f, floor=%.1f)",
            ts_specific,
            gmm_result.ts_specific,
            FLOOR_SPECIFIC,
        )

    if config.ts_housekeeping_threshold is not None:
        ts_housekeeping = config.ts_housekeeping_threshold
    else:
        ts_hk_raw = gmm_result.bg_mean + gmm_result.bg_std
        ts_housekeeping = max(ts_hk_raw, FLOOR_HOUSEKEEPING)
        logger.info(
            "Auto housekeeping threshold: %.3f (bg_mean=%.2f + bg_std=%.2f, floor=%.1f)",
            ts_housekeeping,
            gmm_result.bg_mean,
            gmm_result.bg_std,
            FLOOR_HOUSEKEEPING,
        )

    return ts_enriched, ts_specific, ts_housekeeping, gmm_result


def _build_ts_dataframe(
    ts_matrix: np.ndarray,
    protein_ids: np.ndarray,
    unique_tissues: list[str],
    pop_params: list[PopulationParams],
) -> pd.DataFrame:
    """Assemble TS DataFrame with population params and max tissue."""
    ts_df = pd.DataFrame(ts_matrix, index=protein_ids, columns=unique_tissues)
    ts_df["mu"] = [p.mu for p in pop_params]
    ts_df["sigma"] = [p.sigma for p in pop_params]
    ts_df["pi"] = [p.pi for p in pop_params]

    ts_only = ts_df[unique_tissues]
    ts_df["max_ts"] = ts_only.max(axis=1)
    finite_mask = np.isfinite(ts_df["max_ts"].values)
    ts_df["max_tissue"] = pd.NA
    if finite_mask.any():
        ts_df.loc[finite_mask, "max_tissue"] = ts_only.loc[finite_mask].idxmax(axis=1)
    return ts_df


def _finalize_ts_df(
    ts_df: pd.DataFrame,
    unique_tissues: list[str],
    tissue_labels: np.ndarray,
    log2_matrix: np.ndarray,
    ts_enriched: float,
    ts_specific: float,
    ts_housekeeping: float,
    gmm_result,
) -> None:
    """Add enrichment categories and threshold attrs to ts_df (in-place)."""
    ts_df["enrichment_category"] = _classify_enrichment(
        ts_df[unique_tissues].values,
        tissue_labels_all=tissue_labels,
        unique_tissues=unique_tissues,
        log2_matrix=log2_matrix,
        ts_enriched=ts_enriched,
        ts_specific=ts_specific,
        ts_housekeeping=ts_housekeeping,
    )
    ts_df.attrs["ts_enriched_threshold"] = ts_enriched
    ts_df.attrs["ts_specific_threshold"] = ts_specific
    ts_df.attrs["ts_housekeeping_threshold"] = ts_housekeeping
    if gmm_result is not None:
        ts_df.attrs["gmm"] = gmm_result


def compute_ts_scores(
    log2_matrix: np.ndarray,
    tissue_labels: np.ndarray,
    protein_ids: np.ndarray,
    config: TissueSpecificityConfig,
) -> pd.DataFrame:
    """Compute AdaTiSS tissue specificity scores.

    Returns a DataFrame (proteins × tissues) with TS scores plus columns
    ``mu``, ``sigma``, ``pi``, ``enrichment_category``, ``max_tissue``,
    ``max_ts``.
    """
    unique_tissues = sorted(np.unique(tissue_labels))
    n_proteins = log2_matrix.shape[1]
    sigma_floor = _resolve_sigma_floor(log2_matrix, config.sigma_floor)

    if config.use_pure_mad:
        ts_matrix, pop_params = _compute_ts_vectorized_mad(
            log2_matrix,
            tissue_labels,
            unique_tissues,
            sigma_floor,
        )
    else:
        ts_matrix, pop_params = _compute_ts_loop(
            log2_matrix,
            tissue_labels,
            unique_tissues,
            sigma_floor,
        )

    ts_df = _build_ts_dataframe(
        ts_matrix,
        protein_ids,
        unique_tissues,
        pop_params,
    )
    ts_enriched, ts_specific, ts_housekeeping, gmm_result = _resolve_thresholds(
        config,
        ts_df["max_ts"].dropna().values,
    )
    _finalize_ts_df(
        ts_df,
        unique_tissues,
        tissue_labels,
        log2_matrix,
        ts_enriched,
        ts_specific,
        ts_housekeeping,
        gmm_result,
    )

    cats = ts_df["enrichment_category"].value_counts()
    summary = ", ".join(
        f"{c}: {cats.get(c, 0)}"
        for c in ["tissue-specific", "tissue-enriched", "house-keeping", "other"]
    )
    logger.info("AdaTiSS TS scores (%d proteins): %s", n_proteins, summary)
    return ts_df
