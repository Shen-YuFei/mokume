"""Probabilistic dropout model for differential expression analysis.

Implements a simplified version of the proDA approach (Ahlmann-Eltze & Anders,
2019) which models missing values as probabilistic dropouts rather than
requiring imputation. Each sample gets a sigmoidal dropout curve, and the
likelihood integrates over both observed and missing values.

Reference
---------
Ahlmann-Eltze C, Anders S. proDA: Probabilistic Dropout Analysis for
Identifying Differentially Abundant Proteins in Label-Free Quantitative
Mass Spectrometry. *bioRxiv*. 2019. doi:10.1101/661496

"""

from typing import NamedTuple

import numpy as np
import pandas as pd
from scipy import stats, optimize
from scipy.special import expit  # logistic sigmoid
from statsmodels.stats.multitest import multipletests

from mokume.core.logger import get_logger

logger = get_logger("mokume.analysis.proda")


class DropoutParams(NamedTuple):

    """Per-sample dropout curve parameters for one condition."""

    rho: np.ndarray
    zeta: np.ndarray


def _fit_dropout_curve(intensities: np.ndarray) -> tuple[float, float]:
    """Fit a sigmoidal dropout curve P(miss|mu)=1-sigmoid((mu-rho)/zeta)."""
    observed = intensities[np.isfinite(intensities)]
    n_total = len(intensities)
    n_obs = len(observed)

    if n_obs < 5 or n_obs == n_total:
        return np.nanmedian(observed) - 2.0, 1.0

    obs_min = np.min(observed)

    # Initial guess: rho near where dropout starts
    rho_init = obs_min - 1.0
    zeta_init = 1.0

    def neg_log_lik(params):
        rho, zeta = params
        if zeta <= 0.01:
            return 1e10
        # P(observed) = sigmoid((x - rho) / zeta)
        p_obs = expit((observed - rho) / zeta)
        p_obs = np.clip(p_obs, 1e-10, 1.0 - 1e-10)
        # Log-likelihood for observed values
        ll_obs = np.sum(np.log(p_obs))
        # For missing: P(missing) = 1 - sigmoid((mu - rho) / zeta)
        # Approximate with the expected dropout rate
        n_miss = n_total - n_obs
        if n_miss > 0:
            mean_p_miss = 1.0 - np.mean(p_obs)
            mean_p_miss = np.clip(mean_p_miss, 1e-10, 1.0 - 1e-10)
            ll_miss = n_miss * np.log(mean_p_miss)
        else:
            ll_miss = 0.0
        return -(ll_obs + ll_miss)

    try:
        result = optimize.minimize(
            neg_log_lik,
            x0=[rho_init, zeta_init],
            method="Nelder-Mead",
            options={"maxiter": 500, "xatol": 1e-4, "fatol": 1e-4},
        )
        rho, zeta = result.x
        zeta = max(zeta, 0.1)
    except (RuntimeError, ValueError) as exc:
        logger.debug("Dropout curve fit failed: %s", exc)
        rho = rho_init
        zeta = zeta_init

    return rho, zeta


def _mean_or_rho(
    obs: np.ndarray, dropout_rho: np.ndarray, fallback_shift: float = -3.0,
) -> float:
    """Return mean of observed values, or dropout rho if all missing."""
    if len(obs) > 0:
        return float(np.mean(obs))
    if len(dropout_rho) > 0:
        return float(np.mean(dropout_rho))
    return fallback_shift


def _dropout_se(
    obs: np.ndarray, n_total: int,
) -> float:
    """Compute dropout-inflated standard error for one condition."""
    n_obs = len(obs)
    var = np.var(obs, ddof=1)
    miss_frac = 1.0 - n_obs / n_total if n_total > 0 else 0.0
    return var * (1.0 + miss_frac) / n_obs


def _protein_likelihood_test(
    vals_a: np.ndarray,
    vals_b: np.ndarray,
    dropout_a: DropoutParams,
    dropout_b: DropoutParams,
) -> tuple[float, float, float]:
    """Wald test for one protein with dropout variance inflation."""
    obs_a = vals_a[np.isfinite(vals_a)]
    obs_b = vals_b[np.isfinite(vals_b)]

    if len(obs_a) == 0 and len(obs_b) == 0:
        return np.nan, np.nan, np.nan

    if len(obs_a) == 0:
        return _mean_or_rho(obs_a, dropout_a.rho) - np.mean(obs_b), np.nan, np.nan
    if len(obs_b) == 0:
        return np.mean(obs_a) - _mean_or_rho(obs_b, dropout_b.rho), np.nan, np.nan

    log2fc = np.mean(obs_a) - np.mean(obs_b)
    if len(obs_a) < 2 or len(obs_b) < 2:
        return log2fc, np.nan, np.nan

    se = np.sqrt(_dropout_se(obs_a, len(vals_a)) + _dropout_se(obs_b, len(vals_b)))
    if se < 1e-15:
        return log2fc, np.nan, np.nan

    t_stat = log2fc / se
    df = len(obs_a) + len(obs_b) - 2
    pvalue = 2.0 * stats.t.sf(np.abs(t_stat), df=max(df, 1))
    return log2fc, t_stat, pvalue


def _fit_dropouts(
    log2_matrix: pd.DataFrame,
    samples_a: list[str],
    samples_b: list[str],
) -> tuple[DropoutParams, DropoutParams]:
    """Fit per-sample dropout curves and pack into per-condition params."""
    def _pack_dropouts(sample_names: list[str], dropout_params: dict[str, tuple[float, float]]) -> DropoutParams:
        return DropoutParams(
            rho=np.array([dropout_params[s][0] for s in sample_names]),
            zeta=np.array([dropout_params[s][1] for s in sample_names]),
        )

    all_samples = samples_a + samples_b
    logger.info("Fitting dropout curves for %d samples", len(all_samples))
    dp = {s: _fit_dropout_curve(log2_matrix[s].values) for s in all_samples}
    return _pack_dropouts(samples_a, dp), _pack_dropouts(samples_b, dp)


def _test_all_proteins(
    log2_matrix: pd.DataFrame,
    sample_groups: tuple[list[str], list[str]],
    contrast: tuple[str, str],
    dropouts: tuple[DropoutParams, DropoutParams],
) -> list[dict]:
    """Apply dropout-aware Wald test to every protein."""
    samples_a, samples_b = sample_groups
    cond_a, cond_b = contrast
    dropout_a, dropout_b = dropouts
    results = []
    for protein in log2_matrix.index:
        va = log2_matrix.loc[protein, samples_a].values.astype(float)
        vb = log2_matrix.loc[protein, samples_b].values.astype(float)
        log2fc, t_stat, pval = _protein_likelihood_test(va, vb, dropout_a, dropout_b)
        if np.isnan(pval):
            continue
        oa, ob = va[np.isfinite(va)], vb[np.isfinite(vb)]
        results.append({
            "ProteinName": protein, "log2FC": log2fc, "pvalue": pval,
            "t_stat": t_stat,
            f"mean_{cond_a}": np.mean(oa) if len(oa) else np.nan,
            f"mean_{cond_b}": np.mean(ob) if len(ob) else np.nan,
            "n_a": len(oa), "n_b": len(ob),
            "n_total_a": len(samples_a), "n_total_b": len(samples_b),
        })
    return results


def run_proda(
    log2_matrix: pd.DataFrame,
    samples_a: list[str],
    samples_b: list[str],
    cond_a: str,
    cond_b: str,
) -> pd.DataFrame:
    """Run proDA-style probabilistic dropout DE analysis."""
    dropout_a, dropout_b = _fit_dropouts(log2_matrix, samples_a, samples_b)
    rows = _test_all_proteins(
        log2_matrix,
        (samples_a, samples_b),
        (cond_a, cond_b),
        (dropout_a, dropout_b),
    )
    if not rows:
        logger.warning("No proteins passed DE filtering for proDA")
        return pd.DataFrame()

    de_df = pd.DataFrame(rows)
    _, adj_p, _, _ = multipletests(de_df["pvalue"].values, method="fdr_bh")
    de_df["adj_pvalue"] = adj_p
    return de_df.sort_values("adj_pvalue").reset_index(drop=True)
