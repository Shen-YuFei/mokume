"""Sequential imputation (pure Python).

Reimplements the impSeq algorithm from rrcovNA using iterative regression
on variables ordered by missing-value frequency, without rpy2 or R
dependencies.

Reference
---------
Béguin C, Hulliger B. The BACON-EEM Algorithm for Multivariate Outlier
Detection in Incomplete Survey Data. *Survey Methodology*. 2008.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from mokume.core.logger import get_logger

logger = get_logger("mokume.imputation.impseq")


def _cond_normal_impute(
    row: np.ndarray,
    obs_idx: np.ndarray,
    miss_idx: np.ndarray,
    mu: np.ndarray,
    cov: np.ndarray,
) -> np.ndarray:
    """Impute missing values via conditional normal distribution."""
    sigma_mo = cov[np.ix_(miss_idx, obs_idx)]
    sigma_oo = cov[np.ix_(obs_idx, obs_idx)]
    try:
        sigma_oo_inv = np.linalg.solve(sigma_oo, np.eye(len(obs_idx)))
    except np.linalg.LinAlgError:
        return mu[miss_idx]
    imputed = mu[miss_idx] + sigma_mo @ sigma_oo_inv @ (row[obs_idx] - mu[obs_idx])
    return imputed


def _impseq_core(data: np.ndarray) -> np.ndarray:
    """Sequential imputation matching R rrcovNA::impSeq."""
    _, p = data.shape
    result = data.copy()
    observed = np.isfinite(data)

    complete_mask = observed.all(axis=1)
    n_complete = complete_mask.sum()

    if n_complete < max(2, p):
        col_means = np.nanmean(data, axis=0)
        for j in range(p):
            result[~observed[:, j], j] = col_means[j]
        return result

    complete_data = data[complete_mask]
    mu = np.mean(complete_data, axis=0)
    cov = np.cov(complete_data, rowvar=False, ddof=1)
    if cov.ndim == 0:
        cov = np.array([[cov]])

    incomplete_idx = np.where(~complete_mask)[0]
    n_miss_per_row = (~observed[incomplete_idx]).sum(axis=1)
    order = incomplete_idx[np.argsort(n_miss_per_row)]

    all_used = list(np.where(complete_mask)[0])

    for row_i in order:
        obs_j = np.where(observed[row_i])[0]
        miss_j = np.where(~observed[row_i])[0]
        if len(obs_j) == 0:
            result[row_i, miss_j] = mu[miss_j]
        else:
            result[row_i, miss_j] = _cond_normal_impute(
                result[row_i],
                obs_j,
                miss_j,
                mu,
                cov,
            )

        all_used.append(row_i)
        used_data = result[all_used]
        mu = np.mean(used_data, axis=0)
        if len(all_used) > 1:
            cov = np.cov(used_data, rowvar=False, ddof=1)
            if cov.ndim == 0:
                cov = np.array([[cov]])

    return result


def impute_impseq(
    data: pd.DataFrame,
) -> pd.DataFrame:
    """Impute missing values using sequential regression (pure Python).

    Parameters
    ----------
    data : pd.DataFrame
        Matrix with NaN for missing values.

    Returns
    -------
    pd.DataFrame
        Filled matrix with the same shape and index.
    """
    n_missing = int(data.isna().sum().sum())
    filled = _impseq_core(data.values.copy())
    result = pd.DataFrame(filled, index=data.index, columns=data.columns)
    logger.info("impSeq imputation: %d values imputed", n_missing)
    return result
