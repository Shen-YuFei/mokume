"""
MLE (Maximum Likelihood Estimation) imputation for proteomics data.

Implements the EM algorithm for multivariate normal data with missing
values.  Iteratively estimates the mean vector and covariance matrix,
then imputes each missing entry as the conditional expectation given
the observed entries in the same row.

References
----------
- Little RJA, Rubin DB. Statistical Analysis with Missing Data.
  Wiley, 2002.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from mokume.core.logger import get_logger

logger = get_logger("mokume.imputation.mle")


def _em_step(data: np.ndarray, mu: np.ndarray, sigma: np.ndarray):
    """Single EM iteration: returns updated mu and sigma."""
    n, p = data.shape
    mu_new = np.zeros(p)
    sigma_new = np.zeros((p, p))
    x_filled = data.copy()

    for i in range(n):
        obs = ~np.isnan(data[i])
        mis = np.isnan(data[i])
        if not mis.any():
            x_filled[i] = data[i]
            continue
        if not obs.any():
            x_filled[i] = mu
            continue

        obs_idx = np.where(obs)[0]
        mis_idx = np.where(mis)[0]

        s_oo = sigma[np.ix_(obs_idx, obs_idx)]
        s_mo = sigma[np.ix_(mis_idx, obs_idx)]

        inv_oo = np.linalg.solve(
            s_oo + np.eye(len(obs_idx)) * 1e-8, np.eye(len(obs_idx))
        )
        cond_mean = mu[mis_idx] + s_mo @ inv_oo @ (data[i, obs_idx] - mu[obs_idx])
        x_filled[i, mis_idx] = cond_mean

    mu_new = np.nanmean(x_filled, axis=0)
    centered = x_filled - mu_new
    sigma_new = (centered.T @ centered) / n
    return mu_new, sigma_new, x_filled


def impute_mle(
    data: pd.DataFrame,
    max_iter: int = 50,
    tol: float = 1e-6,
) -> pd.DataFrame:
    """Impute missing values using EM-based MLE.

    Parameters
    ----------
    data : pd.DataFrame
        Protein intensity matrix (rows=proteins, columns=samples) with NaN.
    max_iter : int
        Maximum EM iterations (default 50).
    tol : float
        Convergence tolerance on the mean vector (default 1e-6).

    Returns
    -------
    pd.DataFrame
        Filled matrix with the same shape and index.
    """
    n_missing = int(data.isna().sum().sum())
    values = data.values.copy().astype(float)

    # Initialise mu and sigma from column-wise observed statistics
    mu = np.nanmean(values, axis=0)
    centered = np.where(np.isnan(values), 0, values - mu)
    n_obs = np.sum(~np.isnan(values), axis=0, keepdims=True)
    n_obs = np.maximum(n_obs, 1)
    sigma = (centered.T @ centered) / values.shape[0]
    sigma += np.eye(sigma.shape[0]) * 1e-6

    for iteration in range(max_iter):
        mu_old = mu.copy()
        mu, sigma, filled = _em_step(values, mu, sigma)
        sigma += np.eye(sigma.shape[0]) * 1e-6
        if np.linalg.norm(mu - mu_old) < tol:
            logger.info("MLE-EM converged at iteration %d", iteration + 1)
            break

    result = pd.DataFrame(filled, index=data.index, columns=data.columns)
    logger.info("MLE imputation: %d values imputed", n_missing)
    return result
