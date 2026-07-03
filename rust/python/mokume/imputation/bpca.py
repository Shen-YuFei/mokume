"""Bayesian PCA imputation (pure Python).

Reimplements the BPCA EM algorithm for imputing missing values using a
low-rank probabilistic model, without rpy2 or R dependencies.

Reference
---------
Oba S, Sato M, Takemasa I, et al. A Bayesian missing value estimation
method for gene expression profile data. *Bioinformatics*. 2003;19(16):2088-96.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from mokume.core.logger import get_logger

logger = get_logger("mokume.imputation.bpca")


def _bpca_initmodel(y: np.ndarray, k: int) -> dict:
    """Initialise BPCA model — matches R pcaMethods:::BPCA_initmodel."""
    n, p = y.shape
    nans = np.isnan(y)
    yest = y.copy()
    yest[nans] = 0.0

    row_miss_count = nans.sum(axis=1)
    row_nomiss = np.where(row_miss_count == 0)[0]
    row_miss = np.where(row_miss_count > 0)[0]

    covy = np.cov(yest, rowvar=False, ddof=1)
    if covy.ndim == 0:
        covy = np.array([[float(covy)]])

    U, s, _ = np.linalg.svd(covy)
    U_k = U[:, :k]
    S_k = np.diag(s[:k])
    PA = U_k @ np.sqrt(S_k)

    mean_vec = np.zeros(p)
    for j in range(p):
        obs_j = ~nans[:, j]
        if obs_j.any():
            mean_vec[j] = np.mean(y[obs_j, j])

    denom = np.sum(np.diag(covy)) - np.sum(s[:k])
    tau = 1.0 / max(min(abs(denom), 1e10), 1e-10)
    if denom < 0:
        tau = 1e-10

    galpha0 = 1e-10
    balpha0 = 1.0
    alpha = (2 * galpha0 + p) / (tau * np.diag(PA.T @ PA) + 2 * galpha0 / balpha0)

    return {
        "n": n,
        "p": p,
        "k": k,
        "yest": yest,
        "nans": nans,
        "row_nomiss": row_nomiss,
        "row_miss": row_miss,
        "mean": mean_vec,
        "PA": PA,
        "tau": tau,
        "alpha": alpha,
        "SigW": np.eye(k),
        "galpha0": galpha0,
        "balpha0": balpha0,
        "gmu0": 0.001,
        "btau0": 1.0,
        "gtau0": 1e-10,
        "scores": np.zeros((n, k)),
    }


def _bpca_dostep(M: dict, y: np.ndarray) -> dict:
    """One EM step — matches R pcaMethods:::BPCA_dostep."""
    n, p, k = M["n"], M["p"], M["k"]
    PA = M["PA"]
    tau = M["tau"]
    nans = M["nans"]
    mean_vec = M["mean"]

    scores = np.full((n, k), np.nan)
    Rx = np.eye(k) + tau * (PA.T @ PA) + M["SigW"]
    Rxinv = np.linalg.inv(Rx)

    idx = M["row_nomiss"]
    if len(idx) == 0:
        trS = 0.0
        T = np.zeros((p, k))
    else:
        dy = y[idx] - mean_vec[np.newaxis, :]
        x = tau * Rxinv @ PA.T @ dy.T
        T = dy.T @ x.T
        trS = float(np.sum(dy * dy))
        for ii, row_idx in enumerate(idx):
            scores[row_idx] = x[:, ii]

    for row_idx in M["row_miss"]:
        obs_mask = ~nans[row_idx]
        miss_mask = nans[row_idx]

        dyo = y[row_idx, obs_mask] - mean_vec[obs_mask]
        Wm = PA[miss_mask]
        Wo = PA[obs_mask]

        Rx_local_inv = np.linalg.inv(Rx - tau * (Wm.T @ Wm))
        ex = tau * Wo.T @ dyo
        x = Rx_local_inv @ ex
        dym = Wm @ x

        dy_full = np.zeros(p)
        dy_full[obs_mask] = dyo
        dy_full[miss_mask] = dym

        M["yest"][row_idx] = dy_full + mean_vec

        T += np.outer(dy_full, x)
        T[miss_mask] += Wm @ Rx_local_inv

        n_miss = int(miss_mask.sum())
        trS += (
            float(dy_full @ dy_full)
            + n_miss / tau
            + float(np.trace(Wm @ Rx_local_inv @ Wm.T))
        )
        scores[row_idx] = x

    M["scores"] = scores
    T = T / n
    trS = trS / n

    Rxinv = np.linalg.inv(Rx)
    Dw = Rxinv + tau * T.T @ PA @ Rxinv + np.diag(M["alpha"]) / n
    Dwinv = np.linalg.inv(Dw)

    M["PA"] = T @ Dwinv

    tau_num = p + 2 * M["gtau0"] / n
    tau_denom = (
        trS
        - float(np.trace(T.T @ M["PA"]))
        + float((mean_vec @ mean_vec) * M["gmu0"] + 2 * M["gtau0"] / M["btau0"]) / n
    )
    M["tau"] = float(tau_num / max(abs(tau_denom), 1e-10))
    M["tau"] = max(min(M["tau"], 1e10), 1e-10)

    M["SigW"] = Dwinv * (p / n)
    M["alpha"] = (2 * M["galpha0"] + p) / (
        M["tau"] * np.diag(M["PA"].T @ M["PA"])
        + np.diag(M["SigW"])
        + 2 * M["galpha0"] / M["balpha0"]
    )

    return M


def _bpca_em(
    data: np.ndarray,
    n_pcs: int = 2,
    max_steps: int = 100,
    threshold: float = 1e-4,
) -> np.ndarray:
    """BPCA imputation — matches R pcaMethods::bpca + pca() wrapper."""
    n, p = data.shape
    k = min(n_pcs, min(n, p) - 1)
    if k < 1:
        result = data.copy()
        mu = np.nanmean(data, axis=0)
        nans = np.isnan(data)
        for j in range(p):
            result[nans[:, j], j] = mu[j]
        return result

    center = np.nanmean(data, axis=0)
    y_centered = data.copy() - center
    y_centered[np.isnan(data)] = np.nan

    M = _bpca_initmodel(y_centered, k)

    tau_old = 1000.0
    for step in range(1, max_steps + 1):
        M = _bpca_dostep(M, y_centered)
        if step % 10 == 0:
            dtau = abs(np.log10(max(M["tau"], 1e-300)) - np.log10(max(tau_old, 1e-300)))
            if dtau < threshold:
                break
            tau_old = M["tau"]

    recon = M["scores"] @ M["PA"].T + center
    result = data.copy()
    nans = np.isnan(data)
    result[nans] = recon[nans]
    return result


def impute_bpca(
    data: pd.DataFrame,
    n_pcs: int = 2,
) -> pd.DataFrame:
    """Impute missing values using Bayesian PCA (pure Python).

    Parameters
    ----------
    data : pd.DataFrame
        Protein intensity matrix (rows=proteins, columns=samples) with NaN.
    n_pcs : int
        Number of principal components to use (default 2).

    Returns
    -------
    pd.DataFrame
        Filled matrix with the same shape and index.
    """
    n_missing = int(data.isna().sum().sum())
    filled = _bpca_em(data.values.copy(), n_pcs=n_pcs)
    result = pd.DataFrame(filled, index=data.index, columns=data.columns)
    logger.info("BPCA imputation: %d values imputed", n_missing)
    return result
