"""Robust sequential imputation (pure Python).

Reimplements the impSeqRob algorithm from rrcovNA using iterative
robust regression with MCD covariance estimation, without rpy2 or R
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

logger = get_logger("mokume.imputation.impseqrob")

_MAX_EXTRADIR_DIRECTIONS = 500
_RANDOM_DIRECTION_SEED = 0

# The covariance core is O(p^3) in the number of columns (variables): it
# inverts a p x p (inverse-)covariance per imputed row. Beyond a few thousand
# columns this hangs for hours on wide cohorts, so refuse by default and point
# the caller at column subsetting or a distribution-wise imputer. Callers who
# really mean it can pass ``max_variables=None``.
_MAX_VARIABLES_GUARD = 2000


def _h_alpha_n(alpha: float, n: int, p: int) -> int:
    """Compute h — matches robustbase:::h.alpha.n."""
    n2 = (n + p + 1) // 2
    return int(np.floor(2 * n2 - n + 2 * (n - n2) * alpha))


def _extradir(
    data: np.ndarray,
    max_directions: int = _MAX_EXTRADIR_DIRECTIONS,
) -> np.ndarray:
    """Difference directions for Stahel-Donoho outlyingness.

    rrcovNA uses all pairwise directions only for small complete sets; for
    larger sets it samples at most 500 row-pair directions. Capping this keeps
    the robust covariance step aligned with rrcovNA and avoids O(n^2) work.
    """
    n = data.shape[0]
    n_pairs = n * (n - 1) // 2
    if n_pairs == 0:
        return np.empty((0, data.shape[1]), dtype=data.dtype)
    if n_pairs > max_directions:
        rng = np.random.default_rng(_RANDOM_DIRECTION_SEED)
        first = rng.integers(0, n, size=max_directions)
        second = rng.integers(0, n - 1, size=max_directions)
        second = second + (second >= first)
        return data[first] - data[second]

    idx = np.array([(i, j) for i in range(n) for j in range(i + 1, n)], dtype=int)
    return data[idx[:, 0]] - data[idx[:, 1]]


def _covSD(x: np.ndarray, h: int) -> dict:
    """Stahel-Donoho robust covariance — matches rrcovNA:::.covSD."""
    n, p = x.shape
    B = _extradir(x)
    norms = np.linalg.norm(B, axis=1)
    keep = norms > 1e-12
    B = B[keep]
    norms = norms[keep]
    A = B / norms[:, np.newaxis]

    Y = x @ A.T
    medY = np.median(Y, axis=0)
    madY = np.median(np.abs(Y - medY), axis=0) * 1.4826
    madY = np.where(madY < 1e-12, 1e-12, madY)
    Z = np.abs(Y - medY) / madY
    d = np.max(Z, axis=1)

    ds_ix = np.argsort(d, kind="stable")
    ds_vals = d[ds_ix]
    subset = ds_ix[:h]

    covx = np.cov(x[subset], rowvar=False, ddof=1)
    if covx.ndim == 0:
        covx = np.array([[float(covx)]])
    mx = np.mean(x[subset], axis=0)

    return {
        "cov": covx,
        "mean": mx,
        "d": d,
        "A": A,
        "medY": medY,
        "madY": madY,
        "ds_ix": ds_ix,
        "ds_vals": ds_vals,
    }


def _update_invcov(A: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Sherman-Morrison rank-one inverse update — matches rrcovNA:::.update_invcov."""
    Af = A @ f
    denom = 1.0 + f @ Af
    if abs(denom) > 1e-12:
        return A - np.outer(Af, Af) / denom
    return A


def _impnorm_preimpute(x: np.ndarray, n_needed: int) -> np.ndarray:
    """Pre-impute rows using conditional MVN — matches R rrcovNA:::.impNorm.

    Estimates mu and Sigma from pairwise complete observations, then imputes
    the ``n_needed`` rows with fewest missing values using the conditional
    normal distribution.
    """
    n, p = x.shape
    isnan = np.isnan(x)
    n_miss_per_row = isnan.sum(axis=1)

    # Column means
    mu = np.nanmean(x, axis=0)

    # Pairwise complete-case covariance
    sigma = np.zeros((p, p))
    for j in range(p):
        for k in range(j, p):
            both_obs = (~isnan[:, j]) & (~isnan[:, k])
            if both_obs.sum() > 1:
                xj = x[both_obs, j]
                xk = x[both_obs, k]
                sigma[j, k] = sigma[k, j] = np.cov(xj, xk, ddof=1)[0, 1]
            elif j == k:
                col_obs = x[~isnan[:, j], j]
                sigma[j, j] = (
                    float(np.var(col_obs, ddof=1)) if len(col_obs) > 1 else 1.0
                )

    # Regularise to ensure positive-definiteness
    eigvals = np.linalg.eigvalsh(sigma)
    if eigvals.min() < 1e-10:
        sigma += (abs(eigvals.min()) + 1e-6) * np.eye(p)

    # Select rows with fewest missing values (but >0) — impute these first
    candidate_idx = np.where(n_miss_per_row > 0)[0]
    order = candidate_idx[np.argsort(n_miss_per_row[candidate_idx])]
    impute_rows = order[:n_needed]

    for i in impute_rows:
        obs_mask = ~isnan[i]
        miss_mask = isnan[i]
        obs_idx = np.where(obs_mask)[0]
        miss_idx = np.where(miss_mask)[0]

        if len(obs_idx) == 0:
            x[i, miss_mask] = mu[miss_mask]
            continue

        sigma_mo = sigma[np.ix_(miss_idx, obs_idx)]
        sigma_oo = sigma[np.ix_(obs_idx, obs_idx)]

        try:
            sigma_oo_inv = np.linalg.inv(sigma_oo)
        except np.linalg.LinAlgError:
            sigma_oo_inv = np.linalg.pinv(sigma_oo)

        x[i, miss_mask] = mu[miss_mask] + sigma_mo @ sigma_oo_inv @ (
            x[i, obs_mask] - mu[obs_mask]
        )

    return x


def _impseqrob_core(
    data: np.ndarray,
    alpha: float = 0.9,
) -> np.ndarray:
    """Robust sequential imputation — matches rrcovNA::impSeqRob."""
    n, p = data.shape
    observed = np.isfinite(data)

    n_miss_per_row = (~observed).sum(axis=1)
    sort_order = np.argsort(n_miss_per_row, kind="stable")
    unsort = np.argsort(sort_order, kind="stable")

    x = data[sort_order].copy()
    isnan_x = np.isnan(x)
    risnan = n_miss_per_row[sort_order]

    complobs = np.where(risnan == 0)[0]
    ncomplobs = len(complobs)
    misobs = np.where(risnan > 0)[0]
    nmisobs = len(misobs)

    if nmisobs == 0:
        return data.copy()

    mincomplete = int(np.ceil(p / alpha))
    if ncomplobs < mincomplete:
        # Pre-impute using conditional MVN — matches R rrcovNA:::.impNorm
        x = _impnorm_preimpute(x, mincomplete - ncomplobs)
        isnan_x = np.isnan(x)
        risnan = isnan_x.sum(axis=1)
        complobs = np.where(risnan == 0)[0]
        ncomplobs = len(complobs)
        misobs = np.where(risnan > 0)[0]
        nmisobs = len(misobs)
        if nmisobs == 0:
            return x[unsort]

    h = _h_alpha_n(alpha, ncomplobs, p)
    covm = _covSD(x[complobs], h)
    covx = covm["cov"]
    mx = covm["mean"]
    A = covm["A"]
    medY = covm["medY"]
    madY = covm["madY"]
    d_list = list(covm["d"])
    ds_vals = np.sort(d_list)

    flag = np.ones(n, dtype=int)
    flag[complobs[covm["ds_ix"][:h]]] = 0

    outl = np.zeros(n)
    outl[complobs] = covm["d"]

    nupdate = h
    icovx = None
    mxo = mx.copy()

    for inn in range(nmisobs):
        i = misobs[inn]
        mvar = isnan_x[i]
        obs_mask = ~mvar

        if inn == 0:
            icovx = np.linalg.inv(covx)
        elif flag[misobs[inn - 1]] == 0:
            prev_i = misobs[inn - 1]
            f = (x[prev_i] - mxo) / np.sqrt(max(nupdate - 1, 1))
            scaled = (nupdate - 1) / max(nupdate - 2, 1) * icovx
            icovx = _update_invcov(scaled, f)

        xo = x[i, obs_mask]
        miss_idx = np.where(mvar)[0]
        obs_idx = np.where(obs_mask)[0]

        icov_mm = icovx[np.ix_(miss_idx, miss_idx)]
        icov_mo = icovx[np.ix_(miss_idx, obs_idx)]
        try:
            icov_mm_inv = np.linalg.inv(icov_mm)
        except np.linalg.LinAlgError:
            icov_mm_inv = np.linalg.pinv(icov_mm)

        x[i, mvar] = mx[mvar] - icov_mm_inv @ icov_mo @ (xo - mx[obs_mask])

        Y_row = x[i] @ A.T
        Zx = np.abs(Y_row - medY) / madY
        dx = float(np.max(Zx))

        if dx <= ds_vals[min(h, len(ds_vals) - 1)]:
            flag[i] = 0
            nupdate += 1
            h += 1
            mxo = mx.copy()
            mx = ((nupdate - 1) * mxo + x[i]) / nupdate
            covx = (nupdate - 2) / (nupdate - 1) * covx + np.outer(
                x[i] - mxo, x[i] - mxo
            ) / (nupdate - 1)
        else:
            mxo = mx.copy()

        outl[i] = dx
        d_list.append(dx)
        ds_vals = np.sort(d_list)

    result = x[unsort]
    return result


def impute_impseqrob(
    data: pd.DataFrame,
    alpha: float = 0.9,
    max_variables: int | None = _MAX_VARIABLES_GUARD,
) -> pd.DataFrame:
    """Impute missing values using robust sequential regression (pure Python).

    Parameters
    ----------
    data : pd.DataFrame
        Matrix (rows = features, columns = variables) with NaN for missing
        values.
    alpha : float
        MCD support fraction (default 0.5).
    max_variables : int or None
        Guard against the O(p^3) blow-up: raise if the matrix has more than
        this many columns. Set to ``None`` to disable and run regardless.

    Returns
    -------
    pd.DataFrame
        Filled matrix with the same shape and index.
    """
    n_vars = data.shape[1]
    if max_variables is not None and n_vars > max_variables:
        raise ValueError(
            f"impSeqRob is covariance-based and scales as O(p^3) in the number "
            f"of columns (here p={n_vars} > max_variables={max_variables}); on a "
            f"wide matrix it can run for hours. Subset the columns to the samples "
            f"being compared (e.g. the two contrast groups), or use a "
            f"distribution-wise imputer such as qrilc/mindet/minprob. Pass "
            f"max_variables=None to override and run anyway."
        )
    n_missing = int(data.isna().sum().sum())
    filled = _impseqrob_core(data.values.copy(), alpha=alpha)
    result = pd.DataFrame(filled, index=data.index, columns=data.columns)
    logger.info("impSeqRob imputation: %d values imputed", n_missing)
    return result
