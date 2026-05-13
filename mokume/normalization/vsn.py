"""
VSN (Variance Stabilizing Normalization) — pure Python.

Reimplements the vsn2 algorithm: fits per-sample affine + arsinh
transformation h(x) = arsinh(a + b*x) to stabilise variance across
the intensity range, without rpy2 or R dependencies.

References
----------
- Huber W, von Heydebreck A, Sueltmann H, Poustka A, Vingron M.
  Variance stabilization applied to microarray data calibration and
  to the quantification of differential expression. Bioinformatics.
  2002;18 Suppl 1:S96-104.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from mokume.core.logger import get_logger

logger = get_logger("mokume.normalization.vsn")


def _glog2(x: np.ndarray, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Generalized log2 matching R vsn: (arsinh(a + b*x) - log(2)) / log(2)."""
    z = a[np.newaxis, :] + b[np.newaxis, :] * x
    return (np.arcsinh(z) - np.log(2)) / np.log(2)


def _glog2_deriv(x: np.ndarray, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Derivative of glog2 w.r.t. x: b / (sqrt((a+bx)^2+1) * log(2))."""
    z = a[np.newaxis, :] + b[np.newaxis, :] * x
    return b[np.newaxis, :] / (np.sqrt(z**2 + 1) * np.log(2))


def _vsn_neg_log_lik(
    params: np.ndarray,
    data: np.ndarray,
) -> float:
    """Negative log-likelihood for the VSN model.

    L = -n*K/2*log(sigma2) - 1/(2*sigma2) * sum(h(x)-mu)^2 + sum(log|h'(x)|)
    We maximize L ⟺ minimize -L.
    """
    n_genes, n_samples = data.shape
    a = params[:n_samples]
    b = params[n_samples:]

    h = _glog2(data, a, b)
    h_prime = _glog2_deriv(data, a, b)

    obs = np.isfinite(data)
    h = np.where(obs, h, np.nan)
    h_prime = np.where(obs, h_prime, np.nan)

    mu = np.nanmean(h, axis=1, keepdims=True)
    resid = h - mu
    ss = np.nansum(resid**2)
    n_obs = np.sum(obs)
    sigma2 = max(ss / max(n_obs - n_genes, 1), 1e-12)

    log_jacobian = np.nansum(np.log(np.maximum(h_prime, 1e-300)))

    nll = 0.5 * n_obs * np.log(sigma2) + ss / (2 * sigma2) - log_jacobian
    return nll


def vsn_normalize(data: pd.DataFrame) -> pd.DataFrame:
    """Apply Variance Stabilizing Normalization (pure Python).

    Parameters
    ----------
    data : pd.DataFrame
        Wide matrix (rows=proteins, columns=samples) in **linear**
        intensity space.

    Returns
    -------
    pd.DataFrame
        VSN-transformed matrix (generalised-log scale) with same
        shape and index.
    """
    mat = data.values.astype(float)
    n_samples = mat.shape[1]

    a0 = np.zeros(n_samples)
    b0 = np.ones(n_samples)

    col_medians = np.nanmedian(mat, axis=0)
    ref_median = np.nanmedian(col_medians)
    b0 = np.where(col_medians > 0, ref_median / col_medians, 1.0)

    x0 = np.concatenate([a0, b0])

    bounds = [(None, None)] * n_samples + [(1e-10, None)] * n_samples

    opt_result = minimize(
        _vsn_neg_log_lik,
        x0,
        args=(mat,),
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": 1000, "ftol": 1e-10},
    )

    a_opt = opt_result.x[:n_samples]
    b_opt = opt_result.x[n_samples:]

    transformed = _glog2(mat, a_opt, b_opt)

    result_df = pd.DataFrame(transformed, index=data.index, columns=data.columns)
    logger.info("VSN normalization complete: %d proteins", len(result_df))
    return result_df
