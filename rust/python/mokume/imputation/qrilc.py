"""
QRILC-style left-tail imputation (simplified approximation).

Imputes left-censored missing values by drawing from a normal centered just
below the lowest observed value of each sample, clamping draws to that minimum.
This is a lightweight approximation of QRILC, NOT the canonical method: unlike
``imputeLCMD::impute.QRILC`` (Wei et al. 2018) it does not fit a truncated normal
by quantile regression and does not use the per-column missingness fraction
(``pNA``) to set the truncation quantile, so columns with very different missing
rates get the same fill distribution. Use it as a fast left-shifted MNAR fill,
not as a drop-in for imputeLCMD.

References
----------
- Wei R, et al. Missing Value Imputation Approach for Mass
  Spectrometry-based Metabolomics Data. Sci Rep. 2018;8:663
  (the canonical method this routine approximates).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from mokume.core.logger import get_logger

logger = get_logger("mokume.imputation.qrilc")


def impute_qrilc(
    data: pd.DataFrame,
    tune_sigma: float = 1.0,
    random_state: int = 42,
) -> pd.DataFrame:
    """Impute left-censored missing values with the simplified QRILC heuristic.

    For each sample (column), missing entries are replaced by draws from a
    normal centered at ``min(observed) - sd`` and clamped to ``min(observed)``.
    This is the simplified left-tail approximation described in the module
    docstring; it does not use the per-column missingness fraction and is not a
    drop-in for canonical imputeLCMD QRILC.

    Parameters
    ----------
    data : pd.DataFrame
        Protein intensity matrix (rows=proteins, columns=samples) in
        **log2 space** with NaN for missing values.
    tune_sigma : float
        Multiplier applied to the estimated standard deviation before
        drawing (default 1.0).
    random_state : int
        Seed for the draws. Without it repeated runs on the same matrix
        disagree, and a benchmark reads that spread as a difference between
        imputation methods rather than as noise within one.

    Returns
    -------
    pd.DataFrame
        Filled matrix with the same shape and index.
    """
    result = data.copy()
    rng = np.random.default_rng(random_state)
    n_imputed = 0

    for col in data.columns:
        missing = data[col].isna()
        if not missing.any():
            continue
        observed = data[col].dropna().values
        if len(observed) < 2:
            continue

        sd_obs = np.std(observed, ddof=1)
        if sd_obs < 1e-10:
            sd_obs = 0.1

        # Truncation point: lowest observed value
        trunc = np.min(observed)

        # Shift mean below truncation point
        mu_imp = trunc - sd_obs
        sd_imp = sd_obs * tune_sigma

        draws = rng.normal(mu_imp, sd_imp, size=int(missing.sum()))
        # Clamp: imputed values should not exceed the truncation point
        draws = np.minimum(draws, trunc)

        result.loc[missing, col] = draws
        n_imputed += int(missing.sum())

    logger.info("QRILC imputation: %d values imputed", n_imputed)
    return result
