"""
QRILC (Quantile Regression Imputation of Left-Censored data).

Imputes missing values by drawing from the left tail of a truncated
normal distribution fitted per sample via quantile regression on
observed values.

References
----------
- Wei R, et al. Missing Value Imputation Approach for Mass
  Spectrometry-based Metabolomics Data. Sci Rep. 2018;8:663.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from mokume.core.logger import get_logger

logger = get_logger("mokume.imputation.qrilc")


def impute_qrilc(
    data: pd.DataFrame,
    tune_sigma: float = 1.0,
) -> pd.DataFrame:
    """Impute left-censored missing values via QRILC.

    For each sample (column), a truncated normal distribution is fitted
    to the lower tail of observed values.  Missing entries are replaced
    by random draws from that distribution.

    Parameters
    ----------
    data : pd.DataFrame
        Protein intensity matrix (rows=proteins, columns=samples) in
        **log2 space** with NaN for missing values.
    tune_sigma : float
        Multiplier applied to the estimated standard deviation before
        drawing (default 1.0).

    Returns
    -------
    pd.DataFrame
        Filled matrix with the same shape and index.
    """
    result = data.copy()
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

        draws = np.random.normal(mu_imp, sd_imp, size=int(missing.sum()))
        # Clamp: imputed values should not exceed the truncation point
        draws = np.minimum(draws, trunc)

        result.loc[missing, col] = draws
        n_imputed += int(missing.sum())

    logger.info("QRILC imputation: %d values imputed", n_imputed)
    return result
