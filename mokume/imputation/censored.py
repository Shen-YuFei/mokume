"""
Censored-aware missing value imputation for proteomics data.

Implements multiple imputation strategies that account for the
non-random nature of missingness in mass spectrometry data:

- **MinProb**: Draws replacement values from the low tail of a
  fitted normal distribution (Perseus-style).
- **MinDet**: Replaces missing values with a fixed quantile of the
  per-sample observed values.

References
----------
- Lazar C, et al. Accounting for the Multiple Natures of Missing Values
  in Label-Free Quantitative Proteomics. J Proteome Res. 2016;15(4):1116-25.

"""

import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer

from mokume.core.logger import get_logger

logger = get_logger("mokume.imputation.censored")


def impute_minprob(
    data: pd.DataFrame,
    quantile: float = 0.01,
    shift: float = 1.6,
    scale: float = 0.3,
) -> pd.DataFrame:
    """Perform MinProb imputation using a low-tail normal draw per sample."""
    result = data.copy()
    n_imputed = 0

    for col in data.columns:
        missing = data[col].isna()
        if not missing.any():
            continue
        observed = data[col].dropna().values
        if len(observed) == 0:
            continue

        q_low = np.quantile(observed, quantile)
        sd = np.std(observed) * scale
        if sd < 1e-10:
            sd = 0.1
        mu = q_low - shift * sd

        imputed = np.random.normal(mu, sd, size=missing.sum())
        result.loc[missing, col] = imputed
        n_imputed += missing.sum()

    logger.info("MinProb imputation: %d values imputed", n_imputed)
    return result


def impute_mindet(
    data: pd.DataFrame,
    quantile: float = 0.01,
) -> pd.DataFrame:
    """Apply MinDet imputation using a per-sample observed quantile."""
    result = data.copy()
    n_imputed = 0

    for col in data.columns:
        missing = data[col].isna()
        if not missing.any():
            continue
        observed = data[col].dropna().values
        if len(observed) == 0:
            continue
        fill_val = np.quantile(observed, quantile)
        result.loc[missing, col] = fill_val
        n_imputed += missing.sum()

    logger.info("MinDet imputation: %d values imputed (q=%.3f)", n_imputed, quantile)
    return result


def impute_censored(
    data: pd.DataFrame,
    method: str = "minprob",
    **kwargs,
) -> pd.DataFrame:
    """Dispatch censored-aware imputation by method name."""
    method = method.lower()
    if method == "none":
        logger.info("No imputation applied (method='none')")
        return data.copy()
    if method == "minprob":
        return impute_minprob(data, **kwargs)
    if method == "mindet":
        return impute_mindet(data, **kwargs)
    if method == "knn":
        imputer = KNNImputer(
            n_neighbors=kwargs.get("n_neighbors", 5),
            weights=kwargs.get("weights", "uniform"),
        )
        result = pd.DataFrame(
            imputer.fit_transform(data),
            index=data.index,
            columns=data.columns,
        )
        n_imputed = data.isna().sum().sum()
        logger.info("KNN imputation: %d values imputed", n_imputed)
        return result
    raise ValueError(
        f"Unknown imputation method '{method}'. "
        f"Choose from: 'minprob', 'mindet', 'knn', 'none'"
    )
