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


def classify_missing(
    data: pd.DataFrame,
    quantile_threshold: float = 0.01,
) -> pd.DataFrame:
    """
    Classify each missing value as MNAR or MCAR.

    A missing value is classified as MNAR if the observed values for that
    protein are near or below the per-sample detection threshold.

    Parameters
    ----------
    data : pd.DataFrame
        Protein intensity matrix (proteins × samples), may contain NaN.
    quantile_threshold : float
        Quantile of per-sample observed intensities used as the censoring
        threshold. Default: 0.01 (1st percentile).

    Returns
    -------
    pd.DataFrame
        Boolean matrix of same shape as *data*, True where a value is
        classified as MNAR (left-censored), False for MCAR or observed.
    """
    # Per-sample censoring threshold
    thresholds = data.quantile(quantile_threshold, axis=0)

    is_mnar = pd.DataFrame(False, index=data.index, columns=data.columns)

    for col in data.columns:
        missing_mask = data[col].isna()
        if not missing_mask.any():
            continue
        # For each missing protein, check if its observed values in
        # other samples are near the threshold
        for protein in data.index[missing_mask]:
            obs_vals = data.loc[protein].dropna()
            if len(obs_vals) == 0:
                is_mnar.loc[protein, col] = True
            elif obs_vals.min() <= thresholds[col] * 2.0:
                is_mnar.loc[protein, col] = True

    n_missing = data.isna().sum().sum()
    n_mnar = is_mnar.sum().sum()
    logger.info(
        "Missing value classification: %d total, %d MNAR (%.1f%%), %d MCAR",
        n_missing, n_mnar,
        100 * n_mnar / max(n_missing, 1),
        n_missing - n_mnar,
    )
    return is_mnar


def impute_minprob(
    data: pd.DataFrame,
    quantile: float = 0.01,
    shift: float = 1.6,
    scale: float = 0.3,
) -> pd.DataFrame:
    """MinProb imputation: draw from N(q_low - shift*sd, scale*sd) per sample."""
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
    """
    MinDet (minimum detection) imputation.

    Replaces missing values with the *quantile*-th percentile of observed
    values in the corresponding sample.

    Parameters
    ----------
    data : pd.DataFrame
        Log2 protein intensity matrix (proteins × samples).
    quantile : float
        Quantile to use as the imputation value (default 0.01).

    Returns
    -------
    pd.DataFrame
        Imputed matrix.
    """
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
    """
    Unified interface for censored-aware imputation.

    Parameters
    ----------
    data : pd.DataFrame
        Log2 protein intensity matrix (proteins × samples).
    method : str
        Imputation method: ``"minprob"``, ``"mindet"``,
        ``"knn"``, or ``"none"``.
    **kwargs
        Method-specific keyword arguments.

    Returns
    -------
    pd.DataFrame
        Imputed matrix.
    """
    method = method.lower()
    if method == "none":
        logger.info("No imputation applied (method='none')")
        return data.copy()
    elif method == "minprob":
        return impute_minprob(data, **kwargs)
    elif method == "mindet":
        return impute_mindet(data, **kwargs)
    elif method == "knn":
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
    else:
        raise ValueError(
            f"Unknown imputation method '{method}'. "
            f"Choose from: 'minprob', 'mindet', 'knn', 'none'"
        )
