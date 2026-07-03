"""
GMS (Gaussian Mixture Sampling) imputation for proteomics data.

Fits a two-component Gaussian mixture model per sample (column) to
distinguish the abundance distribution from the missing-not-at-random
(MNAR) component, then draws replacement values from the lower
component.

References
----------
- Lazar C, et al. Accounting for the Multiple Natures of Missing
  Values in Label-Free Quantitative Proteomics Data Sets.
  J Proteome Res. 2016;15(4):1116-25.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from mokume.core.logger import get_logger

logger = get_logger("mokume.imputation.gms")


def impute_gms(
    data: pd.DataFrame,
    n_components: int = 2,
    random_state: int = 42,
) -> pd.DataFrame:
    """Impute missing values using Gaussian Mixture Sampling.

    Parameters
    ----------
    data : pd.DataFrame
        Protein intensity matrix (rows=proteins, columns=samples) with NaN.
    n_components : int
        Number of Gaussian components (default 2).
    random_state : int
        Seed for reproducibility.

    Returns
    -------
    pd.DataFrame
        Filled matrix with the same shape and index.
    """
    from sklearn.mixture import GaussianMixture

    rng = np.random.default_rng(random_state)
    result = data.copy()
    n_imputed = 0

    for col in data.columns:
        missing = data[col].isna()
        if not missing.any():
            continue
        observed = data[col].dropna().values
        if len(observed) < n_components + 1:
            # Fallback: too few observations for GMM
            result.loc[missing, col] = np.nanmean(observed)
            n_imputed += int(missing.sum())
            continue

        gmm = GaussianMixture(
            n_components=n_components,
            random_state=random_state,
        )
        gmm.fit(observed.reshape(-1, 1))

        # Identify the lower-mean component
        lower_idx = int(np.argmin(gmm.means_.ravel()))
        mu = gmm.means_.ravel()[lower_idx]
        sigma = np.sqrt(gmm.covariances_.ravel()[lower_idx])
        if sigma < 1e-10:
            sigma = 0.1

        draws = rng.normal(mu, sigma, size=missing.sum())
        result.loc[missing, col] = draws
        n_imputed += int(missing.sum())

    logger.info("GMS imputation: %d values imputed", n_imputed)
    return result
