"""
Neighbor-averaging (nbavg) imputation for proteomics data.

For each missing value, the average of the K nearest non-missing
neighbors *within the same row* (protein) is used.  If no neighbor
is available, the row mean is used as fallback.

This is a simple, fast heuristic suited to matrices where row-level
(protein-level) correlation across samples is expected.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from mokume.core.logger import get_logger

logger = get_logger("mokume.imputation.nbavg")


def impute_nbavg(
    data: pd.DataFrame,
    k: int = 5,
) -> pd.DataFrame:
    """Impute missing values by averaging K nearest column neighbors.

    Parameters
    ----------
    data : pd.DataFrame
        Protein intensity matrix (rows=proteins, columns=samples) with NaN.
    k : int
        Number of nearest non-missing neighbors to average (default 5).

    Returns
    -------
    pd.DataFrame
        Filled matrix with the same shape and index.
    """
    result = data.copy()
    values = data.values
    n_imputed = 0

    for i in range(values.shape[0]):
        row = values[i]
        missing_idx = np.where(np.isnan(row))[0]
        if len(missing_idx) == 0:
            continue
        observed_idx = np.where(~np.isnan(row))[0]
        if len(observed_idx) == 0:
            continue

        row_mean = np.nanmean(row)
        for j in missing_idx:
            if len(observed_idx) == 0:
                result.iat[i, j] = row_mean
                continue
            # Distance = absolute column-index difference
            dists = np.abs(observed_idx - j)
            nearest = observed_idx[np.argsort(dists)[:k]]
            result.iat[i, j] = np.mean(row[nearest])
            n_imputed += 1

    logger.info("nbavg imputation: %d values imputed (k=%d)", n_imputed, k)
    return result
