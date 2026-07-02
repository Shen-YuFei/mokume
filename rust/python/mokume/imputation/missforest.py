"""
missForest imputation (Stekhoven & Buehlmann 2011) via sklearn.

The original ``missForest`` is an R package by Stekhoven that
iteratively imputes missing values using random forest predictions.
There is no widely-adopted Python re-implementation by the original
author, so we recreate the algorithm with sklearn's
``IterativeImputer`` driven by a ``RandomForestRegressor``. This
mirrors the Stekhoven & Buehlmann (2011) approach for continuous
proteomics intensities while keeping mokume's main dependency set
small (no extra C++ build chains).

References
----------
- Stekhoven, D. J., & Buehlmann, P. (2012). MissForest -
  non-parametric missing value imputation for mixed-type data.
  *Bioinformatics*, 28(1), 112-118.
"""

from typing import Optional

import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.experimental import enable_iterative_imputer  # noqa: F401  pylint: disable=unused-import
from sklearn.impute import IterativeImputer

from mokume.core.logger import get_logger

logger = get_logger("mokume.imputation.missforest")


def impute_missforest(
    data: pd.DataFrame,
    n_estimators: int = 100,
    max_iter: int = 10,
    random_state: Optional[int] = 0,
    n_jobs: int = -1,
) -> pd.DataFrame:
    """Iteratively impute missing values with a Random Forest regressor.

    Parameters mirror sklearn defaults: ``n_estimators`` trees per
    forest, ``max_iter`` imputation passes, ``random_state`` for
    reproducibility, ``n_jobs`` parallel workers (-1 = all cores).
    Returns a same-shape DataFrame with NaNs filled.
    """
    if data.empty or not data.isna().any().any():
        return data.copy()
    estimator = RandomForestRegressor(
        n_estimators=n_estimators,
        n_jobs=n_jobs,
        random_state=random_state,
    )
    imputer = IterativeImputer(
        estimator=estimator,
        max_iter=max_iter,
        random_state=random_state,
    )
    imputed = imputer.fit_transform(data.to_numpy(dtype=float))
    logger.info(
        "missForest: imputed %d values (n_estimators=%d, max_iter=%d)",
        int(data.isna().sum().sum()),
        n_estimators,
        max_iter,
    )
    return pd.DataFrame(imputed, index=data.index, columns=data.columns)
