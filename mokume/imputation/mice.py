"""
MICE (Multiple Imputation by Chained Equations) for proteomics data.

Wraps scikit-learn's ``IterativeImputer`` which implements the MICE
algorithm using BayesianRidge as the default estimator.

References
----------
- van Buuren S, Groothuis-Oudshoorn K. mice: Multivariate
  Imputation by Chained Equations in R. J Stat Softw. 2011;45(3).
"""

from __future__ import annotations

import pandas as pd

from mokume.core.logger import get_logger

logger = get_logger("mokume.imputation.mice")


def impute_mice(
    data: pd.DataFrame,
    max_iter: int = 10,
    n_nearest_features: int | None = None,
    random_state: int = 42,
) -> pd.DataFrame:
    """Impute missing values using MICE (IterativeImputer).

    Parameters
    ----------
    data : pd.DataFrame
        Protein intensity matrix (rows=proteins, columns=samples) with NaN.
    max_iter : int
        Maximum number of imputation rounds (default 10).
    n_nearest_features : int or None
        Number of other features to use for each estimation step.
        ``None`` means use all features.
    random_state : int
        Seed for reproducibility.

    Returns
    -------
    pd.DataFrame
        Filled matrix with the same shape and index.
    """
    from sklearn.experimental import enable_iterative_imputer  # noqa: F401
    from sklearn.impute import IterativeImputer

    n_missing = int(data.isna().sum().sum())
    imputer = IterativeImputer(
        max_iter=max_iter,
        n_nearest_features=n_nearest_features,
        random_state=random_state,
    )
    filled = pd.DataFrame(
        imputer.fit_transform(data),
        index=data.index,
        columns=data.columns,
    )
    logger.info("MICE imputation: %d values imputed", n_missing)
    return filled
