"""
Missing value imputation methods.

This module provides functions for imputing missing values in proteomics data
using various methods including KNN, mean, median, and constant imputation.
"""

from typing import Optional, Union, List

import pandas as pd
from sklearn.impute import KNNImputer, SimpleImputer


def impute_missing_values(
    data: Optional[Union[pd.DataFrame, List[pd.DataFrame], None]],
    method: str = "knn",
    n_neighbors: int = 5,
    weights: str = "uniform",
    metric: str = "nan_euclidean",
    keep_empty_features: bool = True,
    fill_value: float = 0.0,
) -> Union[pd.DataFrame, List[pd.DataFrame], None]:
    """
    Impute missing values in a DataFrame or a list of DataFrames.

    Parameters
    ----------
    data : Optional[Union[pd.DataFrame, List[pd.DataFrame]]]
        A pandas DataFrame or a list of pandas DataFrames containing missing values.
    method : str, optional
        The imputation method to use. Options are:
        - "knn" (default): Use K-Nearest Neighbors imputation.
        - "mean": Impute using the mean of each column.
        - "median": Impute using the median of each column.
        - "most_frequent": Impute using the most frequent value of each column.
        - "constant": Impute using a specific value provided via `fill_value`.
        - "seqknn": Sequential KNN (Kim et al. 2004) - iterative neighbour fill.
        - "missforest": Random Forest iterative imputation (Stekhoven & Buehlmann 2011).
    n_neighbors : int, optional
        The number of neighboring samples to use for KNN imputation. Default is 5.
    weights : str, optional
        The weight function used in KNN prediction. Default is 'uniform'.
    metric : str, optional
        The distance metric used for finding neighbors in KNN. Default is 'nan_euclidean'.
    fill_value : float, optional
        The constant value to use for imputation when `method` is "constant". Default is 0.0.
    keep_empty_features : bool, optional
        Whether to keep features that are entirely empty. Default is True.

    Returns
    -------
    Union[pd.DataFrame, List[pd.DataFrame], None]
        A pandas DataFrame or a list of pandas DataFrames with imputed missing values.
    """
    if data is None:
        return None

    if method not in {
        "knn",
        "mean",
        "median",
        "constant",
        "most_frequent",
        "seqknn",
        "missforest",
    }:
        raise ValueError(
            "Invalid method. Choose from 'knn', 'mean', 'median', 'most_frequent', "
            "'constant', 'seqknn', or 'missforest'."
        )

    if method == "seqknn":
        from mokume.imputation.seqknn import impute_seqknn

        def impute(df: pd.DataFrame) -> pd.DataFrame:
            return impute_seqknn(df, k=n_neighbors)
    elif method == "missforest":
        from mokume.imputation.missforest import impute_missforest

        def impute(df: pd.DataFrame) -> pd.DataFrame:
            return impute_missforest(df)
    else:
        if method == "knn":
            imputer = KNNImputer(
                n_neighbors=n_neighbors,
                weights=weights,
                metric=metric,
                keep_empty_features=keep_empty_features,
            )
        else:
            strategy = method
            imputer = SimpleImputer(strategy=strategy, fill_value=fill_value)

        def impute(df: pd.DataFrame) -> pd.DataFrame:
            imputed_data = imputer.fit_transform(df)
            return pd.DataFrame(imputed_data, columns=df.columns, index=df.index)

    if isinstance(data, pd.DataFrame):
        return impute(data)
    elif isinstance(data, list) and all(isinstance(df, pd.DataFrame) for df in data):
        return [impute(df) for df in data]
    else:
        raise ValueError(
            "The input data must be a pandas DataFrame, a list of DataFrames, or None."
        )
