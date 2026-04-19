"""
Protein-level normalization implementations.

This module provides functions for normalizing protein-level intensities.
"""

import pandas as pd
import numpy as np

from mokume.core.constants import NORM_INTENSITY


def quantile_normalize(
    df: pd.DataFrame, value_column: str = NORM_INTENSITY
) -> pd.DataFrame:
    """
    Apply quantile normalization to protein intensities.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame containing protein intensity data in wide format.
    value_column : str, optional
        The column name for intensity values. Defaults to NORM_INTENSITY.

    Returns
    -------
    pd.DataFrame
        The DataFrame with quantile-normalized intensities.
    """
    # Rank values within each column
    ranked = df.rank(method="average")

    # Get the mean of each rank across all columns
    rank_mean = (
        ranked.stack().groupby(ranked.rank(method="first").stack().values).mean()
    )

    # Map the rank means back to the original positions
    normalized = (
        ranked.stack()
        .map(lambda x: rank_mean[x] if not pd.isna(x) else np.nan)
        .unstack()
    )

    return normalized


def median_center(df: pd.DataFrame, axis: int = 0) -> pd.DataFrame:
    """
    Center data by subtracting the median.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to center.
    axis : int, optional
        Axis along which to compute the median. 0 for columns, 1 for rows.

    Returns
    -------
    pd.DataFrame
        The centered DataFrame.
    """
    median = df.median(axis=axis)
    if axis == 0:
        return df - median
    else:
        return df.sub(median, axis=0)


def total_intensity_normalize(
    df: pd.DataFrame, target_sum: float = 1e9
) -> pd.DataFrame:
    """
    Normalize by total intensity (column sums).

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to normalize.
    target_sum : float, optional
        The target sum for each column. Defaults to 1e9.

    Returns
    -------
    pd.DataFrame
        The normalized DataFrame.
    """
    col_sums = df.sum(axis=0)
    return df * (target_sum / col_sums)
