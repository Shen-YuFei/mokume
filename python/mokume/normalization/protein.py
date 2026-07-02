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
    Apply quantile normalization to a wide protein/peptide x sample matrix.

    Each sample (column) is mapped onto a shared mean reference distribution at
    every value's rank fraction, the standard quantile normalization used by
    limma ``normalizeQuantiles`` / preprocessCore. Columns with differing
    numbers of observed (non-NaN) values are handled by interpolating the rank
    fraction in ``[0, 1]`` into the reference; tied values share their
    average-rank fraction; missing values stay missing. This matches the Rust
    ``mokume`` port cell-for-cell.

    Parameters
    ----------
    df : pd.DataFrame
        Intensity matrix in wide format (rows = proteins/peptides,
        columns = samples).
    value_column : str, optional
        Retained for API compatibility; the wide matrix carries one sample per
        column so no single value column is selected. Defaults to NORM_INTENSITY.

    Returns
    -------
    pd.DataFrame
        The DataFrame with quantile-normalized intensities.
    """
    sorted_columns = [
        np.sort(df[column].dropna().to_numpy(dtype=float)) for column in df.columns
    ]
    grid_size = max((sorted_values.size for sorted_values in sorted_columns), default=0)
    if grid_size == 0:
        return df.copy()

    reference = _mean_reference_distribution(sorted_columns, grid_size)
    normalized = df.copy()
    for column in df.columns:
        normalized[column] = _quantile_map_column(df[column], reference)
    return normalized


def _mean_reference_distribution(sorted_columns: list, grid_size: int) -> np.ndarray:
    """Cross-sample mean of each column's quantile function on a uniform
    ``[0, 1]`` grid of ``grid_size`` points."""
    grid_fractions = (
        np.arange(grid_size, dtype=float) / (grid_size - 1)
        if grid_size > 1
        else np.zeros(1)
    )
    reference = np.zeros(grid_size)
    contributing = 0
    for sorted_values in sorted_columns:
        if sorted_values.size == 0:
            continue
        reference += _interpolate_sorted(sorted_values, grid_fractions)
        contributing += 1
    return reference / contributing


def _quantile_map_column(series: pd.Series, reference: np.ndarray) -> np.ndarray:
    """Map a sample column's observed values onto ``reference`` at their
    average-rank fractions, leaving missing values untouched."""
    observed_count = int(series.notna().sum())
    out = series.to_numpy(dtype=float).copy()
    if observed_count == 0:
        return out
    if observed_count > 1:
        fractions = (
            (series.rank(method="average") - 1.0) / (observed_count - 1)
        ).to_numpy()
    else:
        fractions = np.zeros(series.shape[0])
    mask = series.notna().to_numpy()
    out[mask] = _interpolate_sorted(reference, fractions[mask])
    return out


def _interpolate_sorted(sorted_values: np.ndarray, fractions: np.ndarray) -> np.ndarray:
    """Linear interpolation of an ascending-sorted array at fractions in
    ``[0, 1]``, where a fraction maps to fractional index
    ``fraction * (len - 1)``."""
    fractions = np.asarray(fractions, dtype=float)
    if sorted_values.size == 0:
        return np.full(fractions.shape, np.nan)
    if sorted_values.size == 1:
        return np.full(fractions.shape, float(sorted_values[0]))
    positions = np.clip(fractions, 0.0, 1.0) * (sorted_values.size - 1)
    return np.interp(positions, np.arange(sorted_values.size), sorted_values)


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
