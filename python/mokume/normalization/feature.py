"""
Feature-level normalization implementations.

This module provides functions for normalizing feature intensities
across technical replicates and samples.
"""

import pandas as pd

from mokume.model.normalization import FeatureNormalizationMethod


def normalize_replicates(
    df: pd.DataFrame, method: FeatureNormalizationMethod, *args, **kwargs
) -> pd.Series:
    """
    Normalize the replicate intensities using the specified method.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame containing replicate intensity data.
    method : FeatureNormalizationMethod
        The normalization method to use.
    *args
        Additional positional arguments.
    **kwargs
        Additional keyword arguments.

    Returns
    -------
    pd.Series
        The normalized replicate intensities.
    """
    return method.normalize_replicates(df, *args, **kwargs)


def normalize_sample(
    df: pd.DataFrame, runs: list[str], method: FeatureNormalizationMethod
) -> tuple[dict[str, pd.Series], float]:
    """
    Normalize replicate intensities for a given sample across multiple runs.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame containing replicate intensity data.
    runs : list[str]
        A list of run identifiers for the sample.
    method : FeatureNormalizationMethod
        The normalization method to use.

    Returns
    -------
    tuple[dict[str, pd.Series], float]
        A dictionary mapping each run to its normalized replicate intensities
        and the average metric across all runs.
    """
    return method.normalize_sample(df, runs)


def normalize_runs(
    df: pd.DataFrame,
    technical_replicates: int,
    method: FeatureNormalizationMethod = FeatureNormalizationMethod.Median,
) -> pd.DataFrame:
    """
    Normalize the intensities of runs in the given DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame containing replicate intensity data.
    technical_replicates : int
        The number of technical replicates for each sample.
    method : FeatureNormalizationMethod, optional
        The normalization method to use. Defaults to Median.

    Returns
    -------
    pd.DataFrame
        The DataFrame with normalized replicate intensities.
    """
    return method.normalize_runs(df, technical_replicates)
