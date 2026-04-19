"""
Run/Sample QC preprocessing filters.
"""

from typing import Tuple

import pandas as pd

from mokume.core.logger import get_logger
from mokume.core.constants import (
    INTENSITY,
    NORM_INTENSITY,
    SAMPLE_ID,
    RUN,
    PROTEIN_NAME,
)
from mokume.preprocessing.filters.base import BaseFilter, FilterResult
from mokume.preprocessing.filters.enums import FilterLevel


logger = get_logger("mokume.preprocessing.filters.run_qc")


class RunIntensityFilter(BaseFilter):
    """Filter runs/samples by total intensity."""

    def __init__(
        self,
        min_intensity: float,
        intensity_column: str = NORM_INTENSITY,
        run_column: str = SAMPLE_ID,
    ):
        """
        Initialize the filter.

        Parameters
        ----------
        min_intensity : float
            Minimum total intensity for a run to be included.
        intensity_column : str, optional
            Column name containing intensity values.
        run_column : str, optional
            Column name for run/sample identifiers.
        """
        self.min_intensity = min_intensity
        self.intensity_column = intensity_column
        self.run_column = run_column

    @property
    def name(self) -> str:
        return "RunIntensityFilter"

    @property
    def level(self) -> FilterLevel:
        return FilterLevel.RUN

    def apply(self, df: pd.DataFrame, **kwargs) -> Tuple[pd.DataFrame, FilterResult]:
        input_count = len(df)

        col = self.intensity_column
        if col not in df.columns and NORM_INTENSITY in df.columns:
            col = NORM_INTENSITY
        elif col not in df.columns and INTENSITY in df.columns:
            col = INTENSITY

        run_col = self.run_column
        if run_col not in df.columns and SAMPLE_ID in df.columns:
            run_col = SAMPLE_ID
        elif run_col not in df.columns and RUN in df.columns:
            run_col = RUN

        if col not in df.columns or run_col not in df.columns:
            logger.warning("%s: Required columns not found, skipping filter", self.name)
            return df, self._create_result(input_count, input_count)

        # Calculate total intensity per run
        run_totals = df.groupby(run_col)[col].sum()
        passing_runs = run_totals[run_totals >= self.min_intensity].index
        removed_runs = run_totals[run_totals < self.min_intensity].index.tolist()

        filtered_df = df[df[run_col].isin(passing_runs)].copy()

        output_count = len(filtered_df)

        if removed_runs:
            logger.info(
                "%s: Removed runs with total intensity < %.2e: %s",
                self.name,
                self.min_intensity,
                removed_runs,
            )

        return filtered_df, self._create_result(
            input_count,
            output_count,
            {"min_intensity": self.min_intensity, "removed_runs": removed_runs},
        )


class MinFeaturesFilter(BaseFilter):
    """Filter runs/samples by minimum number of identified features."""

    def __init__(
        self,
        min_features: int = 0,
        min_proteins: int = 0,
        run_column: str = SAMPLE_ID,
        protein_column: str = PROTEIN_NAME,
    ):
        """
        Initialize the filter.

        Parameters
        ----------
        min_features : int, optional
            Minimum number of identified features per run.
        min_proteins : int, optional
            Minimum number of identified proteins per run.
        run_column : str, optional
            Column name for run/sample identifiers.
        protein_column : str, optional
            Column name containing protein identifiers.
        """
        self.min_features = min_features
        self.min_proteins = min_proteins
        self.run_column = run_column
        self.protein_column = protein_column

    @property
    def name(self) -> str:
        return "MinFeaturesFilter"

    @property
    def level(self) -> FilterLevel:
        return FilterLevel.RUN

    def apply(self, df: pd.DataFrame, **kwargs) -> Tuple[pd.DataFrame, FilterResult]:
        input_count = len(df)

        run_col = self.run_column
        if run_col not in df.columns and SAMPLE_ID in df.columns:
            run_col = SAMPLE_ID
        elif run_col not in df.columns and RUN in df.columns:
            run_col = RUN

        if run_col not in df.columns:
            logger.warning("%s: Run column not found, skipping filter", self.name)
            return df, self._create_result(input_count, input_count)

        # Count features per run
        feature_counts = df.groupby(run_col).size()
        passing_runs_features = feature_counts[
            feature_counts >= self.min_features
        ].index

        # Count proteins per run if column exists
        if self.protein_column in df.columns and self.min_proteins > 0:
            protein_counts = df.groupby(run_col)[self.protein_column].nunique()
            passing_runs_proteins = protein_counts[
                protein_counts >= self.min_proteins
            ].index
            passing_runs = passing_runs_features.intersection(passing_runs_proteins)
        else:
            passing_runs = passing_runs_features

        removed_runs = set(df[run_col].unique()) - set(passing_runs)
        filtered_df = df[df[run_col].isin(passing_runs)].copy()

        output_count = len(filtered_df)

        if removed_runs:
            logger.info(
                "%s: Removed runs with < %d features or < %d proteins: %s",
                self.name,
                self.min_features,
                self.min_proteins,
                list(removed_runs),
            )

        return filtered_df, self._create_result(
            input_count,
            output_count,
            {
                "min_features": self.min_features,
                "min_proteins": self.min_proteins,
                "removed_runs": list(removed_runs),
            },
        )


class MissingRateFilter(BaseFilter):
    """Filter runs/samples by missing value rate."""

    def __init__(
        self,
        max_missing_rate: float = 1.0,
        intensity_column: str = NORM_INTENSITY,
        run_column: str = SAMPLE_ID,
    ):
        """
        Initialize the filter.

        Parameters
        ----------
        max_missing_rate : float, optional
            Maximum fraction of missing values allowed per run (0.0-1.0).
        intensity_column : str, optional
            Column name containing intensity values.
        run_column : str, optional
            Column name for run/sample identifiers.
        """
        self.max_missing_rate = max_missing_rate
        self.intensity_column = intensity_column
        self.run_column = run_column

    @property
    def name(self) -> str:
        return "MissingRateFilter"

    @property
    def level(self) -> FilterLevel:
        return FilterLevel.RUN

    def apply(self, df: pd.DataFrame, **kwargs) -> Tuple[pd.DataFrame, FilterResult]:
        input_count = len(df)

        col = self.intensity_column
        if col not in df.columns and NORM_INTENSITY in df.columns:
            col = NORM_INTENSITY

        run_col = self.run_column
        if run_col not in df.columns and SAMPLE_ID in df.columns:
            run_col = SAMPLE_ID

        if col not in df.columns or run_col not in df.columns:
            logger.warning("%s: Required columns not found, skipping filter", self.name)
            return df, self._create_result(input_count, input_count)

        # Calculate missing rate per run
        def missing_rate(group):
            return (group[col].isna() | (group[col] == 0)).mean()

        run_missing = df.groupby(run_col).apply(missing_rate)
        passing_runs = run_missing[run_missing <= self.max_missing_rate].index
        removed_runs = run_missing[run_missing > self.max_missing_rate].index.tolist()

        filtered_df = df[df[run_col].isin(passing_runs)].copy()

        output_count = len(filtered_df)

        if removed_runs:
            logger.info(
                "%s: Removed runs with missing rate > %.1f%%: %s",
                self.name,
                self.max_missing_rate * 100,
                removed_runs,
            )

        return filtered_df, self._create_result(
            input_count,
            output_count,
            {"max_missing_rate": self.max_missing_rate, "removed_runs": removed_runs},
        )


class SampleCorrelationFilter(BaseFilter):
    """Filter samples by pairwise correlation with other samples in the condition."""

    def __init__(
        self,
        min_correlation: float,
        intensity_column: str = NORM_INTENSITY,
        sample_column: str = SAMPLE_ID,
        protein_column: str = PROTEIN_NAME,
    ):
        """
        Initialize the filter.

        Parameters
        ----------
        min_correlation : float
            Minimum pairwise correlation between samples.
        intensity_column : str, optional
            Column name containing intensity values.
        sample_column : str, optional
            Column name for sample identifiers.
        protein_column : str, optional
            Column name containing protein identifiers.
        """
        self.min_correlation = min_correlation
        self.intensity_column = intensity_column
        self.sample_column = sample_column
        self.protein_column = protein_column

    @property
    def name(self) -> str:
        return "SampleCorrelationFilter"

    @property
    def level(self) -> FilterLevel:
        return FilterLevel.SAMPLE

    def apply(self, df: pd.DataFrame, **kwargs) -> Tuple[pd.DataFrame, FilterResult]:
        input_count = len(df)

        col = self.intensity_column
        if col not in df.columns and NORM_INTENSITY in df.columns:
            col = NORM_INTENSITY

        sample_col = self.sample_column
        if sample_col not in df.columns and SAMPLE_ID in df.columns:
            sample_col = SAMPLE_ID

        if (
            col not in df.columns
            or sample_col not in df.columns
            or self.protein_column not in df.columns
        ):
            logger.warning("%s: Required columns not found, skipping filter", self.name)
            return df, self._create_result(input_count, input_count)

        # Pivot to wide format for correlation calculation
        pivot_df = df.pivot_table(
            index=self.protein_column,
            columns=sample_col,
            values=col,
            aggfunc="mean",
        )

        # Calculate pairwise correlations
        corr_matrix = pivot_df.corr()

        # Find samples with low average correlation
        avg_corr = corr_matrix.mean()
        removed_samples = avg_corr[avg_corr < self.min_correlation].index.tolist()

        filtered_df = df[~df[sample_col].isin(removed_samples)].copy()

        output_count = len(filtered_df)

        if removed_samples:
            logger.info(
                "%s: Removed samples with avg correlation < %.2f: %s",
                self.name,
                self.min_correlation,
                removed_samples,
            )

        return filtered_df, self._create_result(
            input_count,
            output_count,
            {
                "min_correlation": self.min_correlation,
                "removed_samples": removed_samples,
            },
        )
