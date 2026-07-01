"""
Intensity-based preprocessing filters.
"""

from typing import Optional, Tuple, List

import pandas as pd
import numpy as np

from mokume.core.logger import get_logger
from mokume.core.constants import INTENSITY, NORM_INTENSITY, SAMPLE_ID, CONDITION
from mokume.preprocessing.filters.base import BaseFilter, FilterResult
from mokume.preprocessing.filters.enums import FilterLevel


logger = get_logger("mokume.preprocessing.filters.intensity")


class MinIntensityFilter(BaseFilter):
    """Filter features below a minimum intensity threshold."""

    def __init__(self, min_intensity: float, intensity_column: str = INTENSITY):
        """
        Initialize the filter.

        Parameters
        ----------
        min_intensity : float
            Minimum intensity threshold. Features below this are removed.
        intensity_column : str, optional
            Column name containing intensity values.
        """
        self.min_intensity = min_intensity
        self.intensity_column = intensity_column

    @property
    def name(self) -> str:
        return "MinIntensityFilter"

    @property
    def level(self) -> FilterLevel:
        return FilterLevel.FEATURE

    def apply(self, df: pd.DataFrame, **kwargs) -> Tuple[pd.DataFrame, FilterResult]:
        input_count = len(df)

        col = self.intensity_column
        if col not in df.columns and NORM_INTENSITY in df.columns:
            col = NORM_INTENSITY
        elif col not in df.columns and INTENSITY in df.columns:
            col = INTENSITY

        if col not in df.columns:
            logger.warning(
                "%s: Intensity column '%s' not found, skipping filter",
                self.name,
                self.intensity_column,
            )
            return df, self._create_result(input_count, input_count)

        mask = df[col] >= self.min_intensity
        filtered_df = df[mask].copy()

        output_count = len(filtered_df)

        logger.debug(
            "%s: Removed %d features below intensity %.2f",
            self.name,
            input_count - output_count,
            self.min_intensity,
        )

        return filtered_df, self._create_result(
            input_count, output_count, {"threshold": self.min_intensity}
        )


class CVThresholdFilter(BaseFilter):
    """Filter features with CV above threshold across replicates."""

    def __init__(
        self,
        cv_threshold: float,
        intensity_column: str = NORM_INTENSITY,
        groupby_columns: Optional[List[str]] = None,
    ):
        """
        Initialize the filter.

        Parameters
        ----------
        cv_threshold : float
            Maximum coefficient of variation (CV) threshold.
        intensity_column : str, optional
            Column name containing intensity values.
        groupby_columns : list[str], optional
            Columns to group by for CV calculation.
        """
        self.cv_threshold = cv_threshold
        self.intensity_column = intensity_column
        self.groupby_columns = groupby_columns or ["ProteinName", "PeptideCanonical"]

    @property
    def name(self) -> str:
        return "CVThresholdFilter"

    @property
    def level(self) -> FilterLevel:
        return FilterLevel.PEPTIDE

    def apply(self, df: pd.DataFrame, **kwargs) -> Tuple[pd.DataFrame, FilterResult]:
        input_count = len(df)

        col = self.intensity_column
        if col not in df.columns and NORM_INTENSITY in df.columns:
            col = NORM_INTENSITY

        # Check if groupby columns exist
        valid_groupby = [c for c in self.groupby_columns if c in df.columns]
        if not valid_groupby:
            logger.warning(
                "%s: No valid groupby columns found, skipping filter", self.name
            )
            return df, self._create_result(input_count, input_count)

        # Calculate CV per peptide
        def calc_cv(x):
            if x.mean() > 0 and len(x) > 1:
                return x.std() / x.mean()
            return np.nan

        cv_df = df.groupby(valid_groupby)[col].agg(calc_cv).reset_index()
        cv_df.columns = list(valid_groupby) + ["cv"]

        # Identify peptides below CV threshold (or with NaN CV - single measurements)
        passing = cv_df[(cv_df["cv"] <= self.cv_threshold) | cv_df["cv"].isna()][
            valid_groupby
        ]

        # Filter original DataFrame
        filtered_df = df.merge(passing, on=valid_groupby, how="inner")

        output_count = len(filtered_df)

        logger.debug(
            "%s: Removed %d features with CV > %.2f",
            self.name,
            input_count - output_count,
            self.cv_threshold,
        )

        return filtered_df, self._create_result(
            input_count, output_count, {"cv_threshold": self.cv_threshold}
        )


class ReplicateAgreementFilter(BaseFilter):
    """Filter features not detected in minimum number of replicates."""

    def __init__(
        self,
        min_replicates: int,
        sample_column: str = SAMPLE_ID,
        condition_column: str = CONDITION,
        groupby_columns: Optional[List[str]] = None,
    ):
        """
        Initialize the filter.

        Parameters
        ----------
        min_replicates : int
            Minimum number of replicates where a feature must be detected.
        sample_column : str, optional
            Column name for sample identifiers.
        condition_column : str, optional
            Column name for condition/group identifiers.
        groupby_columns : list[str], optional
            Columns to group by for replicate counting.
        """
        self.min_replicates = min_replicates
        self.sample_column = sample_column
        self.condition_column = condition_column
        self.groupby_columns = groupby_columns or ["ProteinName", "PeptideCanonical"]

    @property
    def name(self) -> str:
        return "ReplicateAgreementFilter"

    @property
    def level(self) -> FilterLevel:
        return FilterLevel.PEPTIDE

    def apply(self, df: pd.DataFrame, **kwargs) -> Tuple[pd.DataFrame, FilterResult]:
        input_count = len(df)

        # Determine columns for grouping
        valid_groupby = [c for c in self.groupby_columns if c in df.columns]
        if not valid_groupby:
            logger.warning(
                "%s: No valid groupby columns found, skipping filter", self.name
            )
            return df, self._create_result(input_count, input_count)

        # Add condition column if present
        if self.condition_column in df.columns:
            group_cols = [self.condition_column] + valid_groupby
        else:
            group_cols = valid_groupby

        # Count samples per peptide per condition
        if self.sample_column in df.columns:
            rep_counts = (
                df.groupby(group_cols)[self.sample_column].nunique().reset_index()
            )
            rep_counts.columns = list(group_cols) + ["rep_count"]

            # Filter for minimum replicates
            passing = rep_counts[rep_counts["rep_count"] >= self.min_replicates][
                group_cols
            ]
            filtered_df = df.merge(passing, on=group_cols, how="inner")
        else:
            logger.warning(
                "%s: Sample column '%s' not found, skipping filter",
                self.name,
                self.sample_column,
            )
            filtered_df = df

        output_count = len(filtered_df)

        logger.debug(
            "%s: Removed %d features with < %d replicates",
            self.name,
            input_count - output_count,
            self.min_replicates,
        )

        return filtered_df, self._create_result(
            input_count, output_count, {"min_replicates": self.min_replicates}
        )


class QuantileFilter(BaseFilter):
    """Filter features outside specified quantile range."""

    def __init__(
        self,
        lower_quantile: float = 0.0,
        upper_quantile: float = 1.0,
        intensity_column: str = NORM_INTENSITY,
    ):
        """
        Initialize the filter.

        Parameters
        ----------
        lower_quantile : float, optional
            Lower quantile threshold (0.0-1.0).
        upper_quantile : float, optional
            Upper quantile threshold (0.0-1.0).
        intensity_column : str, optional
            Column name containing intensity values.
        """
        self.lower_quantile = lower_quantile
        self.upper_quantile = upper_quantile
        self.intensity_column = intensity_column

    @property
    def name(self) -> str:
        return "QuantileFilter"

    @property
    def level(self) -> FilterLevel:
        return FilterLevel.FEATURE

    def apply(self, df: pd.DataFrame, **kwargs) -> Tuple[pd.DataFrame, FilterResult]:
        input_count = len(df)

        col = self.intensity_column
        if col not in df.columns and NORM_INTENSITY in df.columns:
            col = NORM_INTENSITY

        if col not in df.columns:
            logger.warning(
                "%s: Intensity column '%s' not found, skipping filter",
                self.name,
                self.intensity_column,
            )
            return df, self._create_result(input_count, input_count)

        lower_bound = df[col].quantile(self.lower_quantile)
        upper_bound = df[col].quantile(self.upper_quantile)

        mask = (df[col] >= lower_bound) & (df[col] <= upper_bound)
        filtered_df = df[mask].copy()

        output_count = len(filtered_df)

        logger.debug(
            "%s: Removed %d features outside quantile range [%.3f, %.3f]",
            self.name,
            input_count - output_count,
            self.lower_quantile,
            self.upper_quantile,
        )

        return filtered_df, self._create_result(
            input_count,
            output_count,
            {
                "lower_quantile": self.lower_quantile,
                "upper_quantile": self.upper_quantile,
                "lower_bound": float(lower_bound),
                "upper_bound": float(upper_bound),
            },
        )
