"""
Concrete feature-level (within-run) normalization methods.

Each class registers with the PluginRegistry and implements
``transform_series()`` — the per-run normalization math.
The run-balancing orchestration is handled by the base class.
"""

import pandas as pd

from mokume.core.registry import PluginRegistry
from mokume.normalization.base import FeatureNormalizer


@PluginRegistry.register("normalization.feature", "none")
class NoneFeatureNormalizer(FeatureNormalizer):
    """No-op normalizer — returns data unchanged."""

    @property
    def name(self) -> str:
        return "none"

    def transform_series(self, series: pd.Series) -> float:
        # normalize() is overridden to a no-op; return a harmless scalar so
        # the base-class per-run pathway (if ever reached) does nothing.
        return 1.0

    def normalize(self, df, intensity_col=None, group_col=None, sample_col=None):
        return df


@PluginRegistry.register("normalization.feature", "mean")
class MeanFeatureNormalizer(FeatureNormalizer):
    """Mean normalization: intensity / mean(intensity)."""

    @property
    def name(self) -> str:
        return "mean"

    def transform_series(self, series: pd.Series) -> float:
        return series.mean()


@PluginRegistry.register("normalization.feature", "median")
class MedianFeatureNormalizer(FeatureNormalizer):
    """Median normalization: intensity / median(intensity)."""

    @property
    def name(self) -> str:
        return "median"

    def transform_series(self, series: pd.Series) -> float:
        return series.median()


@PluginRegistry.register("normalization.feature", "max")
class MaxFeatureNormalizer(FeatureNormalizer):
    """Max normalization: intensity / max(intensity)."""

    @property
    def name(self) -> str:
        return "max"

    def transform_series(self, series: pd.Series) -> float:
        return series.max()


@PluginRegistry.register("normalization.feature", "global")
class GlobalFeatureNormalizer(FeatureNormalizer):
    """Global normalization: intensity / sum(intensity)."""

    @property
    def name(self) -> str:
        return "global"

    def transform_series(self, series: pd.Series) -> float:
        return series.sum()


@PluginRegistry.register("normalization.feature", "max_min")
class MaxMinFeatureNormalizer(FeatureNormalizer):
    """Max-Min normalization: (intensity - min) / (max - min)."""

    @property
    def name(self) -> str:
        return "max_min"

    def transform_series(self, series: pd.Series) -> float:
        value_range = series.max() - series.min()
        # Guard against a zero (or NaN) range: a constant run has no spread
        # to balance on, so return a neutral 1.0 to skip scaling and avoid
        # a divide-by-zero in normalize().
        if not value_range or pd.isna(value_range):
            return 1.0
        return value_range


@PluginRegistry.register("normalization.feature", "iqr")
class IQRFeatureNormalizer(FeatureNormalizer):
    """IQR normalization: interquartile range (75th - 25th quantile)."""

    @property
    def name(self) -> str:
        return "iqr"

    def transform_series(self, series: pd.Series) -> float:
        quantiles = series.quantile([0.25, 0.75], interpolation="linear")
        value_range = quantiles.loc[0.75] - quantiles.loc[0.25]
        # Guard against a zero (or NaN) IQR to avoid a divide-by-zero in
        # normalize() when a run has no spread across the quartiles.
        if not value_range or pd.isna(value_range):
            return 1.0
        return value_range
