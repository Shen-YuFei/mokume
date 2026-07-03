"""
Sample-level distribution-alignment normalizers.

This module provides object-oriented wrappers around
``quantile_normalize`` and ``median_center`` from
:mod:`mokume.normalization.protein`, plus a new
``MeanCenterNormalizer`` that subtracts each sample's mean. All three
classes follow the sklearn-style ``fit`` / ``transform`` pattern used
elsewhere in :mod:`mokume.normalization` (see ``RlrNormalizer``,
``IRSNormalizer``).

These methods are dataset-level: they require all samples at once.
"""

from typing import Optional

import pandas as pd

from mokume.normalization.protein import median_center, quantile_normalize


class _StatelessNormalizer:
    """Base class for normalizers without per-sample fitted state."""

    def __init__(self) -> None:
        self._fitted = False

    def fit(self, X: pd.DataFrame, y=None) -> "_StatelessNormalizer":
        """Validate input and mark as fitted."""
        if X.empty:
            raise ValueError("Input DataFrame is empty")
        self._fitted = True
        return self

    def transform(self, X: pd.DataFrame) -> pd.DataFrame:
        """Apply normalization; raises if ``fit`` was not called first."""
        if not self._fitted:
            raise ValueError(
                f"{self.__class__.__name__} has not been fitted. Call fit() first."
            )
        return self._apply(X)

    def fit_transform(self, X: pd.DataFrame, y=None) -> pd.DataFrame:
        """Fit and transform in one step."""
        return self.fit(X, y).transform(X)

    def _apply(self, x: pd.DataFrame) -> pd.DataFrame:  # pragma: no cover - abstract
        raise NotImplementedError


class QuantileNormalizer(_StatelessNormalizer):
    """
    Quantile normalization across samples.

    Wraps :func:`mokume.normalization.protein.quantile_normalize`.

    Parameters
    ----------
    value_column : str, optional
        Column name passed through to the underlying ``quantile_normalize``
        for long-format DataFrames; ignored when ``X`` is wide (proteins
        x samples).
    """

    def __init__(self, value_column: Optional[str] = None):
        super().__init__()
        self.value_column = value_column

    def _apply(self, x: pd.DataFrame) -> pd.DataFrame:
        if self.value_column is not None and self.value_column in x.columns:
            return quantile_normalize(x, value_column=self.value_column)
        return quantile_normalize(x)


class MedianCenterNormalizer(_StatelessNormalizer):
    """
    Sample-level median centering.

    Subtracts each sample's median from its values. Wraps
    :func:`mokume.normalization.protein.median_center` with ``axis=0``
    (the default), which centers each column (sample).
    """

    def _apply(self, x: pd.DataFrame) -> pd.DataFrame:
        return median_center(x, axis=0)


class MeanCenterNormalizer(_StatelessNormalizer):
    """
    Sample-level mean centering.

    Subtracts each sample's mean from its values. Operates on a wide
    DataFrame (proteins x samples). Distinct from
    ``FeatureNormalizationMethod.Mean``, which divides feature-level
    technical replicates by their mean (different scope and operation).
    """

    def _apply(self, x: pd.DataFrame) -> pd.DataFrame:
        return x - x.mean(axis=0)
