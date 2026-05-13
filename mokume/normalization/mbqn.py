"""
MBQN (Mean/Median-Balanced Quantile Normalization).

Performs standard quantile normalization then re-balances each
feature (protein) so that its across-sample mean matches the
mean of the pre-normalized values, preserving rank-invariant
properties while reducing the feature-level bias introduced by
vanilla quantile normalization.

References
----------
- Brombacher E, et al. Tail-robust quantile normalization.
  Proteomics. 2020;20(10):1900068.
"""

from __future__ import annotations

import pandas as pd

from mokume.core.logger import get_logger
from mokume.normalization.protein import quantile_normalize

logger = get_logger("mokume.normalization.mbqn")


def mbqn_normalize(data: pd.DataFrame) -> pd.DataFrame:
    """Apply Mean-Balanced Quantile Normalization.

    Parameters
    ----------
    data : pd.DataFrame
        Wide matrix (rows=proteins, columns=samples) in log2 space.

    Returns
    -------
    pd.DataFrame
        MBQN-normalized matrix with same shape and index.
    """
    row_means_before = data.mean(axis=1)

    qn = quantile_normalize(data)

    row_means_after = qn.mean(axis=1)
    shift = row_means_before - row_means_after

    result = qn.add(shift, axis=0)
    logger.info("MBQN normalization complete: %d proteins", len(result))
    return result


class MBQNNormalizer:
    """Scikit-learn-style wrapper around :func:`mbqn_normalize`."""

    def __init__(self) -> None:
        self._fitted = False

    def fit(self, data: pd.DataFrame, y=None):  # noqa: ARG002
        """Fit (no-op for stateless normalizer)."""
        self._fitted = True
        return self

    def transform(self, data: pd.DataFrame) -> pd.DataFrame:
        """Apply MBQN normalization."""
        return mbqn_normalize(data)

    def fit_transform(self, data: pd.DataFrame, y=None) -> pd.DataFrame:
        """Fit and transform in one step."""
        return self.fit(data, y).transform(data)
