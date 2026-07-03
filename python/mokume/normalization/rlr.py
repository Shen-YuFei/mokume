"""
RLR (Robust Linear Regression) normalization.

This module provides RLR normalization, which applies a robust linear
regression in log space against the per-protein median across all
samples, then rescales each sample to remove the fitted intercept and
slope. RLR originates from NormalyzerDE
(Chawade et al., 2014).

References
----------
- Chawade, A., Alexandersson, E., Levander, F. (2014). NormalyzerDE: a
  data-driven web tool for normalization of proteomics datasets.
  Journal of Proteome Research, 13(6), 3114-3120.
"""

from typing import Optional

import numpy as np
import pandas as pd

from mokume.core.logger import get_logger

logger = get_logger("mokume.normalization.rlr")


class RlrNormalizer:
    """
    Robust Linear Regression (RLR) normalization in log space.

    For each sample, an M-estimator regression is fitted between the
    sample's log2 intensities and the per-protein median (used as a
    pseudo-reference). The fitted intercept and slope are then used
    to rescale the sample so that residuals are centered around zero
    and unit slope.

    Parameters
    ----------
    log_transform : bool, default=False
        If True, input is already log2-transformed.
        If False (default), expects raw intensities and will log2-transform
        internally before fitting (zeros are dropped per sample).

    Attributes
    ----------
    intercepts_ : pd.Series
        Fitted intercept per sample (in log2 space).
    slopes_ : pd.Series
        Fitted slope per sample (in log2 space).
    reference_ : pd.Series
        Per-protein log2 median used as the pseudo-reference.
    """

    def __init__(self, log_transform: bool = False):
        self.log_transform = log_transform
        self.intercepts_: Optional[pd.Series] = None
        self.slopes_: Optional[pd.Series] = None
        self.reference_: Optional[pd.Series] = None

    def _to_log2(self, x: pd.DataFrame) -> pd.DataFrame:
        if self.log_transform:
            return x.replace([np.inf, -np.inf], np.nan)
        positive = x.where(x > 0)
        return np.log2(positive)

    def _fit_sample(
        self, sample_log: pd.Series, reference: pd.Series
    ) -> tuple[float, float]:
        valid = sample_log.notna() & reference.notna()
        if valid.sum() < 5:
            logger.warning(
                "RLR: sample %s has only %d valid proteins; using identity fit",
                sample_log.name,
                int(valid.sum()),
            )
            return 0.0, 1.0
        try:
            import statsmodels.api as sm
        except ImportError as exc:
            raise ImportError(
                "statsmodels is required for RLR normalization. "
                "Install it with: pip install mokume[analysis]"
            ) from exc
        x = sm.add_constant(reference[valid].values)
        y = sample_log[valid].values
        try:
            result = sm.RLM(y, x, M=sm.robust.norms.HuberT()).fit()
        except (np.linalg.LinAlgError, ValueError) as exc:
            logger.warning("RLR fit failed for %s: %s", sample_log.name, exc)
            return 0.0, 1.0
        intercept, slope = float(result.params[0]), float(result.params[1])
        if not np.isfinite(slope) or abs(slope) < 1e-6:
            logger.warning(
                "RLR: degenerate slope for %s (slope=%.3g); using identity",
                sample_log.name,
                slope,
            )
            return 0.0, 1.0
        return intercept, slope

    def fit(self, X: pd.DataFrame, y=None) -> "RlrNormalizer":
        """Compute per-sample RLR intercept and slope."""
        if X.empty:
            raise ValueError("Input DataFrame is empty")
        log_x = self._to_log2(X)
        self.reference_ = log_x.median(axis=1, skipna=True)
        intercepts = {}
        slopes = {}
        for sample in log_x.columns:
            intercepts[sample], slopes[sample] = self._fit_sample(
                log_x[sample], self.reference_
            )
        self.intercepts_ = pd.Series(intercepts)
        self.slopes_ = pd.Series(slopes)
        logger.info("RLR fitted on %d samples", len(log_x.columns))
        return self

    def transform(self, X: pd.DataFrame) -> pd.DataFrame:
        """Apply RLR normalization in log2 space and return the same scale as input."""
        if self.intercepts_ is None or self.slopes_ is None:
            raise ValueError("RlrNormalizer has not been fitted. Call fit() first.")
        log_x = self._to_log2(X)
        result_log = pd.DataFrame(index=log_x.index, columns=log_x.columns, dtype=float)
        for sample in log_x.columns:
            if sample not in self.intercepts_.index:
                logger.warning("Sample %s not in fitted RLR factors", sample)
                result_log[sample] = log_x[sample]
                continue
            alpha = self.intercepts_[sample]
            beta = self.slopes_[sample]
            result_log[sample] = (log_x[sample] - alpha) / beta
        return result_log if self.log_transform else 2**result_log

    def fit_transform(self, X: pd.DataFrame, y=None) -> pd.DataFrame:
        """Fit and transform in one step."""
        return self.fit(X, y).transform(X)


def rlr_normalize(df: pd.DataFrame, log_transform: bool = False) -> pd.DataFrame:
    """Convenience wrapper around :class:`RlrNormalizer.fit_transform`."""
    return RlrNormalizer(log_transform=log_transform).fit_transform(df)
