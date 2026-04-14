"""
LOESS (LOcally Estimated Scatterplot Smoothing) normalization.

Corrects intensity-dependent systematic bias between samples using
local regression on MA plots.  For each sample, the bias curve
M(A) = log2(sample/reference) is fitted via LOWESS and subtracted.

Reference
---------
Cleveland WS. Robust locally weighted regression and smoothing
scatterplots. J Am Stat Assoc. 1979;74(368):829-36.
Yang YH, et al. Normalization for cDNA microarray data: a robust
composite method addressing single and multiple slide systematic
variation. Nucleic Acids Res. 2002;30(4):e15.

"""

import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess

from mokume.core.logger import get_logger

logger = get_logger("mokume.normalization.loess")


class LOESSNormalizer:
    """LOESS regression-based normalization."""

    def __init__(self, frac: float = 0.75, reference: str = "median"):
        """Initialize the LOESS normalizer parameters."""
        self.frac = frac
        self.reference = reference

    def fit_transform(self, data: pd.DataFrame) -> pd.DataFrame:
        """Normalize the data using LOESS."""
        if self.reference == "median":
            ref = data.median(axis=1)
        elif self.reference == "mean":
            ref = data.mean(axis=1)
        else:
            raise ValueError(f"Unknown reference: {self.reference}")

        result = data.copy()
        n_corrected = 0

        for col in data.columns:
            sample = data[col]
            valid = sample.notna() & ref.notna()

            if valid.sum() < 10:
                continue

            m_values = (sample[valid] - ref[valid]).values  # M = log-ratio
            a_values = ((sample[valid] + ref[valid]) / 2.0).values  # A = avg

            # Fit LOWESS: M ~ A
            fitted = lowess(
                m_values, a_values,
                frac=self.frac,
                return_sorted=False,
            )

            # Subtract fitted bias
            result.loc[valid, col] = sample[valid] - fitted
            n_corrected += 1

        logger.info(
            "LOESS normalization: %d samples corrected (frac=%.2f)",
            n_corrected, self.frac,
        )
        return result


def loess_normalize(
    data: pd.DataFrame,
    frac: float = 0.75,
    reference: str = "median",
) -> pd.DataFrame:
    """Apply LOESS normalization through the functional interface."""
    return LOESSNormalizer(frac=frac, reference=reference).fit_transform(data)
