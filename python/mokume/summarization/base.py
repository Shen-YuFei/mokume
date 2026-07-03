"""
Intensity summarization strategies for the mokume package.

This module provides the base class and implementations for various summarization
strategies used when aggregating peptide intensities to protein levels.
"""

from abc import ABC, abstractmethod
import pandas as pd


class SummarizationStrategy(ABC):
    """Abstract base class for intensity summarization strategies."""

    @property
    @abstractmethod
    def name(self) -> str:
        """Return the name of the summarization strategy."""
        pass

    @abstractmethod
    def aggregate(self, values: pd.Series) -> float:
        """
        Aggregate a series of values.

        Parameters
        ----------
        values : pd.Series
            Series of intensity values to aggregate.

        Returns
        -------
        float
            Aggregated value.
        """
        pass

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"


class MedianSummarization(SummarizationStrategy):
    """Use median of peptide intensities."""

    @property
    def name(self) -> str:
        return "Median"

    def aggregate(self, values: pd.Series) -> float:
        return values.median()


class MeanSummarization(SummarizationStrategy):
    """Use mean of peptide intensities."""

    @property
    def name(self) -> str:
        return "Mean"

    def aggregate(self, values: pd.Series) -> float:
        return values.mean()


class SumSummarization(SummarizationStrategy):
    """Use sum of peptide intensities."""

    @property
    def name(self) -> str:
        return "Sum"

    def aggregate(self, values: pd.Series) -> float:
        return values.sum()
