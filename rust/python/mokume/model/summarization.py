"""
Summarization strategy enumerations for the mokume package.

This module provides an enumeration of intensity summarization strategies
used when aggregating peptide intensities to protein levels.
"""

from enum import Enum, auto

import pandas as pd


class SummarizationMethod(Enum):
    """
    Enumeration of intensity summarization methods.

    These methods define how peptide intensities are aggregated when
    summarizing to protein-level quantification.

    Attributes
    ----------
    MEDIAN : auto
        Use the median of peptide intensities.
    MEAN : auto
        Use the mean (average) of peptide intensities.
    SUM : auto
        Use the sum of peptide intensities.
    MAX : auto
        Use the maximum peptide intensity.
    """

    MEDIAN = auto()
    MEAN = auto()
    SUM = auto()
    MAX = auto()

    @classmethod
    def from_str(cls, name: str) -> "SummarizationMethod":
        """
        Convert a string to a SummarizationMethod.

        Parameters
        ----------
        name : str
            The name of the summarization method.

        Returns
        -------
        SummarizationMethod
            The summarization method.

        Raises
        ------
        KeyError
            If the name does not match any summarization method.
        """
        name_ = name.lower()
        for k, v in cls._member_map_.items():
            if k.lower() == name_:
                return v
        raise KeyError(name)

    def aggregate(self, values: pd.Series) -> float:
        """
        Aggregate a series of values using this summarization method.

        Parameters
        ----------
        values : pd.Series
            The values to aggregate.

        Returns
        -------
        float
            The aggregated value.
        """
        if self == SummarizationMethod.MEDIAN:
            return values.median()
        elif self == SummarizationMethod.MEAN:
            return values.mean()
        elif self == SummarizationMethod.SUM:
            return values.sum()
        elif self == SummarizationMethod.MAX:
            return values.max()
        else:
            raise ValueError(f"Unknown summarization method: {self}")

    @property
    def pandas_agg_func(self) -> str:
        """
        Get the pandas aggregation function name for this method.

        Returns
        -------
        str
            The pandas aggregation function name.
        """
        mapping = {
            SummarizationMethod.MEDIAN: "median",
            SummarizationMethod.MEAN: "mean",
            SummarizationMethod.SUM: "sum",
            SummarizationMethod.MAX: "max",
        }
        return mapping.get(self, "sum")

    @property
    def description(self) -> str:
        """
        Get a human-readable description of the summarization method.

        Returns
        -------
        str
            Description of the summarization method.
        """
        descriptions = {
            SummarizationMethod.MEDIAN: "Median of peptide intensities",
            SummarizationMethod.MEAN: "Mean (average) of peptide intensities",
            SummarizationMethod.SUM: "Sum of peptide intensities",
            SummarizationMethod.MAX: "Maximum peptide intensity",
        }
        return descriptions.get(self, "Unknown method")
