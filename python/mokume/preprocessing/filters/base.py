"""
Base class for preprocessing filters.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Optional, Tuple

import pandas as pd

from mokume.core.logger import get_logger
from mokume.preprocessing.filters.enums import FilterLevel


logger = get_logger("mokume.preprocessing.filters")


@dataclass
class FilterResult:
    """
    Result of applying a filter.

    Attributes
    ----------
    input_count : int
        Number of items before filtering.
    output_count : int
        Number of items after filtering.
    removed_count : int
        Number of items removed.
    filter_name : str
        Name of the filter that was applied.
    filter_level : FilterLevel
        Level at which filtering was applied.
    details : dict, optional
        Additional details about the filter operation.
    """

    input_count: int
    output_count: int
    removed_count: int
    filter_name: str
    filter_level: FilterLevel
    details: Optional[dict] = field(default_factory=dict)

    @property
    def removal_rate(self) -> float:
        """Calculate the fraction of items removed."""
        if self.input_count == 0:
            return 0.0
        return self.removed_count / self.input_count

    def __repr__(self) -> str:
        return (
            f"FilterResult({self.filter_name}: "
            f"{self.removed_count}/{self.input_count} removed "
            f"({self.removal_rate:.1%}))"
        )


class BaseFilter(ABC):
    """
    Abstract base class for preprocessing filters.

    All filter implementations should inherit from this class and implement
    the apply() method.
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """Return the filter name."""
        pass

    @property
    @abstractmethod
    def level(self) -> FilterLevel:
        """Return the filter level (FEATURE, PEPTIDE, PROTEIN, RUN)."""
        pass

    @abstractmethod
    def apply(self, df: pd.DataFrame, **kwargs) -> Tuple[pd.DataFrame, FilterResult]:
        """
        Apply the filter to a DataFrame.

        Parameters
        ----------
        df : pd.DataFrame
            Input DataFrame to filter.
        **kwargs
            Additional arguments specific to the filter.

        Returns
        -------
        Tuple[pd.DataFrame, FilterResult]
            Filtered DataFrame and filter result metadata.
        """
        pass

    def _create_result(
        self,
        input_count: int,
        output_count: int,
        details: Optional[dict] = None,
    ) -> FilterResult:
        """Helper to create a FilterResult."""
        return FilterResult(
            input_count=input_count,
            output_count=output_count,
            removed_count=input_count - output_count,
            filter_name=self.name,
            filter_level=self.level,
            details=details or {},
        )

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"
