"""
Base class for protein quantification methods.

This module provides an abstract base class that all quantification methods
should inherit from.
"""

from abc import ABC, abstractmethod
from typing import Optional

import pandas as pd

from mokume.core.constants import (
    PROTEIN_NAME,
    PEPTIDE_CANONICAL,
    NORM_INTENSITY,
    SAMPLE_ID,
)


INTENSITY_COLUMN = "Intensity"
"""Standard output intensity column name for all quantification methods."""


class ProteinQuantificationMethod(ABC):
    """
    Abstract base class for protein quantification methods.

    All quantification methods (iBAQ, Top3, TopN, MaxLFQ, etc.) should
    inherit from this class and implement the required methods.

    All methods produce a standard output column named 'Intensity'.
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """Return the name of the quantification method."""
        pass

    @abstractmethod
    def quantify(
        self,
        peptide_df: pd.DataFrame,
        protein_column: str = PROTEIN_NAME,
        peptide_column: str = PEPTIDE_CANONICAL,
        intensity_column: str = NORM_INTENSITY,
        sample_column: str = SAMPLE_ID,
        run_column: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Quantify proteins from peptide intensities.

        Parameters
        ----------
        peptide_df : pd.DataFrame
            DataFrame containing peptide-level data.
        protein_column : str
            Column name for protein identifiers.
        peptide_column : str
            Column name for peptide sequences.
        intensity_column : str
            Column name for intensity values.
        sample_column : str
            Column name for sample identifiers.
        run_column : str, optional
            Column name for run identifiers. If provided, quantification
            is performed at the run level instead of sample level.
            This enables run-level aggregation similar to DIA-NN's approach.

        Returns
        -------
        pd.DataFrame
            DataFrame containing protein-level quantification values.
            If run_column is provided, the result will contain both
            sample and run identifiers.
        """
        pass

    def _get_grouping_column(
        self, sample_column: str, run_column: Optional[str] = None
    ) -> str:
        """
        Get the column to use for grouping based on aggregation level.

        Parameters
        ----------
        sample_column : str
            Column name for sample identifiers.
        run_column : str, optional
            Column name for run identifiers.

        Returns
        -------
        str
            The column name to use for grouping (run_column if provided,
            otherwise sample_column).
        """
        return run_column if run_column is not None else sample_column

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"
