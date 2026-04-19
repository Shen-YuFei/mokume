"""
All peptides protein quantification method.

This module provides a quantification method that sums all peptide
intensities for each protein.
"""

from typing import Optional

import pandas as pd

from mokume.quantification.base import ProteinQuantificationMethod
from mokume.core.constants import (
    PROTEIN_NAME,
    PEPTIDE_CANONICAL,
    NORM_INTENSITY,
    SAMPLE_ID,
)


class AllPeptidesQuantification(ProteinQuantificationMethod):
    """
    All peptides protein quantification method.

    Calculates protein abundance as the sum of all peptide intensities
    for each protein in each sample (or run if run_column is provided).
    """

    @property
    def name(self) -> str:
        return "AllPeptides"

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
        Quantify proteins using sum of all peptide intensities.

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

        Returns
        -------
        pd.DataFrame
            DataFrame with columns: protein_column, sample_column,
            (run_column if provided), 'Intensity'.
        """
        # Determine grouping columns based on aggregation level
        if run_column is not None and run_column in peptide_df.columns:
            group_cols = [protein_column, sample_column, run_column]
        else:
            group_cols = [protein_column, sample_column]

        result = peptide_df.groupby(group_cols)[intensity_column].sum().reset_index()

        # Rename intensity column
        result = result.rename(columns={intensity_column: "Intensity"})
        return result
