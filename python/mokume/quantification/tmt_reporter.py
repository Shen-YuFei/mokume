"""
TMT reporter intensity protein quantification method.

This module provides the TMT reporter intensity quantification method,
which computes protein abundance as the sum of raw peptide reporter
intensities per (protein, sample). Also known as the ``intensity`` method
and stays in linear space (no log transform).
"""

from typing import Optional

import pandas as pd

from mokume.quantification.base import ProteinQuantificationMethod
from mokume.core.registry import PluginRegistry
from mokume.core.constants import (
    PROTEIN_NAME,
    PEPTIDE_CANONICAL,
    NORM_INTENSITY,
    SAMPLE_ID,
)


@PluginRegistry.register("quantification", "tmt_reporter")
class TMTReporterIntensityQuantification(ProteinQuantificationMethod):
    """
    TMT reporter intensity quantification (``intensity`` method).

    Protein abundance is computed as the sum of raw peptide reporter
    intensities for each protein in each sample. The output stays in
    linear space (no log transform) so downstream code can apply its
    preferred normalization or transformation.
    """

    @property
    def name(self) -> str:
        return "TMTReporterIntensity"

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
        Quantify proteins as the sum of raw reporter intensities.

        Parameters
        ----------
        peptide_df : pd.DataFrame
            DataFrame containing peptide-level data.
        protein_column : str
            Column name for protein identifiers.
        peptide_column : str
            Column name for peptide sequences.
        intensity_column : str
            Column name for intensity values (linear scale).
        sample_column : str
            Column name for sample identifiers.
        run_column : str, optional
            Column name for run identifiers. If provided, quantification
            is performed at the run level instead of sample level.

        Returns
        -------
        pd.DataFrame
            DataFrame with columns: protein_column, sample_column,
            (run_column if provided), 'Intensity'. The 'Intensity'
            column contains the linear-space sum of reporter intensities
            per protein.
        """
        if run_column is not None and run_column in peptide_df.columns:
            group_cols = [protein_column, sample_column, run_column]
        else:
            group_cols = [protein_column, sample_column]

        result = peptide_df.groupby(group_cols)[intensity_column].sum().reset_index()
        return result.rename(columns={intensity_column: "Intensity"})
