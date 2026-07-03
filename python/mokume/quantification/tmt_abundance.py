"""
TMT abundance protein quantification method.

This module provides the TMT abundance quantification method, which
computes protein abundance as the median of log2-transformed peptide
intensities per (protein, sample). Also known as the ``abd`` method,
designed for TMT reporter ion data.
"""

from typing import Optional

import numpy as np
import pandas as pd

from mokume.quantification.base import ProteinQuantificationMethod
from mokume.core.registry import PluginRegistry
from mokume.core.constants import (
    PROTEIN_NAME,
    PEPTIDE_CANONICAL,
    NORM_INTENSITY,
    SAMPLE_ID,
)


@PluginRegistry.register("quantification", "tmt_abundance")
class TMTAbundanceQuantification(ProteinQuantificationMethod):
    """
    TMT abundance quantification (``abd`` method).

    Protein abundance is computed as the median of log2-transformed
    peptide intensities for each protein in each sample. Non-positive
    intensities are treated as missing (log2 undefined) and excluded
    before median aggregation.
    """

    @property
    def name(self) -> str:
        return "TMTAbundance"

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
        Quantify proteins as the median of log2 peptide intensities.

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
            column contains the median log2 abundance per protein.
        """
        if run_column is not None and run_column in peptide_df.columns:
            group_cols = [protein_column, sample_column, run_column]
        else:
            group_cols = [protein_column, sample_column]

        intensities = peptide_df[intensity_column].astype(float)
        log2_intensities = np.log2(intensities.where(intensities > 0))

        log2_df = peptide_df[group_cols].copy()
        log2_df["_log2_intensity"] = log2_intensities

        result = (
            log2_df.dropna(subset=["_log2_intensity"])
            .groupby(group_cols)["_log2_intensity"]
            .median()
            .reset_index()
        )
        return result.rename(columns={"_log2_intensity": "Intensity"})
