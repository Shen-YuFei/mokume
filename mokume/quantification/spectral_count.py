"""
Spectral count protein quantification method.

Protein abundance is measured as the number of identified spectra
(PSMs) per protein per sample.  This is the simplest count-based
quantification and serves as a baseline for label-free experiments.

References
----------
- Liu H, Sadygov RG, Yates JR. A model for random sampling and
  estimation of relative protein abundance in shotgun proteomics.
  Anal Chem. 2004;76(14):4193-201.
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


class SpectralCountQuantification(ProteinQuantificationMethod):
    """
    Spectral count quantification.

    Protein abundance is the number of peptide-spectrum matches (PSMs)
    per protein per sample.
    """

    @property
    def name(self) -> str:
        return "SpectralCount"

    def quantify(
        self,
        peptide_df: pd.DataFrame,
        protein_column: str = PROTEIN_NAME,
        peptide_column: str = PEPTIDE_CANONICAL,
        intensity_column: str = NORM_INTENSITY,
        sample_column: str = SAMPLE_ID,
        run_column: Optional[str] = None,
    ) -> pd.DataFrame:
        """Count spectra per protein per sample.

        Parameters
        ----------
        peptide_df : pd.DataFrame
            DataFrame containing peptide-level data.  Each row is one
            PSM; the method counts rows per (protein, sample).
        protein_column : str
            Column name for protein identifiers.
        peptide_column : str
            Unused (kept for API compatibility).
        intensity_column : str
            Unused (kept for API compatibility).
        sample_column : str
            Column name for sample identifiers.
        run_column : str, optional
            If provided, counting is done per run instead of per sample.

        Returns
        -------
        pd.DataFrame
            DataFrame with columns: *protein_column*, *sample_column*,
            (optionally *run_column*), and ``Intensity`` (the count).
        """
        if run_column is not None and run_column in peptide_df.columns:
            group_cols = [protein_column, sample_column, run_column]
        else:
            group_cols = [protein_column, sample_column]

        result = peptide_df.groupby(group_cols).size().reset_index(name="Intensity")
        return result
