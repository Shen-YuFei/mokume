"""Count-based protein quantification methods.

``SpectralCountQuantification`` counts PSM rows, whereas
``PeptideCountQuantification`` counts distinct canonical peptide sequences.
Keeping both contracts explicit prevents feature-level input from being
misreported as a spectral count.

References
----------
- Liu H, Sadygov RG, Yates JR. A model for random sampling and
  estimation of relative protein abundance in shotgun proteomics.
  Anal Chem. 2004;76(14):4193-201.
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


@PluginRegistry.register("quantification", "spectral_count")
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

        return self._aggregate(peptide_df, group_cols, peptide_column)

    @staticmethod
    def _aggregate(
        peptide_df: pd.DataFrame,
        group_cols: list[str],
        peptide_column: str,
    ) -> pd.DataFrame:
        """Aggregate PSM rows into one count per protein and sample or run."""
        del peptide_column
        return peptide_df.groupby(group_cols).size().reset_index(name="Intensity")


@PluginRegistry.register("quantification", "peptide_count")
class PeptideCountQuantification(SpectralCountQuantification):
    """Count distinct canonical peptides per protein and sample."""

    @property
    def name(self) -> str:
        return "PeptideCount"

    @staticmethod
    def _aggregate(
        peptide_df: pd.DataFrame,
        group_cols: list[str],
        peptide_column: str,
    ) -> pd.DataFrame:
        """Count unique peptide sequences per protein and sample or run."""
        return (
            peptide_df.groupby(group_cols, dropna=False)[peptide_column]
            .nunique()
            .reset_index(name="Intensity")
        )
