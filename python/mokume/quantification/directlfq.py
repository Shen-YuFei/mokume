"""
DirectLFQ protein quantification wrapper.

This module provides a wrapper around the DirectLFQ package for protein
quantification using intensity traces and hierarchical normalization.

DirectLFQ is an optional dependency. Install with:
    pip install mokume[directlfq]

Reference:
    Ammar C, et al. Accurate label-free quantification by directLFQ to compare
    unlimited numbers of proteomes. Mol Cell Proteomics. 2023.
"""

import gc
import os
import tempfile

import pandas as pd
from typing import Optional

from mokume.quantification.base import ProteinQuantificationMethod
from mokume.core.logger import get_logger
from mokume.core.registry import PluginRegistry
from mokume.core.constants import (
    PROTEIN_NAME,
    PEPTIDE_CANONICAL,
    NORM_INTENSITY,
    SAMPLE_ID,
)

logger = get_logger("mokume.quantification.directlfq")

# Lazy import flag
_DIRECTLFQ_AVAILABLE = None


def _check_directlfq_available():
    """Check if directlfq is installed and importable."""
    global _DIRECTLFQ_AVAILABLE
    if _DIRECTLFQ_AVAILABLE is None:
        import importlib.util

        _DIRECTLFQ_AVAILABLE = importlib.util.find_spec("directlfq") is not None
    return _DIRECTLFQ_AVAILABLE


def _import_directlfq():
    """Import directlfq with helpful error message if not available."""
    if not _check_directlfq_available():
        raise ImportError(
            "DirectLFQ support requires the 'directlfq' package.\n"
            "Install with: pip install mokume[directlfq]\n"
            "Or: pip install directlfq"
        )
    import directlfq.lfq_manager as lfq_manager

    return lfq_manager


@PluginRegistry.register("quantification", "directlfq")
class DirectLFQQuantification(ProteinQuantificationMethod):
    """
    DirectLFQ protein quantification using intensity traces.

    DirectLFQ uses a hierarchical normalization approach with variance-guided
    pairwise alignment to produce accurate label-free quantification results.
    It is particularly effective for large-scale proteomics datasets.

    This is a wrapper around the DirectLFQ package by Mann Labs.

    Parameters
    ----------
    min_nonan : int
        Minimum number of non-NaN ion intensities required per protein.
        Default is 1.
    num_cores : int, optional
        Number of CPU cores for parallel processing.
        Default is None (auto-detect).
    deactivate_normalization : bool
        If True, skip the normalization step.
        Default is False.

    Examples
    --------
    >>> from mokume.quantification import DirectLFQQuantification
    >>> directlfq = DirectLFQQuantification(min_nonan=2)
    >>> result = directlfq.quantify(
    ...     peptide_df,
    ...     protein_column="ProteinName",
    ...     peptide_column="PeptideSequence",
    ...     intensity_column="Intensity",
    ...     sample_column="SampleID"
    ... )

    Notes
    -----
    DirectLFQ is an optional dependency. Install with:
        pip install mokume[directlfq]

    References
    ----------
    Ammar C, et al. Accurate label-free quantification by directLFQ to compare
    unlimited numbers of proteomes. Mol Cell Proteomics. 2023.
    """

    def __init__(
        self,
        min_nonan: int = 1,
        num_cores: Optional[int] = None,
        deactivate_normalization: bool = False,
    ):
        """
        Initialize DirectLFQ quantification.

        Parameters
        ----------
        min_nonan : int
            Minimum number of non-NaN ion intensities per protein.
        num_cores : int, optional
            Number of CPU cores for parallel processing.
        deactivate_normalization : bool
            Whether to skip normalization.
        """
        # Validate directlfq is available at init time
        _import_directlfq()

        self.min_nonan = min_nonan
        self.num_cores = num_cores
        self.deactivate_normalization = deactivate_normalization

    @property
    def name(self) -> str:
        return "DirectLFQ"

    @property
    def input_level(self) -> str:
        """DirectLFQ consumes raw per-ion feature intensities.

        Maps to the ``directlfq`` flow in the pipeline runner's
        ``FLOW_DISPATCH``.
        """
        return "peptides_raw"

    def _prepare_input_file(
        self,
        peptide_df: pd.DataFrame,
        protein_column: str,
        peptide_column: str,
        intensity_column: str,
        sample_column: str,
        temp_dir: str,
    ) -> str:
        """
        Convert peptide DataFrame to DirectLFQ input format.

        DirectLFQ expects a matrix format with:
        - 'protein' column
        - 'ion' column (peptide identifier)
        - One column per sample with intensity values

        Parameters
        ----------
        peptide_df : pd.DataFrame
            Input peptide data in long format.
        protein_column : str
            Column name for protein identifiers.
        peptide_column : str
            Column name for peptide sequences.
        intensity_column : str
            Column name for intensity values.
        sample_column : str
            Column name for sample identifiers.
        temp_dir : str
            Directory for temporary file.

        Returns
        -------
        str
            Path to the prepared input file.
        """
        # Pivot to wide format (samples as columns). observed=True keeps
        # pandas from materialising the Cartesian product of every Categorical
        # level (LoadingStage casts SAMPLE_ID/CONDITION to Categorical for
        # legacy dtype parity; on PXD030304 the unobserved product reached
        # ~7.3e12 cells, ~13 TiB of int16, before this guard).
        pivot_df = peptide_df.pivot_table(
            index=[protein_column, peptide_column],
            columns=sample_column,
            values=intensity_column,
            aggfunc="sum",  # Sum if multiple measurements per peptide/sample
            observed=True,
        ).reset_index()

        # Rename columns to DirectLFQ expected format
        pivot_df = pivot_df.rename(
            columns={protein_column: "protein", peptide_column: "ion"}
        )

        # Save to temporary file
        input_path = os.path.join(temp_dir, "directlfq_input.aq_reformat.tsv")
        pivot_df.to_csv(input_path, sep="\t", index=False)

        return input_path

    def _parse_wide_output(
        self,
        result_wide: pd.DataFrame,
        protein_column: str,
        sample_column: str,
    ) -> pd.DataFrame:
        sample_cols = [c for c in result_wide.columns if c != "protein"]
        result_long = result_wide.melt(
            id_vars=["protein"],
            value_vars=sample_cols,
            var_name=sample_column,
            value_name="Intensity",
        )
        result_long = result_long.rename(columns={"protein": protein_column})
        return result_long[
            result_long["Intensity"].notna() & (result_long["Intensity"] > 0)
        ]

    def _parse_output(
        self,
        output_path: str,
        protein_column: str,
        sample_column: str,
    ) -> pd.DataFrame:
        """
        Parse DirectLFQ output and convert to long format.

        Parameters
        ----------
        output_path : str
            Path to DirectLFQ protein intensities output.
        protein_column : str
            Original protein column name to use in output.
        sample_column : str
            Original sample column name to use in output.

        Returns
        -------
        pd.DataFrame
            Results in long format with protein, sample, and intensity columns.
        """
        result_wide = pd.read_csv(output_path, sep="\t")
        return self._parse_wide_output(result_wide, protein_column, sample_column)

    def _run_streamed_directlfq(self, input_path: str) -> pd.DataFrame:
        lfq_manager = _import_directlfq()
        import directlfq.config as lfq_config
        import directlfq.normalization as lfq_norm
        import directlfq.utils as lfq_utils

        from mokume.pipeline.directlfq_streaming import (
            estimate_protein_intensities_streamed,
        )

        lfq_config.set_global_protein_and_ion_id(protein_id="protein", quant_id="ion")
        lfq_config.set_log_processed_proteins(log_processed_proteins=True)
        lfq_config.set_compile_normalized_ion_table(False)
        lfq_config.check_wether_to_copy_numpy_arrays_derived_from_pandas()

        input_file = lfq_manager.prepare_input_filename(input_path)
        input_df = lfq_utils.import_data(
            input_file=input_file,
            input_type_to_use=None,
            filter_dict=None,
        )
        input_df = lfq_utils.sort_input_df_by_protein_and_quant_id(input_df)
        input_df = lfq_utils.remove_potential_quant_id_duplicates(input_df)
        input_df = lfq_utils.index_and_log_transform_input_df(input_df)
        input_df = lfq_utils.remove_allnan_rows_input_df(input_df)

        if not self.deactivate_normalization:
            logger.info("Performing DirectLFQ sample normalization")
            input_df = lfq_norm.NormalizationManagerSamplesOnSelectedProteins(
                input_df,
                num_samples_quadratic=50,
            ).complete_dataframe
            gc.collect()

        logger.info("Estimating DirectLFQ protein intensities")
        return estimate_protein_intensities_streamed(
            input_df,
            min_nonan=self.min_nonan,
            num_samples_quadratic=10,
            num_cores=self.num_cores,
        )

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
        Quantify proteins using DirectLFQ algorithm.

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

        Returns
        -------
        pd.DataFrame
            DataFrame with columns: protein_column, sample_column, 'Intensity'.
        """
        _import_directlfq()

        logger.info("Running DirectLFQ quantification")
        logger.info("Input: %d peptide measurements", len(peptide_df))
        logger.info("Proteins: %d", peptide_df[protein_column].nunique())
        logger.info("Samples: %d", peptide_df[sample_column].nunique())

        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = self._prepare_input_file(
                peptide_df,
                protein_column,
                peptide_column,
                intensity_column,
                sample_column,
                temp_dir,
            )
            logger.info("Prepared DirectLFQ input: %s", input_path)

            try:
                result_wide = self._run_streamed_directlfq(input_path)
            except Exception as e:
                logger.error("DirectLFQ failed: %s", e)
                raise RuntimeError(f"DirectLFQ quantification failed: {e}") from e

            result_df = self._parse_wide_output(
                result_wide,
                protein_column,
                sample_column,
            )

        logger.info(
            "DirectLFQ complete: %d proteins",
            result_df[protein_column].nunique(),
        )

        return result_df


def is_directlfq_available() -> bool:
    """
    Check if DirectLFQ is available for use.

    Returns
    -------
    bool
        True if directlfq package is installed, False otherwise.

    Examples
    --------
    >>> from mokume.quantification.directlfq import is_directlfq_available
    >>> if is_directlfq_available():
    ...     from mokume.quantification import DirectLFQQuantification
    ...     quant = DirectLFQQuantification()
    """
    return _check_directlfq_available()
