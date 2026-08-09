"""
Unified pipeline: features → proteins in one step.

This module provides the main `features_to_proteins` function and
`QuantificationPipeline` class that handle the full proteomics
quantification workflow from feature-level parquet files to
protein intensities.

The pipeline automatically handles:
- Loading and filtering parquet data
- Normalization (run-level and sample-level)
- Protein quantification using various methods
- Optional intermediate exports (peptides, ions)

When DirectLFQ is selected as the quantification method, the pipeline uses the
DirectLFQ package for normalization and a Mokume streaming adapter for protein
estimation to avoid retaining ion-level intermediate results for every protein.
"""

import gc

import pandas as pd
from pathlib import Path
from typing import Optional

from mokume.core.constants import PROTEIN_NAME, load_sdrf
from mokume.core.dataset import QpxDataset
from mokume.core.logger import get_logger
from mokume.model.normalization import parse_normalization_methods
from mokume.pipeline.config import (
    PipelineConfig,
    InputConfig,
    FilterConfig,
    NormalizationConfig,
    QuantificationConfig,
    IRSConfig,
    BatchCorrectionConfig,
    ImputationConfig,
    DEConfig,
    OutputConfig,
    RuntimeConfig,
    validate_de_config,
)
from mokume.pipeline.directlfq_streaming import estimate_protein_intensities_streamed
from mokume.pipeline.stages import (
    LoadingStage,
    NormalizationStage,
    QuantificationStage,
    ImputationStage,
    PostprocessingStage,
)

logger = get_logger("mokume.pipeline")


class QuantificationPipeline:
    """
    Unified pipeline: features → proteins.

    Thin orchestrator that delegates to stage classes:
    - LoadingStage: Data loading and filtering
    - NormalizationStage: Run-level and sample-level normalization
    - QuantificationStage: Protein quantification methods
    - PostprocessingStage: Batch correction, DE, plotting, reports

    Parameters
    ----------
    config : PipelineConfig
        Pipeline configuration object.

    Examples
    --------
    >>> from mokume.pipeline import QuantificationPipeline, PipelineConfig
    >>> from mokume.pipeline.config import InputConfig, QuantificationConfig
    >>>
    >>> config = PipelineConfig(
    ...     input=InputConfig(parquet="data.parquet"),
    ...     quantification=QuantificationConfig(method="maxlfq"),
    ... )
    >>> pipeline = QuantificationPipeline(config)
    >>> proteins = pipeline.run()
    """

    def __init__(self, config: PipelineConfig):
        self.config = config
        self._validate_config()
        self.loading = LoadingStage(config)
        self.normalization = NormalizationStage(config)
        self.quantification = QuantificationStage(config)
        self.imputation = ImputationStage(config)
        self.postprocessing = PostprocessingStage(config)
        # Intermediate tables/flags captured during a run() for run_dataset() to
        # attach to the QpxDataset. Populated as a side effect by the
        # standard-flow / ratio-flow helpers; never changes run()'s contract.
        # Held in a single dict to keep the pipeline's attribute surface small.
        self._run_capture: dict = self._empty_capture()

    @staticmethod
    def _empty_capture() -> dict:
        """Return a fresh, empty per-run capture state."""
        return {
            "peptides": None,
            "psms": None,
            "ratio_refs": None,
            "ratio_sample_to_plex": None,
            "irs_applied": False,
        }

    def _validate_config(self):
        """Validate configuration and check for required parameters."""
        validate_de_config(self.config.de)
        parse_normalization_methods(
            self.config.normalization.run_method,
            self.config.normalization.sample_method,
        )
        source_path = self.config.input.msstats or self.config.input.parquet
        if not Path(source_path).exists():
            source_name = "MSstats" if self.config.input.msstats else "Parquet"
            raise FileNotFoundError(f"{source_name} file not found: {source_path}")

        if (
            self.config.input.msstats
            and self.config.quantification.method.lower() == "ratio"
        ):
            raise ValueError(
                "Ratio quantification requires PSM-level QPX input; "
                "MSstats feature tables do not contain the required PSM evidence"
            )

        if (
            self.config.quantification.method.lower() == "ibaq"
            and not self.config.input.fasta_file
        ):
            raise ValueError("iBAQ quantification requires --fasta-file")

        if (
            self.config.quantification.method.lower() == "ratio"
            and not self.config.input.sdrf
        ):
            raise ValueError("Ratio quantification requires an SDRF file (--sdrf)")

        if (
            self.config.input.fasta_file
            and not Path(self.config.input.fasta_file).exists()
        ):
            raise FileNotFoundError(
                f"FASTA file not found: {self.config.input.fasta_file}"
            )

    def run(self) -> pd.DataFrame:
        """
        Execute the full pipeline.

        Returns
        -------
        pd.DataFrame
            Protein intensities matrix (proteins x samples).
        """
        quant_method = self.config.quantification.method.lower()
        logger.info("Starting pipeline with quant_method=%s", quant_method)

        # Reset per-run capture state so run_dataset() reflects only this run.
        self._run_capture = self._empty_capture()

        if quant_method == "directlfq":
            protein_df = self._run_directlfq_pipeline()
        elif quant_method == "maxlfq" and self._can_run_maxlfq_directlfq_pipeline():
            protein_df = self._run_maxlfq_directlfq_pipeline()
        elif quant_method == "ratio":
            protein_df = self._run_ratio_pipeline()
        else:
            protein_df = self._run_mokume_pipeline()

        # Apply IRS normalization if configured (skip for ratio — ratios already
        # handle cross-plex normalization via per-plex reference division)
        if self.config.irs.enabled and quant_method != "ratio":
            protein_df = self.normalization.apply_irs(protein_df)
            self._run_capture["irs_applied"] = True

        # Apply coverage filter if configured (generic, works with any method)
        if self.config.quantification.coverage_threshold is not None:
            protein_df = self.normalization.apply_coverage_filter(protein_df)

        # Apply imputation if configured (after coverage filter, before batch correction)
        if self.config.imputation.enabled:
            protein_df = self.imputation.impute(protein_df)

        # Apply batch correction if configured
        if self.config.batch.enabled:
            protein_df = self.postprocessing.apply_batch_correction(protein_df)

        # Run differential expression if configured
        de_results = None
        if self.config.de.enabled:
            de_results = self.postprocessing.run_differential_expression(protein_df)

        # Generate plots if configured
        if self.config.output.plot_dir and any(
            [
                self.config.output.plot_volcano,
                self.config.output.plot_heatmap,
                self.config.output.plot_pca,
            ]
        ):
            self.postprocessing.generate_plots(protein_df, de_results)

        # Generate interactive report if configured
        if self.config.output.interactive_report and de_results:
            self.postprocessing.generate_interactive_report(protein_df, de_results)

        return protein_df

    def run_dataset(self) -> QpxDataset:
        """Execute the pipeline and return results as a :class:`QpxDataset`.

        This is a thin wrapper around :meth:`run`: it runs the exact same
        quantification workflow (so protein values are identical to the legacy
        ``run() -> DataFrame`` path), then packages the wide protein matrix and
        any intermediate tables captured during the run into a hierarchical
        :class:`~mokume.core.dataset.QpxDataset` with provenance metadata.

        The dataset levels populated depend on the flow:

        * ``proteins`` -- always set to the wide (protein x sample) matrix
          returned by :meth:`run`.
        * ``peptides`` -- the assembled, normalized peptide table, when the
          standard (mokume-native) flow produced one.
        * ``psms`` -- the loaded PSM table, for the ratio flow.
        * ``sample_info`` -- SDRF metadata when an SDRF file was configured.

        Provenance (``dataset.uns["provenance"]["steps"]``) records at least a
        ``loading`` step and a ``quantification`` step; an
        ``irs_normalization`` step is appended when IRS ran.

        Returns
        -------
        QpxDataset
            Dataset with the protein matrix, intermediate levels, sample
            metadata, and recorded provenance.
        """
        quant_method = self.config.quantification.method

        # Run the legacy pipeline unchanged; this also populates the
        # self._run_capture dict (peptides / psms / ratio refs / irs flag).
        protein_df = self.run()
        capture = self._run_capture

        dataset = QpxDataset()
        dataset.proteins = protein_df
        if capture["peptides"] is not None:
            dataset.peptides = capture["peptides"]
        if capture["psms"] is not None:
            dataset.psms = capture["psms"]

        # Sample metadata from the SDRF, when one is configured. Loading is
        # best-effort: a malformed/absent SDRF must not sink an otherwise
        # successful quantification run, so file/parse errors are logged and
        # sample_info is simply left unset.
        if self.config.input.sdrf:
            try:
                dataset.sample_info = load_sdrf(self.config.input.sdrf)
            except (OSError, ValueError) as exc:  # pragma: no cover - best-effort
                logger.debug("Could not load SDRF sample_info: %s", exc)

        # --- Provenance -------------------------------------------------
        rows_in = None
        if capture["peptides"] is not None:
            rows_in = len(capture["peptides"])
        elif capture["psms"] is not None:
            rows_in = len(capture["psms"])

        dataset.record_step(
            "loading",
            parquet=(
                str(self.config.input.parquet) if self.config.input.parquet else None
            ),
            msstats=(
                str(self.config.input.msstats) if self.config.input.msstats else None
            ),
            sdrf=str(self.config.input.sdrf) if self.config.input.sdrf else None,
            rows_out=rows_in,
        )
        dataset.record_step(
            "quantification",
            method=quant_method,
            rows_out=len(protein_df) if isinstance(protein_df, pd.DataFrame) else None,
        )
        if capture["irs_applied"]:
            dataset.record_step("irs_normalization", stat=self.config.irs.stat)

        # Ratio flow exposes its detected reference samples for downstream
        # consumers (mirrors main's cecilia integration contract).
        if quant_method.lower() == "ratio" and capture["ratio_refs"] is not None:
            sample_to_plex = capture["ratio_sample_to_plex"] or {}
            dataset.uns["ratio_config"] = {
                "reference_samples": list(capture["ratio_refs"]),
                "n_plexes": len(set(sample_to_plex.values())) if sample_to_plex else 0,
            }

        return dataset

    def _run_directlfq_pipeline(self) -> pd.DataFrame:
        """Run pipeline using DirectLFQ package."""
        try:
            import directlfq.normalization as lfq_norm
            import directlfq.config as lfq_config
        except ImportError as exc:
            raise ImportError(
                "DirectLFQ quantification requires the directlfq package.\n"
                "Install with: pip install directlfq\n"
                "Or: pip install mokume[directlfq]"
            ) from exc

        logger.info("Loading and filtering data for DirectLFQ...")
        filtered_table = self.loading.load_for_directlfq()
        logger.info("Filtered data: %d features", filtered_table.num_rows)

        logger.info("Converting to DirectLFQ format...")
        directlfq_input = self.loading.convert_to_directlfq_format(filtered_table)
        logger.info("DirectLFQ input shape: %s", directlfq_input.shape)
        del filtered_table
        gc.collect()

        # Configure DirectLFQ
        if self.config.output.export_ions:
            raise ValueError(
                "--export-ions is not supported by streaming DirectLFQ estimation; "
                "normalized ion tables can be much larger than the protein matrix."
            )
        lfq_config.set_global_protein_and_ion_id(protein_id="protein", quant_id="ion")
        lfq_config.set_compile_normalized_ion_table(False)

        # Run DirectLFQ normalization
        logger.info("Running DirectLFQ sample normalization...")
        normed_df = lfq_norm.NormalizationManagerSamplesOnSelectedProteins(
            directlfq_input,
            num_samples_quadratic=self.config.quantification.directlfq_num_samples_quadratic,
        ).complete_dataframe
        del directlfq_input
        gc.collect()

        logger.info("Running DirectLFQ protein estimation...")
        protein_df = estimate_protein_intensities_streamed(
            normed_df,
            min_nonan=self.config.quantification.directlfq_min_nonan,
            num_samples_quadratic=self.config.quantification.directlfq_num_samples_quadratic,
            num_cores=self.config.quantification.directlfq_num_cores,
        )

        logger.info("DirectLFQ complete: %d proteins", len(protein_df))
        return protein_df

    def _can_run_maxlfq_directlfq_pipeline(self) -> bool:
        """Return whether MaxLFQ can use the low-memory DirectLFQ path."""
        if self.config.output.export_peptides or self.config.output.export_ions:
            logger.info(
                "MaxLFQ DirectLFQ-streaming path disabled because intermediate "
                "peptide/ion export was requested"
            )
            return False
        if self.config.normalization.sample_method.lower() != "none":
            logger.info(
                "MaxLFQ DirectLFQ-streaming path disabled for sample normalization: %s",
                self.config.normalization.sample_method,
            )
            return False
        try:
            from mokume.quantification.directlfq import is_directlfq_available
        except ImportError:
            return False
        return is_directlfq_available()

    def _run_maxlfq_directlfq_pipeline(self) -> pd.DataFrame:
        """Run MaxLFQ through the low-memory DirectLFQ-compatible path."""
        try:
            import directlfq.config as lfq_config
            import directlfq.normalization as lfq_norm
        except ImportError as exc:
            raise ImportError(
                "MaxLFQ DirectLFQ-streaming path requires the directlfq package.\n"
                "Install with: pip install directlfq\n"
                "Or: pip install mokume[directlfq]"
            ) from exc

        logger.info("Loading and filtering data for MaxLFQ DirectLFQ-streaming...")
        filtered_table = self.loading.load_for_maxlfq_directlfq()
        logger.info("Filtered MaxLFQ input rows: %d", filtered_table.num_rows)

        logger.info("Converting MaxLFQ input to DirectLFQ wide format...")
        directlfq_input = self.loading.convert_maxlfq_to_directlfq_format(
            filtered_table
        )
        logger.info("MaxLFQ DirectLFQ input shape: %s", directlfq_input.shape)
        del filtered_table
        gc.collect()

        lfq_config.set_global_protein_and_ion_id(protein_id="protein", quant_id="ion")
        lfq_config.set_compile_normalized_ion_table(False)

        logger.info("Running DirectLFQ sample normalization for MaxLFQ...")
        normed_df = lfq_norm.NormalizationManagerSamplesOnSelectedProteins(
            directlfq_input,
            num_samples_quadratic=self.config.quantification.directlfq_num_samples_quadratic,
        ).complete_dataframe
        del directlfq_input
        gc.collect()

        logger.info("Running DirectLFQ protein estimation for MaxLFQ...")
        protein_df = estimate_protein_intensities_streamed(
            normed_df,
            min_nonan=self.config.quantification.directlfq_min_nonan,
            num_samples_quadratic=self.config.quantification.directlfq_num_samples_quadratic,
            num_cores=self.config.quantification.directlfq_num_cores,
        )
        if "protein" in protein_df.columns:
            protein_df = protein_df.rename(columns={"protein": PROTEIN_NAME})
        sample_cols = [c for c in protein_df.columns if c != PROTEIN_NAME]
        protein_df[sample_cols] = protein_df[sample_cols].replace(0, pd.NA)

        logger.info("MaxLFQ DirectLFQ-streaming complete: %d proteins", len(protein_df))
        return protein_df

    def _run_mokume_pipeline(self) -> pd.DataFrame:
        """Run pipeline using mokume's native implementations."""
        logger.info("Loading and filtering data...")
        peptide_df = self.loading.load_for_mokume()
        logger.info("Processed peptides: %d rows", len(peptide_df))

        # Capture the assembled peptide table so run_dataset() can attach it to
        # the QpxDataset peptides level without reloading the parquet.
        self._run_capture["peptides"] = peptide_df

        # Export peptides if requested
        if self.config.output.export_peptides:
            logger.info("Exporting peptides to %s", self.config.output.export_peptides)
            peptide_df.to_csv(self.config.output.export_peptides, index=False)

        # Quantify proteins
        logger.info(
            "Quantifying proteins with method: %s", self.config.quantification.method
        )
        protein_df = self.quantification.quantify(peptide_df)
        logger.info("Quantification complete: %d proteins", len(protein_df))

        return protein_df

    def _run_ratio_pipeline(self) -> pd.DataFrame:
        """Run ratio-based quantification (PS protocol)."""
        from mokume.quantification.ratio import RatioQuantification

        logger.info("Running ratio-based quantification (PS protocol)...")

        psm_df, ref_samples, sample_to_plex = self.loading.load_for_ratio()

        # Capture the PSM table and detected references so run_dataset() can
        # attach them to the QpxDataset without reloading the parquet.
        self._run_capture["psms"] = psm_df
        self._run_capture["ratio_refs"] = ref_samples
        self._run_capture["ratio_sample_to_plex"] = sample_to_plex

        # Run ratio quantification
        ratio_quant = RatioQuantification(
            reference_samples=ref_samples,
            sample_to_plex=sample_to_plex,
            fraction_merge_method=self.config.quantification.ratio_fraction_merge,
        )
        protein_df = ratio_quant.quantify(psm_df)

        # Remove reference samples from output columns (log2(ref/ref) = 0)
        protein_col = protein_df.columns[0]
        cols_to_keep = [protein_col] + [
            c for c in protein_df.columns if c == protein_col or c not in ref_samples
        ]
        # Deduplicate while preserving order
        seen = set()
        unique_cols = []
        for c in cols_to_keep:
            if c not in seen:
                seen.add(c)
                unique_cols.append(c)
        protein_df = protein_df[unique_cols]

        logger.info("Ratio pipeline complete: %d proteins", len(protein_df))
        return protein_df


def features_to_proteins(
    parquet: Optional[str],
    output: str,
    sdrf: Optional[str] = None,
    quant_method: str = "maxlfq",
    min_aa: int = 7,
    min_unique_peptides: int = 2,
    remove_contaminants: bool = True,
    run_normalization: str = "median",
    sample_normalization: str = "globalMedian",
    normalization_proteins_file: Optional[str] = None,
    fasta_file: Optional[str] = None,
    ion_alignment: Optional[str] = None,
    # piBAQ-specific (paralog-aware iBAQ)
    ibaq_enzyme: str = "Trypsin",
    ibaq_max_aa: int = 50,
    ibaq_min_shared: int = 2,
    ibaq_families_yaml: Optional[str] = None,
    ibaq_min_anchors: int = 1,
    ibaq_high_anchor_threshold: int = 3,
    directlfq_num_cores: Optional[int] = None,
    directlfq_min_nonan: int = 1,
    export_peptides: Optional[str] = None,
    export_ions: Optional[str] = None,
    # Batch correction parameters
    batch_correction: bool = False,
    batch_method: str = "sample_prefix",
    batch_column: Optional[str] = None,
    batch_covariates: Optional[list] = None,
    batch_parametric: bool = True,
    batch_mean_only: bool = False,
    batch_ref: Optional[int] = None,
    # IRS normalization parameters
    irs: bool = False,
    irs_reference_samples: Optional[list] = None,
    irs_sdrf_column: Optional[str] = None,
    irs_sdrf_values: Optional[list] = None,
    irs_reference_regex: str = "pool|powder|ref|reference|bridge",
    irs_stat: str = "median",
    irs_remove_reference: bool = False,
    # Differential expression parameters
    differential_expression: bool = False,
    de_contrasts: Optional[list] = None,
    de_method: str = "auto",
    de_log2fc_threshold: float = 0.5,
    de_fdr_threshold: float = 0.05,
    de_fdr_method: str = "bh",
    de_output: Optional[str] = None,
    de_ensemble_methods: Optional[list] = None,
    de_ensemble_min_k: int = 2,
    # Plotting parameters
    plot_output_dir: Optional[str] = None,
    plot_volcano: bool = False,
    plot_heatmap: bool = False,
    plot_pca: bool = False,
    highlight_genes: Optional[list] = None,
    # Coverage filter
    coverage_threshold: Optional[float] = None,
    # Ratio quantification
    ratio_fraction_merge: str = "mean",
    # Imputation parameters
    impute: bool = False,
    impute_method: str = "none",
    impute_quantile: float = 0.01,
    impute_shift: float = 1.6,
    impute_scale: float = 0.3,
    impute_n_neighbors: int = 5,
    # Interactive report parameters
    interactive_report: bool = False,
    report_output: Optional[str] = None,
    # DuckDB resource limits (cap DuckDB engine only, NOT total process RSS)
    duckdb_memory: Optional[str] = None,
    duckdb_threads: Optional[int] = None,
    *,
    msstats: Optional[str] = None,
) -> pd.DataFrame:
    """
    Quantify proteins from QPX feature parquet or legacy SDRF+MSstats data.

    This is the main entry point for the unified pipeline that handles
    the full workflow from features to proteins in one step.

    Parameters
    ----------
    parquet : str, optional
        Path to the input parquet file (quantms.io/qpx format). Mutually exclusive
        with ``msstats``.
    output : str
        Path for the output protein intensities file.
    sdrf : str, optional
        Path to SDRF file for sample metadata.
    msstats : str, optional
        Path to a quantms ``*_msstats_in.csv`` file. Requires ``sdrf``.
    quant_method : str
        Quantification method. Options:
        - 'directlfq': Uses DirectLFQ package (normalization + quantification)
        - 'ibaq': Intensity-Based Absolute Quantification
        - 'maxlfq': MaxLFQ algorithm
        - 'top3': Top 3 peptides per protein
        - 'top5': Top 5 peptides per protein
        - 'sum': Sum of all peptides
        - 'median': Median of peptides
    min_aa : int
        Minimum amino acid length for peptides. Default: 7.
    min_unique_peptides : int
        Minimum unique peptides per protein. Default: 2.
    remove_contaminants : bool
        Whether to remove contaminants and decoys. Default: True.
    run_normalization : str
        Run/technical replicate normalization method. Options:
        none, median, mean, max, global, max_min, IQR.
        Ignored when quant_method='directlfq'.
    sample_normalization : str
        Sample-to-sample normalization method. Options:
        - 'none': No normalization
        - 'globalMedian': Sample median / global median
        - 'conditionMedian': Condition-specific median
        - 'hierarchical': DirectLFQ-style hierarchical clustering
        Ignored when quant_method='directlfq'.
    normalization_proteins_file : str, optional
        File with protein IDs to use for normalization (one per line).
    fasta_file : str, optional
        FASTA file path. Required for iBAQ quantification.
    ion_alignment : str, optional
        Ion alignment method for MaxLFQ: none, hierarchical.
    ibaq_enzyme : str
        Protease used to digest the FASTA for the piBAQ denominator.
        Default: 'Trypsin'.
    ibaq_max_aa : int
        Maximum peptide length retained from the FASTA digest for piBAQ.
        Default: 50 (min length comes from ``min_aa``).
    ibaq_min_shared : int
        Minimum distinct peptides two proteins must share to co-cluster
        into one piBAQ family. Default: 2.
    ibaq_families_yaml : str, optional
        Path to a YAML file with explicit piBAQ family overrides. When
        ``None`` (default) family discovery is purely data-driven.
    ibaq_min_anchors : int
        Minimum proteotypic ("anchor") peptides a family member needs for
        the family to stay on the per-protein branch; the family rolls up to
        a single iBAQ only when no member reaches it. Default: 1.
    ibaq_high_anchor_threshold : int
        Minimum anchor count (weakest member) for a family to be labelled
        ``EvidenceLevel == "high"``. Default: 3.
    directlfq_num_cores : int, optional
        Number of cores for DirectLFQ parallel processing.
    directlfq_min_nonan : int
        Minimum non-missing values required by DirectLFQ for protein estimation.
    export_peptides : str, optional
        Path to export normalized peptides (for debugging/analysis).
    export_ions : str, optional
        Path to export normalized ions (DirectLFQ only).
    batch_correction : bool
        Whether to apply batch correction after quantification. Default: False.
    batch_method : str
        Batch detection method: sample_prefix, run, column. Default: sample_prefix.
    batch_column : str, optional
        Column name for explicit batch assignment (when batch_method='column').
    batch_covariates : list, optional
        SDRF columns to use as covariates (biological signal to preserve).
    batch_parametric : bool
        Use parametric estimation for ComBat. Default: True.
    batch_mean_only : bool
        Only adjust batch means, not individual effects. Default: False.
    batch_ref : int, optional
        Reference batch ID.

    Returns
    -------
    pd.DataFrame
        Protein intensities matrix.
    """
    config = PipelineConfig(
        input=InputConfig(
            parquet=parquet,
            sdrf=sdrf,
            fasta_file=fasta_file,
            msstats=msstats,
        ),
        filtering=FilterConfig(
            min_aa=min_aa,
            min_unique_peptides=min_unique_peptides,
            remove_contaminants=remove_contaminants,
        ),
        normalization=NormalizationConfig(
            run_method=run_normalization,
            sample_method=sample_normalization,
            proteins_file=normalization_proteins_file,
        ),
        quantification=QuantificationConfig(
            method=quant_method,
            ion_alignment=ion_alignment,
            coverage_threshold=coverage_threshold,
            ratio_fraction_merge=ratio_fraction_merge,
            directlfq_num_cores=directlfq_num_cores,
            directlfq_min_nonan=directlfq_min_nonan,
            ibaq_enzyme=ibaq_enzyme,
            ibaq_max_aa=ibaq_max_aa,
            ibaq_min_shared=ibaq_min_shared,
            ibaq_families_yaml=ibaq_families_yaml,
            ibaq_min_anchors=ibaq_min_anchors,
            ibaq_high_anchor_threshold=ibaq_high_anchor_threshold,
        ),
        irs=IRSConfig(
            enabled=irs,
            reference_samples=irs_reference_samples,
            sdrf_column=irs_sdrf_column,
            sdrf_values=irs_sdrf_values,
            reference_regex=irs_reference_regex,
            stat=irs_stat,
            remove_reference=irs_remove_reference,
        ),
        batch=BatchCorrectionConfig(
            enabled=batch_correction,
            method=batch_method,
            column=batch_column,
            covariates=batch_covariates,
            parametric=batch_parametric,
            mean_only=batch_mean_only,
            ref_batch=batch_ref,
        ),
        imputation=ImputationConfig(
            enabled=impute,
            method=impute_method,
            quantile=impute_quantile,
            shift=impute_shift,
            scale=impute_scale,
            n_neighbors=impute_n_neighbors,
        ),
        de=DEConfig(
            enabled=differential_expression,
            contrasts=de_contrasts,
            method=de_method,
            log2fc_threshold=de_log2fc_threshold,
            fdr_threshold=de_fdr_threshold,
            fdr_method=de_fdr_method,
            output=de_output,
            ensemble_methods=de_ensemble_methods,
            ensemble_min_k=de_ensemble_min_k,
        ),
        output=OutputConfig(
            export_peptides=export_peptides,
            export_ions=export_ions,
            plot_dir=plot_output_dir,
            plot_volcano=plot_volcano,
            plot_heatmap=plot_heatmap,
            plot_pca=plot_pca,
            highlight_genes=highlight_genes,
            interactive_report=interactive_report,
            report_output=report_output,
        ),
        runtime=RuntimeConfig(
            duckdb_memory=duckdb_memory,
            duckdb_threads=duckdb_threads,
        ),
    )

    pipeline = QuantificationPipeline(config)
    protein_df = pipeline.run()

    # Save output
    protein_df.to_csv(output, index=False)
    logger.info("Protein intensities saved to %s", output)

    return protein_df
