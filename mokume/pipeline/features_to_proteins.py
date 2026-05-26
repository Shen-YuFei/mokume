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

When DirectLFQ is selected as the quantification method, the pipeline
delegates ALL processing (normalization + quantification) to the
directlfq package for reproducibility.
"""

import gc

import pandas as pd
from pathlib import Path
from typing import Optional

from mokume.core.logger import get_logger
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
)
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

    def _validate_config(self):
        """Validate configuration and check for required parameters."""
        if not Path(self.config.input.parquet).exists():
            raise FileNotFoundError(
                f"Parquet file not found: {self.config.input.parquet}"
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
        logger.info(f"Starting pipeline with quant_method={quant_method}")

        if quant_method == "directlfq":
            protein_df = self._run_directlfq_pipeline()
        elif quant_method == "ratio":
            protein_df = self._run_ratio_pipeline()
        else:
            protein_df = self._run_mokume_pipeline()

        # Apply IRS normalization if configured (skip for ratio — ratios already
        # handle cross-plex normalization via per-plex reference division)
        if self.config.irs.enabled and quant_method != "ratio":
            protein_df = self.normalization.apply_irs(protein_df)

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

    def _run_directlfq_pipeline(self) -> pd.DataFrame:
        """Run pipeline using DirectLFQ package."""
        try:
            import directlfq.protein_intensity_estimation as lfq_estimation
            import directlfq.normalization as lfq_norm
            import directlfq.config as lfq_config
        except ImportError:
            raise ImportError(
                "DirectLFQ quantification requires the directlfq package.\n"
                "Install with: pip install directlfq\n"
                "Or: pip install mokume[directlfq]"
            )

        logger.info("Loading and filtering data for DirectLFQ...")
        filtered_df = self.loading.load_for_directlfq()
        logger.info(f"Filtered data: {len(filtered_df)} features")

        logger.info("Converting to DirectLFQ format...")
        directlfq_input = self.loading.convert_to_directlfq_format(filtered_df)
        logger.info(f"DirectLFQ input shape: {directlfq_input.shape}")
        del filtered_df
        gc.collect()

        # Configure DirectLFQ
        lfq_config.set_global_protein_and_ion_id(protein_id="protein", quant_id="ion")
        lfq_config.set_compile_normalized_ion_table(
            self.config.output.export_ions is not None
        )

        # Run DirectLFQ normalization
        logger.info("Running DirectLFQ sample normalization...")
        normed_df = lfq_norm.NormalizationManagerSamplesOnSelectedProteins(
            directlfq_input,
            num_samples_quadratic=self.config.quantification.directlfq_num_samples_quadratic,
        ).complete_dataframe
        del directlfq_input
        gc.collect()

        # Run DirectLFQ protein estimation
        logger.info("Running DirectLFQ protein estimation...")
        protein_df, ion_df = lfq_estimation.estimate_protein_intensities(
            normed_df,
            min_nonan=self.config.quantification.directlfq_min_nonan,
            num_samples_quadratic=10,
            num_cores=self.config.quantification.directlfq_num_cores,
        )

        # Export ions if requested
        if self.config.output.export_ions and ion_df is not None:
            logger.info(f"Exporting ions to {self.config.output.export_ions}")
            ion_df.to_csv(self.config.output.export_ions)

        logger.info(f"DirectLFQ complete: {len(protein_df)} proteins")
        return protein_df

    def _run_mokume_pipeline(self) -> pd.DataFrame:
        """Run pipeline using mokume's native implementations."""
        logger.info("Loading and filtering data...")
        peptide_df = self.loading.load_for_mokume()
        logger.info(f"Processed peptides: {len(peptide_df)} rows")

        # Export peptides if requested
        if self.config.output.export_peptides:
            logger.info(f"Exporting peptides to {self.config.output.export_peptides}")
            peptide_df.to_csv(self.config.output.export_peptides, index=False)

        # Quantify proteins
        logger.info(
            f"Quantifying proteins with method: {self.config.quantification.method}"
        )
        protein_df = self.quantification.quantify(peptide_df)
        logger.info(f"Quantification complete: {len(protein_df)} proteins")

        return protein_df

    def _run_ratio_pipeline(self) -> pd.DataFrame:
        """Run ratio-based quantification (PS protocol)."""
        from mokume.quantification.ratio import RatioQuantification

        logger.info("Running ratio-based quantification (PS protocol)...")

        psm_df, ref_samples, sample_to_plex = self.loading.load_for_ratio()

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

        logger.info(f"Ratio pipeline complete: {len(protein_df)} proteins")
        return protein_df


def features_to_proteins(
    parquet: str,
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
) -> pd.DataFrame:
    """
    Quantify proteins directly from feature parquet file.

    This is the main entry point for the unified pipeline that handles
    the full workflow from features to proteins in one step.

    Parameters
    ----------
    parquet : str
        Path to the input parquet file (quantms.io/qpx format).
    output : str
        Path for the output protein intensities file.
    sdrf : str, optional
        Path to SDRF file for sample metadata.
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
    )

    pipeline = QuantificationPipeline(config)
    protein_df = pipeline.run()

    # Save output
    protein_df.to_csv(output, index=False)
    logger.info(f"Protein intensities saved to {output}")

    return protein_df
