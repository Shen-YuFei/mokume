"""
CLI command for unified features to proteins pipeline.

This command provides a single-step workflow from feature-level parquet
files to protein intensities, supporting multiple quantification methods.
"""

import click

from mokume.model.normalization import (
    FeatureNormalizationMethod,
    PeptideNormalizationMethod,
)


# Build choices for sample normalization (including hierarchical)
SAMPLE_NORM_CHOICES = [p.name.lower() for p in PeptideNormalizationMethod]


@click.command(
    "features2proteins", short_help="Quantify proteins from feature parquet file."
)
@click.option(
    "-p",
    "--parquet",
    help="Parquet file (quantms.io/qpx format)",
    required=False,
    type=click.Path(exists=True),
)
@click.option(
    "--msstats",
    help="Legacy quantms *_msstats_in.csv file (requires --sdrf)",
    type=click.Path(exists=True),
)
@click.option(
    "-o",
    "--output",
    help="Output file for protein intensities",
    required=True,
    type=click.Path(),
)
@click.option(
    "-s",
    "--sdrf",
    help="SDRF file for sample metadata",
    default=None,
    type=click.Path(exists=True),
)
@click.option(
    "--quant-method",
    "quant_method",
    help=(
        "Quantification method: directlfq, ibaq, maxlfq, topn, sum, "
        "median, ratio, abd (TMT abundance), intensity (TMT reporter), "
        "spectral_count (PSM-based count)"
    ),
    type=click.Choice(
        [
            "directlfq",
            "ibaq",
            "maxlfq",
            "topn",
            "sum",
            "median",
            "ratio",
            "abd",
            "intensity",
            "spectral_count",
        ],
        case_sensitive=False,
    ),
    default="maxlfq",
    show_default=True,
)
@click.option(
    "--topn",
    "topn_peptides",
    help="Number of top peptides for TopN quantification (used when --quant-method=topn)",
    type=int,
    default=3,
    show_default=True,
)
# Filtering options
@click.option(
    "--min-aa",
    "min_aa",
    help="Minimum number of amino acids for peptides",
    type=int,
    default=7,
    show_default=True,
)
@click.option(
    "--min-unique",
    "min_unique",
    help="Minimum number of unique peptides per protein",
    type=int,
    default=2,
    show_default=True,
)
@click.option(
    "--remove-contaminants/--keep-contaminants",
    "remove_contaminants",
    help="Remove contaminants and decoys",
    default=True,
    show_default=True,
)
# Normalization options (ignored for directlfq)
@click.option(
    "--run-normalization",
    "run_normalization",
    help="Run/technical replicate normalization (ignored for directlfq)",
    type=click.Choice(
        [f.name.lower() for f in FeatureNormalizationMethod], case_sensitive=False
    ),
    default="median",
    show_default=True,
)
@click.option(
    "--sample-normalization",
    "sample_normalization",
    help="Sample normalization method (ignored for directlfq). "
    "Use 'hierarchical' for DirectLFQ-style clustering-based normalization.",
    type=click.Choice(SAMPLE_NORM_CHOICES, case_sensitive=False),
    default="globalmedian",
    show_default=True,
)
@click.option(
    "--normalization-proteins",
    "normalization_proteins",
    help="File with protein IDs to use for normalization (one per line)",
    type=click.Path(exists=True),
    default=None,
)
# Method-specific options
@click.option(
    "--fasta",
    "fasta_file",
    help="FASTA file (required for iBAQ)",
    type=click.Path(exists=True),
    default=None,
)
@click.option(
    "--ion-alignment",
    "ion_alignment",
    help="Ion alignment method for MaxLFQ",
    type=click.Choice(["none", "hierarchical"], case_sensitive=False),
    default=None,
)
# piBAQ-specific options (paralog-aware iBAQ)
@click.option(
    "--ibaq-enzyme",
    "ibaq_enzyme",
    help="Protease used to digest the FASTA for the piBAQ denominator (iBAQ only)",
    default="Trypsin",
    show_default=True,
)
@click.option(
    "--ibaq-max-aa",
    "ibaq_max_aa",
    help="Maximum peptide length from the FASTA digest for piBAQ (iBAQ only)",
    type=int,
    default=50,
    show_default=True,
)
@click.option(
    "--ibaq-min-shared",
    "ibaq_min_shared",
    help="Minimum distinct peptides two proteins must share to co-cluster "
    "into one piBAQ family (iBAQ only)",
    type=int,
    default=2,
    show_default=True,
)
@click.option(
    "--ibaq-families",
    "ibaq_families_yaml",
    help="YAML file with explicit piBAQ family overrides (iBAQ only)",
    type=click.Path(exists=True),
    default=None,
)
@click.option(
    "--ibaq-min-anchors",
    "ibaq_min_anchors",
    help="Minimum proteotypic anchors a family member needs to stay on the "
    "per-protein branch; the family rolls up only when no member reaches it "
    "(iBAQ only)",
    type=int,
    default=1,
    show_default=True,
)
@click.option(
    "--ibaq-high-anchor-threshold",
    "ibaq_high_anchor_threshold",
    help="Minimum anchor count (weakest member) for a family to be labelled "
    "EvidenceLevel='high' (iBAQ only)",
    type=int,
    default=3,
    show_default=True,
)
# DirectLFQ-specific options
@click.option(
    "--directlfq-cores",
    "directlfq_cores",
    help="Number of CPU cores for DirectLFQ",
    type=int,
    default=None,
)
@click.option(
    "--directlfq-min-nonan",
    "directlfq_min_nonan",
    help="Minimum non-NaN values for DirectLFQ",
    type=int,
    default=1,
    show_default=True,
)
# Optional exports
@click.option(
    "--export-peptides",
    "export_peptides",
    help="Export normalized peptides to this file (for debugging/analysis)",
    type=click.Path(),
    default=None,
)
@click.option(
    "--export-ions",
    "export_ions",
    help="Export normalized ions to this file (DirectLFQ only)",
    type=click.Path(),
    default=None,
)
@click.option(
    "--batch-correction",
    "batch_correction",
    help="Enable ComBat batch correction after quantification",
    is_flag=True,
    default=False,
)
@click.option(
    "--batch-method",
    "batch_method",
    help="Batch detection method",
    type=click.Choice(["sample_prefix", "run", "column"], case_sensitive=False),
    default="sample_prefix",
    show_default=True,
)
@click.option(
    "--batch-column",
    "batch_column",
    help="SDRF column to use when --batch-method=column",
    default=None,
)
@click.option(
    "--batch-covariates",
    "batch_covariates",
    help="Comma-separated SDRF columns to preserve as biological covariates",
    default=None,
)
@click.option(
    "--batch-parametric/--batch-nonparametric",
    "batch_parametric",
    help="Use parametric or non-parametric ComBat estimation",
    default=True,
    show_default=True,
)
@click.option(
    "--batch-mean-only",
    "batch_mean_only",
    help="Only adjust batch means, not individual effects",
    is_flag=True,
    default=False,
)
@click.option(
    "--batch-ref",
    "batch_ref",
    help="Reference batch ID for ComBat",
    type=int,
    default=None,
)
# IRS normalization options
@click.option(
    "--irs",
    "irs",
    help="Enable IRS (Internal Reference Scaling) normalization for multi-plex TMT data",
    is_flag=True,
    default=False,
)
@click.option(
    "--irs-reference-samples",
    "irs_reference_samples",
    help="Comma-separated list of reference sample names (source names)",
    default=None,
)
@click.option(
    "--irs-sdrf-column",
    "irs_sdrf_column",
    help="SDRF column to select reference samples from (e.g., 'factor value[disease]')",
    default=None,
)
@click.option(
    "--irs-sdrf-values",
    "irs_sdrf_values",
    help="Comma-separated values in SDRF column that indicate reference samples",
    default=None,
)
@click.option(
    "--irs-reference-regex",
    "irs_reference_regex",
    help="Regex to auto-detect reference samples across SDRF columns",
    default="pool|powder|ref|reference|bridge",
    show_default=True,
)
@click.option(
    "--irs-stat",
    "irs_stat",
    help="Statistic for computing reference intensity per plex",
    type=click.Choice(["median", "mean"], case_sensitive=False),
    default="median",
    show_default=True,
)
@click.option(
    "--irs-remove-reference",
    "irs_remove_reference",
    help="Remove reference samples from the final output",
    is_flag=True,
    default=False,
)
# Coverage filter
@click.option(
    "--coverage-threshold",
    "coverage_threshold",
    help="Minimum fraction of non-missing values per condition to keep a protein (e.g., 0.65)",
    type=float,
    default=None,
)
# Ratio quantification options
@click.option(
    "--ratio-fraction-merge",
    "ratio_fraction_merge",
    help="How to merge fractions in ratio quantification: mean (PS protocol) or max",
    type=click.Choice(["mean", "max"], case_sensitive=False),
    default="mean",
    show_default=True,
)
# Imputation options
@click.option(
    "--impute",
    "impute",
    help="Enable missing-value imputation on the protein matrix",
    is_flag=True,
    default=False,
)
@click.option(
    "--impute-method",
    "impute_method",
    help="Imputation method (operates in log2 space)",
    type=click.Choice(
        [
            "none",
            "knn",
            "minprob",
            "mindet",
            "qrilc",
            "missforest",
            "seqknn",
            "impseq",
        ],
        case_sensitive=False,
    ),
    default="none",
    show_default=True,
)
@click.option(
    "--impute-quantile",
    "impute_quantile",
    help="Quantile for MinProb/MinDet/QRILC low-tail draw",
    type=float,
    default=0.01,
    show_default=True,
)
@click.option(
    "--impute-shift",
    "impute_shift",
    help="MinProb shift in standard deviations",
    type=float,
    default=1.6,
    show_default=True,
)
@click.option(
    "--impute-scale",
    "impute_scale",
    help="MinProb scale factor for the imputation distribution sigma",
    type=float,
    default=0.3,
    show_default=True,
)
@click.option(
    "--impute-n-neighbors",
    "impute_n_neighbors",
    help="Number of neighbours for KNN/SeqKNN imputation",
    type=int,
    default=5,
    show_default=True,
)
# Differential expression options
@click.option(
    "--de",
    "differential_expression",
    help="Enable differential expression analysis",
    is_flag=True,
    default=False,
)
@click.option(
    "--de-contrasts",
    "de_contrasts",
    help="Comma-separated contrasts (e.g., 'A vs B,A vs C' or 'A-B,A-C')",
    default=None,
)
@click.option(
    "--de-contrasts-file",
    "de_contrasts_file",
    help="TSV file with contrasts (columns: group1, group2)",
    type=click.Path(exists=True),
    default=None,
)
@click.option(
    "--de-method",
    "de_method",
    help="DE statistical method (use 'ensemble' for top-k consensus across methods)",
    type=click.Choice(
        [
            "auto",
            "limrots",
            "limma",
            "deqms",
            "proda",
            "rots",
            "ensemble",
        ],
        case_sensitive=False,
    ),
    default="auto",
    show_default=True,
)
@click.option(
    "--de-ensemble-methods",
    "de_ensemble_methods",
    help="Comma-separated DE methods used by --de-method=ensemble (default: limrots,deqms,proda)",
    default=None,
)
@click.option(
    "--de-ensemble-min-k",
    "de_ensemble_min_k",
    help="Minimum number of ensemble members that must agree on direction",
    type=int,
    default=2,
    show_default=True,
)
@click.option(
    "--de-log2fc",
    "de_log2fc_threshold",
    help="Minimum absolute log2 fold change for significance",
    type=float,
    default=0.5,
    show_default=True,
)
@click.option(
    "--de-fdr",
    "de_fdr_threshold",
    help="Maximum adjusted p-value (FDR) for significance",
    type=float,
    default=0.05,
    show_default=True,
)
@click.option(
    "--de-fdr-method",
    "de_fdr_method",
    help="FDR correction method",
    type=click.Choice(["bh", "ihw"], case_sensitive=False),
    default="bh",
    show_default=True,
)
@click.option(
    "--de-output",
    "de_output",
    help="Output file for DE results",
    type=click.Path(),
    default=None,
)
# Plotting options
@click.option(
    "--plot-dir",
    "plot_output_dir",
    help="Output directory for plots",
    type=click.Path(),
    default=None,
)
@click.option(
    "--plot-volcano",
    "plot_volcano",
    help="Generate volcano plot from DE results",
    is_flag=True,
    default=False,
)
@click.option(
    "--plot-heatmap",
    "plot_heatmap",
    help="Generate heatmap of top variable proteins",
    is_flag=True,
    default=False,
)
@click.option(
    "--plot-pca",
    "plot_pca",
    help="Generate PCA plot colored by condition",
    is_flag=True,
    default=False,
)
@click.option(
    "--highlight-genes",
    "highlight_genes",
    help="Comma-separated gene names to highlight on volcano plot",
    default=None,
)
# Interactive report options
@click.option(
    "--interactive-report",
    "interactive_report",
    help="Generate interactive HTML report with plotly (requires mokume[reports])",
    is_flag=True,
    default=False,
)
@click.option(
    "--report-output",
    "report_output",
    help="Output path for interactive HTML report (default: <plot-dir>/report_<contrast>.html)",
    type=click.Path(),
    default=None,
)
# DuckDB resource limits (propagated through to mokume.io.feature.Feature)
@click.option(
    "--duckdb-memory",
    "duckdb_memory",
    help=(
        "DuckDB memory limit (e.g. '80GB', '16384MB'). Caps DuckDB's "
        "internal buffer pool only -- the surrounding Python process "
        "(PyArrow / polars / pandas) is not limited by this value and "
        "peak RSS can be 2-3x larger on wide pivots. Use cgroup MemoryMax, "
        "SLURM --mem, or container resources.limits for a hard cap. "
        "Default: DuckDB autoconfig (~80%% of total RAM)."
    ),
    default=None,
)
@click.option(
    "--duckdb-threads",
    "duckdb_threads",
    help="Number of threads DuckDB may use. Default: all available cores.",
    type=int,
    default=None,
)
@click.pass_context
def features2proteins(
    ctx,
    parquet: str,
    msstats: str,
    output: str,
    sdrf: str,
    quant_method: str,
    topn_peptides: int,
    min_aa: int,
    min_unique: int,
    remove_contaminants: bool,
    run_normalization: str,
    sample_normalization: str,
    normalization_proteins: str,
    fasta_file: str,
    ion_alignment: str,
    ibaq_enzyme: str,
    ibaq_max_aa: int,
    ibaq_min_shared: int,
    ibaq_families_yaml: str,
    ibaq_min_anchors: int,
    ibaq_high_anchor_threshold: int,
    directlfq_cores: int,
    directlfq_min_nonan: int,
    export_peptides: str,
    export_ions: str,
    # Batch correction
    batch_correction: bool,
    batch_method: str,
    batch_column: str,
    batch_covariates: str,
    batch_parametric: bool,
    batch_mean_only: bool,
    batch_ref: int,
    # IRS
    irs: bool,
    irs_reference_samples: str,
    irs_sdrf_column: str,
    irs_sdrf_values: str,
    irs_reference_regex: str,
    irs_stat: str,
    irs_remove_reference: bool,
    # Coverage filter
    coverage_threshold: float,
    # Ratio
    ratio_fraction_merge: str,
    # Imputation
    impute: bool,
    impute_method: str,
    impute_quantile: float,
    impute_shift: float,
    impute_scale: float,
    impute_n_neighbors: int,
    # DE
    differential_expression: bool,
    de_contrasts: str,
    de_contrasts_file: str,
    de_method: str,
    de_ensemble_methods: str,
    de_ensemble_min_k: int,
    de_log2fc_threshold: float,
    de_fdr_threshold: float,
    de_fdr_method: str,
    de_output: str,
    # Plots
    plot_output_dir: str,
    plot_volcano: bool,
    plot_heatmap: bool,
    plot_pca: bool,
    highlight_genes: str,
    # Interactive report
    interactive_report: bool,
    report_output: str,
    # DuckDB resource limits (cap DuckDB engine only, NOT total process RSS)
    duckdb_memory: str,
    duckdb_threads: int,
) -> None:
    """
    Quantify proteins from QPX feature parquet or legacy SDRF+MSstats input.

    This is the recommended unified command that handles the full pipeline
    from features to proteins in one step.

    \b
    QUANTIFICATION METHODS:
      directlfq  - DirectLFQ (uses directlfq package for everything)
      ibaq       - Intensity-Based Absolute Quantification (requires --fasta)
      maxlfq     - MaxLFQ algorithm
      topn       - Top N peptides per protein (use --topn to set N, default 3)
      sum        - Sum of all peptides
      median     - Median of peptides
      ratio      - PS protocol: log2(sample/reference) per plex (requires --sdrf)

    \b
    NORMALIZATION:
      When using directlfq, the DirectLFQ package handles all normalization.
      For other methods, mokume applies normalization:
      - Run normalization: normalizes technical replicates within samples
      - Sample normalization: normalizes samples relative to each other
        Use 'hierarchical' for DirectLFQ-style clustering normalization
        combined with other quantification methods (e.g., iBAQ).

    \b
    IRS NORMALIZATION (multi-plex TMT):
      Use --irs to enable Internal Reference Scaling. Reference channels
      are auto-detected from SDRF 'characteristics[pooled sample]' column,
      or specified via --irs-reference-samples, --irs-sdrf-column/values,
      or --irs-reference-regex.

    \b
    DIFFERENTIAL EXPRESSION:
      Use --de to enable DE analysis. Contrasts must be specified
      via --de-contrasts (e.g., 'NASH vs HL') or --de-contrasts-file (TSV).

    \b
    EXAMPLES:
      # TMT with IRS normalization + DE + volcano plot
      mokume features2proteins -p data.parquet -o proteins.csv -s sdrf.tsv \\
        --quant-method median --irs --irs-remove-reference \\
        --de --de-contrasts NASH-HL --de-output de_results.csv \\
        --plot-dir plots/ --plot-volcano --plot-pca

      # DirectLFQ quantification (uses directlfq package)
      mokume features2proteins -p data.parquet -o proteins.csv --quant-method directlfq

      # iBAQ with hierarchical normalization
      mokume features2proteins -p data.parquet -o proteins.csv \\
        --quant-method ibaq --sample-normalization hierarchical --fasta uniprot.fasta

      # Ratio quantification (PS protocol) with coverage filter + DE
      mokume features2proteins -p data.parquet -o proteins.csv -s sdrf.tsv \\
        --quant-method ratio --coverage-threshold 0.65 \\
        --de --de-method deqms --de-contrasts NASH-HL
    """
    from mokume.pipeline import features_to_proteins as run_pipeline

    if (parquet is None) == (msstats is None):
        raise click.UsageError("Provide exactly one of --parquet or --msstats")
    if msstats and not sdrf:
        raise click.UsageError("--msstats requires --sdrf")

    # Validate iBAQ requires fasta
    if quant_method.lower() == "ibaq" and not fasta_file:
        raise click.UsageError("iBAQ quantification requires --fasta option")

    # Validate ratio requires sdrf
    if quant_method.lower() == "ratio" and not sdrf:
        raise click.UsageError("Ratio quantification requires --sdrf option")

    if batch_correction and batch_method.lower() == "column" and not batch_column:
        raise click.UsageError(
            "Batch correction with method 'column' requires --batch-column option"
        )

    if batch_correction and (batch_column or batch_covariates) and not sdrf:
        raise click.UsageError(
            "Batch correction with --batch-column or --batch-covariates requires --sdrf option"
        )

    if not batch_correction and (
        batch_column
        or batch_covariates
        or batch_method != "sample_prefix"
        or not batch_parametric
        or batch_mean_only
        or batch_ref is not None
    ):
        click.echo(
            "Note: batch correction options are ignored unless --batch-correction is enabled.",
            err=True,
        )

    # Info about DirectLFQ ignoring normalization settings
    if quant_method.lower() == "directlfq":
        if run_normalization != "median" or sample_normalization != "globalmedian":
            click.echo(
                "Note: DirectLFQ handles its own normalization. "
                "--run-normalization and --sample-normalization are ignored.",
                err=True,
            )

    # Info about ratio ignoring normalization/IRS settings
    if quant_method.lower() == "ratio":
        click.echo(
            "Note: Ratio quantification handles cross-plex normalization via "
            "per-plex reference division. --run-normalization, "
            "--sample-normalization, and --irs are ignored.",
            err=True,
        )

    # Handle topn method - construct the method name with N
    effective_quant_method = quant_method
    if quant_method.lower() == "topn":
        effective_quant_method = f"top{topn_peptides}"
        click.echo(f"Using Top{topn_peptides} quantification method")

    # Parse comma-separated CLI values
    parsed_irs_ref_samples = (
        [s.strip() for s in irs_reference_samples.split(",")]
        if irs_reference_samples
        else None
    )
    parsed_irs_sdrf_values = (
        [s.strip() for s in irs_sdrf_values.split(",")] if irs_sdrf_values else None
    )
    parsed_batch_covariates = (
        [s.strip() for s in batch_covariates.split(",")] if batch_covariates else None
    )
    parsed_de_contrasts = (
        [s.strip() for s in de_contrasts.split(",")] if de_contrasts else []
    )
    if de_contrasts_file:
        import csv

        with open(de_contrasts_file, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            if not reader.fieldnames or not {"group1", "group2"}.issubset(
                reader.fieldnames
            ):
                raise click.UsageError(
                    f"Contrasts file '{de_contrasts_file}' must have 'group1' and 'group2' "
                    f"columns. Found: {reader.fieldnames}"
                )
            for row in reader:
                g1 = row.get("group1", "").strip()
                g2 = row.get("group2", "").strip()
                if g1 and g2:
                    parsed_de_contrasts.append(f"{g1} vs {g2}")
        click.echo(
            f"Loaded {len(parsed_de_contrasts)} contrasts "
            f"(inline + file: {de_contrasts_file})"
        )
    parsed_de_contrasts = parsed_de_contrasts or None
    parsed_highlight_genes = (
        [s.strip() for s in highlight_genes.split(",")] if highlight_genes else None
    )
    parsed_de_ensemble_methods = (
        [s.strip() for s in de_ensemble_methods.split(",")]
        if de_ensemble_methods is not None
        else None
    )

    # Run the pipeline
    run_pipeline(
        parquet=parquet,
        msstats=msstats,
        output=output,
        sdrf=sdrf,
        quant_method=effective_quant_method,
        min_aa=min_aa,
        min_unique_peptides=min_unique,
        remove_contaminants=remove_contaminants,
        run_normalization=run_normalization,
        sample_normalization=sample_normalization,
        normalization_proteins_file=normalization_proteins,
        fasta_file=fasta_file,
        ion_alignment=ion_alignment,
        ibaq_enzyme=ibaq_enzyme,
        ibaq_max_aa=ibaq_max_aa,
        ibaq_min_shared=ibaq_min_shared,
        ibaq_families_yaml=ibaq_families_yaml,
        ibaq_min_anchors=ibaq_min_anchors,
        ibaq_high_anchor_threshold=ibaq_high_anchor_threshold,
        directlfq_num_cores=directlfq_cores,
        directlfq_min_nonan=directlfq_min_nonan,
        export_peptides=export_peptides,
        export_ions=export_ions,
        # Batch correction
        batch_correction=batch_correction,
        batch_method=batch_method,
        batch_column=batch_column,
        batch_covariates=parsed_batch_covariates,
        batch_parametric=batch_parametric,
        batch_mean_only=batch_mean_only,
        batch_ref=batch_ref,
        # IRS
        irs=irs,
        irs_reference_samples=parsed_irs_ref_samples,
        irs_sdrf_column=irs_sdrf_column,
        irs_sdrf_values=parsed_irs_sdrf_values,
        irs_reference_regex=irs_reference_regex,
        irs_stat=irs_stat,
        irs_remove_reference=irs_remove_reference,
        # DE
        differential_expression=differential_expression,
        de_contrasts=parsed_de_contrasts,
        de_method=de_method,
        de_log2fc_threshold=de_log2fc_threshold,
        de_fdr_threshold=de_fdr_threshold,
        de_fdr_method=de_fdr_method,
        de_output=de_output,
        de_ensemble_methods=parsed_de_ensemble_methods,
        de_ensemble_min_k=de_ensemble_min_k,
        # Coverage filter
        coverage_threshold=coverage_threshold,
        # Ratio
        ratio_fraction_merge=ratio_fraction_merge,
        # Imputation
        impute=impute,
        impute_method=impute_method,
        impute_quantile=impute_quantile,
        impute_shift=impute_shift,
        impute_scale=impute_scale,
        impute_n_neighbors=impute_n_neighbors,
        # Plots
        plot_output_dir=plot_output_dir,
        plot_volcano=plot_volcano,
        plot_heatmap=plot_heatmap,
        plot_pca=plot_pca,
        highlight_genes=parsed_highlight_genes,
        # Interactive report
        interactive_report=interactive_report,
        report_output=report_output,
        # DuckDB resource limits (cap DuckDB engine only, NOT total process RSS)
        duckdb_memory=duckdb_memory,
        duckdb_threads=duckdb_threads,
    )

    click.echo(f"Protein intensities saved to: {output}")
