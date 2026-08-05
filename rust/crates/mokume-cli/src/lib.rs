use std::{
    ffi::OsString,
    fs::{create_dir_all, read_to_string, write, File},
    path::{Path, PathBuf},
    process::ExitCode,
    sync::Mutex,
};

use clap::{ArgAction, Args, Parser, Subcommand, ValueEnum};
use mokume_core::{
    AggregationLevel, BatchCorrectionConfig, DifferentialExpressionConfig, DirectLfqConfig,
    FeatureToPeptidesConfig, FeatureToProteinsConfig, FilterConfig, IbaqConfig, ImputationConfig,
    InputConfig, IrsChannelConfig, IrsConfig, IrsScope, IrsStat, MaxLfqConfig, MokumeError,
    NormalizationConfig, OutputConfig, OutputFormat, PreprocessingFilterConfig, QuantMethod,
    RatioConfig, RuntimeConfig,
};
use mokume_pipeline::{
    resolve_irs_autodetect_channel, run_features_to_peptides, run_features_to_proteins,
};
use tracing_subscriber::EnvFilter;

mod correct_batches;
mod h5ad;
mod peptides2protein;

#[derive(Debug, Parser)]
#[command(name = "mokume")]
#[command(version)]
#[command(about = "Rust-native proteomics quantification toolkit for quantms/QPX data")]
struct Cli {
    #[arg(
        short = 'v',
        long = "log-level",
        default_value = "debug",
        ignore_case = true
    )]
    log_level: LogLevel,

    #[arg(long = "log-file")]
    log_file: Option<PathBuf>,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum LogLevel {
    Debug,
    Info,
    Warn,
}

impl LogLevel {
    const fn as_filter(self) -> &'static str {
        match self {
            Self::Debug => "debug",
            Self::Info => "info",
            Self::Warn => "warn",
        }
    }
}

#[derive(Debug, Subcommand)]
enum Commands {
    #[command(name = "features2proteins")]
    #[command(about = "Quantify proteins from a QPX feature parquet file")]
    Features2Proteins(Box<Features2ProteinsArgs>),

    #[command(name = "features2peptides")]
    #[command(about = "Convert features to peptide-level output")]
    Features2Peptides(Box<Features2PeptidesArgs>),

    #[command(name = "peptides2protein")]
    #[command(about = "Compute protein quantities from peptide-level input")]
    Peptides2Protein(Box<Peptides2ProteinArgs>),

    #[command(name = "correct-batches")]
    #[command(about = "Correct batch effects in protein quantification output")]
    CorrectBatches(Box<CorrectBatchesArgs>),
}

#[allow(dead_code)]
#[derive(Debug, Args)]
struct Features2PeptidesArgs {
    #[arg(short = 'p', long = "parquet")]
    parquet: PathBuf,

    #[arg(short = 's', long = "sdrf")]
    sdrf: Option<PathBuf>,

    #[arg(long = "min_aa", visible_alias = "min-aa", default_value_t = 7)]
    min_aa: usize,

    #[arg(long = "min_unique", visible_alias = "min-unique", default_value_t = 2)]
    min_unique: usize,

    #[arg(long = "keep-shared-peptides")]
    keep_shared_peptides: bool,

    #[arg(long = "remove_ids", visible_alias = "remove-ids")]
    remove_ids: Option<PathBuf>,

    #[arg(
        long = "remove_decoy_contaminants",
        visible_alias = "remove-decoy-contaminants"
    )]
    remove_decoy_contaminants: bool,

    #[arg(
        long = "remove_low_frequency_peptides",
        visible_alias = "remove-low-frequency-peptides"
    )]
    remove_low_frequency_peptides: bool,

    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    #[arg(long = "skip_normalization", visible_alias = "skip-normalization")]
    skip_normalization: bool,

    #[arg(long = "run-normalization", default_value = "median", value_parser = [
        "none", "mean", "median", "max", "global", "max_min", "iqr",
    ], ignore_case = true)]
    run_normalization: String,

    #[arg(long = "sample-normalization", default_value = "globalmedian", value_parser = [
        "none",
        "globalmedian",
        "conditionmedian",
        "hierarchical",
        "quantile",
        "mediancenter",
        "meancenter",
        "rlr",
        "loess",
        "tmm",
    ], ignore_case = true)]
    sample_normalization: String,

    #[arg(long = "log2")]
    log2: bool,

    #[arg(long = "save_parquet", visible_alias = "save-parquet")]
    save_parquet: bool,

    #[arg(long = "irs_channel", visible_alias = "irs-channel")]
    irs_channel: Option<String>,

    #[arg(long = "irs_autodetect_regex", visible_alias = "irs-autodetect-regex")]
    irs_autodetect_regex: Option<String>,

    #[arg(long = "irs_stat", visible_alias = "irs-stat", default_value = "median", value_parser = ["median", "mean"], ignore_case = true)]
    irs_stat: String,

    #[arg(long = "irs_scope", visible_alias = "irs-scope", default_value = "global", value_parser = ["global", "by_mixture", "two_stage"], ignore_case = true)]
    irs_scope: String,

    #[arg(long = "aggregation_level", visible_alias = "aggregation-level", default_value = "sample", value_parser = ["sample", "run"], ignore_case = true)]
    aggregation_level: String,

    #[arg(long = "filter-config")]
    filter_config: Option<PathBuf>,

    #[arg(long = "generate-filter-config")]
    generate_filter_config: Option<PathBuf>,

    #[arg(long = "filter-min-intensity")]
    filter_min_intensity: Option<f64>,

    #[arg(long = "filter-cv-threshold")]
    filter_cv_threshold: Option<f64>,

    #[arg(long = "filter-charge-states")]
    filter_charge_states: Option<String>,

    #[arg(long = "filter-max-missed-cleavages")]
    filter_max_missed_cleavages: Option<usize>,

    #[arg(long = "filter-exclude-modifications")]
    filter_exclude_modifications: Option<String>,

    #[arg(long = "filter-min-unique-peptides")]
    filter_min_unique_peptides: Option<usize>,

    #[arg(long = "filter-min-features")]
    filter_min_features: Option<usize>,

    #[arg(long = "filter-max-missing-rate")]
    filter_max_missing_rate: Option<f64>,
}

#[allow(dead_code)]
#[derive(Debug, Args)]
struct Peptides2ProteinArgs {
    #[arg(short = 'f', long = "fasta")]
    fasta: Option<PathBuf>,

    #[arg(short = 'p', long = "peptides")]
    peptides: PathBuf,

    #[arg(long = "method", default_value = "ibaq", value_parser = [
        "ibaq", "top3", "topn", "maxlfq", "sum", "directlfq",
    ], ignore_case = true)]
    method: String,

    #[arg(short = 'e', long = "enzyme", default_value = "Trypsin")]
    enzyme: String,

    #[arg(short = 'n', long = "normalize")]
    normalize: bool,

    #[arg(long = "min_aa", visible_alias = "min-aa", default_value_t = 7)]
    min_aa: usize,

    #[arg(long = "max_aa", visible_alias = "max-aa", default_value_t = 30)]
    max_aa: usize,

    #[arg(short = 't', long = "tpa")]
    tpa: bool,

    #[arg(short = 'r', long = "ruler")]
    ruler: bool,

    #[arg(short = 'i', long = "ploidy", default_value_t = 2)]
    ploidy: i32,

    #[arg(short = 'm', long = "organism", default_value = "human")]
    organism: String,

    #[arg(short = 'c', long = "cpc", default_value_t = 200.0)]
    cpc: f64,

    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    #[arg(long = "verbose")]
    verbose: bool,

    #[arg(
        long = "qc_report",
        visible_alias = "qc-report",
        default_value = "QCprofile.pdf"
    )]
    qc_report: PathBuf,

    #[arg(long = "topn_n", visible_alias = "topn-n", default_value_t = 3)]
    topn_n: usize,

    #[arg(long = "threads", default_value_t = -1)]
    threads: i32,

    #[arg(long = "min_nonan", visible_alias = "min-nonan", default_value_t = 1)]
    min_nonan: usize,

    #[arg(long = "families")]
    families_yaml: Option<PathBuf>,

    #[arg(long = "min-shared", default_value_t = 2)]
    min_shared: usize,

    #[arg(long = "min-anchors", default_value_t = 1)]
    min_anchors: usize,

    #[arg(long = "high-anchor-threshold", default_value_t = 3)]
    high_anchor_threshold: usize,
}

#[derive(Debug, Args)]
struct CorrectBatchesArgs {
    #[arg(short = 'f', long = "folder")]
    folder: PathBuf,

    #[arg(short = 'p', long = "pattern", default_value = "*ibaq.tsv")]
    pattern: String,

    #[arg(long = "comment", default_value = "#")]
    comment: String,

    #[arg(long = "sep", default_value = "\t")]
    sep: String,

    #[arg(short = 'o', long = "output")]
    output: PathBuf,

    #[arg(
        long = "sample_id_column",
        visible_aliases = ["sample-id-column", "sid"],
        default_value = "SampleID"
    )]
    sample_id_column: String,

    #[arg(
        long = "protein_id_column",
        visible_aliases = ["protein-id-column", "pid"],
        default_value = "ProteinName"
    )]
    protein_id_column: String,

    #[arg(
        long = "ibaq_raw_column",
        visible_aliases = ["ibaq-raw-column", "ibaq"],
        default_value = "Ibaq"
    )]
    ibaq_raw_column: String,

    #[arg(
        long = "ibaq_corrected_column",
        visible_alias = "ibaq-corrected-column",
        default_value = "IbaqBec"
    )]
    ibaq_corrected_column: String,

    #[arg(long = "export_anndata", visible_alias = "export-anndata")]
    export_anndata: bool,
}

#[derive(Debug, Args)]
struct Features2ProteinsArgs {
    #[arg(short = 'p', long = "parquet")]
    parquet: PathBuf,

    #[arg(short = 'o', long = "output")]
    output: PathBuf,

    #[arg(long = "output-format", default_value = "python-compatible")]
    output_format: OutputFormatArg,

    #[arg(short = 's', long = "sdrf")]
    sdrf: Option<PathBuf>,

    #[arg(
        long = "quant-method",
        alias = "method",
        default_value = "maxlfq",
        ignore_case = true
    )]
    quant_method: QuantMethodArg,

    #[arg(long = "topn", default_value_t = 3)]
    topn_peptides: usize,

    #[arg(long = "min-aa", default_value_t = 7)]
    min_aa: usize,

    #[arg(long = "min-unique", default_value_t = 2)]
    min_unique: usize,

    #[arg(long = "remove-contaminants")]
    remove_contaminants: bool,

    #[arg(long = "keep-contaminants")]
    keep_contaminants: bool,

    #[arg(long = "run-normalization", default_value = "median", value_parser = [
        "none", "mean", "median", "max", "global", "max_min", "iqr",
    ], ignore_case = true)]
    run_normalization: String,

    #[arg(long = "sample-normalization", default_value = "globalmedian", value_parser = [
        "none",
        "globalmedian",
        "conditionmedian",
        "hierarchical",
        "quantile",
        "mediancenter",
        "meancenter",
        "rlr",
        "loess",
        "tmm",
    ], ignore_case = true)]
    sample_normalization: String,

    #[arg(long = "normalization-proteins")]
    normalization_proteins: Option<PathBuf>,

    #[arg(long = "fasta")]
    fasta: Option<PathBuf>,

    #[arg(long = "ion-alignment", value_parser = ["none", "hierarchical"], ignore_case = true)]
    ion_alignment: Option<String>,

    #[arg(long = "ibaq-enzyme", default_value = "Trypsin")]
    ibaq_enzyme: String,

    #[arg(long = "ibaq-max-aa", default_value_t = 50)]
    ibaq_max_aa: usize,

    #[arg(long = "ibaq-min-shared", default_value_t = 2)]
    ibaq_min_shared: usize,

    #[arg(long = "ibaq-families")]
    ibaq_families_yaml: Option<PathBuf>,

    #[arg(long = "ibaq-min-anchors", default_value_t = 1)]
    ibaq_min_anchors: usize,

    #[arg(long = "ibaq-high-anchor-threshold", default_value_t = 3)]
    ibaq_high_anchor_threshold: usize,

    #[arg(long = "directlfq-cores")]
    directlfq_cores: Option<usize>,

    #[arg(long = "directlfq-min-nonan", default_value_t = 1)]
    directlfq_min_nonan: usize,

    #[arg(long = "directlfq-num-samples-quadratic", default_value_t = 50)]
    directlfq_num_samples_quadratic: usize,

    #[arg(long = "export-peptides")]
    export_peptides: Option<PathBuf>,

    #[arg(long = "export-ions")]
    export_ions: Option<PathBuf>,

    #[arg(long = "batch-correction")]
    batch_correction: bool,

    #[arg(long = "batch-method", default_value = "sample_prefix", value_parser = ["sample_prefix", "run", "column"], ignore_case = true)]
    batch_method: String,

    #[arg(long = "batch-column")]
    batch_column: Option<String>,

    #[arg(long = "batch-covariates")]
    batch_covariates: Option<String>,

    #[arg(long = "batch-parametric", action = ArgAction::SetTrue, default_value_t = true)]
    batch_parametric: bool,

    #[arg(long = "batch-nonparametric")]
    batch_nonparametric: bool,

    #[arg(long = "batch-mean-only")]
    batch_mean_only: bool,

    #[arg(long = "batch-ref")]
    batch_ref: Option<i32>,

    #[arg(long = "irs")]
    irs: bool,

    #[arg(
        long = "irs-reference-samples",
        conflicts_with = "irs_reference_sample"
    )]
    irs_reference_samples: Option<String>,

    #[arg(
        long = "irs-reference-sample",
        conflicts_with = "irs_reference_samples"
    )]
    irs_reference_sample: Vec<String>,

    #[arg(long = "irs-sdrf-column")]
    irs_sdrf_column: Option<String>,

    #[arg(long = "irs-sdrf-values")]
    irs_sdrf_values: Option<String>,

    #[arg(
        long = "irs-reference-regex",
        default_value = "pool|powder|ref|reference|bridge"
    )]
    irs_reference_regex: String,

    #[arg(long = "irs-stat", default_value = "median", value_parser = ["median", "mean"], ignore_case = true)]
    irs_stat: String,

    #[arg(long = "irs-remove-reference")]
    irs_remove_reference: bool,

    #[arg(long = "coverage-threshold")]
    coverage_threshold: Option<f64>,

    #[arg(long = "ratio-fraction-merge", default_value = "mean", value_parser = ["mean", "max"], ignore_case = true)]
    ratio_fraction_merge: String,

    #[arg(long = "impute")]
    impute: bool,

    #[arg(long = "impute-method", default_value = "none", value_parser = [
        "none",
        "mean",
        "median",
        "constant",
        "zero",
        "most_frequent",
        "knn",
        "minprob",
        "mindet",
        "qrilc",
        "missforest",
        "seqknn",
        "impseq",
        "gms",
        "bpca",
        "impseqrob",
    ], ignore_case = true)]
    impute_method: String,

    #[arg(long = "impute-quantile", default_value_t = 0.01)]
    impute_quantile: f64,

    #[arg(long = "impute-shift", default_value_t = 1.6)]
    impute_shift: f64,

    #[arg(long = "impute-scale", default_value_t = 0.3)]
    impute_scale: f64,

    #[arg(long = "impute-n-neighbors", default_value_t = 5)]
    impute_n_neighbors: usize,

    #[arg(long = "de")]
    differential_expression: bool,

    #[arg(long = "de-contrasts")]
    de_contrasts: Option<String>,

    #[arg(long = "de-contrasts-file")]
    de_contrasts_file: Option<PathBuf>,

    #[arg(long = "de-method", default_value = "auto", value_parser = [
        "auto",
        "limrots",
        "limma",
        "deqms",
        "proda",
        "rots",
        "ensemble",
    ], ignore_case = true)]
    de_method: String,

    #[arg(long = "de-ensemble-methods")]
    de_ensemble_methods: Option<String>,

    #[arg(long = "de-ensemble-min-k", default_value_t = 2)]
    de_ensemble_min_k: usize,

    #[arg(long = "de-log2fc", default_value_t = 0.5)]
    de_log2fc_threshold: f64,

    #[arg(long = "de-fdr", default_value_t = 0.05)]
    de_fdr_threshold: f64,

    #[arg(long = "de-fdr-method", default_value = "bh", value_parser = ["bh", "ihw"], ignore_case = true)]
    de_fdr_method: String,

    #[arg(long = "de-output")]
    de_output: Option<PathBuf>,

    #[arg(long = "memory", visible_alias = "duckdb-memory")]
    memory: Option<String>,

    #[arg(long = "threads", visible_alias = "duckdb-threads")]
    threads: Option<usize>,
}

impl Features2ProteinsArgs {
    fn into_config(self) -> FeatureToProteinsConfig {
        let quantification = QuantMethod::from(self.quant_method);
        let remove_contaminants = if self.keep_contaminants {
            false
        } else {
            self.remove_contaminants || !self.keep_contaminants
        };

        let topn_peptides = self.topn_peptides;

        FeatureToProteinsConfig {
            input: InputConfig {
                parquet: self.parquet,
                sdrf: self.sdrf,
                fasta: self.fasta,
            },
            output: OutputConfig {
                protein_matrix: self.output,
                export_peptides: self.export_peptides,
                export_ions: self.export_ions,
                format: OutputFormat::from(self.output_format),
            },
            filtering: FilterConfig {
                min_aa: self.min_aa,
                min_unique_peptides: self.min_unique,
                remove_contaminants,
            },
            normalization: NormalizationConfig {
                run_method: self.run_normalization,
                sample_method: self.sample_normalization,
                normalization_proteins: self.normalization_proteins,
            },
            quantification,
            topn_peptides,
            maxlfq: MaxLfqConfig {
                ion_alignment: self.ion_alignment,
                force_builtin: false,
            },
            ibaq: IbaqConfig {
                enzyme: self.ibaq_enzyme,
                max_aa: self.ibaq_max_aa,
                min_shared: self.ibaq_min_shared,
                families_yaml: self.ibaq_families_yaml,
                min_anchors: self.ibaq_min_anchors,
                high_anchor_threshold: self.ibaq_high_anchor_threshold,
            },
            directlfq: DirectLfqConfig {
                cores: self.directlfq_cores,
                min_nonan: self.directlfq_min_nonan,
                num_samples_quadratic: self.directlfq_num_samples_quadratic,
            },
            batch: BatchCorrectionConfig {
                enabled: self.batch_correction,
                method: self.batch_method,
                column: self.batch_column,
                covariates: split_csv_option(self.batch_covariates),
                parametric: if self.batch_nonparametric {
                    false
                } else {
                    self.batch_parametric
                },
                mean_only: self.batch_mean_only,
                ref_batch: self.batch_ref,
            },
            irs: IrsConfig {
                enabled: self.irs,
                reference_samples: if self.irs_reference_sample.is_empty() {
                    split_csv_option(self.irs_reference_samples)
                } else {
                    Some(self.irs_reference_sample)
                },
                sdrf_column: self.irs_sdrf_column,
                sdrf_values: split_csv_option(self.irs_sdrf_values),
                reference_regex: self.irs_reference_regex,
                stat: self.irs_stat,
                remove_reference: self.irs_remove_reference,
            },
            coverage_threshold: self.coverage_threshold,
            ratio: RatioConfig {
                fraction_merge: self.ratio_fraction_merge,
            },
            imputation: ImputationConfig {
                enabled: self.impute,
                method: self.impute_method,
                quantile: self.impute_quantile,
                shift: self.impute_shift,
                scale: self.impute_scale,
                n_neighbors: self.impute_n_neighbors,
            },
            differential_expression: DifferentialExpressionConfig {
                enabled: self.differential_expression,
                contrasts: split_csv_option(self.de_contrasts),
                contrasts_file: self.de_contrasts_file,
                method: self.de_method,
                ensemble_methods: split_ensemble_methods(self.de_ensemble_methods),
                ensemble_min_k: self.de_ensemble_min_k,
                log2fc_threshold: self.de_log2fc_threshold,
                fdr_threshold: self.de_fdr_threshold,
                fdr_method: self.de_fdr_method,
                output: self.de_output,
            },
            runtime: RuntimeConfig {
                memory: self.memory,
                threads: self.threads,
            },
        }
    }
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum OutputFormatArg {
    #[value(name = "python-compatible")]
    PythonCompatible,
    #[value(name = "rust-native")]
    RustNative,
}

impl From<OutputFormatArg> for OutputFormat {
    fn from(value: OutputFormatArg) -> Self {
        match value {
            OutputFormatArg::PythonCompatible => Self::PythonCompatible,
            OutputFormatArg::RustNative => Self::RustNative,
        }
    }
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum QuantMethodArg {
    Directlfq,
    Ibaq,
    Maxlfq,
    Top3,
    Topn,
    Sum,
    Median,
    Ratio,
    Abd,
    Intensity,
    #[value(name = "spectral_count")]
    SpectralCount,
}

impl From<QuantMethodArg> for QuantMethod {
    fn from(value: QuantMethodArg) -> Self {
        match value {
            QuantMethodArg::Directlfq => Self::DirectLfq,
            QuantMethodArg::Ibaq => Self::Ibaq,
            QuantMethodArg::Maxlfq => Self::MaxLfq,
            QuantMethodArg::Top3 | QuantMethodArg::Topn => Self::TopN,
            QuantMethodArg::Sum => Self::Sum,
            QuantMethodArg::Median => Self::Median,
            QuantMethodArg::Ratio => Self::Ratio,
            QuantMethodArg::Abd => Self::Abd,
            QuantMethodArg::Intensity => Self::Intensity,
            QuantMethodArg::SpectralCount => Self::SpectralCount,
        }
    }
}

fn split_csv_option(value: Option<String>) -> Option<Vec<String>> {
    value
        .map(|value| {
            value
                .split(',')
                .map(str::trim)
                .filter(|part| !part.is_empty())
                .map(ToOwned::to_owned)
                .collect::<Vec<_>>()
        })
        .filter(|values| !values.is_empty())
}

/// Split ensemble members without discarding empty entries so validation can
/// report malformed lists instead of silently changing the requested methods.
fn split_ensemble_methods(value: Option<String>) -> Option<Vec<String>> {
    value.map(|value| {
        value
            .split(',')
            .map(str::trim)
            .map(ToOwned::to_owned)
            .collect()
    })
}

/// Dispatch a fully-built [`Cli`] to its subcommand. Shared by the binary entry
/// [`run`] and the library entry [`run_from_args`].
fn dispatch(cli: Cli) -> mokume_core::Result<()> {
    init_logging(cli.log_level, cli.log_file).and_then(|()| match cli.command {
        Commands::Features2Proteins(args) => dispatch_features_to_proteins(args.into_config()),
        Commands::Features2Peptides(args) => dispatch_features_to_peptides(&args),
        Commands::Peptides2Protein(args) => peptides2protein::run_peptides_to_protein(&args),
        Commands::CorrectBatches(args) => correct_batches::run_correct_batches(&args),
    })
}

/// Binary entry point: parse the process arguments (clap prints help / errors and
/// exits on its own) and turn the dispatch result into an [`ExitCode`].
pub fn run() -> ExitCode {
    match dispatch(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("{error}");
            ExitCode::FAILURE
        }
    }
}

/// Library entry point (used by the `mokume-py` PyO3 bindings): parse an explicit
/// argument vector and return the dispatch result. Unlike [`run`], a parse error
/// is returned as `Err` rather than exiting the process, so it never tears down a
/// hosting Python interpreter.
pub fn run_from_args<I, T>(args: I) -> mokume_core::Result<()>
where
    I: IntoIterator<Item = T>,
    T: Into<OsString> + Clone,
{
    let cli = Cli::try_parse_from(args).map_err(|error| MokumeError::InvalidInput {
        message: error.to_string(),
    })?;
    dispatch(cli)
}

/// Console-script entry point for the `mokume` wheel: parse an explicit argument
/// vector and return the process exit code WITHOUT calling `process::exit`, so it
/// never tears down a hosting Python interpreter. clap's help/version are printed
/// to stdout (exit 0) and usage errors to stderr (exit 2), exactly as the
/// standalone binary would; a dispatch failure prints the error and returns 1.
pub fn run_cli_from_args<I, T>(args: I) -> i32
where
    I: IntoIterator<Item = T>,
    T: Into<OsString> + Clone,
{
    match Cli::try_parse_from(args) {
        Ok(cli) => match dispatch(cli) {
            Ok(()) => 0,
            Err(error) => {
                eprintln!("{error}");
                1
            }
        },
        Err(error) => {
            // Reuse clap's own rendering + exit code: help/version go to stdout
            // with code 0, real parse errors to stderr with code 2.
            let _ = error.print();
            error.exit_code()
        }
    }
}

/// Run `features2proteins`: the Rust pipeline is pure-compute. It writes the
/// protein-matrix CSV and, when `--de-output` is set, one differential-expression
/// result CSV per contrast. Plotting, the interactive HTML report, and the other
/// visualization periphery now live in the Python wheel
/// (`python/mokume/commands/`) and are no longer invoked from here.
fn dispatch_features_to_proteins(config: FeatureToProteinsConfig) -> mokume_core::Result<()> {
    run_features_to_proteins(&config)
}

fn dispatch_features_to_peptides(args: &Features2PeptidesArgs) -> mokume_core::Result<()> {
    if !args.parquet.exists() {
        return Err(MokumeError::MissingInput {
            path: args.parquet.clone(),
        });
    }
    if let Some(sdrf) = &args.sdrf {
        if !sdrf.exists() {
            return Err(MokumeError::MissingInput { path: sdrf.clone() });
        }
    }
    // `--generate-filter-config` writes the example config and exits, mirroring
    // the Python command's early return.
    if let Some(path) = &args.generate_filter_config {
        generate_filter_config(path)?;
        return Ok(());
    }

    // `--remove_ids` reads a file of protein IDs to drop; validate it exists up
    // front (Python's `click.Path(exists=True)` does the same).
    if let Some(remove_ids) = &args.remove_ids {
        if !remove_ids.exists() {
            return Err(MokumeError::MissingInput {
                path: remove_ids.clone(),
            });
        }
    }

    // Channel IRS (Python `get_irs_scaling_factors`). All three scopes
    // (`global` / `by_mixture` / `two_stage`) are ported. The reference channel
    // is taken from `--irs_channel`; when that is absent but
    // `--irs_autodetect_regex` and `--sdrf` are given, it is autodetected from
    // the SDRF exactly as Python does (`peptide.py:219-233`). An autodetect that
    // finds no channel leaves IRS disabled, matching Python's warning path.
    // `--save_parquet` and `--aggregation_level run` are implemented below.
    let irs_channel = match (&args.irs_channel, &args.irs_autodetect_regex, &args.sdrf) {
        (Some(channel), _, _) => Some(channel.clone()),
        (None, Some(regex), Some(sdrf)) => resolve_irs_autodetect_channel(sdrf, regex)?,
        _ => None,
    };
    let irs = irs_channel.map(|channel| {
        let stat = if args.irs_stat.eq_ignore_ascii_case("mean") {
            IrsStat::Mean
        } else {
            IrsStat::Median
        };
        let scope = if args.irs_scope.eq_ignore_ascii_case("by_mixture") {
            IrsScope::ByMixture
        } else if args.irs_scope.eq_ignore_ascii_case("two_stage") {
            IrsScope::TwoStage
        } else {
            IrsScope::Global
        };
        IrsChannelConfig {
            channel,
            stat,
            scope,
        }
    });
    let filter_pipeline = build_filter_pipeline(args)?;

    let Some(output) = args.output.clone() else {
        return Err(MokumeError::InvalidInput {
            message: "features2peptides requires --output".to_owned(),
        });
    };

    let aggregation_level = if args.aggregation_level.eq_ignore_ascii_case("run") {
        AggregationLevel::Run
    } else {
        AggregationLevel::Sample
    };

    let config = FeatureToPeptidesConfig {
        input: InputConfig {
            parquet: args.parquet.clone(),
            sdrf: args.sdrf.clone(),
            fasta: None,
        },
        output,
        filtering: FilterConfig {
            min_aa: args.min_aa,
            min_unique_peptides: args.min_unique,
            remove_contaminants: args.remove_decoy_contaminants,
        },
        remove_ids: args.remove_ids.clone(),
        remove_low_frequency_peptides: args.remove_low_frequency_peptides,
        keep_shared_peptides: args.keep_shared_peptides,
        skip_normalization: args.skip_normalization,
        run_normalization: args.run_normalization.clone(),
        sample_normalization: args.sample_normalization.clone(),
        log2: args.log2,
        save_parquet: args.save_parquet,
        aggregation_level,
        filter_pipeline,
        irs,
    };
    run_features_to_peptides(&config)
}

/// Build the optional preprocessing filter pipeline from `--filter-config` and
/// the `--filter-*` overrides, mirroring Python's CLI assembly
/// (`commands/features2peptides.py`): a file (YAML or JSON, by extension) seeds
/// the config, otherwise any override starts from the defaults; then each CLI
/// override is applied. Returns `None` when neither a file nor any override is
/// given (default load-time filtering only).
fn build_filter_pipeline(
    args: &Features2PeptidesArgs,
) -> mokume_core::Result<Option<PreprocessingFilterConfig>> {
    let has_override = args.filter_min_intensity.is_some()
        || args.filter_cv_threshold.is_some()
        || args.filter_charge_states.is_some()
        || args.filter_max_missed_cleavages.is_some()
        || args.filter_exclude_modifications.is_some()
        || args.filter_min_unique_peptides.is_some()
        || args.filter_min_features.is_some()
        || args.filter_max_missing_rate.is_some();

    let mut config = match &args.filter_config {
        Some(path) => load_filter_config(path)?,
        None if has_override => PreprocessingFilterConfig {
            name: "cli_config".to_string(),
            ..PreprocessingFilterConfig::default()
        },
        None => return Ok(None),
    };

    // CLI overrides (Python `apply_overrides`): set only the provided fields.
    if let Some(value) = args.filter_min_intensity {
        config.intensity.min_intensity = value;
    }
    if let Some(value) = args.filter_cv_threshold {
        config.intensity.cv_threshold = Some(value);
    }
    if let Some(value) = &args.filter_charge_states {
        config.peptide.allowed_charge_states = Some(parse_int_csv(value)?);
    }
    if let Some(value) = args.filter_max_missed_cleavages {
        config.peptide.max_missed_cleavages = Some(value);
    }
    if let Some(value) = &args.filter_exclude_modifications {
        config.peptide.exclude_modifications = value
            .split(',')
            .map(|item| item.trim().to_string())
            .filter(|item| !item.is_empty())
            .collect();
    }
    if let Some(value) = args.filter_min_unique_peptides {
        config.protein.min_unique_peptides = value;
    }
    if let Some(value) = args.filter_min_features {
        config.run_qc.min_identified_features = value;
    }
    if let Some(value) = args.filter_max_missing_rate {
        config.run_qc.max_missing_rate = value;
    }

    Ok(Some(config))
}

/// Load a preprocessing filter config from a YAML or JSON file, choosing the
/// parser by extension (`.json` -> JSON, otherwise YAML), mirroring Python's
/// loader (`preprocessing/filters/io.py`).
fn load_filter_config(path: &Path) -> mokume_core::Result<PreprocessingFilterConfig> {
    if !path.exists() {
        return Err(MokumeError::MissingInput {
            path: path.to_path_buf(),
        });
    }
    let text = read_to_string(path).map_err(|source| MokumeError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let is_json = path
        .extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| extension.eq_ignore_ascii_case("json"));
    if is_json {
        serde_json::from_str(&text).map_err(|source| MokumeError::InvalidInput {
            message: format!(
                "failed to parse filter config '{}': {source}",
                path.display()
            ),
        })
    } else {
        serde_yaml::from_str(&text).map_err(|source| MokumeError::InvalidInput {
            message: format!(
                "failed to parse filter config '{}': {source}",
                path.display()
            ),
        })
    }
}

/// Parse a comma-separated integer list (e.g. `"2,3,4"`) into charge states,
/// mirroring Python's `[int(x.strip()) for x in value.split(",")]`.
fn parse_int_csv(value: &str) -> mokume_core::Result<Vec<i32>> {
    value
        .split(',')
        .map(str::trim)
        .filter(|item| !item.is_empty())
        .map(|item| {
            item.parse::<i32>().map_err(|_| MokumeError::InvalidInput {
                message: format!("invalid charge state '{item}': expected an integer"),
            })
        })
        .collect()
}

fn generate_filter_config(path: &Path) -> mokume_core::Result<()> {
    if path
        .extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| extension.eq_ignore_ascii_case("json"))
    {
        write(path, EXAMPLE_FILTER_CONFIG_JSON).map_err(|source| MokumeError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    } else {
        write(path, EXAMPLE_FILTER_CONFIG_YAML).map_err(|source| MokumeError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    }
    Ok(())
}

/// Initialize the global tracing subscriber at most once per process. The binary
/// calls this once; the library entry (`run_from_args`) may be invoked many times
/// from a hosting Python process, so the `Once` guard keeps repeat calls from
/// re-attempting initialization and emitting a spurious "already set" warning.
fn init_logging(level: LogLevel, log_file: Option<PathBuf>) -> mokume_core::Result<()> {
    static INIT: std::sync::Once = std::sync::Once::new();
    let mut outcome: mokume_core::Result<()> = Ok(());
    INIT.call_once(|| outcome = init_logging_once(level, log_file));
    outcome
}

fn init_logging_once(level: LogLevel, log_file: Option<PathBuf>) -> mokume_core::Result<()> {
    let filter = EnvFilter::new(level.as_filter());

    if let Some(path) = log_file {
        if let Some(parent) = path.parent() {
            if !parent.as_os_str().is_empty() {
                create_dir_all(parent).map_err(|source| MokumeError::Io {
                    path: parent.to_path_buf(),
                    source,
                })?;
            }
        }
        let file = File::create(&path).map_err(|source| MokumeError::Io {
            path: path.clone(),
            source,
        })?;
        if let Err(error) = tracing_subscriber::fmt()
            .with_env_filter(filter)
            .with_writer(Mutex::new(file))
            .try_init()
        {
            eprintln!("failed to initialize logging: {error}");
        }
    } else if let Err(error) = tracing_subscriber::fmt().with_env_filter(filter).try_init() {
        eprintln!("failed to initialize logging: {error}");
    }
    Ok(())
}

const EXAMPLE_FILTER_CONFIG_YAML: &str = r#"# Mokume Preprocessing Filter Configuration
# This file defines quality filters applied during peptide normalization

name: example_config

# Global options
enabled: true              # Set to false to disable all filtering
strict_mode: false         # If true, fail on any filter error
log_filtered_counts: true  # Log how many items each filter removes

# Intensity-based filters
intensity:
  min_intensity: 0.0           # Minimum intensity threshold (0 = no filter)
  cv_threshold: null           # Maximum CV across replicates (null = no filter)
  min_replicate_agreement: 1   # Min replicates where feature must be detected
  quantile_lower: 0.0          # Lower quantile for outlier removal (0-1)
  quantile_upper: 1.0          # Upper quantile for outlier removal (0-1)
  remove_zero_intensity: true  # Remove features with zero intensity

# Peptide-level filters
peptide:
  min_search_score: null          # Min search engine score (null = no filter)
  allowed_charge_states: null     # e.g., [2, 3, 4] or null for all charges
  exclude_modifications: []       # Modification names to exclude, e.g., ["Oxidation"]
  max_missed_cleavages: null      # Max missed cleavages (null = no filter)
  fdr_threshold: 0.01             # Peptide FDR threshold (requires q_value column)
  min_peptide_length: 7           # Minimum peptide length in amino acids
  max_peptide_length: 50          # Maximum peptide length in amino acids
  exclude_sequence_patterns: []   # Regex patterns to exclude
  require_unique_peptides: false  # Require peptides unique to one protein

# Protein-level filters
protein:
  fdr_threshold: 0.01         # Protein FDR threshold
  min_coverage: 0.0           # Minimum sequence coverage (0-1)
  min_peptides: 1             # Minimum total peptides per protein
  min_unique_peptides: 2      # Minimum unique peptides per protein
  razor_peptide_handling: keep   # How to handle shared peptides: keep, remove, assign_to_top
  protein_grouping: none         # Grouping strategy: none, subsumption, parsimony
  remove_contaminants: true      # Remove contaminant proteins
  remove_decoys: true            # Remove decoy proteins
  contaminant_patterns:          # Patterns identifying contaminants
    - CONTAMINANT
    - ENTRAP
    - DECOY

# Run/Sample QC filters
run_qc:
  min_total_intensity: 0.0      # Min total intensity per run
  min_identified_features: 0    # Min features per run
  min_identified_proteins: 0    # Min proteins per run
  min_sample_correlation: null  # Min correlation between samples (null = no filter)
  max_missing_rate: 1.0         # Max missing value rate (0-1)
"#;

const EXAMPLE_FILTER_CONFIG_JSON: &str = r#"{
  "name": "example_config",
  "intensity": {
    "min_intensity": 0.0,
    "cv_threshold": null,
    "min_replicate_agreement": 1,
    "quantile_lower": 0.0,
    "quantile_upper": 1.0,
    "remove_zero_intensity": true
  },
  "peptide": {
    "min_search_score": null,
    "allowed_charge_states": null,
    "exclude_modifications": [],
    "max_missed_cleavages": null,
    "fdr_threshold": 0.01,
    "min_peptide_length": 7,
    "max_peptide_length": 50,
    "exclude_sequence_patterns": [],
    "require_unique_peptides": false
  },
  "protein": {
    "fdr_threshold": 0.01,
    "min_coverage": 0.0,
    "min_peptides": 1,
    "min_unique_peptides": 2,
    "razor_peptide_handling": "keep",
    "protein_grouping": "none",
    "remove_contaminants": true,
    "remove_decoys": true,
    "contaminant_patterns": [
      "CONTAMINANT",
      "ENTRAP",
      "DECOY"
    ]
  },
  "run_qc": {
    "min_total_intensity": 0.0,
    "min_identified_features": 0,
    "min_identified_proteins": 0,
    "min_sample_correlation": null,
    "max_missing_rate": 1.0
  },
  "enabled": true,
  "strict_mode": false,
  "log_filtered_counts": true
}
"#;

#[cfg(test)]
mod tests {
    use clap::{CommandFactory, Parser};

    use super::{Cli, Commands};

    #[test]
    fn correct_batches_accepts_short_column_aliases() {
        // Python's click CLI accepts `-sid`/`-pid`/`-ibaq`; clap cannot express
        // single-dash multi-character shorts, so these are offered as `--sid` /
        // `--pid` / `--ibaq` long aliases (closest portable form).
        let cli = Cli::parse_from([
            "mokume",
            "correct-batches",
            "-f",
            "in",
            "-o",
            "out",
            "--sid",
            "MySample",
            "--pid",
            "MyProtein",
            "--ibaq",
            "MyIbaq",
        ]);
        let Commands::CorrectBatches(args) = cli.command else {
            panic!("expected the correct-batches subcommand");
        };
        assert_eq!(args.sample_id_column, "MySample");
        assert_eq!(args.protein_id_column, "MyProtein");
        assert_eq!(args.ibaq_raw_column, "MyIbaq");
    }

    #[test]
    fn parses_python_features2proteins_options() {
        let cli = Cli::parse_from([
            "mokume",
            "features2proteins",
            "-p",
            "input.parquet",
            "-o",
            "protein.csv",
            "-s",
            "input.sdrf.tsv",
            "--quant-method",
            "directlfq",
            "--directlfq-min-nonan",
            "3",
            "--batch-correction",
            "--batch-method",
            "column",
            "--batch-column",
            "characteristics[batch]",
            "--batch-covariates",
            "characteristics[sex],characteristics[organism part]",
            "--batch-nonparametric",
            "--batch-mean-only",
            "--batch-ref",
            "2",
            "--irs",
            "--irs-reference-samples",
            "Pool A,Pool B",
            "--coverage-threshold",
            "0.65",
            "--impute",
            "--impute-method",
            "knn",
            "--de",
            "--de-contrasts",
            "NASH-HL,NASH-Ctrl",
            "--de-fdr-method",
            "ihw",
            "--duckdb-memory",
            "80GB",
            "--duckdb-threads",
            "24",
        ]);

        let Commands::Features2Proteins(args) = cli.command else {
            panic!("expected features2proteins command");
        };
        let config = args.into_config();

        assert_eq!(config.directlfq.min_nonan, 3);
        assert!(config.batch.enabled);
        assert_eq!(config.batch.method, "column");
        assert_eq!(
            config.batch.column.as_deref(),
            Some("characteristics[batch]")
        );
        assert_eq!(
            config.batch.covariates.as_deref(),
            Some(
                &[
                    "characteristics[sex]".to_string(),
                    "characteristics[organism part]".to_string(),
                ][..]
            )
        );
        assert!(!config.batch.parametric);
        assert!(config.batch.mean_only);
        assert_eq!(config.batch.ref_batch, Some(2));
        assert_eq!(
            config.irs.reference_samples.as_deref(),
            Some(&["Pool A".to_string(), "Pool B".to_string()][..])
        );
        assert_eq!(config.coverage_threshold, Some(0.65));
        assert!(config.imputation.enabled);
        assert_eq!(config.imputation.method, "knn");
        assert!(config.differential_expression.enabled);
        assert_eq!(config.differential_expression.fdr_method, "ihw");
        assert_eq!(config.runtime.memory.as_deref(), Some("80GB"));
        assert_eq!(config.runtime.threads, Some(24));
    }

    #[test]
    fn preserves_empty_de_ensemble_members_for_validation() {
        for value in ["", "limma,,deqms"] {
            let cli = Cli::parse_from([
                "mokume",
                "features2proteins",
                "-p",
                "input.parquet",
                "-o",
                "protein.csv",
                "--de-ensemble-methods",
                value,
            ]);
            let Commands::Features2Proteins(args) = cli.command else {
                panic!("expected features2proteins command");
            };
            let methods = args.into_config().differential_expression.ensemble_methods;

            assert!(
                methods
                    .as_ref()
                    .is_some_and(|methods| methods.iter().any(String::is_empty)),
                "empty member from {value:?} was discarded: {methods:?}"
            );
        }
    }

    #[test]
    fn repeatable_reference_sample_preserves_commas() {
        let cli = Cli::parse_from([
            "mokume",
            "features2proteins",
            "-p",
            "input.parquet",
            "-o",
            "protein.csv",
            "--irs-reference-sample",
            "Pool, batch A",
            "--irs-reference-sample",
            "Pool B",
        ]);

        let Commands::Features2Proteins(args) = cli.command else {
            panic!("expected features2proteins command");
        };
        let config = args.into_config();

        assert_eq!(
            config.irs.reference_samples.as_deref(),
            Some(&["Pool, batch A".to_string(), "Pool B".to_string()][..])
        );
    }

    #[test]
    fn rejects_mixed_reference_sample_encodings() {
        let Err(error) = Cli::try_parse_from([
            "mokume",
            "features2proteins",
            "-p",
            "input.parquet",
            "-o",
            "protein.csv",
            "--irs-reference-samples",
            "Pool A,Pool B",
            "--irs-reference-sample",
            "Pool C",
        ]) else {
            panic!("mixed reference-sample encodings should conflict");
        };

        assert_eq!(error.kind(), clap::error::ErrorKind::ArgumentConflict);
    }

    #[test]
    fn top_level_help_lists_compute_commands() {
        let help = render_help(Cli::command());
        for command in [
            "features2proteins",
            "features2peptides",
            "peptides2protein",
            "correct-batches",
        ] {
            assert!(
                help.contains(command),
                "missing command `{command}` in help:\n{help}"
            );
        }
        // The visualization / tissue-map periphery moved to the Python wheel; the
        // Rust CLI is pure-compute and must not expose those subcommands.
        for command in ["tsne_visualization", "tissuemap", "agentic"] {
            assert!(
                !help.contains(command),
                "periphery command `{command}` must stay out of the Rust CLI:\n{help}"
            );
        }
    }

    #[test]
    fn features2proteins_help_lists_python_option_surface() {
        let help = render_subcommand_help("features2proteins");
        for option in [
            "--parquet",
            "--output",
            "--sdrf",
            "--quant-method",
            "--topn",
            "--min-aa",
            "--min-unique",
            "--remove-contaminants",
            "--keep-contaminants",
            "--run-normalization",
            "--sample-normalization",
            "--normalization-proteins",
            "--fasta",
            "--ion-alignment",
            "--ibaq-enzyme",
            "--ibaq-max-aa",
            "--ibaq-min-shared",
            "--ibaq-families",
            "--ibaq-min-anchors",
            "--ibaq-high-anchor-threshold",
            "--directlfq-cores",
            "--directlfq-min-nonan",
            "--directlfq-num-samples-quadratic",
            "--export-peptides",
            "--export-ions",
            "--batch-correction",
            "--batch-method",
            "--batch-column",
            "--batch-covariates",
            "--batch-parametric",
            "--batch-nonparametric",
            "--batch-mean-only",
            "--batch-ref",
            "--irs",
            "--irs-reference-samples",
            "--irs-reference-sample",
            "--irs-sdrf-column",
            "--irs-sdrf-values",
            "--irs-reference-regex",
            "--irs-stat",
            "--irs-remove-reference",
            "--coverage-threshold",
            "--ratio-fraction-merge",
            "--impute",
            "--impute-method",
            "--impute-quantile",
            "--impute-shift",
            "--impute-scale",
            "--impute-n-neighbors",
            "--de",
            "--de-contrasts",
            "--de-contrasts-file",
            "--de-method",
            "--de-ensemble-methods",
            "--de-ensemble-min-k",
            "--de-log2fc",
            "--de-fdr",
            "--de-fdr-method",
            "--de-output",
            "--duckdb-memory",
            "--duckdb-threads",
        ] {
            assert!(
                help.contains(option),
                "missing option `{option}` in help:\n{help}"
            );
        }
        // The plotting / interactive-report flags moved to the Python wheel; clap
        // must reject them now.
        for option in [
            "--plot-dir",
            "--plot-volcano",
            "--plot-heatmap",
            "--plot-pca",
            "--highlight-genes",
            "--interactive-report",
            "--report-output",
        ] {
            assert!(
                !help.contains(option),
                "removed plotting option `{option}` must not appear in help:\n{help}"
            );
        }
    }

    #[test]
    fn features2peptides_help_lists_python_option_surface() {
        let help = render_subcommand_help("features2peptides");
        for option in [
            "--parquet",
            "--sdrf",
            "--min_aa",
            "--min_unique",
            "--keep-shared-peptides",
            "--remove_ids",
            "--remove_decoy_contaminants",
            "--remove_low_frequency_peptides",
            "--output",
            "--skip_normalization",
            "--run-normalization",
            "--sample-normalization",
            "--log2",
            "--save_parquet",
            "--irs_channel",
            "--irs_autodetect_regex",
            "--irs_stat",
            "--irs_scope",
            "--aggregation_level",
            "--filter-config",
            "--generate-filter-config",
            "--filter-min-intensity",
            "--filter-cv-threshold",
            "--filter-charge-states",
            "--filter-max-missed-cleavages",
            "--filter-exclude-modifications",
            "--filter-min-unique-peptides",
            "--filter-min-features",
            "--filter-max-missing-rate",
        ] {
            assert!(
                help.contains(option),
                "missing option `{option}` in help:\n{help}"
            );
        }
    }

    #[test]
    fn peptides2protein_help_lists_python_option_surface() {
        let help = render_subcommand_help("peptides2protein");
        for option in [
            "--fasta",
            "--peptides",
            "--method",
            "--enzyme",
            "--normalize",
            "--min_aa",
            "--max_aa",
            "--tpa",
            "--ruler",
            "--ploidy",
            "--organism",
            "--cpc",
            "--output",
            "--verbose",
            "--qc_report",
            "--topn_n",
            "--threads",
            "--min_nonan",
            "--families",
            "--min-shared",
            "--min-anchors",
            "--high-anchor-threshold",
        ] {
            assert!(
                help.contains(option),
                "missing option `{option}` in help:\n{help}"
            );
        }
    }

    #[test]
    fn correct_batches_help_lists_python_option_surface() {
        let correct_batches_help = render_subcommand_help("correct-batches");
        for option in [
            "--folder",
            "--pattern",
            "--comment",
            "--sep",
            "--output",
            "--sample_id_column",
            "--protein_id_column",
            "--ibaq_raw_column",
            "--ibaq_corrected_column",
            "--export_anndata",
        ] {
            assert!(
                correct_batches_help.contains(option),
                "missing option `{option}` in help:\n{correct_batches_help}"
            );
        }
    }

    #[test]
    fn periphery_subcommands_are_removed() {
        // `tsne_visualization` and `tissuemap` moved to the Python wheel; the Rust
        // CLI no longer registers them as subcommands.
        let mut command = Cli::command();
        for name in ["tsne_visualization", "tsne-visualization", "tissuemap"] {
            assert!(
                command.find_subcommand_mut(name).is_none(),
                "periphery subcommand `{name}` must be removed from the Rust CLI"
            );
        }
    }

    fn render_help(mut command: clap::Command) -> String {
        let mut output = Vec::new();
        if let Err(error) = command.write_long_help(&mut output) {
            panic!("failed to render help: {error}");
        }
        match String::from_utf8(output) {
            Ok(help) => help,
            Err(error) => panic!("help is not valid utf8: {error}"),
        }
    }

    fn render_subcommand_help(name: &str) -> String {
        let mut command = Cli::command();
        let Some(subcommand) = command.find_subcommand_mut(name) else {
            panic!("missing subcommand `{name}`");
        };
        render_help(subcommand.clone())
    }
}
