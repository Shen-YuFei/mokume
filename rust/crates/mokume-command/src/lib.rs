use std::{
    ffi::OsString,
    fs::{create_dir_all, read_to_string, write, File},
    path::{Path, PathBuf},
    str::FromStr,
    sync::Mutex,
};

use clap::{Args, Parser, Subcommand, ValueEnum};
use mokume_core::quant::parse_topn_from_method_name;
use mokume_core::{
    parse_memory_to_bytes, AggregationLevel, BatchCorrectionConfig, DifferentialExpressionConfig,
    DirectLfqConfig, FeatureToPeptidesConfig, FeatureToProteinsConfig, FilterConfig,
    ImputationConfig, InputConfig, IrsChannelConfig, IrsConfig, IrsScope, IrsStat, MaxLfqConfig,
    MokumeError, NamedScoreFilterConfig, NormalizationConfig, OutputConfig, OutputFormat,
    PibaqConfig, PreprocessingFilterConfig, QuantMethod, RatioConfig, RuntimeConfig,
};
use mokume_pipeline::{
    resolve_irs_autodetect_channel, run_features_to_peptides, run_features_to_proteins,
    run_features_to_proteins_with_pibaq_digest, PibaqDigest,
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
        ignore_case = true,
        global = true
    )]
    log_level: LogLevel,

    #[arg(long = "log-file", global = true)]
    log_file: Option<PathBuf>,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
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

/// The FASTA digestion requested by a parsed piBAQ command.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PibaqDigestRequest {
    pub fasta: PathBuf,
    pub enzyme: String,
    pub min_aa: usize,
    pub max_aa: usize,
    pub missed_cleavages: usize,
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
    #[arg(
        short = 'p',
        long = "parquet",
        required_unless_present = "generate_filter_config"
    )]
    parquet: Option<PathBuf>,

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

    #[arg(
        short = 'o',
        long = "output",
        required_unless_present = "generate_filter_config"
    )]
    output: Option<PathBuf>,

    #[arg(
        long = "skip_normalization",
        visible_alias = "skip-normalization",
        conflicts_with_all = ["run_normalization", "sample_normalization"]
    )]
    skip_normalization: bool,

    #[arg(long = "run-normalization", value_parser = [
        "none", "mean", "median", "max", "global", "max_min", "iqr",
    ], ignore_case = true)]
    run_normalization: Option<String>,

    #[arg(long = "sample-normalization", value_parser = [
        "none",
        "globalmedian",
        "conditionmedian",
    ], ignore_case = true)]
    sample_normalization: Option<String>,

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

    #[arg(long = "filter-peptide-fdr", value_parser = parse_fraction)]
    filter_peptide_fdr: Option<f64>,

    #[arg(
        long = "filter-score",
        value_name = "NAME=THRESHOLD",
        value_parser = parse_named_score_filter,
        help = "Filter one exact QPX additional score; comparison direction comes from higher_better"
    )]
    filter_score: Option<NamedScoreFilterConfig>,

    #[arg(long = "filter-exclude-modifications")]
    filter_exclude_modifications: Option<String>,

    #[arg(long = "filter-min-unique-peptides")]
    filter_min_unique_peptides: Option<usize>,

    #[arg(long = "filter-protein-fdr", value_parser = parse_fraction)]
    filter_protein_fdr: Option<f64>,

    #[arg(long = "filter-min-features")]
    filter_min_features: Option<usize>,

    #[arg(long = "filter-max-missing-rate", value_parser = parse_fraction)]
    filter_max_missing_rate: Option<f64>,
}

#[allow(dead_code)]
#[derive(Debug, Args)]
struct Peptides2ProteinArgs {
    #[arg(short = 'f', long = "fasta")]
    fasta: Option<PathBuf>,

    #[arg(short = 'p', long = "peptides")]
    peptides: PathBuf,

    #[arg(
        long = "method",
        default_value = "pibaq",
        value_name = "[pibaq|maxlfq|sum|directlfq|top<N>]",
        value_parser = parse_peptides2protein_method,
        help = "Quantification method: pibaq, maxlfq, sum, directlfq, or top<N> -- the TopN \
family spells its peptide count in the name (e.g. top3, top5)"
    )]
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

    #[arg(short = 'i', long = "ploidy", value_parser = parse_positive_i32)]
    ploidy: Option<i32>,

    #[arg(short = 'm', long = "organism")]
    organism: Option<String>,

    #[arg(short = 'c', long = "cpc", value_parser = parse_positive_f64)]
    cpc: Option<f64>,

    #[arg(short = 'o', long = "output", required = true)]
    output: Option<PathBuf>,

    #[arg(long = "verbose")]
    verbose: bool,

    #[arg(
        long = "qc_report",
        visible_alias = "qc-report",
        default_value = "QCprofile.pdf"
    )]
    qc_report: PathBuf,

    #[arg(
        long = "threads",
        default_value_t = -1,
        value_parser = parse_nonzero_threads,
        help = "DirectLFQ/MaxLFQ worker count; negative values use joblib CPU-relative semantics"
    )]
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

    #[arg(short = 'p', long = "pattern", default_value = "*pibaq.tsv")]
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
        long = "pibaq_raw_column",
        visible_aliases = [
            "pibaq-raw-column",
            "pibaq"
        ],
        default_value = "PiBAQ"
    )]
    pibaq_raw_column: String,

    #[arg(
        long = "pibaq_corrected_column",
        visible_alias = "pibaq-corrected-column",
        default_value = "PiBAQBec"
    )]
    pibaq_corrected_column: String,

    #[arg(long = "export_anndata", visible_alias = "export-anndata")]
    export_anndata: bool,
}

#[derive(Debug, Args)]
struct Features2ProteinsArgs {
    #[arg(
        short = 'p',
        long = "parquet",
        required_unless_present = "msstats",
        conflicts_with = "msstats"
    )]
    parquet: Option<PathBuf>,

    #[arg(
        long = "msstats",
        required_unless_present = "parquet",
        conflicts_with = "parquet",
        requires = "sdrf"
    )]
    msstats: Option<PathBuf>,

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
        value_name = "[directlfq|pibaq|maxlfq|sum|median|ratio|abd|intensity|spectral_count|top<N>]",
        value_parser = parse_quant_method,
        help = "Quantification method: directlfq, pibaq, maxlfq, sum, median, ratio, abd, \
intensity, spectral_count, or top<N> -- the TopN family spells its peptide count in the name \
(e.g. top3, top5)"
    )]
    quant_method: QuantMethodArg,

    #[arg(long = "min-aa", default_value_t = 7)]
    min_aa: usize,

    #[arg(long = "min-unique", default_value_t = 2)]
    min_unique: usize,

    #[arg(long = "keep-contaminants")]
    keep_contaminants: bool,

    #[arg(long = "run-normalization", value_parser = [
        "none", "mean", "median", "max", "global", "max_min", "iqr",
    ], ignore_case = true)]
    run_normalization: Option<String>,

    #[arg(long = "sample-normalization", value_parser = [
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
    sample_normalization: Option<String>,

    #[arg(long = "normalization-proteins")]
    normalization_proteins: Option<PathBuf>,

    #[arg(long = "fasta")]
    fasta: Option<PathBuf>,

    #[arg(long = "pibaq-enzyme", default_value = "Trypsin")]
    pibaq_enzyme: String,

    #[arg(long = "pibaq-max-aa", default_value_t = 50)]
    pibaq_max_aa: usize,

    #[arg(long = "pibaq-min-shared", default_value_t = 2)]
    pibaq_min_shared: usize,

    #[arg(long = "pibaq-families")]
    pibaq_families_yaml: Option<PathBuf>,

    #[arg(long = "pibaq-min-anchors", default_value_t = 1)]
    pibaq_min_anchors: usize,

    #[arg(long = "pibaq-high-anchor-threshold", default_value_t = 3)]
    pibaq_high_anchor_threshold: usize,

    #[arg(long = "directlfq-cores", value_parser = parse_positive_usize)]
    directlfq_cores: Option<usize>,

    #[arg(long = "directlfq-min-nonan", value_parser = parse_positive_usize)]
    directlfq_min_nonan: Option<usize>,

    #[arg(long = "directlfq-num-samples-quadratic", value_parser = parse_positive_usize)]
    directlfq_num_samples_quadratic: Option<usize>,

    #[arg(long = "export-peptides")]
    export_peptides: Option<PathBuf>,

    #[arg(long = "export-ions")]
    export_ions: Option<PathBuf>,

    #[arg(long = "batch-correction")]
    batch_correction: bool,

    #[arg(long = "batch-method", value_parser = ["sample_prefix", "column"], ignore_case = true)]
    batch_method: Option<String>,

    #[arg(long = "batch-column")]
    batch_column: Option<String>,

    #[arg(long = "batch-covariates")]
    batch_covariates: Option<String>,

    #[arg(long = "batch-nonparametric")]
    batch_nonparametric: bool,

    #[arg(long = "batch-mean-only")]
    batch_mean_only: bool,

    #[arg(long = "batch-ref", value_parser = parse_nonnegative_i32)]
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

    #[arg(long = "irs-reference-regex")]
    irs_reference_regex: Option<String>,

    #[arg(long = "irs-stat", value_parser = ["median", "mean"], ignore_case = true)]
    irs_stat: Option<String>,

    #[arg(long = "irs-remove-reference")]
    irs_remove_reference: bool,

    #[arg(long = "coverage-threshold", value_parser = parse_fraction)]
    coverage_threshold: Option<f64>,

    #[arg(
        long = "min-sample-correlation",
        value_parser = parse_correlation,
        requires = "sdrf",
        help = "Drop samples whose mean Pearson correlation to same-condition peers is below this value"
    )]
    min_sample_correlation: Option<f64>,

    #[arg(long = "ratio-fraction-merge", value_parser = ["mean", "max"], ignore_case = true)]
    ratio_fraction_merge: Option<String>,

    #[arg(long = "impute")]
    impute: bool,

    #[arg(long = "impute-method", value_parser = [
        "mean",
        "median",
        "constant",
        "zero",
        "most_frequent",
        "knn",
        "minprob",
        "mindet",
        "qrilc",
        "seqknn",
        "impseq",
        "gms",
        "bpca",
        "impseqrob",
    ], ignore_case = true)]
    impute_method: Option<String>,

    #[arg(long = "impute-quantile", value_parser = parse_fraction)]
    impute_quantile: Option<f64>,

    #[arg(long = "impute-shift", value_parser = parse_finite_f64)]
    impute_shift: Option<f64>,

    #[arg(long = "impute-scale", value_parser = parse_nonnegative_f64)]
    impute_scale: Option<f64>,

    #[arg(long = "impute-n-neighbors", value_parser = parse_positive_usize)]
    impute_n_neighbors: Option<usize>,

    #[arg(long = "de")]
    differential_expression: bool,

    #[arg(long = "de-contrasts")]
    de_contrasts: Option<String>,

    #[arg(long = "de-contrasts-file")]
    de_contrasts_file: Option<PathBuf>,

    #[arg(long = "de-method", value_parser = [
        "auto",
        "limrots",
        "limma",
        "deqms",
        "proda",
        "rots",
        "ensemble",
    ], ignore_case = true)]
    de_method: Option<String>,

    #[arg(long = "de-ensemble-methods")]
    de_ensemble_methods: Option<String>,

    #[arg(
        long = "de-ensemble-min-k",
        value_parser = parse_positive_usize,
        help = "Minimum agreeing ensemble members (default: 2; ensemble only)"
    )]
    de_ensemble_min_k: Option<usize>,

    #[arg(long = "de-log2fc", value_parser = parse_de_log2fc)]
    de_log2fc_threshold: Option<DeLog2FcArg>,

    #[arg(long = "de-effect-size-gate", value_parser = ["mixture", "null_quantile"], ignore_case = true)]
    de_effect_size_gate: Option<String>,

    #[arg(long = "de-fdr", value_parser = parse_fraction)]
    de_fdr_threshold: Option<f64>,

    #[arg(long = "de-fdr-method", value_parser = ["bh", "ihw", "bky", "storey"], ignore_case = true)]
    de_fdr_method: Option<String>,

    #[arg(long = "de-output")]
    de_output: Option<PathBuf>,

    #[arg(
        long = "memory",
        value_parser = parse_memory,
        help = "Linux-only soft process RSS budget (for example 1GB or 512MB); also reduces QPX batch/read-ahead memory"
    )]
    memory: Option<String>,

    #[arg(
        long = "threads",
        visible_alias = "duckdb-threads",
        value_parser = parse_positive_usize
    )]
    threads: Option<usize>,
}

impl Features2ProteinsArgs {
    fn into_config(self) -> mokume_core::Result<FeatureToProteinsConfig> {
        // `top<N>` is the only spelling of a TopN method, so N comes from the
        // method name; the pipeline reads it from `topn_peptides`.
        let QuantMethodArg {
            method: quantification,
            topn,
        } = self.quant_method;
        let topn_peptides = topn.unwrap_or(DEFAULT_TOPN_PEPTIDES);
        let remove_contaminants = !self.keep_contaminants;

        if self.threads.is_some() && self.directlfq_cores.is_some() {
            return Err(MokumeError::InvalidInput {
                message: "choose either --threads or --directlfq-cores, not both".to_owned(),
            });
        }
        if quantification != QuantMethod::DirectLfq
            && (self.directlfq_cores.is_some() || self.directlfq_min_nonan.is_some())
        {
            return Err(MokumeError::InvalidInput {
                message: "--directlfq-cores/--directlfq-min-nonan require --quant-method directlfq"
                    .to_owned(),
            });
        }
        if !matches!(quantification, QuantMethod::DirectLfq | QuantMethod::MaxLfq)
            && self.directlfq_num_samples_quadratic.is_some()
        {
            return Err(MokumeError::InvalidInput {
                message: "--directlfq-num-samples-quadratic only applies to DirectLFQ/MaxLFQ"
                    .to_owned(),
            });
        }

        let quantification_manages_normalization =
            matches!(quantification, QuantMethod::DirectLfq | QuantMethod::Ratio);
        let run_normalization = self.run_normalization.unwrap_or_else(|| {
            if quantification_manages_normalization {
                "none".to_owned()
            } else {
                "median".to_owned()
            }
        });
        let sample_normalization = self.sample_normalization.unwrap_or_else(|| {
            if quantification_manages_normalization {
                "none".to_owned()
            } else {
                "globalmedian".to_owned()
            }
        });

        let batch_method_supplied = self.batch_method.is_some();
        let batch_method = self
            .batch_method
            .unwrap_or_else(|| "sample_prefix".to_owned());
        if !self.batch_correction
            && (batch_method_supplied
                || self.batch_column.is_some()
                || self.batch_covariates.is_some()
                || self.batch_nonparametric
                || self.batch_mean_only
                || self.batch_ref.is_some())
        {
            return Err(MokumeError::InvalidInput {
                message: "batch options require --batch-correction".to_owned(),
            });
        }

        let reference_samples = if self.irs_reference_sample.is_empty() {
            split_csv_option(self.irs_reference_samples)
        } else {
            Some(self.irs_reference_sample)
        };
        if self.irs_sdrf_column.is_some() != self.irs_sdrf_values.is_some() {
            return Err(MokumeError::InvalidInput {
                message: "--irs-sdrf-column and --irs-sdrf-values must be provided together"
                    .to_owned(),
            });
        }
        let selector_count = usize::from(reference_samples.is_some())
            + usize::from(self.irs_sdrf_column.is_some())
            + usize::from(self.irs_reference_regex.is_some());
        if selector_count > 1 {
            return Err(MokumeError::InvalidInput {
                message: "choose one reference selector: samples, SDRF column+values, or regex"
                    .to_owned(),
            });
        }
        if quantification == QuantMethod::Ratio {
            if self.irs {
                return Err(MokumeError::InvalidInput {
                    message: "Ratio quantification cannot also apply IRS".to_owned(),
                });
            }
            if self.irs_sdrf_column.is_some()
                || self.irs_sdrf_values.is_some()
                || self.irs_stat.is_some()
                || self.irs_remove_reference
            {
                return Err(MokumeError::InvalidInput {
                    message: "Ratio accepts --irs-reference-samples or --irs-reference-regex; IRS-only options require --irs"
                        .to_owned(),
                });
            }
        } else if !self.irs
            && (selector_count > 0 || self.irs_stat.is_some() || self.irs_remove_reference)
        {
            return Err(MokumeError::InvalidInput {
                message: "IRS options require --irs".to_owned(),
            });
        }
        let irs_reference_regex = self
            .irs_reference_regex
            .unwrap_or_else(|| "pool|powder|ref|reference|bridge".to_owned());
        let irs_stat = self.irs_stat.unwrap_or_else(|| "median".to_owned());

        if self.ratio_fraction_merge.is_some() && quantification != QuantMethod::Ratio {
            return Err(MokumeError::InvalidInput {
                message: "--ratio-fraction-merge only applies to --quant-method ratio".to_owned(),
            });
        }
        let ratio_fraction_merge = self
            .ratio_fraction_merge
            .unwrap_or_else(|| "mean".to_owned());

        let imputation_options_supplied = self.impute_method.is_some()
            || self.impute_quantile.is_some()
            || self.impute_shift.is_some()
            || self.impute_scale.is_some()
            || self.impute_n_neighbors.is_some();
        if !self.impute && imputation_options_supplied {
            return Err(MokumeError::InvalidInput {
                message: "imputation options require --impute".to_owned(),
            });
        }
        let impute_method = self.impute_method.unwrap_or_else(|| "none".to_owned());
        if self.impute && impute_method.eq_ignore_ascii_case("none") {
            return Err(MokumeError::InvalidInput {
                message: "--impute requires --impute-method".to_owned(),
            });
        }
        let impute_method_lower = impute_method.to_ascii_lowercase();
        if self.impute_quantile.is_some()
            && !matches!(impute_method_lower.as_str(), "mindet" | "minprob")
        {
            return Err(MokumeError::InvalidInput {
                message: "--impute-quantile only applies to mindet/minprob".to_owned(),
            });
        }
        if (self.impute_shift.is_some() || self.impute_scale.is_some())
            && impute_method_lower != "minprob"
        {
            return Err(MokumeError::InvalidInput {
                message: "--impute-shift/--impute-scale only apply to minprob".to_owned(),
            });
        }
        if self.impute_n_neighbors.is_some()
            && !matches!(impute_method_lower.as_str(), "knn" | "seqknn")
        {
            return Err(MokumeError::InvalidInput {
                message: "--impute-n-neighbors only applies to knn/seqknn".to_owned(),
            });
        }

        let de_options_supplied = self.de_contrasts.is_some()
            || self.de_contrasts_file.is_some()
            || self.de_method.is_some()
            || self.de_ensemble_methods.is_some()
            || self.de_ensemble_min_k.is_some()
            || self.de_log2fc_threshold.is_some()
            || self.de_effect_size_gate.is_some()
            || self.de_fdr_threshold.is_some()
            || self.de_fdr_method.is_some()
            || self.de_output.is_some();
        if !self.differential_expression && de_options_supplied {
            return Err(MokumeError::InvalidInput {
                message: "differential-expression options require --de".to_owned(),
            });
        }
        let de_method = self.de_method.unwrap_or_else(|| "auto".to_owned());
        if self.de_ensemble_min_k.is_some() && !de_method.eq_ignore_ascii_case("ensemble") {
            return Err(MokumeError::InvalidInput {
                message: "--de-ensemble-min-k only applies to --de-method ensemble".to_owned(),
            });
        }
        let resolved_de_method = if de_method.eq_ignore_ascii_case("auto") {
            if quantification == QuantMethod::DirectLfq {
                "deqms"
            } else {
                "limrots"
            }
        } else {
            de_method.as_str()
        };
        if self.de_fdr_method.is_some()
            && matches!(
                resolved_de_method.to_ascii_lowercase().as_str(),
                "rots" | "limrots"
            )
        {
            return Err(MokumeError::InvalidInput {
                message: format!(
                    "--de-fdr-method does not apply to {resolved_de_method}, which retains its permutation FDR"
                ),
            });
        }
        let de_fdr_method = self.de_fdr_method.unwrap_or_else(|| "bh".to_owned());
        let (log2fc_threshold, auto_effect_size_gate) = self
            .de_log2fc_threshold
            .unwrap_or(DeLog2FcArg::Fixed(0.5))
            .into_config();
        let effect_size_gate = self.de_effect_size_gate.or(auto_effect_size_gate);

        Ok(FeatureToProteinsConfig {
            input: InputConfig {
                parquet: self.parquet,
                msstats: self.msstats,
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
                run_method: run_normalization,
                sample_method: sample_normalization,
                normalization_proteins: self.normalization_proteins,
            },
            quantification,
            topn_peptides,
            maxlfq: MaxLfqConfig {
                ion_alignment: None,
                force_builtin: false,
            },
            pibaq: PibaqConfig {
                enzyme: self.pibaq_enzyme,
                max_aa: self.pibaq_max_aa,
                min_shared: self.pibaq_min_shared,
                families_yaml: self.pibaq_families_yaml,
                min_anchors: self.pibaq_min_anchors,
                high_anchor_threshold: self.pibaq_high_anchor_threshold,
            },
            directlfq: DirectLfqConfig {
                cores: self.directlfq_cores,
                min_nonan: self.directlfq_min_nonan.unwrap_or(1),
                num_samples_quadratic: self.directlfq_num_samples_quadratic.unwrap_or(50),
            },
            batch: BatchCorrectionConfig {
                enabled: self.batch_correction,
                method: batch_method,
                column: self.batch_column,
                covariates: split_csv_option(self.batch_covariates),
                parametric: !self.batch_nonparametric,
                mean_only: self.batch_mean_only,
                ref_batch: self.batch_ref,
            },
            irs: IrsConfig {
                enabled: self.irs,
                reference_samples,
                sdrf_column: self.irs_sdrf_column,
                sdrf_values: split_csv_option(self.irs_sdrf_values),
                reference_regex: irs_reference_regex,
                stat: irs_stat,
                remove_reference: self.irs_remove_reference,
            },
            coverage_threshold: self.coverage_threshold,
            sample_correlation_threshold: self.min_sample_correlation,
            ratio: RatioConfig {
                fraction_merge: ratio_fraction_merge,
            },
            imputation: ImputationConfig {
                enabled: self.impute,
                method: impute_method,
                quantile: self.impute_quantile.unwrap_or(0.01),
                shift: self.impute_shift.unwrap_or(1.6),
                scale: self.impute_scale.unwrap_or(0.3),
                n_neighbors: self.impute_n_neighbors.unwrap_or(5),
            },
            differential_expression: DifferentialExpressionConfig {
                enabled: self.differential_expression,
                contrasts: split_csv_option(self.de_contrasts),
                contrasts_file: self.de_contrasts_file,
                method: de_method,
                ensemble_methods: split_ensemble_methods(self.de_ensemble_methods),
                ensemble_min_k: self.de_ensemble_min_k.unwrap_or(2),
                log2fc_threshold,
                effect_size_gate,
                fdr_threshold: self.de_fdr_threshold.unwrap_or(0.05),
                fdr_method: de_fdr_method,
                output: self.de_output,
            },
            runtime: RuntimeConfig {
                memory: self.memory,
                threads: self.threads,
            },
        })
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

/// The fixed `--quant-method` names, i.e. every method whose name carries no
/// parameter. The TopN family is spelled `top<N>` and is not listed here.
const FIXED_QUANT_METHODS: &str = "directlfq, pibaq, maxlfq, sum, median, ratio, abd, intensity, \
spectral_count";

/// Default `topn_peptides` for methods outside the TopN family. The pipeline
/// only reads the field when the method is [`QuantMethod::TopN`], in which case
/// N always comes from the `top<N>` name, so this value is inert -- it just
/// keeps the config field populated.
const DEFAULT_TOPN_PEPTIDES: usize = 3;

/// A fixed fold-change threshold or Python-compatible `auto` estimation.
#[derive(Debug, Clone, Copy)]
enum DeLog2FcArg {
    Fixed(f64),
    Auto,
}

impl DeLog2FcArg {
    fn into_config(self) -> (f64, Option<String>) {
        match self {
            Self::Fixed(value) => (value, None),
            Self::Auto => (0.5, Some("mixture".to_string())),
        }
    }
}

impl FromStr for DeLog2FcArg {
    type Err = String;

    fn from_str(value: &str) -> std::result::Result<Self, Self::Err> {
        if value.trim().eq_ignore_ascii_case("auto") {
            return Ok(Self::Auto);
        }
        let threshold = value.parse::<f64>().map_err(|_| {
            format!("invalid log2FC threshold `{value}`: expected `auto` or a non-negative number")
        })?;
        if !threshold.is_finite() || threshold < 0.0 {
            return Err(format!(
                "invalid log2FC threshold `{value}`: expected `auto` or a finite, non-negative number"
            ));
        }
        Ok(Self::Fixed(threshold))
    }
}

fn parse_de_log2fc(value: &str) -> std::result::Result<DeLog2FcArg, String> {
    DeLog2FcArg::from_str(value)
}

/// A validated `--quant-method` value: the parsed [`QuantMethod`] plus, for the
/// `top<N>` family, the N spelled in the name.
///
/// This is a plain `FromStr` newtype rather than a clap `ValueEnum` because
/// `top<N>` is an open-ended family (`top1`, `top3`, `top10`, ...) that no fixed
/// variant list can express. Parsing here (instead of keeping a raw `String` and
/// re-parsing later) means an invalid method is rejected by clap at parse time,
/// with the same exit code as any other bad option.
#[derive(Debug, Clone, Copy)]
struct QuantMethodArg {
    method: QuantMethod,
    /// `Some(N)` only for the `top<N>` family.
    topn: Option<usize>,
}

impl FromStr for QuantMethodArg {
    type Err = String;

    fn from_str(value: &str) -> std::result::Result<Self, Self::Err> {
        let lowered = value.trim().to_ascii_lowercase();
        // `topn` keeps the placeholder letter and means the canonical Top3
        // (Silva 2006), matching what the Python factory does with a `top` name
        // that carries no digits. Normalizing here keeps `top<digits>` as the
        // single internal spelling.
        if lowered == "topn" {
            return Ok(Self {
                method: QuantMethod::TopN,
                topn: Some(3),
            });
        }
        if let Some(topn) = parse_topn_from_method_name(&lowered) {
            return Ok(Self {
                method: QuantMethod::TopN,
                topn: Some(topn),
            });
        }
        // Anything else starting with `top` is a malformed TopN name (`top0`,
        // `topx`, ...); say so rather than reporting a generic unknown method.
        if lowered.starts_with("top") {
            return Err(invalid_topn_message(value));
        }
        let method = QuantMethod::from_str(&lowered).map_err(|_| unknown_method_message(value))?;
        Ok(Self { method, topn: None })
    }
}

/// clap `value_parser` for `--quant-method`. A `fn(&str) -> Result<_, String>`
/// is the form clap accepts directly, and it renders the `String` as the usage
/// error message.
fn parse_quant_method(value: &str) -> std::result::Result<QuantMethodArg, String> {
    QuantMethodArg::from_str(value)
}

/// Methods `peptides2protein` implements. It runs on an already-summarized
/// peptide table, so it offers a smaller set than `features2proteins`.
const PEPTIDES2PROTEIN_METHODS: [&str; 4] = ["pibaq", "maxlfq", "sum", "directlfq"];

/// clap `value_parser` for `peptides2protein --method`.
///
/// Applies the same TopN spelling rules as `--quant-method` over this command's
/// smaller method set, and normalizes bare `topn` to `top<DEFAULT_TOPN_PEPTIDES>`
/// so the runner only ever matches on `top<digits>`.
fn parse_peptides2protein_method(value: &str) -> std::result::Result<String, String> {
    let lowered = value.trim().to_ascii_lowercase();
    if PEPTIDES2PROTEIN_METHODS.contains(&lowered.as_str()) {
        return Ok(lowered);
    }
    if lowered == "topn" {
        return Ok(format!("top{DEFAULT_TOPN_PEPTIDES}"));
    }
    if parse_topn_from_method_name(&lowered).is_some() {
        return Ok(lowered);
    }
    // A `top`-prefixed name that carries no usable N is a malformed TopN request,
    // not an unknown method; say which of the two it is.
    if lowered.starts_with("top") {
        return Err(invalid_topn_message(value));
    }
    Err(format!(
        "unknown peptides2protein method `{value}`: expected one of {}, or `top<N>` (e.g. `top3`)",
        PEPTIDES2PROTEIN_METHODS.join(", ")
    ))
}

fn parse_nonzero_threads(value: &str) -> std::result::Result<i32, String> {
    let threads = value
        .parse::<i32>()
        .map_err(|_| format!("invalid thread count `{value}`: expected a non-zero integer"))?;
    if threads == 0 {
        return Err("invalid thread count `0`: expected a non-zero integer".to_owned());
    }
    Ok(threads)
}

fn parse_memory(value: &str) -> std::result::Result<String, String> {
    parse_memory_to_bytes(value)
        .map(|_| value.to_owned())
        .map_err(|_| {
            format!("invalid memory value `{value}`: expected a positive size such as 1GB or 512MB")
        })
}

fn parse_positive_i32(value: &str) -> std::result::Result<i32, String> {
    let parsed = value
        .parse::<i32>()
        .map_err(|_| format!("invalid positive integer `{value}`"))?;
    if parsed <= 0 {
        return Err(format!("expected a positive integer, got `{value}`"));
    }
    Ok(parsed)
}

fn parse_positive_f64(value: &str) -> std::result::Result<f64, String> {
    let parsed = value
        .parse::<f64>()
        .map_err(|_| format!("invalid positive number `{value}`"))?;
    if !parsed.is_finite() || parsed <= 0.0 {
        return Err(format!("expected a finite positive number, got `{value}`"));
    }
    Ok(parsed)
}

fn parse_nonnegative_i32(value: &str) -> std::result::Result<i32, String> {
    let parsed = value
        .parse::<i32>()
        .map_err(|_| format!("invalid non-negative integer `{value}`"))?;
    if parsed < 0 {
        return Err(format!("expected a non-negative integer, got `{value}`"));
    }
    Ok(parsed)
}

fn parse_finite_f64(value: &str) -> std::result::Result<f64, String> {
    let parsed = value
        .parse::<f64>()
        .map_err(|_| format!("invalid finite number `{value}`"))?;
    if !parsed.is_finite() {
        return Err(format!("expected a finite number, got `{value}`"));
    }
    Ok(parsed)
}

fn parse_nonnegative_f64(value: &str) -> std::result::Result<f64, String> {
    let parsed = parse_finite_f64(value)?;
    if parsed < 0.0 {
        return Err(format!("expected a non-negative number, got `{value}`"));
    }
    Ok(parsed)
}

fn parse_fraction(value: &str) -> std::result::Result<f64, String> {
    let parsed = parse_finite_f64(value)?;
    if !(0.0..=1.0).contains(&parsed) {
        return Err(format!("expected a number between 0 and 1, got `{value}`"));
    }
    Ok(parsed)
}

fn parse_correlation(value: &str) -> std::result::Result<f64, String> {
    let parsed = parse_finite_f64(value)?;
    if !(-1.0..=1.0).contains(&parsed) {
        return Err(format!("expected a number between -1 and 1, got `{value}`"));
    }
    Ok(parsed)
}

fn parse_named_score_filter(value: &str) -> std::result::Result<NamedScoreFilterConfig, String> {
    let (name, threshold) = value
        .rsplit_once('=')
        .ok_or_else(|| format!("invalid score filter `{value}`: expected NAME=THRESHOLD"))?;
    let name = name.trim();
    if name.is_empty() {
        return Err("score filter name cannot be empty".to_owned());
    }
    let threshold = threshold.trim().parse::<f64>().map_err(|_| {
        format!("invalid score filter `{value}`: threshold must be a finite number")
    })?;
    if !threshold.is_finite() {
        return Err(format!(
            "invalid score filter `{value}`: threshold must be a finite number"
        ));
    }
    Ok(NamedScoreFilterConfig {
        name: name.to_owned(),
        threshold,
    })
}

fn parse_positive_usize(value: &str) -> std::result::Result<usize, String> {
    let parsed = value
        .parse::<usize>()
        .map_err(|_| format!("invalid positive integer `{value}`"))?;
    if parsed == 0 {
        return Err("invalid positive integer `0`".to_owned());
    }
    Ok(parsed)
}

/// Error for a `top`-prefixed name that is not a valid `top<N>` (`top0`, `topx`).
fn invalid_topn_message(value: &str) -> String {
    format!(
        "invalid quantification method `{value}`: a TopN method is spelled `top<N>` with N >= 1 \
(e.g. `top1`, `top3`, `top5`)"
    )
}

/// Error for a name that is neither a fixed method nor a `top<N>`.
fn unknown_method_message(value: &str) -> String {
    format!("unknown quantification method `{value}`: expected one of {FIXED_QUANT_METHODS}, or `top<N>` (e.g. `top3`)")
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

/// Dispatch a fully-built [`Cli`] to its subcommand.
fn dispatch(cli: Cli, pibaq_digest: Option<PibaqDigest>) -> mokume_core::Result<()> {
    init_logging(cli.log_level, cli.log_file).and_then(|()| match cli.command {
        Commands::Features2Proteins(args) => args
            .into_config()
            .and_then(|config| dispatch_features_to_proteins(config, pibaq_digest)),
        Commands::Features2Peptides(args) => dispatch_features_to_peptides(&args),
        Commands::Peptides2Protein(args) => {
            peptides2protein::run_peptides_to_protein_with_digest(&args, pibaq_digest)
        }
        Commands::CorrectBatches(args) => correct_batches::run_correct_batches(&args),
    })
}

/// Parse an argv vector and return the runtime pyOpenMS digest request when the
/// selected command is piBAQ. Parse/help errors return `None`; the normal CLI
/// entry point remains authoritative for rendering those errors.
pub fn pibaq_digest_request_from_args<I, T>(args: I) -> Option<PibaqDigestRequest>
where
    I: IntoIterator<Item = T>,
    T: Into<OsString> + Clone,
{
    let cli = Cli::try_parse_from(args).ok()?;
    match cli.command {
        Commands::Features2Proteins(args) if args.quant_method.method == QuantMethod::Pibaq => {
            let fasta = args.fasta?;
            fasta.is_file().then_some(PibaqDigestRequest {
                fasta,
                enzyme: args.pibaq_enzyme,
                min_aa: args.min_aa,
                max_aa: args.pibaq_max_aa,
                missed_cleavages: 0,
            })
        }
        Commands::Peptides2Protein(args)
            if args.method.eq_ignore_ascii_case("pibaq") && args.output.is_some() =>
        {
            let fasta = args.fasta?;
            fasta.is_file().then_some(PibaqDigestRequest {
                fasta,
                enzyme: args.enzyme,
                min_aa: args.min_aa,
                max_aa: args.max_aa,
                missed_cleavages: 0,
            })
        }
        _ => None,
    }
}

/// Library entry point used by the internal `mokume-py` PyO3 crate: parse an
/// explicit argument vector and return the dispatch result. A parse error is
/// returned as `Err` rather than exiting the process, so it never tears down a
/// hosting Python interpreter.
pub fn run_from_args<I, T>(args: I) -> mokume_core::Result<()>
where
    I: IntoIterator<Item = T>,
    T: Into<OsString> + Clone,
{
    let cli = Cli::try_parse_from(args).map_err(|error| MokumeError::InvalidInput {
        message: error.to_string(),
    })?;
    dispatch(cli, None)
}

/// Parse and run a command with a runtime pyOpenMS theoretical-peptide map.
pub fn run_from_args_with_pibaq_digest<I, T>(
    args: I,
    digest: PibaqDigest,
) -> mokume_core::Result<()>
where
    I: IntoIterator<Item = T>,
    T: Into<OsString> + Clone,
{
    let cli = Cli::try_parse_from(args).map_err(|error| MokumeError::InvalidInput {
        message: error.to_string(),
    })?;
    dispatch(cli, Some(digest))
}

/// Console-script entry point for the `mokume` wheel: parse an explicit argument
/// vector and return the process exit code WITHOUT calling `process::exit`, so it
/// never tears down a hosting Python interpreter. clap's help/version are printed
/// to stdout (exit 0) and usage errors to stderr (exit 2); a dispatch failure
/// prints the error and returns 1.
pub fn run_cli_from_args<I, T>(args: I) -> i32
where
    I: IntoIterator<Item = T>,
    T: Into<OsString> + Clone,
{
    match Cli::try_parse_from(args) {
        Ok(cli) => match dispatch(cli, None) {
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

/// Console-script entry point with a runtime pyOpenMS theoretical-peptide map.
pub fn run_cli_from_args_with_pibaq_digest<I, T>(args: I, digest: PibaqDigest) -> i32
where
    I: IntoIterator<Item = T>,
    T: Into<OsString> + Clone,
{
    match Cli::try_parse_from(args) {
        Ok(cli) => match dispatch(cli, Some(digest)) {
            Ok(()) => 0,
            Err(error) => {
                eprintln!("{error}");
                1
            }
        },
        Err(error) => {
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
fn dispatch_features_to_proteins(
    config: FeatureToProteinsConfig,
    pibaq_digest: Option<PibaqDigest>,
) -> mokume_core::Result<()> {
    match pibaq_digest {
        Some(digest) => run_features_to_proteins_with_pibaq_digest(&config, digest),
        None => run_features_to_proteins(&config),
    }
}

fn dispatch_features_to_peptides(args: &Features2PeptidesArgs) -> mokume_core::Result<()> {
    // `--generate-filter-config` writes the example config and exits, mirroring
    // the Python command's early return.
    if let Some(path) = &args.generate_filter_config {
        generate_filter_config(path)?;
        return Ok(());
    }

    let parquet = args
        .parquet
        .as_ref()
        .ok_or_else(|| MokumeError::InvalidInput {
            message: "features2peptides requires --parquet".to_owned(),
        })?;
    if !parquet.exists() {
        return Err(MokumeError::MissingInput {
            path: parquet.clone(),
        });
    }
    if let Some(sdrf) = &args.sdrf {
        if !sdrf.exists() {
            return Err(MokumeError::MissingInput { path: sdrf.clone() });
        }
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

    if args.irs_channel.is_some() && args.irs_autodetect_regex.is_some() {
        return Err(MokumeError::InvalidInput {
            message: "choose either --irs_channel or --irs_autodetect_regex, not both".to_owned(),
        });
    }
    if args.irs_autodetect_regex.is_some() && args.sdrf.is_none() {
        return Err(MokumeError::InvalidInput {
            message: "--irs_autodetect_regex requires --sdrf".to_owned(),
        });
    }
    let irs_requested = args.irs_channel.is_some() || args.irs_autodetect_regex.is_some();
    if !irs_requested
        && (!args.irs_stat.eq_ignore_ascii_case("median")
            || !args.irs_scope.eq_ignore_ascii_case("global"))
    {
        return Err(MokumeError::InvalidInput {
            message: "--irs_stat/--irs_scope require --irs_channel or --irs_autodetect_regex"
                .to_owned(),
        });
    }

    // Channel IRS (Python `get_irs_scaling_factors`). All three scopes
    // (`global` / `by_mixture` / `two_stage`) are ported. The reference channel
    // is taken from `--irs_channel`; when that is absent but
    // `--irs_autodetect_regex` and `--sdrf` are given, it is autodetected from
    // the SDRF exactly as Python does (`peptide.py:219-233`). A requested rule
    // that finds no channel is rejected instead of becoming a successful no-op.
    // `--save_parquet` and `--aggregation_level run` are implemented below.
    let irs_channel = match (&args.irs_channel, &args.irs_autodetect_regex, &args.sdrf) {
        (Some(channel), _, _) => Some(channel.clone()),
        (None, Some(regex), Some(sdrf)) => {
            Some(resolve_irs_autodetect_channel(sdrf, regex)?.ok_or_else(|| {
                MokumeError::InvalidInput {
                    message: format!(
                        "--irs_autodetect_regex `{regex}` matched no reference channel"
                    ),
                }
            })?)
        }
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
            parquet: Some(parquet.clone()),
            msstats: None,
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
        run_normalization: args
            .run_normalization
            .clone()
            .unwrap_or_else(|| "median".to_owned()),
        sample_normalization: args
            .sample_normalization
            .clone()
            .unwrap_or_else(|| "globalmedian".to_owned()),
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
    let has_override = has_filter_override(args);

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
    if let Some(value) = &args.filter_score {
        config.peptide.score = Some(value.clone());
    }
    apply_fdr_overrides(args, &mut config);
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

fn has_filter_override(args: &Features2PeptidesArgs) -> bool {
    args.filter_min_intensity.is_some()
        || args.filter_cv_threshold.is_some()
        || args.filter_charge_states.is_some()
        || args.filter_max_missed_cleavages.is_some()
        || args.filter_peptide_fdr.is_some()
        || args.filter_score.is_some()
        || args.filter_exclude_modifications.is_some()
        || args.filter_min_unique_peptides.is_some()
        || args.filter_protein_fdr.is_some()
        || args.filter_min_features.is_some()
        || args.filter_max_missing_rate.is_some()
}

fn apply_fdr_overrides(args: &Features2PeptidesArgs, config: &mut PreprocessingFilterConfig) {
    if let Some(value) = args.filter_peptide_fdr {
        config.peptide.fdr_threshold = Some(value);
    }
    if let Some(value) = args.filter_protein_fdr {
        config.protein.fdr_threshold = Some(value);
    }
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

/// Initialize the global tracing subscriber once per process.
///
/// Library entry points may be called repeatedly from one Python process. The
/// same logging configuration reuses the installed subscriber; a different
/// configuration is rejected instead of silently pretending that it took
/// effect.
fn init_logging(level: LogLevel, log_file: Option<PathBuf>) -> mokume_core::Result<()> {
    static CONFIG: Mutex<Option<(LogLevel, Option<PathBuf>)>> = Mutex::new(None);
    let requested = (level, log_file);
    let mut active = CONFIG.lock().map_err(|_| MokumeError::InvalidInput {
        message: "logging configuration lock is poisoned".to_owned(),
    })?;

    if let Some(config) = active.as_ref() {
        if config == &requested {
            return Ok(());
        }
        return Err(MokumeError::InvalidInput {
            message: format!(
                "logging is already initialized with log_level={} and log_file={:?}; \
                 requested log_level={} and log_file={:?}",
                config.0.as_filter(),
                config.1,
                requested.0.as_filter(),
                requested.1,
            ),
        });
    }

    init_logging_once(requested.0, requested.1.clone())?;
    *active = Some(requested);
    Ok(())
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
        tracing_subscriber::fmt()
            .with_env_filter(filter)
            .with_ansi(false)
            .with_writer(Mutex::new(file))
            .try_init()
            .map_err(|error| MokumeError::InvalidInput {
                message: format!("failed to initialize logging: {error}"),
            })?;
    } else {
        tracing_subscriber::fmt()
            .with_env_filter(filter)
            .try_init()
            .map_err(|error| MokumeError::InvalidInput {
                message: format!("failed to initialize logging: {error}"),
            })?;
    }
    Ok(())
}

const EXAMPLE_FILTER_CONFIG_YAML: &str = r#"# Mokume Preprocessing Filter Configuration
# This file defines quality filters applied during peptide normalization

name: example_config

# Global options
enabled: true              # Set to false to disable all filtering

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
  allowed_charge_states: null     # e.g., [2, 3, 4] or null for all charges
  exclude_modifications: []       # Modification names to exclude, e.g., ["Oxidation"]
  max_missed_cleavages: null      # Max missed cleavages (null = no filter)
  fdr_threshold: null              # Peptide q-value cutoff (null = no filter)
  min_peptide_length: 7           # Minimum peptide length in amino acids
  max_peptide_length: 50          # Maximum peptide length in amino acids
  exclude_sequence_patterns: []   # Regex patterns to exclude

# Protein-level filters
protein:
  fdr_threshold: null           # Protein-group q-value cutoff (null = no filter)
  min_unique_peptides: 2      # Minimum unique peptides per protein
  razor_peptide_handling: keep   # How to handle shared peptides: keep, remove, assign_to_top
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
  max_missing_rate: 1.0         # Max missing value rate per run (0-1)
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
    "allowed_charge_states": null,
    "exclude_modifications": [],
    "max_missed_cleavages": null,
    "fdr_threshold": null,
    "min_peptide_length": 7,
    "max_peptide_length": 50,
    "exclude_sequence_patterns": []
  },
  "protein": {
    "fdr_threshold": null,
    "min_unique_peptides": 2,
    "razor_peptide_handling": "keep",
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
    "max_missing_rate": 1.0
  },
  "enabled": true
}
"#;

#[cfg(test)]
mod tests {
    use std::path::Path;

    use clap::{CommandFactory, Parser};

    use mokume_core::{MokumeError, PreprocessingFilterConfig};

    use super::{Cli, Commands};

    #[test]
    fn correct_batches_accepts_short_column_aliases() {
        // Python's click CLI accepts `-sid`/`-pid`/`-pibaq`; clap cannot express
        // single-dash multi-character shorts, so these are offered as `--sid` /
        // `--pid` / `--pibaq` long aliases (closest portable form).
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
            "--pibaq",
            "MyPibaq",
        ]);
        let Commands::CorrectBatches(args) = cli.command else {
            panic!("expected the correct-batches subcommand");
        };
        assert_eq!(args.sample_id_column, "MySample");
        assert_eq!(args.protein_id_column, "MyProtein");
        assert_eq!(args.pibaq_raw_column, "MyPibaq");
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
            "--min-sample-correlation",
            "0.90",
            "--impute",
            "--impute-method",
            "knn",
            "--de",
            "--de-method",
            "limma",
            "--de-contrasts",
            "NASH-HL,NASH-Ctrl",
            "--de-fdr-method",
            "ihw",
            "--de-output",
            "de.csv",
            "--memory",
            "1GB",
            "--duckdb-threads",
            "24",
        ]);

        let Commands::Features2Proteins(args) = cli.command else {
            panic!("expected features2proteins command");
        };
        let Ok(config) = args.into_config() else {
            panic!("expected a valid features2proteins config");
        };

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
        assert_eq!(config.sample_correlation_threshold, Some(0.90));
        assert!(config.imputation.enabled);
        assert_eq!(config.imputation.method, "knn");
        assert!(config.differential_expression.enabled);
        assert_eq!(config.differential_expression.fdr_method, "ihw");
        assert_eq!(config.runtime.memory.as_deref(), Some("1GB"));
        assert_eq!(config.runtime.threads, Some(24));
    }

    #[test]
    fn rejects_removed_ibaq_method_name() {
        let features = Cli::try_parse_from([
            "mokume",
            "features2proteins",
            "-p",
            "input.parquet",
            "-o",
            "protein.csv",
            "--quant-method",
            "ibaq",
        ]);
        assert!(features.is_err());

        let peptides = Cli::try_parse_from([
            "mokume",
            "peptides2protein",
            "--peptides",
            "peptides.csv",
            "--method",
            "ibaq",
        ]);
        assert!(peptides.is_err());
    }

    #[test]
    fn rejects_nonpositive_memory_budget() {
        let parsed = Cli::try_parse_from([
            "mokume",
            "features2proteins",
            "--parquet",
            "input.parquet",
            "--output",
            "protein.csv",
            "--memory",
            "0GB",
        ]);

        let Err(error) = parsed else {
            panic!("zero memory budget was accepted");
        };
        assert!(error.to_string().contains("invalid memory value `0GB`"));
    }

    #[test]
    fn parses_adaptive_de_options() {
        let cli = Cli::parse_from([
            "mokume",
            "features2proteins",
            "-p",
            "input.parquet",
            "-o",
            "protein.csv",
            "-s",
            "input.sdrf.tsv",
            "--de",
            "--de-method",
            "limma",
            "--de-contrasts",
            "A vs B",
            "--de-log2fc",
            "auto",
            "--de-fdr-method",
            "storey",
        ]);
        let Commands::Features2Proteins(args) = cli.command else {
            panic!("expected features2proteins command");
        };
        let Ok(config) = args.into_config() else {
            panic!("expected a valid features2proteins config");
        };

        assert_eq!(config.differential_expression.log2fc_threshold, 0.5);
        assert_eq!(
            config.differential_expression.effect_size_gate.as_deref(),
            Some("mixture")
        );
        assert_eq!(config.differential_expression.fdr_method, "storey");
    }

    #[test]
    fn explicit_effect_size_method_uses_numeric_threshold_as_fallback() {
        let cli = Cli::parse_from([
            "mokume",
            "features2proteins",
            "-p",
            "input.parquet",
            "-o",
            "protein.csv",
            "--de",
            "--de-log2fc",
            "0.25",
            "--de-effect-size-gate",
            "null_quantile",
        ]);
        let Commands::Features2Proteins(args) = cli.command else {
            panic!("expected features2proteins command");
        };
        let Ok(config) = args.into_config() else {
            panic!("expected a valid features2proteins config");
        };

        assert_eq!(config.differential_expression.log2fc_threshold, 0.25);
        assert_eq!(
            config.differential_expression.effect_size_gate.as_deref(),
            Some("null_quantile")
        );
    }

    #[test]
    fn parses_native_msstats_input() {
        let cli = Cli::parse_from([
            "mokume",
            "features2proteins",
            "--msstats",
            "input.csv",
            "--sdrf",
            "input.sdrf.tsv",
            "--output",
            "protein.csv",
        ]);
        let Commands::Features2Proteins(args) = cli.command else {
            panic!("expected features2proteins command");
        };
        let Ok(config) = args.into_config() else {
            panic!("expected a valid features2proteins config");
        };

        assert_eq!(
            config.input.msstats.as_deref(),
            Some(Path::new("input.csv"))
        );
        assert!(config.input.parquet.is_none());
    }

    #[test]
    fn native_msstats_input_requires_sdrf() {
        let Err(error) = Cli::try_parse_from([
            "mokume",
            "features2proteins",
            "--msstats",
            "input.csv",
            "--output",
            "protein.csv",
        ]) else {
            panic!("MSstats without SDRF should fail");
        };

        assert_eq!(
            error.kind(),
            clap::error::ErrorKind::MissingRequiredArgument
        );
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
                "--de",
                "--de-method",
                "ensemble",
                "--de-ensemble-methods",
                value,
            ]);
            let Commands::Features2Proteins(args) = cli.command else {
                panic!("expected features2proteins command");
            };
            let Ok(config) = args.into_config() else {
                panic!("expected a valid features2proteins config");
            };
            let methods = config.differential_expression.ensemble_methods;

            assert!(
                methods
                    .as_ref()
                    .is_some_and(|methods| methods.iter().any(String::is_empty)),
                "empty member from {value:?} was discarded: {methods:?}"
            );
        }
    }

    #[test]
    fn rejects_ensemble_min_k_for_single_de_method() {
        let cli = Cli::parse_from([
            "mokume",
            "features2proteins",
            "-p",
            "input.parquet",
            "-o",
            "protein.csv",
            "--de-method",
            "limma",
            "--de-ensemble-min-k",
            "2",
        ]);
        let Commands::Features2Proteins(args) = cli.command else {
            panic!("expected features2proteins command");
        };

        assert!(matches!(
            args.into_config(),
            Err(MokumeError::InvalidInput { .. })
        ));
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
            "--irs",
            "--irs-reference-sample",
            "Pool, batch A",
            "--irs-reference-sample",
            "Pool B",
        ]);

        let Commands::Features2Proteins(args) = cli.command else {
            panic!("expected features2proteins command");
        };
        let Ok(config) = args.into_config() else {
            panic!("expected a valid features2proteins config");
        };

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
            "--min-aa",
            "--min-unique",
            "--keep-contaminants",
            "--run-normalization",
            "--sample-normalization",
            "--normalization-proteins",
            "--fasta",
            "--pibaq-enzyme",
            "--pibaq-max-aa",
            "--pibaq-min-shared",
            "--pibaq-families",
            "--pibaq-min-anchors",
            "--pibaq-high-anchor-threshold",
            "--directlfq-cores",
            "--directlfq-min-nonan",
            "--directlfq-num-samples-quadratic",
            "--export-peptides",
            "--export-ions",
            "--batch-correction",
            "--batch-method",
            "--batch-column",
            "--batch-covariates",
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
            "--min-sample-correlation",
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
            "--memory",
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
            "--ion-alignment",
            "--duckdb-memory",
            "--remove-contaminants",
            "--batch-parametric",
            "--impute-method missforest",
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
        // N is spelled in the method name (`top5`), so the companion option is
        // gone from both CLIs and must stay gone.
        assert!(
            !help.contains("--topn"),
            "removed option `--topn` must not appear in help:\n{help}"
        );
        assert!(!help.contains("--ibaq-"));
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
            "--filter-peptide-fdr",
            "--filter-score",
            "--filter-exclude-modifications",
            "--filter-min-unique-peptides",
            "--filter-protein-fdr",
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
    fn filter_config_rejects_unknown_keys() {
        for yaml in [
            "unknown_option: true\n",
            "intensity:\n  min_intensitty: 5\n",
        ] {
            let parsed = serde_yaml::from_str::<PreprocessingFilterConfig>(yaml);
            assert!(parsed.is_err(), "unknown filter key was accepted: {yaml}");
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
        // N is spelled in the method name (`--method top5`), so the companion
        // option is gone from both CLIs and must stay gone.
        assert!(
            !help.contains("--topn_n"),
            "removed option `--topn_n` must not appear in help:\n{help}"
        );
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
            "--pibaq_raw_column",
            "--pibaq_corrected_column",
            "--export_anndata",
        ] {
            assert!(
                correct_batches_help.contains(option),
                "missing option `{option}` in help:\n{correct_batches_help}"
            );
        }
        assert!(!correct_batches_help.contains("--ibaq_"));
        assert!(!correct_batches_help.contains("--ibaq-"));
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
