use std::{
    ffi::OsString,
    fs::{create_dir_all, read_to_string, write, File},
    path::{Path, PathBuf},
    str::FromStr,
    sync::Mutex,
};

use clap::{ArgAction, Args, Parser, Subcommand, ValueEnum};
use mokume_core::quant::parse_topn_from_method_name;
use mokume_core::{
    AggregationLevel, BatchCorrectionConfig, DifferentialExpressionConfig, DirectLfqConfig,
    FeatureToPeptidesConfig, FeatureToProteinsConfig, FilterConfig, ImputationConfig, InputConfig,
    IrsChannelConfig, IrsConfig, IrsScope, IrsStat, MaxLfqConfig, MokumeError, NormalizationConfig,
    OutputConfig, OutputFormat, PibaqConfig, PreprocessingFilterConfig, QuantMethod, RatioConfig,
    RuntimeConfig,
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

    #[arg(long = "de-log2fc", default_value = "0.5", value_parser = parse_de_log2fc)]
    de_log2fc_threshold: DeLog2FcArg,

    #[arg(long = "de-effect-size-gate", value_parser = ["mixture", "null_quantile"], ignore_case = true)]
    de_effect_size_gate: Option<String>,

    #[arg(long = "de-fdr", default_value_t = 0.05)]
    de_fdr_threshold: f64,

    #[arg(long = "de-fdr-method", default_value = "bh", value_parser = ["bh", "ihw", "bky", "storey"], ignore_case = true)]
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
        let remove_contaminants = if self.keep_contaminants {
            false
        } else {
            self.remove_contaminants || !self.keep_contaminants
        };

        // `top<N>` is the only spelling of a TopN method, so N comes from the
        // method name; the pipeline reads it from `topn_peptides`.
        let QuantMethodArg {
            method: quantification,
            topn,
        } = self.quant_method;
        let topn_peptides = topn.unwrap_or(DEFAULT_TOPN_PEPTIDES);
        let (log2fc_threshold, auto_effect_size_gate) = self.de_log2fc_threshold.into_config();
        let effect_size_gate = self.de_effect_size_gate.or(auto_effect_size_gate);

        FeatureToProteinsConfig {
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
                log2fc_threshold,
                effect_size_gate,
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
        Commands::Features2Proteins(args) => {
            dispatch_features_to_proteins(args.into_config(), pibaq_digest)
        }
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
            parquet: Some(args.parquet.clone()),
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
    use std::path::Path;

    use clap::{CommandFactory, Parser};

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
        let config = args.into_config();

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
            "--de-log2fc",
            "0.25",
            "--de-effect-size-gate",
            "null_quantile",
        ]);
        let Commands::Features2Proteins(args) = cli.command else {
            panic!("expected features2proteins command");
        };
        let config = args.into_config();

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
        let config = args.into_config();

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
            "--min-aa",
            "--min-unique",
            "--remove-contaminants",
            "--keep-contaminants",
            "--run-normalization",
            "--sample-normalization",
            "--normalization-proteins",
            "--fasta",
            "--ion-alignment",
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
