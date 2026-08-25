use std::path::PathBuf;

use clap::Args;
use mokume_core::QuantMethod;

use crate::parsers::{
    parse_correlation, parse_de_log2fc, parse_finite_f64, parse_fraction, parse_memory,
    parse_nonnegative_f64, parse_nonnegative_i32, parse_positive_usize, parse_quant_method,
    DeLog2FcArg, OutputFormatArg, QuantMethodArg,
};
use crate::PibaqDigestRequest;

mod config;
#[derive(Debug, Args)]
pub(crate) struct Features2ProteinsArgs {
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
    pub(crate) fn into_config(self) -> mokume_core::Result<mokume_core::FeatureToProteinsConfig> {
        config::into_config(self)
    }

    pub(crate) fn into_pibaq_digest_request(self) -> Option<PibaqDigestRequest> {
        if self.quant_method.method != QuantMethod::Pibaq {
            return None;
        }
        let fasta = self.fasta?;
        fasta.is_file().then_some(PibaqDigestRequest {
            fasta,
            enzyme: self.pibaq_enzyme,
            min_aa: self.min_aa,
            max_aa: self.pibaq_max_aa,
            missed_cleavages: 0,
        })
    }
}
