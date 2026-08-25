use std::path::PathBuf;

use clap::Args;

use crate::parsers::{
    parse_nonzero_threads, parse_peptides2protein_method, parse_positive_f64, parse_positive_i32,
};

#[derive(Debug, Args)]
pub(crate) struct Peptides2ProteinArgs {
    #[arg(short = 'f', long = "fasta")]
    pub(crate) fasta: Option<PathBuf>,

    #[arg(short = 'p', long = "peptides")]
    pub(crate) peptides: PathBuf,

    #[arg(
        long = "method",
        default_value = "pibaq",
        value_name = "[pibaq|maxlfq|sum|directlfq|top<N>]",
        value_parser = parse_peptides2protein_method,
        help = "Quantification method: pibaq, maxlfq, sum, directlfq, or top<N> -- the TopN \
family spells its peptide count in the name (e.g. top3, top5)"
    )]
    pub(crate) method: String,

    #[arg(short = 'e', long = "enzyme", default_value = "Trypsin")]
    pub(crate) enzyme: String,

    #[arg(short = 'n', long = "normalize")]
    pub(crate) normalize: bool,

    #[arg(long = "min_aa", visible_alias = "min-aa", default_value_t = 7)]
    pub(crate) min_aa: usize,

    #[arg(long = "max_aa", visible_alias = "max-aa", default_value_t = 30)]
    pub(crate) max_aa: usize,

    #[arg(short = 't', long = "tpa")]
    pub(crate) tpa: bool,

    #[arg(short = 'r', long = "ruler")]
    pub(crate) ruler: bool,

    #[arg(short = 'i', long = "ploidy", value_parser = parse_positive_i32)]
    pub(crate) ploidy: Option<i32>,

    #[arg(short = 'm', long = "organism")]
    pub(crate) organism: Option<String>,

    #[arg(short = 'c', long = "cpc", value_parser = parse_positive_f64)]
    pub(crate) cpc: Option<f64>,

    #[arg(short = 'o', long = "output", required = true)]
    pub(crate) output: Option<PathBuf>,

    #[arg(long = "verbose")]
    pub(crate) verbose: bool,

    #[arg(
        long = "qc_report",
        visible_alias = "qc-report",
        default_value = "QCprofile.pdf"
    )]
    pub(crate) qc_report: PathBuf,

    #[arg(
        long = "threads",
        default_value_t = -1,
        value_parser = parse_nonzero_threads,
        help = "DirectLFQ/MaxLFQ worker count; negative values use joblib CPU-relative semantics"
    )]
    pub(crate) threads: i32,

    #[arg(long = "min_nonan", visible_alias = "min-nonan", default_value_t = 1)]
    pub(crate) min_nonan: usize,

    #[arg(long = "families")]
    pub(crate) families_yaml: Option<PathBuf>,

    #[arg(long = "min-shared", default_value_t = 2)]
    pub(crate) min_shared: usize,

    #[arg(long = "min-anchors", default_value_t = 1)]
    pub(crate) min_anchors: usize,

    #[arg(long = "high-anchor-threshold", default_value_t = 3)]
    pub(crate) high_anchor_threshold: usize,
}

#[derive(Debug, Args)]
pub(crate) struct CorrectBatchesArgs {
    #[arg(short = 'f', long = "folder")]
    pub(crate) folder: PathBuf,

    #[arg(short = 'p', long = "pattern", default_value = "*pibaq.tsv")]
    pub(crate) pattern: String,

    #[arg(long = "comment", default_value = "#")]
    pub(crate) comment: String,

    #[arg(long = "sep", default_value = "\t")]
    pub(crate) sep: String,

    #[arg(short = 'o', long = "output")]
    pub(crate) output: PathBuf,

    #[arg(
        long = "sample_id_column",
        visible_aliases = ["sample-id-column", "sid"],
        default_value = "SampleID"
    )]
    pub(crate) sample_id_column: String,

    #[arg(
        long = "protein_id_column",
        visible_aliases = ["protein-id-column", "pid"],
        default_value = "ProteinName"
    )]
    pub(crate) protein_id_column: String,

    #[arg(
        long = "pibaq_raw_column",
        visible_aliases = [
            "pibaq-raw-column",
            "pibaq"
        ],
        default_value = "PiBAQ"
    )]
    pub(crate) pibaq_raw_column: String,

    #[arg(
        long = "pibaq_corrected_column",
        visible_alias = "pibaq-corrected-column",
        default_value = "PiBAQBec"
    )]
    pub(crate) pibaq_corrected_column: String,

    #[arg(long = "export_anndata", visible_alias = "export-anndata")]
    pub(crate) export_anndata: bool,
}
