use std::path::PathBuf;

use clap::Args;

use crate::parsers::{
    parse_peptides2protein_method, parse_positive_f64, parse_positive_i32, parse_positive_usize,
};

#[derive(Debug, Args)]
pub(crate) struct Peptides2ProteinArgs {
    #[arg(short = 'f', long = "fasta", value_name = "FILE")]
    pub(crate) fasta: Option<PathBuf>,

    #[arg(short = 'p', long = "peptides", value_name = "FILE")]
    pub(crate) peptides: PathBuf,

    #[arg(
        long = "quant-method",
        default_value = "pibaq",
        value_name = "METHOD",
        value_parser = parse_peptides2protein_method,
        help = "[possible values: pibaq, maxlfq, sum, directlfq, top<N> (e.g. top3)]"
    )]
    pub(crate) quant_method: String,

    #[arg(long = "enzyme", value_name = "NAME", default_value = "Trypsin")]
    pub(crate) enzyme: String,

    #[arg(long = "normalize")]
    pub(crate) normalize: bool,

    #[arg(long = "min-aa", value_name = "N", default_value_t = 7)]
    pub(crate) min_aa: usize,

    #[arg(long = "max-aa", value_name = "N", default_value_t = 30)]
    pub(crate) max_aa: usize,

    #[arg(long = "tpa")]
    pub(crate) tpa: bool,

    #[arg(long = "ruler")]
    pub(crate) ruler: bool,

    #[arg(long = "ploidy", value_name = "N", value_parser = parse_positive_i32)]
    pub(crate) ploidy: Option<i32>,

    #[arg(long = "organism", value_name = "NAME")]
    pub(crate) organism: Option<String>,

    #[arg(long = "cpc", value_name = "VALUE", value_parser = parse_positive_f64)]
    pub(crate) cpc: Option<f64>,

    #[arg(short = 'o', long = "output", value_name = "FILE", required = true)]
    pub(crate) output: Option<PathBuf>,

    #[arg(long = "qc-report", value_name = "FILE")]
    pub(crate) qc_report: Option<PathBuf>,

    #[arg(
        short = 't',
        long = "threads",
        value_name = "N",
        value_parser = parse_positive_usize,
        help = "DirectLFQ/MaxLFQ only"
    )]
    pub(crate) threads: Option<usize>,

    #[arg(
        long = "directlfq-min-nonan",
        value_name = "N",
        value_parser = parse_positive_usize
    )]
    pub(crate) directlfq_min_nonan: Option<usize>,

    #[arg(long = "families", value_name = "FILE")]
    pub(crate) families_yaml: Option<PathBuf>,

    #[arg(long = "min-shared", value_name = "N", default_value_t = 2)]
    pub(crate) min_shared: usize,

    #[arg(long = "min-anchors", value_name = "N", default_value_t = 1)]
    pub(crate) min_anchors: usize,

    #[arg(long = "high-anchor-threshold", value_name = "N", default_value_t = 3)]
    pub(crate) high_anchor_threshold: usize,
}

#[derive(Debug, Args)]
pub(crate) struct CorrectBatchesArgs {
    #[arg(short = 'i', long = "input", value_name = "DIR")]
    pub(crate) input: PathBuf,

    #[arg(
        short = 'p',
        long = "pattern",
        value_name = "GLOB",
        default_value = "*pibaq.tsv"
    )]
    pub(crate) pattern: String,

    #[arg(long = "comment", value_name = "PREFIX", default_value = "#")]
    pub(crate) comment: String,

    #[arg(long = "sep", value_name = "TEXT", default_value = "\t")]
    pub(crate) sep: String,

    #[arg(short = 'o', long = "output", value_name = "FILE")]
    pub(crate) output: PathBuf,

    #[arg(
        long = "sample-id-column",
        value_name = "COLUMN",
        default_value = "SampleID"
    )]
    pub(crate) sample_id_column: String,

    #[arg(
        long = "protein-id-column",
        value_name = "COLUMN",
        default_value = "ProteinName"
    )]
    pub(crate) protein_id_column: String,

    #[arg(
        long = "pibaq-raw-column",
        value_name = "COLUMN",
        default_value = "PiBAQ"
    )]
    pub(crate) pibaq_raw_column: String,

    #[arg(
        long = "pibaq-corrected-column",
        value_name = "COLUMN",
        default_value = "PiBAQBec"
    )]
    pub(crate) pibaq_corrected_column: String,

    #[arg(long = "export-anndata")]
    pub(crate) export_anndata: bool,
}
