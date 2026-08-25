use clap::CommandFactory;

use crate::Cli;

const REQUIRED_FEATURES_TO_PROTEINS_OPTIONS: &[&str] = &[
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
];

const REMOVED_FEATURES_TO_PROTEINS_OPTIONS: &[&str] = &[
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
];

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
    for option in REQUIRED_FEATURES_TO_PROTEINS_OPTIONS.iter().copied() {
        assert!(
            help.contains(option),
            "missing option `{option}` in help:\n{help}"
        );
    }
    // The plotting / interactive-report flags moved to the Python wheel; clap
    // must reject them now.
    for option in REMOVED_FEATURES_TO_PROTEINS_OPTIONS.iter().copied() {
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
