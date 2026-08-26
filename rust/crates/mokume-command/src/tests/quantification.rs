use std::path::Path;

use clap::Parser;

use crate::{Cli, Commands, Features2ProteinsArgs, QuantifyCommands};

#[test]
fn parses_true_spectral_count_psm_input() {
    let cli = Cli::parse_from([
        "mokume",
        "quantify",
        "features2proteins",
        "--psm",
        "input.psm.parquet",
        "--parquet",
        "input.feature.parquet",
        "--sdrf",
        "input.sdrf.tsv",
        "--output",
        "counts.csv",
        "--quant-method",
        "spectral-count",
    ]);
    let args = features_to_proteins_args(cli);
    let Ok(config) = args.into_config() else {
        panic!("expected a valid spectral-count config");
    };
    assert_eq!(
        config.input.psm.as_deref(),
        Some(Path::new("input.psm.parquet"))
    );
    assert_eq!(
        config.input.parquet.as_deref(),
        Some(Path::new("input.feature.parquet"))
    );
    assert_eq!(
        config.quantification,
        mokume_core::QuantMethod::SpectralCount
    );
    assert_eq!(config.normalization.run_method, "none");
    assert_eq!(config.normalization.sample_method, "none");
}

#[test]
fn spectral_count_requires_paired_qpx_inputs() {
    let Err(error) = Cli::try_parse_from([
        "mokume",
        "quantify",
        "features2proteins",
        "--psm",
        "input.psm.parquet",
        "--sdrf",
        "input.sdrf.tsv",
        "--output",
        "counts.csv",
        "--quant-method",
        "spectral-count",
    ]) else {
        panic!("spectral-count accepted a PSM without its feature QPX");
    };
    assert!(error.to_string().contains("--parquet"), "{error}");

    let cli = Cli::parse_from([
        "mokume",
        "quantify",
        "features2proteins",
        "--parquet",
        "input.feature.parquet",
        "--sdrf",
        "input.sdrf.tsv",
        "--output",
        "counts.csv",
        "--quant-method",
        "spectral-count",
    ]);
    let Err(error) = features_to_proteins_args(cli).into_config() else {
        panic!("spectral-count accepted a feature QPX without its PSM QPX");
    };
    assert!(error.to_string().contains("--psm and --parquet"), "{error}");
}

#[test]
fn feature_input_uses_explicit_peptide_count_name() {
    let cli = Cli::parse_from([
        "mokume",
        "quantify",
        "features2proteins",
        "--parquet",
        "input.feature.parquet",
        "--output",
        "counts.csv",
        "--quant-method",
        "peptide-count",
    ]);
    let args = features_to_proteins_args(cli);
    let Ok(config) = args.into_config() else {
        panic!("expected a valid peptide-count config");
    };
    assert_eq!(
        config.quantification,
        mokume_core::QuantMethod::PeptideCount
    );
    assert_eq!(config.normalization.run_method, "none");
    assert_eq!(config.normalization.sample_method, "none");
}

#[test]
fn count_methods_reject_irs() {
    for (method, inputs) in [
        ("peptide-count", vec!["--parquet", "input.feature.parquet"]),
        (
            "spectral-count",
            vec![
                "--psm",
                "input.psm.parquet",
                "--parquet",
                "input.feature.parquet",
            ],
        ),
    ] {
        let mut argv = vec!["mokume", "quantify", "features2proteins"];
        argv.extend(inputs);
        argv.extend([
            "--sdrf",
            "input.sdrf.tsv",
            "--output",
            "counts.csv",
            "--quant-method",
            method,
            "--irs",
            "--irs-reference-sample",
            "Pool",
        ]);
        let cli = Cli::parse_from(argv);
        let args = features_to_proteins_args(cli);
        let Err(error) = args.into_config() else {
            panic!("IRS was accepted for a count method");
        };
        assert!(error.to_string().contains("cannot apply IRS"), "{error}");
    }
}

#[test]
fn pibaq_uses_thirty_aa_and_method_specific_min_unique() {
    let cli = Cli::parse_from([
        "mokume",
        "quantify",
        "features2proteins",
        "--parquet",
        "input.parquet",
        "--output",
        "protein.csv",
        "--quant-method",
        "pibaq",
        "--fasta",
        "proteins.fasta",
    ]);
    let args = features_to_proteins_args(cli);
    let Ok(config) = args.into_config() else {
        panic!("expected a valid piBAQ config");
    };
    assert_eq!(config.pibaq.max_aa, 30);
    assert_eq!(config.filtering.min_unique_peptides, 0);

    let explicit = Cli::parse_from([
        "mokume",
        "quantify",
        "features2proteins",
        "--parquet",
        "input.parquet",
        "--output",
        "protein.csv",
        "--quant-method",
        "pibaq",
        "--fasta",
        "proteins.fasta",
        "--min-unique",
        "0",
    ]);
    let args = features_to_proteins_args(explicit);
    assert!(args.into_config().is_err());
}

fn features_to_proteins_args(cli: Cli) -> Box<Features2ProteinsArgs> {
    let Commands::Quantify(quantify) = cli.command else {
        panic!("expected quantify command");
    };
    let QuantifyCommands::Features2Proteins(args) = quantify.command else {
        panic!("expected features2proteins command");
    };
    args
}
