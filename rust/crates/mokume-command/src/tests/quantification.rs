use std::path::Path;

use clap::Parser;

use crate::{Cli, Commands};

#[test]
fn parses_true_spectral_count_psm_input() {
    let cli = Cli::parse_from([
        "mokume",
        "features2proteins",
        "--psm",
        "input.psm.parquet",
        "--sdrf",
        "input.sdrf.tsv",
        "--output",
        "counts.csv",
        "--quant-method",
        "spectral_count",
    ]);
    let Commands::Features2Proteins(args) = cli.command else {
        panic!("expected features2proteins command");
    };
    let Ok(config) = args.into_config() else {
        panic!("expected a valid spectral-count config");
    };
    assert_eq!(
        config.input.psm.as_deref(),
        Some(Path::new("input.psm.parquet"))
    );
    assert_eq!(
        config.quantification,
        mokume_core::QuantMethod::SpectralCount
    );
    assert_eq!(config.normalization.run_method, "none");
    assert_eq!(config.normalization.sample_method, "none");
}

#[test]
fn feature_input_uses_explicit_peptide_count_name() {
    let cli = Cli::parse_from([
        "mokume",
        "features2proteins",
        "--parquet",
        "input.feature.parquet",
        "--output",
        "counts.csv",
        "--quant-method",
        "peptide_count",
    ]);
    let Commands::Features2Proteins(args) = cli.command else {
        panic!("expected features2proteins command");
    };
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
    for (method, input_flag, input) in [
        ("peptide_count", "--parquet", "input.feature.parquet"),
        ("spectral_count", "--psm", "input.psm.parquet"),
    ] {
        let cli = Cli::parse_from([
            "mokume",
            "features2proteins",
            input_flag,
            input,
            "--sdrf",
            "input.sdrf.tsv",
            "--output",
            "counts.csv",
            "--quant-method",
            method,
            "--irs",
            "--irs-reference-samples",
            "Pool",
        ]);
        let Commands::Features2Proteins(args) = cli.command else {
            panic!("expected features2proteins command");
        };
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
    let Commands::Features2Proteins(args) = cli.command else {
        panic!("expected features2proteins command");
    };
    let Ok(config) = args.into_config() else {
        panic!("expected a valid piBAQ config");
    };
    assert_eq!(config.pibaq.max_aa, 30);
    assert_eq!(config.filtering.min_unique_peptides, 0);

    let explicit = Cli::parse_from([
        "mokume",
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
    let Commands::Features2Proteins(args) = explicit.command else {
        panic!("expected features2proteins command");
    };
    assert!(args.into_config().is_err());
}
