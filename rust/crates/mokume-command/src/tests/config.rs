use std::path::Path;

use clap::Parser;
use mokume_core::{MokumeError, PreprocessingFilterConfig};

use crate::{Cli, Commands};

const PYTHON_FEATURES_TO_PROTEINS_ARGS: &[&str] = &[
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
];

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
    let cli = Cli::parse_from(PYTHON_FEATURES_TO_PROTEINS_ARGS.iter().copied());

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
    assert_eq!(config.batch.ref_batch.as_deref(), Some("2"));
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
fn impute_method_enables_imputation_and_constant_maps_to_zero() {
    let cli = Cli::parse_from([
        "mokume",
        "features2proteins",
        "--parquet",
        "input.parquet",
        "--output",
        "protein.csv",
        "--impute-method",
        "constant",
    ]);
    let Commands::Features2Proteins(args) = cli.command else {
        panic!("expected features2proteins command");
    };
    let Ok(config) = args.into_config() else {
        panic!("expected a valid imputation config");
    };
    assert!(config.imputation.enabled);
    assert_eq!(config.imputation.method, "zero");
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
fn filter_config_rejects_unknown_keys() {
    for yaml in [
        "unknown_option: true\n",
        "intensity:\n  min_intensitty: 5\n",
    ] {
        let parsed = serde_yaml::from_str::<PreprocessingFilterConfig>(yaml);
        assert!(parsed.is_err(), "unknown filter key was accepted: {yaml}");
    }
}
