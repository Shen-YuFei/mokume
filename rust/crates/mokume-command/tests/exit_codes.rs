use std::{
    fs::{create_dir_all, write},
    path::{Path, PathBuf},
    time::{SystemTime, UNIX_EPOCH},
};

use mokume_command::run_from_args;

#[test]
fn features2proteins_missing_input_returns_error() {
    let error = match run(&[
        "features2proteins",
        "--parquet",
        "definitely-missing.feature.parquet",
        "--output",
        "protein.csv",
    ]) {
        Ok(()) => panic!("command must fail for a missing input"),
        Err(error) => error,
    };
    let message = error.to_string();
    assert!(
        message.contains("input file does not exist: definitely-missing.feature.parquet"),
        "unexpected error: {message}"
    );
}

#[test]
fn peptides2protein_command_writes_protein_tsv() -> Result<(), Box<dyn std::error::Error>> {
    let root = temp_root()?;
    create_dir_all(&root)?;
    let peptides = root.join("peptides.csv");
    let output = root.join("protein.tsv");
    write(
        &peptides,
        "ProteinName,PeptideCanonical,SampleID,Condition,NormIntensity\n\
P1,PEPTIDEAK,S1,A,100.0\n\
P1,APEPTIDECK,S1,A,300.0\n\
P2,ALYAAEK,S1,A,500.0\n",
    )?;

    run(&[
        "peptides2protein",
        "--method",
        "sum",
        "--peptides",
        path_str(&peptides)?,
        "--output",
        path_str(&output)?,
    ])?;
    let written = std::fs::read_to_string(&output)?;
    let header = written.lines().next().unwrap_or_default();
    assert_eq!(header, "ProteinName\tSampleID\tIntensity\tCondition");
    assert!(
        written.lines().any(|line| line.starts_with("P1\tS1\t400")),
        "P1/S1 sum must be 400:\n{written}"
    );
    Ok(())
}

#[test]
fn correct_batches_command_writes_corrected_tsv() -> Result<(), Box<dyn std::error::Error>> {
    let root = temp_root()?;
    let input = root.join("input");
    create_dir_all(&input)?;
    write(
        input.join("batchA_pibaq.tsv"),
        "# header comment\nProteinName\tSampleID\tPiBAQ\n\
P1\tB1-s1\t10.0\nP2\tB1-s1\t5.0\nP1\tB1-s2\t11.0\nP2\tB1-s2\t6.0\n",
    )?;
    write(
        input.join("batchB_pibaq.tsv"),
        "# header comment\nProteinName\tSampleID\tPiBAQ\n\
P1\tB2-s1\t20.0\nP2\tB2-s1\t8.0\nP1\tB2-s2\t21.0\nP2\tB2-s2\t7.5\n",
    )?;
    let output = root.join("corrected.tsv");

    run(&[
        "correct-batches",
        "--folder",
        path_str(&input)?,
        "--output",
        path_str(&output)?,
    ])?;
    let written = std::fs::read_to_string(&output)?;
    let header = written.lines().next().unwrap_or_default();
    assert!(
        header.ends_with("PiBAQBec"),
        "output header must end with the corrected column:\n{header}"
    );
    assert_eq!(
        written.lines().count(),
        9,
        "output must contain a header and 8 data rows:\n{written}"
    );
    Ok(())
}

#[test]
fn correct_batches_rejects_input_output_collision() -> Result<(), Box<dyn std::error::Error>> {
    let root = temp_root()?;
    let input = root.join("input");
    create_dir_all(&input)?;
    let input_file = input.join("batchA_pibaq.tsv");
    write(
        &input_file,
        "ProteinName\tSampleID\tPiBAQ\n\
P1\tB1-s1\t10.0\nP1\tB1-s2\t11.0\n",
    )?;
    write(
        input.join("batchB_pibaq.tsv"),
        "ProteinName\tSampleID\tPiBAQ\n\
P1\tB2-s1\t20.0\nP1\tB2-s2\t21.0\n",
    )?;
    let before = std::fs::read(&input_file)?;

    let error = match run(&[
        "correct-batches",
        "--folder",
        path_str(&input)?,
        "--output",
        path_str(&input_file)?,
    ]) {
        Ok(()) => panic!("an output/input collision must fail"),
        Err(error) => error,
    };
    assert!(
        error.to_string().contains("also an input file"),
        "unexpected error: {error}"
    );
    assert_eq!(std::fs::read(&input_file)?, before);
    Ok(())
}

#[test]
fn correct_batches_export_anndata_writes_h5ad() -> Result<(), Box<dyn std::error::Error>> {
    let root = temp_root()?;
    let input = root.join("input");
    create_dir_all(&input)?;
    // Two batches (B1, B2) of two samples each, satisfying ComBat's minimum.
    write(
        input.join("batchA_pibaq.tsv"),
        "ProteinName\tSampleID\tPiBAQ\n\
P1\tB1-s1\t10.0\nP2\tB1-s1\t5.0\nP3\tB1-s1\t1.0\n\
P1\tB1-s2\t11.0\nP2\tB1-s2\t6.0\nP3\tB1-s2\t2.0\n",
    )?;
    write(
        input.join("batchB_pibaq.tsv"),
        "ProteinName\tSampleID\tPiBAQ\n\
P1\tB2-s1\t20.0\nP2\tB2-s1\t8.0\nP3\tB2-s1\t3.0\n\
P1\tB2-s2\t21.0\nP2\tB2-s2\t7.5\nP3\tB2-s2\t2.5\n",
    )?;
    let output = root.join("corrected.tsv");

    run(&[
        "correct-batches",
        "--folder",
        path_str(&input)?,
        "--output",
        path_str(&output)?,
        "--export_anndata",
    ])?;

    // The corrected TSV and the AnnData `.h5ad` are both written, the latter at
    // the output path with its extension replaced by `.h5ad`.
    assert!(output.is_file(), "expected corrected TSV at {output:?}");
    let h5ad = root.join("corrected.h5ad");
    assert!(h5ad.is_file(), "expected AnnData file at {h5ad:?}");

    // An `.h5ad` is an HDF5 container; its first eight bytes are the HDF5
    // signature. Checking it confirms a real HDF5 file was produced.
    let bytes = std::fs::read(&h5ad)?;
    assert!(
        bytes.starts_with(&[0x89, b'H', b'D', b'F', b'\r', b'\n', 0x1a, b'\n']),
        "output is not an HDF5 file"
    );
    Ok(())
}

#[test]
fn features2peptides_generates_filter_config() -> Result<(), Box<dyn std::error::Error>> {
    let root = temp_root()?;
    create_dir_all(&root)?;
    let yaml = root.join("filters.yaml");
    let json = root.join("filters.json");

    for output_path in [&yaml, &json] {
        run(&[
            "features2peptides",
            "--generate-filter-config",
            path_str(output_path)?,
        ])?;
    }

    let yaml = std::fs::read_to_string(yaml)?;
    assert!(
        yaml.contains("name: example_config"),
        "yaml config is missing the example name:\n{yaml}"
    );
    assert!(
        yaml.contains("min_unique_peptides: 2"),
        "yaml config is missing protein min_unique_peptides:\n{yaml}"
    );
    assert!(
        yaml.contains("max_missing_rate: 1.0"),
        "yaml config is missing run-QC max_missing_rate:\n{yaml}"
    );

    let json = std::fs::read_to_string(json)?;
    assert!(
        json.contains(r#""name": "example_config""#),
        "json config is missing the example name:\n{json}"
    );
    assert!(
        json.contains(r#""min_unique_peptides": 2"#),
        "json config is missing protein min_unique_peptides:\n{json}"
    );
    assert!(
        json.contains(r#""max_missing_rate": 1.0"#),
        "json config is missing run-QC max_missing_rate:\n{json}"
    );
    Ok(())
}

#[test]
fn unimplemented_features2proteins_options_return_stable_error(
) -> Result<(), Box<dyn std::error::Error>> {
    let root = temp_root()?;
    create_dir_all(&root)?;
    let parquet = root.join("input.parquet");
    let sdrf = root.join("input.sdrf.tsv");
    write(&parquet, [])?;
    write(
        &sdrf,
        "source name\tcomment[data file]\nSample1\trun1.raw\n",
    )?;
    let parquet = path_str(&parquet)?;
    let sdrf = path_str(&sdrf)?;

    for case in [
        FeatureToProteinsCase {
            // Cell-based methods (default sum) now carry dataset-level
            // normalization into the export; non-cell-based methods (maxlfq) do
            // not and stay rejected.
            args: &[
                "--quant-method",
                "maxlfq",
                "--sample-normalization",
                "quantile",
                "--export-peptides",
                "peptides.csv",
            ],
            stage: "dataset-normalization-export-peptides",
        },
        // The `--plot-*` / `--interactive-report` / `--report-output` flags were
        // removed from the Rust command interface (plotting / reports moved to
        // the Python periphery), so they are no longer exercised here.
        FeatureToProteinsCase {
            args: &["--quant-method", "sum", "--export-ions", "ions.csv"],
            stage: "export-ions",
        },
    ] {
        let mut args = vec![
            "features2proteins",
            "--parquet",
            parquet,
            "--sdrf",
            sdrf,
            "--output",
            "protein.csv",
        ];
        args.extend_from_slice(case.args);
        let error = match run(&args) {
            Ok(()) => panic!("unimplemented options must fail"),
            Err(error) => error,
        };
        let message = error.to_string();
        let expected = format!("stage `{}` is not implemented yet", case.stage);
        assert!(
            message.contains(&expected),
            "unexpected error for {:?}: {message}",
            case.args
        );
    }
    let error = match run(&[
        "features2proteins",
        "--parquet",
        parquet,
        "--sdrf",
        sdrf,
        "--output",
        "protein.csv",
        "--quant-method",
        "directlfq",
        "--export-peptides",
        "peptides.csv",
    ]) {
        Ok(()) => panic!("DirectLFQ peptide export must fail"),
        Err(error) => error,
    };
    assert!(error
        .to_string()
        .contains("export-peptides is not supported by directlfq"));
    Ok(())
}

struct FeatureToProteinsCase<'a> {
    args: &'a [&'a str],
    stage: &'a str,
}

fn run(args: &[&str]) -> mokume_core::Result<()> {
    let mut argv = Vec::with_capacity(args.len() + 1);
    argv.push("mokume");
    argv.extend_from_slice(args);
    run_from_args(argv)
}

fn temp_root() -> Result<PathBuf, Box<dyn std::error::Error>> {
    let timestamp = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
    Ok(tempfile::Builder::new()
        .prefix(&format!("mokume-command-test-{timestamp}-"))
        .tempdir()?
        .keep())
}

fn path_str(path: &Path) -> Result<&str, Box<dyn std::error::Error>> {
    path.to_str()
        .ok_or_else(|| format!("path is not UTF-8: {}", path.display()).into())
}
