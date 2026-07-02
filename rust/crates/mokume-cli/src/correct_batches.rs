//! `correct-batches` command: parametric ComBat batch-effect correction for
//! long-format iBAQ TSV files.
//!
//! This mirrors `mokume.commands.batch_correct.run_batch_correction` in the
//! Python package. The flow:
//!   1. glob `folder` for `pattern`, concatenate the long-format TSV files
//!      (honoring the comment character and the separator);
//!   2. pivot to a protein (row) x sample (column) matrix of raw iBAQ, filling
//!      missing cells with 0.0 (matching `pivot_wider(..., fillna=True)`);
//!   3. validate sample IDs and derive integer batch labels from the prefix
//!      before the first '-' using first-seen (`pd.factorize`) order;
//!   4. run parametric ComBat (`combat_parametric` with default `ComBatParams`,
//!      equivalent to `apply_batch_correction(df, batch, kwargs={})`);
//!   5. merge the corrected values back onto the raw long rows by
//!      (sample, protein) and write the long-format TSV with the corrected
//!      column appended.
//!
//! When `--export-anndata` is set, an AnnData `.h5ad` file is written alongside
//! the TSV (the output path with its extension replaced by `.h5ad`), mirroring
//! `mokume.io.parquet.create_anndata` + `adata.write`. See [`export_anndata`].
//!
//! Determinism note: Python's output row order follows `glob.glob` (filesystem
//! order). Here the matched files are sorted by path and rows are emitted in
//! file/append order, which is deterministic. The corrected values are keyed by
//! (sample, protein), so they are independent of that ordering. Both the matrix
//! columns (samples) and rows (proteins) are sorted before ComBat, matching the
//! pandas `pivot_table` ordering that drives the Python batch labels.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use mokume_core::{MokumeError, Result};
use mokume_stats::batch::{combat_parametric, ComBatParams};

use crate::CorrectBatchesArgs;

/// A single parsed long-format input row, retaining every column so the output
/// reproduces the input plus the appended corrected column.
struct LongRow {
    values: Vec<String>,
}

/// The combined long-format table read from all matched files.
struct LongTable {
    headers: Vec<String>,
    rows: Vec<LongRow>,
    sample_index: usize,
    protein_index: usize,
    ibaq_index: usize,
}

/// Entry point for the `correct-batches` command.
pub fn run_correct_batches(args: &CorrectBatchesArgs) -> Result<()> {
    let separator = single_byte(&args.sep, "sep")?;
    let comment = optional_single_byte(&args.comment, "comment")?;

    let table = load_long_table(
        &args.folder,
        &args.pattern,
        separator,
        comment,
        &args.sample_id_column,
        &args.protein_id_column,
        &args.ibaq_raw_column,
    )?;

    // Sorted, unique proteins (rows) and samples (columns), matching the pandas
    // `pivot_table` ordering the Python code relies on.
    let mut proteins = BTreeSet::new();
    let mut samples = BTreeSet::new();
    for row in &table.rows {
        proteins.insert(cell(row, table.protein_index)?.to_owned());
        samples.insert(cell(row, table.sample_index)?.to_owned());
    }
    let proteins = proteins.into_iter().collect::<Vec<_>>();
    let samples = samples.into_iter().collect::<Vec<_>>();

    validate_sample_ids(&samples)?;
    let batch = batch_labels(&samples)?;
    ensure_batch_sizes(&batch)?;

    let matrix = build_matrix(&table, &proteins, &samples)?;
    let corrected = combat_parametric(&matrix, &batch, ComBatParams::default());

    // Key corrected values by (sample, protein) so the merge is independent of
    // row order.
    let mut corrected_by_key: HashMap<(&str, &str), f64> = HashMap::new();
    for (protein_idx, protein) in proteins.iter().enumerate() {
        for (sample_idx, sample) in samples.iter().enumerate() {
            corrected_by_key.insert(
                (sample.as_str(), protein.as_str()),
                corrected[protein_idx][sample_idx],
            );
        }
    }

    write_output(&args.output, &table, &args.ibaq_corrected_column, |row| {
        let sample = cell(row, table.sample_index)?;
        let protein = cell(row, table.protein_index)?;
        Ok(corrected_by_key.get(&(sample, protein)).copied())
    })?;

    if args.export_anndata {
        export_anndata(args, &table, &proteins, &samples, &matrix, &corrected)?;
    }

    Ok(())
}

/// Write the AnnData `.h5ad` export alongside the TSV, mirroring the Python
/// `create_anndata(df_ibaq, obs_col=SampleID, var_col=ProteinName,
/// value_col=Ibaq, layer_cols=[IbaqBec])` path.
///
/// The matrices supplied here are protein (row) x sample (column); AnnData wants
/// observation (sample) x variable (protein), so both are transposed. `X` holds
/// the raw iBAQ (missing cells already 0 from `build_matrix`). The `IbaqBec`
/// layer is built from the *merged long table*: a corrected value only lands in
/// a cell when that (sample, protein) pair was present in the input; cells that
/// were absent stay 0 (the Python layer pivot fills them with 0 rather than the
/// ComBat output, because the absent row never receives an `IbaqBec` value).
fn export_anndata(
    args: &CorrectBatchesArgs,
    table: &LongTable,
    proteins: &[String],
    samples: &[String],
    raw: &[Vec<f64>],
    corrected: &[Vec<f64>],
) -> Result<()> {
    // Which (protein, sample) pairs appeared in the input long table.
    let mut present: BTreeSet<(usize, usize)> = BTreeSet::new();
    let protein_pos = index_map(proteins);
    let sample_pos = index_map(samples);
    for row in &table.rows {
        let protein = cell(row, table.protein_index)?;
        let sample = cell(row, table.sample_index)?;
        if let (Some(&p), Some(&s)) = (protein_pos.get(protein), sample_pos.get(sample)) {
            present.insert((p, s));
        }
    }

    // Transpose protein x sample into sample (obs) x protein (var) layout.
    let n_obs = samples.len();
    let n_var = proteins.len();
    let mut x = vec![vec![0.0_f64; n_var]; n_obs];
    let mut layer = vec![vec![0.0_f64; n_var]; n_obs];
    for (protein_idx, _) in proteins.iter().enumerate() {
        for (sample_idx, _) in samples.iter().enumerate() {
            x[sample_idx][protein_idx] = raw[protein_idx][sample_idx];
            if present.contains(&(protein_idx, sample_idx)) {
                layer[sample_idx][protein_idx] = corrected[protein_idx][sample_idx];
            }
        }
    }

    let output_path = anndata_path(&args.output);
    let layers = vec![(args.ibaq_corrected_column.clone(), layer)];
    let export = crate::h5ad::AnnDataExport {
        obs_names: samples,
        obs_index_name: &args.sample_id_column,
        var_names: proteins,
        var_index_name: &args.protein_id_column,
        x: &x,
        layers: &layers,
    };
    crate::h5ad::write_h5ad(&output_path, &export)
}

/// Replace the output file's extension with `.h5ad`, matching Python's
/// `Path(output).with_suffix(".h5ad")`.
fn anndata_path(output: &Path) -> PathBuf {
    output.with_extension("h5ad")
}

/// Glob `folder/pattern`, read every matched TSV, and concatenate the rows.
/// Files are visited in sorted path order for determinism. Schemas must match
/// the first file (matching `combine_ibaq_tsv_files`).
fn load_long_table(
    folder: &Path,
    pattern: &str,
    separator: u8,
    comment: Option<u8>,
    sample_column: &str,
    protein_column: &str,
    ibaq_column: &str,
) -> Result<LongTable> {
    let paths = matched_files(folder, pattern)?;
    if paths.is_empty() {
        return Err(MokumeError::InvalidInput {
            message: format!(
                "No files found in the directory '{}' matching the pattern '{pattern}'.",
                folder.display()
            ),
        });
    }

    let mut headers: Option<Vec<String>> = None;
    let mut rows = Vec::new();

    for path in &paths {
        let (file_headers, file_rows) = read_long_file(path, separator, comment)?;
        match &headers {
            None => headers = Some(file_headers),
            Some(expected) => {
                if !same_columns(expected, &file_headers) {
                    return Err(MokumeError::InvalidInput {
                        message: format!(
                            "Schema mismatch in file '{}'. Expected columns: {expected:?}, got: {file_headers:?}",
                            path.display()
                        ),
                    });
                }
            }
        }
        rows.extend(file_rows.into_iter().map(|values| LongRow { values }));
    }

    let headers = headers.unwrap_or_default();
    let sample_index = column_index(&headers, sample_column)?;
    let protein_index = column_index(&headers, protein_column)?;
    let ibaq_index = column_index(&headers, ibaq_column)?;

    Ok(LongTable {
        headers,
        rows,
        sample_index,
        protein_index,
        ibaq_index,
    })
}

/// Read a single long-format TSV, returning `(headers, rows)`.
fn read_long_file(
    path: &Path,
    separator: u8,
    comment: Option<u8>,
) -> Result<(Vec<String>, Vec<Vec<String>>)> {
    let mut builder = csv::ReaderBuilder::new();
    builder.delimiter(separator).flexible(false);
    if let Some(comment) = comment {
        builder.comment(Some(comment));
    }
    let mut reader = builder
        .from_path(path)
        .map_err(|source| csv_error(path, source))?;

    let headers = reader
        .headers()
        .map_err(|source| csv_error(path, source))?
        .iter()
        .map(ToOwned::to_owned)
        .collect::<Vec<_>>();

    let mut rows = Vec::new();
    for record in reader.records() {
        let record = record.map_err(|source| csv_error(path, source))?;
        rows.push(record.iter().map(ToOwned::to_owned).collect::<Vec<_>>());
    }
    Ok((headers, rows))
}

/// Build the protein (row) x sample (column) matrix of raw iBAQ values, filling
/// missing cells with 0.0 (matching `pivot_wider(..., fillna=True)`). Duplicate
/// (protein, sample) combinations are rejected, matching `pivot_wider`.
fn build_matrix(
    table: &LongTable,
    proteins: &[String],
    samples: &[String],
) -> Result<Vec<Vec<f64>>> {
    let protein_pos = index_map(proteins);
    let sample_pos = index_map(samples);
    let mut matrix = vec![vec![0.0_f64; samples.len()]; proteins.len()];
    let mut seen = BTreeSet::new();

    for row in &table.rows {
        let protein = cell(row, table.protein_index)?;
        let sample = cell(row, table.sample_index)?;
        if !seen.insert((protein.to_owned(), sample.to_owned())) {
            return Err(MokumeError::InvalidInput {
                message: format!(
                    "Found duplicate combination of protein '{protein}' and sample '{sample}'."
                ),
            });
        }
        let raw = cell(row, table.ibaq_index)?;
        let value = parse_value(raw, table.ibaq_index)?;
        let (Some(&r), Some(&c)) = (protein_pos.get(protein), sample_pos.get(sample)) else {
            // Unreachable: proteins/samples were derived from these same rows.
            continue;
        };
        matrix[r][c] = value;
    }
    Ok(matrix)
}

/// Write the long-format output: every original column followed by the corrected
/// column. `lookup` returns the corrected value for a row (None leaves the cell
/// empty, matching a left-merge miss).
fn write_output(
    output: &Path,
    table: &LongTable,
    corrected_column: &str,
    lookup: impl Fn(&LongRow) -> Result<Option<f64>>,
) -> Result<()> {
    let file = File::create(output).map_err(|source| MokumeError::Io {
        path: output.to_path_buf(),
        source,
    })?;
    let mut writer = BufWriter::new(file);

    let mut header_line = table.headers.join("\t");
    header_line.push('\t');
    header_line.push_str(corrected_column);
    writeln!(writer, "{header_line}").map_err(|source| MokumeError::Io {
        path: output.to_path_buf(),
        source,
    })?;

    for row in &table.rows {
        let corrected = lookup(row)?;
        let mut line = row.values.join("\t");
        line.push('\t');
        if let Some(value) = corrected {
            line.push_str(&format_value(value));
        }
        writeln!(writer, "{line}").map_err(|source| MokumeError::Io {
            path: output.to_path_buf(),
            source,
        })?;
    }
    writer.flush().map_err(|source| MokumeError::Io {
        path: output.to_path_buf(),
        source,
    })?;
    Ok(())
}

/// Derive integer batch labels from sample names: the prefix before the first
/// '-' is the batch ID; labels follow first-seen order (`pd.factorize`).
fn batch_labels(samples: &[String]) -> Result<Vec<usize>> {
    let mut order: Vec<String> = Vec::new();
    let mut seen: HashMap<String, usize> = HashMap::new();
    let mut labels = Vec::with_capacity(samples.len());
    for sample in samples {
        let prefix = batch_prefix(sample)?;
        let label = match seen.get(&prefix) {
            Some(&label) => label,
            None => {
                let label = order.len();
                order.push(prefix.clone());
                seen.insert(prefix, label);
                label
            }
        };
        labels.push(label);
    }
    Ok(labels)
}

/// Extract and validate the batch prefix (text before the first '-').
fn batch_prefix(sample: &str) -> Result<String> {
    let prefix = sample.split('-').next().unwrap_or("");
    if prefix.is_empty() {
        return Err(MokumeError::InvalidInput {
            message: format!("Invalid sample name format: {sample}"),
        });
    }
    if !prefix.chars().all(|c| c.is_ascii_alphanumeric()) {
        return Err(MokumeError::InvalidInput {
            message: format!("Invalid batch ID format: {prefix}"),
        });
    }
    Ok(prefix.to_owned())
}

/// Validate sample IDs against `SAMPLE_ID_REGEX`:
/// `^[A-Za-z0-9]+(?:-[A-Za-z0-9]+)*$` (alphanumeric segments joined by '-').
fn validate_sample_ids(samples: &[String]) -> Result<()> {
    let invalid = samples
        .iter()
        .filter(|sample| !is_valid_sample_id(sample))
        .cloned()
        .collect::<Vec<_>>();
    if invalid.is_empty() {
        Ok(())
    } else {
        Err(MokumeError::InvalidInput {
            message: format!("Invalid sample IDs found in the data: {invalid:?}"),
        })
    }
}

/// True when `sample` matches `^[A-Za-z0-9]+(?:-[A-Za-z0-9]+)*$`.
fn is_valid_sample_id(sample: &str) -> bool {
    if sample.is_empty() {
        return false;
    }
    sample
        .split('-')
        .all(|segment| !segment.is_empty() && segment.chars().all(|c| c.is_ascii_alphanumeric()))
}

/// Every batch must contain at least two samples (matching
/// `TooFewSamplesInBatch` in the Python `apply_batch_correction`).
fn ensure_batch_sizes(batch: &[usize]) -> Result<()> {
    let mut counts: BTreeMap<usize, usize> = BTreeMap::new();
    for &label in batch {
        *counts.entry(label).or_insert(0) += 1;
    }
    let short = counts
        .iter()
        .filter(|&(_, &count)| count < 2)
        .map(|(&label, _)| label)
        .collect::<Vec<_>>();
    if short.is_empty() {
        Ok(())
    } else {
        Err(MokumeError::InvalidInput {
            message: format!(
                "Batches must contain at least two samples, the following batch factors did not: {short:?}"
            ),
        })
    }
}

/// Return the matched files for `folder/pattern`, sorted by path. Supports the
/// shell-style wildcards `*` and `?` used by mokume's default `*ibaq.tsv`.
fn matched_files(folder: &Path, pattern: &str) -> Result<Vec<PathBuf>> {
    let entries = std::fs::read_dir(folder).map_err(|source| MokumeError::Io {
        path: folder.to_path_buf(),
        source,
    })?;
    let mut matches = Vec::new();
    for entry in entries {
        let entry = entry.map_err(|source| MokumeError::Io {
            path: folder.to_path_buf(),
            source,
        })?;
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        if let Some(name) = path.file_name().and_then(|name| name.to_str()) {
            if glob_match(pattern, name) {
                matches.push(path);
            }
        }
    }
    matches.sort();
    Ok(matches)
}

/// Minimal glob matcher for `*` (any run, including empty) and `?` (single
/// character). Sufficient for the `*ibaq.tsv` style patterns mokume uses.
fn glob_match(pattern: &str, name: &str) -> bool {
    let pattern = pattern.chars().collect::<Vec<_>>();
    let name = name.chars().collect::<Vec<_>>();
    glob_match_at(&pattern, &name)
}

fn glob_match_at(pattern: &[char], name: &[char]) -> bool {
    match pattern.first() {
        None => name.is_empty(),
        Some('*') => {
            // Match zero or more characters.
            glob_match_at(&pattern[1..], name)
                || (!name.is_empty() && glob_match_at(pattern, &name[1..]))
        }
        Some('?') => !name.is_empty() && glob_match_at(&pattern[1..], &name[1..]),
        Some(&c) => name.first() == Some(&c) && glob_match_at(&pattern[1..], &name[1..]),
    }
}

/// Map values to their positions for fast lookup.
fn index_map(values: &[String]) -> HashMap<&str, usize> {
    values
        .iter()
        .enumerate()
        .map(|(index, value)| (value.as_str(), index))
        .collect()
}

/// Two header sets are equal when they contain the same column names (order is
/// irrelevant, matching the Python set comparison in `combine_ibaq_tsv_files`).
fn same_columns(left: &[String], right: &[String]) -> bool {
    let left: BTreeSet<&String> = left.iter().collect();
    let right: BTreeSet<&String> = right.iter().collect();
    left == right
}

fn column_index(headers: &[String], column: &str) -> Result<usize> {
    headers
        .iter()
        .position(|header| header == column)
        .ok_or_else(|| MokumeError::InvalidInput {
            message: format!("Column '{column}' not found in the input files."),
        })
}

fn cell(row: &LongRow, index: usize) -> Result<&str> {
    row.values
        .get(index)
        .map(String::as_str)
        .ok_or_else(|| MokumeError::InvalidInput {
            message: format!("Row is missing column {index}."),
        })
}

fn parse_value(raw: &str, index: usize) -> Result<f64> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Ok(0.0);
    }
    trimmed
        .parse::<f64>()
        .map_err(|_| MokumeError::InvalidInput {
            message: format!("Column {index} value '{raw}' is not numeric."),
        })
}

/// Format a corrected value. `{}` on `f64` gives Rust's shortest round-trip
/// representation, matching Python's `repr(float)` for the precision the golden
/// test compares (relative 1e-6).
fn format_value(value: f64) -> String {
    format!("{value}")
}

fn single_byte(value: &str, name: &str) -> Result<u8> {
    let bytes = value.as_bytes();
    if bytes.len() == 1 {
        Ok(bytes[0])
    } else {
        Err(MokumeError::InvalidInput {
            message: format!("--{name} must be a single byte, got '{value}'."),
        })
    }
}

fn optional_single_byte(value: &str, name: &str) -> Result<Option<u8>> {
    if value.is_empty() {
        Ok(None)
    } else {
        single_byte(value, name).map(Some)
    }
}

fn csv_error(path: &Path, source: csv::Error) -> MokumeError {
    MokumeError::InvalidInput {
        message: format!("Error reading file '{}': {source}", path.display()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::error::Error;
    use std::time::{SystemTime, UNIX_EPOCH};

    // `super::*` brings the crate's single-parameter `Result` alias into scope,
    // so the test helpers spell out the two-parameter standard `Result`.
    type TestResult<T> = std::result::Result<T, Box<dyn Error>>;

    fn temp_dir(tag: &str) -> TestResult<PathBuf> {
        let nanos = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
        Ok(tempfile::Builder::new()
            .prefix(&format!("mokume-correct-batches-{tag}-{nanos}-"))
            .tempdir()?
            .keep())
    }

    fn write_file(dir: &Path, name: &str, contents: &str) -> TestResult<PathBuf> {
        let path = dir.join(name);
        let mut file = File::create(&path)?;
        file.write_all(contents.as_bytes())?;
        Ok(path)
    }

    #[test]
    fn glob_match_handles_wildcards() {
        assert!(glob_match("*ibaq.tsv", "batchA_ibaq.tsv"));
        assert!(glob_match("*ibaq.tsv", "ibaq.tsv"));
        assert!(!glob_match("*ibaq.tsv", "proteins.tsv"));
        assert!(glob_match("sample?.tsv", "sample1.tsv"));
        assert!(!glob_match("sample?.tsv", "sample12.tsv"));
        assert!(glob_match("*", "anything"));
    }

    #[test]
    fn batch_labels_use_first_seen_order() -> TestResult<()> {
        let samples = vec![
            "B2-s1".to_string(),
            "B1-s1".to_string(),
            "B2-s2".to_string(),
            "B1-s2".to_string(),
        ];
        assert_eq!(batch_labels(&samples)?, vec![0, 1, 0, 1]);
        Ok(())
    }

    #[test]
    fn invalid_batch_prefix_is_rejected() {
        assert!(batch_prefix("good-s1").is_ok());
        assert!(batch_prefix("-s1").is_err());
    }

    #[test]
    fn sample_id_validation_matches_regex() {
        assert!(is_valid_sample_id("B1-s1"));
        assert!(is_valid_sample_id("PXD000001"));
        assert!(!is_valid_sample_id("B1--s1"));
        assert!(!is_valid_sample_id("B1-"));
        assert!(!is_valid_sample_id("B1_s1"));
        assert!(!is_valid_sample_id("bad sample"));
    }

    fn args_for(folder: &Path, output: &Path) -> CorrectBatchesArgs {
        CorrectBatchesArgs {
            folder: folder.to_path_buf(),
            pattern: "*ibaq.tsv".to_string(),
            comment: "#".to_string(),
            sep: "\t".to_string(),
            output: output.to_path_buf(),
            sample_id_column: "SampleID".to_string(),
            protein_id_column: "ProteinName".to_string(),
            ibaq_raw_column: "Ibaq".to_string(),
            ibaq_corrected_column: "IbaqBec".to_string(),
            export_anndata: false,
        }
    }

    const BATCH_A: &str = "# synthetic ibaq fixture for correct-batches golden test\n\
ProteinName\tSampleID\tCondition\tIbaq\n\
P1\tB1-s1\tctrl\t10.0\n\
P2\tB1-s1\tctrl\t5.0\n\
P3\tB1-s1\tctrl\t1.0\n\
P4\tB1-s1\tctrl\t50.0\n\
P5\tB1-s1\tctrl\t3.0\n\
P1\tB1-s2\tcase\t11.0\n\
P2\tB1-s2\tcase\t6.0\n\
P3\tB1-s2\tcase\t2.0\n\
P4\tB1-s2\tcase\t52.0\n\
P5\tB1-s2\tcase\t3.5\n\
P1\tB1-s3\tcase\t9.5\n\
P2\tB1-s3\tcase\t4.0\n\
P3\tB1-s3\tcase\t1.5\n\
P4\tB1-s3\tcase\t48.0\n";

    const BATCH_B: &str = "# synthetic ibaq fixture for correct-batches golden test\n\
ProteinName\tSampleID\tCondition\tIbaq\n\
P1\tB2-s1\tctrl\t20.0\n\
P2\tB2-s1\tctrl\t8.0\n\
P3\tB2-s1\tctrl\t3.0\n\
P4\tB2-s1\tctrl\t30.0\n\
P5\tB2-s1\tctrl\t6.0\n\
P1\tB2-s2\tcase\t21.0\n\
P2\tB2-s2\tcase\t7.5\n\
P3\tB2-s2\tcase\t2.5\n\
P4\tB2-s2\tcase\t31.0\n\
P5\tB2-s2\tcase\t5.5\n\
P1\tB2-s3\tcase\t19.0\n\
P2\tB2-s3\tcase\t9.0\n\
P3\tB2-s3\tcase\t4.0\n\
P4\tB2-s3\tcase\t29.0\n\
P5\tB2-s3\tcase\t6.5\n";

    // Expected IbaqBec values keyed by (SampleID, ProteinName), captured from the
    // Python oracle:
    //   conda run -n Bigbio python -m mokume.mokume_cli correct-batches \
    //     --folder <fixture-dir> --pattern "*ibaq.tsv" --output <out.tsv>
    // inmoose 0.9.1 / pandas. B1-s3/P5 is intentionally absent in the input, so
    // it never appears in the output (left-merge miss).
    const EXPECTED: &[(&str, &str, f64)] = &[
        ("B1-s1", "P1", 14.86027101136503),
        ("B1-s1", "P2", 6.568250615964994),
        ("B1-s1", "P3", 1.8528067846912328),
        ("B1-s1", "P4", 40.18479019681545),
        ("B1-s1", "P5", 4.791576963794292),
        ("B1-s2", "P1", 15.799571789576207),
        ("B1-s2", "P2", 7.472688758460274),
        ("B1-s2", "P3", 2.802714985157337),
        ("B1-s2", "P4", 41.94917557277646),
        ("B1-s2", "P5", 5.225650064373009),
        ("B1-s3", "P1", 14.390620622259444),
        ("B1-s3", "P2", 5.663812473469715),
        ("B1-s3", "P3", 2.327760884924285),
        ("B1-s3", "P4", 38.42040482085444),
        ("B2-s1", "P1", 15.14370361017618),
        ("B2-s1", "P2", 6.405668897416727),
        ("B2-s1", "P3", 2.169016234610666),
        ("B2-s1", "P4", 39.86488123350397),
        ("B2-s1", "P5", 4.093616827669412),
        ("B2-s2", "P1", 16.181696236203024),
        ("B2-s2", "P2", 5.836487008977049),
        ("B2-s2", "P3", 1.6604208788502293),
        ("B2-s2", "P4", 41.08220056732911),
        ("B2-s2", "P5", 3.4445768562464525),
        ("B2-s3", "P1", 14.105710984149336),
        ("B2-s3", "P2", 7.544032674296082),
        ("B2-s3", "P3", 3.186206946131539),
        ("B2-s3", "P4", 38.647561899678834),
        ("B2-s3", "P5", 4.742656799092371),
    ];

    /// Golden test: the Rust command reproduces the Python oracle's IbaqBec
    /// values to relative 1e-6 on a synthetic 2-batch / 6-sample / 5-protein
    /// dataset (one cell intentionally missing).
    #[test]
    fn correct_batches_matches_python_oracle() -> TestResult<()> {
        let dir = temp_dir("oracle")?;
        write_file(&dir, "batchA_ibaq.tsv", BATCH_A)?;
        write_file(&dir, "batchB_ibaq.tsv", BATCH_B)?;
        let output = dir.join("corrected.tsv");

        run_correct_batches(&args_for(&dir, &output))?;

        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(&output)?;
        let headers = reader
            .headers()?
            .iter()
            .map(ToOwned::to_owned)
            .collect::<Vec<_>>();
        let sample_col = column_index(&headers, "SampleID")?;
        let protein_col = column_index(&headers, "ProteinName")?;
        let bec_col = column_index(&headers, "IbaqBec")?;

        let mut actual: HashMap<(String, String), f64> = HashMap::new();
        for record in reader.records() {
            let record = record?;
            let sample = field(&record, sample_col)?.to_string();
            let protein = field(&record, protein_col)?.to_string();
            let value = field(&record, bec_col)?.parse::<f64>()?;
            actual.insert((sample, protein), value);
        }

        assert_eq!(
            actual.len(),
            EXPECTED.len(),
            "row count must match the oracle"
        );
        for &(sample, protein, expected) in EXPECTED {
            let Some(&got) = actual.get(&(sample.to_string(), protein.to_string())) else {
                panic!("missing corrected value for ({sample}, {protein})");
            };
            let tolerance = expected.abs() * 1e-6;
            assert!(
                (got - expected).abs() <= tolerance,
                "({sample}, {protein}): got {got}, expected {expected}"
            );
        }
        Ok(())
    }

    fn field(record: &csv::StringRecord, index: usize) -> TestResult<&str> {
        record
            .get(index)
            .ok_or_else(|| format!("record missing column {index}").into())
    }

    #[test]
    fn export_anndata_writes_h5ad_file() -> TestResult<()> {
        let dir = temp_dir("anndata")?;
        write_file(&dir, "batchA_ibaq.tsv", BATCH_A)?;
        write_file(&dir, "batchB_ibaq.tsv", BATCH_B)?;
        let output = dir.join("corrected.tsv");
        let mut args = args_for(&dir, &output);
        args.export_anndata = true;

        run_correct_batches(&args)?;

        // The `.h5ad` lands next to the TSV with the suffix replaced.
        let h5ad = dir.join("corrected.h5ad");
        assert!(
            h5ad.is_file(),
            "expected AnnData file at {}",
            h5ad.display()
        );

        // Read X / IbaqBec back with the hdf5 crate and confirm the corrected
        // layer matches the oracle (samples x proteins) and that the absent
        // input cell B1-s3/P5 carries 0 in both X and the layer.
        let file = hdf5_metno::File::open(&h5ad)?;
        let obs: Vec<hdf5_metno::types::VarLenUnicode> =
            file.group("obs")?.dataset("SampleID")?.read_raw()?;
        let obs: Vec<String> = obs.iter().map(|s| s.as_str().to_string()).collect();
        let var: Vec<hdf5_metno::types::VarLenUnicode> =
            file.group("var")?.dataset("ProteinName")?.read_raw()?;
        let var: Vec<String> = var.iter().map(|s| s.as_str().to_string()).collect();

        let x: ndarray::Array2<f64> = file.dataset("X")?.read_2d()?;
        let layer: ndarray::Array2<f64> = file.dataset("layers/IbaqBec")?.read_2d()?;
        assert_eq!(x.shape(), &[obs.len(), var.len()]);
        assert_eq!(layer.shape(), &[obs.len(), var.len()]);

        let obs_idx = |name: &str| obs.iter().position(|s| s == name);
        let var_idx = |name: &str| var.iter().position(|s| s == name);

        // Corrected layer values match the Python oracle where present.
        for &(sample, protein, expected) in EXPECTED {
            let (Some(i), Some(j)) = (obs_idx(sample), var_idx(protein)) else {
                panic!("missing ({sample}, {protein}) in AnnData index");
            };
            let got = layer[[i, j]];
            let tolerance = expected.abs() * 1e-6;
            assert!(
                (got - expected).abs() <= tolerance,
                "layer ({sample}, {protein}): got {got}, expected {expected}"
            );
        }

        // The absent input cell keeps 0 in both X and the layer (it never
        // received an IbaqBec value in the merged long table).
        let (Some(i), Some(j)) = (obs_idx("B1-s3"), var_idx("P5")) else {
            panic!("missing B1-s3/P5 in AnnData index");
        };
        assert_eq!(x[[i, j]], 0.0, "absent X cell must be 0");
        assert_eq!(layer[[i, j]], 0.0, "absent layer cell must be 0");

        // Raw X reflects the input iBAQ values.
        let (Some(i), Some(j)) = (obs_idx("B1-s1"), var_idx("P4")) else {
            panic!("missing B1-s1/P4 in AnnData index");
        };
        assert_eq!(x[[i, j]], 50.0);
        Ok(())
    }
}
