//! Batch-label detection.
//!
//! Port of `mokume.model.batch_correction.BatchDetectionMethod` and
//! `mokume.postprocessing.batch_correction.detect_batches`
//! (`batch_correction.py:116`). Produces 0-indexed batch labels for a list of
//! samples under one of five detection methods, mirroring Python's
//! `pd.factorize` first-occurrence encoding.
//!
//! The methods:
//!   - `sample_prefix`: batch is the sample-name prefix (text before the first
//!     `-`, else the `[A-Za-z]+[0-9]+` run immediately before `_`, else the whole
//!     name);
//!   - `run`: batch is the run name looked up in `run_info` (errors when
//!     `run_info` is absent, matching Python's `ValueError`);
//!   - `fraction` / `techreplicate`: batch is the fraction / tech-replicate id
//!     looked up in `run_info`, defaulting missing samples to `"1"`; with no
//!     `run_info` both fall back to `sample_prefix` (Python logs a warning and
//!     recurses);
//!   - `column`: batch is an explicit per-sample value (errors when the value
//!     list is absent or length-mismatched, matching Python).

use std::collections::HashMap;

/// Batch-detection method, mirroring `BatchDetectionMethod`. The string forms
/// match the Python enum *values* (`sample_prefix`/`run`/`fraction`/
/// `techreplicate`/`column`), which [`BatchDetectionMethod::parse_name`]
/// normalises case-insensitively with `-`/space mapped to `_`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BatchDetectionMethod {
    SamplePrefix,
    RunName,
    Fraction,
    TechReplicate,
    ExplicitColumn,
}

impl BatchDetectionMethod {
    /// The canonical lowercase value, equal to the Python enum value.
    pub fn as_value(self) -> &'static str {
        match self {
            BatchDetectionMethod::SamplePrefix => "sample_prefix",
            BatchDetectionMethod::RunName => "run",
            BatchDetectionMethod::Fraction => "fraction",
            BatchDetectionMethod::TechReplicate => "techreplicate",
            BatchDetectionMethod::ExplicitColumn => "column",
        }
    }

    /// Parse a method name, mirroring `BatchDetectionMethod.from_str`: lowercase,
    /// then `-` and space become `_`, then match an enum value. Returns `None`
    /// for unknown names (Python raises `ValueError`; the pipeline caller maps an
    /// unknown name to `sample_prefix` with a warning, which it does on `None`).
    /// Named `parse_name` rather than `from_str` to avoid shadowing the standard
    /// `FromStr::from_str` trait method.
    pub fn parse_name(name: &str) -> Option<Self> {
        let normalized = name.to_ascii_lowercase().replace(['-', ' '], "_");
        match normalized.as_str() {
            "sample_prefix" => Some(BatchDetectionMethod::SamplePrefix),
            "run" => Some(BatchDetectionMethod::RunName),
            "fraction" => Some(BatchDetectionMethod::Fraction),
            "techreplicate" => Some(BatchDetectionMethod::TechReplicate),
            "column" => Some(BatchDetectionMethod::ExplicitColumn),
            _ => None,
        }
    }
}

/// Error returned when a method's required input is missing or inconsistent,
/// mirroring the `ValueError`s `detect_batches`/`_detect_explicit_batches` raise.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DetectBatchesError {
    /// `run` method requested without a `run_info` mapping.
    RunInfoRequired,
    /// `column` method requested without explicit batch values.
    ColumnValuesRequired,
    /// `column` value count does not match the sample count.
    ColumnLengthMismatch { values: usize, samples: usize },
}

/// Extract the batch-like prefix of a sample id, mirroring `_sample_prefix`
/// (`batch_correction.py:85`): the text before the first `-`, else the
/// `[A-Za-z]+[0-9]+` run immediately followed by `_`, else the whole id.
pub fn sample_prefix(sample_id: &str) -> &str {
    if let Some((prefix, _)) = sample_id.split_once('-') {
        return prefix;
    }
    if let Some(prefix) = leading_alpha_digit_prefix(sample_id) {
        return prefix;
    }
    sample_id
}

/// Match Python's `^([a-zA-Z]+\d+)_` capture: one or more ASCII letters, then one
/// or more ASCII digits, immediately followed by `_`. Returns the letters+digits
/// run (without the trailing `_`) when it matches.
fn leading_alpha_digit_prefix(name: &str) -> Option<&str> {
    let bytes = name.as_bytes();
    let mut end = 0;
    while end < bytes.len() && bytes[end].is_ascii_alphabetic() {
        end += 1;
    }
    if end == 0 {
        return None;
    }
    let digits_start = end;
    while end < bytes.len() && bytes[end].is_ascii_digit() {
        end += 1;
    }
    if end == digits_start {
        return None;
    }
    (bytes.get(end) == Some(&b'_')).then(|| &name[..end])
}

/// `pd.factorize` first-occurrence, 0-indexed labels for `values`.
pub fn factorize(values: &[String]) -> Vec<usize> {
    let mut labels: HashMap<&str, usize> = HashMap::new();
    let mut result = Vec::with_capacity(values.len());
    for value in values {
        let next = labels.len();
        result.push(*labels.entry(value.as_str()).or_insert(next));
    }
    result
}

/// Detect 0-indexed batch labels for `sample_ids`, mirroring `detect_batches`
/// (`batch_correction.py:116`). `run_info` maps a sample id to its run /
/// fraction / tech-replicate value (the same map the Python `run_info` argument
/// carries). `column_values` supplies the explicit per-sample batch values for
/// the `column` method. Missing required inputs return [`DetectBatchesError`],
/// matching Python's `ValueError`s.
pub fn detect_batches(
    sample_ids: &[String],
    method: BatchDetectionMethod,
    run_info: Option<&HashMap<String, String>>,
    column_values: Option<&[String]>,
) -> Result<Vec<usize>, DetectBatchesError> {
    match method {
        BatchDetectionMethod::SamplePrefix => Ok(factorize(
            &sample_ids
                .iter()
                .map(|id| sample_prefix(id).to_owned())
                .collect::<Vec<_>>(),
        )),
        BatchDetectionMethod::RunName => {
            let run_info = run_info.ok_or(DetectBatchesError::RunInfoRequired)?;
            // Python `run_info.get(s, s)`: default to the sample id itself.
            let values = sample_ids
                .iter()
                .map(|id| run_info.get(id).cloned().unwrap_or_else(|| id.clone()))
                .collect::<Vec<_>>();
            Ok(factorize(&values))
        }
        BatchDetectionMethod::Fraction | BatchDetectionMethod::TechReplicate => match run_info {
            None => {
                // Python `_fallback_prefix_batches`: warn + recurse on prefix.
                Ok(factorize(
                    &sample_ids
                        .iter()
                        .map(|id| sample_prefix(id).to_owned())
                        .collect::<Vec<_>>(),
                ))
            }
            Some(run_info) => {
                // Python `run_info.get(s, "1")`: default missing samples to "1".
                let values = sample_ids
                    .iter()
                    .map(|id| run_info.get(id).cloned().unwrap_or_else(|| "1".to_owned()))
                    .collect::<Vec<_>>();
                Ok(factorize(&values))
            }
        },
        BatchDetectionMethod::ExplicitColumn => {
            let values = column_values.ok_or(DetectBatchesError::ColumnValuesRequired)?;
            if values.len() != sample_ids.len() {
                return Err(DetectBatchesError::ColumnLengthMismatch {
                    values: values.len(),
                    samples: sample_ids.len(),
                });
            }
            Ok(factorize(values))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn samples(names: &[&str]) -> Vec<String> {
        names.iter().map(|s| (*s).to_owned()).collect()
    }

    #[test]
    fn from_str_normalizes_like_python() {
        assert_eq!(
            BatchDetectionMethod::parse_name("Sample-Prefix"),
            Some(BatchDetectionMethod::SamplePrefix)
        );
        assert_eq!(
            BatchDetectionMethod::parse_name("RUN"),
            Some(BatchDetectionMethod::RunName)
        );
        assert_eq!(
            BatchDetectionMethod::parse_name("TechReplicate"),
            Some(BatchDetectionMethod::TechReplicate)
        );
        // The Python enum value is `techreplicate` (no separator), so a spaced or
        // hyphenated form normalises to `tech_replicate` and does not match,
        // matching Python's `from_str` raising `ValueError` on those inputs.
        assert_eq!(BatchDetectionMethod::parse_name("Tech Replicate"), None);
        assert_eq!(BatchDetectionMethod::parse_name("tech-replicate"), None);
        assert_eq!(BatchDetectionMethod::parse_name("unknown"), None);
    }

    #[test]
    fn sample_prefix_matches_python() {
        // Python `_sample_prefix` oracle.
        assert_eq!(sample_prefix("PXD001-S1"), "PXD001");
        assert_eq!(sample_prefix("a-b-c"), "a");
        assert_eq!(sample_prefix("PXD001_S1"), "PXD001");
        assert_eq!(sample_prefix("sampleA"), "sampleA");
        assert_eq!(sample_prefix("PXD001S2_S1"), "PXD001S2_S1");
    }

    #[test]
    fn detect_sample_prefix_factorizes() {
        let ids = samples(&["PXD001-S1", "PXD001-S2", "PXD002-S1", "PXD002-S2"]);
        let labels = detect_batches(&ids, BatchDetectionMethod::SamplePrefix, None, None);
        assert_eq!(labels, Ok(vec![0, 0, 1, 1]));
    }

    #[test]
    fn detect_run_requires_run_info() {
        let ids = samples(&["S1", "S2"]);
        assert_eq!(
            detect_batches(&ids, BatchDetectionMethod::RunName, None, None),
            Err(DetectBatchesError::RunInfoRequired)
        );
        let mut run_info = HashMap::new();
        run_info.insert("S1".to_owned(), "runA".to_owned());
        run_info.insert("S2".to_owned(), "runB".to_owned());
        assert_eq!(
            detect_batches(&ids, BatchDetectionMethod::RunName, Some(&run_info), None),
            Ok(vec![0, 1])
        );
        // Missing sample defaults to its own id (a distinct factorize code).
        let ids3 = samples(&["S1", "S2", "S3"]);
        assert_eq!(
            detect_batches(&ids3, BatchDetectionMethod::RunName, Some(&run_info), None),
            Ok(vec![0, 1, 2])
        );
    }

    #[test]
    fn detect_fraction_falls_back_to_prefix_without_run_info() {
        let ids = samples(&["PXD001-S1", "PXD001-S2", "PXD002-S1"]);
        // No run_info -> prefix fallback.
        assert_eq!(
            detect_batches(&ids, BatchDetectionMethod::Fraction, None, None),
            Ok(vec![0, 0, 1])
        );
        assert_eq!(
            detect_batches(&ids, BatchDetectionMethod::TechReplicate, None, None),
            Ok(vec![0, 0, 1])
        );
    }

    #[test]
    fn detect_fraction_defaults_missing_to_one() {
        let ids = samples(&["S1", "S2", "S3"]);
        let mut run_info = HashMap::new();
        run_info.insert("S1".to_owned(), "2".to_owned());
        run_info.insert("S2".to_owned(), "2".to_owned());
        // S3 missing -> "1"; S1/S2 -> "2". Factorize first-occurrence: 2->0, 1->1.
        assert_eq!(
            detect_batches(&ids, BatchDetectionMethod::Fraction, Some(&run_info), None),
            Ok(vec![0, 0, 1])
        );
    }

    #[test]
    fn detect_column_validates_length() {
        let ids = samples(&["S1", "S2", "S3"]);
        assert_eq!(
            detect_batches(&ids, BatchDetectionMethod::ExplicitColumn, None, None),
            Err(DetectBatchesError::ColumnValuesRequired)
        );
        let values = samples(&["A", "B"]);
        assert_eq!(
            detect_batches(
                &ids,
                BatchDetectionMethod::ExplicitColumn,
                None,
                Some(&values)
            ),
            Err(DetectBatchesError::ColumnLengthMismatch {
                values: 2,
                samples: 3
            })
        );
        let values_ok = samples(&["A", "A", "B"]);
        assert_eq!(
            detect_batches(
                &ids,
                BatchDetectionMethod::ExplicitColumn,
                None,
                Some(&values_ok)
            ),
            Ok(vec![0, 0, 1])
        );
    }
}
