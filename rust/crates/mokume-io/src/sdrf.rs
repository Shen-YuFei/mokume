use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::path::Path;

use csv::{ReaderBuilder, StringRecord, Trim};
use mokume_core::{MokumeError, Result};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SdrfRecord {
    pub sample_accession: String,
    pub run_accession: Option<String>,
    pub data_file: String,
    pub file_uri: Option<String>,
    pub label: Option<String>,
    pub fraction: Option<String>,
    pub biological_replicate: Option<u32>,
    pub technical_replicate: Option<u32>,
    pub condition: Option<String>,
    /// Raw value of `characteristics[pooled sample]` (SDRF-Proteomics spec). The
    /// ratio reference detection treats `pooled` (or an `SN=...` member list) as a
    /// reference channel, mirroring Python's `detect_pooled_from_sdrf`.
    pub pooled_sample: Option<String>,
}

#[derive(Debug, Clone, Default)]
pub struct SdrfTable {
    records: Vec<SdrfRecord>,
    by_run: HashMap<String, usize>,
    by_run_label: HashMap<(String, String), usize>,
}

impl SdrfTable {
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        if !path.exists() {
            return Err(MokumeError::MissingInput { path });
        }

        let file = File::open(&path).map_err(|source| MokumeError::Io {
            path: path.clone(),
            source,
        })?;
        Self::from_reader(file)
    }

    pub fn from_reader(reader: impl Read) -> Result<Self> {
        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .flexible(true)
            .trim(Trim::All)
            .from_reader(reader);
        let headers = reader
            .headers()
            .map_err(|source| invalid_input(format!("failed to read SDRF header: {source}")))?
            .clone();
        let columns = SdrfColumns::from_headers(&headers)?;

        let mut records = Vec::new();
        for record in reader.records() {
            let record = record
                .map_err(|source| invalid_input(format!("failed to read SDRF row: {source}")))?;
            records.push(columns.to_record(&record)?);
        }

        Ok(Self::from_records(records))
    }

    pub fn from_records(records: Vec<SdrfRecord>) -> Self {
        let mut table = Self {
            records,
            by_run: HashMap::new(),
            by_run_label: HashMap::new(),
        };
        table.rebuild_index();
        table
    }

    pub fn len(&self) -> usize {
        self.records.len()
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    pub fn records(&self) -> &[SdrfRecord] {
        &self.records
    }

    pub fn lookup(&self, run_file_name: &str, label: Option<&str>) -> Option<&SdrfRecord> {
        let run_key = normalize_file_key(run_file_name);
        if let Some(label) = label {
            let label_key = normalize_label_key(label);
            if let Some(index) = self.by_run_label.get(&(run_key.clone(), label_key)) {
                return self.records.get(*index);
            }
        }

        self.by_run
            .get(&run_key)
            .and_then(|index| self.records.get(*index))
    }

    fn rebuild_index(&mut self) {
        for (index, record) in self.records.iter().enumerate() {
            for key in record.run_keys() {
                self.by_run.entry(key.clone()).or_insert(index);
                if let Some(label) = record.label.as_deref() {
                    self.by_run_label
                        .entry((key, normalize_label_key(label)))
                        .or_insert(index);
                }
            }
        }
    }
}

/// Raw SDRF rows for covariate / explicit-batch-column extraction, mirroring
/// Python's `load_sdrf` (`pd.read_csv(sep="\t")` then lowercasing every header).
/// The typed [`SdrfTable`] only keeps the spec columns it needs for run/label
/// lookup, but batch correction reads arbitrary `characteristics[*]` /
/// `factor value[*]` / custom batch columns by name, so it parses the file
/// through this verbatim view instead. Headers are lowercased and trimmed to
/// match Python's column addressing; cell values are kept case-sensitive
/// (factorize encodes by value identity).
#[derive(Debug, Clone, Default)]
pub struct SdrfRawTable {
    headers: Vec<String>,
    rows: Vec<Vec<String>>,
}

impl SdrfRawTable {
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        if !path.exists() {
            return Err(MokumeError::MissingInput { path });
        }
        let file = File::open(&path).map_err(|source| MokumeError::Io {
            path: path.clone(),
            source,
        })?;
        Self::from_reader(file)
    }

    pub fn from_reader(reader: impl Read) -> Result<Self> {
        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .flexible(true)
            .trim(Trim::All)
            .from_reader(reader);
        let headers = reader
            .headers()
            .map_err(|source| invalid_input(format!("failed to read SDRF header: {source}")))?
            .iter()
            .map(normalize_header)
            .collect();
        let mut rows = Vec::new();
        for record in reader.records() {
            let record = record
                .map_err(|source| invalid_input(format!("failed to read SDRF row: {source}")))?;
            rows.push(record.iter().map(ToOwned::to_owned).collect());
        }
        Ok(Self { headers, rows })
    }

    pub fn headers(&self) -> &[String] {
        &self.headers
    }

    /// Index of the column whose lowercased header equals `name` (which the
    /// caller lowercases), or `None` when absent.
    pub fn column_index(&self, name: &str) -> Option<usize> {
        self.headers.iter().position(|header| header == name)
    }

    pub fn row_count(&self) -> usize {
        self.rows.len()
    }

    /// Cell at `row`/`col`, or `""` when the (flexible) row is short, matching
    /// how a missing trailing SDRF cell reads as empty.
    pub fn cell(&self, row: usize, col: usize) -> &str {
        self.rows
            .get(row)
            .and_then(|cells| cells.get(col))
            .map_or("", String::as_str)
    }
}

impl SdrfRecord {
    fn run_keys(&self) -> Vec<String> {
        let mut keys = Vec::with_capacity(3);
        push_key(&mut keys, &self.data_file);
        if let Some(file_uri) = self.file_uri.as_deref() {
            push_key(&mut keys, file_uri);
        }
        if let Some(run_accession) = self.run_accession.as_deref() {
            push_key(&mut keys, run_accession);
        }
        keys
    }
}

#[derive(Debug, Clone)]
struct SdrfColumns {
    source_name: Option<usize>,
    assay_name: Option<usize>,
    data_file: Option<usize>,
    file_uri: Option<usize>,
    label: Option<usize>,
    fraction: Option<usize>,
    biological_replicate: Option<usize>,
    technical_replicate: Option<usize>,
    condition: Option<usize>,
    pooled_sample: Option<usize>,
}

impl SdrfColumns {
    fn from_headers(headers: &StringRecord) -> Result<Self> {
        let data_file = header_index(
            headers,
            &["comment[data file]", "comment[raw file]", "data file"],
        );
        let file_uri = header_index(
            headers,
            &["comment[file uri]", "comment[associated file uri]"],
        );
        if data_file.is_none() && file_uri.is_none() {
            return Err(invalid_input(
                "SDRF must contain `comment[data file]` or `comment[file uri]`",
            ));
        }

        Ok(Self {
            source_name: header_index(headers, &["source name", "sample accession"]),
            assay_name: header_index(headers, &["assay name", "run accession"]),
            data_file,
            file_uri,
            label: header_index(headers, &["comment[label]", "label"]),
            fraction: header_index(headers, &["comment[fraction identifier]", "fraction"]),
            biological_replicate: header_index(
                headers,
                &[
                    "characteristics[biological replicate]",
                    "comment[biological replicate]",
                    "biological replicate",
                ],
            ),
            technical_replicate: header_index(
                headers,
                &[
                    "comment[technical replicate]",
                    "technical replicate",
                    "technical_replica",
                ],
            ),
            condition: first_factor_value(headers).or_else(|| {
                header_index(
                    headers,
                    &[
                        "factor value[cell line]",
                        "characteristics[organism part]",
                        "characteristics[cell line]",
                        "characteristics[disease]",
                    ],
                )
            }),
            pooled_sample: header_index(headers, &["characteristics[pooled sample]"]),
        })
    }

    fn to_record(&self, record: &StringRecord) -> Result<SdrfRecord> {
        let file_uri = optional_cell(record, self.file_uri).map(ToOwned::to_owned);
        let data_file = optional_cell(record, self.data_file)
            .map(ToOwned::to_owned)
            .or_else(|| file_uri.as_deref().map(basename).map(ToOwned::to_owned))
            .ok_or_else(|| invalid_input("SDRF row has no data file"))?;
        let sample_accession = optional_cell(record, self.source_name)
            .map(ToOwned::to_owned)
            .unwrap_or_else(|| data_file.clone());

        Ok(SdrfRecord {
            sample_accession,
            run_accession: optional_cell(record, self.assay_name).map(ToOwned::to_owned),
            data_file,
            file_uri,
            label: optional_cell(record, self.label).map(extract_label_value),
            fraction: optional_cell(record, self.fraction).map(ToOwned::to_owned),
            biological_replicate: optional_cell(record, self.biological_replicate)
                .and_then(parse_positive_u32),
            technical_replicate: optional_cell(record, self.technical_replicate)
                .and_then(parse_positive_u32),
            condition: optional_cell(record, self.condition).map(ToOwned::to_owned),
            pooled_sample: optional_cell(record, self.pooled_sample).map(ToOwned::to_owned),
        })
    }
}

fn push_key(keys: &mut Vec<String>, value: &str) {
    let key = normalize_file_key(value);
    if !key.is_empty() && !keys.iter().any(|known| known == &key) {
        keys.push(key);
    }
}

pub fn normalize_file_key(value: &str) -> String {
    let key = basename(value)
        .trim()
        .trim_start_matches("file://")
        .to_ascii_lowercase();
    // Strip a trailing run extension so the SDRF `comment[data file]` key matches
    // the QPX `run_file_name`, which is stored extension-less. Without this the
    // per-feature SDRF lookup misses and the Condition column falls back to the
    // run filename. Mirrors `mokume-pipeline` de.rs::strip_run_extension.
    for ext in [".raw", ".mzml", ".d", ".wiff"] {
        if let Some(stem) = key.strip_suffix(ext) {
            return stem.to_owned();
        }
    }
    key
}

pub fn normalize_label_key(value: &str) -> String {
    extract_label_value(value).trim().to_ascii_lowercase()
}

fn extract_label_value(value: &str) -> String {
    let trimmed = value.trim();
    for part in trimmed.split(';') {
        let part = part.trim();
        if let Some(value) = part.strip_prefix("NT=") {
            return value.trim().to_owned();
        }
    }
    trimmed.to_owned()
}

fn basename(value: &str) -> &str {
    let trimmed = value.trim();
    let without_forward = trimmed
        .rsplit_once('/')
        .map_or(trimmed, |(_, basename)| basename);
    without_forward
        .rsplit_once('\\')
        .map_or(without_forward, |(_, basename)| basename)
}

fn header_index(headers: &StringRecord, candidates: &[&str]) -> Option<usize> {
    headers.iter().enumerate().find_map(|(index, header)| {
        let header = normalize_header(header);
        candidates
            .iter()
            .any(|candidate| header == normalize_header(candidate))
            .then_some(index)
    })
}

fn first_factor_value(headers: &StringRecord) -> Option<usize> {
    headers.iter().enumerate().find_map(|(index, header)| {
        normalize_header(header)
            .starts_with("factor value[")
            .then_some(index)
    })
}

fn normalize_header(value: &str) -> String {
    value.trim().to_ascii_lowercase()
}

fn optional_cell(record: &StringRecord, index: Option<usize>) -> Option<&str> {
    index
        .and_then(|index| record.get(index))
        .map(str::trim)
        .filter(|value| !value.is_empty())
}

fn parse_positive_u32(value: &str) -> Option<u32> {
    value.trim().parse::<u32>().ok()
}

fn invalid_input(message: impl Into<String>) -> MokumeError {
    MokumeError::InvalidInput {
        message: message.into(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_sdrf_and_uses_label_aware_lookup() -> Result<()> {
        let input = concat!(
            "source name\tassay name\tcomment[data file]\tcomment[label]\tcomment[fraction identifier]\tcomment[technical replicate]\tfactor value[cell line]\n",
            "sample-1\trun 1\trun1.raw\tAC=MS:1002038;NT=label free sample\t1\t2\tHeLa\n",
            "sample-2\trun 2\trun2.raw\tTMT126\t1\t1\tK562\n",
        );

        let table = SdrfTable::from_reader(input.as_bytes())?;
        let label_free = table
            .lookup("/tmp/run1.raw", Some("label free sample"))
            .ok_or_else(|| invalid_input("missing label-free lookup"))?;
        let tmt = table
            .lookup("run2.raw", Some("TMT126"))
            .ok_or_else(|| invalid_input("missing TMT lookup"))?;

        assert_eq!(table.len(), 2);
        assert_eq!(label_free.sample_accession, "sample-1");
        assert_eq!(label_free.condition.as_deref(), Some("HeLa"));
        assert_eq!(label_free.technical_replicate, Some(2));
        assert_eq!(tmt.sample_accession, "sample-2");
        Ok(())
    }

    #[test]
    fn raw_table_lowercases_headers_and_addresses_cells_by_name() -> Result<()> {
        let input = concat!(
            "Source Name\tcharacteristics[sex]\tbatch\n",
            "S1\tmale\tb1\n",
            "S2\tfemale\tb2\n",
        );
        let raw = SdrfRawTable::from_reader(input.as_bytes())?;

        // Headers are lowercased (Python's `load_sdrf`); cells stay verbatim.
        assert_eq!(
            raw.headers(),
            ["source name", "characteristics[sex]", "batch"]
        );
        assert_eq!(raw.column_index("source name"), Some(0));
        assert_eq!(raw.column_index("Source Name"), None);
        assert_eq!(raw.column_index("missing"), None);
        assert_eq!(raw.row_count(), 2);

        let sex = raw
            .column_index("characteristics[sex]")
            .ok_or_else(|| invalid_input("missing sex column"))?;
        assert_eq!(raw.cell(0, sex), "male");
        assert_eq!(raw.cell(1, sex), "female");
        // Out-of-range row/column reads as empty (short flexible rows).
        assert_eq!(raw.cell(9, sex), "");
        Ok(())
    }

    #[test]
    fn parses_pooled_sample_column() -> Result<()> {
        let input = concat!(
            "source name\tcomment[data file]\tcharacteristics[pooled sample]\n",
            "p1_1\tp1_ref.raw\tpooled\n",
            "p1_2\tp1_a.raw\tnot pooled\n",
        );

        let table = SdrfTable::from_reader(input.as_bytes())?;
        let pooled = table
            .lookup("p1_ref.raw", None)
            .ok_or_else(|| invalid_input("missing pooled lookup"))?;
        let other = table
            .lookup("p1_a.raw", None)
            .ok_or_else(|| invalid_input("missing non-pooled lookup"))?;

        assert_eq!(pooled.pooled_sample.as_deref(), Some("pooled"));
        assert_eq!(other.pooled_sample.as_deref(), Some("not pooled"));
        Ok(())
    }

    #[test]
    fn normalize_file_key_strips_run_extensions() {
        // The QPX `run_file_name` is stored extension-less while the SDRF
        // `comment[data file]` carries a run extension; both must normalize to
        // the same key so the per-feature lookup matches.
        for name in [
            "Run_Condition_A_01.raw",
            "Run_Condition_A_01.mzML",
            "Run_Condition_A_01.d",
            "Run_Condition_A_01.wiff",
            "/data/Run_Condition_A_01.raw",
            "Run_Condition_A_01",
        ] {
            assert_eq!(
                normalize_file_key(name),
                "run_condition_a_01",
                "key for {name}"
            );
        }
    }

    #[test]
    fn lookup_matches_extensionless_run_key() -> Result<()> {
        // Regression: SDRF data files carry `.raw`, but the QPX run key is
        // extension-less. Without extension stripping the lookup misses and the
        // Condition column falls back to the run filename.
        let input = concat!(
            "source name\tcomment[data file]\tfactor value[group]\n",
            "A1\tRun_Condition_A_01.raw\tA\n",
            "B1\tRun_Condition_B_01.raw\tB\n",
        );
        let table = SdrfTable::from_reader(input.as_bytes())?;
        let a = table
            .lookup("Run_Condition_A_01", None)
            .ok_or_else(|| invalid_input("missing extension-less lookup"))?;
        assert_eq!(a.condition.as_deref(), Some("A"));
        Ok(())
    }
}
