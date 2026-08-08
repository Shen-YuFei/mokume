use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::path::Path;

use csv::{Reader, ReaderBuilder, StringRecord, Trim};
use mokume_core::{MokumeError, Result};

use crate::{normalize_label_key, QpxFeatureRecord, SdrfTable};

const TMT6: &[&str] = &["TMT126", "TMT127", "TMT128", "TMT129", "TMT130", "TMT131"];
const TMT10: &[&str] = &[
    "TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C", "TMT129N", "TMT129C", "TMT130N",
    "TMT130C", "TMT131",
];
const TMT11: &[&str] = &[
    "TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C", "TMT129N", "TMT129C", "TMT130N",
    "TMT130C", "TMT131N", "TMT131C",
];
const TMT16: &[&str] = &[
    "TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C", "TMT129N", "TMT129C", "TMT130N",
    "TMT130C", "TMT131N", "TMT131C", "TMT132N", "TMT132C", "TMT133N", "TMT133C", "TMT134N",
];
const ITRAQ4: &[&str] = &["ITRAQ114", "ITRAQ115", "ITRAQ116", "ITRAQ117"];
const ITRAQ8: &[&str] = &[
    "ITRAQ113", "ITRAQ114", "ITRAQ115", "ITRAQ116", "ITRAQ117", "ITRAQ118", "ITRAQ119", "ITRAQ121",
];

#[derive(Debug, Clone, Copy)]
struct Columns {
    protein: usize,
    peptide: usize,
    charge: usize,
    intensity: usize,
    run: Option<usize>,
    reference: Option<usize>,
    channel: Option<usize>,
}

impl Columns {
    fn from_headers(headers: &StringRecord) -> Result<Self> {
        let columns = Self {
            protein: required_column(headers, &["ProteinName"])?,
            peptide: required_column(headers, &["PeptideSequence"])?,
            charge: required_column(headers, &["Charge", "PrecursorCharge"])?,
            intensity: required_column(headers, &["Intensity"])?,
            run: column(headers, &["Run"]),
            reference: column(headers, &["Reference"]),
            channel: column(headers, &["Channel"]),
        };
        if columns.run.is_none() && columns.reference.is_none() {
            return Err(invalid_input(
                "MSstats input requires a Run or Reference column",
            ));
        }
        Ok(columns)
    }
}

#[derive(Debug)]
enum Labeling {
    LabelFree,
    Isobaric {
        channels: &'static [&'static str],
        declared: HashMap<String, String>,
    },
}

impl Labeling {
    fn from_sdrf(sdrf: &SdrfTable) -> Result<Self> {
        let labels = sdrf
            .records()
            .iter()
            .filter_map(|record| record.label.as_deref())
            .map(normalize_label_key)
            .filter(|label| !label.is_empty())
            .collect::<HashSet<_>>();

        if labels.len() == 1 && labels.contains("label free sample") {
            return Ok(Self::LabelFree);
        }

        let channels = if labels.iter().any(|label| label.contains("tmt")) {
            if labels.len() > 11
                || ["tmt132n", "tmt132c", "tmt133n", "tmt133c", "tmt134n"]
                    .iter()
                    .any(|label| labels.contains(*label))
            {
                TMT16
            } else if labels.len() == 11 || labels.contains("tmt131c") {
                TMT11
            } else if labels.len() > 6 {
                TMT10
            } else {
                TMT6
            }
        } else if labels.iter().any(|label| label.contains("itraq")) {
            if labels.len() > 4 {
                ITRAQ8
            } else {
                ITRAQ4
            }
        } else {
            return Err(invalid_input(format!(
                "cannot infer LFQ, TMT, or iTRAQ labeling from SDRF labels: {labels:?}"
            )));
        };

        let canonical = channels
            .iter()
            .map(|label| (normalize_label_key(label), (*label).to_owned()))
            .collect::<HashMap<_, _>>();
        let declared = labels
            .into_iter()
            .filter_map(|label| canonical.get(&label).cloned().map(|value| (label, value)))
            .collect();
        Ok(Self::Isobaric { channels, declared })
    }

    fn resolve_channel(&self, value: Option<&str>) -> Result<Option<String>> {
        let Self::Isobaric { channels, declared } = self else {
            return Ok(None);
        };
        let text = value
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .ok_or_else(|| {
                invalid_input("MSstats input is missing a channel value for an isobaric row")
            })?;
        let canonical = text
            .parse::<f64>()
            .ok()
            .filter(|number| number.is_finite() && number.fract() == 0.0 && *number >= 1.0)
            .and_then(|number| channels.get(number as usize - 1))
            .map(|label| (*label).to_owned())
            .or_else(|| declared.get(&normalize_label_key(text)).cloned())
            .ok_or_else(|| {
                invalid_input(format!(
                    "MSstats channel `{text}` is not declared by the SDRF plex"
                ))
            })?;
        Ok(Some(canonical))
    }
}

pub struct MsstatsReader<'a> {
    reader: Reader<File>,
    columns: Columns,
    sdrf: &'a SdrfTable,
    labeling: Labeling,
    sequence_proteins: HashMap<String, String>,
    non_unique_sequences: HashSet<String>,
    row_number: usize,
}

impl<'a> MsstatsReader<'a> {
    pub fn open(path: impl AsRef<Path>, sdrf: &'a SdrfTable) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        if !path.exists() {
            return Err(MokumeError::MissingInput { path });
        }

        let mut first_pass = open_csv(&path)?;
        let headers = csv_headers(&mut first_pass, &path)?;
        let columns = Columns::from_headers(&headers)?;
        let (sequence_proteins, non_unique_sequences, rows) =
            collect_sequence_proteins(&mut first_pass, columns, &path)?;
        if rows == 0 {
            return Err(invalid_input("MSstats input contains no rows"));
        }

        let labeling = Labeling::from_sdrf(sdrf)?;
        if matches!(labeling, Labeling::Isobaric { .. }) && columns.channel.is_none() {
            return Err(invalid_input(
                "MSstats input is missing required column: Channel",
            ));
        }

        let mut reader = open_csv(&path)?;
        csv_headers(&mut reader, &path)?;
        Ok(Self {
            reader,
            columns,
            sdrf,
            labeling,
            sequence_proteins,
            non_unique_sequences,
            row_number: 1,
        })
    }

    fn parse_record(&self, record: &StringRecord, row: usize) -> Result<QpxFeatureRecord> {
        let protein = required_value(record, self.columns.protein, "ProteinName", row)?;
        let peptidoform = required_value(record, self.columns.peptide, "PeptideSequence", row)?;
        let sequence = canonical_sequence(peptidoform);
        if sequence.is_empty() {
            return Err(row_error(
                row,
                "PeptideSequence does not contain an amino-acid sequence",
            ));
        }
        let charge = parse_charge(required_value(record, self.columns.charge, "Charge", row)?)
            .ok_or_else(|| row_error(row, "Charge is not a valid 16-bit integer"))?;
        let intensity = required_value(record, self.columns.intensity, "Intensity", row)?
            .parse::<f32>()
            .ok()
            .filter(|value| value.is_finite())
            .map(f64::from)
            .ok_or_else(|| row_error(row, "Intensity is not a finite number"))?;

        let run_file_name = self.resolve_run_file(record, row)?;

        let label = self
            .labeling
            .resolve_channel(optional_value(record, self.columns.channel))?;
        self.sdrf.lookup(&run_file_name, label.as_deref())?;

        let protein_accessions = protein
            .split(';')
            .map(str::trim)
            .map(ToOwned::to_owned)
            .collect::<Vec<_>>();
        let unique = protein_accessions.len() == 1
            && !self.non_unique_sequences.contains(&sequence)
            && self.sequence_proteins.get(&sequence).map(String::as_str) == Some(protein);

        Ok(QpxFeatureRecord {
            sequence,
            peptidoform: peptidoform.to_owned(),
            charge,
            run_file_name,
            sample_accession: None,
            protein_accessions,
            anchor_protein: None,
            unique: Some(unique),
            is_decoy: None,
            pg_global_qvalue: None,
            label,
            intensity,
        })
    }

    fn resolve_run_file(&self, record: &StringRecord, row: usize) -> Result<String> {
        let reference_file = optional_value(record, self.columns.reference)
            .and_then(|value| self.sdrf.resolve_run_file(value));
        let run_file = optional_value(record, self.columns.run)
            .and_then(|value| self.sdrf.resolve_run_file(value));
        if reference_file.is_some() && run_file.is_some() && reference_file != run_file {
            return Err(row_error(
                row,
                "Run and Reference map to different SDRF data files",
            ));
        }
        reference_file.or(run_file).ok_or_else(|| {
            row_error(
                row,
                "Run and Reference cannot be mapped to an SDRF data file",
            )
        })
    }
}

impl Iterator for MsstatsReader<'_> {
    type Item = Result<QpxFeatureRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        let record = self.reader.records().next()?;
        self.row_number += 1;
        Some(
            record
                .map_err(|source| {
                    row_error(self.row_number, format!("failed to read CSV row: {source}"))
                })
                .and_then(|record| self.parse_record(&record, self.row_number)),
        )
    }
}

fn open_csv(path: &Path) -> Result<Reader<File>> {
    ReaderBuilder::new()
        .flexible(true)
        .trim(Trim::All)
        .from_path(path)
        .map_err(|source| csv_error(path, source))
}

fn csv_headers(reader: &mut Reader<File>, path: &Path) -> Result<StringRecord> {
    reader
        .headers()
        .cloned()
        .map_err(|source| csv_error(path, source))
}

fn collect_sequence_proteins(
    reader: &mut Reader<File>,
    columns: Columns,
    path: &Path,
) -> Result<(HashMap<String, String>, HashSet<String>, usize)> {
    let mut proteins = HashMap::new();
    let mut non_unique = HashSet::new();
    let mut rows = 0;
    for record in reader.records() {
        rows += 1;
        let record = record.map_err(|source| csv_error(path, source))?;
        let sequence = canonical_sequence(record.get(columns.peptide).unwrap_or(""));
        let protein = record.get(columns.protein).unwrap_or("").trim();
        if sequence.is_empty() || protein.is_empty() {
            continue;
        }
        match proteins.get(&sequence) {
            Some(known) if known != protein => {
                non_unique.insert(sequence);
            }
            None => {
                proteins.insert(sequence, protein.to_owned());
            }
            Some(_) => {}
        }
    }
    Ok((proteins, non_unique, rows))
}

fn canonical_sequence(value: &str) -> String {
    let mut sequence = String::new();
    let mut parentheses = 0_u32;
    let mut brackets = 0_u32;
    for character in value.trim().chars() {
        match character {
            '(' => parentheses += 1,
            ')' => parentheses = parentheses.saturating_sub(1),
            '[' => brackets += 1,
            ']' => brackets = brackets.saturating_sub(1),
            value if parentheses == 0 && brackets == 0 && value.is_ascii_alphabetic() => {
                sequence.push(value.to_ascii_uppercase());
            }
            _ => {}
        }
    }
    sequence
}

fn parse_charge(value: &str) -> Option<i32> {
    let value = value.trim().parse::<f64>().ok()?;
    (value.is_finite()
        && value.fract() == 0.0
        && value >= f64::from(i16::MIN)
        && value <= f64::from(i16::MAX))
    .then_some(value as i32)
}

fn column(headers: &StringRecord, aliases: &[&str]) -> Option<usize> {
    aliases.iter().find_map(|alias| {
        headers
            .iter()
            .position(|header| header.trim().eq_ignore_ascii_case(alias))
    })
}

fn required_column(headers: &StringRecord, aliases: &[&str]) -> Result<usize> {
    column(headers, aliases).ok_or_else(|| {
        invalid_input(format!(
            "MSstats input is missing required column: {}",
            aliases[0]
        ))
    })
}

fn optional_value(record: &StringRecord, column: Option<usize>) -> Option<&str> {
    column
        .and_then(|index| record.get(index))
        .map(str::trim)
        .filter(|value| !value.is_empty())
}

fn required_value<'a>(
    record: &'a StringRecord,
    column: usize,
    name: &str,
    row: usize,
) -> Result<&'a str> {
    record
        .get(column)
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .ok_or_else(|| row_error(row, format!("{name} is empty")))
}

fn csv_error(path: &Path, source: csv::Error) -> MokumeError {
    invalid_input(format!(
        "failed to read MSstats CSV `{}`: {source}",
        path.display()
    ))
}

fn row_error(row: usize, message: impl Into<String>) -> MokumeError {
    invalid_input(format!("invalid MSstats row {row}: {}", message.into()))
}

fn invalid_input(message: impl Into<String>) -> MokumeError {
    MokumeError::InvalidInput {
        message: message.into(),
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use tempfile::NamedTempFile;

    use super::*;

    fn msstats_file(contents: &str) -> Result<NamedTempFile> {
        let mut file = NamedTempFile::new().map_err(|source| MokumeError::Io {
            path: "temporary MSstats.csv".into(),
            source,
        })?;
        file.write_all(contents.as_bytes())
            .map_err(|source| MokumeError::Io {
                path: file.path().to_path_buf(),
                source,
            })?;
        Ok(file)
    }

    #[test]
    fn reads_lfq_reference_alias_and_canonicalizes_sequence() -> Result<()> {
        let sdrf = SdrfTable::from_reader(
            concat!(
                "source name\tassay name\tcomment[data file]\tcomment[label]\n",
                "sample-1\trun-1\tsample.raw\tlabel free sample\n",
            )
            .as_bytes(),
        )?;
        let input = msstats_file(concat!(
            "ProteinName,PeptideSequence,Charge,Reference,Intensity\n",
            "P12345,K.AC(ox)DE[+1]F.R,2,C:\\\\data\\\\sample.raw,123.125\n",
        ))?;

        let rows = MsstatsReader::open(input.path(), &sdrf)?.collect::<Result<Vec<_>>>()?;

        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].sequence, "KACDEFR");
        assert_eq!(rows[0].run_file_name, "sample.raw");
        assert_eq!(rows[0].label, None);
        assert_eq!(rows[0].unique, Some(true));
        Ok(())
    }

    #[test]
    fn maps_tmt10_channel_position_to_canonical_label() -> Result<()> {
        let labels = TMT10
            .iter()
            .enumerate()
            .map(|(index, label)| format!("sample-{}\tplex.raw\t{label}\n", index + 1))
            .collect::<String>();
        let sdrf = SdrfTable::from_reader(
            format!("source name\tcomment[data file]\tcomment[label]\n{labels}").as_bytes(),
        )?;
        let input = msstats_file(concat!(
            "ProteinName,PeptideSequence,PrecursorCharge,Run,Channel,Intensity\n",
            "P12345,PEPTIDE,3,plex,10,42\n",
        ))?;

        let row = MsstatsReader::open(input.path(), &sdrf)?
            .next()
            .ok_or_else(|| invalid_input("missing MSstats row"))??;

        assert_eq!(row.run_file_name, "plex.raw");
        assert_eq!(row.label.as_deref(), Some("TMT131"));
        Ok(())
    }

    #[test]
    fn rejects_conflicting_run_and_reference_aliases() -> Result<()> {
        let sdrf = SdrfTable::from_reader(
            concat!(
                "source name\tcomment[data file]\tcomment[label]\n",
                "sample-1\trun-1.raw\tlabel free sample\n",
                "sample-2\trun-2.raw\tlabel free sample\n",
            )
            .as_bytes(),
        )?;
        let input = msstats_file(concat!(
            "ProteinName,PeptideSequence,Charge,Run,Reference,Intensity\n",
            "P12345,PEPTIDE,2,run-1,run-2,42\n",
        ))?;

        let message = MsstatsReader::open(input.path(), &sdrf)?
            .next()
            .ok_or_else(|| invalid_input("missing MSstats row"))?
            .err()
            .ok_or_else(|| invalid_input("conflicting aliases unexpectedly succeeded"))?
            .to_string();

        assert!(message.contains("Run and Reference map to different"));
        Ok(())
    }
}
