use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

use mokume_core::{FeatureToProteinsConfig, Result};
use mokume_io::{QpxPsmParquetReader, QpxPsmRecord, SdrfTable};

use crate::memory::MemoryPlan;

#[derive(Debug)]
pub(crate) struct SpectralCountCell {
    pub protein_group: String,
    pub sample: String,
    pub spectra: usize,
    pub sequences: HashSet<String>,
}

#[derive(Debug)]
pub(crate) struct SpectralCountResult {
    pub cells: Vec<SpectralCountCell>,
    pub target_psms: usize,
    pub unique_spectra: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct SpectrumKey {
    run: String,
    scan: Vec<i64>,
}

#[derive(Debug)]
struct SpectrumAssignment {
    proteins: BTreeSet<String>,
    sample: String,
    sequences: HashSet<String>,
}

pub(crate) fn count(
    config: &FeatureToProteinsConfig,
    sdrf: &SdrfTable,
    memory: &MemoryPlan,
) -> Result<SpectralCountResult> {
    let psm = config
        .input
        .psm
        .as_deref()
        .ok_or_else(|| super::invalid_input("spectral_count requires --psm input"))?;
    let reader = QpxPsmParquetReader::open(psm, memory.qpx_batch_size())?;
    let (assignments, target_psms) = collect_assignments(config, sdrf, memory, reader)?;
    if assignments.is_empty() {
        return Err(super::invalid_input(
            "spectral_count found no usable target PSM after decoy, peptide-length, contaminant, and protein-accession filtering",
        ));
    }
    let unique_spectra = assignments.len();
    let cells = collapse_assignments(assignments, config.filtering.min_unique_peptides);
    Ok(SpectralCountResult {
        cells,
        target_psms,
        unique_spectra,
    })
}

fn collect_assignments(
    config: &FeatureToProteinsConfig,
    sdrf: &SdrfTable,
    memory: &MemoryPlan,
    mut reader: QpxPsmParquetReader,
) -> Result<(HashMap<SpectrumKey, SpectrumAssignment>, usize)> {
    let mut assignments = HashMap::new();
    let mut run_samples = HashMap::new();
    let mut target_psms = 0;
    for batch in &mut reader {
        for record in batch? {
            target_psms += usize::from(add_record(
                config,
                sdrf,
                &mut run_samples,
                &mut assignments,
                record,
            )?);
        }
        memory.check("PSM spectral-count streaming")?;
    }
    Ok((assignments, target_psms))
}

fn add_record(
    config: &FeatureToProteinsConfig,
    sdrf: &SdrfTable,
    run_samples: &mut HashMap<String, String>,
    assignments: &mut HashMap<SpectrumKey, SpectrumAssignment>,
    record: QpxPsmRecord,
) -> Result<bool> {
    if record.is_decoy || record.sequence.len() < config.filtering.min_aa {
        return Ok(false);
    }
    let Some(proteins) = protein_set(&record.protein_accessions) else {
        return Ok(false);
    };
    if config.filtering.remove_contaminants
        && record
            .protein_accessions
            .iter()
            .any(|accession| super::is_contaminant(accession))
    {
        return Ok(false);
    }
    if record.scan.is_empty() {
        return Err(super::invalid_input(format!(
            "PSM in run `{}` has no scan identity; spectral_count cannot deduplicate it",
            record.run_file_name
        )));
    }
    let sample = sample_for_run(sdrf, run_samples, &record.run_file_name)?;
    insert_assignment(assignments, record, proteins, sample)?;
    Ok(true)
}

fn sample_for_run(
    sdrf: &SdrfTable,
    cache: &mut HashMap<String, String>,
    run: &str,
) -> Result<String> {
    if let Some(sample) = cache.get(run) {
        return Ok(sample.clone());
    }
    let sample = sdrf.lookup(run, None)?.sample_accession.clone();
    cache.insert(run.to_owned(), sample.clone());
    Ok(sample)
}

fn insert_assignment(
    assignments: &mut HashMap<SpectrumKey, SpectrumAssignment>,
    record: QpxPsmRecord,
    proteins: BTreeSet<String>,
    sample: String,
) -> Result<()> {
    let key = SpectrumKey {
        run: record.run_file_name,
        scan: record.scan,
    };
    let assignment = assignments
        .entry(key)
        .or_insert_with(|| SpectrumAssignment {
            proteins: BTreeSet::new(),
            sample: sample.clone(),
            sequences: HashSet::new(),
        });
    if assignment.sample != sample {
        return Err(super::invalid_input(
            "one spectrum maps to multiple samples in the SDRF",
        ));
    }
    assignment.proteins.extend(proteins);
    assignment.sequences.insert(record.sequence);
    Ok(())
}

fn collapse_assignments(
    assignments: HashMap<SpectrumKey, SpectrumAssignment>,
    min_unique: usize,
) -> Vec<SpectralCountCell> {
    let mut cells = BTreeMap::<(String, String), (usize, HashSet<String>)>::new();
    for assignment in assignments.into_values() {
        let protein_group = assignment
            .proteins
            .into_iter()
            .collect::<Vec<_>>()
            .join(";");
        let cell = cells.entry((protein_group, assignment.sample)).or_default();
        cell.0 += 1;
        cell.1.extend(assignment.sequences);
    }
    cells
        .into_iter()
        .filter(|(_, (_, sequences))| sequences.len() >= min_unique)
        .map(
            |((protein_group, sample), (spectra, sequences))| SpectralCountCell {
                protein_group,
                sample,
                spectra,
                sequences,
            },
        )
        .collect()
}

fn protein_set(accessions: &[String]) -> Option<BTreeSet<String>> {
    let proteins = accessions
        .iter()
        .map(|accession| super::parse_protein_accession(accession))
        .filter(|accession| !accession.is_empty())
        .collect::<BTreeSet<_>>();
    (!proteins.is_empty()).then_some(proteins)
}

#[cfg(test)]
mod tests {
    use std::collections::{BTreeSet, HashMap, HashSet};

    use mokume_io::QpxPsmRecord;

    use super::{
        collapse_assignments, insert_assignment, protein_set, SpectrumAssignment, SpectrumKey,
    };

    fn assignment(protein: &str, sample: &str, sequence: &str) -> SpectrumAssignment {
        SpectrumAssignment {
            proteins: BTreeSet::from([protein.to_owned()]),
            sample: sample.to_owned(),
            sequences: HashSet::from([sequence.to_owned()]),
        }
    }

    #[test]
    fn counts_unique_spectra_not_unique_peptides() {
        let assignments = HashMap::from([
            (
                SpectrumKey {
                    run: "run-a".to_owned(),
                    scan: vec![1],
                },
                assignment("P1", "S1", "PEPTIDEA"),
            ),
            (
                SpectrumKey {
                    run: "run-a".to_owned(),
                    scan: vec![2],
                },
                assignment("P1", "S1", "PEPTIDEA"),
            ),
            (
                SpectrumKey {
                    run: "run-a".to_owned(),
                    scan: vec![3],
                },
                assignment("P1", "S1", "PEPTIDEB"),
            ),
        ]);
        let cells = collapse_assignments(assignments, 2);
        assert_eq!(cells.len(), 1);
        assert_eq!(cells[0].spectra, 3);
        assert_eq!(cells[0].sequences.len(), 2);
    }

    #[test]
    fn merges_all_assignments_for_one_spectrum() -> Result<(), Box<dyn std::error::Error>> {
        let mut assignments = HashMap::new();
        let first = QpxPsmRecord {
            sequence: "PEPTIDEA".to_owned(),
            run_file_name: "run-a".to_owned(),
            scan: vec![1],
            protein_accessions: vec!["P1".to_owned()],
            is_decoy: false,
        };
        insert_assignment(
            &mut assignments,
            first.clone(),
            BTreeSet::from(["P1".to_owned()]),
            "S1".to_owned(),
        )?;
        insert_assignment(
            &mut assignments,
            first,
            BTreeSet::from(["P2".to_owned()]),
            "S1".to_owned(),
        )?;
        let cells = collapse_assignments(assignments, 1);
        assert_eq!(cells.len(), 1);
        assert_eq!(cells[0].protein_group, "P1;P2");
        assert_eq!(cells[0].spectra, 1);
        Ok(())
    }

    #[test]
    fn shared_protein_set_is_sorted_and_deduplicated() {
        let accessions = ["sp|P2|two", "P1", "P2"].map(str::to_owned);
        assert_eq!(
            protein_set(&accessions),
            Some(BTreeSet::from(["P1".to_owned(), "P2".to_owned()]))
        );
    }
}
