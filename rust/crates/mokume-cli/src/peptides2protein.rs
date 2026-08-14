//! `peptides2protein` command: compute per-protein quantities from a
//! peptide-level intensity matrix.
//!
//! This mirrors `mokume.commands.peptides2protein.peptides2protein` in the
//! Python package. Two code paths exist:
//!
//!   * iBAQ (`--method ibaq`, the default) reuses the existing piBAQ core via
//!     `mokume_pipeline::run_ibaq_from_peptides`. A FASTA is required to derive
//!     theoretical peptide counts. The output is the Python long-format table
//!     `ProteinName, SampleID, Condition, NormIntensity, Ibaq, FamilyId,
//!     EvidenceLevel, FamilySize`.
//!
//!   * The deterministic generic methods `sum`, `top3`, and `topn` are computed
//!     directly from the long-format rows, faithfully reproducing the Python
//!     `AllPeptidesQuantification` / `TopNQuantification.quantify` group-bys
//!     (raw per-row aggregation; no canonical/charge collapse). The output is
//!     `ProteinName, SampleID, Intensity` plus `Condition` (when the input has
//!     it) and `IntensityNorm` (when `--normalize` is set).
//!
//!   * `maxlfq` and `directlfq` roll the peptide matrix up with the DirectLFQ
//!     estimator (canonical peptides as ions) via
//!     `mokume_pipeline::run_lfq_from_peptides`, mirroring Python's
//!     `DirectLFQQuantification`; `maxlfq` delegates to DirectLFQ with
//!     `min_nonan = 2` (its `min_peptides`). The output is the same `Intensity`
//!     long-format table, keeping only `> 0` rows (Python's `_parse_wide_output`).
//!
//! iBAQ extras (P3) are computed as deterministic post-processing on the piBAQ
//! result, wired into `--tpa`, `--ruler`, and `--normalize` (see the
//! [`ibaq_extras`] module):
//!   * `--tpa`: Total Protein Approach -- `MolecularWeight` + `TPA` columns.
//!   * `--normalize`: the Python `normalize_ibaq` PRIDE/ProteomicsDB transforms
//!     (`IbaqNorm`, `IbaqLog`, `IbaqPpb`).
//!   * `--ruler`: the ProteomicRuler copy number / moles / weight /
//!     concentration columns (requires `--tpa`).
//!
//! iBAQ digestion computes natively in Rust for the pyOpenMS enzymes whose
//! cleavage rules are ported (Trypsin[/P], Lys-C[/P], Arg-C[/P], Chymotrypsin[/P],
//! Glu-C, Asp-N, Lys-N, PepsinA, and the other rules `supports_ibaq_enzyme`
//! accepts -- that function is the authoritative ported set), all with zero missed
//! cleavages; the per-enzyme digests are oracle-locked against pyOpenMS in the
//! pipeline crate. For any other enzyme pyOpenMS knows (CNBr, V8-DE, unspecific
//! cleavage, ...) the Rust kernel has no cleavage rule, so iBAQ for that enzyme is
//! not supported here and the command fails with a clear error pointing to the
//! Python wheel; the default `Trypsin` path and every other ported enzyme are
//! never affected.
//!
//! Input format: a comma-separated CSV (matching Python's `pd.read_csv`), a
//! tab-separated `.tsv`, or a parquet file (matching Python's `is_parquet`
//! magic-byte check + `pd.read_parquet`), with at minimum `ProteinName`,
//! `SampleID`, and `NormIntensity`. The peptide column is `PeptideCanonical`
//! (preferred) or `PeptideSequence`; it is required for iBAQ. `Condition` is
//! optional. Parquet `SampleID` / `Condition` columns may be dictionary-encoded
//! (pandas `Categorical`), as written by `features2peptides --save_parquet`.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, ErrorKind, Read, Write};
use std::path::Path;

use mokume_core::{MokumeError, Result};
use mokume_io::read_peptide_parquet;
use mokume_pipeline::{
    run_ibaq_from_peptides, run_lfq_from_peptides, supports_ibaq_enzyme, IbaqFromPeptidesParams,
    LfqPeptideObservation, PeptideObservation,
};

use crate::Peptides2ProteinArgs;

/// Column header for the protein identifier (Python `PROTEIN_NAME`).
const PROTEIN_NAME: &str = "ProteinName";
/// Column header for the sample identifier (Python `SAMPLE_ID`).
const SAMPLE_ID: &str = "SampleID";
/// Column header for the normalized intensity (Python `NORM_INTENSITY`).
const NORM_INTENSITY: &str = "NormIntensity";
/// Optional condition column (Python `CONDITION`).
const CONDITION: &str = "Condition";
/// Preferred peptide-sequence column (Python `PEPTIDE_CANONICAL`).
const PEPTIDE_CANONICAL: &str = "PeptideCanonical";
/// Fallback peptide-sequence column (Python `PEPTIDE_SEQUENCE`).
const PEPTIDE_SEQUENCE: &str = "PeptideSequence";

/// A parsed peptide-level input table, retaining the per-row values needed by
/// every supported method.
struct PeptideTable {
    /// `(protein, sample, peptide, condition, intensity)` per accepted row.
    rows: Vec<PeptideRow>,
    /// Whether the input carried a `Condition` column.
    has_condition: bool,
    /// Whether the input carried a peptide-sequence column.
    has_peptide: bool,
}

#[derive(Debug, Clone, PartialEq)]
struct PeptideRow {
    protein: String,
    sample: String,
    peptide: String,
    condition: String,
    intensity: f64,
}

/// Entry point for the `peptides2protein` command.
pub fn run_peptides_to_protein(args: &Peptides2ProteinArgs) -> Result<()> {
    let method = args.method.to_ascii_lowercase();

    if !args.peptides.exists() {
        return Err(MokumeError::MissingInput {
            path: args.peptides.clone(),
        });
    }
    let Some(output) = args.output.as_ref() else {
        return Err(MokumeError::InvalidInput {
            message: "peptides2protein requires --output".to_owned(),
        });
    };

    // `--verbose` in Python writes the `--qc_report` QC PDF, and this happens
    // *only* on the iBAQ path (`ibaq.py:969`); the generic / LFQ branches in
    // `commands/peptides2protein.py` never touch `verbose`/`qc_report`. The QC
    // plotting itself moved to the Python wheel, so the iBAQ path prints a pointer
    // to the wheel command instead of drawing a PDF (see `run_ibaq`); `--verbose`
    // on any other method stays the same no-op it is in Python.

    match method.as_str() {
        "ibaq" => run_ibaq(args, output),
        "sum" | "top3" | "topn" => run_generic(args, &method, output),
        "maxlfq" | "directlfq" => run_lfq(args, &method, output),
        other => Err(MokumeError::InvalidInput {
            message: format!("unknown peptides2protein method '{other}'"),
        }),
    }
}

/// DirectLFQ default `num_samples_quadratic` (the global-stage knob). The Python
/// `DirectLFQQuantification` uses directlfq's default of 50; `peptides2protein`
/// does not expose it, so it is fixed here too.
const DIRECTLFQ_NUM_SAMPLES_QUADRATIC: usize = 50;

/// DirectLFQ / MaxLFQ path: roll the peptide matrix up with the DirectLFQ
/// estimator (Python's `DirectLFQQuantification`; `maxlfq` delegates to it with
/// `min_nonan = 2`). Emits the same long-format table as the deterministic
/// methods, keeping only `Intensity > 0` rows (Python's `_parse_wide_output`).
fn run_lfq(args: &Peptides2ProteinArgs, method: &str, output: &Path) -> Result<()> {
    let table = load_peptide_table(&args.peptides)?;
    if !table.has_peptide {
        return Err(MokumeError::InvalidInput {
            message: format!(
                "{method} requires a peptide column ('{PEPTIDE_CANONICAL}' or '{PEPTIDE_SEQUENCE}')"
            ),
        });
    }

    // Python's maxlfq delegates to DirectLFQ with min_nonan = 2 (its min_peptides);
    // the directlfq method uses the configured --min_nonan.
    let min_nonan = if method == "maxlfq" {
        2
    } else {
        args.min_nonan
    };

    let mut condition_by_sample: HashMap<String, String> = HashMap::new();
    let mut observations = Vec::with_capacity(table.rows.len());
    for row in &table.rows {
        condition_by_sample
            .entry(row.sample.clone())
            .or_insert_with(|| row.condition.clone());
        observations.push(LfqPeptideObservation {
            protein: row.protein.clone(),
            peptide: row.peptide.clone(),
            sample: row.sample.clone(),
            intensity: row.intensity,
        });
    }

    let mut results: Vec<GenericRow> =
        run_lfq_from_peptides(&observations, min_nonan, DIRECTLFQ_NUM_SAMPLES_QUADRATIC)
            .into_iter()
            .map(|row| GenericRow {
                protein: row.protein,
                sample: row.sample,
                intensity: row.intensity,
            })
            .collect();
    results.sort_by(|left, right| {
        left.protein
            .cmp(&right.protein)
            .then_with(|| left.sample.cmp(&right.sample))
    });

    // `--normalize` divides each protein's intensity by its sample total, exactly
    // as the Python generic path does after `quantify`.
    let norm = if args.normalize {
        let mut totals: HashMap<&str, f64> = HashMap::new();
        for row in &results {
            *totals.entry(row.sample.as_str()).or_insert(0.0) += row.intensity;
        }
        Some(totals)
    } else {
        None
    };

    write_generic_output(
        output,
        &results,
        &condition_by_sample,
        table.has_condition,
        norm.as_ref(),
    )
}

/// iBAQ path: reuse the existing piBAQ core; emit the Python long-format table
/// plus any requested extras (`--tpa`, `--normalize`, `--ruler`).
///
/// For an enzyme whose cleavage rule is not ported to Rust (anything outside the
/// set [`supports_ibaq_enzyme`] accepts), the Rust kernel cannot digest the
/// proteome, so iBAQ for that enzyme is not supported here: the command fails with
/// a clear `InvalidInput` that points to the Python wheel. The `--fasta`
/// precondition is still enforced first, matching the native path. The default
/// `Trypsin` path and every other ported enzyme stay on the native Rust kernel
/// below.
fn run_ibaq(args: &Peptides2ProteinArgs, output: &Path) -> Result<()> {
    // Resolve the organism up front (matching Python `OrganismDescription.get`,
    // which raises when the name is unknown). The default is `human`.
    let organism = if args.organism.is_empty() {
        None
    } else {
        Some(ibaq_extras::resolve_organism(&args.organism)?)
    };

    // The proteomic ruler requires TPA + non-zero ploidy/cpc + organism, exactly
    // as the Python `peptides_to_protein` guard demands.
    if args.ruler && (!args.tpa || args.ploidy == 0 || args.cpc == 0.0 || args.organism.is_empty())
    {
        return Err(MokumeError::InvalidInput {
            message:
                "`ploidy`, `cpc`, `organism` and `tpa` are required to calculate protein weight and concentration"
                    .to_owned(),
        });
    }

    let Some(fasta) = args.fasta.as_ref() else {
        return Err(MokumeError::InvalidInput {
            message: "the --fasta option is required for the iBAQ method".to_owned(),
        });
    };
    if !fasta.exists() {
        return Err(MokumeError::MissingInput {
            path: fasta.clone(),
        });
    }

    // Enzymes whose cleavage rule is not ported to the Rust kernel cannot be
    // digested here. Point the user to the Python wheel rather than silently
    // producing wrong numbers. The `--fasta` precondition above runs first so an
    // unported enzyme without a FASTA fails on the missing FASTA, as before.
    if !supports_ibaq_enzyme(&args.enzyme) {
        return Err(MokumeError::InvalidInput {
            message: format!(
                "iBAQ digestion enzyme '{}' is not ported to the Rust kernel. \
Install the Rust wheel and run it there: pip install mokume-rs[ibaq]; \
python -m mokume.commands.peptides2protein_ibaq --enzyme '{}' ...",
                args.enzyme, args.enzyme
            ),
        });
    }

    let table = load_peptide_table(&args.peptides)?;
    if !table.has_peptide {
        return Err(MokumeError::InvalidInput {
            message: format!(
                "iBAQ requires a peptide column ('{PEPTIDE_CANONICAL}' or '{PEPTIDE_SEQUENCE}')"
            ),
        });
    }

    // The piBAQ core keys on (peptide, sample); carry a condition lookup so the
    // long-format output can report each sample's condition.
    let mut condition_by_sample: HashMap<String, String> = HashMap::new();
    let mut observations = Vec::with_capacity(table.rows.len());
    for row in &table.rows {
        condition_by_sample
            .entry(row.sample.clone())
            .or_insert_with(|| row.condition.clone());
        observations.push(PeptideObservation {
            peptide: row.peptide.clone(),
            sample: row.sample.clone(),
            intensity: row.intensity,
        });
    }

    let params = IbaqFromPeptidesParams {
        fasta: fasta.clone(),
        min_aa: args.min_aa,
        max_aa: args.max_aa,
        min_shared: args.min_shared,
        min_anchors: args.min_anchors,
        high_anchor_threshold: args.high_anchor_threshold,
        families_yaml: args.families_yaml.clone(),
        enzyme: args.enzyme.clone(),
        tpa: args.tpa,
    };
    let rows = run_ibaq_from_peptides(&observations, &params)?;

    // Lift the piBAQ rows into the extra-aware records, attach the requested
    // post-processing columns, then emit. The piBAQ allocation math is never
    // re-touched: the extras only read `NormIntensity`, `Ibaq`, and (for
    // TPA/ruler) the molecular weight the core already computed.
    let mut records = ibaq_extras::IbaqExtraRow::lift(rows, &condition_by_sample);

    if args.normalize {
        ibaq_extras::normalize_ibaq(&mut records);
    }
    if args.ruler {
        let Some(organism) = organism.as_ref() else {
            return Err(MokumeError::InvalidInput {
                message: "the --organism option is required for the proteomic ruler".to_owned(),
            });
        };
        ibaq_extras::apply_ruler(&mut records, organism, args.ploidy, args.cpc);
    }

    // Stable, deterministic output ordering: protein then sample.
    records.sort_by(|left, right| {
        left.protein
            .cmp(&right.protein)
            .then_with(|| left.sample.cmp(&right.sample))
    });

    write_ibaq_output(
        output,
        &records,
        table.has_condition,
        args.tpa,
        args.normalize,
        args.ruler,
    )?;

    // The Rust kernel owns the numbers (the table just written). The QC report
    // PDF is plotting periphery and now lives in the Python wheel, so `--verbose`
    // writes no PDF here; it only prints a one-line pointer to the wheel command
    // that draws the same density + box plots from the table.
    if args.verbose {
        eprintln!(
            "note: QC report generation moved to the Python wheel: \
pip install mokume-rs[plotting]; \
mokume.peptides2protein_qc(protein_table=\"{}\", qc_report=\"{}\")",
            output.display(),
            args.qc_report.display()
        );
    }

    Ok(())
}

/// Generic deterministic path (`sum` / `top3` / `topn`): faithful re-creation of
/// the Python `quantify` group-bys.
fn run_generic(args: &Peptides2ProteinArgs, method: &str, output: &Path) -> Result<()> {
    let topn = match method {
        "sum" => None,
        "top3" => Some(3),
        "topn" => Some(args.topn_n),
        _ => unreachable!("caller restricts the method set"),
    };

    let table = load_peptide_table(&args.peptides)?;

    // Group raw rows by (protein, sample); Python does not collapse charges or
    // peptidoforms here, so neither do we.
    let mut groups: HashMap<(String, String), Vec<f64>> = HashMap::new();
    let mut condition_by_sample: HashMap<String, String> = HashMap::new();
    for row in &table.rows {
        condition_by_sample
            .entry(row.sample.clone())
            .or_insert_with(|| row.condition.clone());
        groups
            .entry((row.protein.clone(), row.sample.clone()))
            .or_default()
            .push(row.intensity);
    }

    let mut results: Vec<GenericRow> = groups
        .into_iter()
        .map(|((protein, sample), intensities)| {
            let intensity = aggregate(intensities, topn);
            GenericRow {
                protein,
                sample,
                intensity,
            }
        })
        .collect();
    results.sort_by(|left, right| {
        left.protein
            .cmp(&right.protein)
            .then_with(|| left.sample.cmp(&right.sample))
    });

    // `--normalize` divides each protein's intensity by its sample total, the
    // generic-path transform in the Python command (the iBAQ path is deferred).
    let norm = if args.normalize {
        let mut totals: HashMap<&str, f64> = HashMap::new();
        for row in &results {
            *totals.entry(row.sample.as_str()).or_insert(0.0) += row.intensity;
        }
        Some(totals)
    } else {
        None
    };

    write_generic_output(
        output,
        &results,
        &condition_by_sample,
        table.has_condition,
        norm.as_ref(),
    )
}

/// Aggregate one (protein, sample) group: sum (when `topn` is `None`) or the
/// mean of the `topn` most intense rows. Mirrors Python's `groupby.sum` and
/// `sort_values(desc).head(n).mean`.
fn aggregate(mut intensities: Vec<f64>, topn: Option<usize>) -> f64 {
    match topn {
        None => intensities.iter().sum(),
        Some(n) => {
            intensities.sort_by(|left, right| right.total_cmp(left));
            let selected = &intensities[..intensities.len().min(n.max(1))];
            if selected.is_empty() {
                0.0
            } else {
                selected.iter().sum::<f64>() / selected.len() as f64
            }
        }
    }
}

struct GenericRow {
    protein: String,
    sample: String,
    intensity: f64,
}

/// Read a peptide-level table from CSV / `.tsv` / parquet, dispatching on the
/// parquet magic bytes (`PAR1`) exactly like Python's `is_parquet`, so the same
/// `--peptides` path resolves to the same loader regardless of file extension.
fn load_peptide_table(path: &Path) -> Result<PeptideTable> {
    if looks_like_parquet(path)? {
        load_peptide_table_parquet(path)
    } else {
        load_peptide_table_csv(path)
    }
}

/// Detect a parquet file by its leading `PAR1` magic bytes, mirroring Python's
/// `mokume.core.constants.is_parquet` (which reads the first four bytes). A file
/// shorter than four bytes is treated as non-parquet, leaving the CSV/TSV path
/// to surface any parse error.
fn looks_like_parquet(path: &Path) -> Result<bool> {
    let mut file = File::open(path).map_err(|source| MokumeError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let mut header = [0_u8; 4];
    match file.read_exact(&mut header) {
        Ok(()) => Ok(&header == b"PAR1"),
        Err(error) if error.kind() == ErrorKind::UnexpectedEof => Ok(false),
        Err(source) => Err(MokumeError::Io {
            path: path.to_path_buf(),
            source,
        }),
    }
}

/// Load a peptide table from a parquet file, applying the same intensity filter
/// (`dropna` + `> 0`) and `"Empty"` condition default as the CSV path so both
/// inputs yield identical [`PeptideTable`]s. Mirrors Python's `pd.read_parquet`
/// followed by the generic / iBAQ numeric coercion.
fn load_peptide_table_parquet(path: &Path) -> Result<PeptideTable> {
    let raw = read_peptide_parquet(path)?;
    let mut rows = Vec::with_capacity(raw.rows.len());
    for row in raw.rows {
        let Some(intensity) = row.intensity else {
            continue;
        };
        if !intensity.is_finite() || intensity <= 0.0 {
            continue;
        }
        let condition = if raw.has_condition {
            row.condition.unwrap_or_default()
        } else {
            "Empty".to_owned()
        };
        rows.push(PeptideRow {
            protein: row.protein,
            sample: row.sample,
            peptide: row.peptide.unwrap_or_default(),
            condition,
            intensity,
        });
    }

    Ok(PeptideTable {
        rows,
        has_condition: raw.has_condition,
        has_peptide: raw.has_peptide,
    })
}

/// Read a peptide-level table from CSV (comma) or `.tsv` (tab). Rows with a
/// missing or non-positive intensity are dropped, matching the Python
/// `dropna` + `> 0` filter applied on the iBAQ path and the implicit numeric
/// coercion on the generic path.
fn load_peptide_table_csv(path: &Path) -> Result<PeptideTable> {
    let delimiter = if path
        .extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| extension.eq_ignore_ascii_case("tsv"))
    {
        b'\t'
    } else {
        b','
    };

    let mut reader = csv::ReaderBuilder::new()
        .delimiter(delimiter)
        .from_path(path)
        .map_err(|source| csv_error(path, source))?;

    let headers = reader
        .headers()
        .map_err(|source| csv_error(path, source))?
        .iter()
        .map(ToOwned::to_owned)
        .collect::<Vec<_>>();

    let protein_index = column_index(&headers, PROTEIN_NAME)?;
    let sample_index = column_index(&headers, SAMPLE_ID)?;
    let intensity_index = column_index(&headers, NORM_INTENSITY)?;
    let condition_index = optional_column_index(&headers, CONDITION);
    let peptide_index = optional_column_index(&headers, PEPTIDE_CANONICAL)
        .or_else(|| optional_column_index(&headers, PEPTIDE_SEQUENCE));

    let mut rows = Vec::new();
    for record in reader.records() {
        let record = record.map_err(|source| csv_error(path, source))?;
        let raw_intensity = field(&record, intensity_index, path)?;
        let trimmed = raw_intensity.trim();
        if trimmed.is_empty() {
            continue;
        }
        let intensity = trimmed
            .parse::<f64>()
            .map_err(|_| MokumeError::InvalidInput {
                message: format!("'{trimmed}' in column '{NORM_INTENSITY}' is not numeric"),
            })?;
        if !intensity.is_finite() || intensity <= 0.0 {
            continue;
        }
        let condition = condition_index
            .map(|index| field(&record, index, path))
            .transpose()?
            .map_or_else(|| "Empty".to_owned(), ToOwned::to_owned);
        let peptide = peptide_index
            .map(|index| field(&record, index, path))
            .transpose()?
            .unwrap_or("")
            .to_owned();
        rows.push(PeptideRow {
            protein: field(&record, protein_index, path)?.to_owned(),
            sample: field(&record, sample_index, path)?.to_owned(),
            peptide,
            condition,
            intensity,
        });
    }

    Ok(PeptideTable {
        rows,
        has_condition: condition_index.is_some(),
        has_peptide: peptide_index.is_some(),
    })
}

/// Write the iBAQ long-format output (tab-separated, matching Python). Columns
/// are appended in the Python `peptides_to_protein` order: the base eight, then
/// the TPA pair, then the normalize triple, then the ruler quartet.
fn write_ibaq_output(
    output: &Path,
    rows: &[ibaq_extras::IbaqExtraRow],
    has_condition: bool,
    tpa: bool,
    normalize: bool,
    ruler: bool,
) -> Result<()> {
    let mut writer = create_writer(output)?;

    let mut header = vec![PROTEIN_NAME, SAMPLE_ID];
    if has_condition {
        header.push(CONDITION);
    }
    header.extend([
        NORM_INTENSITY,
        "Ibaq",
        "FamilyId",
        "EvidenceLevel",
        "FamilySize",
    ]);
    if tpa {
        header.extend([ibaq_extras::MOLECULAR_WEIGHT, ibaq_extras::TPA]);
    }
    if normalize {
        header.extend([
            ibaq_extras::IBAQ_NORM,
            ibaq_extras::IBAQ_LOG,
            ibaq_extras::IBAQ_PPB,
        ]);
    }
    if ruler {
        header.extend([
            ibaq_extras::COPY_NUMBER,
            ibaq_extras::MOLES_NMOL,
            ibaq_extras::WEIGHT_NG,
            ibaq_extras::CONCENTRATION_NM,
        ]);
    }
    write_line(&mut writer, output, &header.join("\t"))?;

    for row in rows {
        let mut fields = vec![row.protein.clone(), row.sample.clone()];
        if has_condition {
            fields.push(row.condition.clone());
        }
        fields.push(format_float(row.norm_intensity));
        fields.push(format_float(row.ibaq));
        fields.push(row.family_id.clone());
        fields.push(row.evidence_level.to_owned());
        fields.push(row.family_size.to_string());
        if tpa {
            fields.push(format_optional(row.molecular_weight));
            fields.push(format_optional(row.tpa));
        }
        if normalize {
            fields.push(format_optional(row.ibaq_norm));
            fields.push(format_optional(row.ibaq_log));
            fields.push(format_optional(row.ibaq_ppb));
        }
        if ruler {
            fields.push(format_optional(row.copy_number));
            fields.push(format_optional(row.moles_nmol));
            fields.push(format_optional(row.weight_ng));
            fields.push(format_optional(row.concentration_nm));
        }
        write_line(&mut writer, output, &fields.join("\t"))?;
    }

    flush(&mut writer, output)
}

/// Write the generic-method long-format output (tab-separated, matching Python).
fn write_generic_output(
    output: &Path,
    rows: &[GenericRow],
    condition_by_sample: &HashMap<String, String>,
    has_condition: bool,
    norm: Option<&HashMap<&str, f64>>,
) -> Result<()> {
    let mut writer = create_writer(output)?;

    let mut header = vec![PROTEIN_NAME, SAMPLE_ID, "Intensity"];
    if has_condition {
        header.push(CONDITION);
    }
    if norm.is_some() {
        header.push("IntensityNorm");
    }
    write_line(&mut writer, output, &header.join("\t"))?;

    for row in rows {
        let mut fields = vec![
            row.protein.clone(),
            row.sample.clone(),
            format_float(row.intensity),
        ];
        if has_condition {
            fields.push(
                condition_by_sample
                    .get(&row.sample)
                    .cloned()
                    .unwrap_or_else(|| "Empty".to_owned()),
            );
        }
        if let Some(totals) = norm {
            let total = totals.get(row.sample.as_str()).copied().unwrap_or(0.0);
            let normalized = if total > 0.0 {
                row.intensity / total
            } else {
                0.0
            };
            fields.push(format_float(normalized));
        }
        write_line(&mut writer, output, &fields.join("\t"))?;
    }

    flush(&mut writer, output)
}

fn create_writer(output: &Path) -> Result<BufWriter<File>> {
    if let Some(parent) = output.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(|source| MokumeError::Io {
                path: parent.to_path_buf(),
                source,
            })?;
        }
    }
    let file = File::create(output).map_err(|source| MokumeError::Io {
        path: output.to_path_buf(),
        source,
    })?;
    Ok(BufWriter::new(file))
}

fn write_line(writer: &mut BufWriter<File>, output: &Path, line: &str) -> Result<()> {
    writeln!(writer, "{line}").map_err(|source| MokumeError::Io {
        path: output.to_path_buf(),
        source,
    })
}

fn flush(writer: &mut BufWriter<File>, output: &Path) -> Result<()> {
    writer.flush().map_err(|source| MokumeError::Io {
        path: output.to_path_buf(),
        source,
    })
}

/// `{}` on `f64` gives Rust's shortest round-trip representation, matching the
/// precision the golden test compares (relative 1e-9).
fn format_float(value: f64) -> String {
    format!("{value}")
}

/// Render an optional extra-column value; absent values become an empty cell
/// (the extras are always populated when their flag is set, so this only
/// guards against a missing post-processing pass).
fn format_optional(value: Option<f64>) -> String {
    value.map(format_float).unwrap_or_default()
}

fn column_index(headers: &[String], column: &str) -> Result<usize> {
    optional_column_index(headers, column).ok_or_else(|| MokumeError::InvalidInput {
        message: format!("could not find required column '{column}' in the peptide file"),
    })
}

fn optional_column_index(headers: &[String], column: &str) -> Option<usize> {
    headers.iter().position(|header| header == column)
}

fn field<'a>(record: &'a csv::StringRecord, index: usize, path: &Path) -> Result<&'a str> {
    record.get(index).ok_or_else(|| MokumeError::InvalidInput {
        message: format!("row in '{}' is missing column {index}", path.display()),
    })
}

fn csv_error(path: &Path, source: csv::Error) -> MokumeError {
    MokumeError::InvalidInput {
        message: format!("error reading file '{}': {source}", path.display()),
    }
}

/// iBAQ post-processing extras: TPA (carried over from the piBAQ core),
/// `normalize_ibaq` (rIBAQ / `IbaqLog` / `IbaqPpb`), and the ProteomicRuler
/// copy-number / moles / weight / concentration columns.
///
/// Every transform is a deterministic function of the piBAQ `NormIntensity` /
/// `Ibaq` (and, for TPA/ruler, the `MolecularWeight`) the verified core already
/// produced -- the allocation math is never re-touched. Each formula mirrors
/// `mokume.quantification.ibaq` exactly so the output cells match the Python
/// oracle to 1e-9.
mod ibaq_extras {
    use std::collections::HashMap;

    use mokume_core::{MokumeError, Result};
    use serde::Deserialize;

    /// `MolecularWeight` column header (Python `MOLECULARWEIGHT`).
    pub(super) const MOLECULAR_WEIGHT: &str = "MolecularWeight";
    /// `TPA` column header (Python `TPA`).
    pub(super) const TPA: &str = "TPA";
    /// Relative iBAQ column header (Python `IBAQ_NORMALIZED`).
    pub(super) const IBAQ_NORM: &str = "IbaqNorm";
    /// Log-shifted iBAQ column header (Python `IBAQ_LOG`).
    pub(super) const IBAQ_LOG: &str = "IbaqLog";
    /// Parts-per-billion iBAQ column header (Python `IBAQ_PPB`).
    pub(super) const IBAQ_PPB: &str = "IbaqPpb";
    /// Copy-number column header (Python `COPYNUMBER`).
    pub(super) const COPY_NUMBER: &str = "CopyNumber";
    /// Moles column header (Python `MOLES_NMOL`).
    pub(super) const MOLES_NMOL: &str = "Moles[nmol]";
    /// Weight column header (Python `WEIGHT_NG`).
    pub(super) const WEIGHT_NG: &str = "Weight[ng]";
    /// Concentration column header (Python `CONCENTRATION_NM`).
    pub(super) const CONCENTRATION_NM: &str = "Concentration[nM]";

    /// Avogadro's number, the Python `AVAGADRO` constant (kept spelled as the
    /// physical constant; the Python source has the historical misspelling).
    const AVOGADRO: f64 = 6.02214129e23;
    /// Average base-pair mass, the Python `AVERAGE_BASE_PAIR_MASS` constant.
    const AVERAGE_BASE_PAIR_MASS: f64 = 617.96;
    /// PRIDE parts-per-billion multiplier (`* 100_000_000`).
    const IBAQ_PPB_FACTOR: f64 = 100_000_000.0;
    /// ProteomicsDB log shift (`log10(rIBAQ) + 10`).
    const IBAQ_LOG_SHIFT: f64 = 10.0;
    /// Nanomole-per-mole scale (`1e9 / AVOGADRO`) used for `Moles[nmol]`.
    const NMOL_PER_MOLE: f64 = 1e9;
    /// Nanogram divisor (`/ 1e-9`) used when deriving the ruler volume.
    const WEIGHT_NG_TO_GRAMS: f64 = 1e-9;

    /// Embedded copy of `mokume/data/organisms.json`, parsed at runtime to look
    /// up genome size and histone entries for the proteomic ruler.
    const ORGANISMS_JSON: &str = include_str!("data/organisms.json");

    /// One iBAQ output row, extended with the optional extra columns. The base
    /// columns mirror [`mokume_pipeline::IbaqProteinRow`]; the `Option` fields
    /// are populated by the corresponding post-processing pass.
    #[derive(Debug, Clone)]
    pub(super) struct IbaqExtraRow {
        pub(super) protein: String,
        pub(super) sample: String,
        pub(super) condition: String,
        pub(super) norm_intensity: f64,
        pub(super) ibaq: f64,
        pub(super) family_id: String,
        pub(super) evidence_level: &'static str,
        pub(super) family_size: usize,
        pub(super) molecular_weight: Option<f64>,
        pub(super) tpa: Option<f64>,
        pub(super) ibaq_norm: Option<f64>,
        pub(super) ibaq_log: Option<f64>,
        pub(super) ibaq_ppb: Option<f64>,
        pub(super) copy_number: Option<f64>,
        pub(super) moles_nmol: Option<f64>,
        pub(super) weight_ng: Option<f64>,
        pub(super) concentration_nm: Option<f64>,
    }

    impl IbaqExtraRow {
        /// Lift the piBAQ core rows into extra-aware records, resolving each
        /// sample's condition (defaulting to `Empty`, as the Python loader fills
        /// a missing `Condition`).
        pub(super) fn lift(
            rows: Vec<mokume_pipeline::IbaqProteinRow>,
            condition_by_sample: &HashMap<String, String>,
        ) -> Vec<Self> {
            rows.into_iter()
                .map(|row| {
                    let condition = condition_by_sample
                        .get(&row.sample)
                        .cloned()
                        .unwrap_or_else(|| "Empty".to_owned());
                    Self {
                        protein: row.protein,
                        sample: row.sample,
                        condition,
                        norm_intensity: row.norm_intensity,
                        ibaq: row.ibaq,
                        family_id: row.family_id,
                        evidence_level: row.evidence_level,
                        family_size: row.family_size,
                        molecular_weight: row.molecular_weight,
                        tpa: row.tpa,
                        ibaq_norm: None,
                        ibaq_log: None,
                        ibaq_ppb: None,
                        copy_number: None,
                        moles_nmol: None,
                        weight_ng: None,
                        concentration_nm: None,
                    }
                })
                .collect()
        }
    }

    /// `normalize_ibaq`: per (SampleID, Condition) group, divide each protein's
    /// iBAQ by the group's total iBAQ (rIBAQ), then derive the ProteomicsDB log
    /// (`10 + log10(rIBAQ)` where positive, else 0) and the PRIDE parts-per-
    /// billion (`rIBAQ * 1e8`). Mirrors `mokume.quantification.ibaq.normalize_ibaq`.
    pub(super) fn normalize_ibaq(rows: &mut [IbaqExtraRow]) {
        let mut totals: HashMap<(&str, &str), f64> = HashMap::new();
        for row in rows.iter() {
            *totals
                .entry((row.sample.as_str(), row.condition.as_str()))
                .or_insert(0.0) += row.ibaq;
        }
        // Snapshot the totals so the borrow ends before the mutable pass.
        let totals: HashMap<(String, String), f64> = totals
            .into_iter()
            .map(|((sample, condition), total)| ((sample.to_owned(), condition.to_owned()), total))
            .collect();
        for row in rows.iter_mut() {
            let total = totals
                .get(&(row.sample.clone(), row.condition.clone()))
                .copied()
                .unwrap_or(0.0);
            let norm = row.ibaq / total;
            row.ibaq_norm = Some(norm);
            row.ibaq_log = Some(if norm > 0.0 {
                norm.log10() + IBAQ_LOG_SHIFT
            } else {
                0.0
            });
            row.ibaq_ppb = Some(norm * IBAQ_PPB_FACTOR);
        }
    }

    /// ProteomicRuler description for a single organism, deserialized from the
    /// embedded `organisms.json`. Only the fields the ruler reads are kept.
    #[derive(Debug, Clone, Deserialize)]
    pub(super) struct Organism {
        name: String,
        genome_size: f64,
        #[serde(default)]
        histone_entries: Vec<String>,
    }

    /// Resolve an organism name case-insensitively, matching the Python
    /// `OrganismDescription.get(name.upper())`. Returns an error for an unknown
    /// organism, mirroring the Python `KeyError`.
    pub(super) fn resolve_organism(name: &str) -> Result<Organism> {
        let registry: HashMap<String, Organism> =
            serde_json::from_str(ORGANISMS_JSON).map_err(|source| MokumeError::InvalidInput {
                message: format!("could not parse the embedded organism table: {source}"),
            })?;
        let key = name.to_ascii_uppercase();
        if let Some(organism) = registry.get(&key) {
            return Ok(organism.clone());
        }
        // The registry is keyed by the uppercased `name` field in Python; the
        // JSON keys already match, but fall back to a name scan for robustness.
        registry
            .into_values()
            .find(|organism| organism.name.eq_ignore_ascii_case(name))
            .ok_or_else(|| MokumeError::InvalidInput {
                message: format!("could not resolve organism description for {name}"),
            })
    }

    /// Apply the proteomic ruler per Condition group, mirroring
    /// `ConcentrationWeightByProteomicRuler.apply_by_condition`. For each
    /// condition independently: histone intensity = max(sum of `NormIntensity`
    /// over histone proteins, 1.0); then per protein
    /// `CopyNumber = NormIntensity / histone * dna_mass * AVOGADRO / MW`,
    /// `Moles[nmol] = CopyNumber * 1e9 / AVOGADRO`,
    /// `Weight[ng] = Moles[nmol] * MW`; finally the group volume
    /// `= sum(Weight[ng]) / 1e-9 / cpc` and
    /// `Concentration[nM] = volume * Moles[nmol]`.
    pub(super) fn apply_ruler(
        rows: &mut [IbaqExtraRow],
        organism: &Organism,
        ploidy: i32,
        cpc: f64,
    ) {
        let dna_mass = f64::from(ploidy) * organism.genome_size * AVERAGE_BASE_PAIR_MASS / AVOGADRO;
        let histones: std::collections::HashSet<&str> = organism
            .histone_entries
            .iter()
            .map(String::as_str)
            .collect();

        // Group row indices by condition so each group's histone intensity and
        // volume are computed over exactly that condition's proteins.
        let mut groups: HashMap<String, Vec<usize>> = HashMap::new();
        for (index, row) in rows.iter().enumerate() {
            groups.entry(row.condition.clone()).or_default().push(index);
        }

        for indices in groups.values() {
            let histone_intensity = indices
                .iter()
                .filter(|&&index| histones.contains(rows[index].protein.as_str()))
                .map(|&index| rows[index].norm_intensity)
                .sum::<f64>()
                .max(1.0);

            // First pass: copy number, moles, weight; accumulate the group's
            // total weight for the volume term.
            let mut total_weight = 0.0;
            for &index in indices {
                let row = &mut rows[index];
                let mw = row.molecular_weight.unwrap_or(1.0);
                let copy_number = row.norm_intensity / histone_intensity * dna_mass * AVOGADRO / mw;
                let moles = copy_number * (NMOL_PER_MOLE / AVOGADRO);
                let weight = moles * mw;
                row.copy_number = Some(copy_number);
                row.moles_nmol = Some(moles);
                row.weight_ng = Some(weight);
                total_weight += weight;
            }

            // Second pass: the group volume depends on the summed weight, so the
            // concentration is derived only once every weight is known.
            let volume = total_weight / WEIGHT_NG_TO_GRAMS / cpc;
            for &index in indices {
                let row = &mut rows[index];
                let moles = row.moles_nmol.unwrap_or(0.0);
                row.concentration_nm = Some(volume * moles);
            }
        }
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

    fn temp_dir(tag: &str) -> TestResult<std::path::PathBuf> {
        let nanos = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
        Ok(tempfile::Builder::new()
            .prefix(&format!("mokume-peptides2protein-{tag}-{nanos}-"))
            .tempdir()?
            .keep())
    }

    fn write_file(path: &Path, contents: &str) -> TestResult<()> {
        let mut file = File::create(path)?;
        file.write_all(contents.as_bytes())?;
        Ok(())
    }

    // Synthetic peptide matrix shared by the golden tests. Each
    // (protein, sample, peptide) is unique so the Python generic group-bys and
    // the Rust path coincide exactly. P3's `THIDPECK` is intentionally absent
    // from the FASTA digest (the digest yields `ATHIDPECK`), exercising the
    // FASTA-driven peptide filtering on the iBAQ path.
    const PEPTIDES_CSV: &str = "ProteinName,PeptideCanonical,SampleID,Condition,NormIntensity\n\
P1,PEPTIDEAK,S1,A,100.0\n\
P1,APEPTIDECK,S1,A,300.0\n\
P1,ASHAEDPEPK,S1,A,50.0\n\
P1,PEPTIDEAK,S2,B,40.0\n\
P1,APEPTIDECK,S2,B,60.0\n\
P2,ALYAAEK,S1,A,500.0\n\
P3,THIDPEAK,S1,A,700.0\n\
P3,THIDPECK,S1,A,900.0\n";

    const PROTEOME_FASTA: &str =
        ">P1\nPEPTIDEAKAPEPTIDECKASHAEDPEPK\n>P2\nALYAAEK\n>P3\nTHIDPEAKATHIDPECK\n";

    fn base_args(peptides: &Path, output: &Path) -> Peptides2ProteinArgs {
        Peptides2ProteinArgs {
            fasta: None,
            peptides: peptides.to_path_buf(),
            method: "ibaq".to_owned(),
            enzyme: "Trypsin".to_owned(),
            normalize: false,
            min_aa: 7,
            max_aa: 30,
            tpa: false,
            ruler: false,
            ploidy: 2,
            organism: "human".to_owned(),
            cpc: 200.0,
            output: Some(output.to_path_buf()),
            verbose: false,
            qc_report: std::path::PathBuf::from("QCprofile.pdf"),
            topn_n: 3,
            threads: -1,
            min_nonan: 1,
            families_yaml: None,
            min_shared: 2,
            min_anchors: 1,
            high_anchor_threshold: 3,
        }
    }

    fn read_table(path: &Path) -> TestResult<(Vec<String>, Vec<Vec<String>>)> {
        let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?;
        let headers = reader.headers()?.iter().map(ToOwned::to_owned).collect();
        let rows = reader
            .records()
            .map(|record| record.map(|record| record.iter().map(ToOwned::to_owned).collect()))
            .collect::<std::result::Result<Vec<Vec<String>>, _>>()?;
        Ok((headers, rows))
    }

    /// Look up the cell value at (protein, sample, column), returning an error
    /// (never panicking via unwrap/expect) when the row or column is absent.
    fn cell(
        headers: &[String],
        rows: &[Vec<String>],
        protein: &str,
        sample: &str,
        column: &str,
    ) -> TestResult<String> {
        let position = |name: &str| {
            headers
                .iter()
                .position(|header| header == name)
                .ok_or_else(|| -> Box<dyn Error> { format!("column '{name}' is missing").into() })
        };
        let protein_idx = position(PROTEIN_NAME)?;
        let sample_idx = position(SAMPLE_ID)?;
        let column_idx = position(column)?;
        let row = rows
            .iter()
            .find(|row| {
                row.get(protein_idx).is_some_and(|p| p == protein)
                    && row.get(sample_idx).is_some_and(|s| s == sample)
            })
            .ok_or_else(|| -> Box<dyn Error> {
                format!("row ({protein}, {sample}) is missing").into()
            })?;
        row.get(column_idx)
            .cloned()
            .ok_or_else(|| format!("cell ({protein}, {sample}, {column}) is missing").into())
    }

    /// Assert the numeric value at (protein, sample, column) is within 1e-9 of
    /// `expected`, without unwrap/expect (forbidden by the lint policy).
    fn assert_cell_close(
        headers: &[String],
        rows: &[Vec<String>],
        protein: &str,
        sample: &str,
        column: &str,
        expected: f64,
    ) -> TestResult<()> {
        let value = cell(headers, rows, protein, sample, column)?;
        let got = value
            .parse::<f64>()
            .map_err(|_| -> Box<dyn Error> { format!("'{value}' is not numeric").into() })?;
        assert!(
            (got - expected).abs() <= 1e-9,
            "({protein}, {sample}, {column}): got {got}, expected {expected}"
        );
        Ok(())
    }

    // Golden oracle (Python), captured with:
    //   conda run -n Bigbio python -m mokume.mokume_cli peptides2protein \
    //     --method sum -p peptides.csv -o out.tsv
    // Columns: ProteinName SampleID Intensity Condition.
    #[test]
    fn peptides2protein_sum_matches_python_oracle() -> TestResult<()> {
        let dir = temp_dir("sum")?;
        let peptides = dir.join("peptides.csv");
        let output = dir.join("out.tsv");
        write_file(&peptides, PEPTIDES_CSV)?;

        let mut args = base_args(&peptides, &output);
        args.method = "sum".to_owned();
        run_peptides_to_protein(&args)?;

        let (headers, rows) = read_table(&output)?;
        assert_eq!(
            headers,
            vec!["ProteinName", "SampleID", "Intensity", "Condition"]
        );
        assert_cell_close(&headers, &rows, "P1", "S1", "Intensity", 450.0)?;
        assert_cell_close(&headers, &rows, "P1", "S2", "Intensity", 100.0)?;
        assert_cell_close(&headers, &rows, "P2", "S1", "Intensity", 500.0)?;
        assert_cell_close(&headers, &rows, "P3", "S1", "Intensity", 1600.0)?;
        assert_eq!(cell(&headers, &rows, "P1", "S2", "Condition")?, "B");
        Ok(())
    }

    /// Build the parquet form of [`PEPTIDES_CSV`] via the shared
    /// `write_peptide_parquet` schema (dictionary-encoded `SampleID`/`Condition`,
    /// float32 `NormIntensity`), so the parity test exercises the real reader.
    fn peptide_parquet_rows_from_csv(csv: &str) -> TestResult<Vec<mokume_io::PeptideParquetRow>> {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b',')
            .from_reader(csv.as_bytes());
        let mut rows = Vec::new();
        for record in reader.records() {
            let record = record?;
            let get = |index: usize| -> TestResult<String> {
                record
                    .get(index)
                    .map(ToOwned::to_owned)
                    .ok_or_else(|| -> Box<dyn Error> { format!("missing field {index}").into() })
            };
            let intensity = get(4)?
                .parse::<f64>()
                .map_err(|_| -> Box<dyn Error> { "non-numeric intensity".into() })?;
            rows.push(mokume_io::PeptideParquetRow {
                protein_name: get(0)?,
                peptide_canonical: get(1)?,
                sample_id: get(2)?,
                bio_replicate: 1,
                condition: get(3)?,
                run: None,
                tech_replicate: None,
                norm_intensity: intensity,
            });
        }
        Ok(rows)
    }

    // The parquet and CSV loaders must produce byte-identical protein output for
    // the same data (Python reads either via `pd.read_parquet`/`pd.read_csv`).
    #[test]
    fn peptides2protein_parquet_input_matches_csv_input() -> TestResult<()> {
        let dir = temp_dir("parquet-parity")?;
        let csv_path = dir.join("peptides.csv");
        let parquet_path = dir.join("peptides.parquet");
        write_file(&csv_path, PEPTIDES_CSV)?;
        let parquet_rows = peptide_parquet_rows_from_csv(PEPTIDES_CSV)?;
        mokume_io::write_peptide_parquet(&parquet_path, &parquet_rows, false)?;

        let csv_out = dir.join("from_csv.tsv");
        let parquet_out = dir.join("from_parquet.tsv");

        let mut csv_args = base_args(&csv_path, &csv_out);
        csv_args.method = "sum".to_owned();
        run_peptides_to_protein(&csv_args)?;

        let mut parquet_args = base_args(&parquet_path, &parquet_out);
        parquet_args.method = "sum".to_owned();
        run_peptides_to_protein(&parquet_args)?;

        let (csv_headers, csv_rows) = read_table(&csv_out)?;
        let (parquet_headers, parquet_rows_out) = read_table(&parquet_out)?;
        assert_eq!(csv_headers, parquet_headers);
        assert_eq!(csv_rows, parquet_rows_out);
        Ok(())
    }

    // Oracle:
    //   ... peptides2protein --method topn --topn_n 2 -p peptides.csv -o out.tsv
    #[test]
    fn peptides2protein_topn_matches_python_oracle() -> TestResult<()> {
        let dir = temp_dir("topn")?;
        let peptides = dir.join("peptides.csv");
        let output = dir.join("out.tsv");
        write_file(&peptides, PEPTIDES_CSV)?;

        let mut args = base_args(&peptides, &output);
        args.method = "topn".to_owned();
        args.topn_n = 2;
        run_peptides_to_protein(&args)?;

        let (headers, rows) = read_table(&output)?;
        assert_cell_close(&headers, &rows, "P1", "S1", "Intensity", 200.0)?;
        assert_cell_close(&headers, &rows, "P1", "S2", "Intensity", 50.0)?;
        assert_cell_close(&headers, &rows, "P2", "S1", "Intensity", 500.0)?;
        assert_cell_close(&headers, &rows, "P3", "S1", "Intensity", 800.0)?;
        Ok(())
    }

    // Oracle:
    //   ... peptides2protein --method top3 -p peptides.csv -o out.tsv
    #[test]
    fn peptides2protein_top3_matches_python_oracle() -> TestResult<()> {
        let dir = temp_dir("top3")?;
        let peptides = dir.join("peptides.csv");
        let output = dir.join("out.tsv");
        write_file(&peptides, PEPTIDES_CSV)?;

        let mut args = base_args(&peptides, &output);
        args.method = "top3".to_owned();
        run_peptides_to_protein(&args)?;

        let (headers, rows) = read_table(&output)?;
        assert_cell_close(&headers, &rows, "P1", "S1", "Intensity", 150.0)?;
        assert_cell_close(&headers, &rows, "P3", "S1", "Intensity", 800.0)?;
        Ok(())
    }

    // Oracle:
    //   ... peptides2protein --method sum -n -p peptides.csv -o out.tsv
    // Adds IntensityNorm = Intensity / sum(Intensity per SampleID).
    #[test]
    fn peptides2protein_sum_normalize_matches_python_oracle() -> TestResult<()> {
        let dir = temp_dir("sumnorm")?;
        let peptides = dir.join("peptides.csv");
        let output = dir.join("out.tsv");
        write_file(&peptides, PEPTIDES_CSV)?;

        let mut args = base_args(&peptides, &output);
        args.method = "sum".to_owned();
        args.normalize = true;
        run_peptides_to_protein(&args)?;

        let (headers, rows) = read_table(&output)?;
        assert!(headers.contains(&"IntensityNorm".to_owned()));
        // S1 total = 450 + 500 + 1600 = 2550; S2 total = 100.
        assert_cell_close(&headers, &rows, "P1", "S1", "IntensityNorm", 450.0 / 2550.0)?;
        assert_cell_close(
            &headers,
            &rows,
            "P3",
            "S1",
            "IntensityNorm",
            1600.0 / 2550.0,
        )?;
        assert_cell_close(&headers, &rows, "P1", "S2", "IntensityNorm", 1.0)?;
        Ok(())
    }

    // Oracle:
    //   ... peptides2protein --method ibaq -f proteome.fasta -p peptides.csv -o out.tsv
    // Columns: ProteinName SampleID Condition NormIntensity Ibaq FamilyId EvidenceLevel FamilySize.
    #[test]
    fn peptides2protein_ibaq_matches_python_oracle() -> TestResult<()> {
        let dir = temp_dir("ibaq")?;
        let peptides = dir.join("peptides.csv");
        let fasta = dir.join("proteome.fasta");
        let output = dir.join("out.tsv");
        write_file(&peptides, PEPTIDES_CSV)?;
        write_file(&fasta, PROTEOME_FASTA)?;

        let mut args = base_args(&peptides, &output);
        args.method = "ibaq".to_owned();
        args.fasta = Some(fasta);
        run_peptides_to_protein(&args)?;

        let (headers, rows) = read_table(&output)?;
        assert_eq!(
            headers,
            vec![
                "ProteinName",
                "SampleID",
                "Condition",
                "NormIntensity",
                "Ibaq",
                "FamilyId",
                "EvidenceLevel",
                "FamilySize",
            ]
        );
        // P1/S1: NormIntensity 450, Ibaq 450/3 = 150 (3 proteotypic peptides).
        assert_cell_close(&headers, &rows, "P1", "S1", "NormIntensity", 450.0)?;
        assert_cell_close(&headers, &rows, "P1", "S1", "Ibaq", 150.0)?;
        assert_cell_close(&headers, &rows, "P1", "S2", "Ibaq", 100.0 / 3.0)?;
        assert_cell_close(&headers, &rows, "P2", "S1", "Ibaq", 500.0)?;
        // P3/S1: only THIDPEAK is in the digest, so NormIntensity 700, Ibaq 350.
        assert_cell_close(&headers, &rows, "P3", "S1", "NormIntensity", 700.0)?;
        assert_cell_close(&headers, &rows, "P3", "S1", "Ibaq", 350.0)?;
        assert_eq!(cell(&headers, &rows, "P1", "S1", "EvidenceLevel")?, "high");
        assert_eq!(
            cell(&headers, &rows, "P2", "S1", "EvidenceLevel")?,
            "medium"
        );
        Ok(())
    }

    // `--high-anchor-threshold` only re-buckets the `EvidenceLevel` annotation
    // (Python `_classify_evidence`); it must never change the `Ibaq` values.
    // P1 has 3 anchors, P2 has 1 (min_anchors stays 1). Lowering the threshold to
    // 1 promotes P2 medium -> high; raising it to 4 demotes P1 high -> medium.
    #[test]
    fn peptides2protein_ibaq_high_anchor_threshold_only_changes_evidence() -> TestResult<()> {
        let dir = temp_dir("ibaq-threshold")?;
        let peptides = dir.join("peptides.csv");
        let fasta = dir.join("proteome.fasta");
        write_file(&peptides, PEPTIDES_CSV)?;
        write_file(&fasta, PROTEOME_FASTA)?;

        let run = |threshold: usize, out: &Path| -> TestResult<()> {
            let mut args = base_args(&peptides, out);
            args.method = "ibaq".to_owned();
            args.fasta = Some(fasta.clone());
            args.high_anchor_threshold = threshold;
            run_peptides_to_protein(&args)?;
            Ok(())
        };

        let low_out = dir.join("low.tsv");
        let high_out = dir.join("high.tsv");
        run(1, &low_out)?;
        run(4, &high_out)?;

        let (low_h, low_rows) = read_table(&low_out)?;
        let (high_h, high_rows) = read_table(&high_out)?;

        // Evidence re-buckets with the threshold.
        assert_eq!(
            cell(&low_h, &low_rows, "P1", "S1", "EvidenceLevel")?,
            "high"
        );
        assert_eq!(
            cell(&low_h, &low_rows, "P2", "S1", "EvidenceLevel")?,
            "high"
        );
        assert_eq!(
            cell(&high_h, &high_rows, "P1", "S1", "EvidenceLevel")?,
            "medium"
        );
        assert_eq!(
            cell(&high_h, &high_rows, "P2", "S1", "EvidenceLevel")?,
            "medium"
        );

        // Quantification values are invariant across thresholds.
        for (protein, sample, expected) in [
            ("P1", "S1", 150.0),
            ("P2", "S1", 500.0),
            ("P3", "S1", 350.0),
        ] {
            assert_cell_close(&low_h, &low_rows, protein, sample, "Ibaq", expected)?;
            assert_cell_close(&high_h, &high_rows, protein, sample, "Ibaq", expected)?;
        }
        Ok(())
    }

    #[test]
    fn peptides2protein_ibaq_requires_fasta() -> TestResult<()> {
        let dir = temp_dir("nofasta")?;
        let peptides = dir.join("peptides.csv");
        let output = dir.join("out.tsv");
        write_file(&peptides, PEPTIDES_CSV)?;

        let args = base_args(&peptides, &output);
        let Err(error) = run_peptides_to_protein(&args) else {
            panic!("iBAQ without --fasta must fail");
        };
        assert!(matches!(error, MokumeError::InvalidInput { .. }));
        Ok(())
    }

    // Oracle (TPA columns MolecularWeight + TPA appended):
    //   ... peptides2protein --method ibaq -f proteome.fasta -p peptides.csv \
    //       --tpa -o out.tsv
    // MolecularWeight is the pyOpenMS getMonoWeight of the canonical protein;
    // TPA = NormIntensity / MolecularWeight.
    #[test]
    fn peptides2protein_ibaq_tpa_matches_python_oracle() -> TestResult<()> {
        let dir = temp_dir("ibaqtpa")?;
        let peptides = dir.join("peptides.csv");
        let fasta = dir.join("proteome.fasta");
        let output = dir.join("out.tsv");
        write_file(&peptides, PEPTIDES_CSV)?;
        write_file(&fasta, PROTEOME_FASTA)?;

        let mut args = base_args(&peptides, &output);
        args.fasta = Some(fasta);
        args.tpa = true;
        run_peptides_to_protein(&args)?;

        let (headers, rows) = read_table(&output)?;
        assert!(headers.contains(&"MolecularWeight".to_owned()));
        assert!(headers.contains(&"TPA".to_owned()));
        // pyOpenMS getMonoWeight: P1=3143.460508429001, P2=764.4068587863999,
        // P3=1903.9098218451002.
        assert_cell_close(
            &headers,
            &rows,
            "P1",
            "S1",
            "MolecularWeight",
            3143.460508429001,
        )?;
        assert_cell_close(
            &headers,
            &rows,
            "P2",
            "S1",
            "MolecularWeight",
            764.4068587863999,
        )?;
        // TPA = NormIntensity / MolecularWeight.
        assert_cell_close(
            &headers,
            &rows,
            "P1",
            "S1",
            "TPA",
            450.0 / 3143.460508429001,
        )?;
        assert_cell_close(
            &headers,
            &rows,
            "P1",
            "S2",
            "TPA",
            100.0 / 3143.460508429001,
        )?;
        assert_cell_close(
            &headers,
            &rows,
            "P2",
            "S1",
            "TPA",
            500.0 / 764.4068587863999,
        )?;
        assert_cell_close(
            &headers,
            &rows,
            "P3",
            "S1",
            "TPA",
            700.0 / 1903.9098218451002,
        )?;
        Ok(())
    }

    // Oracle (normalize_ibaq columns IbaqNorm + IbaqLog + IbaqPpb appended):
    //   ... peptides2protein --method ibaq -f proteome.fasta -p peptides.csv \
    //       -n -o out.tsv
    // rIBAQ = Ibaq / sum(Ibaq per SampleID,Condition);
    // IbaqLog = 10 + log10(rIBAQ); IbaqPpb = rIBAQ * 1e8.
    #[test]
    fn peptides2protein_ibaq_normalize_matches_python_oracle() -> TestResult<()> {
        let dir = temp_dir("ibaqnorm")?;
        let peptides = dir.join("peptides.csv");
        let fasta = dir.join("proteome.fasta");
        let output = dir.join("out.tsv");
        write_file(&peptides, PEPTIDES_CSV)?;
        write_file(&fasta, PROTEOME_FASTA)?;

        let mut args = base_args(&peptides, &output);
        args.fasta = Some(fasta);
        args.normalize = true;
        run_peptides_to_protein(&args)?;

        let (headers, rows) = read_table(&output)?;
        assert!(headers.contains(&"IbaqNorm".to_owned()));
        assert!(headers.contains(&"IbaqLog".to_owned()));
        assert!(headers.contains(&"IbaqPpb".to_owned()));
        // Condition A total Ibaq = 150 + 500 + 350 = 1000; Condition B = 33.333...
        assert_cell_close(&headers, &rows, "P1", "S1", "IbaqNorm", 0.15)?;
        assert_cell_close(&headers, &rows, "P2", "S1", "IbaqNorm", 0.5)?;
        assert_cell_close(&headers, &rows, "P3", "S1", "IbaqNorm", 0.35)?;
        assert_cell_close(&headers, &rows, "P1", "S2", "IbaqNorm", 1.0)?;
        assert_cell_close(
            &headers,
            &rows,
            "P1",
            "S1",
            "IbaqLog",
            0.15_f64.log10() + 10.0,
        )?;
        assert_cell_close(&headers, &rows, "P1", "S2", "IbaqLog", 10.0)?;
        assert_cell_close(&headers, &rows, "P2", "S1", "IbaqPpb", 50_000_000.0)?;
        assert_cell_close(&headers, &rows, "P1", "S2", "IbaqPpb", 100_000_000.0)?;
        Ok(())
    }

    // Oracle (ProteomicRuler columns; ruler requires --tpa):
    //   ... peptides2protein --method ibaq -f proteome.fasta -p peptides.csv \
    //       --tpa --ruler --organism human --ploidy 2 --cpc 200 -o out.tsv
    // P1/P2/P3 are not human histones, so histone_intensity = max(0, 1) = 1.
    // dna_mass = 2 * 3.22e9 * 617.96 / 6.02214129e23.
    #[test]
    fn peptides2protein_ibaq_ruler_matches_python_oracle() -> TestResult<()> {
        let dir = temp_dir("ibaqruler")?;
        let peptides = dir.join("peptides.csv");
        let fasta = dir.join("proteome.fasta");
        let output = dir.join("out.tsv");
        write_file(&peptides, PEPTIDES_CSV)?;
        write_file(&fasta, PROTEOME_FASTA)?;

        let mut args = base_args(&peptides, &output);
        args.fasta = Some(fasta);
        args.tpa = true;
        args.ruler = true;
        args.organism = "human".to_owned();
        args.ploidy = 2;
        args.cpc = 200.0;
        run_peptides_to_protein(&args)?;

        let (headers, rows) = read_table(&output)?;
        for column in [
            "CopyNumber",
            "Moles[nmol]",
            "Weight[ng]",
            "Concentration[nM]",
        ] {
            assert!(headers.contains(&column.to_owned()), "missing {column}");
        }
        // Captured from the Python oracle (see module-level command above).
        assert_cell_close(&headers, &rows, "P1", "S1", "CopyNumber", 569705926063.9502)?;
        assert_cell_close(&headers, &rows, "P2", "S1", "CopyNumber", 2603104848063.672)?;
        assert_cell_close(
            &headers,
            &rows,
            "P3",
            "S1",
            "CopyNumber",
            1463180476321.2395,
        )?;
        assert_cell_close(
            &headers,
            &rows,
            "P1",
            "S2",
            "CopyNumber",
            126601316903.10007,
        )?;
        assert_cell_close(
            &headers,
            &rows,
            "P1",
            "S1",
            "Moles[nmol]",
            0.0009460188637718765,
        )?;
        assert_cell_close(
            &headers,
            &rows,
            "P1",
            "S1",
            "Weight[ng]",
            2.9737729384957685,
        )?;
        // Concentration is a per-Condition volume term; condition A and B differ.
        assert_cell_close(
            &headers,
            &rows,
            "P1",
            "S1",
            "Concentration[nM]",
            51576.16376717422,
        )?;
        assert_cell_close(
            &headers,
            &rows,
            "P2",
            "S1",
            "Concentration[nM]",
            235662.2176540015,
        )?;
        assert_cell_close(
            &headers,
            &rows,
            "P3",
            "S1",
            "Concentration[nM]",
            132463.49110155678,
        )?;
        assert_cell_close(
            &headers,
            &rows,
            "P1",
            "S2",
            "Concentration[nM]",
            694.6284682447709,
        )?;
        Ok(())
    }

    // The proteomic ruler requires --tpa (mirrors the Python guard); without it
    // the command must fail with an InvalidInput, not silently skip the ruler.
    #[test]
    fn peptides2protein_ruler_requires_tpa() -> TestResult<()> {
        let dir = temp_dir("rulernotpa")?;
        let peptides = dir.join("peptides.csv");
        let fasta = dir.join("proteome.fasta");
        let output = dir.join("out.tsv");
        write_file(&peptides, PEPTIDES_CSV)?;
        write_file(&fasta, PROTEOME_FASTA)?;

        let mut args = base_args(&peptides, &output);
        args.fasta = Some(fasta);
        args.ruler = true;
        assert!(matches!(
            run_peptides_to_protein(&args),
            Err(MokumeError::InvalidInput { .. })
        ));
        Ok(())
    }

    // An unknown organism must fail (mirrors the Python KeyError), even when the
    // ruler is not requested -- the organism is always resolved on the iBAQ path.
    #[test]
    fn peptides2protein_rejects_unknown_organism() -> TestResult<()> {
        let dir = temp_dir("badorg")?;
        let peptides = dir.join("peptides.csv");
        let fasta = dir.join("proteome.fasta");
        let output = dir.join("out.tsv");
        write_file(&peptides, PEPTIDES_CSV)?;
        write_file(&fasta, PROTEOME_FASTA)?;

        let mut args = base_args(&peptides, &output);
        args.fasta = Some(fasta);
        args.organism = "martian".to_owned();
        assert!(matches!(
            run_peptides_to_protein(&args),
            Err(MokumeError::InvalidInput { .. })
        ));
        Ok(())
    }

    #[test]
    fn peptides2protein_runs_lfq_methods() -> TestResult<()> {
        // maxlfq and directlfq both roll the peptide matrix up via the DirectLFQ
        // estimator and write the long-format `Intensity` table.
        let dir = temp_dir("lfq")?;
        let peptides = dir.join("peptides.csv");
        write_file(&peptides, PEPTIDES_CSV)?;

        for method in ["maxlfq", "directlfq"] {
            let output = dir.join(format!("{method}.tsv"));
            let mut args = base_args(&peptides, &output);
            args.method = method.to_owned();
            run_peptides_to_protein(&args)?;

            let (headers, rows) = read_table(&output)?;
            assert_eq!(headers.first().map(String::as_str), Some(PROTEIN_NAME));
            assert_eq!(headers.get(1).map(String::as_str), Some(SAMPLE_ID));
            assert_eq!(headers.get(2).map(String::as_str), Some("Intensity"));
            assert!(!rows.is_empty(), "{method} produced rows");
        }
        Ok(())
    }

    #[test]
    fn peptides2protein_unsupported_enzyme_errors_pointing_to_wheel() -> TestResult<()> {
        // `CNBr` is a real protease pyOpenMS knows but the Rust port has not wired
        // a cleavage rule for. The iBAQ path does not digest it natively and no
        // longer delegates to Python; it fails with a clear `InvalidInput` that
        // names the enzyme and points to the Python wheel (`mokume-rs[ibaq]`).
        // Supported non-Trypsin enzymes (Lys-C, Chymotrypsin, ...) still compute
        // natively in Rust; their digests are oracle-locked in the pipeline crate.
        let dir = temp_dir("enzyme")?;
        let peptides = dir.join("peptides.csv");
        let fasta = dir.join("proteome.fasta");
        let output = dir.join("out.tsv");
        write_file(&peptides, PEPTIDES_CSV)?;
        write_file(&fasta, PROTEOME_FASTA)?;

        let mut args = base_args(&peptides, &output);
        args.fasta = Some(fasta);
        args.enzyme = "CNBr".to_owned();

        match run_peptides_to_protein(&args) {
            Err(MokumeError::InvalidInput { message }) => {
                assert!(
                    message.contains("CNBr") && message.contains("mokume-rs[ibaq]"),
                    "error must name the enzyme and the wheel: {message}"
                );
            }
            other => panic!("expected InvalidInput pointing to the wheel, got {other:?}"),
        }
        Ok(())
    }

    #[test]
    fn peptides2protein_unsupported_enzyme_still_requires_fasta() -> TestResult<()> {
        // The unported-enzyme path keeps the same `--fasta` precondition the native
        // iBAQ path enforces: an unported enzyme without a FASTA still fails with an
        // InvalidInput error.
        let dir = temp_dir("enzyme-nofasta")?;
        let peptides = dir.join("peptides.csv");
        let output = dir.join("out.tsv");
        write_file(&peptides, PEPTIDES_CSV)?;

        let mut args = base_args(&peptides, &output);
        args.enzyme = "CNBr".to_owned();
        let Err(error) = run_peptides_to_protein(&args) else {
            panic!("iBAQ with an unported enzyme but no --fasta must fail");
        };
        assert!(matches!(error, MokumeError::InvalidInput { .. }));
        Ok(())
    }

    #[test]
    fn peptides2protein_verbose_writes_table_no_qc_pdf() -> TestResult<()> {
        // `--verbose` on the iBAQ path computes the protein table in Rust and, since
        // the QC plotting moved to the Python wheel, writes no QC PDF: it only
        // prints a pointer to the wheel command. The run must still succeed and the
        // protein table must contain the `Ibaq` column.
        let dir = temp_dir("verbose")?;
        let peptides = dir.join("peptides.csv");
        let fasta = dir.join("proteome.fasta");
        let output = dir.join("out.tsv");
        let qc = dir.join("qc.pdf");
        write_file(&peptides, PEPTIDES_CSV)?;
        write_file(&fasta, PROTEOME_FASTA)?;

        let mut args = base_args(&peptides, &output);
        args.method = "ibaq".to_owned();
        args.fasta = Some(fasta);
        args.verbose = true;
        args.qc_report = qc.clone();
        run_peptides_to_protein(&args)?;

        // The protein table is written and carries the iBAQ column.
        let (headers, rows) = read_table(&output)?;
        assert!(headers.iter().any(|header| header == "Ibaq"));
        assert!(!rows.is_empty(), "the protein table must contain rows");

        // No QC PDF is produced; that now lives in the Python wheel.
        assert!(
            !qc.exists(),
            "no QC PDF must be written by the Rust compute path"
        );
        Ok(())
    }
}
