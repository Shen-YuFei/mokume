#!/usr/bin/env python3
"""Generate an auditable proteotypic-only iBAQ benchmark baseline.

This benchmark-only implementation deliberately does not call a production
``ibaq`` quantification method. It consumes the same normalized peptide table
used by the piBAQ comparison and uses Mokume's runtime pyOpenMS FASTA digest.
The numerator and denominator both contain only proteotypic peptides. Peptides
mapping to multiple canonical FASTA accessions are recorded in the assignment
audit and discarded. This symmetric baseline isolates piBAQ's shared-peptide
handling; it is not Mokume's production piBAQ implementation or an exact
reproduction of ibaqpy or MaxQuant protein-group inference.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd

from mokume.core.constants import get_accession
from mokume.io.fasta import canonicalize_isoform, digest_fasta_full

ENZYME = "Trypsin"
MIN_AA = 7
MAX_AA = 30
DEFAULT_OUTPUT_DIR = (
    Path(__file__).resolve().parent.parent / "data" / "proteotypic-ibaq"
)
PROTEIN = "ProteinName"
PEPTIDE = "PeptideCanonical"
SAMPLE = "SampleID"
CONDITION = "Condition"
INTENSITY = "NormIntensity"
STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class DatasetInput:
    """One peptide table and the technology label used for its output folder."""

    technology: str
    path: Path


def parse_args() -> argparse.Namespace:
    """Parse explicit local inputs; at least one technology must be supplied."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tmt-peptides", type=Path)
    parser.add_argument("--lfq-peptides", type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Output root (default: ignored benchmark data/proteotypic-ibaq directory)",
    )
    args = parser.parse_args()
    if args.tmt_peptides is None and args.lfq_peptides is None:
        parser.error("at least one of --tmt-peptides or --lfq-peptides is required")
    return args


def resolve_inputs(args: argparse.Namespace) -> tuple[list[DatasetInput], Path, Path]:
    """Resolve paths and report every missing input in one error."""
    datasets = [
        DatasetInput(technology, path.resolve())
        for technology, path in (
            ("tmt", args.tmt_peptides),
            ("lfq", args.lfq_peptides),
        )
        if path is not None
    ]
    fasta = args.fasta.resolve()
    missing = [str(dataset.path) for dataset in datasets if not dataset.path.is_file()]
    if not fasta.is_file():
        missing.append(str(fasta))
    if missing:
        raise FileNotFoundError(f"Missing benchmark input(s): {', '.join(missing)}")
    return datasets, fasta, args.output_dir.resolve()


def validate_intensities(peptides: pd.DataFrame, path: Path) -> pd.Series:
    """Validate non-negative, finite intensities without silently dropping rows."""
    raw_intensities = peptides[INTENSITY]
    intensities = pd.to_numeric(raw_intensities, errors="coerce")
    invalid_intensities = (
        intensities.isna() | ~np.isfinite(intensities) | intensities.lt(0)
    )
    if invalid_intensities.any():
        examples = raw_intensities.loc[invalid_intensities].astype(str).unique()[:3]
        raise ValueError(
            f"{path} contains invalid {INTENSITY} values, for example: "
            f"{examples.tolist()}"
        )
    if not intensities.gt(0).any():
        raise ValueError(f"{path} contains no positive peptide measurements")
    return intensities


def validate_conditions(peptides: pd.DataFrame, path: Path) -> None:
    """Require an optional condition column to map one-to-one from samples."""
    if CONDITION not in peptides:
        return
    missing_condition = peptides[CONDITION].isna() | peptides[CONDITION].astype(
        str
    ).str.strip().eq("")
    if missing_condition.any():
        raise ValueError(f"{path} contains missing {CONDITION} values")
    conditions_per_sample = peptides.groupby(SAMPLE, observed=True)[CONDITION].nunique()
    if conditions_per_sample.gt(1).any():
        samples = sorted(conditions_per_sample[conditions_per_sample.gt(1)].index)
        raise ValueError(
            f"{path} maps a sample to multiple conditions, for example: {samples[:3]}"
        )


def validate_sequences_and_keys(peptides: pd.DataFrame, path: Path) -> None:
    """Reject modified sequences and duplicate peptide-level measurement keys."""
    invalid_peptides = ~peptides[PEPTIDE].map(
        lambda peptide: bool(peptide) and set(peptide) <= STANDARD_AMINO_ACIDS
    )
    if invalid_peptides.any():
        examples = sorted(peptides.loc[invalid_peptides, PEPTIDE].unique())[:3]
        raise ValueError(
            f"{path} contains modified or non-standard peptide sequences, "
            f"for example: {examples}"
        )
    duplicate_keys = [PROTEIN, PEPTIDE, SAMPLE]
    if peptides.duplicated(duplicate_keys).any():
        raise ValueError(
            f"{path} contains duplicate protein/peptide/sample measurements; "
            "provide a peptide-level Mokume intermediate"
        )


def read_peptides(path: Path) -> pd.DataFrame:
    """Load and validate one Mokume peptide-level intermediate."""
    if path.suffix.lower() == ".parquet":
        peptides = pd.read_parquet(path)
    else:
        separator = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
        peptides = pd.read_csv(path, sep=separator)
    if PEPTIDE not in peptides and "PeptideSequence" in peptides:
        peptides = peptides.rename(columns={"PeptideSequence": PEPTIDE})
    required = {PROTEIN, PEPTIDE, SAMPLE, INTENSITY}
    missing = required.difference(peptides.columns)
    if missing:
        raise ValueError(f"{path} is missing peptide columns: {sorted(missing)}")
    peptides = peptides.copy()
    if peptides[list(required)].isna().any(axis=None):
        raise ValueError(f"{path} contains missing required peptide values")
    for column in (PROTEIN, PEPTIDE, SAMPLE):
        peptides[column] = peptides[column].astype(str).str.strip()
    if peptides[[PROTEIN, PEPTIDE, SAMPLE]].eq("").any(axis=None):
        raise ValueError(f"{path} contains blank required peptide values")
    peptides[INTENSITY] = validate_intensities(peptides, path)
    peptides[PEPTIDE] = peptides[PEPTIDE].str.upper()
    validate_conditions(peptides, path)
    validate_sequences_and_keys(peptides, path)
    return peptides


def group_accessions(protein_group: str) -> tuple[str, ...]:
    """Return all canonical accessions in a reported group without choosing a lead."""
    accessions = {
        canonicalize_isoform(get_accession(token.strip()))
        for token in str(protein_group).split(";")
        if token.strip()
    }
    return tuple(sorted(accession for accession in accessions if accession))


def assignment_row(
    protein_group: str,
    peptide: str,
    peptide_to_accessions: Mapping[str, set[str]],
) -> dict[str, object]:
    """Classify one reported group/peptide pair using its complete FASTA mapping."""
    reported = group_accessions(protein_group)
    fasta_accessions = tuple(sorted(peptide_to_accessions.get(str(peptide), set())))
    assigned: str | None = None
    if not fasta_accessions:
        status = "peptide_not_in_fasta"
    elif len(fasta_accessions) > 1:
        status = "shared_in_fasta"
    elif fasta_accessions[0] not in reported:
        status = "reported_group_mismatch"
    else:
        status = "assigned_proteotypic"
        assigned = fasta_accessions[0]
    return {
        "ReportedProteinGroup": protein_group,
        PEPTIDE: peptide,
        "ReportedAccessions": ";".join(reported),
        "FastaAccessions": ";".join(fasta_accessions),
        "AssignedProtein": assigned,
        "AssignmentStatus": status,
    }


def build_assignment_audit(
    peptides: pd.DataFrame,
    peptide_to_accessions: Mapping[str, set[str]],
) -> pd.DataFrame:
    """Build one deterministic audit row per reported group/peptide pair."""
    pairs = peptides[[PROTEIN, PEPTIDE]].drop_duplicates()
    records = [
        assignment_row(protein_group, peptide, peptide_to_accessions)
        for protein_group, peptide in pairs.itertuples(index=False, name=None)
    ]
    return pd.DataFrame.from_records(records)


def proteotypic_theoretical_counts(
    peptide_to_accessions: Mapping[str, set[str]],
) -> dict[str, int]:
    """Count theoretical peptides owned by exactly one canonical FASTA accession."""
    counts: dict[str, int] = {}
    for owners in peptide_to_accessions.values():
        if len(owners) == 1:
            owner = next(iter(owners))
            counts[owner] = counts.get(owner, 0) + 1
    return counts


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of one provenance input or output."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def installed_version(distribution: str) -> str:
    """Return an installed version, including a clear source-tree fallback."""
    try:
        return version(distribution)
    except PackageNotFoundError:
        return "uninstalled-source-tree"


def provenance_record(path: Path) -> dict[str, str]:
    """Describe one provenance file with its resolved path and SHA-256."""
    return {"path": str(path.resolve()), "sha256": sha256_file(path)}


def write_provenance(
    dataset: DatasetInput,
    fasta: Path,
    output_dir: Path,
    counts: dict[str, object],
) -> Path:
    """Write the complete runtime and file provenance for one baseline."""
    result_path = output_dir / "proteotypic_ibaq.parquet"
    audit_path = output_dir / "proteotypic_ibaq_assignment_audit.parquet"
    provenance = {
        "technology": dataset.technology,
        "definition": "proteotypic-only numerator / proteotypic-only denominator",
        "parameters": {
            "enzyme": ENZYME,
            "min_aa": MIN_AA,
            "max_aa": MAX_AA,
            "missed_cleavages": 0,
        },
        "software": {
            "mokume": installed_version("mokume"),
            "pyopenms": installed_version("pyopenms"),
        },
        "inputs": {
            "script": provenance_record(Path(__file__)),
            "peptides": provenance_record(dataset.path),
            "fasta": provenance_record(fasta),
        },
        "outputs": {
            "baseline": provenance_record(result_path),
            "assignment_audit": provenance_record(audit_path),
        },
        "counts": counts,
    }
    path = output_dir / "provenance.json"
    path.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def compute_proteotypic_ibaq(
    peptides: pd.DataFrame,
    audit: pd.DataFrame,
    theoretical_counts: Mapping[str, int],
) -> pd.DataFrame:
    """Compute iBAQ with a symmetric proteotypic-only denominator."""
    assignments = audit[["ReportedProteinGroup", PEPTIDE, "AssignedProtein"]]
    assigned = peptides.merge(
        assignments,
        left_on=[PROTEIN, PEPTIDE],
        right_on=["ReportedProteinGroup", PEPTIDE],
        how="left",
        validate="many_to_one",
    ).dropna(subset=["AssignedProtein"])
    assigned = assigned[assigned[INTENSITY] > 0]
    group_columns = ["AssignedProtein", SAMPLE]
    if CONDITION in assigned.columns:
        group_columns.append(CONDITION)
    assigned_keys = [*group_columns, PEPTIDE]
    if assigned.duplicated(assigned_keys).any():
        raise ValueError(
            "multiple reported groups assign the same peptide/sample measurement "
            "to one FASTA accession"
        )
    result = assigned.groupby(group_columns, observed=True, as_index=False).agg(
        ProteotypicIntensity=(INTENSITY, "sum"),
        ObservedProteotypicPeptides=(PEPTIDE, "nunique"),
    )
    result = result.rename(columns={"AssignedProtein": PROTEIN})
    result["TheoreticalProteotypicPeptides"] = result[PROTEIN].map(theoretical_counts)
    result = result[result["TheoreticalProteotypicPeptides"].fillna(0) > 0].copy()
    result["Intensity"] = (
        result["ProteotypicIntensity"] / result["TheoreticalProteotypicPeptides"]
    )
    return result.sort_values([PROTEIN, SAMPLE]).reset_index(drop=True)


def write_dataset(
    dataset: DatasetInput,
    fasta: Path,
    output_root: Path,
    peptide_to_accessions: Mapping[str, set[str]],
    theoretical_counts: Mapping[str, int],
) -> None:
    """Compute and write one technology's baseline plus assignment audit."""
    peptides = read_peptides(dataset.path)
    audit = build_assignment_audit(peptides, peptide_to_accessions)
    result = compute_proteotypic_ibaq(peptides, audit, theoretical_counts)
    if result.empty:
        raise ValueError(f"{dataset.path} produced no proteotypic-only iBAQ values")
    output_dir = output_root / dataset.technology
    output_dir.mkdir(parents=True, exist_ok=True)
    result_path = output_dir / "proteotypic_ibaq.parquet"
    audit_path = output_dir / "proteotypic_ibaq_assignment_audit.parquet"
    result.to_parquet(result_path, index=False)
    audit.to_parquet(audit_path, index=False)
    status_counts = audit["AssignmentStatus"].value_counts().to_dict()
    counts = {
        "input_rows": len(peptides),
        "positive_input_rows": int(peptides[INTENSITY].gt(0).sum()),
        "zero_input_rows": int(peptides[INTENSITY].eq(0).sum()),
        "result_rows": len(result),
        "result_proteins": int(result[PROTEIN].nunique()),
        "assignment_status": status_counts,
    }
    provenance_path = write_provenance(
        dataset,
        fasta,
        output_dir,
        counts,
    )
    print(f"{dataset.technology.upper()}: {len(result):,} protein-sample values")
    print(f"  assignment status: {status_counts}")
    print(f"  baseline: {result_path}")
    print(f"  audit: {audit_path}")
    print(f"  provenance: {provenance_path}")


def main() -> None:
    """Generate requested proteotypic-only iBAQ baselines."""
    datasets, fasta, output_root = resolve_inputs(parse_args())
    _, peptide_to_accessions, _ = digest_fasta_full(
        str(fasta), ENZYME, MIN_AA, MAX_AA, canonicalize_isoforms=True
    )
    theoretical_counts = proteotypic_theoretical_counts(peptide_to_accessions)
    for dataset in datasets:
        write_dataset(
            dataset,
            fasta,
            output_root,
            peptide_to_accessions,
            theoretical_counts,
        )


if __name__ == "__main__":
    main()
