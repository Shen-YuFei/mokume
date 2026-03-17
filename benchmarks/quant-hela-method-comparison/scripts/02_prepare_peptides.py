#!/usr/bin/env python3
"""
Phase 2: Data Preparation

Convert downloaded raw data to standardized peptide format for quantification.
Handles parquet, TSV, and various tool output formats.
"""

import sys
from typing import Optional

import pandas as pd
import numpy as np

from config import (
    ALL_DATASETS,
    HELA_DATASETS,
    RAW_DATA_DIR,
    PROCESSED_DIR,
    QUANTMS_PARQUET_MAPPING,
    DatasetInfo,
)


def detect_data_format(df: pd.DataFrame, file_format: str = None) -> str:
    """
    Detect the data format based on column names or explicit format.

    Returns one of: 'quantms', 'diann', 'maxquant', 'msstats', 'unknown'
    """
    # Use explicit format if provided
    if file_format == "msstats_csv":
        return "msstats"
    if file_format == "parquet":
        # Check column names for quantms parquet format
        columns = set(df.columns.str.lower())
        if "pg" in columns or "peptidoform" in columns:
            return "quantms"

    columns = set(df.columns.str.lower())

    # quantms/ibaqpy format (parquet files)
    if "pg" in columns and "peptidoform" in columns:
        return "quantms"

    # MSstats format (CSV from proteomicslfq)
    if "proteinname" in columns and "peptidesequence" in columns:
        return "msstats"

    # DIA-NN format
    if "protein.group" in [c.lower() for c in df.columns]:
        return "diann"

    # MaxQuant peptides.txt format
    if "leading razor protein" in [c.lower() for c in df.columns]:
        return "maxquant"

    return "unknown"


def handle_ibaqpy_format(df: pd.DataFrame) -> pd.DataFrame:
    """
    Handle ibaqpy-research format with nested intensities.

    The intensities column contains arrays of dictionaries:
    [{'sample_accession': 'X', 'channel': 'LFQ', 'intensity': 123.0}, ...]

    This function explodes the data to one row per sample.
    """
    # Check if this is ibaqpy format
    if "intensities" not in df.columns:
        return df

    # Check if intensities contains nested dicts
    sample_val = df["intensities"].iloc[0]
    if not isinstance(sample_val, (list, np.ndarray)):
        return df
    if len(sample_val) == 0:
        return df
    if not isinstance(sample_val[0], dict):
        return df

    print("  Detected ibaqpy nested format - expanding intensities...")

    # Handle pg_accessions first (take first protein)
    if "pg_accessions" in df.columns:
        df["ProteinName"] = df["pg_accessions"].apply(
            lambda x: x[0] if isinstance(x, (list, np.ndarray)) and len(x) > 0 else str(x) if pd.notna(x) else ""
        )

    # Rename sequence to PeptideSequence
    if "sequence" in df.columns:
        df["PeptideSequence"] = df["sequence"]

    # Explode intensities
    rows = []
    for idx, row in df.iterrows():
        intensities = row["intensities"]
        if isinstance(intensities, (list, np.ndarray)):
            for entry in intensities:
                if isinstance(entry, dict) and "intensity" in entry:
                    new_row = {
                        "ProteinName": row.get("ProteinName", ""),
                        "PeptideSequence": row.get("PeptideSequence", row.get("sequence", "")),
                        "SampleID": entry.get("sample_accession", row.get("reference_file_name", "")),
                        "NormIntensity": entry.get("intensity", 0.0),
                    }
                    rows.append(new_row)

    result = pd.DataFrame(rows)
    print(f"  Expanded to {len(result)} rows")
    return result


def handle_array_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Handle simple array-type columns (non-nested).

    For columns like pg_accessions that are arrays of strings.
    """
    # Handle pg_accessions (take first protein)
    if "pg_accessions" in df.columns and "ProteinName" not in df.columns:
        df["pg_accessions"] = df["pg_accessions"].apply(
            lambda x: x[0] if isinstance(x, (list, np.ndarray)) and len(x) > 0 else str(x) if pd.notna(x) else ""
        )

    return df


def normalize_columns(df: pd.DataFrame, data_format: str) -> pd.DataFrame:
    """
    Normalize column names to standard mokume format.

    Handles multiple possible source column names for each target.
    """
    # First check for ibaqpy nested format and expand
    df = handle_ibaqpy_format(df)

    # If already expanded, return early
    required = ["ProteinName", "PeptideSequence", "SampleID", "NormIntensity"]
    if all(col in df.columns for col in required):
        return df

    # Handle simple array-type columns
    df = handle_array_columns(df)

    # Handle MSstats format - create SampleID from BioReplicate
    if data_format == "msstats":
        if "SampleID" not in df.columns:
            if "BioReplicate" in df.columns:
                df["SampleID"] = df["BioReplicate"].astype(str)
            elif "Run" in df.columns:
                df["SampleID"] = df["Run"].astype(str)
        if "NormIntensity" not in df.columns and "Intensity" in df.columns:
            df["NormIntensity"] = df["Intensity"]

    # Target columns we need
    targets = ["ProteinName", "PeptideSequence", "SampleID", "NormIntensity"]

    # Build a reverse mapping: target -> list of possible source names
    # Use QUANTMS_PARQUET_MAPPING for all formats (ibaqpy-research uses this)
    source_mapping = QUANTMS_PARQUET_MAPPING

    # Create actual column rename mapping
    rename_map = {}
    existing_cols = set(df.columns)

    for source, target in source_mapping.items():
        if source in existing_cols and target not in existing_cols:
            # Only rename if source exists and target doesn't
            if target not in rename_map.values():
                rename_map[source] = target

    # Fallback: try to auto-detect by partial match
    existing_after = existing_cols - set(rename_map.keys())
    assigned_targets = set(rename_map.values())

    for col in list(existing_after):
        col_lower = col.lower()
        if "ProteinName" not in assigned_targets:
            if "protein" in col_lower and "accession" in col_lower:
                rename_map[col] = "ProteinName"
                assigned_targets.add("ProteinName")
                continue
        if "PeptideSequence" not in assigned_targets:
            if col == "sequence" or (("peptide" in col_lower or "sequence" in col_lower) and "protein" not in col_lower):
                rename_map[col] = "PeptideSequence"
                assigned_targets.add("PeptideSequence")
                continue
        if "SampleID" not in assigned_targets:
            if "sample" in col_lower and "id" in col_lower:
                rename_map[col] = "SampleID"
                assigned_targets.add("SampleID")
                continue
        if "NormIntensity" not in assigned_targets:
            if col_lower == "intensity" or "normintensity" in col_lower:
                rename_map[col] = "NormIntensity"
                assigned_targets.add("NormIntensity")
                continue

    # Apply mapping
    df = df.rename(columns=rename_map)
    return df


def extract_protein_accession(protein_val) -> str:
    """
    Extract UniProt accession from protein identifier.

    Handles formats like:
    - "P12345"
    - "sp|P12345|PROTEIN_NAME"
    - "P12345;P67890" (takes first)
    - ["P12345", "P67890"] (array, takes first)
    - np.array(["P12345", "P67890"]) (numpy array, takes first)
    """
    # Handle None/NaN
    if protein_val is None:
        return ""

    # Handle numpy arrays or lists
    if isinstance(protein_val, (list, np.ndarray)):
        if len(protein_val) == 0:
            return ""
        protein_str = str(protein_val[0])
    elif pd.isna(protein_val):
        return ""
    else:
        protein_str = str(protein_val)

    # Handle semicolon-separated groups (take first)
    if ";" in protein_str:
        protein_str = protein_str.split(";")[0]

    # Handle UniProt format (sp|P12345|NAME)
    if "|" in protein_str:
        parts = protein_str.split("|")
        if len(parts) >= 2:
            return parts[1]

    # Already an accession
    return protein_str.strip()


def prepare_peptide_data(
    df: pd.DataFrame,
    dataset_id: str,
    file_format: str = None
) -> pd.DataFrame:
    """
    Prepare peptide-level data for quantification.

    Parameters
    ----------
    df : pd.DataFrame
        Raw data with intensity values
    dataset_id : str
        Dataset identifier for logging
    file_format : str, optional
        Explicit file format hint (parquet, msstats_csv, etc.)

    Returns
    -------
    pd.DataFrame
        Standardized peptide data with columns:
        - ProteinName: UniProt accession
        - PeptideSequence: Peptide identifier
        - SampleID: Sample identifier
        - NormIntensity: Intensity value
    """
    print(f"  Input shape: {df.shape}")

    # Detect format
    data_format = detect_data_format(df, file_format)
    print(f"  Detected format: {data_format}")

    # Normalize columns
    df = normalize_columns(df, data_format)

    # Check required columns
    required = ["ProteinName", "PeptideSequence", "SampleID", "NormIntensity"]
    missing = [c for c in required if c not in df.columns]

    if missing:
        print(f"  WARNING: Missing columns: {missing}")
        print(f"  Available columns: {list(df.columns)}")

        # Try to find alternatives
        for col in df.columns:
            print(f"    - {col}")

        return pd.DataFrame()

    # Extract accessions
    df["ProteinName"] = df["ProteinName"].apply(extract_protein_accession)

    # Remove rows with missing values
    original_len = len(df)
    df = df[df["NormIntensity"].notna()]
    df = df[df["NormIntensity"] > 0]
    df = df[df["ProteinName"].notna()]
    df = df[df["ProteinName"] != ""]
    df = df[df["PeptideSequence"].notna()]

    # Remove decoys and contaminants
    df = df[~df["ProteinName"].str.startswith("DECOY_", na=False)]
    df = df[~df["ProteinName"].str.startswith("CONTAM_", na=False)]
    df = df[~df["ProteinName"].str.startswith("CON__", na=False)]
    df = df[~df["ProteinName"].str.startswith("REV__", na=False)]

    print(f"  Filtered: {original_len} -> {len(df)} rows")

    # Select and order columns
    result = df[["ProteinName", "PeptideSequence", "SampleID", "NormIntensity"]].copy()

    # Add optional columns if available
    optional_cols = ["Condition", "BioReplicate", "TechReplicate", "Fraction", "Run"]
    for col in optional_cols:
        if col in df.columns:
            result[col] = df[col]

    # Add default Condition if missing (required by iBAQ)
    if "Condition" not in result.columns:
        result["Condition"] = "default"

    # Summary statistics
    n_proteins = result["ProteinName"].nunique()
    n_peptides = result["PeptideSequence"].nunique()
    n_samples = result["SampleID"].nunique()

    print(f"  Proteins: {n_proteins}, Peptides: {n_peptides}, Samples: {n_samples}")

    return result


def load_raw_data(dataset: DatasetInfo) -> Optional[pd.DataFrame]:
    """
    Load raw data for a dataset.

    Tries the expected file path based on file_format, then alternatives.
    """
    feature_path = dataset.local_feature_path

    # Try the expected path first
    if feature_path.exists():
        try:
            if feature_path.suffix == ".parquet":
                print(f"  Loading parquet: {feature_path.name}")
                return pd.read_parquet(feature_path)
            elif feature_path.suffix == ".csv":
                print(f"  Loading CSV: {feature_path.name}")
                return pd.read_csv(feature_path)
            elif feature_path.suffix == ".tsv":
                print(f"  Loading TSV: {feature_path.name}")
                return pd.read_csv(feature_path, sep="\t")
        except Exception as e:
            print(f"  ERROR loading file: {e}")
            return None

    # Try alternative extensions
    base_name = feature_path.stem.replace("_feature", "")
    for ext in [".parquet", ".csv", ".tsv"]:
        alt_path = RAW_DATA_DIR / f"{dataset.project_id}_feature{ext}"
        if alt_path.exists():
            if ext == ".parquet":
                print(f"  Loading parquet: {alt_path.name}")
                return pd.read_parquet(alt_path)
            elif ext == ".csv":
                print(f"  Loading CSV: {alt_path.name}")
                return pd.read_csv(alt_path)
            else:
                print(f"  Loading TSV: {alt_path.name}")
                return pd.read_csv(alt_path, sep="\t")

    # Check for pattern matches
    for pattern in [f"{dataset.project_id}*feature*", f"*{dataset.project_id}*"]:
        matches = list(RAW_DATA_DIR.glob(pattern))
        for match in matches:
            if match.suffix in [".parquet", ".tsv", ".csv"]:
                print(f"  Found alternative: {match.name}")
                if match.suffix == ".parquet":
                    return pd.read_parquet(match)
                elif match.suffix == ".tsv":
                    return pd.read_csv(match, sep="\t")
                else:
                    return pd.read_csv(match)

    print(f"  ERROR: No data file found for {dataset.project_id}")
    print(f"  Expected: {feature_path}")
    return None


def process_dataset(dataset: DatasetInfo, force: bool = False) -> bool:
    """
    Process a single dataset.

    Returns True if successful.
    """
    print(f"\n{dataset.project_id}: {dataset.name}")

    # Check if already processed
    output_path = dataset.peptide_path
    if output_path.exists() and not force:
        print(f"  Already processed: {output_path.name}")
        return True

    # Load raw data
    df = load_raw_data(dataset)
    if df is None:
        return False

    # Prepare peptide data
    peptide_df = prepare_peptide_data(df, dataset.project_id, dataset.file_format)
    if peptide_df.empty:
        return False

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    peptide_df.to_parquet(output_path, index=False)
    print(f"  Saved: {output_path.name}")

    return True


def process_all_datasets(
    datasets: dict = None,
    force: bool = False,
) -> dict:
    """
    Process all datasets.

    Returns dict with processing status.
    """
    if datasets is None:
        datasets = ALL_DATASETS

    print("=" * 70)
    print("Peptide Data Preparation")
    print("=" * 70)
    print(f"\nDatasets to process: {len(datasets)}")
    print(f"Output directory: {PROCESSED_DIR}")

    results = {}

    for dataset_id, dataset in datasets.items():
        success = process_dataset(dataset, force=force)
        results[dataset_id] = {"success": success}

    # Summary
    print("\n" + "=" * 70)
    print("Processing Summary")
    print("=" * 70)

    success_count = sum(1 for r in results.values() if r["success"])
    print(f"\nSuccessfully processed: {success_count}/{len(datasets)} datasets")

    for dataset_id, result in results.items():
        status = "OK" if result["success"] else "FAILED"
        print(f"  {dataset_id}: {status}")

    return results


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Prepare peptide data for quantification"
    )
    parser.add_argument(
        "--dataset",
        type=str,
        help="Process specific dataset (project ID)"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Reprocess even if output exists"
    )
    parser.add_argument(
        "--hela-only",
        action="store_true",
        help="Only process HeLa datasets"
    )

    args = parser.parse_args()

    # Select datasets
    if args.dataset:
        if args.dataset in ALL_DATASETS:
            datasets = {args.dataset: ALL_DATASETS[args.dataset]}
        else:
            print(f"Unknown dataset: {args.dataset}")
            print(f"Available: {', '.join(ALL_DATASETS.keys())}")
            sys.exit(1)
    elif args.hela_only:
        datasets = HELA_DATASETS
    else:
        datasets = ALL_DATASETS

    # Run
    results = process_all_datasets(datasets=datasets, force=args.force)

    # Exit code
    success_count = sum(1 for r in results.values() if r["success"])
    sys.exit(0 if success_count > 0 else 1)


if __name__ == "__main__":
    main()
