#!/usr/bin/env python3
"""
Generate iBAQ quantification for PXD007683 benchmark.

This script:
1. Downloads the FASTA file if not present
2. Loads peptide-level data
3. Computes iBAQ values using theoretical tryptic peptide counts
4. Saves results to parquet files
"""

import shutil
import urllib.request
import urllib.parse
import warnings

import pandas as pd
import numpy as np

from config import (
    DATA_DIR, LOCAL_DATA_DIR, LOCAL_QUANTIFIED_DIR,
    FASTA_URL, FASTA_PATH,
)
from benchmark_utils import (
    parse_fasta, count_theoretical_peptides, get_accession,
)

warnings.filterwarnings("ignore")

_opener = urllib.request.build_opener()


def download_fasta():
    """Download FASTA file if not present."""
    if FASTA_PATH.exists():
        print(f"FASTA file already exists: {FASTA_PATH}")
        return True

    print(f"Downloading FASTA from: {FASTA_URL}")
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    try:
        if urllib.parse.urlparse(FASTA_URL).scheme not in ("http", "https"):
            raise ValueError(f"URL scheme not allowed: {FASTA_URL}")
        with _opener.open(FASTA_URL) as response, open(FASTA_PATH, "wb") as out_file:
            shutil.copyfileobj(response, out_file)
        print(f"Downloaded to: {FASTA_PATH}")
        return True
    except Exception as e:
        print(f"Failed to download FASTA: {e}")
        return False


def compute_ibaq(technology: str) -> pd.DataFrame:
    """
    Compute iBAQ values for a technology.

    iBAQ = Sum of peptide intensities / Number of theoretical tryptic peptides

    Parameters
    ----------
    technology : str
        Either "tmt" or "lfq"

    Returns
    -------
    pd.DataFrame
        DataFrame with iBAQ values per protein per sample
    """
    PROTEIN_NAME = "ProteinName"
    SAMPLE_ID = "SampleID"
    NORM_INTENSITY = "NormIntensity"

    # Load peptide-level data
    peptide_path = LOCAL_DATA_DIR / f"{technology}_peptides.parquet"
    if not peptide_path.exists():
        raise FileNotFoundError(f"Peptide data not found: {peptide_path}")

    print(f"  Loading peptide data from: {peptide_path}")
    peptides_df = pd.read_parquet(peptide_path)
    print(f"  Loaded {len(peptides_df)} peptide entries")

    # Get unique protein accessions (extract from protein group format)
    all_proteins = peptides_df[PROTEIN_NAME].unique()

    # Extract individual accessions from protein groups (e.g., "sp|P12345|NAME;sp|P67890|NAME")
    protein_accessions = set()
    for prot in all_proteins:
        for p in str(prot).split(";"):
            acc = get_accession(p.strip())
            if acc:
                protein_accessions.add(acc)

    print(f"  Found {len(protein_accessions)} unique protein accessions")

    # Parse FASTA and compute theoretical peptide counts
    print("  Computing theoretical peptides from FASTA...")
    fasta_proteins = parse_fasta(str(FASTA_PATH))

    peptide_counts = {}
    found_proteins = set()
    for acc in protein_accessions:
        if acc in fasta_proteins:
            seq = fasta_proteins[acc]
            n_pep = count_theoretical_peptides(seq, enzyme="Trypsin", min_aa=7, max_aa=30)
            peptide_counts[acc] = n_pep
            found_proteins.add(acc)

    print(f"  Found {len(found_proteins)} proteins in FASTA ({len(peptide_counts)} with peptide counts)")

    # Create mapping from full protein name to accession
    protein_to_accession = {}
    for prot in all_proteins:
        for p in str(prot).split(";"):
            acc = get_accession(p.strip())
            if acc and acc in found_proteins:
                protein_to_accession[prot] = acc
                break

    print(f"  Mapped {len(protein_to_accession)} protein groups to FASTA entries")

    # Compute iBAQ for each sample
    samples = peptides_df[SAMPLE_ID].unique()
    print(f"  Computing iBAQ for {len(samples)} samples...")

    results = []
    for sample in samples:
        sample_df = peptides_df[peptides_df[SAMPLE_ID] == sample].copy()

        # Group by protein and sum intensities
        protein_intensities = sample_df.groupby(PROTEIN_NAME)[NORM_INTENSITY].sum()

        for protein, intensity in protein_intensities.items():
            accession = protein_to_accession.get(protein)

            if accession and accession in peptide_counts:
                n_peptides = peptide_counts[accession]
                if n_peptides > 0:
                    ibaq_value = intensity / n_peptides
                else:
                    ibaq_value = np.nan
            else:
                ibaq_value = np.nan
                n_peptides = None

            results.append({
                "ProteinName": protein,
                "SampleID": sample,
                "Intensity": ibaq_value,
                "TheoreticalPeptides": n_peptides,
            })

    result_df = pd.DataFrame(results)

    # Remove proteins with no FASTA match
    n_before = len(result_df)
    result_df = result_df.dropna(subset=["Intensity"])
    n_after = len(result_df)
    print(f"  Removed {n_before - n_after} entries without FASTA match ({n_after} remaining)")

    return result_df


def main():
    """Generate iBAQ quantification for both technologies."""
    print("=" * 60)
    print("PXD007683: Generate iBAQ Quantification")
    print("=" * 60)

    # Download FASTA
    print("\n1. Checking FASTA file...")
    if not download_fasta():
        print("ERROR: Could not download FASTA file")
        return

    # Process each technology
    for technology in ["tmt", "lfq"]:
        print(f"\n2. Processing {technology.upper()}...")

        try:
            ibaq_df = compute_ibaq(technology)
            print(f"  Computed iBAQ for {ibaq_df['ProteinName'].nunique()} proteins")

            # Save to parquet
            output_dir = LOCAL_QUANTIFIED_DIR / technology
            output_dir.mkdir(parents=True, exist_ok=True)
            output_path = output_dir / "ibaq.parquet"

            ibaq_df.to_parquet(output_path, index=False)
            print(f"  Saved to: {output_path}")

            # Summary stats
            pivot = ibaq_df.pivot(index="ProteinName", columns="SampleID", values="Intensity")
            print(f"  Proteins: {len(pivot)}")
            print(f"  Samples: {len(pivot.columns)}")
            print(f"  Missing rate: {pivot.isnull().mean().mean()*100:.1f}%")

        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()

    print("\n" + "=" * 60)
    print("Done!")


if __name__ == "__main__":
    main()
