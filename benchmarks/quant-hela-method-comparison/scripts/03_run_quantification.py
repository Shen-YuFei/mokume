#!/usr/bin/env python3
"""
Phase 3: Run Quantification Methods

Run all protein quantification methods on prepared peptide data.
Methods: iBAQ, rIBAQ, DirectLFQ, Top3, TopN, Sum
"""

import sys
import time
from pathlib import Path
from typing import Optional, Dict, List

import pandas as pd
import numpy as np

# Add parent directory to path for mokume imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from config import (
    ALL_DATASETS,
    HELA_DATASETS,
    PROCESSED_DIR,
    PROTEIN_QUANT_DIR,
    FASTA_FILE,
    QUANTIFICATION_METHODS,
    TOPN_VALUES,
    IBAQ_PARAMS,
    DatasetInfo,
)


def run_ibaq(
    peptide_df: pd.DataFrame,
    fasta_path: Path,
    dataset_id: str,
) -> Optional[pd.DataFrame]:
    """
    Run iBAQ quantification.

    Requires FASTA file for theoretical peptide calculation.
    Uses the file-based peptides_to_protein function.
    """
    import tempfile
    import os
    import gzip
    import shutil

    if not fasta_path.exists():
        print(f"  WARNING: FASTA file not found: {fasta_path}")
        return None

    try:
        from mokume.quantification import peptides_to_protein

        # Create temporary files for input/output
        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = os.path.join(tmpdir, "peptides.parquet")
            output_file = os.path.join(tmpdir, "ibaq_output.tsv")

            # Handle gzipped FASTA files
            fasta_file = str(fasta_path)
            if fasta_path.suffix == ".gz":
                decompressed_fasta = os.path.join(tmpdir, "reference.fasta")
                with gzip.open(fasta_path, 'rb') as f_in:
                    with open(decompressed_fasta, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                fasta_file = decompressed_fasta

            # Save peptide data
            peptide_df.to_parquet(input_file, index=False)

            # Run iBAQ (file-based function)
            peptides_to_protein(
                fasta=fasta_file,
                peptides=input_file,
                enzyme=IBAQ_PARAMS["enzyme"],
                normalize=True,
                min_aa=IBAQ_PARAMS["min_aa"],
                max_aa=IBAQ_PARAMS["max_aa"],
                tpa=False,
                ruler=False,
                ploidy=2,
                cpc=200.0,
                organism="human",
                output=output_file,
                verbose=False,
                qc_report="",
            )

            # Read results
            if os.path.exists(output_file):
                result = pd.read_csv(output_file, sep="\t")
                return result
            else:
                print(f"  ERROR: iBAQ output file not generated")
                return None

    except Exception as e:
        print(f"  ERROR running iBAQ: {e}")
        import traceback
        traceback.print_exc()
        return None


def run_maxlfq(
    peptide_df: pd.DataFrame,
    min_peptides: int = 2,
) -> Optional[pd.DataFrame]:
    """
    Run MaxLFQ quantification.

    Uses DirectLFQ if available, otherwise built-in implementation.
    """
    try:
        from mokume.quantification import MaxLFQQuantification

        maxlfq = MaxLFQQuantification(min_peptides=min_peptides, threads=-1)
        result = maxlfq.quantify(
            peptide_df,
            protein_column="ProteinName",
            peptide_column="PeptideSequence",
            intensity_column="NormIntensity",
            sample_column="SampleID",
        )
        return result

    except Exception as e:
        print(f"  ERROR running MaxLFQ: {e}")
        import traceback
        traceback.print_exc()
        return None


def run_top3(peptide_df: pd.DataFrame) -> Optional[pd.DataFrame]:
    """Run Top3 quantification."""
    try:
        from mokume.quantification import Top3Quantification

        top3 = Top3Quantification()
        result = top3.quantify(
            peptide_df,
            protein_column="ProteinName",
            peptide_column="PeptideSequence",
            intensity_column="NormIntensity",
            sample_column="SampleID",
        )
        return result

    except Exception as e:
        print(f"  ERROR running Top3: {e}")
        return None


def run_topn(peptide_df: pd.DataFrame, n: int = 10) -> Optional[pd.DataFrame]:
    """Run TopN quantification."""
    try:
        from mokume.quantification import TopNQuantification

        topn = TopNQuantification(n=n)
        result = topn.quantify(
            peptide_df,
            protein_column="ProteinName",
            peptide_column="PeptideSequence",
            intensity_column="NormIntensity",
            sample_column="SampleID",
        )
        return result

    except Exception as e:
        print(f"  ERROR running Top{n}: {e}")
        return None


def run_sum(peptide_df: pd.DataFrame) -> Optional[pd.DataFrame]:
    """Run Sum (all peptides) quantification."""
    try:
        from mokume.quantification import AllPeptidesQuantification

        sum_quant = AllPeptidesQuantification()
        result = sum_quant.quantify(
            peptide_df,
            protein_column="ProteinName",
            peptide_column="PeptideSequence",
            intensity_column="NormIntensity",
            sample_column="SampleID",
        )
        return result

    except Exception as e:
        print(f"  ERROR running Sum: {e}")
        return None


def run_directlfq(peptide_df: pd.DataFrame) -> Optional[pd.DataFrame]:
    """
    Run DirectLFQ quantification directly.

    DirectLFQ uses hierarchical normalization with variance-guided
    pairwise alignment for accurate label-free quantification.
    """
    try:
        from mokume.quantification import DirectLFQQuantification

        directlfq = DirectLFQQuantification(min_nonan=1)
        result = directlfq.quantify(
            peptide_df,
            protein_column="ProteinName",
            peptide_column="PeptideSequence",
            intensity_column="NormIntensity",
            sample_column="SampleID",
        )
        return result

    except Exception as e:
        print(f"  ERROR running DirectLFQ: {e}")
        import traceback
        traceback.print_exc()
        return None


def quantify_dataset(
    dataset: DatasetInfo,
    methods: List[str] = None,
    force: bool = False,
) -> Dict[str, bool]:
    """
    Run all quantification methods on a dataset.

    Returns dict mapping method name to success status.
    """
    if methods is None:
        methods = QUANTIFICATION_METHODS

    print(f"\n{dataset.project_id}: {dataset.name}")

    # Load peptide data
    peptide_path = dataset.peptide_path
    if not peptide_path.exists():
        print(f"  ERROR: Peptide file not found: {peptide_path}")
        print("  Run 02_prepare_peptides.py first")
        return {m: False for m in methods}

    peptide_df = pd.read_parquet(peptide_path)
    print(f"  Loaded {len(peptide_df):,} peptide measurements")
    print(f"  Proteins: {peptide_df['ProteinName'].nunique()}, Samples: {peptide_df['SampleID'].nunique()}")

    results = {}

    # Create output directory
    output_dir = PROTEIN_QUANT_DIR / dataset.project_id
    output_dir.mkdir(parents=True, exist_ok=True)

    for method in methods:
        output_path = output_dir / f"{method}.parquet"

        if output_path.exists() and not force:
            print(f"  {method}: Already computed")
            results[method] = True
            continue

        print(f"  {method}: Running...", end=" ", flush=True)
        start_time = time.time()

        if method == "ibaq":
            result_df = run_ibaq(peptide_df, FASTA_FILE, dataset.project_id)
        elif method == "ibaq_raw":
            # ibaq_raw uses the same iBAQ computation - reuse if already computed
            ibaq_path = output_dir / "ibaq.parquet"
            if ibaq_path.exists():
                result_df = pd.read_parquet(ibaq_path)
            else:
                result_df = run_ibaq(peptide_df, FASTA_FILE, dataset.project_id)
                # Also save as ibaq
                if result_df is not None:
                    result_df.to_parquet(ibaq_path, index=False)
        elif method == "maxlfq":
            result_df = run_maxlfq(peptide_df)
        elif method == "directlfq":
            result_df = run_directlfq(peptide_df)
        elif method == "top3":
            result_df = run_top3(peptide_df)
        elif method == "topn":
            # Run multiple TopN values
            for n in TOPN_VALUES:
                topn_output = output_dir / f"top{n}.parquet"
                if topn_output.exists() and not force:
                    continue
                topn_result = run_topn(peptide_df, n=n)
                if topn_result is not None:
                    topn_result.to_parquet(topn_output, index=False)
            result_df = run_topn(peptide_df, n=10)  # Default TopN
        elif method == "sum":
            result_df = run_sum(peptide_df)
        else:
            print(f"Unknown method: {method}")
            results[method] = False
            continue

        elapsed = time.time() - start_time

        if result_df is not None and len(result_df) > 0:
            result_df.to_parquet(output_path, index=False)
            n_proteins = result_df["ProteinName"].nunique()
            print(f"Done ({elapsed:.1f}s, {n_proteins} proteins)")
            results[method] = True
        else:
            print(f"Failed ({elapsed:.1f}s)")
            results[method] = False

    return results


def quantify_all_datasets(
    datasets: dict = None,
    methods: List[str] = None,
    force: bool = False,
) -> dict:
    """
    Run quantification on all datasets.

    Returns nested dict: {dataset_id: {method: success}}
    """
    if datasets is None:
        datasets = ALL_DATASETS
    if methods is None:
        methods = QUANTIFICATION_METHODS

    print("=" * 70)
    print("Protein Quantification")
    print("=" * 70)
    print(f"\nDatasets: {len(datasets)}")
    print(f"Methods: {', '.join(methods)}")
    print(f"Output directory: {PROTEIN_QUANT_DIR}")
    print(f"FASTA file: {FASTA_FILE}")

    # Check FASTA
    if "ibaq" in methods and not FASTA_FILE.exists():
        print(f"\nWARNING: FASTA file not found. iBAQ will be skipped.")

    all_results = {}

    for dataset_id, dataset in datasets.items():
        results = quantify_dataset(dataset, methods=methods, force=force)
        all_results[dataset_id] = results

    # Summary
    print("\n" + "=" * 70)
    print("Quantification Summary")
    print("=" * 70)

    # Create summary table
    print(f"\n{'Dataset':<20}", end="")
    for method in methods:
        print(f"{method:>10}", end="")
    print()
    print("-" * (20 + 10 * len(methods)))

    for dataset_id, results in all_results.items():
        print(f"{dataset_id:<20}", end="")
        for method in methods:
            status = "OK" if results.get(method, False) else "FAIL"
            print(f"{status:>10}", end="")
        print()

    return all_results


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Run protein quantification methods"
    )
    parser.add_argument(
        "--dataset",
        type=str,
        help="Process specific dataset (project ID)"
    )
    parser.add_argument(
        "--method",
        type=str,
        choices=QUANTIFICATION_METHODS,
        help="Run specific method only"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Recompute even if output exists"
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

    # Select methods
    if args.method:
        methods = [args.method]
    else:
        methods = QUANTIFICATION_METHODS

    # Run
    results = quantify_all_datasets(
        datasets=datasets,
        methods=methods,
        force=args.force
    )

    # Count successes
    total_success = sum(
        sum(1 for v in r.values() if v)
        for r in results.values()
    )
    total_attempts = sum(len(r) for r in results.values())

    print(f"\nTotal: {total_success}/{total_attempts} successful")
    sys.exit(0 if total_success > 0 else 1)


if __name__ == "__main__":
    main()
