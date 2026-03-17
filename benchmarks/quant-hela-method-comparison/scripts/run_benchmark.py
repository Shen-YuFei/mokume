#!/usr/bin/env python3
"""
HeLa Benchmark Runner

Master script to run the complete benchmarking pipeline:
1. Download data from PRIDE
2. Prepare peptide data
3. Run quantification methods
4. Compute evaluation metrics
5. Generate visualizations

Usage:
    python run_benchmark.py              # Run all phases
    python run_benchmark.py --phase 3    # Run from phase 3 onwards
    python run_benchmark.py --phase 3 --stop 3  # Run only phase 3
"""

import argparse
from pathlib import Path

# Phase modules
from config import ALL_DATASETS, HELA_DATASETS, TMT_LFQ_COMPARISON

# Use importlib to handle numbered module names
import importlib.util


def import_module(name, path):
    """Import a module from a file path."""
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_phase_1(datasets, force=False):
    """Phase 1: Download data."""
    print("\n" + "=" * 70)
    print("PHASE 1: Data Acquisition")
    print("=" * 70)

    script_dir = Path(__file__).parent
    mod = import_module("download_data", script_dir / "01_download_data.py")
    results = mod.download_all_datasets(datasets=datasets)
    return results


def run_phase_2(datasets, force=False):
    """Phase 2: Prepare peptides."""
    print("\n" + "=" * 70)
    print("PHASE 2: Peptide Preparation")
    print("=" * 70)

    script_dir = Path(__file__).parent
    mod = import_module("prepare_peptides", script_dir / "02_prepare_peptides.py")
    results = mod.process_all_datasets(datasets=datasets, force=force)
    return results


def run_phase_3(datasets, force=False):
    """Phase 3: Run quantification."""
    print("\n" + "=" * 70)
    print("PHASE 3: Protein Quantification")
    print("=" * 70)

    script_dir = Path(__file__).parent
    mod = import_module("run_quantification", script_dir / "03_run_quantification.py")
    results = mod.quantify_all_datasets(datasets=datasets, force=force)
    return results


def run_phase_4(datasets, force=False):
    """Phase 4: Compute metrics."""
    print("\n" + "=" * 70)
    print("PHASE 4: Evaluation Metrics")
    print("=" * 70)

    script_dir = Path(__file__).parent
    mod = import_module("compute_metrics", script_dir / "04_compute_metrics.py")
    results = mod.run_all_metrics(datasets=datasets)
    mod.save_results(results)
    return results


def run_phase_5(force=False):
    """Phase 5: Generate plots."""
    print("\n" + "=" * 70)
    print("PHASE 5: Visualizations")
    print("=" * 70)

    script_dir = Path(__file__).parent
    mod = import_module("generate_plots", script_dir / "05_generate_plots.py")
    mod.generate_all_plots()
    return {}


def main():
    parser = argparse.ArgumentParser(
        description="Run HeLa benchmarking pipeline"
    )
    parser.add_argument(
        "--phase",
        type=int,
        choices=[1, 2, 3, 4, 5],
        default=1,
        help="Start from this phase (default: 1)"
    )
    parser.add_argument(
        "--stop",
        type=int,
        choices=[1, 2, 3, 4, 5],
        default=5,
        help="Stop at this phase (default: 5)"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Reprocess even if output exists"
    )
    parser.add_argument(
        "--hela-only",
        action="store_true",
        help="Only use HeLa datasets (skip PXD007683)"
    )
    parser.add_argument(
        "--comparison-only",
        action="store_true",
        help="Only use PXD007683 TMT/LFQ comparison"
    )

    args = parser.parse_args()

    # Select datasets
    if args.hela_only:
        datasets = HELA_DATASETS
    elif args.comparison_only:
        datasets = TMT_LFQ_COMPARISON
    else:
        datasets = ALL_DATASETS

    print("=" * 70)
    print("HeLa Protein Quantification Benchmark")
    print("=" * 70)
    print(f"\nDatasets: {len(datasets)}")
    print(f"Phases: {args.phase} to {args.stop}")
    print(f"Force recompute: {args.force}")

    # Run phases
    phases = {
        1: ("Data Acquisition", run_phase_1),
        2: ("Peptide Preparation", run_phase_2),
        3: ("Protein Quantification", run_phase_3),
        4: ("Evaluation Metrics", run_phase_4),
        5: ("Visualizations", run_phase_5),
    }

    for phase_num in range(args.phase, args.stop + 1):
        name, func = phases[phase_num]

        if phase_num == 5:
            func(force=args.force)
        else:
            func(datasets, force=args.force)

    print("\n" + "=" * 70)
    print("BENCHMARK COMPLETE")
    print("=" * 70)
    print("\nResults saved to:")
    print("  - Analysis: hela/analysis/")
    print("  - Plots: hela/plots/")


if __name__ == "__main__":
    # Change to script directory for relative imports
    import os
    os.chdir(Path(__file__).parent)

    main()
