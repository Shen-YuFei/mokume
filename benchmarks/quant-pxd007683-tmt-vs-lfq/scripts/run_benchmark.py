#!/usr/bin/env python3
"""
Run the complete PXD007683 benchmark pipeline.

This script runs all benchmark analyses in sequence:
1. Download data (if needed)
2. Grid search over methods
3. Variance decomposition
4. Fold-change accuracy
5. Stability metrics
6. Cross-technology correlation
7. Generate report

Usage:
    python run_benchmark.py           # Run all steps
    python run_benchmark.py --quick   # Skip grid search (faster)
    python run_benchmark.py --step 3  # Run only step 3
"""

import subprocess
import sys
import argparse
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent

STEPS = [
    ("00_download_data.py", "Download benchmark data"),
    ("01_grid_search_methods.py", "Grid search over methods"),
    ("02_variance_decomposition.py", "Variance decomposition (PCA)"),
    ("03_fold_change_accuracy.py", "Fold-change accuracy analysis"),
    ("04_stability_metrics.py", "Stability metrics (CV analysis)"),
    ("05_cross_technology_correlation.py", "Cross-technology correlation"),
    ("06_generate_report.py", "Generate comprehensive report"),
]


def run_script(script_name: str, description: str) -> bool:
    """Run a single benchmark script."""
    script_path = SCRIPT_DIR / script_name

    if not script_path.exists():
        print(f"  Script not found: {script_path}")
        return False

    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"Script: {script_name}")
    print("=" * 60)

    try:
        result = subprocess.run(
            [sys.executable, str(script_path)],
            cwd=str(SCRIPT_DIR),
            check=True,
        )
        return result.returncode == 0
    except subprocess.CalledProcessError as e:
        print(f"\nERROR: Script failed with return code {e.returncode}")
        return False
    except Exception as e:
        print(f"\nERROR: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(description="Run PXD007683 benchmark pipeline")
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Skip grid search (faster execution)",
    )
    parser.add_argument(
        "--step",
        type=int,
        choices=range(len(STEPS)),
        help="Run only a specific step (0-indexed)",
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Skip data download step",
    )
    args = parser.parse_args()

    print("=" * 60)
    print("PXD007683 Benchmark Pipeline")
    print("=" * 60)
    print("\nSteps:")
    for i, (script, desc) in enumerate(STEPS):
        print(f"  {i}. {desc}")

    # Determine which steps to run
    if args.step is not None:
        steps_to_run = [STEPS[args.step]]
        print(f"\nRunning only step {args.step}")
    else:
        steps_to_run = STEPS.copy()

        if args.skip_download:
            steps_to_run = [(s, d) for s, d in steps_to_run if "download" not in s.lower()]
            print("\nSkipping download step")

        if args.quick:
            steps_to_run = [(s, d) for s, d in steps_to_run if "grid_search" not in s.lower()]
            print("\nQuick mode: skipping grid search")

    # Run steps
    results = []
    for script, description in steps_to_run:
        success = run_script(script, description)
        results.append((script, success))

        if not success:
            print(f"\nStep failed: {script}")
            if input("Continue with remaining steps? (y/n): ").lower() != "y":
                break

    # Summary
    print("\n" + "=" * 60)
    print("Pipeline Summary")
    print("=" * 60)

    for script, success in results:
        status = "OK" if success else "FAILED"
        print(f"  [{status}] {script}")

    failed = sum(1 for _, s in results if not s)
    if failed == 0:
        print("\nAll steps completed successfully!")
        print(f"\nResults are in: {SCRIPT_DIR.parent / 'results'}")
        print(f"Figures are in: {SCRIPT_DIR.parent / 'figures'}")
    else:
        print(f"\n{failed} step(s) failed.")
        sys.exit(1)


if __name__ == "__main__":
    main()
