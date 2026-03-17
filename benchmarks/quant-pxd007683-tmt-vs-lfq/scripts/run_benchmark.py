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

import runpy
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


def _validate_step(step_index: int):
    """Validate step index and return (script_name, description, resolved_path) or None."""
    if not 0 <= step_index < len(STEPS):
        print(f"  Invalid step index: {step_index}")
        return None

    script_name, description = STEPS[step_index]
    script_path = (SCRIPT_DIR / script_name).resolve()

    if not script_path.is_relative_to(SCRIPT_DIR.resolve()):
        print(f"  Script path escapes allowed directory: {script_path}")
        return None

    if not script_path.exists():
        print(f"  Script not found: {script_path}")
        return None

    return script_name, description, script_path


def run_step(step_index: int) -> bool:
    """Run a single benchmark script by its index in STEPS."""
    validated = _validate_step(step_index)
    if validated is None:
        return False

    script_name, description, script_path = validated

    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"Script: {script_name}")
    print("=" * 60)

    try:
        runpy.run_path(str(script_path), run_name="__main__")
        return True
    except SystemExit as e:
        if e.code is not None and e.code != 0:
            print(f"\nERROR: Script exited with code {e.code}")
            return False
        return True
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

    # Determine which step indices to run
    if args.step is not None:
        indices_to_run = [args.step]
        print(f"\nRunning only step {args.step}")
    else:
        indices_to_run = list(range(len(STEPS)))

        if args.skip_download:
            indices_to_run = [i for i in indices_to_run if "download" not in STEPS[i][0].lower()]
            print("\nSkipping download step")

        if args.quick:
            indices_to_run = [i for i in indices_to_run if "grid_search" not in STEPS[i][0].lower()]
            print("\nQuick mode: skipping grid search")

    # Run steps
    results = []
    for idx in indices_to_run:
        success = run_step(idx)
        results.append((STEPS[idx][0], success))

        if not success:
            print(f"\nStep failed: {STEPS[idx][0]}")
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
