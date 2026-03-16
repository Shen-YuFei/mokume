#!/usr/bin/env python3
"""
Stability metrics analysis (CV analysis).

This script answers Q1: What is the best absolute expression algorithm?

Metrics:
- CV within conditions (intra-condition variability)
- CV across all samples
- CV distribution by quantification method
"""

import sys
from pathlib import Path
import warnings

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Add parent to path for config import
sys.path.insert(0, str(Path(__file__).parent))
from config import (
    RESULTS_DIR, FIGURES_DIR,
    SAMPLE_CONDITIONS, CONDITION_COLORS, METHOD_COLORS,
    CV_THRESHOLDS, LOCAL_QUANTIFIED_DIR,
    FIGURE_DPI, FIGURE_FORMAT,
)

warnings.filterwarnings("ignore")


def load_quantified_data(technology: str, method: str = "maxlfq") -> pd.DataFrame:
    """Load pre-quantified data from benchmarks-local."""
    parquet_path = LOCAL_QUANTIFIED_DIR / technology / f"{method}.parquet"

    if not parquet_path.exists():
        raise FileNotFoundError(f"Data not found: {parquet_path}")

    df = pd.read_parquet(parquet_path)

    # Pivot to wide format (proteins x samples)
    df_wide = df.pivot(index="ProteinName", columns="SampleID", values="Intensity")

    return df_wide


def get_condition(sample_id: str) -> str:
    """Get condition for a sample ID."""
    return SAMPLE_CONDITIONS.get(str(sample_id), "Unknown")


def compute_cv_metrics(df_wide: pd.DataFrame) -> dict:
    """Compute CV metrics for a quantified dataset."""
    # Map columns to conditions
    col_to_condition = {str(col): get_condition(col) for col in df_wide.columns}

    # Overall CV (across all samples)
    overall_cv = df_wide.std(axis=1) / df_wide.mean(axis=1)
    overall_cv = overall_cv.replace([np.inf, -np.inf], np.nan).dropna()

    # CV within conditions
    condition_cvs = {}
    for condition in set(col_to_condition.values()):
        if condition == "Unknown":
            continue
        cols = [c for c in df_wide.columns if col_to_condition.get(str(c)) == condition]
        if len(cols) >= 2:
            subset = df_wide[cols]
            cv = subset.std(axis=1) / subset.mean(axis=1)
            cv = cv.replace([np.inf, -np.inf], np.nan).dropna()
            condition_cvs[condition] = cv

    # Combine all within-condition CVs
    all_within_cvs = pd.concat(condition_cvs.values()) if condition_cvs else pd.Series()

    # Quality assessment
    def assess_quality(cv_series):
        if len(cv_series) == 0:
            return {}
        return {
            "excellent": (cv_series < CV_THRESHOLDS["excellent"]).mean(),
            "good": (cv_series < CV_THRESHOLDS["good"]).mean(),
            "acceptable": (cv_series < CV_THRESHOLDS["acceptable"]).mean(),
            "poor": (cv_series >= CV_THRESHOLDS["poor"]).mean(),
        }

    return {
        "n_proteins": len(df_wide),
        "overall_cv_median": overall_cv.median(),
        "overall_cv_mean": overall_cv.mean(),
        "overall_cv_std": overall_cv.std(),
        "within_condition_cv_median": all_within_cvs.median() if len(all_within_cvs) > 0 else np.nan,
        "within_condition_cv_mean": all_within_cvs.mean() if len(all_within_cvs) > 0 else np.nan,
        "overall_cv_distribution": overall_cv,
        "within_cv_distribution": all_within_cvs,
        "condition_cv_distributions": condition_cvs,
        "quality_overall": assess_quality(overall_cv),
        "quality_within": assess_quality(all_within_cvs),
    }


def compare_methods(technology: str) -> pd.DataFrame:
    """Compare CV metrics across different quantification methods."""
    # Available methods in the pre-quantified data
    methods = ["sum", "top3", "top5", "top10", "maxlfq", "directlfq"]

    results = []

    for method in methods:
        try:
            df_wide = load_quantified_data(technology, method=method)
            metrics = compute_cv_metrics(df_wide)

            results.append({
                "technology": technology,
                "method": method,
                "n_proteins": metrics["n_proteins"],
                "overall_cv_median": metrics["overall_cv_median"],
                "overall_cv_mean": metrics["overall_cv_mean"],
                "within_cv_median": metrics["within_condition_cv_median"],
                "within_cv_mean": metrics["within_condition_cv_mean"],
                "pct_excellent": metrics["quality_within"].get("excellent", 0) * 100,
                "pct_good": metrics["quality_within"].get("good", 0) * 100,
            })
            print(f"  {method}: CV={metrics['within_condition_cv_median']:.3f}")

        except FileNotFoundError:
            print(f"  {method}: not available")
            continue
        except Exception as e:
            print(f"  {method}: ERROR - {e}")
            continue

    return pd.DataFrame(results)


def plot_cv_comparison(results: pd.DataFrame, output_path: Path):
    """Plot CV comparison across methods."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, tech in zip(axes, ["tmt", "lfq"]):
        tech_data = results[results["technology"] == tech]

        if len(tech_data) == 0:
            continue

        methods = tech_data["method"].tolist()
        x = np.arange(len(methods))
        width = 0.35

        # Overall CV
        ax.bar(x - width/2, tech_data["overall_cv_median"], width,
               label="Overall CV", color="steelblue", alpha=0.7)

        # Within-condition CV
        ax.bar(x + width/2, tech_data["within_cv_median"], width,
               label="Within-condition CV", color="coral", alpha=0.7)

        ax.set_xlabel("Method")
        ax.set_ylabel("Median CV")
        ax.set_title(f"{tech.upper()}")
        ax.set_xticks(x)
        ax.set_xticklabels(methods, rotation=45, ha="right")
        ax.legend()
        ax.grid(True, alpha=0.3, axis="y")

        # Add threshold lines
        ax.axhline(CV_THRESHOLDS["excellent"], color="green", linestyle="--", alpha=0.5, label="_nolegend_")
        ax.axhline(CV_THRESHOLDS["good"], color="orange", linestyle="--", alpha=0.5, label="_nolegend_")

    plt.suptitle("CV Comparison by Quantification Method", fontsize=14)
    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()


def plot_cv_distribution(metrics: dict, technology: str, method: str, output_path: Path):
    """Plot CV distribution histogram."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Overall CV distribution
    ax = axes[0]
    cv_data = metrics["overall_cv_distribution"]
    ax.hist(cv_data, bins=50, alpha=0.7, color="steelblue", edgecolor="white")
    ax.axvline(cv_data.median(), color="red", linestyle="--", linewidth=2,
               label=f"Median: {cv_data.median():.3f}")
    ax.axvline(CV_THRESHOLDS["good"], color="orange", linestyle="--",
               label=f"Good threshold: {CV_THRESHOLDS['good']}")
    ax.set_xlabel("Coefficient of Variation")
    ax.set_ylabel("Count")
    ax.set_title("Overall CV Distribution")
    ax.legend()
    ax.set_xlim(0, 1)

    # Within-condition CV distribution
    ax = axes[1]
    cv_data = metrics["within_cv_distribution"]
    if len(cv_data) > 0:
        ax.hist(cv_data, bins=50, alpha=0.7, color="coral", edgecolor="white")
        ax.axvline(cv_data.median(), color="red", linestyle="--", linewidth=2,
                   label=f"Median: {cv_data.median():.3f}")
        ax.axvline(CV_THRESHOLDS["good"], color="orange", linestyle="--",
                   label=f"Good threshold: {CV_THRESHOLDS['good']}")
    ax.set_xlabel("Coefficient of Variation")
    ax.set_ylabel("Count")
    ax.set_title("Within-Condition CV Distribution")
    ax.legend()
    ax.set_xlim(0, 1)

    plt.suptitle(f"{technology.upper()} - {method}: CV Distributions", fontsize=12)
    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()


def plot_cv_by_condition(metrics: dict, technology: str, output_path: Path):
    """Plot CV distribution by condition."""
    condition_cvs = metrics["condition_cv_distributions"]

    if not condition_cvs:
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    positions = []
    labels = []

    for i, (condition, cvs) in enumerate(sorted(condition_cvs.items())):
        bp = ax.boxplot(cvs.dropna(), positions=[i], widths=0.6,
                       patch_artist=True)
        color = CONDITION_COLORS.get(condition, "gray")
        for patch in bp["boxes"]:
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        positions.append(i)
        labels.append(f"{condition}\n(n={len(cvs)})")

    ax.set_xticks(positions)
    ax.set_xticklabels(labels)
    ax.set_ylabel("CV")
    ax.set_title(f"{technology.upper()}: CV Distribution by Condition")
    ax.grid(True, alpha=0.3, axis="y")

    # Add threshold lines
    ax.axhline(CV_THRESHOLDS["excellent"], color="green", linestyle="--", alpha=0.5, label="Excellent (<10%)")
    ax.axhline(CV_THRESHOLDS["good"], color="orange", linestyle="--", alpha=0.5, label="Good (<20%)")
    ax.legend()

    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()


def main():
    """Run stability metrics analysis."""
    print("=" * 60)
    print("PXD007683 Benchmark: Stability Metrics (CV Analysis)")
    print("=" * 60)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    all_results = []

    for technology in ["tmt", "lfq"]:
        print(f"\n{technology.upper()}:")
        print("-" * 40)

        try:
            # Compare all methods
            print("  Comparing methods:")
            results = compare_methods(technology)
            all_results.append(results)

            # Detailed analysis for best method (maxlfq)
            print("\n  Detailed analysis (MaxLFQ):")
            df_wide = load_quantified_data(technology, method="maxlfq")
            metrics = compute_cv_metrics(df_wide)

            print(f"    Proteins: {metrics['n_proteins']}")
            print(f"    Overall CV (median): {metrics['overall_cv_median']:.3f}")
            print(f"    Within-condition CV (median): {metrics['within_condition_cv_median']:.3f}")
            print(f"    Quality: {metrics['quality_within']['good']*100:.1f}% proteins with CV < 20%")

            # Generate plots
            plot_cv_distribution(
                metrics, technology, "maxlfq",
                FIGURES_DIR / f"cv_distribution_{technology}.{FIGURE_FORMAT}"
            )

            plot_cv_by_condition(
                metrics, technology,
                FIGURES_DIR / f"cv_by_condition_{technology}.{FIGURE_FORMAT}"
            )

        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
            continue

    # Save and plot combined results
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined.to_csv(RESULTS_DIR / "stability_metrics.csv", index=False)
        print(f"\nResults saved to: {RESULTS_DIR / 'stability_metrics.csv'}")

        plot_cv_comparison(
            combined,
            FIGURES_DIR / f"cv_method_comparison.{FIGURE_FORMAT}"
        )
        print(f"Comparison plot saved to: {FIGURES_DIR / 'cv_method_comparison.png'}")

    # Summary
    print("\n" + "=" * 60)
    print("Summary: Best Methods by CV")
    print("=" * 60)

    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        for tech in ["tmt", "lfq"]:
            tech_data = combined[combined["technology"] == tech]
            if len(tech_data) > 0:
                best = tech_data.loc[tech_data["within_cv_median"].idxmin()]
                print(f"\n{tech.upper()}:")
                print(f"  Best method: {best['method']}")
                print(f"  Within-condition CV: {best['within_cv_median']:.3f}")
                print(f"  Proteins with CV < 20%: {best['pct_good']:.1f}%")


if __name__ == "__main__":
    main()
