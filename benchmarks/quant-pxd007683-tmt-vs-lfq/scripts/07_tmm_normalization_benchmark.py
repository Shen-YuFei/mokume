#!/usr/bin/env python3
"""
TMM Normalization Benchmark.

This script evaluates TMM (Trimmed Mean of M-values) normalization:
- Compares CV before and after TMM normalization
- Tests effect on different quantification methods
- Evaluates robustness to composition bias

TMM is particularly useful for datasets with spike-ins (like PXD007683 with yeast)
where a subset of proteins has different abundance patterns.
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
    SAMPLE_CONDITIONS, CONDITION_COLORS,
    CV_THRESHOLDS, LOCAL_QUANTIFIED_DIR,
    FIGURE_DPI, FIGURE_FORMAT,
    get_methods_for_technology,
)

# Add mokume to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))
from mokume.normalization.tmm import TMMNormalizer, tmm_normalize

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

    return {
        "n_proteins": len(df_wide),
        "overall_cv_median": overall_cv.median(),
        "overall_cv_mean": overall_cv.mean(),
        "within_cv_median": all_within_cvs.median() if len(all_within_cvs) > 0 else np.nan,
        "within_cv_mean": all_within_cvs.mean() if len(all_within_cvs) > 0 else np.nan,
        "pct_good": (all_within_cvs < CV_THRESHOLDS["good"]).mean() * 100 if len(all_within_cvs) > 0 else 0,
        "pct_excellent": (all_within_cvs < CV_THRESHOLDS["excellent"]).mean() * 100 if len(all_within_cvs) > 0 else 0,
        "overall_cv_distribution": overall_cv,
        "within_cv_distribution": all_within_cvs,
    }


def benchmark_tmm(technology: str, method: str = "maxlfq") -> dict:
    """
    Benchmark TMM normalization for a given technology and method.

    Returns metrics before and after TMM normalization.
    """
    # Load data
    df_wide = load_quantified_data(technology, method=method)
    print(f"    Loaded {len(df_wide)} proteins, {len(df_wide.columns)} samples")

    # Replace zeros with NaN for normalization
    df_clean = df_wide.replace(0, np.nan)

    # Drop rows with too many missing values (need at least 3 samples)
    min_samples = 3
    df_clean = df_clean.dropna(thresh=min_samples)
    print(f"    After filtering: {len(df_clean)} proteins")

    # Metrics before TMM
    metrics_before = compute_cv_metrics(df_clean)
    print(f"    Before TMM: CV={metrics_before['within_cv_median']:.4f}")

    # Apply TMM normalization
    try:
        tmm = TMMNormalizer(m_trim=0.3, a_trim=0.05)
        df_tmm = tmm.fit_transform(df_clean)

        # Get normalization factors for reporting
        norm_factors = tmm.get_norm_factors()
        ref_sample = tmm.ref_sample_

        # Metrics after TMM
        metrics_after = compute_cv_metrics(df_tmm)
        print(f"    After TMM: CV={metrics_after['within_cv_median']:.4f}")
        print(f"    Reference sample: {ref_sample}")
        print(f"    Factor range: {norm_factors.min():.3f} - {norm_factors.max():.3f}")

        # Compute improvement
        cv_improvement = (metrics_before['within_cv_median'] - metrics_after['within_cv_median']) / metrics_before['within_cv_median'] * 100

        return {
            "technology": technology,
            "method": method,
            "n_proteins": len(df_clean),
            "ref_sample": ref_sample,
            "factor_min": norm_factors.min(),
            "factor_max": norm_factors.max(),
            "factor_std": norm_factors.std(),
            "cv_before": metrics_before['within_cv_median'],
            "cv_after": metrics_after['within_cv_median'],
            "cv_improvement_pct": cv_improvement,
            "pct_good_before": metrics_before['pct_good'],
            "pct_good_after": metrics_after['pct_good'],
            "overall_cv_before": metrics_before['overall_cv_median'],
            "overall_cv_after": metrics_after['overall_cv_median'],
            "cv_dist_before": metrics_before['within_cv_distribution'],
            "cv_dist_after": metrics_after['within_cv_distribution'],
        }

    except Exception as e:
        print(f"    ERROR applying TMM: {e}")
        import traceback
        traceback.print_exc()
        return None


def plot_tmm_comparison(results: list, output_path: Path):
    """Plot TMM normalization effect across technologies and methods."""
    df = pd.DataFrame([r for r in results if r is not None])

    if len(df) == 0:
        print("No results to plot")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, tech in zip(axes, ["tmt", "lfq"]):
        tech_data = df[df["technology"] == tech]

        if len(tech_data) == 0:
            ax.set_title(f"{tech.upper()}: No data")
            continue

        methods = tech_data["method"].tolist()
        x = np.arange(len(methods))
        width = 0.35

        # Before TMM
        ax.bar(x - width/2, tech_data["cv_before"], width,
               label="Before TMM", color="lightcoral", alpha=0.8)

        # After TMM
        ax.bar(x + width/2, tech_data["cv_after"], width,
               label="After TMM", color="steelblue", alpha=0.8)

        ax.set_xlabel("Quantification Method")
        ax.set_ylabel("Median Within-Condition CV")
        ax.set_title(f"{tech.upper()}: TMM Normalization Effect")
        ax.set_xticks(x)
        ax.set_xticklabels(methods, rotation=45, ha="right")
        ax.legend()
        ax.grid(True, alpha=0.3, axis="y")

        # Add improvement annotations
        for i, row in tech_data.iterrows():
            idx = methods.index(row["method"])
            improvement = row["cv_improvement_pct"]
            color = "green" if improvement > 0 else "red"
            sign = "+" if improvement > 0 else ""
            ax.annotate(f"{sign}{improvement:.1f}%",
                       xy=(idx, max(row["cv_before"], row["cv_after"]) + 0.005),
                       ha="center", fontsize=8, color=color)

    plt.suptitle("Effect of TMM Normalization on CV", fontsize=14)
    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()
    print(f"  Comparison plot saved: {output_path}")


def plot_cv_distribution_comparison(result: dict, output_path: Path):
    """Plot CV distribution before and after TMM."""
    if result is None:
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    cv_before = result["cv_dist_before"]
    cv_after = result["cv_dist_after"]

    # Before TMM
    ax = axes[0]
    ax.hist(cv_before, bins=50, alpha=0.7, color="lightcoral", edgecolor="white")
    ax.axvline(cv_before.median(), color="red", linestyle="--", linewidth=2,
               label=f"Median: {cv_before.median():.3f}")
    ax.axvline(CV_THRESHOLDS["good"], color="orange", linestyle="--",
               label=f"Good threshold: {CV_THRESHOLDS['good']}")
    ax.set_xlabel("CV")
    ax.set_ylabel("Count")
    ax.set_title("Before TMM")
    ax.legend()
    ax.set_xlim(0, 0.5)

    # After TMM
    ax = axes[1]
    ax.hist(cv_after, bins=50, alpha=0.7, color="steelblue", edgecolor="white")
    ax.axvline(cv_after.median(), color="red", linestyle="--", linewidth=2,
               label=f"Median: {cv_after.median():.3f}")
    ax.axvline(CV_THRESHOLDS["good"], color="orange", linestyle="--",
               label=f"Good threshold: {CV_THRESHOLDS['good']}")
    ax.set_xlabel("CV")
    ax.set_ylabel("Count")
    ax.set_title("After TMM")
    ax.legend()
    ax.set_xlim(0, 0.5)

    tech = result["technology"].upper()
    method = result["method"]
    plt.suptitle(f"{tech} ({method}): CV Distribution Before/After TMM", fontsize=12)
    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()


def plot_normalization_factors(results: list, output_path: Path):
    """Plot TMM normalization factors across samples."""
    # Collect norm factors for maxlfq from each technology
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax, tech in zip(axes, ["tmt", "lfq"]):
        # Re-run TMM to get the factors
        try:
            df_wide = load_quantified_data(tech, method="maxlfq")
            df_clean = df_wide.replace(0, np.nan).dropna(thresh=3)

            tmm = TMMNormalizer(m_trim=0.3, a_trim=0.05)
            tmm.fit(df_clean)
            factors = tmm.get_norm_factors()

            # Color by condition
            colors = [CONDITION_COLORS.get(get_condition(s), "gray") for s in factors.index]

            ax.bar(range(len(factors)), factors.values, color=colors, alpha=0.8)
            ax.axhline(1.0, color="black", linestyle="--", linewidth=1, label="No adjustment")
            ax.set_xlabel("Sample")
            ax.set_ylabel("TMM Normalization Factor")
            ax.set_title(f"{tech.upper()}: TMM Factors by Sample")
            ax.set_xticks(range(len(factors)))
            ax.set_xticklabels(factors.index, rotation=45, ha="right")

            # Add reference sample annotation
            ref_idx = list(factors.index).index(tmm.ref_sample_)
            ax.annotate("ref", xy=(ref_idx, factors.iloc[ref_idx]),
                       xytext=(ref_idx, factors.iloc[ref_idx] + 0.05),
                       ha="center", fontsize=8)

        except Exception as e:
            ax.set_title(f"{tech.upper()}: Error - {e}")

    # Add legend for conditions
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=color, label=cond)
                      for cond, color in CONDITION_COLORS.items()]
    axes[0].legend(handles=legend_elements, loc="upper right", fontsize=8)

    plt.suptitle("TMM Normalization Factors", fontsize=14)
    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()
    print(f"  Factors plot saved: {output_path}")


def main():
    """Run TMM normalization benchmark."""
    print("=" * 60)
    print("PXD007683 Benchmark: TMM Normalization")
    print("=" * 60)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    results = []

    for technology in ["tmt", "lfq"]:
        print(f"\n{technology.upper()}:")
        print("-" * 40)

        # Get appropriate methods for this technology
        # DirectLFQ should NOT be used for TMT
        methods = get_methods_for_technology(technology)
        print(f"  Methods: {methods}")

        for method in methods:
            print(f"\n  Method: {method}")
            try:
                result = benchmark_tmm(technology, method)
                if result is not None:
                    results.append(result)

                    # Plot CV distribution for maxlfq
                    if method == "maxlfq":
                        plot_cv_distribution_comparison(
                            result,
                            FIGURES_DIR / f"tmm_cv_dist_{technology}.{FIGURE_FORMAT}"
                        )

            except FileNotFoundError as e:
                print(f"    Skipped: {e}")
            except Exception as e:
                print(f"    ERROR: {e}")
                import traceback
                traceback.print_exc()

    # Save results
    if results:
        # Create summary DataFrame (without distribution columns)
        summary_cols = ["technology", "method", "n_proteins", "ref_sample",
                       "factor_min", "factor_max", "factor_std",
                       "cv_before", "cv_after", "cv_improvement_pct",
                       "pct_good_before", "pct_good_after",
                       "overall_cv_before", "overall_cv_after"]
        summary_df = pd.DataFrame([{k: v for k, v in r.items() if k in summary_cols}
                                   for r in results])
        summary_df.to_csv(RESULTS_DIR / "tmm_benchmark.csv", index=False)
        print(f"\nResults saved to: {RESULTS_DIR / 'tmm_benchmark.csv'}")

        # Generate plots
        plot_tmm_comparison(results, FIGURES_DIR / f"tmm_comparison.{FIGURE_FORMAT}")
        plot_normalization_factors(results, FIGURES_DIR / f"tmm_factors.{FIGURE_FORMAT}")

    # Summary
    print("\n" + "=" * 60)
    print("TMM Normalization Summary")
    print("=" * 60)

    if results:
        summary_df = pd.DataFrame([{k: v for k, v in r.items()
                                    if k not in ["cv_dist_before", "cv_dist_after"]}
                                   for r in results])

        for tech in ["tmt", "lfq"]:
            tech_data = summary_df[summary_df["technology"] == tech]
            if len(tech_data) > 0:
                print(f"\n{tech.upper()}:")
                print(f"  Average CV improvement: {tech_data['cv_improvement_pct'].mean():.1f}%")
                best = tech_data.loc[tech_data["cv_improvement_pct"].idxmax()]
                print(f"  Best improvement: {best['method']} ({best['cv_improvement_pct']:.1f}%)")
                print(f"  Factor variation (maxlfq): {tech_data[tech_data['method']=='maxlfq']['factor_std'].values[0]:.3f}")

        # Overall recommendation
        print("\n" + "-" * 40)
        print("Recommendation:")
        avg_improvement = summary_df["cv_improvement_pct"].mean()
        if avg_improvement > 5:
            print(f"  TMM normalization provides {avg_improvement:.1f}% average CV improvement.")
            print("  Recommended for this dataset with yeast spike-in.")
        elif avg_improvement > 0:
            print(f"  TMM provides modest {avg_improvement:.1f}% CV improvement.")
            print("  May be beneficial for datasets with composition bias.")
        else:
            print(f"  TMM shows {avg_improvement:.1f}% change (no improvement).")
            print("  GlobalMedian or Hierarchical normalization may be more suitable.")


if __name__ == "__main__":
    main()
