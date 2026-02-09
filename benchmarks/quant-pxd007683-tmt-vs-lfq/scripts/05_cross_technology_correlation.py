#!/usr/bin/env python3
"""
Cross-technology correlation analysis: LFQ vs TMT.

This script analyzes the correlation between LFQ and TMT quantification
for the same biological samples.

Key analyses:
- Overall correlation (all proteins, all samples)
- Per-protein correlation (same protein across samples)
- Per-sample correlation (all proteins in matched samples)
"""

import sys
from pathlib import Path
import warnings

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# Add parent to path for config import
sys.path.insert(0, str(Path(__file__).parent))
from config import (
    RESULTS_DIR, FIGURES_DIR,
    SAMPLE_CONDITIONS, CONDITION_COLORS, TECHNOLOGY_COLORS,
    LOCAL_QUANTIFIED_DIR,
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


def compute_overall_correlation(tmt_wide: pd.DataFrame, lfq_wide: pd.DataFrame) -> dict:
    """Compute overall correlation between TMT and LFQ."""
    # Find common proteins
    common_proteins = tmt_wide.index.intersection(lfq_wide.index)

    if len(common_proteins) == 0:
        return {"error": "No common proteins found"}

    # Find common samples (by sample number)
    tmt_samples = set(str(s) for s in tmt_wide.columns)
    lfq_samples = set(str(s) for s in lfq_wide.columns)
    common_samples = sorted(tmt_samples & lfq_samples, key=lambda x: int(x) if x.isdigit() else 0)

    if len(common_samples) == 0:
        return {"error": "No matching samples found"}

    # Create aligned matrices
    tmt_aligned = tmt_wide.loc[common_proteins, common_samples]
    lfq_aligned = lfq_wide.loc[common_proteins, common_samples]

    # Log transform for correlation
    tmt_log = np.log10(tmt_aligned.replace(0, np.nan))
    lfq_log = np.log10(lfq_aligned.replace(0, np.nan))

    # Flatten and remove NaN pairs
    tmt_flat = tmt_log.values.flatten()
    lfq_flat = lfq_log.values.flatten()

    valid = ~(np.isnan(tmt_flat) | np.isnan(lfq_flat))
    tmt_valid = tmt_flat[valid]
    lfq_valid = lfq_flat[valid]

    if len(tmt_valid) < 10:
        return {"error": "Not enough valid values for correlation"}

    # Compute correlations
    pearson_r, pearson_p = stats.pearsonr(tmt_valid, lfq_valid)
    spearman_r, spearman_p = stats.spearmanr(tmt_valid, lfq_valid)

    return {
        "n_proteins": len(common_proteins),
        "n_samples": len(common_samples),
        "n_values": len(tmt_valid),
        "pearson_r": pearson_r,
        "pearson_p": pearson_p,
        "spearman_r": spearman_r,
        "spearman_p": spearman_p,
        "tmt_log": tmt_log,
        "lfq_log": lfq_log,
        "common_samples": common_samples,
    }


def compute_per_protein_correlation(tmt_wide: pd.DataFrame, lfq_wide: pd.DataFrame) -> pd.DataFrame:
    """Compute correlation for each protein across samples."""
    common_proteins = tmt_wide.index.intersection(lfq_wide.index)

    tmt_samples = set(str(s) for s in tmt_wide.columns)
    lfq_samples = set(str(s) for s in lfq_wide.columns)
    common_samples = sorted(tmt_samples & lfq_samples, key=lambda x: int(x) if x.isdigit() else 0)

    if len(common_samples) < 3:
        return pd.DataFrame()

    results = []
    for protein in common_proteins:
        tmt_vals = np.log10(tmt_wide.loc[protein, common_samples].replace(0, np.nan))
        lfq_vals = np.log10(lfq_wide.loc[protein, common_samples].replace(0, np.nan))

        valid = ~(tmt_vals.isna() | lfq_vals.isna())
        n_valid = valid.sum()

        if n_valid >= 3:
            r, p = stats.pearsonr(tmt_vals[valid], lfq_vals[valid])
            results.append({
                "protein": protein,
                "pearson_r": r,
                "pearson_p": p,
                "n_samples": n_valid,
            })

    return pd.DataFrame(results)


def compute_per_sample_correlation(tmt_wide: pd.DataFrame, lfq_wide: pd.DataFrame) -> pd.DataFrame:
    """Compute correlation for each sample across proteins."""
    common_proteins = tmt_wide.index.intersection(lfq_wide.index)

    tmt_samples = set(str(s) for s in tmt_wide.columns)
    lfq_samples = set(str(s) for s in lfq_wide.columns)
    common_samples = sorted(tmt_samples & lfq_samples, key=lambda x: int(x) if x.isdigit() else 0)

    results = []
    for sample in common_samples:
        tmt_vals = np.log10(tmt_wide.loc[common_proteins, sample].replace(0, np.nan))
        lfq_vals = np.log10(lfq_wide.loc[common_proteins, sample].replace(0, np.nan))

        valid = ~(tmt_vals.isna() | lfq_vals.isna())
        n_valid = valid.sum()

        if n_valid >= 10:
            r, p = stats.pearsonr(tmt_vals[valid], lfq_vals[valid])
            condition = get_condition(sample)
            results.append({
                "sample": sample,
                "condition": condition,
                "pearson_r": r,
                "pearson_p": p,
                "n_proteins": n_valid,
            })

    return pd.DataFrame(results)


def plot_overall_correlation(correlation_data: dict, output_path: Path):
    """Create scatter plot of overall TMT vs LFQ correlation."""
    tmt_log = correlation_data["tmt_log"]
    lfq_log = correlation_data["lfq_log"]

    fig, ax = plt.subplots(figsize=(8, 8))

    # Flatten and plot
    tmt_flat = tmt_log.values.flatten()
    lfq_flat = lfq_log.values.flatten()

    valid = ~(np.isnan(tmt_flat) | np.isnan(lfq_flat))

    # 2D histogram for density
    h = ax.hist2d(
        tmt_flat[valid], lfq_flat[valid],
        bins=100, cmap="viridis", cmin=1
    )
    plt.colorbar(h[3], ax=ax, label="Count")

    # Diagonal line
    lims = [
        min(np.nanmin(tmt_flat), np.nanmin(lfq_flat)),
        max(np.nanmax(tmt_flat), np.nanmax(lfq_flat))
    ]
    ax.plot(lims, lims, "r--", alpha=0.8, label="y=x")

    # Add correlation text
    r = correlation_data["pearson_r"]
    ax.text(0.05, 0.95, f"Pearson r = {r:.3f}\nn = {correlation_data['n_values']:,}",
            transform=ax.transAxes, fontsize=12, verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

    ax.set_xlabel("TMT (log10 intensity)")
    ax.set_ylabel("LFQ (log10 intensity)")
    ax.set_title("Cross-Technology Correlation: TMT vs LFQ")
    ax.legend(loc="lower right")

    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()


def plot_per_protein_correlation_distribution(per_protein_df: pd.DataFrame, output_path: Path):
    """Plot distribution of per-protein correlations."""
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.hist(per_protein_df["pearson_r"], bins=50, alpha=0.7, color="steelblue", edgecolor="white")

    median_r = per_protein_df["pearson_r"].median()
    ax.axvline(median_r, color="red", linestyle="--", linewidth=2,
               label=f"Median r = {median_r:.3f}")

    # Count proteins with high correlation
    high_corr = (per_protein_df["pearson_r"] > 0.8).sum()
    pct_high = high_corr / len(per_protein_df) * 100

    ax.text(0.05, 0.95, f"Proteins with r > 0.8: {high_corr} ({pct_high:.1f}%)",
            transform=ax.transAxes, fontsize=11, verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

    ax.set_xlabel("Pearson Correlation (per protein)")
    ax.set_ylabel("Count")
    ax.set_title("Per-Protein TMT-LFQ Correlation Distribution")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()


def plot_per_sample_correlation(per_sample_df: pd.DataFrame, output_path: Path):
    """Plot per-sample correlations colored by condition."""
    fig, ax = plt.subplots(figsize=(10, 6))

    # Sort by sample number
    per_sample_df = per_sample_df.copy()
    per_sample_df["sample_num"] = per_sample_df["sample"].apply(lambda x: int(x) if str(x).isdigit() else 0)
    per_sample_df = per_sample_df.sort_values("sample_num")

    for condition, color in CONDITION_COLORS.items():
        mask = per_sample_df["condition"] == condition
        if mask.any():
            subset = per_sample_df[mask]
            ax.bar(
                [str(s) for s in subset["sample"]],
                subset["pearson_r"],
                color=color,
                alpha=0.7,
                label=condition,
                edgecolor="white",
            )

    ax.axhline(0.9, color="green", linestyle="--", alpha=0.5, label="r = 0.9")

    ax.set_xlabel("Sample")
    ax.set_ylabel("Pearson Correlation")
    ax.set_title("Per-Sample TMT-LFQ Correlation")
    ax.legend()
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()


def main():
    """Run cross-technology correlation analysis."""
    print("=" * 60)
    print("PXD007683 Benchmark: Cross-Technology Correlation")
    print("=" * 60)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    try:
        # Load data
        print("\nLoading TMT data (MaxLFQ)...")
        tmt_wide = load_quantified_data("tmt", method="maxlfq")
        print(f"  Proteins: {len(tmt_wide)}, Samples: {len(tmt_wide.columns)}")

        print("\nLoading LFQ data (MaxLFQ)...")
        lfq_wide = load_quantified_data("lfq", method="maxlfq")
        print(f"  Proteins: {len(lfq_wide)}, Samples: {len(lfq_wide.columns)}")

        # Overall correlation
        print("\nComputing overall correlation...")
        overall = compute_overall_correlation(tmt_wide, lfq_wide)

        if "error" in overall:
            print(f"  ERROR: {overall['error']}")
            return

        print(f"  Common proteins: {overall['n_proteins']}")
        print(f"  Matched samples: {overall['n_samples']}")
        print(f"  Pearson r: {overall['pearson_r']:.4f}")
        print(f"  Spearman r: {overall['spearman_r']:.4f}")

        # Per-protein correlation
        print("\nComputing per-protein correlations...")
        per_protein = compute_per_protein_correlation(tmt_wide, lfq_wide)
        print(f"  Proteins analyzed: {len(per_protein)}")
        print(f"  Median correlation: {per_protein['pearson_r'].median():.4f}")
        print(f"  Proteins with r > 0.8: {(per_protein['pearson_r'] > 0.8).sum()} "
              f"({(per_protein['pearson_r'] > 0.8).mean()*100:.1f}%)")

        # Per-sample correlation
        print("\nComputing per-sample correlations...")
        per_sample = compute_per_sample_correlation(tmt_wide, lfq_wide)
        print(f"  Samples analyzed: {len(per_sample)}")
        print(f"  Mean correlation: {per_sample['pearson_r'].mean():.4f}")
        print(f"  Range: {per_sample['pearson_r'].min():.4f} - {per_sample['pearson_r'].max():.4f}")

        # Save results
        results_summary = pd.DataFrame([{
            "metric": "overall",
            "n_proteins": overall["n_proteins"],
            "n_samples": overall["n_samples"],
            "pearson_r": overall["pearson_r"],
            "spearman_r": overall["spearman_r"],
        }])
        results_summary.to_csv(RESULTS_DIR / "cross_technology_overall.csv", index=False)

        per_protein.to_csv(RESULTS_DIR / "cross_technology_per_protein.csv", index=False)
        per_sample.to_csv(RESULTS_DIR / "cross_technology_per_sample.csv", index=False)

        print(f"\nResults saved to: {RESULTS_DIR}")

        # Generate plots
        print("\nGenerating plots...")

        plot_overall_correlation(
            overall,
            FIGURES_DIR / f"cross_tech_overall.{FIGURE_FORMAT}"
        )

        plot_per_protein_correlation_distribution(
            per_protein,
            FIGURES_DIR / f"cross_tech_per_protein.{FIGURE_FORMAT}"
        )

        plot_per_sample_correlation(
            per_sample,
            FIGURES_DIR / f"cross_tech_per_sample.{FIGURE_FORMAT}"
        )

        print(f"Plots saved to: {FIGURES_DIR}")

        # Summary
        print("\n" + "=" * 60)
        print("Summary: Cross-Technology Correlation")
        print("=" * 60)
        print(f"\nOverall Pearson correlation: {overall['pearson_r']:.4f}")
        print(f"Overall Spearman correlation: {overall['spearman_r']:.4f}")
        print(f"\nPer-protein median correlation: {per_protein['pearson_r'].median():.4f}")
        print(f"Per-sample mean correlation: {per_sample['pearson_r'].mean():.4f}")

        # Interpretation
        print("\nInterpretation:")
        if overall['pearson_r'] > 0.9:
            print("  Excellent agreement between TMT and LFQ")
        elif overall['pearson_r'] > 0.8:
            print("  Good agreement between TMT and LFQ")
        elif overall['pearson_r'] > 0.7:
            print("  Moderate agreement between TMT and LFQ")
        else:
            print("  Poor agreement between TMT and LFQ - further investigation needed")

    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
