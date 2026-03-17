#!/usr/bin/env python3
"""
Fold-change accuracy analysis using yeast spike-in ground truth.

This script answers Q3: Can we accurately detect expected fold-changes?

Ground truth:
- Yeast proteins spiked at 10%, 5%, 3.3% → expected fold-changes: 3x, 2x, 1.5x
- Human proteins should show no change (fold-change = 1.0)

Metrics:
- Observed vs expected fold-change for yeast proteins
- False positive rate for human proteins
- Fold-change accuracy (RMSE, bias)
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
    SAMPLE_CONDITIONS, EXPECTED_FOLD_CHANGES, SPECIES_PATTERNS,
    CONDITION_COLORS, TECHNOLOGY_COLORS,
    LOCAL_QUANTIFIED_DIR, LOCAL_DATA_DIR,
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


def load_species_info(technology: str) -> pd.DataFrame:
    """Load species information from peptide-level data."""
    peptide_path = LOCAL_DATA_DIR / f"{technology}_peptides.parquet"

    if not peptide_path.exists():
        return None

    df = pd.read_parquet(peptide_path)

    # Get protein-level species info (if any peptide is yeast, protein is yeast)
    species = df.groupby("ProteinName").agg({
        "IsYeast": "max",
        "IsHuman": "max",
    }).reset_index()

    return species


def identify_species(protein_name: str) -> str:
    """Identify species from protein name."""
    protein_upper = str(protein_name).upper()

    for species, patterns in SPECIES_PATTERNS.items():
        for pattern in patterns:
            if pattern.upper() in protein_upper:
                return species

    return "unknown"


def get_condition(sample_id: str) -> str:
    """Get condition for a sample ID."""
    return SAMPLE_CONDITIONS.get(str(sample_id), "Unknown")


def compute_fold_changes(df_wide: pd.DataFrame, species_df: pd.DataFrame = None) -> pd.DataFrame:
    """Compute observed fold-changes between conditions."""
    # Map columns to conditions
    col_to_condition = {col: get_condition(col) for col in df_wide.columns}

    # Compute mean per condition
    condition_means = {}
    for condition in set(col_to_condition.values()):
        if condition == "Unknown":
            continue
        cols = [c for c, cond in col_to_condition.items() if cond == condition]
        if cols:
            condition_means[condition] = df_wide[cols].mean(axis=1)

    # Merge species info if available
    if species_df is not None:
        species_map = dict(zip(species_df["ProteinName"], species_df.apply(
            lambda r: "yeast" if r["IsYeast"] else ("human" if r["IsHuman"] else "unknown"),
            axis=1
        )))
    else:
        species_map = None

    # Compute fold-changes for each expected comparison
    results = []
    for (cond1, cond2), expected_fc in EXPECTED_FOLD_CHANGES.items():
        if cond1 not in condition_means or cond2 not in condition_means:
            continue

        mean1 = condition_means[cond1]
        mean2 = condition_means[cond2]

        # Compute fold-change (cond1 / cond2)
        # Use log2 for better interpretation
        log2_fc = np.log2(mean1 / mean2)
        expected_log2_fc = np.log2(expected_fc)

        for protein in df_wide.index:
            if pd.notna(log2_fc.get(protein)):
                # Determine species
                if species_map and protein in species_map:
                    species = species_map[protein]
                else:
                    species = identify_species(protein)

                results.append({
                    "protein": protein,
                    "species": species,
                    "comparison": f"{cond1}_vs_{cond2}",
                    "expected_fc": expected_fc,
                    "expected_log2_fc": expected_log2_fc,
                    "observed_log2_fc": log2_fc[protein],
                    "fc_error": log2_fc[protein] - expected_log2_fc,
                })

    return pd.DataFrame(results)


def compute_accuracy_metrics(fc_df: pd.DataFrame) -> dict:
    """Compute fold-change accuracy metrics."""
    metrics = {}

    for comparison in fc_df["comparison"].unique():
        comp_df = fc_df[fc_df["comparison"] == comparison]

        # Yeast proteins (should show expected fold-change)
        yeast = comp_df[comp_df["species"] == "yeast"]
        if len(yeast) > 0:
            expected = yeast["expected_log2_fc"].iloc[0]

            # RMSE
            rmse = np.sqrt(np.mean((yeast["observed_log2_fc"] - expected) ** 2))

            # Bias
            bias = np.mean(yeast["observed_log2_fc"] - expected)

            # Observed median
            observed_median = yeast["observed_log2_fc"].median()

            # Compression ratio
            compression_ratio = observed_median / expected if expected != 0 else np.nan

            metrics[f"{comparison}_yeast"] = {
                "n_proteins": len(yeast),
                "expected_log2_fc": expected,
                "observed_median_log2_fc": observed_median,
                "rmse": rmse,
                "bias": bias,
                "compression_ratio": compression_ratio,
            }

        # Human proteins (should show no change, FC = 1, log2FC = 0)
        human = comp_df[comp_df["species"] == "human"]
        if len(human) > 0:
            # False positive analysis
            # How many human proteins show significant change?
            # Using a threshold of |log2FC| > 1 (2-fold change)
            false_positives = (np.abs(human["observed_log2_fc"]) > 1).sum()
            fp_rate = false_positives / len(human)

            # RMSE from expected (0)
            rmse = np.sqrt(np.mean(human["observed_log2_fc"] ** 2))

            metrics[f"{comparison}_human"] = {
                "n_proteins": len(human),
                "expected_log2_fc": 0.0,
                "observed_median_log2_fc": human["observed_log2_fc"].median(),
                "rmse": rmse,
                "false_positive_rate": fp_rate,
                "false_positives": false_positives,
            }

    return metrics


def plot_fold_change_accuracy(fc_df: pd.DataFrame, technology: str, output_path: Path):
    """Plot observed vs expected fold-changes."""
    # Filter to major comparison (10% vs 3.3%, expected 3x)
    comparison = "QY_10pct_vs_QY_3pct"
    comp_df = fc_df[fc_df["comparison"] == comparison]

    if len(comp_df) == 0:
        print(f"  No data for comparison {comparison}")
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left: Histogram of fold-changes by species
    ax = axes[0]
    for species, color in [("yeast", CONDITION_COLORS["QY_10pct"]),
                           ("human", CONDITION_COLORS["QY_3pct"])]:
        subset = comp_df[comp_df["species"] == species]
        if len(subset) > 0:
            ax.hist(
                subset["observed_log2_fc"].dropna(),
                bins=50,
                alpha=0.6,
                label=f"{species.capitalize()} (n={len(subset)})",
                color=color,
            )

    # Add expected values
    expected_yeast = np.log2(3.0)  # Expected for 10% vs 3.3%
    ax.axvline(expected_yeast, color="red", linestyle="--", label=f"Expected yeast (log2={expected_yeast:.2f})")
    ax.axvline(0, color="green", linestyle="--", label="Expected human (log2=0)")

    ax.set_xlabel("Observed log2(Fold-Change)")
    ax.set_ylabel("Count")
    ax.set_title(f"{technology.upper()}: 10% vs 3.3% Yeast")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Right: Box plot by species
    ax = axes[1]
    species_data = []
    species_labels = []

    for species in ["yeast", "human"]:
        subset = comp_df[comp_df["species"] == species]["observed_log2_fc"].dropna()
        if len(subset) > 0:
            species_data.append(subset.values)
            species_labels.append(f"{species.capitalize()}\n(n={len(subset)})")

    if species_data:
        bp = ax.boxplot(species_data, labels=species_labels, patch_artist=True)

        colors = [CONDITION_COLORS["QY_10pct"], CONDITION_COLORS["QY_3pct"]]
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)

        ax.axhline(expected_yeast, color="red", linestyle="--", label=f"Expected yeast")
        ax.axhline(0, color="green", linestyle="--", label="Expected human")

    ax.set_ylabel("Observed log2(Fold-Change)")
    ax.set_title(f"{technology.upper()}: Fold-Change Distribution")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()


def plot_fc_comparison(tmt_fc: pd.DataFrame, lfq_fc: pd.DataFrame, output_path: Path):
    """Compare fold-change accuracy between TMT and LFQ."""
    comparison = "QY_10pct_vs_QY_3pct"

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax, (fc_df, tech, color) in zip(
        axes,
        [(tmt_fc, "TMT", TECHNOLOGY_COLORS["TMT"]),
         (lfq_fc, "LFQ", TECHNOLOGY_COLORS["LFQ"])]
    ):
        comp_df = fc_df[(fc_df["comparison"] == comparison) & (fc_df["species"] == "yeast")]

        if len(comp_df) == 0:
            continue

        expected = np.log2(3.0)

        ax.hist(comp_df["observed_log2_fc"].dropna(), bins=50, alpha=0.7, color=color)
        ax.axvline(expected, color="red", linestyle="--", linewidth=2, label=f"Expected ({expected:.2f})")

        median_obs = comp_df["observed_log2_fc"].median()
        ax.axvline(median_obs, color="blue", linestyle="-", linewidth=2,
                   label=f"Observed median ({median_obs:.2f})")

        # Compute compression
        compression = median_obs / expected if expected != 0 else np.nan
        ax.set_title(f"{tech}: Yeast proteins\n(compression: {compression:.2f})")
        ax.set_xlabel("Observed log2(Fold-Change)")
        ax.set_ylabel("Count")
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.suptitle("Fold-Change Accuracy: TMT vs LFQ (10% vs 3.3% yeast)", fontsize=12)
    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()


def main():
    """Run fold-change accuracy analysis."""
    print("=" * 60)
    print("PXD007683 Benchmark: Fold-Change Accuracy")
    print("=" * 60)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    all_fc_data = {}
    all_metrics = {}

    for technology in ["tmt", "lfq"]:
        print(f"\n{technology.upper()}:")
        print("-" * 40)

        try:
            # Load pre-quantified data
            print("  Loading pre-quantified data (MaxLFQ)...")
            df_wide = load_quantified_data(technology, method="maxlfq")
            print(f"  Proteins: {len(df_wide)}")

            # Load species info
            print("  Loading species information...")
            species_df = load_species_info(technology)
            if species_df is not None:
                print(f"  Yeast proteins: {species_df['IsYeast'].sum()}")
                print(f"  Human proteins: {species_df['IsHuman'].sum()}")

            # Compute fold-changes
            print("  Computing fold-changes...")
            fc_df = compute_fold_changes(df_wide, species_df)
            fc_df["technology"] = technology
            all_fc_data[technology] = fc_df

            # Species breakdown
            print(f"  Species in results:")
            for species, count in fc_df["species"].value_counts().items():
                print(f"    {species}: {count}")

            # Compute accuracy metrics
            print("  Computing accuracy metrics...")
            metrics = compute_accuracy_metrics(fc_df)
            all_metrics[technology] = metrics

            # Print summary
            for key, vals in metrics.items():
                print(f"\n  {key}:")
                for k, v in vals.items():
                    if isinstance(v, float):
                        print(f"    {k}: {v:.3f}")
                    else:
                        print(f"    {k}: {v}")

            # Generate plots
            plot_fold_change_accuracy(
                fc_df, technology,
                FIGURES_DIR / f"fold_change_{technology}.{FIGURE_FORMAT}"
            )

        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
            continue

    # Save results
    if all_fc_data:
        combined_fc = pd.concat(all_fc_data.values(), ignore_index=True)
        combined_fc.to_csv(RESULTS_DIR / "fold_change_data.csv", index=False)
        print(f"\nFold-change data saved to: {RESULTS_DIR / 'fold_change_data.csv'}")

    if all_metrics:
        # Flatten metrics for CSV
        metrics_rows = []
        for tech, tech_metrics in all_metrics.items():
            for comparison, vals in tech_metrics.items():
                row = {"technology": tech, "comparison": comparison}
                row.update(vals)
                metrics_rows.append(row)

        metrics_df = pd.DataFrame(metrics_rows)
        metrics_df.to_csv(RESULTS_DIR / "fold_change_metrics.csv", index=False)
        print(f"Metrics saved to: {RESULTS_DIR / 'fold_change_metrics.csv'}")

    # Combined comparison plot
    if len(all_fc_data) == 2:
        plot_fc_comparison(
            all_fc_data["tmt"],
            all_fc_data["lfq"],
            FIGURES_DIR / f"fold_change_comparison.{FIGURE_FORMAT}"
        )
        print(f"Comparison plot saved to: {FIGURES_DIR / 'fold_change_comparison.png'}")

    # Summary
    print("\n" + "=" * 60)
    print("Summary: Fold-Change Accuracy (10% vs 3.3% yeast)")
    print("=" * 60)

    for tech, metrics in all_metrics.items():
        yeast_key = "QY_10pct_vs_QY_3pct_yeast"
        human_key = "QY_10pct_vs_QY_3pct_human"

        if yeast_key in metrics:
            m = metrics[yeast_key]
            print(f"\n{tech.upper()} - Yeast proteins:")
            print(f"  Expected log2FC: {m['expected_log2_fc']:.2f}")
            print(f"  Observed median log2FC: {m['observed_median_log2_fc']:.2f}")
            print(f"  Compression ratio: {m['compression_ratio']:.2f}")
            print(f"  RMSE: {m['rmse']:.3f}")

        if human_key in metrics:
            m = metrics[human_key]
            print(f"\n{tech.upper()} - Human proteins (negative control):")
            print(f"  Expected log2FC: 0.0")
            print(f"  Observed median log2FC: {m['observed_median_log2_fc']:.2f}")
            print(f"  False positive rate (|log2FC| > 1): {m['false_positive_rate']:.1%}")


if __name__ == "__main__":
    main()
