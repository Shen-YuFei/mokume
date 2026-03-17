#!/usr/bin/env python3
"""
Phase 4: Preprocessing Filters Benchmark.

Tests filtering stringency at peptide and protein levels.

Filters tested:
- Peptide length: 7 (default), 8, 9 amino acids minimum
- Missed cleavages: 0, 1, 2 (trypsin)
- Minimum unique peptides per protein: 1, 2, 3

Note: FDR/score filters not available as data is already pre-filtered by quantms.
"""

import sys
from pathlib import Path
import warnings
from itertools import product

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Add parent to path for config import
sys.path.insert(0, str(Path(__file__).parent))
from config import (
    RESULTS_DIR, FIGURES_DIR,
    SAMPLE_CONDITIONS,
    CV_THRESHOLDS,
    FIGURE_DPI, FIGURE_FORMAT,
    EXPECTED_FOLD_CHANGES, SPECIES_PATTERNS,
)

# Add mokume to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

# Import mokume filters
from mokume.preprocessing.filters.peptide import (
    PeptideLengthFilter,
    MissedCleavageFilter,
)
from mokume.preprocessing.filters.protein import MinPeptideFilter

warnings.filterwarnings("ignore")

# =============================================================================
# Paths
# =============================================================================

LOCAL_DIR = Path(__file__).parent.parent.parent.parent / "benchmarks-local" / "PXD007683"
LOCAL_DATA_DIR = LOCAL_DIR / "data" / "processed"

# =============================================================================
# Filter Configurations to Test
# =============================================================================

# Peptide length filters
PEPTIDE_LENGTH_CONFIGS = [
    {"min_length": 7, "max_length": 50, "name": "len≥7"},   # Default
    {"min_length": 8, "max_length": 50, "name": "len≥8"},   # Stricter
    {"min_length": 9, "max_length": 50, "name": "len≥9"},   # Most strict
]

# Missed cleavage filters
MISSED_CLEAVAGE_CONFIGS = [
    {"max_missed_cleavages": 2, "name": "mc≤2"},  # Permissive
    {"max_missed_cleavages": 1, "name": "mc≤1"},  # Default
    {"max_missed_cleavages": 0, "name": "mc=0"},  # Strict
]

# Minimum peptides per protein
MIN_PEPTIDE_CONFIGS = [
    {"min_unique_peptides": 1, "name": "pep≥1"},  # Permissive
    {"min_unique_peptides": 2, "name": "pep≥2"},  # Default
    {"min_unique_peptides": 3, "name": "pep≥3"},  # Strict
]

# Default quantification method per technology
DEFAULT_QUANT = {
    "tmt": "sum",
    "lfq": "directlfq"
}


def load_peptide_data(technology: str) -> pd.DataFrame:
    """Load peptide-level data."""
    parquet_path = LOCAL_DATA_DIR / f"{technology}_peptides.parquet"

    if not parquet_path.exists():
        raise FileNotFoundError(f"Data not found: {parquet_path}")

    return pd.read_parquet(parquet_path)


def get_condition(sample_id: str) -> str:
    """Get condition for a sample ID."""
    return SAMPLE_CONDITIONS.get(str(sample_id), "Unknown")


def get_species(protein_name: str) -> str:
    """Identify species from protein name."""
    protein_upper = protein_name.upper()
    for species, patterns in SPECIES_PATTERNS.items():
        for pattern in patterns:
            if pattern.upper() in protein_upper:
                return species
    return "unknown"


# =============================================================================
# Filter Application
# =============================================================================

def apply_filters(
    df: pd.DataFrame,
    min_length: int = 7,
    max_missed_cleavages: int = 2,
    min_unique_peptides: int = 1
) -> pd.DataFrame:
    """Apply preprocessing filters to peptide data."""
    df_filtered = df.copy()

    # Peptide length filter
    length_filter = PeptideLengthFilter(
        min_length=min_length,
        max_length=50,
        sequence_column="PeptideSequence"
    )
    df_filtered, _ = length_filter.apply(df_filtered)

    # Missed cleavage filter
    mc_filter = MissedCleavageFilter(
        max_missed_cleavages=max_missed_cleavages,
        sequence_column="PeptideSequence",
        enzyme="trypsin"
    )
    df_filtered, _ = mc_filter.apply(df_filtered)

    # Minimum peptides filter
    pep_filter = MinPeptideFilter(
        min_peptides=1,
        min_unique_peptides=min_unique_peptides,
        protein_column="ProteinName",
        peptide_column="PeptideSequence"
    )
    df_filtered, _ = pep_filter.apply(df_filtered)

    return df_filtered


def quantify_sum(df: pd.DataFrame) -> pd.DataFrame:
    """Simple sum quantification for filtered peptides."""
    protein_intensities = df.groupby(["ProteinName", "SampleID"])["NormIntensity"].sum().reset_index()
    protein_intensities.columns = ["ProteinName", "SampleID", "Intensity"]
    df_wide = protein_intensities.pivot(index="ProteinName", columns="SampleID", values="Intensity")
    return df_wide


# =============================================================================
# Evaluation Metrics
# =============================================================================

def compute_metrics(df: pd.DataFrame) -> dict:
    """Compute all evaluation metrics."""
    col_to_condition = {str(col): get_condition(col) for col in df.columns}

    # Q1: Within-condition CV
    condition_cvs = []
    for condition in set(col_to_condition.values()):
        if condition == "Unknown":
            continue
        cols = [c for c in df.columns if col_to_condition.get(str(c)) == condition]
        if len(cols) >= 2:
            subset = df[cols]
            cv = subset.std(axis=1) / subset.mean(axis=1)
            cv = cv.replace([np.inf, -np.inf], np.nan).dropna()
            condition_cvs.extend(cv.tolist())

    within_cv = np.median(condition_cvs) if condition_cvs else np.nan

    # Q3: Fold-change accuracy for yeast proteins
    yeast_proteins = [p for p in df.index if get_species(p) == "yeast"]

    fc_rmse = np.nan
    if yeast_proteins:
        yeast_df = df.loc[yeast_proteins]
        condition_intensities = {}
        for condition in ["QY_10pct", "QY_5pct", "QY_3pct"]:
            cols = [c for c in yeast_df.columns if col_to_condition.get(str(c)) == condition]
            if cols:
                condition_intensities[condition] = yeast_df[cols].mean(axis=1)

        if len(condition_intensities) >= 2:
            observed_fcs = []
            expected_fcs = []
            for (cond1, cond2), expected_fc in EXPECTED_FOLD_CHANGES.items():
                if cond1 in condition_intensities and cond2 in condition_intensities:
                    obs_fc = (condition_intensities[cond1] / condition_intensities[cond2]).dropna()
                    obs_log2_fc = np.log2(obs_fc)
                    expected_log2_fc = np.log2(expected_fc)
                    observed_fcs.extend(obs_log2_fc.tolist())
                    expected_fcs.extend([expected_log2_fc] * len(obs_log2_fc))

            if observed_fcs:
                observed_fcs = np.array(observed_fcs)
                expected_fcs = np.array(expected_fcs)
                fc_rmse = np.sqrt(np.mean((observed_fcs - expected_fcs) ** 2))

    # Quality metrics
    n_proteins = len(df)
    pct_good_cv = (np.array(condition_cvs) < CV_THRESHOLDS["good"]).mean() * 100 if condition_cvs else 0

    return {
        "n_proteins": n_proteins,
        "within_cv": within_cv,
        "fc_rmse": fc_rmse,
        "pct_good_cv": pct_good_cv,
    }


# =============================================================================
# Main Benchmark
# =============================================================================

def run_benchmark(technology: str) -> list:
    """Run preprocessing filter benchmark for a technology."""
    results = []

    print("\n  Loading peptide data...")
    try:
        df_peptides = load_peptide_data(technology)
        n_peptides_orig = len(df_peptides)
        n_proteins_orig = df_peptides["ProteinName"].nunique()
        print(f"  Original: {n_peptides_orig:,} peptide rows, {n_proteins_orig:,} proteins")
    except FileNotFoundError as e:
        print(f"  Skipping {technology}: {e}")
        return results

    # Test all combinations
    total = len(PEPTIDE_LENGTH_CONFIGS) * len(MISSED_CLEAVAGE_CONFIGS) * len(MIN_PEPTIDE_CONFIGS)
    count = 0

    for len_cfg, mc_cfg, pep_cfg in product(
        PEPTIDE_LENGTH_CONFIGS, MISSED_CLEAVAGE_CONFIGS, MIN_PEPTIDE_CONFIGS
    ):
        count += 1
        combo_name = f"{len_cfg['name']}_{mc_cfg['name']}_{pep_cfg['name']}"
        print(f"  [{count}/{total}] {combo_name}...", end=" ")

        try:
            # Apply filters
            df_filtered = apply_filters(
                df_peptides,
                min_length=len_cfg["min_length"],
                max_missed_cleavages=mc_cfg["max_missed_cleavages"],
                min_unique_peptides=pep_cfg["min_unique_peptides"]
            )

            n_peptides = len(df_filtered)
            n_proteins = df_filtered["ProteinName"].nunique()

            if n_peptides == 0:
                print("No peptides left, skipping")
                continue

            # Quantify
            df_quant = quantify_sum(df_filtered)
            df_quant = df_quant.replace(0, np.nan).dropna(thresh=3)

            # Normalize (global median)
            global_median = df_quant.median().median()
            sample_medians = df_quant.median()
            factors = global_median / sample_medians
            df_norm = df_quant * factors

            # Compute metrics
            metrics = compute_metrics(df_norm)

            result = {
                "technology": technology,
                "min_length": len_cfg["min_length"],
                "max_missed_cleavages": mc_cfg["max_missed_cleavages"],
                "min_unique_peptides": pep_cfg["min_unique_peptides"],
                "filter_combo": combo_name,
                "n_peptides": n_peptides,
                "n_peptides_orig": n_peptides_orig,
                "pct_peptides_kept": n_peptides / n_peptides_orig * 100,
                **metrics
            }
            results.append(result)

            print(f"pep={n_peptides:,}, prot={metrics['n_proteins']}, CV={metrics['within_cv']:.4f}")

        except Exception as e:
            print(f"ERROR: {e}")
            import traceback
            traceback.print_exc()
            continue

    return results


def plot_filter_heatmap(results_df: pd.DataFrame, output_path: Path):
    """Plot heatmap of filter combinations."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    for tech_idx, tech in enumerate(["tmt", "lfq"]):
        tech_df = results_df[results_df["technology"] == tech]

        if len(tech_df) == 0:
            continue

        # Create pivot table for length × missed cleavages (for min_pep=2)
        for ax_idx, (metric, title, cmap) in enumerate([
            ("within_cv", "Within-Condition CV", "RdYlGn_r"),
            ("n_proteins", "Proteins Quantified", "RdYlGn"),
        ]):
            ax = axes[tech_idx, ax_idx]

            # Filter for default min_pep=2
            subset = tech_df[tech_df["min_unique_peptides"] == 2].copy()

            if len(subset) == 0:
                continue

            # Pivot: length × missed cleavages
            pivot = subset.pivot(
                index="min_length",
                columns="max_missed_cleavages",
                values=metric
            )

            sns.heatmap(pivot, annot=True, fmt=".3f" if metric == "within_cv" else ".0f",
                       cmap=cmap, ax=ax)
            ax.set_title(f"{tech.upper()}: {title} (min_pep=2)")
            ax.set_xlabel("Max Missed Cleavages")
            ax.set_ylabel("Min Peptide Length")

    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()


def plot_coverage_quality_tradeoff(results_df: pd.DataFrame, output_path: Path):
    """Plot coverage vs quality trade-off."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax, tech in zip(axes, ["tmt", "lfq"]):
        tech_df = results_df[results_df["technology"] == tech].copy()

        if len(tech_df) == 0:
            continue

        # Color by stringency (number of strict filters)
        tech_df["stringency"] = (
            (tech_df["min_length"] - 7) +
            (2 - tech_df["max_missed_cleavages"]) +
            (tech_df["min_unique_peptides"] - 1)
        )

        scatter = ax.scatter(
            tech_df["n_proteins"],
            tech_df["within_cv"],
            c=tech_df["stringency"],
            cmap="coolwarm",
            s=100,
            alpha=0.7,
            edgecolor="black"
        )

        # Label extremes
        best_cv = tech_df.loc[tech_df["within_cv"].idxmin()]
        best_n = tech_df.loc[tech_df["n_proteins"].idxmax()]

        ax.annotate(f"Best CV\n{best_cv['filter_combo']}",
                   (best_cv["n_proteins"], best_cv["within_cv"]),
                   xytext=(10, 10), textcoords="offset points", fontsize=8)
        ax.annotate(f"Most proteins\n{best_n['filter_combo']}",
                   (best_n["n_proteins"], best_n["within_cv"]),
                   xytext=(10, -20), textcoords="offset points", fontsize=8)

        ax.set_xlabel("Proteins Quantified")
        ax.set_ylabel("Within-Condition CV")
        ax.set_title(f"{tech.upper()}: Coverage vs Quality")
        ax.grid(True, alpha=0.3)

        plt.colorbar(scatter, ax=ax, label="Filter Stringency")

    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()


def plot_individual_filter_impact(results_df: pd.DataFrame, output_path: Path):
    """Plot impact of individual filters."""
    fig, axes = plt.subplots(3, 2, figsize=(14, 12))

    filters = [
        ("min_length", [7, 8, 9], "Peptide Length"),
        ("max_missed_cleavages", [2, 1, 0], "Max Missed Cleavages"),
        ("min_unique_peptides", [1, 2, 3], "Min Unique Peptides"),
    ]

    for row, (filter_col, values, title) in enumerate(filters):
        for col, metric in enumerate(["within_cv", "n_proteins"]):
            ax = axes[row, col]

            for tech in ["tmt", "lfq"]:
                tech_df = results_df[results_df["technology"] == tech]

                if len(tech_df) == 0:
                    continue

                # Get mean metric for each filter value
                # (averaging over other filter combinations)
                means = tech_df.groupby(filter_col)[metric].mean()

                ax.plot(means.index, means.values, "o-",
                       label=tech.upper(), linewidth=2, markersize=8)

            ax.set_xlabel(title)
            ax.set_ylabel(metric)
            ax.set_title(f"{title}: {metric}")
            ax.legend()
            ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()


def main():
    """Run Phase 4: Preprocessing Filters Benchmark."""
    print("=" * 60)
    print("Phase 4: Preprocessing Filters Benchmark")
    print("=" * 60)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    all_results = []

    for technology in ["tmt", "lfq"]:
        print(f"\n{'='*60}")
        print(f"{technology.upper()}")
        print("=" * 60)

        results = run_benchmark(technology)
        all_results.extend(results)

    # Save results
    if all_results:
        results_df = pd.DataFrame(all_results)
        results_df.to_csv(RESULTS_DIR / "preprocessing_filters.csv", index=False)
        print(f"\nResults saved to: {RESULTS_DIR / 'preprocessing_filters.csv'}")

        # Generate plots
        print("\nGenerating plots...")

        plot_filter_heatmap(
            results_df,
            FIGURES_DIR / f"filter_heatmap.{FIGURE_FORMAT}"
        )

        plot_coverage_quality_tradeoff(
            results_df,
            FIGURES_DIR / f"filter_coverage_quality.{FIGURE_FORMAT}"
        )

        plot_individual_filter_impact(
            results_df,
            FIGURES_DIR / f"filter_individual_impact.{FIGURE_FORMAT}"
        )

        # Summary
        print("\n" + "=" * 60)
        print("SUMMARY: Filter Impact")
        print("=" * 60)

        for tech in ["tmt", "lfq"]:
            print(f"\n{tech.upper()}:")
            tech_df = results_df[results_df["technology"] == tech]

            if len(tech_df) == 0:
                continue

            # Best by CV
            best_cv = tech_df.loc[tech_df["within_cv"].idxmin()]
            print(f"  Best CV: {best_cv['filter_combo']} "
                  f"(CV={best_cv['within_cv']:.4f}, n={best_cv['n_proteins']})")

            # Best coverage
            best_n = tech_df.loc[tech_df["n_proteins"].idxmax()]
            print(f"  Most proteins: {best_n['filter_combo']} "
                  f"(CV={best_n['within_cv']:.4f}, n={best_n['n_proteins']})")

            # Best RMSE
            valid_rmse = tech_df[tech_df["fc_rmse"].notna()]
            if len(valid_rmse) > 0:
                best_rmse = valid_rmse.loc[valid_rmse["fc_rmse"].idxmin()]
                print(f"  Best RMSE: {best_rmse['filter_combo']} "
                      f"(RMSE={best_rmse['fc_rmse']:.3f})")

            # Default comparison (len≥7, mc≤2, pep≥2)
            default = tech_df[
                (tech_df["min_length"] == 7) &
                (tech_df["max_missed_cleavages"] == 2) &
                (tech_df["min_unique_peptides"] == 2)
            ]
            if len(default) > 0:
                d = default.iloc[0]
                print("\n  Default (len≥7, mc≤2, pep≥2):")
                print(f"    CV={d['within_cv']:.4f}, n={d['n_proteins']}, "
                      f"RMSE={d['fc_rmse']:.3f}")

            # Individual filter impact
            print("\n  Filter impact (mean change from baseline):")
            baseline = tech_df[
                (tech_df["min_length"] == 7) &
                (tech_df["max_missed_cleavages"] == 2) &
                (tech_df["min_unique_peptides"] == 1)
            ]
            if len(baseline) > 0:
                base_cv = baseline.iloc[0]["within_cv"]
                base_n = baseline.iloc[0]["n_proteins"]

                # Length impact
                for l in [8, 9]:
                    subset = tech_df[(tech_df["min_length"] == l) &
                                    (tech_df["max_missed_cleavages"] == 2) &
                                    (tech_df["min_unique_peptides"] == 1)]
                    if len(subset) > 0:
                        cv_change = (subset.iloc[0]["within_cv"] - base_cv) / base_cv * 100
                        n_change = subset.iloc[0]["n_proteins"] - base_n
                        print(f"    len≥{l}: CV {cv_change:+.1f}%, proteins {n_change:+d}")


if __name__ == "__main__":
    main()
