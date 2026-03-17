#!/usr/bin/env python3
"""
Phase 2: Quantification × Normalization Interaction Benchmark.

Tests whether the best normalization method varies by quantification method.

Combinations tested:
- TMT Quantification: sum, top3, top10, maxlfq, ibaq (NO directlfq - designed for LFQ)
- LFQ Quantification: sum, top3, top10, maxlfq, directlfq, ibaq
- Sample Normalization: none, globalmedian, hierarchical, tmm

Note: MaxLFQ is included for comparison but is NOT recommended for TMT.
"""

import sys
from pathlib import Path
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Add parent to path for config import
sys.path.insert(0, str(Path(__file__).parent))
from config import (
    RESULTS_DIR, FIGURES_DIR,
    SAMPLE_CONDITIONS,
    CV_THRESHOLDS, LOCAL_QUANTIFIED_DIR,
    FIGURE_DPI, FIGURE_FORMAT,
    EXPECTED_FOLD_CHANGES, SPECIES_PATTERNS,
)

# Add mokume to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

warnings.filterwarnings("ignore")

# =============================================================================
# Methods to Test
# =============================================================================

# Quantification methods appropriate for each technology
QUANT_METHODS_TMT = ["sum", "top3", "top10", "ibaq"]  # NO directlfq, maxlfq not recommended
QUANT_METHODS_LFQ = ["sum", "top3", "top10", "maxlfq", "directlfq", "ibaq"]

# Sample normalizations to test
SAMPLE_NORMALIZATIONS = ["none", "globalmedian", "hierarchical", "tmm"]


def load_quantified_data(technology: str, method: str) -> pd.DataFrame:
    """Load pre-quantified data from benchmarks-local."""
    parquet_path = LOCAL_QUANTIFIED_DIR / technology / f"{method}.parquet"

    if not parquet_path.exists():
        raise FileNotFoundError(f"Data not found: {parquet_path}")

    df = pd.read_parquet(parquet_path)
    df_wide = df.pivot(index="ProteinName", columns="SampleID", values="Intensity")
    return df_wide


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
# Normalization Functions
# =============================================================================

def apply_sample_normalization(df: pd.DataFrame, method: str) -> pd.DataFrame:
    """Apply sample-level normalization."""
    if method == "none":
        return df.copy()

    elif method == "globalmedian":
        global_median = df.median().median()
        sample_medians = df.median()
        factors = global_median / sample_medians
        return df * factors

    elif method == "hierarchical":
        from mokume.normalization.hierarchical import HierarchicalSampleNormalizer
        df_log = np.log2(df.replace(0, np.nan))
        normalizer = HierarchicalSampleNormalizer(min_overlap=10)
        df_normalized = normalizer.fit_transform(df_log)
        return 2 ** df_normalized

    elif method == "tmm":
        from mokume.normalization.tmm import TMMNormalizer
        df_clean = df.replace(0, np.nan)
        tmm = TMMNormalizer()
        return tmm.fit_transform(df_clean)

    else:
        raise ValueError(f"Unknown sample normalization: {method}")


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
    """Run quantification × normalization grid for a technology."""
    results = []

    # Get appropriate quantification methods
    if technology == "tmt":
        quant_methods = QUANT_METHODS_TMT
    else:
        quant_methods = QUANT_METHODS_LFQ

    print(f"\n  Quantification methods: {quant_methods}")
    print(f"  Normalization methods: {SAMPLE_NORMALIZATIONS}")

    total = len(quant_methods) * len(SAMPLE_NORMALIZATIONS)
    count = 0

    for quant_method in quant_methods:
        # Load quantified data
        try:
            df_base = load_quantified_data(technology, quant_method)
            df_base = df_base.replace(0, np.nan).dropna(thresh=3)
            n_proteins = len(df_base)
        except FileNotFoundError as e:
            print(f"  Skipping {quant_method}: {e}")
            continue

        for sample_norm in SAMPLE_NORMALIZATIONS:
            count += 1
            combo_name = f"{quant_method}+{sample_norm}"
            print(f"  [{count}/{total}] {combo_name}...", end=" ")

            try:
                # Apply normalization
                df_norm = apply_sample_normalization(df_base, sample_norm)

                # Compute metrics
                metrics = compute_metrics(df_norm)

                result = {
                    "technology": technology,
                    "quant_method": quant_method,
                    "sample_norm": sample_norm,
                    **metrics
                }
                results.append(result)

                print(f"CV={metrics['within_cv']:.4f}, RMSE={metrics['fc_rmse']:.3f}")

            except Exception as e:
                print(f"ERROR: {e}")
                continue

    return results


def plot_interaction_heatmap(results_df: pd.DataFrame, metric: str, output_path: Path):
    """Plot heatmap showing quant × norm interaction."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, tech in zip(axes, ["tmt", "lfq"]):
        tech_df = results_df[results_df["technology"] == tech]

        if len(tech_df) == 0:
            continue

        # Pivot for heatmap
        pivot = tech_df.pivot(
            index="quant_method",
            columns="sample_norm",
            values=metric
        )

        # Color scheme
        if metric in ["within_cv", "fc_rmse"]:
            cmap = "RdYlGn_r"
            fmt = ".3f"
        else:
            cmap = "RdYlGn"
            fmt = ".1f"

        sns.heatmap(pivot, annot=True, fmt=fmt, cmap=cmap, ax=ax,
                    cbar_kws={"label": metric})
        ax.set_title(f"{tech.upper()}: {metric}")
        ax.set_xlabel("Sample Normalization")
        ax.set_ylabel("Quantification Method")

    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()


def plot_best_combinations(results_df: pd.DataFrame, output_path: Path):
    """Plot best method combinations."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, tech in zip(axes, ["tmt", "lfq"]):
        tech_df = results_df[results_df["technology"] == tech].copy()

        if len(tech_df) == 0:
            continue

        tech_df["combo"] = tech_df["quant_method"] + " + " + tech_df["sample_norm"]

        # Sort by CV
        tech_df = tech_df.sort_values("within_cv")

        # Plot
        colors = plt.cm.RdYlGn_r(np.linspace(0.2, 0.8, len(tech_df)))
        bars = ax.barh(range(len(tech_df)), tech_df["within_cv"], color=colors)

        ax.set_yticks(range(len(tech_df)))
        ax.set_yticklabels(tech_df["combo"], fontsize=8)
        ax.set_xlabel("Within-Condition CV (lower = better)")
        ax.set_title(f"{tech.upper()}: All Combinations Ranked by CV")
        ax.grid(True, alpha=0.3, axis="x")

        # Highlight best
        bars[0].set_color("#2ecc71")
        bars[0].set_edgecolor("black")
        bars[0].set_linewidth(2)

    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()


def plot_quant_comparison(results_df: pd.DataFrame, output_path: Path):
    """Plot quantification method comparison (best norm for each)."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    metrics = [("within_cv", "Within-Condition CV", True),
               ("fc_rmse", "Fold-Change RMSE", True),
               ("n_proteins", "Proteins Quantified", False),
               ("pct_good_cv", "% Proteins CV < 20%", False)]

    for ax, (metric, title, lower_better) in zip(axes.flat, metrics):
        for i, tech in enumerate(["tmt", "lfq"]):
            tech_df = results_df[results_df["technology"] == tech]

            if len(tech_df) == 0:
                continue

            # Get best normalization for each quant method
            if lower_better:
                best = tech_df.loc[tech_df.groupby("quant_method")[metric].idxmin()]
            else:
                best = tech_df.loc[tech_df.groupby("quant_method")[metric].idxmax()]

            quant_methods = best["quant_method"].tolist()
            values = best[metric].tolist()
            norms = best["sample_norm"].tolist()

            x = np.arange(len(quant_methods))
            width = 0.35
            offset = (i - 0.5) * width

            bars = ax.bar(x + offset, values, width, label=tech.upper(),
                         color=["#ff7f00", "#984ea3"][i], alpha=0.8)

            # Annotate with best norm
            for j, (bar, norm) in enumerate(zip(bars, norms)):
                ax.annotate(norm, xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                           xytext=(0, 3), textcoords="offset points",
                           ha="center", va="bottom", fontsize=7, rotation=45)

        ax.set_xlabel("Quantification Method")
        ax.set_ylabel(metric)
        ax.set_title(title)
        ax.set_xticks(x)
        ax.set_xticklabels(quant_methods, rotation=45, ha="right")
        ax.legend()
        ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()


def main():
    """Run Phase 2: Quantification × Normalization Interaction Benchmark."""
    print("=" * 60)
    print("Phase 2: Quantification × Normalization Interaction")
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
        results_df.to_csv(RESULTS_DIR / "quant_norm_interaction.csv", index=False)
        print(f"\nResults saved to: {RESULTS_DIR / 'quant_norm_interaction.csv'}")

        # Generate plots
        print("\nGenerating plots...")

        for metric in ["within_cv", "fc_rmse"]:
            plot_interaction_heatmap(
                results_df, metric,
                FIGURES_DIR / f"quant_norm_heatmap_{metric}.{FIGURE_FORMAT}"
            )

        plot_best_combinations(
            results_df,
            FIGURES_DIR / f"quant_norm_ranking.{FIGURE_FORMAT}"
        )

        plot_quant_comparison(
            results_df,
            FIGURES_DIR / f"quant_method_comparison.{FIGURE_FORMAT}"
        )

        # Summary
        print("\n" + "=" * 60)
        print("SUMMARY: Best Combinations")
        print("=" * 60)

        for tech in ["tmt", "lfq"]:
            print(f"\n{tech.upper()}:")
            tech_df = results_df[results_df["technology"] == tech]

            if len(tech_df) == 0:
                continue

            # Best by CV
            best_cv = tech_df.loc[tech_df["within_cv"].idxmin()]
            print(f"  Best CV: {best_cv['quant_method']} + {best_cv['sample_norm']} "
                  f"(CV={best_cv['within_cv']:.4f})")

            # Best by RMSE
            best_rmse = tech_df.loc[tech_df["fc_rmse"].idxmin()]
            print(f"  Best RMSE: {best_rmse['quant_method']} + {best_rmse['sample_norm']} "
                  f"(RMSE={best_rmse['fc_rmse']:.3f})")

            # Best per quantification method
            print("\n  Best normalization per quant method (by CV):")
            for quant in tech_df["quant_method"].unique():
                quant_df = tech_df[tech_df["quant_method"] == quant]
                best = quant_df.loc[quant_df["within_cv"].idxmin()]
                print(f"    {quant}: {best['sample_norm']} (CV={best['within_cv']:.4f})")


if __name__ == "__main__":
    main()
