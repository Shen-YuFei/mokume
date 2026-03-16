#!/usr/bin/env python3
"""
Phase 1: Normalization Grid Benchmark.

Tests all sample-level and post-quantification normalization methods
with MaxLFQ quantification (the best performer from previous benchmarks).

Combinations tested:
- Sample Normalization: none, globalmedian, conditionmedian, hierarchical, tmm
- Post Normalization: none, quantile, median_center

Total: 5 × 3 = 15 combinations per technology = 30 total
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
    SAMPLE_CONDITIONS, CONDITION_COLORS,
    CV_THRESHOLDS, LOCAL_QUANTIFIED_DIR,
    FIGURE_DPI, FIGURE_FORMAT,
    EXPECTED_FOLD_CHANGES, SPECIES_PATTERNS,
)

# Add mokume to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

warnings.filterwarnings("ignore")

# =============================================================================
# Normalization Methods to Test
# =============================================================================

SAMPLE_NORMALIZATIONS = ["none", "globalmedian", "conditionmedian", "hierarchical", "tmm"]
POST_NORMALIZATIONS = ["none", "quantile", "median_center"]


def get_default_quant_method(technology: str) -> str:
    """
    Get the appropriate default quantification method for each technology.

    - TMT: Use 'sum' (reporter ion intensities, all samples in same run)
    - LFQ: Use 'directlfq' (MS1 intensities, run-to-run variation)

    MaxLFQ is NOT appropriate for TMT because it uses peptide ratio alignment
    which is designed for LFQ's run-to-run variation.
    """
    if technology.lower() == "tmt":
        return "sum"  # Simple aggregation for reporter ions
    else:
        return "directlfq"  # DirectLFQ for label-free


def load_quantified_data(technology: str, method: str = None) -> pd.DataFrame:
    """Load pre-quantified data from benchmarks-local."""
    if method is None:
        method = get_default_quant_method(technology)

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
        # Normalize each sample to global median
        global_median = df.median().median()
        sample_medians = df.median()
        factors = global_median / sample_medians
        return df * factors

    elif method == "conditionmedian":
        # Normalize within each condition
        result = df.copy()
        col_to_condition = {str(col): get_condition(col) for col in df.columns}

        for condition in set(col_to_condition.values()):
            if condition == "Unknown":
                continue
            cols = [c for c in df.columns if col_to_condition.get(str(c)) == condition]
            if len(cols) > 1:
                condition_median = df[cols].median().median()
                for col in cols:
                    sample_median = df[col].median()
                    if sample_median > 0:
                        result[col] = df[col] * (condition_median / sample_median)
        return result

    elif method == "hierarchical":
        # Use mokume's hierarchical normalizer
        from mokume.normalization.hierarchical import HierarchicalSampleNormalizer

        # Work in log2 space
        df_log = np.log2(df.replace(0, np.nan))
        normalizer = HierarchicalSampleNormalizer(min_overlap=10)
        df_normalized = normalizer.fit_transform(df_log)
        # Convert back to linear space
        return 2 ** df_normalized

    elif method == "tmm":
        # Use mokume's TMM normalizer
        from mokume.normalization.tmm import TMMNormalizer

        df_clean = df.replace(0, np.nan)
        tmm = TMMNormalizer()
        return tmm.fit_transform(df_clean)

    else:
        raise ValueError(f"Unknown sample normalization: {method}")


def apply_post_normalization(df: pd.DataFrame, method: str) -> pd.DataFrame:
    """Apply post-quantification normalization."""
    if method == "none":
        return df.copy()

    elif method == "quantile":
        # Quantile normalization
        from scipy.stats import rankdata

        result = df.copy()
        # Get ranks for each column
        ranks = df.apply(lambda x: rankdata(x, method='average', nan_policy='omit'), axis=0)
        # Get sorted values (ignoring NaN)
        sorted_means = df.apply(lambda x: np.sort(x.dropna().values), axis=0)

        # Use mean of sorted values at each rank
        n_rows = len(df)
        for col in df.columns:
            col_ranks = ranks[col].values
            col_sorted = np.sort(df[col].dropna().values)

            # Interpolate for each rank
            for i, rank in enumerate(col_ranks):
                if not np.isnan(rank) and not np.isnan(df[col].iloc[i]):
                    # Map rank to quantile value
                    idx = int(rank) - 1
                    if idx < len(col_sorted):
                        result.iloc[i, result.columns.get_loc(col)] = col_sorted[idx]

        return result

    elif method == "median_center":
        # Center each sample by subtracting median (in log space)
        df_log = np.log2(df.replace(0, np.nan))
        medians = df_log.median()
        global_median = medians.median()
        df_centered = df_log - medians + global_median
        return 2 ** df_centered

    else:
        raise ValueError(f"Unknown post normalization: {method}")


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

    # Overall CV
    overall_cv = (df.std(axis=1) / df.mean(axis=1)).replace([np.inf, -np.inf], np.nan).median()

    # Q2: Variance explained by condition (using ANOVA-like approach)
    # Compute between-group and within-group variance
    condition_means = {}
    for condition in set(col_to_condition.values()):
        if condition == "Unknown":
            continue
        cols = [c for c in df.columns if col_to_condition.get(str(c)) == condition]
        if cols:
            condition_means[condition] = df[cols].mean(axis=1)

    if len(condition_means) > 1:
        # Between-group variance
        grand_mean = df.mean(axis=1)
        between_var = sum(
            len([c for c in df.columns if col_to_condition.get(str(c)) == cond]) *
            ((condition_means[cond] - grand_mean) ** 2).mean()
            for cond in condition_means
        )
        # Total variance
        total_var = df.var(axis=1).mean() * (len(df.columns) - 1)
        var_explained = (between_var / total_var * 100) if total_var > 0 else 0
    else:
        var_explained = 0

    # Q3: Fold-change accuracy for yeast proteins
    yeast_proteins = [p for p in df.index if get_species(p) == "yeast"]
    human_proteins = [p for p in df.index if get_species(p) == "human"]

    fc_rmse = np.nan
    fc_bias = np.nan

    if yeast_proteins:
        yeast_df = df.loc[yeast_proteins]

        # Compute mean per condition
        condition_intensities = {}
        for condition in ["QY_10pct", "QY_5pct", "QY_3pct"]:
            cols = [c for c in yeast_df.columns if col_to_condition.get(str(c)) == condition]
            if cols:
                condition_intensities[condition] = yeast_df[cols].mean(axis=1)

        if len(condition_intensities) >= 2:
            # Compute fold-changes
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
                fc_bias = np.mean(observed_fcs - expected_fcs)

    # Proteins quantified
    n_proteins = len(df)
    n_yeast = len(yeast_proteins)
    n_human = len(human_proteins)

    # Quality: % proteins with CV < 20%
    pct_good_cv = (np.array(condition_cvs) < CV_THRESHOLDS["good"]).mean() * 100 if condition_cvs else 0

    return {
        "n_proteins": n_proteins,
        "n_yeast": n_yeast,
        "n_human": n_human,
        "within_cv": within_cv,
        "overall_cv": overall_cv,
        "var_explained": var_explained,
        "fc_rmse": fc_rmse,
        "fc_bias": fc_bias,
        "pct_good_cv": pct_good_cv,
    }


# =============================================================================
# Main Benchmark
# =============================================================================

def run_benchmark(technology: str) -> list:
    """Run normalization grid for a technology."""
    results = []

    # Get appropriate quantification method for this technology
    quant_method = get_default_quant_method(technology)
    print(f"\n  Using quantification method: {quant_method}")
    print(f"  (TMT uses 'sum' for reporter ions, LFQ uses 'directlfq' for MS1)")

    # Load base data
    print(f"  Loading {quant_method} data...")
    try:
        df_base = load_quantified_data(technology, method=quant_method)
        print(f"  Loaded {len(df_base)} proteins, {len(df_base.columns)} samples")
    except FileNotFoundError as e:
        print(f"  ERROR: {e}")
        return results

    # Clean data
    df_base = df_base.replace(0, np.nan)
    df_base = df_base.dropna(thresh=3)  # At least 3 non-NA values
    print(f"  After filtering: {len(df_base)} proteins")

    # Test all combinations
    total = len(SAMPLE_NORMALIZATIONS) * len(POST_NORMALIZATIONS)
    count = 0

    for sample_norm, post_norm in product(SAMPLE_NORMALIZATIONS, POST_NORMALIZATIONS):
        count += 1
        combo_name = f"{sample_norm}+{post_norm}"
        print(f"  [{count}/{total}] {combo_name}...", end=" ")

        try:
            # Apply normalizations
            df_norm = apply_sample_normalization(df_base, sample_norm)
            df_norm = apply_post_normalization(df_norm, post_norm)

            # Compute metrics
            metrics = compute_metrics(df_norm)

            result = {
                "technology": technology,
                "quant_method": quant_method,
                "sample_norm": sample_norm,
                "post_norm": post_norm,
                **metrics
            }
            results.append(result)

            print(f"CV={metrics['within_cv']:.4f}, RMSE={metrics['fc_rmse']:.3f}")

        except Exception as e:
            print(f"ERROR: {e}")
            continue

    return results


def plot_heatmap(results_df: pd.DataFrame, metric: str, output_path: Path):
    """Plot heatmap of metric across normalizations."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, tech in zip(axes, ["tmt", "lfq"]):
        tech_df = results_df[results_df["technology"] == tech]

        if len(tech_df) == 0:
            continue

        # Pivot for heatmap
        pivot = tech_df.pivot(
            index="sample_norm",
            columns="post_norm",
            values=metric
        )

        # Determine color scheme (lower is better for CV/RMSE, higher for var_explained)
        if metric in ["within_cv", "overall_cv", "fc_rmse"]:
            cmap = "RdYlGn_r"  # Red = bad (high), Green = good (low)
            fmt = ".3f"
        else:
            cmap = "RdYlGn"  # Green = good (high)
            fmt = ".1f"

        sns.heatmap(pivot, annot=True, fmt=fmt, cmap=cmap, ax=ax,
                    cbar_kws={"label": metric})
        ax.set_title(f"{tech.upper()}: {metric}")
        ax.set_xlabel("Post Normalization")
        ax.set_ylabel("Sample Normalization")

    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()


def plot_comparison_bars(results_df: pd.DataFrame, output_path: Path):
    """Plot bar comparison of best methods."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    metrics = ["within_cv", "fc_rmse", "var_explained", "pct_good_cv"]
    titles = ["Within-Condition CV (lower=better)",
              "Fold-Change RMSE (lower=better)",
              "Variance Explained by Condition % (higher=better)",
              "% Proteins with CV < 20% (higher=better)"]

    for ax, metric, title in zip(axes.flat, metrics, titles):
        # Get best and worst for each technology
        for i, tech in enumerate(["tmt", "lfq"]):
            tech_df = results_df[results_df["technology"] == tech].copy()
            tech_df["combo"] = tech_df["sample_norm"] + "+" + tech_df["post_norm"]

            # Sort by metric
            ascending = metric in ["within_cv", "fc_rmse"]
            tech_df = tech_df.sort_values(metric, ascending=ascending)

            # Plot top 5
            top5 = tech_df.head(5)
            x = np.arange(5) + i * 6
            colors = ["#2ecc71" if j == 0 else "#3498db" for j in range(5)]

            ax.barh(x, top5[metric], color=colors, alpha=0.8)
            ax.set_yticks(x)
            ax.set_yticklabels([f"{tech.upper()}: {c}" for c in top5["combo"]])

        ax.set_xlabel(metric)
        ax.set_title(title)
        ax.grid(True, alpha=0.3, axis="x")

    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()


def main():
    """Run Phase 1: Normalization Grid Benchmark."""
    print("=" * 60)
    print("Phase 1: Normalization Grid Benchmark")
    print("=" * 60)
    print(f"\nSample normalizations: {SAMPLE_NORMALIZATIONS}")
    print(f"Post normalizations: {POST_NORMALIZATIONS}")
    print(f"Total combinations per technology: {len(SAMPLE_NORMALIZATIONS) * len(POST_NORMALIZATIONS)}")

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
        results_df.to_csv(RESULTS_DIR / "normalization_grid.csv", index=False)
        print(f"\nResults saved to: {RESULTS_DIR / 'normalization_grid.csv'}")

        # Generate plots
        print("\nGenerating plots...")

        for metric in ["within_cv", "fc_rmse", "var_explained", "pct_good_cv"]:
            plot_heatmap(
                results_df, metric,
                FIGURES_DIR / f"norm_grid_heatmap_{metric}.{FIGURE_FORMAT}"
            )

        plot_comparison_bars(
            results_df,
            FIGURES_DIR / f"norm_grid_comparison.{FIGURE_FORMAT}"
        )

        # Summary
        print("\n" + "=" * 60)
        print("SUMMARY: Best Normalization Combinations")
        print("=" * 60)

        for tech in ["tmt", "lfq"]:
            print(f"\n{tech.upper()}:")
            tech_df = results_df[results_df["technology"] == tech]

            # Best by CV
            best_cv = tech_df.loc[tech_df["within_cv"].idxmin()]
            print(f"  Best CV: {best_cv['sample_norm']}+{best_cv['post_norm']} "
                  f"(CV={best_cv['within_cv']:.4f})")

            # Best by RMSE
            best_rmse = tech_df.loc[tech_df["fc_rmse"].idxmin()]
            print(f"  Best RMSE: {best_rmse['sample_norm']}+{best_rmse['post_norm']} "
                  f"(RMSE={best_rmse['fc_rmse']:.3f})")

            # Best by variance explained
            best_var = tech_df.loc[tech_df["var_explained"].idxmax()]
            print(f"  Best Var: {best_var['sample_norm']}+{best_var['post_norm']} "
                  f"(Var={best_var['var_explained']:.1f}%)")

            # Baseline (no normalization)
            baseline = tech_df[(tech_df["sample_norm"] == "none") &
                              (tech_df["post_norm"] == "none")]
            if len(baseline) > 0:
                baseline = baseline.iloc[0]
                print(f"  Baseline: CV={baseline['within_cv']:.4f}, "
                      f"RMSE={baseline['fc_rmse']:.3f}")


if __name__ == "__main__":
    main()
