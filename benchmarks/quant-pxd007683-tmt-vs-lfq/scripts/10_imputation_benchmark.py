#!/usr/bin/env python3
"""
Phase 3: Imputation Method Benchmark.

Tests imputation methods for handling missing values.
This is mainly relevant for LFQ data (~17% missing) vs TMT (~0% missing).

Imputation methods tested:
- none: Complete cases only (no imputation)
- knn: K-Nearest Neighbors imputation (good for MAR - Missing At Random)
- min: Minimum value imputation (good for MNAR - Missing Not At Random/below LOD)
- median: Median imputation (simple baseline)
- mean: Mean imputation (simple baseline)

Note: Imputation is expected to have minimal impact on TMT due to low missingness.
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
    SAMPLE_CONDITIONS,
    CV_THRESHOLDS, LOCAL_QUANTIFIED_DIR,
    FIGURE_DPI, FIGURE_FORMAT,
    EXPECTED_FOLD_CHANGES, SPECIES_PATTERNS,
)

# Add mokume to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))
from mokume.imputation import impute_missing_values

warnings.filterwarnings("ignore")

# =============================================================================
# Imputation Methods to Test
# =============================================================================

IMPUTATION_METHODS = ["none", "knn", "min", "median", "mean"]

# Default quantification method per technology (from Phase 2 findings)
DEFAULT_QUANT = {
    "tmt": "sum",       # Simple aggregation for TMT reporter ions
    "lfq": "directlfq"  # DirectLFQ for label-free
}


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
# Imputation Functions
# =============================================================================

def apply_imputation(df: pd.DataFrame, method: str) -> pd.DataFrame:
    """Apply imputation method to handle missing values."""
    if method == "none":
        # No imputation - return data as-is (complete cases will be used in metrics)
        return df.copy()

    elif method == "knn":
        # K-Nearest Neighbors imputation
        df_clean = df.replace(0, np.nan)
        return impute_missing_values(df_clean, method="knn", n_neighbors=5)

    elif method == "min":
        # Minimum value imputation (for MNAR - below limit of detection)
        df_clean = df.replace(0, np.nan)
        # Use global minimum or per-sample minimum
        global_min = df_clean.min().min()
        # Shift down slightly to represent "below detection"
        fill_value = global_min * 0.5 if global_min > 0 else 1.0
        return impute_missing_values(df_clean, method="constant", fill_value=fill_value)

    elif method == "median":
        # Median imputation
        df_clean = df.replace(0, np.nan)
        return impute_missing_values(df_clean, method="median")

    elif method == "mean":
        # Mean imputation
        df_clean = df.replace(0, np.nan)
        return impute_missing_values(df_clean, method="mean")

    else:
        raise ValueError(f"Unknown imputation method: {method}")


def apply_sample_normalization(df: pd.DataFrame, method: str = "globalmedian") -> pd.DataFrame:
    """Apply sample-level normalization."""
    if method == "none":
        return df.copy()
    elif method == "globalmedian":
        global_median = df.median().median()
        sample_medians = df.median()
        factors = global_median / sample_medians
        return df * factors
    else:
        raise ValueError(f"Unknown normalization: {method}")


# =============================================================================
# Evaluation Metrics
# =============================================================================

def compute_metrics(df: pd.DataFrame, use_complete_cases: bool = True) -> dict:
    """Compute all evaluation metrics."""
    col_to_condition = {str(col): get_condition(col) for col in df.columns}

    # For "none" imputation, use complete cases only
    if use_complete_cases:
        df_eval = df.dropna()
    else:
        df_eval = df

    # Q1: Within-condition CV
    condition_cvs = []
    for condition in set(col_to_condition.values()):
        if condition == "Unknown":
            continue
        cols = [c for c in df_eval.columns if col_to_condition.get(str(c)) == condition]
        if len(cols) >= 2:
            subset = df_eval[cols]
            cv = subset.std(axis=1) / subset.mean(axis=1)
            cv = cv.replace([np.inf, -np.inf], np.nan).dropna()
            condition_cvs.extend(cv.tolist())

    within_cv = np.median(condition_cvs) if condition_cvs else np.nan

    # Q3: Fold-change accuracy for yeast proteins
    yeast_proteins = [p for p in df_eval.index if get_species(p) == "yeast"]

    fc_rmse = np.nan
    if yeast_proteins:
        yeast_df = df_eval.loc[yeast_proteins]
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
    n_proteins = len(df_eval)
    n_proteins_total = len(df)
    pct_complete = (n_proteins / n_proteins_total * 100) if n_proteins_total > 0 else 0
    pct_good_cv = (np.array(condition_cvs) < CV_THRESHOLDS["good"]).mean() * 100 if condition_cvs else 0

    # Missing value statistics (on original data before imputation)
    pct_missing = df.isna().sum().sum() / (df.shape[0] * df.shape[1]) * 100

    return {
        "n_proteins": n_proteins,
        "n_proteins_total": n_proteins_total,
        "pct_complete": pct_complete,
        "pct_missing": pct_missing,
        "within_cv": within_cv,
        "fc_rmse": fc_rmse,
        "pct_good_cv": pct_good_cv,
    }


# =============================================================================
# Main Benchmark
# =============================================================================

def run_benchmark(technology: str) -> list:
    """Run imputation benchmark for a technology."""
    results = []

    quant_method = DEFAULT_QUANT[technology]
    print(f"\n  Quantification method: {quant_method}")
    print(f"  Imputation methods: {IMPUTATION_METHODS}")

    # Load quantified data
    try:
        df_base = load_quantified_data(technology, quant_method)
        df_base = df_base.replace(0, np.nan)

        # Report missing values
        total_values = df_base.shape[0] * df_base.shape[1]
        missing_values = df_base.isna().sum().sum()
        pct_missing = missing_values / total_values * 100
        print(f"  Missing values: {pct_missing:.2f}% ({missing_values}/{total_values})")

    except FileNotFoundError as e:
        print(f"  Skipping {technology}: {e}")
        return results

    for imp_method in IMPUTATION_METHODS:
        print(f"  Testing {imp_method}...", end=" ")

        try:
            # Apply imputation
            if imp_method == "none":
                df_imputed = df_base.copy()
                use_complete_cases = True
            else:
                df_imputed = apply_imputation(df_base, imp_method)
                use_complete_cases = False

            # Apply standard normalization
            df_norm = apply_sample_normalization(df_imputed, "globalmedian")

            # Compute metrics
            metrics = compute_metrics(df_norm, use_complete_cases=use_complete_cases)

            result = {
                "technology": technology,
                "quant_method": quant_method,
                "imputation": imp_method,
                **metrics
            }
            results.append(result)

            print(f"n={metrics['n_proteins']}, CV={metrics['within_cv']:.4f}, RMSE={metrics['fc_rmse']:.3f}")

        except Exception as e:
            print(f"ERROR: {e}")
            import traceback
            traceback.print_exc()
            continue

    return results


def plot_imputation_comparison(results_df: pd.DataFrame, output_path: Path):
    """Plot imputation method comparison."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    metrics = [
        ("within_cv", "Within-Condition CV", True),
        ("fc_rmse", "Fold-Change RMSE", True),
        ("n_proteins", "Proteins Quantified", False),
        ("pct_good_cv", "% Proteins CV < 20%", False)
    ]

    for ax, (metric, title, lower_better) in zip(axes.flat, metrics):
        for i, tech in enumerate(["tmt", "lfq"]):
            tech_df = results_df[results_df["technology"] == tech]

            if len(tech_df) == 0:
                continue

            methods = tech_df["imputation"].tolist()
            values = tech_df[metric].tolist()

            x = np.arange(len(methods))
            width = 0.35
            offset = (i - 0.5) * width

            ax.bar(x + offset, values, width, label=tech.upper(),
                   color=["#ff7f00", "#984ea3"][i], alpha=0.8)

        ax.set_xlabel("Imputation Method")
        ax.set_ylabel(metric)
        ax.set_title(title)
        ax.set_xticks(x)
        ax.set_xticklabels(methods, rotation=45, ha="right")
        ax.legend()
        ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()


def plot_missing_value_impact(results_df: pd.DataFrame, output_path: Path):
    """Plot the impact of imputation on coverage and quality."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax, tech in zip(axes, ["tmt", "lfq"]):
        tech_df = results_df[results_df["technology"] == tech].copy()

        if len(tech_df) == 0:
            continue

        # Scatter plot: n_proteins vs CV
        colors = plt.cm.Set2(np.linspace(0, 1, len(tech_df)))

        for i, (_, row) in enumerate(tech_df.iterrows()):
            ax.scatter(row["n_proteins"], row["within_cv"],
                      s=200, c=[colors[i]], label=row["imputation"],
                      edgecolor="black", linewidth=1.5)
            ax.annotate(row["imputation"],
                       (row["n_proteins"], row["within_cv"]),
                       xytext=(5, 5), textcoords="offset points",
                       fontsize=9)

        ax.set_xlabel("Proteins Quantified")
        ax.set_ylabel("Within-Condition CV")
        ax.set_title(f"{tech.upper()}: Coverage vs Quality Trade-off")
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()


def main():
    """Run Phase 3: Imputation Method Benchmark."""
    print("=" * 60)
    print("Phase 3: Imputation Method Benchmark")
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
        results_df.to_csv(RESULTS_DIR / "imputation_benchmark.csv", index=False)
        print(f"\nResults saved to: {RESULTS_DIR / 'imputation_benchmark.csv'}")

        # Generate plots
        print("\nGenerating plots...")

        plot_imputation_comparison(
            results_df,
            FIGURES_DIR / f"imputation_comparison.{FIGURE_FORMAT}"
        )

        plot_missing_value_impact(
            results_df,
            FIGURES_DIR / f"imputation_coverage_quality.{FIGURE_FORMAT}"
        )

        # Summary
        print("\n" + "=" * 60)
        print("SUMMARY: Imputation Impact")
        print("=" * 60)

        for tech in ["tmt", "lfq"]:
            print(f"\n{tech.upper()}:")
            tech_df = results_df[results_df["technology"] == tech]

            if len(tech_df) == 0:
                continue

            # Report missing values
            pct_missing = tech_df["pct_missing"].iloc[0]
            print(f"  Missing values: {pct_missing:.2f}%")

            # Best by coverage
            best_coverage = tech_df.loc[tech_df["n_proteins"].idxmax()]
            print(f"  Best coverage: {best_coverage['imputation']} "
                  f"(n={best_coverage['n_proteins']}, CV={best_coverage['within_cv']:.4f})")

            # Best by CV
            best_cv = tech_df.loc[tech_df["within_cv"].idxmin()]
            print(f"  Best CV: {best_cv['imputation']} "
                  f"(n={best_cv['n_proteins']}, CV={best_cv['within_cv']:.4f})")

            # Best by RMSE
            valid_rmse = tech_df[tech_df["fc_rmse"].notna()]
            if len(valid_rmse) > 0:
                best_rmse = valid_rmse.loc[valid_rmse["fc_rmse"].idxmin()]
                print(f"  Best RMSE: {best_rmse['imputation']} "
                      f"(RMSE={best_rmse['fc_rmse']:.3f})")

            # Impact summary
            none_result = tech_df[tech_df["imputation"] == "none"]
            if len(none_result) > 0:
                none_cv = none_result["within_cv"].iloc[0]
                none_n = none_result["n_proteins"].iloc[0]

                print(f"\n  Impact vs no imputation (complete cases, n={none_n}):")
                for _, row in tech_df.iterrows():
                    if row["imputation"] != "none":
                        cv_change = (row["within_cv"] - none_cv) / none_cv * 100
                        n_change = row["n_proteins"] - none_n
                        print(f"    {row['imputation']}: "
                              f"CV {cv_change:+.1f}%, proteins {n_change:+d}")


if __name__ == "__main__":
    main()
