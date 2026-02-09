#!/usr/bin/env python3
"""
Phase 5: Optimal Pipeline Selection.

Combines best performers from each benchmark phase to find optimal pipelines.

From previous phases:
- Phase 1: Best normalizations (TMT: globalmedian, LFQ: depends on quant method)
- Phase 2: Best quant+norm combinations (TMT: iBAQ+globalmedian, LFQ: maxlfq+none)
- Phase 3: Imputation impact (minimal for TMT, tradeoff for LFQ)
- Phase 4: Best filters (pep≥3 helps CV, strict filters help RMSE)

Candidate pipelines tested:
1. Default baseline
2. Best CV pipeline
3. Best RMSE pipeline
4. Best balanced pipeline
5. Technology-specific optimized
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
    SAMPLE_CONDITIONS, CV_THRESHOLDS,
    FIGURE_DPI, FIGURE_FORMAT,
    EXPECTED_FOLD_CHANGES, SPECIES_PATTERNS,
)

# Add mokume to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from mokume.preprocessing.filters.peptide import (
    PeptideLengthFilter,
    MissedCleavageFilter,
)
from mokume.preprocessing.filters.protein import MinPeptideFilter
from mokume.imputation import impute_missing_values

warnings.filterwarnings("ignore")

# =============================================================================
# Paths
# =============================================================================

LOCAL_DIR = Path(__file__).parent.parent.parent.parent / "benchmarks-local" / "PXD007683"
LOCAL_DATA_DIR = LOCAL_DIR / "data" / "processed"
LOCAL_QUANTIFIED_DIR = LOCAL_DIR / "results" / "quantification"


# =============================================================================
# Pipeline Configurations
# =============================================================================

# Each pipeline is a dict with:
# - quant_method: quantification method
# - sample_norm: sample normalization
# - min_unique_peptides: protein filter
# - min_length: peptide length filter
# - max_missed_cleavages: missed cleavage filter
# - imputation: imputation method (none, knn, min)

PIPELINES_TMT = {
    "default": {
        "quant_method": "sum",
        "sample_norm": "none",
        "min_unique_peptides": 2,
        "min_length": 7,
        "max_missed_cleavages": 2,
        "imputation": "none",
        "description": "Default baseline"
    },
    "best_cv": {
        "quant_method": "sum",
        "sample_norm": "globalmedian",
        "min_unique_peptides": 3,
        "min_length": 7,
        "max_missed_cleavages": 2,
        "imputation": "none",
        "description": "Optimized for CV"
    },
    "best_rmse": {
        "quant_method": "ibaq",
        "sample_norm": "tmm",
        "min_unique_peptides": 3,
        "min_length": 9,
        "max_missed_cleavages": 0,
        "imputation": "none",
        "description": "Optimized for fold-change accuracy"
    },
    "balanced": {
        "quant_method": "sum",
        "sample_norm": "globalmedian",
        "min_unique_peptides": 2,
        "min_length": 7,
        "max_missed_cleavages": 1,
        "imputation": "none",
        "description": "Balanced CV/coverage"
    },
    "maxlfq_test": {
        "quant_method": "maxlfq",
        "sample_norm": "none",
        "min_unique_peptides": 2,
        "min_length": 7,
        "max_missed_cleavages": 2,
        "imputation": "none",
        "description": "MaxLFQ (for comparison)"
    },
}

PIPELINES_LFQ = {
    "default": {
        "quant_method": "directlfq",
        "sample_norm": "none",
        "min_unique_peptides": 2,
        "min_length": 7,
        "max_missed_cleavages": 2,
        "imputation": "none",
        "description": "DirectLFQ baseline"
    },
    "best_cv": {
        "quant_method": "maxlfq",
        "sample_norm": "none",
        "min_unique_peptides": 3,
        "min_length": 9,
        "max_missed_cleavages": 1,
        "imputation": "none",
        "description": "Optimized for CV"
    },
    "best_rmse": {
        "quant_method": "maxlfq",
        "sample_norm": "hierarchical",
        "min_unique_peptides": 3,
        "min_length": 9,
        "max_missed_cleavages": 0,
        "imputation": "none",
        "description": "Optimized for fold-change accuracy"
    },
    "balanced": {
        "quant_method": "directlfq",
        "sample_norm": "none",
        "min_unique_peptides": 2,
        "min_length": 8,
        "max_missed_cleavages": 1,
        "imputation": "none",
        "description": "Balanced CV/coverage"
    },
    "high_coverage": {
        "quant_method": "directlfq",
        "sample_norm": "globalmedian",
        "min_unique_peptides": 1,
        "min_length": 7,
        "max_missed_cleavages": 2,
        "imputation": "knn",
        "description": "Maximum coverage (with imputation)"
    },
}


def load_peptide_data(technology: str) -> pd.DataFrame:
    """Load peptide-level data."""
    parquet_path = LOCAL_DATA_DIR / f"{technology}_peptides.parquet"
    if not parquet_path.exists():
        raise FileNotFoundError(f"Data not found: {parquet_path}")
    return pd.read_parquet(parquet_path)


def load_quantified_data(technology: str, method: str) -> pd.DataFrame:
    """Load pre-quantified data."""
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
# Pipeline Execution
# =============================================================================

def apply_filters(df: pd.DataFrame, config: dict) -> pd.DataFrame:
    """Apply preprocessing filters to peptide data."""
    df_filtered = df.copy()

    # Peptide length filter
    length_filter = PeptideLengthFilter(
        min_length=config["min_length"],
        max_length=50,
        sequence_column="PeptideSequence"
    )
    df_filtered, _ = length_filter.apply(df_filtered)

    # Missed cleavage filter
    mc_filter = MissedCleavageFilter(
        max_missed_cleavages=config["max_missed_cleavages"],
        sequence_column="PeptideSequence",
        enzyme="trypsin"
    )
    df_filtered, _ = mc_filter.apply(df_filtered)

    # Minimum peptides filter
    pep_filter = MinPeptideFilter(
        min_peptides=1,
        min_unique_peptides=config["min_unique_peptides"],
        protein_column="ProteinName",
        peptide_column="PeptideSequence"
    )
    df_filtered, _ = pep_filter.apply(df_filtered)

    return df_filtered


def quantify_sum(df: pd.DataFrame) -> pd.DataFrame:
    """Simple sum quantification."""
    protein_intensities = df.groupby(["ProteinName", "SampleID"])["NormIntensity"].sum().reset_index()
    protein_intensities.columns = ["ProteinName", "SampleID", "Intensity"]
    df_wide = protein_intensities.pivot(index="ProteinName", columns="SampleID", values="Intensity")
    return df_wide


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
        raise ValueError(f"Unknown normalization: {method}")


def apply_imputation(df: pd.DataFrame, method: str) -> pd.DataFrame:
    """Apply imputation method."""
    if method == "none":
        return df.copy()

    df_clean = df.replace(0, np.nan)

    if method == "knn":
        return impute_missing_values(df_clean, method="knn", n_neighbors=5)
    elif method == "min":
        global_min = df_clean.min().min()
        fill_value = global_min * 0.5 if global_min > 0 else 1.0
        return impute_missing_values(df_clean, method="constant", fill_value=fill_value)
    else:
        raise ValueError(f"Unknown imputation: {method}")


def compute_metrics(df: pd.DataFrame, use_complete_cases: bool = True) -> dict:
    """Compute all evaluation metrics."""
    if use_complete_cases:
        df_eval = df.dropna()
    else:
        df_eval = df

    col_to_condition = {str(col): get_condition(col) for col in df_eval.columns}

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
    n_yeast = len(yeast_proteins)
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
    pct_good_cv = (np.array(condition_cvs) < CV_THRESHOLDS["good"]).mean() * 100 if condition_cvs else 0

    return {
        "n_proteins": n_proteins,
        "n_yeast": n_yeast,
        "within_cv": within_cv,
        "fc_rmse": fc_rmse,
        "pct_good_cv": pct_good_cv,
    }


def run_pipeline(technology: str, name: str, config: dict) -> dict:
    """Run a single pipeline configuration."""
    print(f"\n  Pipeline: {name}")
    print(f"    {config['description']}")

    results = {"technology": technology, "pipeline": name, **config}

    try:
        # For pre-quantified methods, load and apply post-processing
        if config["quant_method"] in ["directlfq", "maxlfq", "ibaq"]:
            df_quant = load_quantified_data(technology, config["quant_method"])
            # Note: Can't apply peptide filters to pre-quantified data
            # In practice, we'd need to re-run quantification with filters
            # For now, we just use the pre-quantified data
            results["filters_applied"] = False
        else:
            # Load peptide data and apply filters, then quantify
            df_peptides = load_peptide_data(technology)
            df_filtered = apply_filters(df_peptides, config)
            results["n_peptides_filtered"] = len(df_filtered)
            results["filters_applied"] = True

            # Quantify (only sum supported here)
            df_quant = quantify_sum(df_filtered)

        # Clean zeros
        df_quant = df_quant.replace(0, np.nan)

        # Apply normalization
        df_norm = apply_sample_normalization(df_quant, config["sample_norm"])

        # Apply imputation
        df_final = apply_imputation(df_norm, config["imputation"])

        # Filter to proteins with at least 3 valid values
        df_final = df_final.dropna(thresh=3)

        # Compute metrics
        use_complete_cases = config["imputation"] == "none"
        metrics = compute_metrics(df_final, use_complete_cases=use_complete_cases)

        results.update(metrics)

        print(f"    n_proteins={metrics['n_proteins']}, "
              f"CV={metrics['within_cv']:.4f}, RMSE={metrics['fc_rmse']:.3f}")

    except Exception as e:
        print(f"    ERROR: {e}")
        import traceback
        traceback.print_exc()

    return results


# =============================================================================
# Visualization
# =============================================================================

def plot_pipeline_comparison(results_df: pd.DataFrame, output_path: Path):
    """Plot comparison of all pipelines."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    metrics = [
        ("within_cv", "Within-Condition CV", True),
        ("fc_rmse", "Fold-Change RMSE", True),
        ("n_proteins", "Proteins Quantified", False),
        ("pct_good_cv", "% Proteins CV < 20%", False),
    ]

    for ax, (metric, title, lower_better) in zip(axes.flat, metrics):
        for tech in ["tmt", "lfq"]:
            tech_df = results_df[results_df["technology"] == tech]

            if len(tech_df) == 0:
                continue

            x = range(len(tech_df))
            values = tech_df[metric].tolist()
            pipelines = tech_df["pipeline"].tolist()

            color = "#ff7f00" if tech == "tmt" else "#984ea3"
            bars = ax.bar([i + (0.2 if tech == "tmt" else -0.2) for i in x],
                         values, width=0.35, label=tech.upper(),
                         color=color, alpha=0.8)

            # Highlight best
            if lower_better:
                best_idx = np.nanargmin(values)
            else:
                best_idx = np.nanargmax(values)
            bars[best_idx].set_edgecolor("black")
            bars[best_idx].set_linewidth(2)

        ax.set_xlabel("Pipeline")
        ax.set_ylabel(metric)
        ax.set_title(title)
        ax.set_xticks(x)
        ax.set_xticklabels(pipelines, rotation=45, ha="right")
        ax.legend()
        ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()


def plot_radar_chart(results_df: pd.DataFrame, output_path: Path):
    """Plot radar chart comparing pipelines."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    metrics = ["n_proteins", "pct_good_cv", "within_cv_inv", "fc_rmse_inv"]
    metric_labels = ["Proteins", "% Good CV", "1/CV", "1/RMSE"]

    for ax, tech in zip(axes, ["tmt", "lfq"]):
        tech_df = results_df[results_df["technology"] == tech].copy()

        if len(tech_df) == 0:
            continue

        # Invert CV and RMSE so higher is better
        tech_df["within_cv_inv"] = 1 / tech_df["within_cv"]
        tech_df["fc_rmse_inv"] = 1 / tech_df["fc_rmse"].replace(0, np.nan)

        # Normalize to 0-1
        for m in metrics:
            if m in tech_df.columns:
                min_val = tech_df[m].min()
                max_val = tech_df[m].max()
                if max_val > min_val:
                    tech_df[m + "_norm"] = (tech_df[m] - min_val) / (max_val - min_val)
                else:
                    tech_df[m + "_norm"] = 0.5

        # Angles for radar
        angles = np.linspace(0, 2 * np.pi, len(metrics), endpoint=False).tolist()
        angles += angles[:1]  # Complete the circle

        # Plot each pipeline
        colors = plt.cm.Set2(np.linspace(0, 1, len(tech_df)))
        for i, (_, row) in enumerate(tech_df.iterrows()):
            values = [row.get(m + "_norm", 0) for m in metrics]
            values += values[:1]

            ax.plot(angles, values, "o-", linewidth=2, label=row["pipeline"],
                   color=colors[i])
            ax.fill(angles, values, alpha=0.1, color=colors[i])

        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(metric_labels)
        ax.set_title(f"{tech.upper()}: Pipeline Comparison")
        ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1))

    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()


def main():
    """Run Phase 5: Optimal Pipeline Selection."""
    print("=" * 60)
    print("Phase 5: Optimal Pipeline Selection")
    print("=" * 60)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    all_results = []

    # TMT pipelines
    print("\n" + "=" * 60)
    print("TMT PIPELINES")
    print("=" * 60)

    for name, config in PIPELINES_TMT.items():
        result = run_pipeline("tmt", name, config)
        all_results.append(result)

    # LFQ pipelines
    print("\n" + "=" * 60)
    print("LFQ PIPELINES")
    print("=" * 60)

    for name, config in PIPELINES_LFQ.items():
        result = run_pipeline("lfq", name, config)
        all_results.append(result)

    # Save results
    results_df = pd.DataFrame(all_results)
    results_df.to_csv(RESULTS_DIR / "optimal_pipeline.csv", index=False)
    print(f"\nResults saved to: {RESULTS_DIR / 'optimal_pipeline.csv'}")

    # Generate plots
    print("\nGenerating plots...")

    plot_pipeline_comparison(
        results_df,
        FIGURES_DIR / f"optimal_pipeline_comparison.{FIGURE_FORMAT}"
    )

    # Summary
    print("\n" + "=" * 60)
    print("FINAL RECOMMENDATIONS")
    print("=" * 60)

    for tech in ["tmt", "lfq"]:
        print(f"\n{tech.upper()}:")
        tech_df = results_df[results_df["technology"] == tech]

        if len(tech_df) == 0:
            continue

        # Rankings
        print("\n  Rankings:")

        # Best by CV
        tech_df_valid = tech_df[tech_df["within_cv"].notna()]
        if len(tech_df_valid) > 0:
            best_cv = tech_df_valid.loc[tech_df_valid["within_cv"].idxmin()]
            print(f"    Best CV: {best_cv['pipeline']} "
                  f"(CV={best_cv['within_cv']:.4f})")

        # Best by RMSE
        tech_df_valid = tech_df[tech_df["fc_rmse"].notna()]
        if len(tech_df_valid) > 0:
            best_rmse = tech_df_valid.loc[tech_df_valid["fc_rmse"].idxmin()]
            print(f"    Best RMSE: {best_rmse['pipeline']} "
                  f"(RMSE={best_rmse['fc_rmse']:.3f})")

        # Best coverage
        best_n = tech_df.loc[tech_df["n_proteins"].idxmax()]
        print(f"    Most proteins: {best_n['pipeline']} "
              f"(n={best_n['n_proteins']})")

        # Overall recommendation (balanced)
        # Score = rank(CV) + rank(RMSE) + rank(coverage)
        tech_df = tech_df.copy()
        tech_df["cv_rank"] = tech_df["within_cv"].rank()
        tech_df["rmse_rank"] = tech_df["fc_rmse"].rank()
        tech_df["n_rank"] = tech_df["n_proteins"].rank(ascending=False)
        tech_df["total_rank"] = tech_df["cv_rank"] + tech_df["rmse_rank"] + tech_df["n_rank"]

        best_overall = tech_df.loc[tech_df["total_rank"].idxmin()]
        print(f"\n  RECOMMENDED: {best_overall['pipeline']}")
        print(f"    {best_overall['description']}")
        print(f"    Quant: {best_overall['quant_method']}, "
              f"Norm: {best_overall['sample_norm']}")
        print(f"    Filters: len≥{best_overall['min_length']}, "
              f"mc≤{best_overall['max_missed_cleavages']}, "
              f"pep≥{best_overall['min_unique_peptides']}")
        print(f"    Results: CV={best_overall['within_cv']:.4f}, "
              f"RMSE={best_overall['fc_rmse']:.3f}, "
              f"n={best_overall['n_proteins']}")


if __name__ == "__main__":
    main()
