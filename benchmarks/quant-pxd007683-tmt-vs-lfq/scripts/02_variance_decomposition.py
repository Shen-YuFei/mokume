#!/usr/bin/env python3
"""
Variance decomposition analysis: Technology vs Condition vs Residual.

This script answers Q2: What proportion of variance is explained by:
- Condition (yeast spike-in level) - this is biological variance we want to PRESERVE
- Residual (unexplained)

Uses PCA and ANOVA-style variance decomposition.
"""

import logging
import sys
from pathlib import Path
import warnings

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt

# Add parent to path for config import
sys.path.insert(0, str(Path(__file__).parent))
from config import (
    RESULTS_DIR, FIGURES_DIR,
    SAMPLE_CONDITIONS, CONDITION_COLORS,
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


def compute_variance_decomposition(df_wide: pd.DataFrame) -> dict:
    """
    Compute variance decomposition using ANOVA-style approach.

    Returns percentage of variance explained by each factor.
    """
    # Log transform for variance stabilization
    log_data = np.log10(df_wide.replace(0, np.nan))

    # Remove rows with too many NaN
    valid_rows = log_data.notna().sum(axis=1) >= log_data.shape[1] * 0.5
    log_data = log_data.loc[valid_rows]

    # Get complete cases for variance decomposition
    complete_data = log_data.dropna()

    if len(complete_data) < 10:
        return {"error": "Not enough complete cases"}

    # Map columns to conditions
    col_conditions = {col: get_condition(col) for col in complete_data.columns}

    # Total variance
    total_var = complete_data.values.var()

    # Compute mean per condition
    condition_means = {}
    for condition in set(col_conditions.values()):
        if condition == "Unknown":
            continue
        cols = [c for c, cond in col_conditions.items() if cond == condition]
        if cols:
            condition_means[condition] = complete_data[cols].mean(axis=1)

    if len(condition_means) < 2:
        return {"error": "Not enough conditions"}

    # Between-condition variance
    grand_mean = complete_data.mean(axis=1)
    n_samples_per_cond = {
        cond: len([c for c, co in col_conditions.items() if co == cond])
        for cond in condition_means
    }
    between_cond_var = sum(
        n_samples_per_cond[cond] * ((mean - grand_mean) ** 2).mean()
        for cond, mean in condition_means.items()
    ) / complete_data.shape[1]

    # Within-condition variance (residual)
    within_cond_var = total_var - between_cond_var

    return {
        "total_variance": total_var,
        "condition_variance": between_cond_var,
        "residual_variance": within_cond_var,
        "pct_condition": between_cond_var / total_var * 100 if total_var > 0 else 0,
        "pct_residual": within_cond_var / total_var * 100 if total_var > 0 else 0,
        "n_proteins": len(complete_data),
    }


def run_pca_analysis(df_wide: pd.DataFrame) -> dict:
    """Run PCA and compute variance explained by condition."""
    # Log transform
    log_data = np.log10(df_wide.replace(0, np.nan))

    # Get complete cases
    complete_data = log_data.dropna()

    if len(complete_data) < 10:
        return {"error": "Not enough complete cases"}

    # Transpose for PCA (samples as rows)
    X = complete_data.T.values

    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # PCA
    n_components = min(5, X.shape[0], X.shape[1])
    pca = PCA(n_components=n_components)
    pca_coords = pca.fit_transform(X_scaled)

    # Create DataFrame with PCA coordinates
    pca_df = pd.DataFrame(
        pca_coords,
        index=complete_data.columns,
        columns=[f"PC{i+1}" for i in range(pca_coords.shape[1])]
    )

    # Add condition info
    pca_df["condition"] = [get_condition(idx) for idx in pca_df.index]
    pca_df["sample_id"] = pca_df.index

    # Compute silhouette score for condition clustering
    valid_mask = pca_df["condition"] != "Unknown"
    valid_df = pca_df[valid_mask]

    silhouette = np.nan
    if len(valid_df["condition"].unique()) >= 2 and len(valid_df) >= 4:
        try:
            labels = pd.factorize(valid_df["condition"])[0]
            coords = valid_df[["PC1", "PC2"]].values
            silhouette = silhouette_score(coords, labels)
        except Exception as exc:
            logging.warning("Silhouette score computation failed: %s", exc)

    return {
        "pca_df": pca_df,
        "explained_variance": pca.explained_variance_ratio_,
        "silhouette_condition": silhouette,
        "n_proteins": len(complete_data),
    }


def plot_pca_by_condition(pca_results: dict, technology: str, output_path: Path):
    """Create PCA plot colored by condition."""
    pca_df = pca_results["pca_df"]
    var_explained = pca_results["explained_variance"]

    fig, ax = plt.subplots(figsize=(8, 6))

    for condition, color in CONDITION_COLORS.items():
        mask = pca_df["condition"] == condition
        if mask.any():
            ax.scatter(
                pca_df.loc[mask, "PC1"],
                pca_df.loc[mask, "PC2"],
                c=color,
                label=condition,
                s=100,
                alpha=0.8,
                edgecolor="white",
            )
            # Add sample labels
            for idx in pca_df[mask].index:
                ax.annotate(
                    str(idx),
                    (pca_df.loc[idx, "PC1"], pca_df.loc[idx, "PC2"]),
                    fontsize=8,
                    alpha=0.7,
                )

    ax.set_xlabel(f"PC1 ({var_explained[0]*100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({var_explained[1]*100:.1f}% variance)")
    ax.set_title(f"PCA: {technology.upper()} - Colored by Condition")
    ax.legend(title="Condition")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()


def plot_combined_pca(tmt_pca: dict, lfq_pca: dict, output_path: Path):
    """Create combined PCA plot with both technologies."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, (pca_results, tech) in zip(axes, [(tmt_pca, "TMT"), (lfq_pca, "LFQ")]):
        pca_df = pca_results["pca_df"]
        var_explained = pca_results["explained_variance"]

        for condition, color in CONDITION_COLORS.items():
            mask = pca_df["condition"] == condition
            if mask.any():
                ax.scatter(
                    pca_df.loc[mask, "PC1"],
                    pca_df.loc[mask, "PC2"],
                    c=color,
                    label=condition,
                    s=100,
                    alpha=0.8,
                    edgecolor="white",
                )

        ax.set_xlabel(f"PC1 ({var_explained[0]*100:.1f}%)")
        ax.set_ylabel(f"PC2 ({var_explained[1]*100:.1f}%)")
        ax.set_title(f"{tech} (Silhouette: {pca_results['silhouette_condition']:.3f})")
        ax.legend(title="Condition")
        ax.grid(True, alpha=0.3)

    plt.suptitle("PCA Analysis: Condition Separation", fontsize=14)
    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI)
    plt.close()


def main():
    """Run variance decomposition analysis."""
    print("=" * 60)
    print("PXD007683 Benchmark: Variance Decomposition")
    print("=" * 60)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    results = []
    pca_results = {}

    for technology in ["tmt", "lfq"]:
        print(f"\n{technology.upper()}:")
        print("-" * 40)

        try:
            # Load pre-quantified data
            print("  Loading pre-quantified data (MaxLFQ)...")
            df_wide = load_quantified_data(technology, method="maxlfq")
            print(f"  Proteins: {len(df_wide)}, Samples: {len(df_wide.columns)}")

            # Variance decomposition
            print("  Computing variance decomposition...")
            var_results = compute_variance_decomposition(df_wide)

            if "error" not in var_results:
                print(f"    Condition explains: {var_results['pct_condition']:.1f}% of variance")
                print(f"    Residual: {var_results['pct_residual']:.1f}% of variance")
                results.append({
                    "technology": technology,
                    **var_results,
                })

            # PCA analysis
            print("  Running PCA analysis...")
            pca_result = run_pca_analysis(df_wide)

            if "error" not in pca_result:
                pca_results[technology] = pca_result
                print(f"    PC1 explains: {pca_result['explained_variance'][0]*100:.1f}% of variance")
                print(f"    PC2 explains: {pca_result['explained_variance'][1]*100:.1f}% of variance")
                print(f"    Silhouette score (condition): {pca_result['silhouette_condition']:.3f}")

                # Save individual PCA plot
                plot_pca_by_condition(
                    pca_result,
                    technology,
                    FIGURES_DIR / f"pca_{technology}.{FIGURE_FORMAT}"
                )

        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
            continue

    # Save results
    if results:
        results_df = pd.DataFrame(results)
        results_df.to_csv(RESULTS_DIR / "variance_decomposition.csv", index=False)
        print(f"\nResults saved to: {RESULTS_DIR / 'variance_decomposition.csv'}")

    # Combined PCA plot
    if len(pca_results) == 2:
        plot_combined_pca(
            pca_results["tmt"],
            pca_results["lfq"],
            FIGURES_DIR / f"pca_combined.{FIGURE_FORMAT}"
        )
        print(f"Combined PCA plot saved to: {FIGURES_DIR / 'pca_combined.png'}")

    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)

    if results:
        for r in results:
            print(f"\n{r['technology'].upper()}:")
            print(f"  Variance explained by condition: {r['pct_condition']:.1f}%")
            print(f"  Residual variance: {r['pct_residual']:.1f}%")
            print(f"  Proteins analyzed: {r['n_proteins']}")


if __name__ == "__main__":
    main()
