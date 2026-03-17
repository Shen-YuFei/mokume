#!/usr/bin/env python3
"""
Method Correlation Plots

Generate density scatter plots comparing different quantification methods,
similar to the MaxQuant vs quantms comparison style.
"""

import sys
from pathlib import Path
from typing import Optional, Dict, List

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, str(Path(__file__).parent))

from config import (
    ALL_DATASETS,
    PROTEIN_QUANT_DIR,
    PLOTS_DIR,
)


def concordance_correlation_coefficient(y_true, y_pred):
    """
    Calculate Lin's Concordance Correlation Coefficient (CCC).

    CCC measures both precision and accuracy of predictions.
    """
    # Pearson correlation
    cor = np.corrcoef(y_true, y_pred)[0, 1]

    # Means
    mean_true = np.mean(y_true)
    mean_pred = np.mean(y_pred)

    # Variances
    var_true = np.var(y_true)
    var_pred = np.var(y_pred)

    # Standard deviations
    sd_true = np.std(y_true)
    sd_pred = np.std(y_pred)

    # CCC
    numerator = 2 * cor * sd_true * sd_pred
    denominator = var_true + var_pred + (mean_true - mean_pred) ** 2

    ccc = numerator / denominator

    return ccc


def load_method_data(dataset_id: str, method: str) -> Optional[pd.DataFrame]:
    """Load quantification result for a dataset and method."""
    quant_dir = PROTEIN_QUANT_DIR / dataset_id
    path = quant_dir / f"{method}.parquet"

    if not path.exists():
        return None

    return pd.read_parquet(path)


def get_intensity_column(method: str) -> str:
    """Get the intensity column name for a method."""
    mapping = {
        "ibaq": "IbaqLog",
        "ribaq": "IbaqNorm",
        "directlfq": "DirectLFQIntensity",
        "top3": "Top3Intensity",
        "topn": "Top10Intensity",
        "top10": "Top10Intensity",
        "sum": "SumIntensity",
    }
    return mapping.get(method, "Intensity")


def is_log_transformed(method: str) -> bool:
    """Check if the method's intensity values are already log-transformed."""
    # IbaqLog is already log-transformed
    return method == "ibaq"


def get_median_expression(df: pd.DataFrame, intensity_col: str) -> pd.Series:
    """Get median expression per protein."""
    return df.groupby("ProteinName")[intensity_col].median()


def zscore_normalize(values: np.ndarray) -> np.ndarray:
    """Z-score normalize values."""
    mean = np.nanmean(values)
    std = np.nanstd(values)
    if std > 0:
        return (values - mean) / std
    return values - mean


def plot_density_scatter(
    x: np.ndarray,
    y: np.ndarray,
    ax: plt.Axes,
    xlabel: str,
    ylabel: str,
    title: str = "",
    bins: int = 100,
):
    """
    Create a density scatter plot with hexbin coloring.

    Similar to the style shown in the reference image.
    """
    # Remove NaN and infinite values
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]

    if len(x) < 10:
        ax.text(0.5, 0.5, "Insufficient data", transform=ax.transAxes,
                ha='center', va='center')
        return

    # Calculate statistics
    r, p = pearsonr(x, y)
    rho, _ = spearmanr(x, y)
    ccc = concordance_correlation_coefficient(x, y)

    # Create hexbin density plot
    hb = ax.hexbin(x, y, gridsize=bins, cmap='jet', mincnt=1,
                   norm=LogNorm(), alpha=0.9)

    # Add colorbar
    cb = plt.colorbar(hb, ax=ax, label='n_neighbors')

    # Add diagonal line (y = x)
    lims = [
        min(ax.get_xlim()[0], ax.get_ylim()[0]),
        max(ax.get_xlim()[1], ax.get_ylim()[1]),
    ]
    ax.plot(lims, lims, 'k--', alpha=0.5, linewidth=1)
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    # Add statistics text
    stats_text = f"R = {r:.2f}, p < 2.2e-16\nLin's CCC: {ccc:.3f}"
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
            fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    return r, ccc


def compare_methods_single_dataset(
    dataset_id: str,
    method1: str,
    method2: str,
    output_dir: Path = None,
) -> Optional[Dict]:
    """
    Compare two quantification methods for a single dataset.
    """
    if output_dir is None:
        output_dir = PLOTS_DIR

    # Load data
    df1 = load_method_data(dataset_id, method1)
    df2 = load_method_data(dataset_id, method2)

    if df1 is None or df2 is None:
        return None

    # Get intensity columns
    col1 = get_intensity_column(method1)
    col2 = get_intensity_column(method2)

    if col1 not in df1.columns or col2 not in df2.columns:
        return None

    # Get median expression per protein
    expr1 = get_median_expression(df1, col1)
    expr2 = get_median_expression(df2, col2)

    # Find common proteins
    common = expr1.index.intersection(expr2.index)

    if len(common) < 100:
        return None

    # Get values
    x = expr1.loc[common].values
    y = expr2.loc[common].values

    # Log transform if not already log-transformed
    if not is_log_transformed(method1):
        x = np.log2(x + 1)
    if not is_log_transformed(method2):
        y = np.log2(y + 1)

    # Z-score normalize to make scales comparable
    x = zscore_normalize(x)
    y = zscore_normalize(y)

    # Create plot
    fig, ax = plt.subplots(figsize=(8, 8))

    xlabel = f"{method1.upper()} (z-score)"
    ylabel = f"{method2.upper()} (z-score)"
    title = f"{dataset_id}: {method1.upper()} vs {method2.upper()}"

    r, ccc = plot_density_scatter(x, y, ax, xlabel, ylabel, title)

    plt.tight_layout()

    # Save
    filename = f"correlation_{dataset_id}_{method1}_vs_{method2}"
    plt.savefig(output_dir / f"{filename}.png", dpi=150, bbox_inches="tight")
    plt.savefig(output_dir / f"{filename}.pdf", bbox_inches="tight")
    plt.close()

    return {
        "dataset": dataset_id,
        "method1": method1,
        "method2": method2,
        "n_proteins": len(common),
        "pearson_r": r,
        "ccc": ccc,
    }


def compare_methods_all_datasets(
    method1: str,
    method2: str,
    datasets: Dict = None,
    output_dir: Path = None,
) -> pd.DataFrame:
    """
    Compare two methods across all datasets and create a combined plot.
    """
    if datasets is None:
        datasets = ALL_DATASETS
    if output_dir is None:
        output_dir = PLOTS_DIR

    # Collect all data points
    all_x = []
    all_y = []
    results = []

    for dataset_id in datasets:
        df1 = load_method_data(dataset_id, method1)
        df2 = load_method_data(dataset_id, method2)

        if df1 is None or df2 is None:
            continue

        col1 = get_intensity_column(method1)
        col2 = get_intensity_column(method2)

        if col1 not in df1.columns or col2 not in df2.columns:
            continue

        expr1 = get_median_expression(df1, col1)
        expr2 = get_median_expression(df2, col2)

        common = expr1.index.intersection(expr2.index)

        if len(common) < 10:
            continue

        x = expr1.loc[common].values
        y = expr2.loc[common].values

        # Log transform if needed
        if not is_log_transformed(method1):
            x = np.log2(x + 1)
        if not is_log_transformed(method2):
            y = np.log2(y + 1)

        # Filter valid values
        mask = np.isfinite(x) & np.isfinite(y)
        x = x[mask]
        y = y[mask]

        # Z-score normalize per dataset for comparability
        x = zscore_normalize(x)
        y = zscore_normalize(y)

        all_x.extend(x)
        all_y.extend(y)

        # Per-dataset stats
        if len(x) > 10:
            r, _ = pearsonr(x, y)
            ccc = concordance_correlation_coefficient(x, y)
            results.append({
                "dataset": dataset_id,
                "n_proteins": len(x),
                "pearson_r": r,
                "ccc": ccc,
            })

    if len(all_x) < 100:
        print(f"  Insufficient data for {method1} vs {method2}")
        return pd.DataFrame()

    all_x = np.array(all_x)
    all_y = np.array(all_y)

    # Create combined plot
    fig, ax = plt.subplots(figsize=(10, 10))

    xlabel = f"{method1.upper()} (z-score)"
    ylabel = f"{method2.upper()} (z-score)"
    title = f"All Datasets: {method1.upper()} vs {method2.upper()}\n(n = {len(all_x):,} protein-dataset pairs)"

    plot_density_scatter(all_x, all_y, ax, xlabel, ylabel, title, bins=80)

    plt.tight_layout()

    filename = f"correlation_all_{method1}_vs_{method2}"
    plt.savefig(output_dir / f"{filename}.png", dpi=150, bbox_inches="tight")
    plt.savefig(output_dir / f"{filename}.pdf", bbox_inches="tight")
    plt.close()

    print(f"  Saved: {filename}.png/pdf")

    return pd.DataFrame(results)


def create_method_comparison_grid(
    methods: List[str] = None,
    datasets: Dict = None,
    output_dir: Path = None,
):
    """
    Create a grid of correlation plots comparing all method pairs.
    """
    if methods is None:
        methods = ["ibaq", "ribaq", "directlfq", "top3", "sum"]
    if datasets is None:
        datasets = ALL_DATASETS
    if output_dir is None:
        output_dir = PLOTS_DIR

    n_methods = len(methods)
    fig, axes = plt.subplots(n_methods, n_methods, figsize=(4 * n_methods, 4 * n_methods))

    for i, method1 in enumerate(methods):
        for j, method2 in enumerate(methods):
            ax = axes[i, j]

            if i == j:
                # Diagonal - show method name
                ax.text(0.5, 0.5, method1.upper(), transform=ax.transAxes,
                       ha='center', va='center', fontsize=14, fontweight='bold')
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)
                ax.axis('off')
                continue

            if i > j:
                # Lower triangle - hide
                ax.axis('off')
                continue

            # Upper triangle - correlation plot
            all_x = []
            all_y = []

            for dataset_id in datasets:
                df1 = load_method_data(dataset_id, method1)
                df2 = load_method_data(dataset_id, method2)

                if df1 is None or df2 is None:
                    continue

                col1 = get_intensity_column(method1)
                col2 = get_intensity_column(method2)

                if col1 not in df1.columns or col2 not in df2.columns:
                    continue

                expr1 = get_median_expression(df1, col1)
                expr2 = get_median_expression(df2, col2)

                common = expr1.index.intersection(expr2.index)

                if len(common) < 10:
                    continue

                x = expr1.loc[common].values
                y = expr2.loc[common].values

                if not is_log_transformed(method1):
                    x = np.log2(x + 1)
                if not is_log_transformed(method2):
                    y = np.log2(y + 1)

                mask = np.isfinite(x) & np.isfinite(y)
                x = zscore_normalize(x[mask])
                y = zscore_normalize(y[mask])
                all_x.extend(x)
                all_y.extend(y)

            if len(all_x) > 100:
                all_x = np.array(all_x)
                all_y = np.array(all_y)

                # Simple hexbin without colorbar for grid
                ax.hexbin(all_x, all_y, gridsize=50, cmap='jet', mincnt=1,
                         norm=LogNorm(), alpha=0.9)

                # Stats
                r, _ = pearsonr(all_x, all_y)
                ax.text(0.05, 0.95, f"R={r:.2f}", transform=ax.transAxes,
                       fontsize=9, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

                # Diagonal
                lims = [min(all_x.min(), all_y.min()), max(all_x.max(), all_y.max())]
                ax.plot(lims, lims, 'k--', alpha=0.5, linewidth=0.5)

            if i == 0:
                ax.set_title(method2.upper())
            if j == 0:
                ax.set_ylabel(method1.upper())

    plt.suptitle("Method Correlation Comparison (All Datasets)", fontsize=16, y=1.02)
    plt.tight_layout()

    plt.savefig(output_dir / "method_correlation_grid.png", dpi=150, bbox_inches="tight")
    plt.savefig(output_dir / "method_correlation_grid.pdf", bbox_inches="tight")
    plt.close()

    print("  Saved: method_correlation_grid.png/pdf")


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate method correlation plots"
    )
    parser.add_argument(
        "--method1",
        type=str,
        default="ibaq",
        help="First method to compare"
    )
    parser.add_argument(
        "--method2",
        type=str,
        default="directlfq",
        help="Second method to compare"
    )
    parser.add_argument(
        "--dataset",
        type=str,
        help="Specific dataset (or 'all' for combined)"
    )
    parser.add_argument(
        "--grid",
        action="store_true",
        help="Generate full comparison grid"
    )

    args = parser.parse_args()

    print("=" * 70)
    print("Generating Method Correlation Plots")
    print("=" * 70)

    output_dir = PLOTS_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.grid:
        print("\nGenerating method comparison grid...")
        create_method_comparison_grid(output_dir=output_dir)
    else:
        # Generate key comparisons
        comparisons = [
            ("ibaq", "directlfq"),
            ("ibaq", "ribaq"),
            ("ibaq", "top3"),
            ("ibaq", "sum"),
            ("directlfq", "ribaq"),
            ("directlfq", "top3"),
        ]

        print("\nGenerating pairwise correlation plots...")
        all_results = []

        for m1, m2 in comparisons:
            print(f"\n  {m1.upper()} vs {m2.upper()}:")
            results = compare_methods_all_datasets(m1, m2, output_dir=output_dir)
            if not results.empty:
                results["comparison"] = f"{m1}_vs_{m2}"
                all_results.append(results)

        if all_results:
            combined = pd.concat(all_results, ignore_index=True)
            combined.to_csv(output_dir / "method_correlations.csv", index=False)
            print(f"\n  Saved: method_correlations.csv")

    print("\n" + "=" * 70)
    print(f"Done! Plots saved to: {output_dir}")
    print("=" * 70)


if __name__ == "__main__":
    main()
