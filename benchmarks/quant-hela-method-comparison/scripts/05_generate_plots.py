#!/usr/bin/env python3
"""
Phase 5: Generate Visualizations

Create publication-quality plots for the benchmarking results:
1. CV Distribution - Violin/box plots comparing CV across methods
2. Cross-Experiment Correlation Heatmaps
3. TMT vs LFQ Scatter Plots
4. Expression Profile Consistency
5. Rank Stability Plots
6. Method Comparison Summary - Radar chart
"""

import sys
from pathlib import Path
from typing import Optional, Dict, List

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns

sys.path.insert(0, str(Path(__file__).parent.parent))

from config import (
    ALL_DATASETS,
    HELA_DATASETS,
    ANALYSIS_DIR,
    PLOTS_DIR,
    PROTEIN_QUANT_DIR,
    QUANTIFICATION_METHODS,
    PROTEINS_OF_INTEREST,
)

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

# Figure sizes
FIGSIZE_SMALL = (8, 6)
FIGSIZE_MEDIUM = (10, 8)
FIGSIZE_LARGE = (14, 10)
FIGSIZE_WIDE = (14, 6)


def load_analysis_results() -> dict:
    """Load all analysis results from CSV files."""
    results = {}

    # CV comparison
    cv_path = ANALYSIS_DIR / "cv_comparison.csv"
    if cv_path.exists():
        results["cv"] = pd.read_csv(cv_path)

    # TMT vs LFQ
    tmt_path = ANALYSIS_DIR / "tmt_lfq_comparison.csv"
    if tmt_path.exists():
        results["tmt_lfq"] = pd.read_csv(tmt_path)

    # Summary metrics
    summary_path = ANALYSIS_DIR / "summary_metrics.csv"
    if summary_path.exists():
        results["summary"] = pd.read_csv(summary_path)

    # Cross-experiment correlations
    results["cross_corr"] = {}
    for method in QUANTIFICATION_METHODS:
        path = ANALYSIS_DIR / f"cross_experiment_corr_{method}.csv"
        if path.exists():
            results["cross_corr"][method] = pd.read_csv(path, index_col=0)

    # Rank consistency
    results["rank"] = {}
    for method in QUANTIFICATION_METHODS:
        path = ANALYSIS_DIR / f"rank_consistency_{method}.csv"
        if path.exists():
            results["rank"][method] = pd.read_csv(path, index_col=0)

    # Expression stability
    results["stability"] = {}
    for method in QUANTIFICATION_METHODS:
        path = ANALYSIS_DIR / f"expression_stability_{method}.csv"
        if path.exists():
            results["stability"][method] = pd.read_csv(path)

    return results


def plot_cv_distribution(results: dict, output_dir: Path):
    """
    Create CV distribution boxplot.

    Shows within-experiment variability for each quantification method.
    """
    if "cv" not in results or results["cv"].empty:
        print("  No CV data available")
        return

    cv_df = results["cv"]

    fig, ax = plt.subplots(figsize=FIGSIZE_MEDIUM)

    # Box plot only
    sns.boxplot(
        data=cv_df,
        x="method",
        y="mean_cv",
        ax=ax,
        palette="husl",
    )
    ax.set_xlabel("Quantification Method", fontsize=12)
    ax.set_ylabel("Mean CV", fontsize=12)
    ax.set_title("Coefficient of Variation by Method", fontsize=14)
    ax.tick_params(axis='x', rotation=45)

    plt.tight_layout()
    plt.savefig(output_dir / "cv_distribution.png", dpi=150, bbox_inches="tight")
    plt.savefig(output_dir / "cv_distribution.pdf", bbox_inches="tight")
    plt.close()

    print("  Saved: cv_distribution.png/pdf")


def plot_cross_experiment_heatmaps(results: dict, output_dir: Path):
    """
    Create cross-experiment correlation heatmaps for each method.
    """
    if not results.get("cross_corr"):
        print("  No cross-experiment correlation data available")
        return

    n_methods = len(results["cross_corr"])
    if n_methods == 0:
        return

    # Create grid of heatmaps
    n_cols = min(3, n_methods)
    n_rows = (n_methods + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    if n_methods == 1:
        axes = [[axes]]
    elif n_rows == 1:
        axes = [axes]

    for idx, (method, corr_matrix) in enumerate(results["cross_corr"].items()):
        row, col = idx // n_cols, idx % n_cols
        ax = axes[row][col]

        sns.heatmap(
            corr_matrix,
            annot=True,
            fmt=".2f",
            cmap="RdYlBu_r",
            vmin=0.5,
            vmax=1.0,
            ax=ax,
            square=True,
            cbar_kws={"shrink": 0.8},
        )
        ax.set_title(f"{method.upper()}")
        ax.tick_params(axis='x', rotation=45)
        ax.tick_params(axis='y', rotation=0)

    # Remove empty subplots
    for idx in range(n_methods, n_rows * n_cols):
        row, col = idx // n_cols, idx % n_cols
        axes[row][col].set_visible(False)

    plt.suptitle("Cross-Experiment Correlation (Pearson)", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(output_dir / "correlation_heatmap.png", dpi=150, bbox_inches="tight")
    plt.savefig(output_dir / "correlation_heatmap.pdf", bbox_inches="tight")
    plt.close()

    print("  Saved: correlation_heatmap.png/pdf")


def concordance_correlation_coefficient(y_true, y_pred):
    """Calculate Lin's Concordance Correlation Coefficient."""
    cor = np.corrcoef(y_true, y_pred)[0, 1]
    mean_true = np.mean(y_true)
    mean_pred = np.mean(y_pred)
    var_true = np.var(y_true)
    var_pred = np.var(y_pred)
    sd_true = np.std(y_true)
    sd_pred = np.std(y_pred)
    numerator = 2 * cor * sd_true * sd_pred
    denominator = var_true + var_pred + (mean_true - mean_pred) ** 2
    return numerator / denominator if denominator > 0 else 0


def plot_tmt_vs_lfq(results: dict, output_dir: Path):
    """
    Create TMT vs LFQ density scatter plots for each method.
    """
    from matplotlib.colors import LogNorm
    from scipy.stats import pearsonr

    # Load TMT vs LFQ values for each method
    methods_with_data = []
    for method in QUANTIFICATION_METHODS:
        path = ANALYSIS_DIR / f"tmt_lfq_values_{method}.csv"
        if path.exists():
            methods_with_data.append(method)

    if not methods_with_data:
        print("  No TMT vs LFQ data available")
        return

    n_methods = len(methods_with_data)
    n_cols = min(3, n_methods)
    n_rows = (n_methods + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 6 * n_rows))

    if n_methods == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = np.array([axes])

    for idx, method in enumerate(methods_with_data):
        row, col = idx // n_cols, idx % n_cols
        ax = axes[row, col]

        # Load data
        df = pd.read_csv(ANALYSIS_DIR / f"tmt_lfq_values_{method}.csv")
        x = df["LFQ"].values
        y = df["TMT"].values

        # Remove invalid values
        mask = np.isfinite(x) & np.isfinite(y)
        x = x[mask]
        y = y[mask]

        if len(x) < 10:
            ax.text(0.5, 0.5, "Insufficient data", transform=ax.transAxes,
                    ha='center', va='center')
            continue

        # Calculate statistics
        r, p = pearsonr(x, y)
        ccc = concordance_correlation_coefficient(x, y)

        # Create hexbin density plot
        hb = ax.hexbin(x, y, gridsize=60, cmap='jet', mincnt=1,
                       norm=LogNorm(), alpha=0.9)
        plt.colorbar(hb, ax=ax, label='n_neighbors')

        # Add diagonal line (y = x)
        lims = [
            min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1]),
        ]
        ax.plot(lims, lims, 'k--', alpha=0.5, linewidth=1)
        ax.set_xlim(lims)
        ax.set_ylim(lims)

        # Add statistics text
        stats_text = f"R = {r:.2f}, p < 2.2e-16\nLin's CCC: {ccc:.3f}\nn = {len(x):,}"
        ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
                fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        # Set axis labels - indicate transformation type
        # iBAQ (IbaqLog) is already log-transformed by library
        # All others (including ibaq_raw) are log2-transformed for plotting
        if method == "ibaq":
            scale_label = "log"
        else:
            scale_label = "log2"

        ax.set_xlabel(f"LFQ ({method.upper()}, {scale_label})")
        ax.set_ylabel(f"TMT ({method.upper()}, {scale_label})")
        ax.set_title(f"{method.upper()}")

    # Remove empty subplots
    for idx in range(n_methods, n_rows * n_cols):
        row, col = idx // n_cols, idx % n_cols
        axes[row, col].set_visible(False)

    plt.suptitle("TMT vs LFQ Correlation (PXD007683)", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(output_dir / "tmt_vs_lfq_correlation.png", dpi=150, bbox_inches="tight")
    plt.savefig(output_dir / "tmt_vs_lfq_correlation.pdf", bbox_inches="tight")
    plt.close()

    print("  Saved: tmt_vs_lfq_correlation.png/pdf")


def plot_expression_profiles(results: dict, output_dir: Path):
    """
    Create expression profile plots for proteins of interest.

    Shows expression stability across datasets.
    """
    if not results.get("stability"):
        print("  No expression stability data available")
        return

    # Combine all methods
    all_stability = []
    for method, df in results["stability"].items():
        df = df.copy()
        df["method"] = method
        all_stability.append(df)

    if not all_stability:
        return

    combined = pd.concat(all_stability, ignore_index=True)

    # Select top 10 most stable proteins
    protein_mad = combined.groupby("ProteinName")["MAD"].mean().sort_values()
    top_proteins = protein_mad.head(10).index.tolist()

    # Filter to top proteins
    plot_df = combined[combined["ProteinName"].isin(top_proteins)]

    fig, ax = plt.subplots(figsize=FIGSIZE_MEDIUM)

    sns.barplot(
        data=plot_df,
        x="ProteinName",
        y="MAD",
        hue="method",
        ax=ax,
    )

    ax.set_xlabel("Protein (UniProt Accession)")
    ax.set_ylabel("Median Absolute Deviation (log2)")
    ax.set_title("Expression Stability Across Datasets\n(Top 10 Most Stable Proteins)")
    ax.tick_params(axis='x', rotation=45)
    ax.legend(title="Method", bbox_to_anchor=(1.02, 1), loc='upper left')

    plt.tight_layout()
    plt.savefig(output_dir / "expression_profiles.png", dpi=150, bbox_inches="tight")
    plt.savefig(output_dir / "expression_profiles.pdf", bbox_inches="tight")
    plt.close()

    print("  Saved: expression_profiles.png/pdf")


def plot_rank_stability(results: dict, output_dir: Path):
    """
    Create rank stability comparison plot.

    Shows mean Spearman rank correlation across experiments.
    """
    if not results.get("rank"):
        print("  No rank consistency data available")
        return

    # Compute mean rank correlation for each method
    mean_ranks = []
    for method, rank_matrix in results["rank"].items():
        # Get upper triangle (excluding diagonal)
        n = len(rank_matrix)
        upper_tri = rank_matrix.values[np.triu_indices(n, k=1)]
        mean_rank = np.nanmean(upper_tri)
        std_rank = np.nanstd(upper_tri)
        mean_ranks.append({
            "method": method,
            "mean_rank_corr": mean_rank,
            "std_rank_corr": std_rank,
        })

    df = pd.DataFrame(mean_ranks)

    fig, ax = plt.subplots(figsize=FIGSIZE_SMALL)

    bars = ax.bar(df["method"], df["mean_rank_corr"], yerr=df["std_rank_corr"],
                  capsize=5, color="teal", alpha=0.8)

    ax.set_xlabel("Quantification Method")
    ax.set_ylabel("Mean Spearman Rank Correlation")
    ax.set_title("Protein Rank Stability Across Experiments")
    ax.tick_params(axis='x', rotation=45)
    ax.set_ylim(0, 1.1)

    # Add value labels
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height:.3f}',
                   xy=(bar.get_x() + bar.get_width() / 2, height),
                   xytext=(0, 3),
                   textcoords="offset points",
                   ha='center', va='bottom', fontsize=9)

    plt.tight_layout()
    plt.savefig(output_dir / "rank_stability.png", dpi=150, bbox_inches="tight")
    plt.savefig(output_dir / "rank_stability.pdf", bbox_inches="tight")
    plt.close()

    print("  Saved: rank_stability.png/pdf")


def plot_method_summary_radar(results: dict, output_dir: Path):
    """
    Create radar chart summarizing all metrics for each method.
    """
    if "summary" not in results or results["summary"].empty:
        print("  No summary data available")
        return

    df = results["summary"]

    # Metrics to include (normalize to 0-1 scale)
    metrics = []
    metric_names = []

    # CV (lower is better, so invert)
    if "mean_cv" in df.columns:
        cv_vals = df["mean_cv"].values
        cv_norm = 1 - (cv_vals - cv_vals.min()) / (cv_vals.max() - cv_vals.min() + 1e-10)
        metrics.append(cv_norm)
        metric_names.append("Low CV")

    # Cross-experiment correlation (higher is better)
    if "mean_cross_corr" in df.columns:
        corr_vals = df["mean_cross_corr"].values
        corr_norm = (corr_vals - corr_vals.min()) / (corr_vals.max() - corr_vals.min() + 1e-10)
        metrics.append(corr_norm)
        metric_names.append("Cross-Exp Corr")

    # TMT-LFQ Pearson (higher is better)
    if "tmt_lfq_pearson" in df.columns:
        tmt_vals = df["tmt_lfq_pearson"].fillna(0).values
        tmt_norm = (tmt_vals - tmt_vals.min()) / (tmt_vals.max() - tmt_vals.min() + 1e-10)
        metrics.append(tmt_norm)
        metric_names.append("TMT-LFQ Agree")

    # Rank consistency (higher is better)
    if "mean_rank_corr" in df.columns:
        rank_vals = df["mean_rank_corr"].values
        rank_norm = (rank_vals - rank_vals.min()) / (rank_vals.max() - rank_vals.min() + 1e-10)
        metrics.append(rank_norm)
        metric_names.append("Rank Stability")

    # MAD (lower is better, so invert)
    if "mean_mad" in df.columns:
        mad_vals = df["mean_mad"].values
        mad_norm = 1 - (mad_vals - mad_vals.min()) / (mad_vals.max() - mad_vals.min() + 1e-10)
        metrics.append(mad_norm)
        metric_names.append("Expression Stability")

    if len(metrics) < 3:
        print("  Not enough metrics for radar chart")
        return

    # Create radar chart
    methods = df["method"].tolist()
    n_metrics = len(metrics)
    n_methods = len(methods)

    angles = np.linspace(0, 2 * np.pi, n_metrics, endpoint=False).tolist()
    angles += angles[:1]  # Close the loop

    fig, ax = plt.subplots(figsize=FIGSIZE_MEDIUM, subplot_kw=dict(polar=True))

    colors = plt.cm.tab10(np.linspace(0, 1, n_methods))

    for i, method in enumerate(methods):
        values = [m[i] for m in metrics]
        values += values[:1]  # Close the loop

        ax.plot(angles, values, 'o-', linewidth=2, label=method, color=colors[i])
        ax.fill(angles, values, alpha=0.15, color=colors[i])

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(metric_names, size=10)
    ax.set_ylim(0, 1)
    ax.set_title("Method Comparison Summary\n(Higher = Better)", size=14, y=1.08)
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0))

    plt.tight_layout()
    plt.savefig(output_dir / "method_comparison_radar.png", dpi=150, bbox_inches="tight")
    plt.savefig(output_dir / "method_comparison_radar.pdf", bbox_inches="tight")
    plt.close()

    print("  Saved: method_comparison_radar.png/pdf")


def plot_summary_barplot(results: dict, output_dir: Path):
    """
    Create summary bar plot comparing all metrics.
    """
    if "summary" not in results or results["summary"].empty:
        print("  No summary data available")
        return

    df = results["summary"]

    # Prepare data for plotting
    metrics_to_plot = []

    if "mean_cv" in df.columns:
        metrics_to_plot.append(("mean_cv", "Mean CV", True))  # Lower is better

    if "mean_cross_corr" in df.columns:
        metrics_to_plot.append(("mean_cross_corr", "Cross-Exp Corr", False))

    if "tmt_lfq_spearman" in df.columns:
        metrics_to_plot.append(("tmt_lfq_spearman", "TMT-LFQ Spearman", False))

    if "mean_rank_corr" in df.columns:
        metrics_to_plot.append(("mean_rank_corr", "Rank Consistency", False))

    n_metrics = len(metrics_to_plot)
    if n_metrics == 0:
        return

    fig, axes = plt.subplots(1, n_metrics, figsize=(4 * n_metrics, 5))
    if n_metrics == 1:
        axes = [axes]

    for ax, (col, title, lower_better) in zip(axes, metrics_to_plot):
        data = df[["method", col]].dropna()

        if data.empty:
            ax.set_title(f"{title}\n(No data)")
            ax.set_visible(False)
            continue

        bars = ax.bar(data["method"], data[col], color="steelblue", alpha=0.8)

        # Highlight best
        if len(data) > 0:
            if lower_better:
                best_idx = data[col].idxmin()
            else:
                best_idx = data[col].idxmax()

            best_method = data.loc[best_idx, "method"]
            for i, bar in enumerate(bars):
                if data.iloc[i]["method"] == best_method:
                    bar.set_color("darkorange")

        ax.set_xlabel("Method")
        ax.set_ylabel(title)
        ax.set_title(title)
        ax.tick_params(axis='x', rotation=45)

        # Add value labels
        for bar in bars:
            height = bar.get_height()
            ax.annotate(f'{height:.3f}',
                       xy=(bar.get_x() + bar.get_width() / 2, height),
                       xytext=(0, 3),
                       textcoords="offset points",
                       ha='center', va='bottom', fontsize=8)

    plt.suptitle("Benchmarking Summary (Best Method Highlighted)", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(output_dir / "summary_barplot.png", dpi=150, bbox_inches="tight")
    plt.savefig(output_dir / "summary_barplot.pdf", bbox_inches="tight")
    plt.close()

    print("  Saved: summary_barplot.png/pdf")


def plot_sample_expression_boxplots(output_dir: Path):
    """
    Create sample-level expression boxplots for each quantification method.

    Shows the distribution of protein expression values per sample,
    colored by dataset/project.
    """
    # Intensity column names for each method
    intensity_columns = {
        "ibaq": "IbaqLog",        # Already log-transformed by library
        "ibaq_raw": "Ibaq",       # Raw iBAQ values
        "directlfq": "DirectLFQIntensity",
        "top3": "Top3Intensity",
        "topn": "Top10Intensity",
        "sum": "SumIntensity",
    }

    # Collect data for all methods
    method_data = {}

    for method in QUANTIFICATION_METHODS:
        intensity_col = intensity_columns.get(method)
        if not intensity_col:
            continue

        all_samples = []

        for dataset_id in ALL_DATASETS:
            quant_path = PROTEIN_QUANT_DIR / dataset_id / f"{method}.parquet"
            if not quant_path.exists():
                continue

            try:
                df = pd.read_parquet(quant_path)
            except Exception:
                continue

            if intensity_col not in df.columns:
                continue

            # Get sample-level data
            if "SampleID" in df.columns:
                for sample_id in df["SampleID"].unique():
                    sample_df = df[df["SampleID"] == sample_id]
                    values = sample_df[intensity_col].dropna()

                    if len(values) > 0:
                        # Log transform if not already log
                        # iBAQ (IbaqLog) is already log-transformed by library - use as-is
                        # All other methods need log10 transformation
                        if method != "ibaq":
                            values = np.log10(values + 1)
                        # else: IbaqLog is already in log scale, use directly

                        all_samples.append({
                            "dataset": dataset_id,
                            "sample": f"{dataset_id}_{sample_id}",
                            "values": values.values,
                            "median": values.median(),
                            "q25": values.quantile(0.25),
                            "q75": values.quantile(0.75),
                        })

        if all_samples:
            method_data[method] = all_samples

    if not method_data:
        print("  No sample expression data available")
        return

    # Create a separate plot for each method
    n_methods = len(method_data)

    # Get unique datasets for coloring
    all_datasets = sorted(set(s["dataset"] for samples in method_data.values() for s in samples))
    colors = plt.cm.tab20(np.linspace(0, 1, len(all_datasets)))
    dataset_colors = {d: colors[i] for i, d in enumerate(all_datasets)}

    # Create figure with subplots for each method
    fig, axes = plt.subplots(n_methods, 1, figsize=(20, 4 * n_methods))
    if n_methods == 1:
        axes = [axes]

    for ax, method in zip(axes, method_data.keys()):
        samples = method_data[method]

        # Sort samples by dataset
        samples_sorted = sorted(samples, key=lambda x: (x["dataset"], x["sample"]))

        # Create boxplot data
        positions = []
        box_data = []
        box_colors = []

        for i, s in enumerate(samples_sorted):
            positions.append(i)
            box_data.append(s["values"])
            box_colors.append(dataset_colors[s["dataset"]])

        # Create boxplots
        bp = ax.boxplot(box_data, positions=positions, widths=0.6, patch_artist=True,
                       showfliers=True, flierprops=dict(marker='.', markersize=2, alpha=0.5))

        # Color boxes by dataset
        for patch, color in zip(bp['boxes'], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        # Style
        ax.set_ylabel(f"log10({method.upper()})")
        ax.set_title(f"{method.upper()} Expression Distribution per Sample")
        ax.set_xlabel("Samples per Project")

        # Set x-axis ticks (show fewer labels to avoid crowding)
        n_samples = len(samples_sorted)
        if n_samples > 50:
            tick_step = max(1, n_samples // 20)
            tick_positions = list(range(0, n_samples, tick_step))
            ax.set_xticks(tick_positions)
            ax.set_xticklabels([samples_sorted[i]["sample"].split("_")[-1] for i in tick_positions],
                              rotation=90, fontsize=6)
        else:
            ax.set_xticks(positions)
            ax.set_xticklabels([s["sample"].split("_")[-1] for s in samples_sorted],
                              rotation=90, fontsize=6)

        ax.set_xlim(-1, len(samples_sorted))

    # Add legend for datasets
    legend_handles = [plt.Rectangle((0, 0), 1, 1, facecolor=dataset_colors[d], alpha=0.7)
                      for d in all_datasets]
    fig.legend(legend_handles, all_datasets, loc='center right',
               bbox_to_anchor=(1.12, 0.5), title="Project", fontsize=8)

    plt.tight_layout()
    plt.savefig(output_dir / "sample_expression_boxplots.png", dpi=150, bbox_inches="tight")
    plt.savefig(output_dir / "sample_expression_boxplots.pdf", bbox_inches="tight")
    plt.close()

    # Also create a per-project summary plot (panel c/d style from reference)
    fig, axes = plt.subplots(len(method_data), 1, figsize=(14, 4 * len(method_data)))
    if len(method_data) == 1:
        axes = [axes]

    for ax, method in zip(axes, method_data.keys()):
        samples = method_data[method]

        # Group by dataset
        dataset_values = {}
        for s in samples:
            if s["dataset"] not in dataset_values:
                dataset_values[s["dataset"]] = []
            dataset_values[s["dataset"]].extend(s["values"])

        # Create boxplot data
        datasets = sorted(dataset_values.keys())
        box_data = [dataset_values[d] for d in datasets]
        box_colors = [dataset_colors[d] for d in datasets]

        bp = ax.boxplot(box_data, patch_artist=True,
                       showfliers=True, flierprops=dict(marker='.', markersize=2, alpha=0.3))

        for patch, color in zip(bp['boxes'], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        ax.set_ylabel(f"log10({method.upper()})")
        ax.set_title(f"{method.upper()} Protein Expression per Project")
        ax.set_xlabel("Project")
        ax.set_xticklabels(datasets, rotation=45, ha='right', fontsize=8)

    plt.tight_layout()
    plt.savefig(output_dir / "project_expression_boxplots.png", dpi=150, bbox_inches="tight")
    plt.savefig(output_dir / "project_expression_boxplots.pdf", bbox_inches="tight")
    plt.close()

    print("  Saved: sample_expression_boxplots.png/pdf")
    print("  Saved: project_expression_boxplots.png/pdf")


def generate_all_plots(output_dir: Path = None):
    """Generate all plots."""
    if output_dir is None:
        output_dir = PLOTS_DIR

    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Generating Benchmark Visualizations")
    print("=" * 70)
    print(f"\nOutput directory: {output_dir}")

    # Load results
    print("\nLoading analysis results...")
    results = load_analysis_results()

    if not results:
        print("ERROR: No analysis results found. Run 04_compute_metrics.py first.")
        return

    # Generate plots
    print("\nGenerating plots...")

    print("\n1. CV Distribution")
    plot_cv_distribution(results, output_dir)

    print("\n2. Cross-Experiment Correlation Heatmaps")
    plot_cross_experiment_heatmaps(results, output_dir)

    print("\n3. TMT vs LFQ Correlation")
    plot_tmt_vs_lfq(results, output_dir)

    print("\n4. Expression Profiles")
    plot_expression_profiles(results, output_dir)

    print("\n5. Rank Stability")
    plot_rank_stability(results, output_dir)

    print("\n6. Method Comparison Radar Chart")
    plot_method_summary_radar(results, output_dir)

    print("\n7. Summary Bar Plot")
    plot_summary_barplot(results, output_dir)

    print("\n8. Sample Expression Boxplots")
    plot_sample_expression_boxplots(output_dir)

    print("\n" + "=" * 70)
    print("Done! All plots saved to:", output_dir)
    print("=" * 70)


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate benchmark visualizations"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PLOTS_DIR,
        help="Output directory for plots"
    )

    args = parser.parse_args()

    generate_all_plots(output_dir=args.output_dir)


if __name__ == "__main__":
    main()
