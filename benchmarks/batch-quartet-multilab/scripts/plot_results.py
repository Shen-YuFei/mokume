#!/usr/bin/env python
"""
Generate visualizations for the batch effect correction benchmark.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# Paths
BENCHMARK_DIR = Path(__file__).parent.parent
RESULTS_DIR = BENCHMARK_DIR / "results"
FIGURES_DIR = BENCHMARK_DIR / "figures"

FIGURES_DIR.mkdir(exist_ok=True)


def plot_metrics_comparison():
    """Plot CV and SNR comparison."""
    metrics = pd.read_csv(RESULTS_DIR / "benchmark_metrics.csv")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # CV comparison
    ax = axes[0]
    methods = metrics["Quantification"].unique()
    x = np.arange(len(methods))
    width = 0.35

    # Get values in correct order (None is stored as NaN in CSV)
    raw_cv = []
    combat_cv = []
    for method in methods:
        raw_val = metrics[(metrics["Quantification"] == method) & (metrics["Batch_Correction"].isna())]["Mean_CV"].values
        combat_val = metrics[(metrics["Quantification"] == method) & (metrics["Batch_Correction"] == "ComBat")]["Mean_CV"].values
        raw_cv.append(raw_val[0] if len(raw_val) > 0 else 0)
        combat_cv.append(combat_val[0] if len(combat_val) > 0 else 0)

    raw_cv = np.array(raw_cv)
    combat_cv = np.array(combat_cv)

    bars1 = ax.bar(x - width/2, raw_cv, width, label="No Correction", color="#1f77b4", alpha=0.8)
    bars2 = ax.bar(x + width/2, combat_cv, width, label="ComBat", color="#ff7f0e", alpha=0.8)

    ax.set_ylabel("Mean CV (lower is better)")
    ax.set_xlabel("Quantification Method")
    ax.set_title("Coefficient of Variation")
    ax.set_xticks(x)
    ax.set_xticklabels(methods)
    ax.legend()
    ax.grid(axis="y", alpha=0.3)

    # Add value labels
    for bar in bars1:
        height = bar.get_height()
        ax.annotate(f'{height:.3f}', xy=(bar.get_x() + bar.get_width()/2, height),
                    xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=9)
    for bar in bars2:
        height = bar.get_height()
        ax.annotate(f'{height:.3f}', xy=(bar.get_x() + bar.get_width()/2, height),
                    xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=9)

    # SNR comparison
    ax = axes[1]
    raw_snr = []
    combat_snr = []
    for method in methods:
        raw_val = metrics[(metrics["Quantification"] == method) & (metrics["Batch_Correction"].isna())]["SNR"].values
        combat_val = metrics[(metrics["Quantification"] == method) & (metrics["Batch_Correction"] == "ComBat")]["SNR"].values
        raw_snr.append(raw_val[0] if len(raw_val) > 0 else 0)
        combat_snr.append(combat_val[0] if len(combat_val) > 0 else 0)

    raw_snr = np.array(raw_snr)
    combat_snr = np.array(combat_snr)

    bars1 = ax.bar(x - width/2, raw_snr, width, label="No Correction", color="#1f77b4", alpha=0.8)
    bars2 = ax.bar(x + width/2, combat_snr, width, label="ComBat", color="#ff7f0e", alpha=0.8)

    ax.set_ylabel("SNR (dB, higher is better)")
    ax.set_xlabel("Quantification Method")
    ax.set_title("Signal-to-Noise Ratio")
    ax.set_xticks(x)
    ax.set_xticklabels(methods)
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)

    # Add value labels
    for bar in bars1:
        height = bar.get_height()
        ax.annotate(f'{height:.2f}', xy=(bar.get_x() + bar.get_width()/2, height),
                    xytext=(0, 3 if height >= 0 else -12), textcoords="offset points",
                    ha='center', va='bottom' if height >= 0 else 'top', fontsize=9)
    for bar in bars2:
        height = bar.get_height()
        ax.annotate(f'{height:.2f}', xy=(bar.get_x() + bar.get_width()/2, height),
                    xytext=(0, 3 if height >= 0 else -12), textcoords="offset points",
                    ha='center', va='bottom' if height >= 0 else 'top', fontsize=9)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "metrics_comparison.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {FIGURES_DIR / 'metrics_comparison.png'}")


def plot_pca_comparison():
    """Plot PCA before and after batch correction."""
    meta = pd.read_csv(RESULTS_DIR / "metadata.csv")

    methods = ["maxlfq", "top3", "ibaq", "directlfq"]
    titles = ["MaxLFQ", "Top3", "iBAQ", "DirectLFQ"]

    # Filter to only methods that have data
    available_methods = []
    available_titles = []
    for method, title in zip(methods, titles):
        if (RESULTS_DIR / f"{method}_raw.csv").exists():
            available_methods.append(method)
            available_titles.append(title)

    n_methods = len(available_methods)
    fig, axes = plt.subplots(2, n_methods, figsize=(5*n_methods, 10))

    methods = available_methods
    titles = available_titles

    for i, (method, title) in enumerate(zip(methods, titles)):
        for j, (suffix, correction_title) in enumerate([("raw", "No Correction"), ("combat", "ComBat")]):
            ax = axes[j, i]

            try:
                matrix = pd.read_csv(RESULTS_DIR / f"{method}_{suffix}.csv", index_col=0)
            except FileNotFoundError:
                ax.text(0.5, 0.5, "Data not available", ha='center', va='center')
                ax.set_title(f"{title} - {correction_title}")
                continue

            # Get complete cases
            matrix_complete = matrix.dropna()

            if len(matrix_complete) < 10:
                ax.text(0.5, 0.5, "Insufficient data", ha='center', va='center')
                ax.set_title(f"{title} - {correction_title}")
                continue

            # PCA
            pca = PCA(n_components=2)
            pca_coords = pca.fit_transform(matrix_complete.T)

            # Get sample metadata
            samples = matrix_complete.columns.tolist()
            sample_meta = meta[meta["run_id"].isin(samples)].set_index("run_id").loc[samples]

            # Color by sample type
            colors = {"D5": "#e41a1c", "D6": "#377eb8", "F7": "#4daf4a", "M8": "#984ea3"}
            markers = {"DDA": "o", "DIA": "s"}

            for sample_type in ["D5", "D6", "F7", "M8"]:
                for mode in ["DDA", "DIA"]:
                    mask = (sample_meta["sample"] == sample_type) & (sample_meta["mode"] == mode)
                    if mask.sum() > 0:
                        ax.scatter(
                            pca_coords[mask.values, 0],
                            pca_coords[mask.values, 1],
                            c=colors[sample_type],
                            marker=markers[mode],
                            s=60,
                            alpha=0.7,
                            label=f"{sample_type} ({mode})" if j == 0 and i == 0 else ""
                        )

            ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
            ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
            ax.set_title(f"{title} - {correction_title}")
            ax.grid(alpha=0.3)

    # Add legend
    axes[0, 0].legend(loc="upper left", fontsize=8, ncol=2)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "pca_comparison.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {FIGURES_DIR / 'pca_comparison.png'}")


def plot_batch_effect_diagnosis():
    """Plot batch effect before and after correction."""
    meta = pd.read_csv(RESULTS_DIR / "metadata.csv")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Use MaxLFQ as example
    try:
        raw_matrix = pd.read_csv(RESULTS_DIR / "maxlfq_raw.csv", index_col=0)
        combat_matrix = pd.read_csv(RESULTS_DIR / "maxlfq_combat.csv", index_col=0)
    except FileNotFoundError:
        print("MaxLFQ results not found")
        return

    for ax, (matrix, title) in zip(axes, [(raw_matrix, "Before Correction"), (combat_matrix, "After ComBat")]):
        # Calculate median intensity per sample
        sample_medians = matrix.median()

        # Get batch info
        samples = sample_medians.index.tolist()
        sample_meta = meta[meta["run_id"].isin(samples)].set_index("run_id").loc[samples]

        # Color by batch
        batches = sample_meta["batch"].unique()
        colors = plt.cm.tab10(np.linspace(0, 1, len(batches)))
        batch_colors = dict(zip(batches, colors))

        x = np.arange(len(samples))
        batch_cols = [batch_colors[b] for b in sample_meta["batch"]]

        ax.bar(x, sample_medians.values, color=batch_cols, alpha=0.8)
        ax.set_xlabel("Samples")
        ax.set_ylabel("Median Log Intensity")
        ax.set_title(f"MaxLFQ - {title}")
        ax.set_xticks([])

        # Add batch legend
        for batch, color in batch_colors.items():
            ax.bar([], [], color=color, label=batch)
        ax.legend(loc="upper right", fontsize=8, title="Batch")

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "batch_effect_diagnosis.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {FIGURES_DIR / 'batch_effect_diagnosis.png'}")


def plot_improvement_summary():
    """Plot summary of improvement from batch correction."""
    metrics = pd.read_csv(RESULTS_DIR / "benchmark_metrics.csv")

    fig, ax = plt.subplots(figsize=(10, 6))

    methods = metrics["Quantification"].unique()

    # Calculate improvement
    improvements = []
    for method in methods:
        raw_cv = metrics[(metrics["Quantification"] == method) & (metrics["Batch_Correction"].isna())]["Mean_CV"].values
        combat_cv = metrics[(metrics["Quantification"] == method) & (metrics["Batch_Correction"] == "ComBat")]["Mean_CV"].values
        raw_snr = metrics[(metrics["Quantification"] == method) & (metrics["Batch_Correction"].isna())]["SNR"].values
        combat_snr = metrics[(metrics["Quantification"] == method) & (metrics["Batch_Correction"] == "ComBat")]["SNR"].values

        if len(raw_cv) > 0 and len(combat_cv) > 0:
            cv_improvement = (raw_cv[0] - combat_cv[0]) / raw_cv[0] * 100
        else:
            cv_improvement = 0

        if len(raw_snr) > 0 and len(combat_snr) > 0:
            snr_improvement = combat_snr[0] - raw_snr[0]
        else:
            snr_improvement = 0

        improvements.append({
            "Method": method,
            "CV_Reduction_%": cv_improvement,
            "SNR_Gain_dB": snr_improvement
        })

    imp_df = pd.DataFrame(improvements)

    x = np.arange(len(methods))
    width = 0.35

    bars1 = ax.bar(x - width/2, imp_df["CV_Reduction_%"], width, label="CV Reduction (%)", color="#2ca02c", alpha=0.8)
    bars2 = ax.bar(x + width/2, imp_df["SNR_Gain_dB"], width, label="SNR Gain (dB)", color="#9467bd", alpha=0.8)

    ax.set_ylabel("Improvement")
    ax.set_xlabel("Quantification Method")
    ax.set_title("Improvement from Protein-Level Batch Correction (ComBat)")
    ax.set_xticks(x)
    ax.set_xticklabels(methods)
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    ax.axhline(y=0, color='k', linestyle='-', alpha=0.3)

    # Add value labels
    for bar in bars1:
        height = bar.get_height()
        ax.annotate(f'{height:.1f}%', xy=(bar.get_x() + bar.get_width()/2, height),
                    xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=10)
    for bar in bars2:
        height = bar.get_height()
        ax.annotate(f'{height:.2f}', xy=(bar.get_x() + bar.get_width()/2, height),
                    xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=10)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "improvement_summary.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {FIGURES_DIR / 'improvement_summary.png'}")

    # Print summary
    print("\n" + "=" * 70)
    print("IMPROVEMENT SUMMARY")
    print("=" * 70)
    print(imp_df.to_string(index=False))


def main():
    print("=" * 70)
    print("Generating benchmark visualizations...")
    print("=" * 70)

    plot_metrics_comparison()
    plot_pca_comparison()
    plot_batch_effect_diagnosis()
    plot_improvement_summary()

    print("\nAll figures saved to:", FIGURES_DIR)


if __name__ == "__main__":
    main()
