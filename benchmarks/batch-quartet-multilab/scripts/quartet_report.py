"""Metrics and figures for the refreshed Quartet Rust benchmark."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


plt.switch_backend("Agg")

METHODS = ("pibaq", "maxlfq", "directlfq", "top3")
METHOD_LABELS = dict(zip(METHODS, ("piBAQ", "MaxLFQ", "DirectLFQ", "Top3")))
SAMPLE_TYPES = ("D5", "D6", "F7", "M8")
EXPECTED_BATCHES = (
    "DDA_APT",
    "DIA_APT",
    "DIA_BGI",
    "DDA_FDU",
    "DIA_FDU",
    "DDA_NVG",
)
SAMPLE_COLORS = {
    "D5": "#DC2626",
    "D6": "#2563EB",
    "F7": "#16A34A",
    "M8": "#9333EA",
}
MODE_MARKERS = {"DDA": "o", "DIA": "s"}


def technical_cv(matrix_log2: pd.DataFrame, metadata: pd.DataFrame) -> float:
    """Return the mean sample-group median CV on the linear scale."""
    linear = np.exp2(matrix_log2)
    group_medians = []
    for sample_type in SAMPLE_TYPES:
        columns = metadata.loc[metadata["sample"] == sample_type, "run_id"]
        values = linear.loc[:, columns]
        cv = values.std(axis=1, ddof=1) / values.mean(axis=1)
        group_medians.append(float(cv[np.isfinite(cv)].median()))
    return float(np.mean(group_medians))


def biological_snr(matrix_log2: pd.DataFrame, metadata: pd.DataFrame) -> float:
    """Measure PCA separation between groups relative to within-group spread."""
    coordinates = PCA(n_components=2).fit_transform(matrix_log2.T)
    labels = metadata.set_index("run_id").loc[matrix_log2.columns, "sample"].to_numpy()
    within = []
    between = []
    for left, right in combinations(range(len(coordinates)), 2):
        distance = float(np.linalg.norm(coordinates[left] - coordinates[right]))
        target = within if labels[left] == labels[right] else between
        target.append(distance)
    noise = np.sqrt(np.mean(np.square(within)))
    signal = np.sqrt(np.mean(np.square(between)))
    return float(20 * np.log10(signal / noise))


def inter_batch_correlation(matrix_log2: pd.DataFrame, metadata: pd.DataFrame) -> float:
    """Correlate protein-by-sample-type profiles between acquisition batches."""
    profiles = {}
    for batch in EXPECTED_BATCHES:
        batch_meta = metadata.loc[metadata["batch"] == batch]
        values = []
        for sample_type in SAMPLE_TYPES:
            columns = batch_meta.loc[batch_meta["sample"] == sample_type, "run_id"]
            values.append(matrix_log2.loc[:, columns].median(axis=1).to_numpy())
        profiles[batch] = np.concatenate(values)
    correlations = [
        np.corrcoef(profiles[left], profiles[right])[0, 1]
        for left, right in combinations(EXPECTED_BATCHES, 2)
    ]
    return float(np.mean(correlations))


def batch_rmse(matrix_log2: pd.DataFrame, metadata: pd.DataFrame) -> float:
    """Measure batch-centroid residuals after removing sample-type means."""
    sample_means = {}
    for sample_type in SAMPLE_TYPES:
        columns = metadata.loc[metadata["sample"] == sample_type, "run_id"]
        sample_means[sample_type] = matrix_log2.loc[:, columns].mean(axis=1)
    residuals = []
    for batch in EXPECTED_BATCHES:
        batch_meta = metadata.loc[metadata["batch"] == batch]
        for sample_type in SAMPLE_TYPES:
            columns = batch_meta.loc[batch_meta["sample"] == sample_type, "run_id"]
            centroid = matrix_log2.loc[:, columns].mean(axis=1)
            residuals.extend((centroid - sample_means[sample_type]).to_numpy())
    return float(np.sqrt(np.mean(np.square(residuals))))


def calculate_metrics(
    raw: dict[str, pd.DataFrame],
    corrected: dict[str, pd.DataFrame],
    coverage: dict[str, int],
    complete_coverage: dict[str, int],
    metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Compute matched-universe diagnostics without declaring a winner."""
    records = []
    for method in METHODS:
        for correction, matrix in (("Raw", raw[method]), ("ComBat", corrected[method])):
            records.append(
                {
                    "Quantification": METHOD_LABELS[method],
                    "Batch_Correction": correction,
                    "Total_Proteins": coverage[method],
                    "Complete_Proteins": complete_coverage[method],
                    "Common_Proteins": len(matrix),
                    "Mean_CV": technical_cv(matrix, metadata),
                    "SNR": biological_snr(matrix, metadata),
                    "Inter_Batch_Correlation": inter_batch_correlation(
                        matrix, metadata
                    ),
                    "Batch_RMSE": batch_rmse(matrix, metadata),
                }
            )
    metrics = pd.DataFrame.from_records(records)
    finite_columns = ["Mean_CV", "SNR", "Inter_Batch_Correlation", "Batch_RMSE"]
    if not np.isfinite(metrics[finite_columns]).all().all():
        raise ValueError("Quartet metrics contain non-finite values")
    return metrics


def save_results(
    results_dir: Path,
    metadata: pd.DataFrame,
    raw: dict[str, pd.DataFrame],
    corrected: dict[str, pd.DataFrame],
    metrics: pd.DataFrame,
) -> None:
    """Write the refreshed tracked result contract."""
    metadata.drop(columns="source_folder").to_csv(
        results_dir / "metadata.csv", index=False
    )
    metrics.to_csv(results_dir / "benchmark_metrics.csv", index=False)
    for method in METHODS:
        raw[method].rename_axis("ProteinName").to_csv(results_dir / f"{method}_raw.csv")
        corrected[method].rename_axis("ProteinName").to_csv(
            results_dir / f"{method}_combat.csv"
        )


def plot_batch_effect_diagnosis(
    raw: dict[str, pd.DataFrame],
    corrected: dict[str, pd.DataFrame],
    metadata: pd.DataFrame,
    output: Path,
) -> None:
    """Plot the original two-panel MaxLFQ sample-median diagnosis."""
    figure, axes = plt.subplots(1, 2, figsize=(14, 6))
    colors = plt.get_cmap("tab10")(np.linspace(0, 1, len(EXPECTED_BATCHES)))
    batch_colors = dict(zip(EXPECTED_BATCHES, colors))
    sample_metadata = metadata.set_index("run_id")
    panels = (
        (raw["maxlfq"], "Before correction"),
        (corrected["maxlfq"], "After Rust ComBat"),
    )
    for axis, (matrix, title) in zip(axes, panels):
        plot_sample_medians(axis, matrix, title, sample_metadata, batch_colors)
    figure.tight_layout()
    figure.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(figure)


def plot_sample_medians(
    axis: plt.Axes,
    matrix: pd.DataFrame,
    title: str,
    sample_metadata: pd.DataFrame,
    batch_colors: dict[str, object],
) -> None:
    """Draw one raw/corrected sample-median panel."""
    sample_medians = matrix.median(axis=0)
    ordered_metadata = sample_metadata.loc[sample_medians.index]
    bar_colors = [batch_colors[batch] for batch in ordered_metadata["batch"]]
    axis.bar(
        np.arange(len(sample_medians)),
        sample_medians.to_numpy(),
        color=bar_colors,
        alpha=0.8,
    )
    axis.set_xlabel("Samples")
    axis.set_ylabel("Median log2 intensity")
    axis.set_title(f"MaxLFQ - {title}")
    axis.set_xticks([])
    for batch, color in batch_colors.items():
        axis.bar([], [], color=color, label=batch)
    axis.legend(loc="upper right", fontsize=8, title="Batch")


def plot_pca_panel(
    axis: plt.Axes,
    matrix: pd.DataFrame,
    method: str,
    correction: str,
    sample_metadata: pd.DataFrame,
) -> None:
    """Draw one method/correction PCA panel."""
    pca = PCA(n_components=2)
    coordinates = pca.fit_transform(matrix.T)
    ordered_metadata = sample_metadata.loc[matrix.columns]
    for sample_type in SAMPLE_TYPES:
        for mode, marker in MODE_MARKERS.items():
            mask = (ordered_metadata["sample"] == sample_type) & (
                ordered_metadata["mode"] == mode
            )
            axis.scatter(
                coordinates[mask.to_numpy(), 0],
                coordinates[mask.to_numpy(), 1],
                color=SAMPLE_COLORS[sample_type],
                marker=marker,
                s=42,
                alpha=0.72,
                label=(
                    f"{sample_type} {mode}"
                    if correction == "Raw" and method == METHODS[0]
                    else None
                ),
            )
    axis.set_title(f"{METHOD_LABELS[method]} — {correction}")
    axis.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
    axis.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)")
    axis.grid(alpha=0.2)


def plot_pca(
    raw: dict[str, pd.DataFrame],
    corrected: dict[str, pd.DataFrame],
    metadata: pd.DataFrame,
    output: Path,
) -> None:
    """Plot raw and corrected PCA for each current Rust quantification method."""
    figure, axes = plt.subplots(2, len(METHODS), figsize=(20, 10))
    sample_metadata = metadata.set_index("run_id")
    for column, method in enumerate(METHODS):
        for row, (correction, matrix) in enumerate(
            (("Raw", raw[method]), ("Rust ComBat", corrected[method]))
        ):
            plot_pca_panel(
                axes[row, column], matrix, method, correction, sample_metadata
            )
    axes[0, 0].legend(frameon=False, fontsize=8, ncol=2)
    figure.suptitle("Quartet PCA: sample type color, acquisition mode shape")
    figure.tight_layout()
    figure.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(figure)
