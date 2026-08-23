"""Metrics and figures for the refreshed HeLa quantification benchmark."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from config import (
    IBAQPY_DATASETS,
    PROTEINS_OF_INTEREST,
    QUANTIFICATION_METHODS,
    TMT_LFQ_DATASETS,
)


plt.switch_backend("Agg")

METHODS = tuple(QUANTIFICATION_METHODS)
METHOD_LABELS = {
    "pibaq": "piBAQ",
    "maxlfq": "MaxLFQ",
    "directlfq": "DirectLFQ",
    "top3": "Top3",
    "top10": "Top10",
    "sum": "Sum",
}
MAIN_DATASETS = tuple(IBAQPY_DATASETS)
COMPARISON_DATASETS = tuple(TMT_LFQ_DATASETS)


@dataclass(frozen=True)
class BenchmarkResults:
    """Computed metric tables written by the benchmark refresh."""

    cv_summary: pd.DataFrame
    pearson: dict[str, pd.DataFrame]
    spearman: dict[str, pd.DataFrame]
    stability: dict[str, pd.DataFrame]
    tmt_summary: pd.DataFrame
    tmt_values: dict[str, pd.DataFrame]


def method_medians(
    matrices: dict[str, dict[str, pd.DataFrame]],
    datasets: tuple[str, ...],
) -> dict[str, dict[str, pd.Series]]:
    """Return positive per-dataset protein medians for every method."""
    medians = {method: {} for method in METHODS}
    for dataset_id in datasets:
        for method in METHODS:
            values = matrices[dataset_id][method].median(axis=1, skipna=True)
            medians[method][dataset_id] = values.where(values > 0).dropna()
    return medians


def common_cv_summary(
    matrices: dict[str, dict[str, pd.DataFrame]],
) -> pd.DataFrame:
    """Describe linear-scale sample dispersion on one universe per dataset."""
    records = []
    for dataset_id, method_matrices in matrices.items():
        if any(len(matrix.columns) < 3 for matrix in method_matrices.values()):
            print(f"skipping within-dataset CV for {dataset_id}: fewer than 3 samples")
            continue
        common = set.intersection(
            *(
                set(matrix.index[matrix.notna().sum(axis=1) >= 3])
                for matrix in method_matrices.values()
            )
        )
        if not common:
            raise ValueError(f"No common CV universe for {dataset_id}")
        proteins = sorted(common)
        for method, matrix in method_matrices.items():
            values = matrix.loc[proteins]
            cv = values.std(axis=1, ddof=1) / values.mean(axis=1)
            cv = cv[np.isfinite(cv)]
            records.append(
                {
                    "dataset": dataset_id,
                    "method": method,
                    "mean_cv": cv.mean(),
                    "median_cv": cv.median(),
                    "cv_q25": cv.quantile(0.25),
                    "cv_q75": cv.quantile(0.75),
                    "n_proteins": len(cv),
                }
            )
    return pd.DataFrame.from_records(records)


def pairwise_common_proteins(
    medians: dict[str, dict[str, pd.Series]],
    left: str,
    right: str,
) -> list[str]:
    """Return one shared protein universe for every method in a dataset pair."""
    common = set.intersection(
        *(
            set(medians[method][left].index).intersection(medians[method][right].index)
            for method in METHODS
        )
    )
    return sorted(common)


def cross_experiment_correlations(
    medians: dict[str, dict[str, pd.Series]],
    datasets: tuple[str, ...],
) -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame]]:
    """Compute Pearson and Spearman matrices on pairwise common proteins."""
    pearson = {
        method: pd.DataFrame(np.nan, index=datasets, columns=datasets)
        for method in METHODS
    }
    for matrix in pearson.values():
        np.fill_diagonal(matrix.values, 1.0)
    spearman = {method: matrix.copy() for method, matrix in pearson.items()}
    for left_index, left in enumerate(datasets):
        for right in datasets[left_index + 1 :]:
            proteins = pairwise_common_proteins(medians, left, right)
            if len(proteins) < 10:
                continue
            for method in METHODS:
                left_values = np.log2(medians[method][left].loc[proteins])
                right_values = np.log2(medians[method][right].loc[proteins])
                pearson_value = left_values.corr(right_values, method="pearson")
                spearman_value = left_values.corr(right_values, method="spearman")
                pearson[method].loc[left, right] = pearson_value
                pearson[method].loc[right, left] = pearson_value
                spearman[method].loc[left, right] = spearman_value
                spearman[method].loc[right, left] = spearman_value
    return pearson, spearman


def expression_stability(
    medians: dict[str, dict[str, pd.Series]],
    datasets: tuple[str, ...],
) -> dict[str, pd.DataFrame]:
    """Compute MAD on matched dataset observations for named proteins."""
    records = {method: [] for method in METHODS}
    for protein in PROTEINS_OF_INTEREST:
        eligible = [
            dataset_id
            for dataset_id in datasets
            if all(protein in medians[method][dataset_id].index for method in METHODS)
        ]
        if len(eligible) < 2:
            continue
        for method in METHODS:
            values = np.log2(
                np.array(
                    [
                        medians[method][dataset_id].loc[protein]
                        for dataset_id in eligible
                    ]
                )
            )
            records[method].append(
                {
                    "ProteinName": protein,
                    "MAD": np.median(np.abs(values - np.median(values))),
                    "IQR": np.quantile(values, 0.75) - np.quantile(values, 0.25),
                    "Range": values.max() - values.min(),
                    "Mean": values.mean(),
                    "Std": values.std(),
                    "N_datasets": len(values),
                }
            )
    return {
        method: pd.DataFrame.from_records(method_records)
        for method, method_records in records.items()
    }


def tmt_lfq_comparison(
    matrices: dict[str, dict[str, pd.DataFrame]],
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    """Compare TMT and LFQ on one protein universe shared by all methods."""
    medians = method_medians(matrices, COMPARISON_DATASETS)
    proteins = pairwise_common_proteins(medians, "PXD007683-TMT", "PXD007683-LFQ")
    if len(proteins) < 10:
        raise ValueError("PXD007683 TMT/LFQ common universe has fewer than 10 proteins")
    summaries = []
    values = {}
    for method in METHODS:
        tmt = np.log2(medians[method]["PXD007683-TMT"].loc[proteins])
        lfq = np.log2(medians[method]["PXD007683-LFQ"].loc[proteins])
        frame = pd.DataFrame(
            {"ProteinName": proteins, "TMT": tmt.to_numpy(), "LFQ": lfq.to_numpy()}
        )
        values[method] = frame
        summaries.append(
            {
                "method": method,
                "n_proteins": len(frame),
                "pearson_r": frame["TMT"].corr(frame["LFQ"], method="pearson"),
                "spearman_r": frame["TMT"].corr(frame["LFQ"], method="spearman"),
            }
        )
    return pd.DataFrame.from_records(summaries), values


def upper_triangle_mean(matrix: pd.DataFrame) -> float:
    """Return the finite mean above a square matrix's diagonal."""
    values = matrix.to_numpy()[np.triu_indices(len(matrix), k=1)]
    return float(np.nanmean(values))


def save_results(results_dir: Path, results: BenchmarkResults) -> pd.DataFrame:
    """Write the complete tracked result contract for the refreshed benchmark."""
    results.cv_summary.to_csv(results_dir / "cv_comparison.csv", index=False)
    for method in METHODS:
        results.pearson[method].to_csv(
            results_dir / f"cross_experiment_corr_{method}.csv"
        )
        results.spearman[method].to_csv(results_dir / f"rank_consistency_{method}.csv")
        results.stability[method].to_csv(
            results_dir / f"expression_stability_{method}.csv", index=False
        )
        results.tmt_values[method].to_csv(
            results_dir / f"tmt_lfq_values_{method}.csv", index=False
        )
    results.tmt_summary.to_csv(results_dir / "tmt_lfq_comparison.csv", index=False)
    summary = []
    for method in METHODS:
        method_cv = results.cv_summary.loc[
            results.cv_summary["method"] == method, "mean_cv"
        ]
        method_stability = results.stability[method]
        summary.append(
            {
                "method": method,
                "mean_cv": method_cv.mean(),
                "mean_cross_corr": upper_triangle_mean(results.pearson[method]),
                "mean_rank_corr": upper_triangle_mean(results.spearman[method]),
                "mean_mad": method_stability["MAD"].mean(),
            }
        )
    summary_frame = pd.DataFrame.from_records(summary)
    summary_frame.to_csv(results_dir / "summary_metrics.csv", index=False)
    return summary_frame


def plot_cv(cv_summary: pd.DataFrame, figures_dir: Path) -> None:
    """Render descriptive per-dataset CV on the common protein universe."""
    frame = cv_summary.copy()
    frame["method"] = pd.Categorical(frame["method"], METHODS, ordered=True)
    figure, axis = plt.subplots(figsize=(10, 6))
    sns.boxplot(data=frame, x="method", y="mean_cv", ax=axis, color="#93C5FD")
    axis.set_xticks(range(len(METHODS)), [METHOD_LABELS[method] for method in METHODS])
    axis.set_xlabel("Quantification method")
    axis.set_ylabel("Mean within-dataset protein CV")
    axis.set_title("Within-dataset sample dispersion on a matched protein universe")
    figure.tight_layout()
    figure.savefig(figures_dir / "cv_distribution.png", dpi=180)
    plt.close(figure)


def plot_correlations(pearson: dict[str, pd.DataFrame], figures_dir: Path) -> None:
    """Render one cross-experiment Pearson heatmap per current method."""
    figure, axes = plt.subplots(2, 3, figsize=(18, 11))
    for axis, method in zip(axes.flat, METHODS):
        sns.heatmap(
            pearson[method],
            ax=axis,
            cmap="viridis",
            vmin=0.0,
            vmax=1.0,
            square=True,
            xticklabels=False,
            yticklabels=False,
            cbar_kws={"shrink": 0.7},
        )
        axis.set_title(METHOD_LABELS[method])
        axis.set_xlabel("20 human datasets")
        axis.set_ylabel("20 human datasets")
    figure.suptitle(
        "Cross-experiment Pearson correlation on pairwise matched proteins",
        fontsize=15,
    )
    figure.tight_layout()
    figure.savefig(figures_dir / "correlation_heatmap.png", dpi=180)
    plt.close(figure)
