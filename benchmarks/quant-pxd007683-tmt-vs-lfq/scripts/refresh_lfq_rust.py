#!/usr/bin/env python3
"""Recompute and plot the PXD007683 LFQ benchmark with the Rust kernel."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

import mokume

plt.switch_backend("Agg")


METHODS = ("pibaq", "maxlfq", "directlfq", "sum", "top3", "top5", "top10")
METHOD_LABELS = {
    "pibaq": "piBAQ",
    "maxlfq": "MaxLFQ",
    "directlfq": "DirectLFQ",
    "sum": "Sum",
    "top3": "Top3",
    "top5": "Top5",
    "top10": "Top10",
}
METHOD_COLORS = {
    "pibaq": "#7C3AED",
    "maxlfq": "#0F766E",
    "directlfq": "#2563EB",
    "sum": "#64748B",
    "top3": "#D97706",
    "top5": "#EA580C",
    "top10": "#DC2626",
}
CONDITION_ORDER = ("Y10", "Y5", "Y3.3")
CONDITION_COLORS = {"Y10": "#C2410C", "Y5": "#D97706", "Y3.3": "#0F766E"}
CONTRASTS = (
    ("Y10", "Y5", 1.0),
    ("Y10", "Y3.3", float(np.log2(3.0))),
    ("Y5", "Y3.3", float(np.log2(1.5))),
)


@dataclass(frozen=True)
class BenchmarkPaths:
    """Resolved inputs and outputs for one benchmark refresh."""

    feature: Path
    sdrf: Path
    fasta: Path
    matrix_dir: Path
    results_dir: Path
    figures_dir: Path


def parse_args() -> argparse.Namespace:
    """Parse explicit data paths so the script stays portable."""
    benchmark_dir = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feature", required=True, type=Path)
    parser.add_argument("--sdrf", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument(
        "--matrix-dir",
        type=Path,
        default=benchmark_dir / "data" / "current-rust" / "lfq",
    )
    parser.add_argument("--results-dir", type=Path, default=benchmark_dir / "results")
    parser.add_argument("--figures-dir", type=Path, default=benchmark_dir / "figures")
    parser.add_argument("--threads", required=True, type=int)
    return parser.parse_args()


def resolve_paths(args: argparse.Namespace) -> BenchmarkPaths:
    """Resolve and validate all input and output paths."""
    inputs = (args.feature, args.sdrf, args.fasta)
    missing = [str(path) for path in inputs if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing benchmark input(s): {', '.join(missing)}")
    if args.threads < 1:
        raise ValueError("--threads must be greater than zero")
    paths = BenchmarkPaths(
        feature=args.feature.resolve(),
        sdrf=args.sdrf.resolve(),
        fasta=args.fasta.resolve(),
        matrix_dir=args.matrix_dir.resolve(),
        results_dir=args.results_dir.resolve(),
        figures_dir=args.figures_dir.resolve(),
    )
    for output_dir in (paths.matrix_dir, paths.results_dir, paths.figures_dir):
        output_dir.mkdir(parents=True, exist_ok=True)
    return paths


def load_design(
    sdrf_path: Path,
) -> tuple[pd.DataFrame, list[str], dict[str, list[str]]]:
    """Load the 11-sample LFQ design and validate its spike-in groups."""
    design = pd.read_csv(sdrf_path, sep="\t")
    required = {"source name", "factor value[spiked compound]", "comment[data file]"}
    missing = required.difference(design.columns)
    if missing:
        raise ValueError(f"SDRF is missing columns: {sorted(missing)}")
    samples = design["source name"].astype(str).tolist()
    groups = {
        condition: design.loc[
            design["factor value[spiked compound]"] == condition, "source name"
        ]
        .astype(str)
        .tolist()
        for condition in CONDITION_ORDER
    }
    observed = {condition: len(group) for condition, group in groups.items()}
    if observed != {"Y10": 3, "Y5": 4, "Y3.3": 4}:
        raise ValueError(f"Unexpected PXD007683 LFQ group sizes: {observed}")
    return design, samples, groups


def run_quantification(paths: BenchmarkPaths, threads: int) -> None:
    """Run every protein roll-up through ``mokume.features2proteins``."""
    for method in METHODS:
        output = paths.matrix_dir / f"{method}.csv"
        options = {
            "parquet": str(paths.feature),
            "sdrf": str(paths.sdrf),
            "quant_method": method,
            "run_normalization": "none",
            "sample_normalization": "none",
            "min_aa": 7,
            "threads": threads,
            "output": str(output),
        }
        if method == "pibaq":
            options.update(
                fasta=str(paths.fasta),
                pibaq_enzyme="Trypsin",
                pibaq_max_aa=30,
            )
        else:
            options["min_unique"] = 2
        mokume.features2proteins(**options)
        if not output.is_file():
            raise RuntimeError(f"Rust kernel did not write {output}")


def parse_fasta_species(fasta_path: Path) -> dict[str, str]:
    """Map UniProt accessions to the human or yeast benchmark component."""
    species = {}
    with fasta_path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            token = line[1:].split(maxsplit=1)[0]
            parts = token.split("|")
            accession = parts[1] if len(parts) >= 2 else parts[0]
            entry_name = parts[2] if len(parts) >= 3 else token
            if "_YEAST" in entry_name:
                species[accession] = "yeast"
            elif "_HUMAN" in entry_name:
                species[accession] = "human"
    return species


def load_matrices(matrix_dir: Path, samples: list[str]) -> dict[str, pd.DataFrame]:
    """Load Rust matrices and represent non-positive cells as missing evidence."""
    matrices = {}
    expected_columns = ["ProteinName", *samples]
    for method in METHODS:
        path = matrix_dir / f"{method}.csv"
        matrix = pd.read_csv(path)
        if matrix.columns.tolist() != expected_columns:
            raise ValueError(f"Unexpected columns in {path}: {matrix.columns.tolist()}")
        if matrix["ProteinName"].duplicated().any():
            raise ValueError(f"Duplicate proteins in {path}")
        values = matrix.set_index("ProteinName").apply(pd.to_numeric, errors="coerce")
        matrices[method] = values.where(values > 0)
    return matrices


def classified_proteins(
    matrices: dict[str, pd.DataFrame], species: dict[str, str]
) -> set[str]:
    """Return proteins shared by all method outputs and assigned to a species."""
    common = set.intersection(*(set(matrix.index) for matrix in matrices.values()))
    return {protein for protein in common if species.get(protein) in {"human", "yeast"}}


def common_cv_values(
    matrices: dict[str, pd.DataFrame],
    groups: dict[str, list[str]],
    proteins: set[str],
) -> pd.DataFrame:
    """Compute CV on the same complete protein-condition universe for every method."""
    records = []
    for condition in CONDITION_ORDER:
        columns = groups[condition]
        complete = set(proteins)
        for matrix in matrices.values():
            complete.intersection_update(
                matrix.index[matrix[columns].notna().all(axis=1)]
            )
        ordered = sorted(complete)
        for method, matrix in matrices.items():
            values = matrix.loc[ordered, columns]
            cv = values.std(axis=1, ddof=1) / values.mean(axis=1)
            records.extend(
                {
                    "method": method,
                    "condition": condition,
                    "protein": protein,
                    "cv": value,
                }
                for protein, value in cv.items()
                if np.isfinite(value)
            )
    return pd.DataFrame.from_records(records)


def contrast_values(
    matrix: pd.DataFrame,
    group_a: list[str],
    group_b: list[str],
) -> tuple[pd.Series, pd.Series]:
    """Return log2 fold changes and the two-replicate evidence mask."""
    values_a = np.log2(matrix[group_a])
    values_b = np.log2(matrix[group_b])
    fold_change = values_a.median(axis=1) - values_b.median(axis=1)
    enough = (values_a.notna().sum(axis=1) >= 2) & (values_b.notna().sum(axis=1) >= 2)
    return fold_change, enough


def measurable_fold_changes(
    matrices: dict[str, pd.DataFrame],
    group_a: list[str],
    group_b: list[str],
    proteins: set[str],
) -> tuple[dict[str, pd.Series], set[str]]:
    """Return method fold changes and their common measurable proteins."""
    fold_changes = {}
    common = set(proteins)
    for method, matrix in matrices.items():
        fold_change, enough = contrast_values(matrix, group_a, group_b)
        fold_changes[method] = fold_change
        common.intersection_update(matrix.index[enough])
    return fold_changes, common


def organism_fold_change_records(
    fold_changes: dict[str, pd.Series],
    proteins: set[str],
    species: dict[str, str],
    contrast: str,
    organism_truth: tuple[str, float],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Build detailed and summary records for one organism and contrast."""
    organism, expected = organism_truth
    organism_proteins = sorted(
        protein for protein in proteins if species.get(protein) == organism
    )
    value_records = []
    summary_records = []
    for method, fold_change in fold_changes.items():
        observed = fold_change.loc[organism_proteins].dropna()
        error = observed - expected
        value_records.extend(
            {
                "method": method,
                "contrast": contrast,
                "species": organism,
                "expected_log2fc": expected,
                "protein": protein,
                "observed_log2fc": value,
            }
            for protein, value in observed.items()
        )
        summary_records.append(
            {
                "method": method,
                "contrast": contrast,
                "species": organism,
                "expected_log2fc": expected,
                "n_proteins": len(observed),
                "median_log2fc": observed.median(),
                "median_bias": observed.median() - expected,
                "rmse": float(np.sqrt(np.mean(np.square(error)))),
            }
        )
    return value_records, summary_records


def common_fold_changes(
    matrices: dict[str, pd.DataFrame],
    groups: dict[str, list[str]],
    proteins: set[str],
    species: dict[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute fold changes on a common measurable universe per contrast."""
    values_records = []
    summary_records = []
    for group_a, group_b, yeast_expected in CONTRASTS:
        fold_changes, common = measurable_fold_changes(
            matrices, groups[group_a], groups[group_b], proteins
        )
        contrast = f"{group_a}/{group_b}"
        for organism_truth in (("yeast", yeast_expected), ("human", 0.0)):
            records = organism_fold_change_records(
                fold_changes,
                common,
                species,
                contrast,
                organism_truth,
            )
            values_records.extend(records[0])
            summary_records.extend(records[1])
    return (
        pd.DataFrame.from_records(values_records),
        pd.DataFrame.from_records(summary_records),
    )


def method_metrics(
    matrices: dict[str, pd.DataFrame],
    species: dict[str, str],
    cv_values: pd.DataFrame,
) -> pd.DataFrame:
    """Summarize coverage, completeness, and common-universe CV."""
    records = []
    for method, matrix in matrices.items():
        assigned = matrix.index.map(species.get).isin(["human", "yeast"])
        values = matrix.loc[assigned]
        method_cv = cv_values.loc[cv_values["method"] == method, "cv"]
        records.append(
            {
                "method": method,
                "n_proteins": int(values.notna().any(axis=1).sum()),
                "matrix_completeness": float(values.notna().mean().mean()),
                "common_cv_n": len(method_cv),
                "common_cv_median": method_cv.median(),
                "common_cv_q75": method_cv.quantile(0.75),
            }
        )
    return pd.DataFrame.from_records(records)


def peptide_observation_counts(
    feature_path: Path, design: pd.DataFrame
) -> pd.DataFrame:
    """Count observed and absent canonical peptide sequences per LFQ sample."""
    features = pd.read_parquet(feature_path, columns=["sequence", "run_file_name"])
    features = features.dropna(subset=["sequence", "run_file_name"])
    total_peptides = int(features["sequence"].nunique())
    counts = features.groupby("run_file_name")["sequence"].nunique()
    run_to_sample = {
        Path(str(data_file)).stem: str(sample)
        for data_file, sample in zip(
            design["comment[data file]"], design["source name"]
        )
    }
    records = []
    for run, count in counts.items():
        sample = run_to_sample.get(Path(str(run)).stem)
        if sample is None:
            raise ValueError(f"Feature run {run!r} is absent from the SDRF")
        condition = design.loc[
            design["source name"] == sample, "factor value[spiked compound]"
        ]
        records.append(
            {
                "sample": sample,
                "condition": str(condition.iloc[0]),
                "observed_peptides": int(count),
                "missing_peptides": total_peptides - int(count),
                "missing_rate": 1.0 - int(count) / total_peptides,
            }
        )
    result = pd.DataFrame.from_records(records).sort_values("sample")
    if len(result) != len(design):
        raise ValueError(
            f"Matched {len(result)} feature runs to {len(design)} SDRF samples"
        )
    return result


def style_axis(axis: plt.Axes) -> None:
    """Apply the shared benchmark figure style."""
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(axis="y", color="#CBD5E1", alpha=0.45, linewidth=0.7)
    axis.set_axisbelow(True)


def plot_mean_cv(metrics: pd.DataFrame, output: Path) -> None:
    """Plot median CV using the common complete protein-condition universe."""
    ordered = metrics.sort_values("common_cv_median")
    labels = [METHOD_LABELS[method] for method in ordered["method"]]
    colors = [METHOD_COLORS[method] for method in ordered["method"]]
    figure, axis = plt.subplots(figsize=(9.2, 5.2))
    bars = axis.barh(labels, ordered["common_cv_median"] * 100, color=colors)
    axis.bar_label(bars, fmt="%.1f%%", padding=4, fontsize=9)
    axis.set_xlabel("Median within-condition CV (%)")
    axis.set_title(
        "PXD007683 LFQ — current Rust quantification",
        loc="left",
        weight="bold",
        y=1.075,
    )
    axis.text(
        0,
        1.015,
        "Same complete protein-condition universe across all methods; "
        "no normalization or imputation",
        transform=axis.transAxes,
        color="#475569",
        fontsize=9,
    )
    style_axis(axis)
    figure.tight_layout()
    figure.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(figure)


def plot_cv_distribution(cv_values: pd.DataFrame, output: Path) -> None:
    """Plot common-universe CV distributions for all Rust methods."""
    data = [
        cv_values.loc[cv_values["method"] == method, "cv"] * 100 for method in METHODS
    ]
    figure, axis = plt.subplots(figsize=(11.5, 5.6))
    boxes = axis.boxplot(
        data,
        tick_labels=[METHOD_LABELS[m] for m in METHODS],
        showfliers=False,
        patch_artist=True,
    )
    for patch, method in zip(boxes["boxes"], METHODS):
        patch.set_facecolor(METHOD_COLORS[method])
        patch.set_edgecolor(METHOD_COLORS[method])
        patch.set_alpha(0.3)
        patch.set_linewidth(2)
    for median in boxes["medians"]:
        median.set_color("#0F172A")
        median.set_linewidth(1.4)
    axis.set_ylabel("Within-condition CV (%)")
    axis.set_title("Protein-condition CV distributions", loc="left", weight="bold")
    axis.set_ylim(0, float(cv_values["cv"].quantile(0.95) * 100 * 1.05))
    axis.tick_params(axis="x", rotation=20)
    style_axis(axis)
    figure.tight_layout()
    figure.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(figure)


def plot_fold_changes(values: pd.DataFrame, output: Path) -> None:
    """Plot yeast fold-change distributions against the declared spike-in truth."""
    yeast = values.loc[values["species"] == "yeast"]
    figure, axes = plt.subplots(1, 3, figsize=(16, 5.4), sharey=True)
    for axis, (group_a, group_b, expected) in zip(axes, CONTRASTS):
        contrast = f"{group_a}/{group_b}"
        distributions = [
            yeast.loc[
                (yeast["method"] == method) & (yeast["contrast"] == contrast),
                "observed_log2fc",
            ]
            for method in METHODS
        ]
        violin = axis.violinplot(distributions, showmedians=True, showextrema=False)
        for body, method in zip(violin["bodies"], METHODS):
            body.set_facecolor(METHOD_COLORS[method])
            body.set_edgecolor("none")
            body.set_alpha(0.72)
        violin["cmedians"].set_color("#0F172A")
        violin["cmedians"].set_linewidth(1.5)
        axis.axhline(expected, color="#0F172A", linestyle="--", linewidth=1.2)
        axis.set_title(f"{contrast}  truth={expected:.3f}", weight="bold")
        axis.set_xticks(range(1, len(METHODS) + 1), [METHOD_LABELS[m] for m in METHODS])
        axis.tick_params(axis="x", rotation=45)
        style_axis(axis)
    axes[0].set_ylabel("Observed yeast log2 fold change")
    figure.suptitle(
        "PXD007683 LFQ spike-in recovery — common measurable yeast proteins",
        x=0.06,
        ha="left",
        weight="bold",
    )
    figure.tight_layout()
    figure.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(figure)


def plot_missing_peptides(counts: pd.DataFrame, output: Path) -> None:
    """Plot observed versus absent peptides in the across-sample union."""
    figure, axis = plt.subplots(figsize=(11, 5.5))
    positions = np.arange(len(counts))
    colors = [CONDITION_COLORS[condition] for condition in counts["condition"]]
    axis.bar(positions, counts["observed_peptides"], color=colors)
    axis.bar(
        positions,
        counts["missing_peptides"],
        bottom=counts["observed_peptides"],
        color="#E2E8F0",
    )
    axis.set_xticks(
        positions, counts["sample"].str.replace("HYE LFQ ", "S", regex=False)
    )
    axis.set_ylabel("Unique canonical peptides")
    axis.set_title(
        "Peptide observation completeness by LFQ sample", loc="left", weight="bold"
    )
    handles = [
        *(
            Patch(color=CONDITION_COLORS[condition], label=condition)
            for condition in CONDITION_ORDER
        ),
        Patch(color="#E2E8F0", label="Absent from sample"),
    ]
    axis.legend(handles=handles, frameon=False, ncol=4)
    style_axis(axis)
    figure.tight_layout()
    figure.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    """Run quantification, metrics, and plots as one reproducible refresh."""
    args = parse_args()
    paths = resolve_paths(args)
    design, samples, groups = load_design(paths.sdrf)
    run_quantification(paths, args.threads)
    species = parse_fasta_species(paths.fasta)
    matrices = load_matrices(paths.matrix_dir, samples)
    proteins = classified_proteins(matrices, species)
    cv_values = common_cv_values(matrices, groups, proteins)
    fold_values, fold_summary = common_fold_changes(matrices, groups, proteins, species)
    metrics = method_metrics(matrices, species, cv_values)
    peptide_counts = peptide_observation_counts(paths.feature, design)

    metrics.to_csv(
        paths.results_dir / "pxd007683_lfq_rust_method_metrics.csv", index=False
    )
    fold_summary.to_csv(
        paths.results_dir / "pxd007683_lfq_rust_fold_change.csv", index=False
    )
    peptide_counts.to_csv(
        paths.results_dir / "pxd007683_lfq_peptide_observations.csv", index=False
    )

    plot_mean_cv(metrics, paths.figures_dir / "method_mean_cv_lfq.png")
    plot_cv_distribution(cv_values, paths.figures_dir / "method_per_p_cv_lfq.png")
    plot_fold_changes(fold_values, paths.figures_dir / "fold_change_lfq.png")
    plot_missing_peptides(
        peptide_counts, paths.figures_dir / "missing_peptides_by_sample.png"
    )
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
