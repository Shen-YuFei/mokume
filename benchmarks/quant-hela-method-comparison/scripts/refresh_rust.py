#!/usr/bin/env python3
"""Recompute the HeLa/human quantification benchmark with the Rust kernel."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import duckdb
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from config import (
    IBAQPY_DATASETS,
    PIBAQ_PARAMS,
    PROTEINS_OF_INTEREST,
    QUANTIFICATION_METHODS,
    TMT_LFQ_DATASETS,
)
import mokume

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
ALL_DATASETS = {**IBAQPY_DATASETS, **TMT_LFQ_DATASETS}


@dataclass(frozen=True)
class RefreshPaths:
    """Resolved benchmark inputs and temporary/current output directories."""

    raw_dir: Path
    fasta: Path
    work_dir: Path
    results_dir: Path
    figures_dir: Path


@dataclass(frozen=True)
class BenchmarkResults:
    """Computed metric tables written by the benchmark refresh."""

    cv_summary: pd.DataFrame
    pearson: dict[str, pd.DataFrame]
    spearman: dict[str, pd.DataFrame]
    stability: dict[str, pd.DataFrame]
    tmt_summary: pd.DataFrame
    tmt_values: dict[str, pd.DataFrame]


def parse_args() -> argparse.Namespace:
    """Parse explicit paths so validation can precede every tracked overwrite."""
    benchmark_dir = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--results-dir", type=Path, default=benchmark_dir / "results")
    parser.add_argument("--figures-dir", type=Path, default=benchmark_dir / "figures")
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def input_candidates(raw_dir: Path, dataset_id: str) -> tuple[Path, ...]:
    """Return current and legacy local names for one public source file."""
    dataset = ALL_DATASETS[dataset_id]
    if dataset.file_format == "msstats_csv":
        return (
            raw_dir / dataset.feature_file,
            raw_dir / f"{dataset_id}_feature.csv",
        )
    return (
        raw_dir / f"{dataset_id}_feature.parquet",
        raw_dir / dataset.feature_file,
    )


def resolve_paths(args: argparse.Namespace) -> tuple[RefreshPaths, dict[str, Path]]:
    """Validate every input before creating any benchmark output."""
    if args.threads < 1:
        raise ValueError("--threads must be greater than zero")
    if not args.fasta.is_file():
        raise FileNotFoundError(f"Missing FASTA: {args.fasta}")
    sources = {}
    missing = []
    for dataset_id in ALL_DATASETS:
        source = next(
            (
                path
                for path in input_candidates(args.raw_dir, dataset_id)
                if path.is_file()
            ),
            None,
        )
        if source is None:
            missing.append(dataset_id)
        else:
            sources[dataset_id] = source.resolve()
    if missing:
        raise FileNotFoundError(f"Missing benchmark inputs: {', '.join(missing)}")
    paths = RefreshPaths(
        raw_dir=args.raw_dir.resolve(),
        fasta=args.fasta.resolve(),
        work_dir=args.work_dir.resolve(),
        results_dir=args.results_dir.resolve(),
        figures_dir=args.figures_dir.resolve(),
    )
    for output_dir in (paths.work_dir, paths.results_dir, paths.figures_dir):
        output_dir.mkdir(parents=True, exist_ok=True)
    (paths.work_dir / "duckdb-tmp").mkdir(parents=True, exist_ok=True)
    return paths, sources


def prepare_parquet_source(
    connection: duckdb.DuckDBPyConnection,
    source: Path,
    output: Path,
) -> None:
    """Stream a QPX feature parquet into canonical peptide-level rows."""
    relation = connection.sql(
        """
        WITH expanded AS (
            SELECT pg_accessions[1] AS protein,
                   sequence AS peptide,
                   observation.sample_accession AS sample,
                   CAST(observation.intensity AS DOUBLE) AS intensity
            FROM read_parquet(?)
            CROSS JOIN UNNEST(intensities) AS nested(observation)
        )
        SELECT protein AS ProteinName,
               peptide AS PeptideSequence,
               sample AS SampleID,
               SUM(intensity) AS NormIntensity
        FROM expanded
        WHERE protein IS NOT NULL AND protein <> ''
          AND peptide IS NOT NULL AND peptide <> ''
          AND sample IS NOT NULL AND sample <> ''
          AND isfinite(intensity) AND intensity > 0
          AND NOT starts_with(protein, 'DECOY_')
          AND NOT starts_with(protein, 'CONTAM_')
          AND NOT starts_with(protein, 'CON__')
          AND NOT starts_with(protein, 'REV__')
        GROUP BY protein, peptide, sample
        """,
        params=[str(source)],
    )
    relation.write_parquet(str(output), compression="zstd")


def prepare_msstats_source(
    connection: duckdb.DuckDBPyConnection,
    source: Path,
    output: Path,
) -> None:
    """Collapse public MSstats rows to canonical peptide/sample intensities."""
    relation = connection.sql(
        """
        WITH source_rows AS (
            SELECT CAST(ProteinName AS VARCHAR) AS protein_raw,
                   regexp_replace(
                       CAST(PeptideSequence AS VARCHAR),
                       '\\([^)]*\\)|\\[[^]]*\\]|[.\\-]',
                       '',
                       'g'
                   ) AS peptide,
                   CAST(BioReplicate AS VARCHAR) AS sample,
                   CAST(Condition AS VARCHAR) AS condition,
                   CAST(Intensity AS DOUBLE) AS intensity
            FROM read_csv_auto(?, header = true, sample_size = -1)
        ), canonical AS (
            SELECT CASE
                       WHEN contains(split_part(protein_raw, ';', 1), '|')
                       THEN split_part(split_part(protein_raw, ';', 1), '|', 2)
                       ELSE split_part(protein_raw, ';', 1)
                   END AS protein,
                   peptide,
                   sample,
                   condition,
                   intensity
            FROM source_rows
        )
        SELECT protein AS ProteinName,
               peptide AS PeptideSequence,
               sample AS SampleID,
               SUM(intensity) AS NormIntensity,
               condition AS Condition
        FROM canonical
        WHERE protein IS NOT NULL AND protein <> ''
          AND peptide IS NOT NULL AND peptide <> ''
          AND sample IS NOT NULL AND sample <> ''
          AND isfinite(intensity) AND intensity > 0
          AND NOT starts_with(protein, 'DECOY_')
          AND NOT starts_with(protein, 'CONTAM_')
          AND NOT starts_with(protein, 'CON__')
          AND NOT starts_with(protein, 'REV__')
        GROUP BY protein, peptide, sample, condition
        """,
        params=[str(source)],
    )
    relation.write_parquet(str(output), compression="zstd")


def validate_peptide_table(
    connection: duckdb.DuckDBPyConnection,
    path: Path,
    dataset_id: str,
) -> dict[str, int]:
    """Fail on an empty or malformed prepared peptide table."""
    row = connection.execute(
        """
        SELECT COUNT(*), COUNT(DISTINCT ProteinName),
               COUNT(DISTINCT PeptideSequence), COUNT(DISTINCT SampleID),
               COUNT(*) FILTER (
                   WHERE NOT isfinite(NormIntensity) OR NormIntensity <= 0
               )
        FROM read_parquet(?)
        """,
        [str(path)],
    ).fetchone()
    rows, proteins, peptides, samples, invalid = map(int, row)
    if rows == 0 or proteins == 0 or peptides == 0 or samples < 1 or invalid:
        raise ValueError(
            f"Invalid prepared table for {dataset_id}: "
            f"rows={rows}, proteins={proteins}, peptides={peptides}, "
            f"samples={samples}, invalid={invalid}"
        )
    if dataset_id in COMPARISON_DATASETS and samples != 11:
        raise ValueError(f"{dataset_id} has {samples} samples; expected 11")
    return {
        "rows": rows,
        "proteins": proteins,
        "peptides": peptides,
        "samples": samples,
    }


def prepare_inputs(
    paths: RefreshPaths,
    sources: dict[str, Path],
    threads: int,
    force: bool,
) -> dict[str, Path]:
    """Prepare every source with DuckDB using bounded, streaming aggregation."""
    peptide_dir = paths.work_dir / "peptides"
    peptide_dir.mkdir(parents=True, exist_ok=True)
    connection = duckdb.connect()
    prepared = {}
    try:
        connection.execute("SET threads=?", [threads])
        connection.execute("SET memory_limit='80GB'")
        connection.execute("SET temp_directory=?", [str(paths.work_dir / "duckdb-tmp")])
        for dataset_id, source in sources.items():
            output = peptide_dir / f"{dataset_id}.parquet"
            if force and output.exists():
                output.unlink()
            if not output.exists():
                if ALL_DATASETS[dataset_id].file_format == "msstats_csv":
                    prepare_msstats_source(connection, source, output)
                else:
                    prepare_parquet_source(connection, source, output)
            stats = validate_peptide_table(connection, output, dataset_id)
            print(
                f"prepared {dataset_id}: {stats['rows']:,} rows, "
                f"{stats['proteins']:,} proteins, "
                f"{stats['peptides']:,} peptides, {stats['samples']} samples"
            )
            prepared[dataset_id] = output
    finally:
        connection.close()
    return prepared


def run_method(
    paths: RefreshPaths,
    dataset_id: str,
    peptide_path: Path,
    method: str,
    threads: int,
) -> Path:
    """Run one current quantification method through the wheel's Rust kernel."""
    output = paths.work_dir / "protein-quant" / dataset_id / f"{method}.tsv"
    options = {
        "peptides": str(peptide_path),
        "method": method,
        "output": str(output),
        "threads": threads,
    }
    if method == "pibaq":
        options.update(
            fasta=str(paths.fasta),
            enzyme=ALL_DATASETS[dataset_id].enzyme,
            min_aa=PIBAQ_PARAMS["min_aa"],
            max_aa=PIBAQ_PARAMS["max_aa"],
        )
    elif method == "directlfq":
        options["min_nonan"] = 1
    mokume.peptides2protein(**options)
    return output


def load_matrix(path: Path, method: str) -> pd.DataFrame:
    """Load one positive linear-intensity result as protein x sample."""
    result = pd.read_csv(path, sep="\t")
    value_column = "PiBAQ" if method == "pibaq" else "Intensity"
    required = {"ProteinName", "SampleID", value_column}
    missing = required.difference(result.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    result[value_column] = pd.to_numeric(result[value_column], errors="coerce")
    result = result.loc[
        np.isfinite(result[value_column]) & (result[value_column] > 0),
        ["ProteinName", "SampleID", value_column],
    ].copy()
    if result.empty:
        raise ValueError(f"{path} contains no positive finite quantities")
    if result.duplicated(["ProteinName", "SampleID"]).any():
        raise ValueError(f"{path} contains duplicate protein/sample cells")
    matrix = result.pivot(
        index="ProteinName", columns="SampleID", values=value_column
    ).sort_index()
    matrix.columns = matrix.columns.astype(str)
    return matrix.reindex(sorted(matrix.columns), axis=1)


def quantify_all(
    paths: RefreshPaths,
    prepared: dict[str, Path],
    threads: int,
    force: bool,
) -> dict[str, dict[str, pd.DataFrame]]:
    """Run and validate all six methods on every prepared dataset."""
    matrix_root = paths.work_dir / "protein-quant"
    matrices = {}
    for dataset_id, peptide_path in prepared.items():
        dataset_dir = matrix_root / dataset_id
        dataset_dir.mkdir(parents=True, exist_ok=True)
        method_matrices = {}
        for method in METHODS:
            output = dataset_dir / f"{method}.tsv"
            if force or not output.exists():
                run_method(
                    paths,
                    dataset_id,
                    peptide_path,
                    method,
                    threads,
                )
            matrix = load_matrix(output, method)
            print(
                f"quantified {dataset_id}/{method}: "
                f"{len(matrix):,} proteins x {len(matrix.columns)} samples"
            )
            method_matrices[method] = matrix
        matrices[dataset_id] = method_matrices
    return matrices


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


def save_results(
    results_dir: Path,
    results: BenchmarkResults,
) -> pd.DataFrame:
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


def main() -> None:
    """Run preparation, Rust quantification, matched metrics, and plotting."""
    args = parse_args()
    paths, sources = resolve_paths(args)
    prepared = prepare_inputs(paths, sources, args.threads, args.force)
    matrices = quantify_all(paths, prepared, args.threads, args.force)
    cv_summary = common_cv_summary(matrices)
    medians = method_medians(matrices, MAIN_DATASETS)
    pearson, spearman = cross_experiment_correlations(medians, MAIN_DATASETS)
    stability = expression_stability(medians, MAIN_DATASETS)
    tmt_summary, tmt_values = tmt_lfq_comparison(matrices)
    summary = save_results(
        paths.results_dir,
        BenchmarkResults(
            cv_summary=cv_summary,
            pearson=pearson,
            spearman=spearman,
            stability=stability,
            tmt_summary=tmt_summary,
            tmt_values=tmt_values,
        ),
    )
    plot_cv(cv_summary, paths.figures_dir)
    plot_correlations(pearson, paths.figures_dir)
    print("\nRefreshed summary")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
