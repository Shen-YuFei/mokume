#!/usr/bin/env python3
"""Recompute the HeLa/human quantification benchmark with the Rust kernel."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

from config import (
    IBAQPY_DATASETS,
    PIBAQ_PARAMS,
    TMT_LFQ_DATASETS,
)
from hela_report import (
    BenchmarkResults,
    COMPARISON_DATASETS,
    MAIN_DATASETS,
    METHODS,
    common_cv_summary,
    cross_experiment_correlations,
    expression_stability,
    method_medians,
    plot_correlations,
    plot_cv,
    save_results,
    tmt_lfq_comparison,
)
import mokume

ALL_DATASETS = {**IBAQPY_DATASETS, **TMT_LFQ_DATASETS}


@dataclass(frozen=True)
class RefreshPaths:
    """Resolved benchmark inputs and temporary/current output directories."""

    raw_dir: Path
    fasta: Path
    work_dir: Path
    results_dir: Path
    figures_dir: Path


def parse_args() -> argparse.Namespace:
    """Parse explicit paths so validation can precede every tracked overwrite."""
    benchmark_dir = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    for name in ("results", "figures"):
        parser.add_argument(f"--{name}-dir", type=Path, default=benchmark_dir / name)
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
        "quant_method": method,
        "output": str(output),
    }
    if method in {"directlfq", "maxlfq"}:
        options["threads"] = threads
    if method == "pibaq":
        options.update(
            fasta=str(paths.fasta),
            enzyme=ALL_DATASETS[dataset_id].enzyme,
            min_aa=PIBAQ_PARAMS["min_aa"],
            max_aa=PIBAQ_PARAMS["max_aa"],
        )
    elif method == "directlfq":
        options["directlfq_min_nonan"] = 1
    mokume.peptides2protein(**options)
    return output


def load_matrix(path: Path, method: str) -> pd.DataFrame:
    """Load one positive linear-intensity result as protein x sample."""
    result = pd.read_csv(path, sep="\t")
    value_column = "PiBAQ" if method == "pibaq" else "Intensity"
    columns = ["ProteinName", "SampleID", value_column]
    missing = set(columns).difference(result.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    quantities = pd.to_numeric(result[value_column], errors="coerce")
    positive = np.isfinite(quantities) & quantities.gt(0)
    result = result.loc[positive, columns].copy()
    result[value_column] = quantities.loc[positive]
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
