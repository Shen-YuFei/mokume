#!/usr/bin/env python3
"""Recompute the Quartet batch benchmark with the current Rust kernel."""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import duckdb
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

import mokume


plt.switch_backend("Agg")

METHODS = ("pibaq", "maxlfq", "directlfq", "top3")
METHOD_LABELS = {
    "pibaq": "piBAQ",
    "maxlfq": "MaxLFQ",
    "directlfq": "DirectLFQ",
    "top3": "Top3",
}
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


@dataclass(frozen=True)
class RefreshPaths:
    """Resolved sources and isolated output locations."""

    source_dir: Path
    fasta: Path
    work_dir: Path
    results_dir: Path
    figures_dir: Path


@dataclass(frozen=True)
class QuartetAnalysis:
    """Raw, complete, and matched matrices used by the tracked outputs."""

    raw_full: dict[str, pd.DataFrame]
    corrected_complete: dict[str, pd.DataFrame]
    raw_matched: dict[str, pd.DataFrame]
    corrected_matched: dict[str, pd.DataFrame]
    coverage: dict[str, int]
    complete_coverage: dict[str, int]


def parse_args() -> argparse.Namespace:
    """Parse explicit paths so tracked outputs need not be touched during a run."""
    benchmark_dir = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--results-dir", type=Path, default=benchmark_dir / "results")
    parser.add_argument("--figures-dir", type=Path, default=benchmark_dir / "figures")
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def resolve_paths(args: argparse.Namespace) -> tuple[RefreshPaths, dict[str, Path]]:
    """Validate the six source files and FASTA before creating outputs."""
    if args.threads < 1:
        raise ValueError("--threads must be greater than zero")
    if not args.fasta.is_file():
        raise FileNotFoundError(f"Missing FASTA: {args.fasta}")
    source_dir = args.source_dir.resolve()
    evidence = {}
    missing = []
    for folder in ("APT_DDA", "APT_DIA", "BGI_DIA", "FDU_DDA", "FDU_DIA", "NVG_DDA"):
        path = source_dir / folder / "evidence.txt"
        if path.is_file():
            evidence[folder] = path
        else:
            missing.append(folder)
    if missing:
        raise FileNotFoundError(f"Missing Quartet evidence files: {', '.join(missing)}")
    paths = RefreshPaths(
        source_dir=source_dir,
        fasta=args.fasta.resolve(),
        work_dir=args.work_dir.resolve(),
        results_dir=args.results_dir.resolve(),
        figures_dir=args.figures_dir.resolve(),
    )
    for output_dir in (paths.work_dir, paths.results_dir, paths.figures_dir):
        output_dir.mkdir(parents=True, exist_ok=True)
    (paths.work_dir / "duckdb-tmp").mkdir(parents=True, exist_ok=True)
    return paths, evidence


def batch_name(folder: str) -> str:
    """Convert source folder names to the upstream mode-lab batch convention."""
    lab, mode = folder.split("_", maxsplit=1)
    return f"{mode}_{lab}"


def ordered_runs(raw_files: list[str]) -> list[tuple[str, str, int]]:
    """Map fixed injection order, while honoring explicit Quartet labels."""
    explicit = []
    for raw_file in raw_files:
        match = re.search(r"_(D5|D6|F7|M8)_([123])$", raw_file)
        if match:
            explicit.append((raw_file, match.group(1), int(match.group(2))))
    if explicit:
        if len(explicit) != 12:
            raise ValueError("Only part of one batch has explicit Quartet labels")
        order = {sample: position for position, sample in enumerate(SAMPLE_TYPES)}
        return sorted(explicit, key=lambda item: (item[2], order[item[1]]))

    numbered = []
    for raw_file in raw_files:
        if raw_file.startswith("Exp"):
            match = re.match(r"Exp(\d+)", raw_file)
        else:
            match = re.search(r"(\d+)$", raw_file)
        if match is None:
            raise ValueError(f"Cannot resolve Quartet injection order: {raw_file}")
        numbered.append((int(match.group(1)), raw_file))
    ordered = [raw_file for _, raw_file in sorted(numbered)]
    return [
        (raw_file, SAMPLE_TYPES[index % 4], index // 4 + 1)
        for index, raw_file in enumerate(ordered)
    ]


def read_raw_files(connection: duckdb.DuckDBPyConnection, path: Path) -> list[str]:
    """Read the 12 run identifiers from one MaxQuant evidence file."""
    rows = connection.execute(
        """
        SELECT DISTINCT "Raw file"
        FROM read_csv(?, delim='\t', header=true, all_varchar=true)
        WHERE "Raw file" IS NOT NULL AND "Raw file" <> ''
        """,
        [str(path)],
    ).fetchall()
    raw_files = [str(row[0]) for row in rows]
    if len(raw_files) != 12:
        raise ValueError(f"{path} has {len(raw_files)} runs; expected 12")
    return raw_files


def build_metadata(
    connection: duckdb.DuckDBPyConnection,
    evidence: dict[str, Path],
) -> pd.DataFrame:
    """Build the balanced 6-batch, 72-sample design from source run names."""
    records = []
    for folder, path in evidence.items():
        lab, mode = folder.split("_", maxsplit=1)
        batch = batch_name(folder)
        for order, (raw_file, sample, tube) in enumerate(
            ordered_runs(read_raw_files(connection, path)), start=1
        ):
            records.append(
                {
                    "run_id": f"{batch.replace('_', '')}-{sample}-{tube}",
                    "raw_file": raw_file,
                    "mode": mode,
                    "lab": lab,
                    "batch": batch,
                    "sample": sample,
                    "tube": tube,
                    "order": order,
                    "source_folder": folder,
                }
            )
    metadata = pd.DataFrame.from_records(records)
    if len(metadata) != 72 or metadata["run_id"].nunique() != 72:
        raise ValueError("Quartet metadata is not a unique 72-sample design")
    observed_batches = tuple(metadata["batch"].drop_duplicates())
    if observed_batches != EXPECTED_BATCHES:
        raise ValueError(f"Unexpected batch order: {observed_batches}")
    group_sizes = metadata.groupby(["batch", "sample"], observed=True).size()
    if not (group_sizes == 3).all():
        raise ValueError("Every Quartet batch/sample group must have three replicates")
    return metadata


def prepare_batch(
    connection: duckdb.DuckDBPyConnection,
    evidence_path: Path,
    mapping: pd.DataFrame,
    output: Path,
) -> None:
    """Stream one evidence file into canonical positive peptide quantities."""
    connection.register("run_mapping", mapping[["raw_file", "run_id"]])
    relation = connection.sql(
        """
        SELECT CAST(source.Proteins AS VARCHAR) AS ProteinName,
               CAST(source.Sequence AS VARCHAR) AS PeptideSequence,
               mapping.run_id AS SampleID,
               SUM(TRY_CAST(source.Intensity AS DOUBLE)) AS NormIntensity
        FROM read_csv(?, delim='\t', header=true, all_varchar=true) AS source
        INNER JOIN run_mapping AS mapping
            ON source."Raw file" = mapping.raw_file
        WHERE TRY_CAST(source.Intensity AS DOUBLE) > 0
          AND source.Proteins IS NOT NULL AND source.Proteins <> ''
          AND NOT contains(source.Proteins, ';')
          AND source.Sequence IS NOT NULL AND source.Sequence <> ''
          AND COALESCE(source.Reverse, '') <> '+'
          AND COALESCE(source."Potential contaminant", '') <> '+'
          AND NOT starts_with(source.Proteins, 'REV__')
          AND NOT starts_with(source.Proteins, 'CON__')
        GROUP BY source.Proteins, source.Sequence, mapping.run_id
        """,
        params=[str(evidence_path)],
    )
    relation.write_parquet(str(output), compression="zstd")
    connection.unregister("run_mapping")


def prepare_batch_files(
    connection: duckdb.DuckDBPyConnection,
    prepared_dir: Path,
    evidence: dict[str, Path],
    metadata: pd.DataFrame,
    force: bool,
) -> list[Path]:
    """Prepare each source batch and return its parquet path."""
    batch_paths = []
    for folder, source in evidence.items():
        output = prepared_dir / f"{folder}.parquet"
        if force and output.exists():
            output.unlink()
        if not output.exists():
            mapping = metadata.loc[metadata["source_folder"] == folder]
            prepare_batch(connection, source, mapping, output)
        batch_paths.append(output)
    return batch_paths


def combine_batches(
    connection: duckdb.DuckDBPyConnection,
    batch_paths: list[Path],
    combined: Path,
    force: bool,
) -> None:
    """Combine prepared batches without interpolating paths into SQL."""
    if force and combined.exists():
        combined.unlink()
    if combined.exists():
        return
    relation = connection.read_parquet([str(path) for path in batch_paths]).project(
        "ProteinName, PeptideSequence, SampleID, NormIntensity"
    )
    relation.write_parquet(str(combined), compression="zstd")


def peptide_table_stats(
    connection: duckdb.DuckDBPyConnection, combined: Path
) -> tuple[int, int, int, int, int]:
    """Return row, protein, peptide, sample, and invalid-value counts."""
    stats = connection.execute(
        """
        SELECT COUNT(*), COUNT(DISTINCT ProteinName),
               COUNT(DISTINCT PeptideSequence), COUNT(DISTINCT SampleID),
               COUNT(*) FILTER (
                   WHERE NOT isfinite(NormIntensity) OR NormIntensity <= 0
               )
        FROM read_parquet(?)
        """,
        [str(combined)],
    ).fetchone()
    return tuple(map(int, stats))


def prepare_peptides(
    paths: RefreshPaths,
    evidence: dict[str, Path],
    metadata: pd.DataFrame,
    threads: int,
    force: bool,
) -> Path:
    """Prepare all six batches and combine them into one peptide table."""
    prepared_dir = paths.work_dir / "prepared-batches"
    prepared_dir.mkdir(parents=True, exist_ok=True)
    combined = paths.work_dir / "quartet-peptides.parquet"
    connection = duckdb.connect()
    try:
        connection.execute("SET threads=?", [threads])
        connection.execute("SET memory_limit='80GB'")
        connection.execute("SET temp_directory=?", [str(paths.work_dir / "duckdb-tmp")])
        batch_paths = prepare_batch_files(
            connection, prepared_dir, evidence, metadata, force
        )
        combine_batches(connection, batch_paths, combined, force)
        stats = peptide_table_stats(connection, combined)
    finally:
        connection.close()
    rows, proteins, peptides, samples, invalid = stats
    if rows == 0 or proteins == 0 or peptides == 0 or samples != 72 or invalid:
        raise ValueError(
            "Invalid Quartet peptide table: "
            f"rows={rows}, proteins={proteins}, peptides={peptides}, "
            f"samples={samples}, invalid={invalid}"
        )
    print(
        f"prepared: {rows:,} rows, {proteins:,} proteins, "
        f"{peptides:,} peptides, {samples} samples"
    )
    return combined


def quantify_method(
    peptides: Path,
    output: Path,
    method: str,
    fasta: Path,
    threads: int,
) -> None:
    """Run one method through the current wheel's Rust kernel."""
    options = {
        "peptides": str(peptides),
        "method": method,
        "output": str(output),
        "threads": threads,
    }
    if method == "pibaq":
        options.update(fasta=str(fasta), enzyme="Trypsin", min_aa=7, max_aa=50)
    elif method == "directlfq":
        options["min_nonan"] = 1
    mokume.peptides2protein(**options)


def load_matrix(path: Path, method: str, samples: list[str]) -> pd.DataFrame:
    """Load a positive linear protein matrix and verify its sample contract."""
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
    if result.duplicated(["ProteinName", "SampleID"]).any():
        raise ValueError(f"{path} contains duplicate protein/sample cells")
    matrix = result.pivot(index="ProteinName", columns="SampleID", values=value_column)
    unknown = set(matrix.columns).difference(samples)
    if unknown:
        raise ValueError(f"{path} contains unknown samples: {sorted(unknown)}")
    if set(samples).difference(matrix.columns):
        raise ValueError(f"{path} does not contain all 72 samples")
    return matrix.reindex(columns=samples).sort_index()


def quantify_all(
    paths: RefreshPaths,
    peptides: Path,
    metadata: pd.DataFrame,
    threads: int,
    force: bool,
) -> tuple[dict[str, pd.DataFrame], dict[str, int]]:
    """Quantify all four methods and retain their pre-matching coverage."""
    output_dir = paths.work_dir / "protein-quant"
    output_dir.mkdir(parents=True, exist_ok=True)
    samples = metadata["run_id"].tolist()
    matrices = {}
    coverage = {}
    for method in METHODS:
        output = output_dir / f"{method}.tsv"
        if force or not output.exists():
            quantify_method(peptides, output, method, paths.fasta, threads)
        matrix = load_matrix(output, method, samples)
        matrices[method] = matrix
        coverage[method] = len(matrix)
        print(f"quantified {method}: {len(matrix):,} proteins x 72 samples")
    return matrices, coverage


def complete_and_matched_matrices(
    matrices: dict[str, pd.DataFrame],
) -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame]]:
    """Retain each complete universe and one fair universe across methods."""
    complete = {method: matrix.dropna() for method, matrix in matrices.items()}
    common = sorted(
        set.intersection(*(set(matrix.index) for matrix in complete.values()))
    )
    if len(common) < 10:
        raise ValueError(f"Only {len(common)} matched complete Quartet proteins")
    for method, matrix in complete.items():
        print(f"complete {method}: {len(matrix):,} proteins")
    print(f"matched complete universe: {len(common):,} proteins")
    matched = {method: matrix.loc[common].copy() for method, matrix in complete.items()}
    return complete, matched


def correct_one_method(
    paths: RefreshPaths,
    method: str,
    raw_log2: pd.DataFrame,
) -> pd.DataFrame:
    """Apply native Rust ComBat to a complete log2 protein matrix."""
    method_dir = paths.work_dir / "combat" / method
    method_dir.mkdir(parents=True, exist_ok=True)
    input_path = method_dir / "input.tsv"
    output_path = method_dir / "corrected.tsv"
    long_frame = (
        raw_log2.rename_axis("ProteinName")
        .reset_index()
        .melt(id_vars="ProteinName", var_name="SampleID", value_name="Log2Intensity")
    )
    long_frame.to_csv(input_path, sep="\t", index=False)
    if output_path.exists():
        output_path.unlink()
    mokume.correct_batches(
        folder=str(method_dir),
        pattern="input.tsv",
        output=str(output_path),
        sample_id_column="SampleID",
        protein_id_column="ProteinName",
        pibaq_raw_column="Log2Intensity",
        pibaq_corrected_column="CorrectedIntensity",
        sep="\t",
        comment="",
    )
    corrected = pd.read_csv(output_path, sep="\t")
    corrected["CorrectedIntensity"] = pd.to_numeric(
        corrected["CorrectedIntensity"], errors="coerce"
    )
    matrix = corrected.pivot(
        index="ProteinName", columns="SampleID", values="CorrectedIntensity"
    ).reindex(index=raw_log2.index, columns=raw_log2.columns)
    if matrix.shape != raw_log2.shape or not np.isfinite(matrix.to_numpy()).all():
        raise ValueError(f"Invalid ComBat output for {method}: {matrix.shape}")
    return matrix


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
    if (
        not np.isfinite(
            metrics[["Mean_CV", "SNR", "Inter_Batch_Correlation", "Batch_RMSE"]]
        )
        .all()
        .all()
    ):
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


def analyze_matrices(
    paths: RefreshPaths,
    matrices: dict[str, pd.DataFrame],
    coverage: dict[str, int],
) -> QuartetAnalysis:
    """Build complete, matched, and corrected matrix views."""
    complete, matched = complete_and_matched_matrices(matrices)
    raw_full = {method: np.log2(matrix) for method, matrix in matrices.items()}
    raw_complete = {method: np.log2(matrix) for method, matrix in complete.items()}
    corrected_complete = {
        method: correct_one_method(paths, method, raw_complete[method])
        for method in METHODS
    }
    raw_matched = {method: np.log2(matrix) for method, matrix in matched.items()}
    corrected_matched = {
        method: corrected_complete[method].loc[raw_matched[method].index]
        for method in METHODS
    }
    complete_coverage = {method: len(matrix) for method, matrix in complete.items()}
    return QuartetAnalysis(
        raw_full=raw_full,
        corrected_complete=corrected_complete,
        raw_matched=raw_matched,
        corrected_matched=corrected_matched,
        coverage=coverage,
        complete_coverage=complete_coverage,
    )


def main() -> None:
    """Run source preparation, Rust kernels, metrics, and plots."""
    args = parse_args()
    paths, evidence = resolve_paths(args)
    connection = duckdb.connect()
    try:
        connection.execute("SET threads=?", [args.threads])
        metadata = build_metadata(connection, evidence)
    finally:
        connection.close()
    peptides = prepare_peptides(paths, evidence, metadata, args.threads, args.force)
    matrices, coverage = quantify_all(
        paths, peptides, metadata, args.threads, args.force
    )
    analysis = analyze_matrices(paths, matrices, coverage)
    metrics = calculate_metrics(
        analysis.raw_matched,
        analysis.corrected_matched,
        analysis.coverage,
        analysis.complete_coverage,
        metadata,
    )
    save_results(
        paths.results_dir,
        metadata,
        analysis.raw_full,
        analysis.corrected_complete,
        metrics,
    )
    plot_batch_effect_diagnosis(
        analysis.raw_matched,
        analysis.corrected_matched,
        metadata,
        paths.figures_dir / "batch_effect_diagnosis.png",
    )
    plot_pca(
        analysis.raw_matched,
        analysis.corrected_matched,
        metadata,
        paths.figures_dir / "pca_comparison.png",
    )
    print("\nRefreshed Quartet metrics")
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
