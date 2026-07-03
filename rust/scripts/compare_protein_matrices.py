#!/usr/bin/env python3
"""Compare Python and Rust mokume protein matrix outputs.

Example real-data parity command:

  python scripts/compare_protein_matrices.py \
    --python /tmp/mokume-parity/PXD003539/python/sum_none.csv \
    --rust /tmp/mokume-parity/PXD003539/rust/sum_none.csv \
    --output /tmp/mokume-parity/PXD003539/sum_none_report.json

The JSON report includes protein/sample counts, set differences, shared cell
counts, absolute and relative error summaries, per-sample Spearman correlation,
and per-sample total intensity differences.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Iterable


@dataclass(frozen=True)
class Matrix:
    protein_column: str
    samples: list[str]
    values: dict[str, dict[str, float | None]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare Python and Rust mokume protein matrices.",
        epilog=(
            "Reports contain protein_counts, sample_counts, protein_set_differences, "
            "sample_set_differences, shared_matrix_cells, error_metrics, "
            "sample_metrics with Spearman and total intensity deltas, and "
            "largest_cell_differences. Write reports under /tmp/mokume-parity/<PXD>/ "
            "for real-data parity runs."
        ),
    )
    parser.add_argument("--python", required=True, type=Path, help="Python output CSV")
    parser.add_argument("--rust", required=True, type=Path, help="Rust output CSV")
    parser.add_argument("--output", type=Path, help="Optional JSON report path")
    parser.add_argument(
        "--top-n",
        type=int,
        default=25,
        help="Number of largest cell differences to include",
    )
    return parser.parse_args()


def load_matrix(path: Path) -> Matrix:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"{path} has no header")
        protein_column = reader.fieldnames[0]
        samples = reader.fieldnames[1:]
        values: dict[str, dict[str, float | None]] = {}
        for row in reader:
            protein = (row.get(protein_column) or "").strip()
            if not protein:
                continue
            values[protein] = {
                sample: parse_float(row.get(sample, "")) for sample in samples
            }
    return Matrix(protein_column=protein_column, samples=samples, values=values)


def parse_float(value: str | None) -> float | None:
    if value is None:
        return None
    value = value.strip()
    if not value:
        return None
    parsed = float(value)
    return parsed if math.isfinite(parsed) else None


def percentile(values: list[float], quantile: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    position = (len(ordered) - 1) * quantile
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def ranks(values: list[float]) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda item: item[1])
    output = [0.0] * len(values)
    index = 0
    while index < len(indexed):
        end = index + 1
        while end < len(indexed) and indexed[end][1] == indexed[index][1]:
            end += 1
        rank = (index + end - 1) / 2.0 + 1.0
        for original, _ in indexed[index:end]:
            output[original] = rank
        index = end
    return output


def pearson(left: list[float], right: list[float]) -> float | None:
    if len(left) < 2 or len(right) < 2:
        return None
    left_mean = sum(left) / len(left)
    right_mean = sum(right) / len(right)
    numerator = sum((a - left_mean) * (b - right_mean) for a, b in zip(left, right))
    left_den = math.sqrt(sum((a - left_mean) ** 2 for a in left))
    right_den = math.sqrt(sum((b - right_mean) ** 2 for b in right))
    if left_den == 0.0 or right_den == 0.0:
        return None
    return numerator / (left_den * right_den)


def spearman(left: list[float], right: list[float]) -> float | None:
    if len(left) < 2:
        return None
    return pearson(ranks(left), ranks(right))


def finite_pairs(
    python_matrix: Matrix,
    rust_matrix: Matrix,
    proteins: Iterable[str],
    samples: Iterable[str],
) -> Iterable[tuple[str, str, float, float]]:
    for protein in proteins:
        py_row = python_matrix.values.get(protein, {})
        rs_row = rust_matrix.values.get(protein, {})
        for sample in samples:
            py_value = py_row.get(sample)
            rs_value = rs_row.get(sample)
            if py_value is not None and rs_value is not None:
                yield protein, sample, py_value, rs_value


def compare(python_matrix: Matrix, rust_matrix: Matrix, top_n: int) -> dict:
    python_proteins = set(python_matrix.values)
    rust_proteins = set(rust_matrix.values)
    python_samples = set(python_matrix.samples)
    rust_samples = set(rust_matrix.samples)
    shared_proteins = sorted(python_proteins & rust_proteins)
    shared_samples = sorted(python_samples & rust_samples)

    abs_errors: list[float] = []
    rel_errors: list[float] = []
    largest: list[dict] = []
    for protein, sample, py_value, rs_value in finite_pairs(
        python_matrix, rust_matrix, shared_proteins, shared_samples
    ):
        abs_error = abs(py_value - rs_value)
        denominator = max(abs(py_value), 1e-12)
        rel_error = abs_error / denominator
        abs_errors.append(abs_error)
        rel_errors.append(rel_error)
        largest.append(
            {
                "protein": protein,
                "sample": sample,
                "python": py_value,
                "rust": rs_value,
                "absolute_error": abs_error,
                "relative_error": rel_error,
            }
        )

    largest.sort(key=lambda item: item["absolute_error"], reverse=True)

    sample_metrics = {}
    for sample in shared_samples:
        py_values: list[float] = []
        rs_values: list[float] = []
        for protein in shared_proteins:
            py_value = python_matrix.values.get(protein, {}).get(sample)
            rs_value = rust_matrix.values.get(protein, {}).get(sample)
            if py_value is not None and rs_value is not None:
                py_values.append(py_value)
                rs_values.append(rs_value)
        py_total = sum(py_values)
        rs_total = sum(rs_values)
        sample_metrics[sample] = {
            "shared_cells": len(py_values),
            "spearman": spearman(py_values, rs_values),
            "python_total_intensity": py_total,
            "rust_total_intensity": rs_total,
            "total_intensity_difference": rs_total - py_total,
            "total_intensity_relative_difference": relative_difference(
                rs_total, py_total
            ),
        }

    error_metrics = {
        "absolute": {
            "median": median(abs_errors) if abs_errors else None,
            "p95": percentile(abs_errors, 0.95),
            "max": max(abs_errors) if abs_errors else None,
        },
        "relative": {
            "median": median(rel_errors) if rel_errors else None,
            "p95": percentile(rel_errors, 0.95),
            "max": max(rel_errors) if rel_errors else None,
        },
    }

    return {
        "protein_counts": {
            "python": len(python_proteins),
            "rust": len(rust_proteins),
            "shared": len(shared_proteins),
            "python_only": len(python_proteins - rust_proteins),
            "rust_only": len(rust_proteins - python_proteins),
        },
        "sample_counts": {
            "python": len(python_samples),
            "rust": len(rust_samples),
            "shared": len(shared_samples),
            "python_only": len(python_samples - rust_samples),
            "rust_only": len(rust_samples - python_samples),
        },
        "protein_set_differences": {
            "python_only": sorted(python_proteins - rust_proteins)[:top_n],
            "rust_only": sorted(rust_proteins - python_proteins)[:top_n],
        },
        "sample_set_differences": {
            "python_only": sorted(python_samples - rust_samples),
            "rust_only": sorted(rust_samples - python_samples),
        },
        "shared_matrix_cells": len(abs_errors),
        "error_metrics": error_metrics,
        "median_absolute_error": error_metrics["absolute"]["median"],
        "p95_absolute_error": error_metrics["absolute"]["p95"],
        "median_relative_error": error_metrics["relative"]["median"],
        "sample_metrics": sample_metrics,
        "largest_cell_differences": largest[:top_n],
    }


def relative_difference(value: float, reference: float) -> float | None:
    denominator = max(abs(reference), 1e-12)
    result = (value - reference) / denominator
    return result if math.isfinite(result) else None


def main() -> int:
    args = parse_args()
    report = compare(load_matrix(args.python), load_matrix(args.rust), args.top_n)
    rendered = json.dumps(report, indent=2, sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
