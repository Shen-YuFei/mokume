"""Lightweight post-run matrix summaries for Mokume Studio."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from mokume.studio.catalog import command_schema
from mokume.studio.models import ProjectRecord, RunRecord, RunStatus
from mokume.studio.paths import ProjectPaths


def qc_summary(
    record: RunRecord,
    artifacts: list[dict[str, Any]],
    project: ProjectRecord,
) -> dict[str, Any]:
    """Build a bounded summary from registered tabular artifacts."""
    if record.status is not RunStatus.SUCCEEDED:
        return {"available": False, "reason": "QC is available after a successful run"}
    candidates = _qc_tables(artifacts, project)
    if not candidates:
        return {
            "available": False,
            "reason": "No tabular quantification artifact was produced",
        }
    summaries = []
    errors = []
    for path in candidates:
        try:
            summaries.append((path, _summarize_matrix(path)))
        except (OSError, RuntimeError, TypeError, ValueError) as exc:
            errors.append(str(exc))
    if not summaries:
        reason = errors[0] if errors else "no readable table"
        return {"available": False, "reason": f"QC summary unavailable: {reason}"}
    table, summary = max(summaries, key=lambda item: item[0].stat().st_size)
    entity_counts: dict[str, int] = {}
    for _path, item in summaries:
        label = item["entity_label"]
        entity_counts[label] = max(entity_counts.get(label, 0), item["entity_count"])
    summary.update(
        available=True,
        artifact_path=str(table),
        entity_counts=entity_counts,
        normalization=_normalization_summary(record.argv),
        reports=_report_artifacts(artifacts),
    )
    return summary


def _qc_tables(artifacts: list[dict[str, Any]], project: ProjectRecord) -> list[Path]:
    guard = ProjectPaths(project.root)
    candidates = []
    for artifact in artifacts:
        path = Path(artifact["path"])
        if path.suffix.casefold() not in {".csv", ".tsv", ".txt"}:
            continue
        try:
            candidates.append(guard.resolve_existing(path))
        except (OSError, RuntimeError, ValueError):
            continue
    return candidates


def _summarize_matrix(path: Path) -> dict[str, Any]:
    frame, numeric, positive, log_values = _matrix_data(path)
    missing, quartiles, median_cv = _matrix_statistics(positive, log_values)
    correlations, pca = _sample_diagnostics(log_values, list(numeric.columns))
    return {
        "entity_count": int(len(frame)),
        "entity_label": _entity_label(frame.columns[0]),
        "sample_count": int(numeric.shape[1]),
        "missing_percent": missing,
        "median_cv_percent": median_cv,
        "log2_intensity_quartiles": quartiles,
        "correlation": correlations,
        "pca": pca,
        "truncated": len(frame) == 50000,
    }


def _matrix_data(path: Path) -> tuple[Any, Any, Any, Any]:
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        first_line = stream.readline()
    delimiter = "\t" if first_line.count("\t") > first_line.count(",") else ","
    frame = pd.read_csv(path, sep=delimiter, nrows=50000)
    numeric = frame.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    numeric = numeric.loc[:, numeric.notna().any(axis=0)]
    if frame.empty or numeric.empty:
        raise ValueError("table has no numeric sample matrix")
    values = numeric.to_numpy(dtype=float)
    observed = np.isfinite(values) & (values > 0)
    positive = np.where(observed, values, np.nan)
    return frame, numeric, positive, np.log2(positive)


def _matrix_statistics(
    positive: Any, log_values: Any
) -> tuple[float, list[float], float | None]:
    observed = np.isfinite(positive)
    missing = 100.0 * (1.0 - float(observed.sum()) / observed.size)
    means = np.nanmean(positive, axis=1)
    counts = np.sum(observed, axis=1)
    deviations = np.nanstd(positive, axis=1, ddof=1)
    cvs = np.divide(
        deviations,
        means,
        out=np.full_like(means, np.nan),
        where=(means > 0) & (counts > 1),
    )
    quartiles = [
        round(float(value), 3) for value in np.nanpercentile(log_values, [25, 50, 75])
    ]
    median_cv = (
        round(float(np.nanmedian(cvs) * 100), 2) if np.isfinite(cvs).any() else None
    )
    return round(missing, 2), quartiles, median_cv


def _sample_diagnostics(
    values: Any, samples: list[str]
) -> tuple[dict[str, Any], dict[str, Any]]:
    empty_correlation = {"median": None, "lowest_sample": None, "outliers": []}
    empty_pca = {"variance_percent": [], "outliers": []}
    if len(samples) < 2 or not np.isfinite(values).any():
        return empty_correlation, empty_pca
    with np.errstate(all="ignore"):
        row_medians = np.nanmedian(values, axis=1)
    global_median = float(np.nanmedian(values))
    filled = np.where(np.isfinite(values), values, row_medians[:, None])
    filled = np.nan_to_num(filled, nan=global_median)
    return _correlation_diagnostics(filled, samples), _pca_diagnostics(filled, samples)


def _correlation_diagnostics(filled: Any, samples: list[str]) -> dict[str, Any]:
    correlation = np.corrcoef(filled.T)
    means = (correlation.sum(axis=1) - 1.0) / max(1, len(samples) - 1)
    if not np.isfinite(means).any():
        return {"median": None, "lowest_sample": None, "outliers": []}
    median_correlation = float(np.nanmedian(means))
    spread = float(np.nanmedian(np.abs(means - median_correlation)))
    cutoff = min(0.7, median_correlation - 3 * spread)
    outliers = [sample for sample, value in zip(samples, means) if value < cutoff]
    return {
        "median": round(median_correlation, 4),
        "lowest_sample": samples[int(np.nanargmin(means))],
        "outliers": outliers,
    }


def _pca_diagnostics(filled: Any, samples: list[str]) -> dict[str, Any]:
    centered = filled.T - filled.T.mean(axis=0, keepdims=True)
    u_matrix, singular, _vectors = np.linalg.svd(centered, full_matrices=False)
    variance = singular**2
    ratio = variance[:2] / variance.sum() * 100 if variance.sum() else np.array([])
    scores = u_matrix[:, :2] * singular[:2]
    distances = np.linalg.norm(scores - np.median(scores, axis=0), axis=1)
    median_distance = float(np.median(distances))
    mad = float(np.median(np.abs(distances - median_distance)))
    pca_outliers = [
        sample
        for sample, value in zip(samples, distances)
        if mad and value > median_distance + 3 * mad
    ]
    return {
        "variance_percent": [round(float(value), 2) for value in ratio],
        "outliers": pca_outliers,
    }


def _entity_label(identifier: Any) -> str:
    name = str(identifier).casefold()
    if "protein" in name:
        return "proteins"
    if "peptide" in name or "sequence" in name:
        return "peptides"
    return "features"


def _normalization_summary(argv: list[str]) -> dict[str, Any]:
    options = _command_options(argv)
    names = ("run-normalization", "sample-normalization", "impute-method")
    methods = {name: options[name] for name in names if name in options}
    if "normalize" in options:
        methods["normalize"] = True
    return {
        "methods": methods,
        "comparison_available": False,
        "message": "This workflow did not emit a matched pre-normalization matrix",
    }


def _command_options(argv: list[str]) -> dict[str, str | bool]:
    path_length = 1
    for command in command_schema():
        if argv[: len(command["path"])] == command["path"]:
            path_length = len(command["path"])
            break
    options: dict[str, str | bool] = {}
    index = path_length
    while index < len(argv):
        token = argv[index]
        if not token.startswith("--"):
            index += 1
            continue
        name = token[2:]
        if "=" in name:
            name, value = name.split("=", 1)
            options[name] = value
            index += 1
        elif index + 1 < len(argv) and not argv[index + 1].startswith("-"):
            options[name] = argv[index + 1]
            index += 2
        else:
            options[name] = True
            index += 1
    return options


def _report_artifacts(artifacts: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {"id": item["id"], "path": item["path"], "media_type": item["media_type"]}
        for item in artifacts
        if item["media_type"] in {"text/html", "application/pdf"}
        or Path(item["path"]).suffix.casefold() in {".html", ".pdf"}
    ]
