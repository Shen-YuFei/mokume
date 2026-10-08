"""Structured preflight review and persisted run insights for Studio."""

from __future__ import annotations

import json
import re
import shlex
import shutil
from collections.abc import Callable, Iterable
from datetime import datetime
from pathlib import Path
from typing import Any

import psutil
from fastapi import Depends, FastAPI, HTTPException, status

from mokume.studio.catalog import (
    OUTPUT_ARGUMENTS,
    command_schema,
    validate_and_canonicalize,
)
from mokume.studio.command_insights import (
    parse_occurrences,
    planned_workflow_steps,
    template_from_argv,
)
from mokume.studio.models import ProjectRecord, RunRecord, RunStatus, ValidationRequest
from mokume.studio.paths import ProjectPaths
from mokume.studio.qc import qc_summary
from mokume.studio.sdrf import read_sdrf


Dependency = Callable[..., Any]
_OPTION_PATTERN = re.compile(r"--([a-z0-9-]+)")
_FILE_SUFFIXES = {
    "parquet": (".parquet",),
    "msstats": (".csv", ".tsv"),
    "sdrf": (".tsv", ".txt"),
    "fasta": (".fa", ".faa", ".fasta"),
    "filter-config": (".json", ".yaml", ".yml"),
    "config": (".json", ".yaml", ".yml"),
    "families": (".json", ".yaml", ".yml"),
    "pibaq-families": (".json", ".yaml", ".yml"),
}


def install_insight_routes(
    app: FastAPI,
    runtime: Any,
    require_session: Dependency,
    require_mutation: Dependency,
) -> None:
    """Install review and run-detail endpoints on the authenticated app."""

    @app.post("/api/commands/review")
    async def command_review(
        payload: ValidationRequest,
        _session=Depends(require_mutation),
    ) -> dict:
        project = runtime.store.active_project()
        if project is None:
            raise HTTPException(status.HTTP_409_CONFLICT, "Open a project folder first")
        return review_command(payload.argv, project.root)

    @app.get("/api/runs/{run_id}/details")
    async def run_details(run_id: str, _session=Depends(require_session)) -> dict:
        project = runtime.store.active_project()
        record = runtime.store.get_run(run_id)
        if project is None or record is None or record.project_id != project.id:
            raise HTTPException(status.HTTP_404_NOT_FOUND, "Run not found")
        return build_run_details(runtime.store, record, project)


def review_command(argv: list[str], project_root: str | Path) -> dict:
    """Return a non-mutating, structured review even when validation fails."""
    spec = _find_command(argv)
    guard = ProjectPaths(project_root)
    canonical, validation_error = _try_canonicalize(argv, guard.root)
    inspected_argv = canonical or argv
    occurrences = parse_occurrences(inspected_argv, spec) if spec else []
    inspection = _review_inspection(guard, occurrences)
    checks = _review_checks(validation_error, occurrences, inspection)
    valid = canonical is not None and not any(
        item["status"] == "error" for item in checks
    )
    return {
        "valid": valid,
        "error": validation_error,
        "error_parameter": _error_parameter(validation_error),
        "argv": canonical or argv,
        "command": shlex.join(canonical or argv),
        "workflow": list(spec["path"]) if spec else argv[:1],
        "display_name": _workflow_display_name(spec, argv),
        "template": template_from_argv(canonical or argv, spec),
        "checks": checks,
        **inspection,
        "planned_steps": planned_workflow_steps(canonical or argv),
    }


def _workflow_display_name(spec: dict[str, Any] | None, argv: list[str]) -> str:
    if spec is None:
        return "Unknown"
    return spec.get("display_name", " ".join(spec["path"] or argv[:1]))


def _review_inspection(
    guard: ProjectPaths, occurrences: list[dict[str, Any]]
) -> dict[str, Any]:
    inputs = _inspect_paths(guard, occurrences, output=False)
    return {
        "inputs": inputs,
        "outputs": _inspect_paths(guard, occurrences, output=True),
        "resources": _resource_estimate(occurrences, inputs, guard.root),
        "sdrf": _sdrf_preview(guard, occurrences),
    }


def build_run_details(store: Any, record: RunRecord, project: ProjectRecord) -> dict:
    """Assemble one workspace-scoped run view from immutable local records."""
    directory = _run_directory(record, project)
    parameters = _read_json(directory / "parameters.json") if directory else {}
    provenance = _read_json(directory / "provenance.json") if directory else {}
    artifacts = [
        item.model_dump(mode="json") for item in store.list_artifacts(record.id)
    ]
    events = store.events_after(record.id)
    return {
        "run": record.model_dump(mode="json"),
        "command": shlex.join(record.argv),
        "duration_seconds": _duration_seconds(record),
        "template": template_from_argv(record.argv, _find_command(record.argv)),
        "parameters": parameters,
        "provenance": provenance,
        "artifacts": artifacts,
        "logs": _run_logs(directory),
        "stages": _stage_summary(record, events),
        "planned_steps": planned_workflow_steps(record.argv),
        "qc": qc_summary(record, artifacts, project),
    }


def _find_command(argv: list[str]) -> dict[str, Any] | None:
    for spec in command_schema():
        path = list(spec["path"])
        if argv[: len(path)] == path:
            return spec
    return None


def _try_canonicalize(
    argv: list[str], root: Path
) -> tuple[list[str] | None, str | None]:
    try:
        return validate_and_canonicalize(argv, root), None
    except (OSError, RuntimeError, ValueError) as exc:
        return None, str(exc)


def _path_values(occurrence: dict[str, Any]) -> Iterable[tuple[str, str, str]]:
    flag = occurrence["flag"]
    names = flag.get("value_names") or []
    parameter = str(flag.get("long") or flag.get("id"))
    for position, value in enumerate(occurrence["values"]):
        name = names[min(position, len(names) - 1)] if names else ""
        if name in {"FILE", "DIR"}:
            yield parameter, name.lower(), value


def _inspect_paths(
    guard: ProjectPaths,
    occurrences: list[dict[str, Any]],
    *,
    output: bool,
) -> list[dict[str, Any]]:
    records = []
    for occurrence in occurrences:
        parameter = str(occurrence["flag"].get("long") or "")
        if (parameter in OUTPUT_ARGUMENTS) is not output:
            continue
        for name, expected_kind, requested in _path_values(occurrence):
            records.append(_inspect_path(guard, name, expected_kind, requested, output))
    return records


def _inspect_path(
    guard: ProjectPaths,
    parameter: str,
    expected_kind: str,
    requested: str,
    output: bool,
) -> dict[str, Any]:
    record = {
        "parameter": parameter,
        "path": requested,
        "kind": expected_kind,
        "status": "ok",
        "message": "Output path is available" if output else "Input is readable",
        "size_bytes": None,
    }
    try:
        path = (
            guard.resolve_output(requested, allow_existing=True)
            if output
            else guard.resolve_existing(requested)
        )
        record["path"] = str(path)
        if output:
            record["exists"] = path.exists()
            if path.exists():
                record.update(status="error", message="Output already exists")
            return record
        actual_kind = "dir" if path.is_dir() else "file"
        if actual_kind != expected_kind:
            record.update(
                status="error", message=f"Expected {expected_kind}, found {actual_kind}"
            )
        elif expected_kind == "file" and not _suffix_matches(parameter, path):
            record.update(
                status="warning", message="File type is unusual for this parameter"
            )
        record["size_bytes"] = _path_size(path)
        return record
    except (OSError, RuntimeError, ValueError) as exc:
        record.update(
            status="error", message=str(exc), exists=False if output else None
        )
        return record


def _suffix_matches(parameter: str, path: Path) -> bool:
    suffixes = _FILE_SUFFIXES.get(parameter)
    return suffixes is None or path.name.casefold().endswith(suffixes)


def _path_size(path: Path) -> int:
    if path.is_file():
        return path.stat().st_size
    total = 0
    for position, candidate in enumerate(path.rglob("*")):
        if position >= 5000:
            break
        if candidate.is_file():
            total += candidate.stat().st_size
    return total


def _review_checks(
    error: str | None,
    occurrences: list[dict[str, Any]],
    inspection: dict[str, Any],
) -> list[dict[str, Any]]:
    parameter = _error_parameter(error)
    inputs = inspection["inputs"]
    outputs = inspection["outputs"]
    metadata = inspection["sdrf"]
    resources = inspection["resources"]
    checks = [
        {
            "id": "contract",
            "label": "Workflow contract",
            "status": "error" if error else "ok",
            "message": error or "Parameters and dependencies are valid",
            "parameter": parameter,
        }
    ]
    checks.extend(_path_check("input", item) for item in inputs)
    checks.extend(_path_check("output", item) for item in outputs)
    if metadata:
        issues = metadata.get("issues", [])
        checks.append(
            {
                "id": "sdrf",
                "label": "SDRF sample mapping",
                "status": "warning" if issues else "ok",
                "message": f"{metadata['row_count']} rows; {len(issues)} mapping warnings",
                "parameter": "sdrf",
            }
        )
        checks.extend(_metadata_dependency_checks(occurrences, metadata))
    constrained = []
    if resources["suggested_peak_memory_bytes"] > resources["available_memory_bytes"]:
        constrained.append("estimated peak memory exceeds currently available memory")
    if resources["suggested_free_disk_bytes"] > resources["free_disk_bytes"]:
        constrained.append("estimated output space exceeds currently free disk")
    checks.append(
        {
            "id": "resources",
            "label": "Resource availability",
            "status": "warning" if constrained else "ok",
            "message": "; ".join(constrained)
            or "Estimated resources are currently available",
            "parameter": "memory" if constrained else None,
        }
    )
    return checks


def _metadata_dependency_checks(
    occurrences: list[dict[str, Any]], metadata: dict[str, Any]
) -> list[dict[str, Any]]:
    checks = []
    options = _options_by_name(occurrences)
    condition_required = options.get("sample-normalization") == ["condition-median"]
    condition_required |= any(
        name in options
        for name in (
            "de-contrast",
            "de-contrast-file",
            "coverage-threshold",
            "min-sample-correlation",
        )
    )
    if condition_required:
        checks.append(
            _metadata_column_check(
                "condition", metadata["columns"].get("condition"), "sdrf"
            )
        )
    requested_columns = [
        ("batch-column", options.get("batch-column", [])),
        ("batch-covariate", options.get("batch-covariate", [])),
        ("irs-sdrf-column", options.get("irs-sdrf-column", [])),
    ]
    headers = set(metadata.get("headers", []))
    for parameter, columns in requested_columns:
        for column in columns:
            checks.append(_metadata_column_check(column, column in headers, parameter))
    reference_check = _default_reference_check(options, metadata)
    if reference_check:
        checks.append(reference_check)
    return checks


def _options_by_name(
    occurrences: list[dict[str, Any]],
) -> dict[str, list[str]]:
    options: dict[str, list[str]] = {}
    for item in occurrences:
        name = str(item["flag"].get("long"))
        options.setdefault(name, []).extend(item["values"])
    return options


def _default_reference_check(
    options: dict[str, list[str]], metadata: dict[str, Any]
) -> dict[str, Any] | None:
    default_reference = (
        options.get("quant-method") == ["ratio"] or "irs" in options
    ) and not any(
        selector in options
        for selector in (
            "irs-reference-sample",
            "irs-reference-regex",
            "irs-sdrf-column",
        )
    )
    if not default_reference:
        return None
    found = metadata.get("reference_count", 0) > 0
    return {
        "id": "sdrf:reference",
        "label": "SDRF reference samples",
        "status": "ok" if found else "error",
        "message": (
            f"{metadata['reference_count']} reference rows detected"
            if found
            else "No default IRS reference samples were detected"
        ),
        "parameter": "sdrf",
    }


def _metadata_column_check(
    label: str, available: str | bool | None, parameter: str
) -> dict[str, Any]:
    return {
        "id": f"sdrf:column:{label}",
        "label": f"SDRF column: {label}",
        "status": "ok" if available else "error",
        "message": "Column is available"
        if available
        else "Required column was not found",
        "parameter": parameter,
    }


def _path_check(kind: str, item: dict[str, Any]) -> dict[str, Any]:
    return {
        "id": f"{kind}:{item['parameter']}:{item['path']}",
        "label": f"--{item['parameter']}",
        "status": item["status"],
        "message": item["message"],
        "parameter": item["parameter"],
    }


def _error_parameter(error: str | None) -> str | None:
    if not error:
        return None
    match = _OPTION_PATTERN.search(error)
    return match.group(1) if match else None


def _occurrence_value(occurrences: list[dict[str, Any]], name: str) -> str | None:
    for occurrence in reversed(occurrences):
        flag_name = occurrence["flag"].get("long")
        if flag_name == name and occurrence["values"]:
            return occurrence["values"][0]
    return None


def _resource_estimate(
    occurrences: list[dict[str, Any]],
    inputs: list[dict[str, Any]],
    project_root: Path,
) -> dict[str, Any]:
    total = sum(item.get("size_bytes") or 0 for item in inputs)
    threads = _occurrence_value(occurrences, "threads") or "24"
    try:
        thread_count = int(threads)
    except ValueError:
        thread_count = 24
    return {
        "threads": thread_count,
        "configured_memory": _occurrence_value(occurrences, "memory"),
        "input_bytes": total,
        "suggested_peak_memory_bytes": max(512 * 1024**2, total * 3),
        "suggested_free_disk_bytes": max(256 * 1024**2, total * 2),
        "available_memory_bytes": psutil.virtual_memory().available,
        "free_disk_bytes": shutil.disk_usage(project_root).free,
        "estimate_only": True,
    }


def _sdrf_preview(
    guard: ProjectPaths, occurrences: list[dict[str, Any]]
) -> dict[str, Any] | None:
    requested = _occurrence_value(occurrences, "sdrf")
    if requested is None:
        return None
    try:
        path = guard.resolve_existing(requested)
        return read_sdrf(path)
    except (OSError, RuntimeError, ValueError) as exc:
        return {
            "path": requested,
            "row_count": 0,
            "columns": {},
            "rows": [],
            "issues": [str(exc)],
        }


def _run_directory(record: RunRecord, project: ProjectRecord) -> Path | None:
    try:
        path = ProjectPaths(project.root).resolve_existing(record.run_directory)
    except (OSError, RuntimeError, ValueError):
        return None
    return path if path.is_dir() else None


def _read_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return payload if isinstance(payload, dict) else {}


def _run_logs(directory: Path | None) -> dict[str, list[str]]:
    if directory is None:
        return {"stdout": [], "stderr": []}
    return {
        name: _tail_lines(directory / f"{name}.log") for name in ("stdout", "stderr")
    }


def _tail_lines(path: Path, limit: int = 300) -> list[str]:
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return []
    return lines[-limit:]


def _duration_seconds(record: RunRecord) -> float | None:
    if not record.started_at or not record.finished_at:
        return None
    started = datetime.fromisoformat(record.started_at)
    finished = datetime.fromisoformat(record.finished_at)
    return max(0.0, (finished - started).total_seconds())


def _stage_summary(
    record: RunRecord, events: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    definitions = (
        ("inputs", "Read & verify inputs"),
        ("workflow", "Run configured workflow"),
        ("artifacts", "Collect outputs"),
        ("provenance", "Publish provenance"),
    )
    stage_events = [event for event in events if event["type"] == "stage"]
    stages = []
    for stage_id, label in definitions:
        matching = [
            event["payload"]
            for event in stage_events
            if event["payload"].get("stage") == stage_id
        ]
        latest = matching[-1] if matching else {}
        stages.append(
            {
                "id": stage_id,
                "label": label,
                "status": latest.get(
                    "status", _derived_stage_status(record.status, stage_id)
                ),
                "elapsed_seconds": latest.get("elapsed_seconds"),
                "error": latest.get("error"),
            }
        )
    return stages


def _derived_stage_status(status_value: RunStatus, stage_id: str) -> str:
    if status_value is RunStatus.SUCCEEDED:
        return "succeeded"
    if status_value in {RunStatus.FAILED, RunStatus.CANCELLED, RunStatus.INTERRUPTED}:
        return "failed" if stage_id == "workflow" else "pending"
    if status_value is RunStatus.RUNNING:
        return "running" if stage_id == "workflow" else "pending"
    return "pending"
