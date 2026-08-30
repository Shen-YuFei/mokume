"""Isolated worker entry point for one immutable Mokume Studio job."""

from __future__ import annotations

import argparse
import hashlib
import json
import mimetypes
import os
import platform
import sys
import uuid
from pathlib import Path

import mokume

from mokume.agentic.knowledge import load_knowledge_graph
from mokume.agentic.service import (
    EvaluationRequest,
    InspectionRequest,
    RecommendationService,
)
from mokume.studio.catalog import command_paths, validate_and_canonicalize
from mokume.studio.jobs import (
    path_snapshot,
    validate_spec_integrity,
    write_terminal_files,
)
from mokume.studio.models import (
    ArtifactRecord,
    JobOperation,
    JobSpec,
    RunStatus,
    utc_now,
)
from mokume.studio.paths import ProjectPaths
from mokume.studio.science import DatasetStatus, ScienceStore
from mokume.studio.state import StateStore


def build_parser() -> argparse.ArgumentParser:
    """Build the private worker parser."""
    parser = argparse.ArgumentParser(prog="python -m mokume.studio.worker")
    parser.add_argument("--state-directory", required=True)
    parser.add_argument("--run-id", required=True)
    return parser


def file_identity(path: Path, guard: ProjectPaths) -> dict:
    """Record a streaming content identity without loading files into memory."""
    identity = path_snapshot(path, guard)
    if path.is_file():
        identity["sha256"] = _sha256(path)
    else:
        for entry in identity.get("entries", []):
            entry["sha256"] = _sha256(Path(entry["resolved_path"]))
    return identity


def execute(store: StateStore, spec: JobSpec) -> int:
    """Revalidate and execute one native or scientific operation."""
    _validate_spec(spec)
    guard, resolved_inputs = _validated_inputs(spec)
    input_identities = [file_identity(path, guard) for path in resolved_inputs]
    started_at = utc_now()
    store.update_run(spec.run_id, RunStatus.RUNNING, worker_pid=os.getpid())
    output_paths, knowledge_fingerprint = _execute_operation(store, spec, guard)
    artifacts = _register_outputs(store, spec.run_id, guard, output_paths)
    provenance = _provenance(
        spec,
        artifacts,
        input_identities,
        started_at,
        knowledge_fingerprint,
    )
    store.update_run(spec.run_id, RunStatus.SUCCEEDED)
    write_terminal_files(store, spec.run_id, spec.run_directory, provenance)
    return 0


def _validate_spec(spec: JobSpec) -> None:
    validate_spec_integrity(spec)
    if spec.operation is JobOperation.NATIVE:
        canonical = validate_and_canonicalize(spec.argv, spec.project_root)
        if canonical != spec.argv:
            raise RuntimeError("job argv changed during worker validation")


def _validated_inputs(spec: JobSpec) -> tuple[ProjectPaths, list[Path]]:
    guard = ProjectPaths(spec.project_root)
    guard.resolve_existing(spec.run_directory)
    snapshots = spec.parameters.get("input_snapshots", [])
    resolved_inputs = [guard.resolve_existing(item["path"]) for item in snapshots]
    input_snapshots = [path_snapshot(path, guard) for path in resolved_inputs]
    if input_snapshots != snapshots:
        raise RuntimeError("input files changed after the run was approved")
    return guard, resolved_inputs


def _execute_operation(
    store: StateStore,
    spec: JobSpec,
    guard: ProjectPaths,
) -> tuple[list[Path], str | None]:
    if spec.operation is JobOperation.NATIVE:
        _inputs, outputs = command_paths(spec.argv)
        mokume.run(spec.argv)
        return outputs, None
    if spec.operation is JobOperation.INSPECT_DATASET:
        return _inspect_dataset(store, spec)
    if spec.operation is JobOperation.EVALUATE_RECOMMENDATION:
        return _evaluate_recommendation(spec, guard)
    raise RuntimeError(f"unsupported worker operation: {spec.operation.value}")


def _agentic_service(spec: JobSpec):
    fingerprint = load_knowledge_graph().fingerprint
    if spec.payload.get("knowledge_fingerprint") != fingerprint:
        raise RuntimeError("knowledge snapshot changed after approval")
    return RecommendationService(), fingerprint


def _inspect_dataset(
    store: StateStore,
    spec: JobSpec,
) -> tuple[list[Path], str]:
    service, fingerprint = _agentic_service(spec)
    dataset_id = spec.payload["dataset_id"]
    science = ScienceStore(store)
    science.update_dataset(dataset_id, DatasetStatus.RUNNING)
    result = service.inspect_dataset(InspectionRequest(**spec.payload["request"]))
    destination = Path(spec.run_directory) / "inspection.json"
    destination.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    science.update_dataset(dataset_id, DatasetStatus.READY, result=result, error=None)
    return [destination], fingerprint


def _evaluate_recommendation(
    spec: JobSpec,
    guard: ProjectPaths,
) -> tuple[list[Path], str]:
    service, fingerprint = _agentic_service(spec)
    request = EvaluationRequest(**spec.payload["request"])
    service.evaluate_recommendation(request)
    return [guard.resolve_existing(request.output_dir)], fingerprint


def _provenance(
    spec: JobSpec,
    artifacts: list[ArtifactRecord],
    input_identities: list[dict],
    started_at: str,
    knowledge_fingerprint: str | None,
) -> dict:
    return {
        "run_id": spec.run_id,
        "command": spec.argv if spec.operation is JobOperation.NATIVE else [],
        "operation": spec.operation.value,
        "contract_version": 1,
        "mokume_version": mokume.__version__,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "knowledge_fingerprint": knowledge_fingerprint,
        "approved_hash": spec.approved_hash,
        "parameters": spec.parameters,
        "threads": spec.threads,
        "inputs": input_identities,
        "artifacts": [artifact.model_dump() for artifact in artifacts],
        "plan_source": spec.payload.get("plan_source"),
        "started_at": started_at,
        "finished_at": utc_now(),
        "status": RunStatus.SUCCEEDED.value,
    }


def main(argv: list[str] | None = None) -> int:
    """Load one spec by ID and report a stable terminal status."""
    args = build_parser().parse_args(argv)
    store = StateStore(args.state_directory)
    spec = None
    try:
        spec = _load_spec(store, args.run_id)
        return execute(store, spec)
    except KeyboardInterrupt:
        _fail_inspection(store, spec, "worker interrupted")
        record = store.get_run(args.run_id)
        if record and record.status is not RunStatus.CANCELLING:
            store.update_run(
                args.run_id, RunStatus.INTERRUPTED, error="worker interrupted"
            )
            write_terminal_files(store, args.run_id, record.run_directory)
        return 130
    except (KeyError, OSError, RuntimeError, TypeError, ValueError) as exc:
        _fail_inspection(store, spec, str(exc))
        store.update_run(args.run_id, RunStatus.FAILED, error=str(exc))
        record = store.get_run(args.run_id)
        if record:
            write_terminal_files(store, args.run_id, record.run_directory)
        print(f"Mokume Studio worker failed: {exc}", file=sys.stderr)
        return 1


def _fail_inspection(
    store: StateStore,
    spec: JobSpec | None,
    error: str,
) -> None:
    if spec is None or spec.operation is not JobOperation.INSPECT_DATASET:
        return
    dataset_id = spec.payload.get("dataset_id")
    if isinstance(dataset_id, str):
        ScienceStore(store).update_dataset(
            dataset_id,
            DatasetStatus.FAILED,
            error=error,
        )


def _load_spec(store: StateStore, run_id: str) -> JobSpec:
    path = store.spec_directory / f"{run_id}.json"
    try:
        spec = JobSpec.model_validate_json(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"cannot load job spec {run_id}") from exc
    if spec.run_id != run_id:
        raise RuntimeError("job spec ID mismatch")
    return spec


def _register_outputs(
    store: StateStore,
    run_id: str,
    guard: ProjectPaths,
    output_paths: list[Path],
) -> list[ArtifactRecord]:
    artifacts = []
    for requested in output_paths:
        try:
            output = guard.resolve_existing(requested)
        except ValueError:
            continue
        candidates = [output] if output.is_file() else sorted(output.rglob("*"))
        for path in (candidate for candidate in candidates if candidate.is_file()):
            media_type = (
                mimetypes.guess_type(path.name)[0] or "application/octet-stream"
            )
            artifact = ArtifactRecord(
                id=str(uuid.uuid4()),
                run_id=run_id,
                path=str(path),
                media_type=media_type,
                size=path.stat().st_size,
                sha256=_sha256(path),
            )
            store.register_artifact(artifact)
            store.add_event(run_id, "artifact", artifact.model_dump())
            artifacts.append(artifact)
    return artifacts


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


if __name__ == "__main__":
    raise SystemExit(main())
