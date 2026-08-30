"""Isolated worker entry point for one immutable Mokume Studio job."""

from __future__ import annotations

import argparse
import hashlib
import mimetypes
import os
import platform
import sys
import uuid
from pathlib import Path

import mokume

from mokume.studio.catalog import command_paths, validate_and_canonicalize
from mokume.studio.jobs import canonical_hash, path_snapshot, write_terminal_files
from mokume.studio.models import ArtifactRecord, JobSpec, RunStatus, utc_now
from mokume.studio.paths import ProjectPaths
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
    """Revalidate and execute a single Rust-backed command."""
    if canonical_hash(spec.parameters) != spec.approved_hash:
        raise RuntimeError("job parameters no longer match the approved hash")
    if spec.parameters.get("argv") != spec.argv:
        raise RuntimeError("job argv does not match immutable parameters")
    canonical = validate_and_canonicalize(spec.argv, spec.project_root)
    if canonical != spec.argv:
        raise RuntimeError("job argv changed during worker validation")
    guard = ProjectPaths(spec.project_root)
    run_directory = guard.resolve_existing(spec.run_directory)
    inputs, outputs = command_paths(spec.argv)
    resolved_inputs = [guard.resolve_existing(path) for path in inputs]
    input_snapshots = [path_snapshot(path, guard) for path in resolved_inputs]
    if input_snapshots != spec.parameters.get("input_snapshots"):
        raise RuntimeError("input files changed after the run was approved")
    input_identities = [file_identity(path, guard) for path in resolved_inputs]
    started_at = utc_now()
    store.update_run(spec.run_id, RunStatus.RUNNING, worker_pid=os.getpid())

    mokume.run(spec.argv)
    artifacts = _register_outputs(store, spec.run_id, guard, outputs)
    provenance = {
        "run_id": spec.run_id,
        "command": spec.argv,
        "contract_version": 1,
        "mokume_version": mokume.__version__,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "knowledge_fingerprint": None,
        "approved_hash": spec.approved_hash,
        "parameters": spec.parameters,
        "threads": spec.threads,
        "inputs": input_identities,
        "artifacts": [artifact.model_dump() for artifact in artifacts],
        "started_at": started_at,
        "finished_at": utc_now(),
        "status": RunStatus.SUCCEEDED.value,
    }
    store.update_run(spec.run_id, RunStatus.SUCCEEDED)
    write_terminal_files(store, spec.run_id, run_directory, provenance)
    return 0


def main(argv: list[str] | None = None) -> int:
    """Load one spec by ID and report a stable terminal status."""
    args = build_parser().parse_args(argv)
    store = StateStore(args.state_directory)
    try:
        spec = _load_spec(store, args.run_id)
        return execute(store, spec)
    except KeyboardInterrupt:
        record = store.get_run(args.run_id)
        if record and record.status is not RunStatus.CANCELLING:
            store.update_run(
                args.run_id, RunStatus.INTERRUPTED, error="worker interrupted"
            )
            write_terminal_files(store, args.run_id, record.run_directory)
        return 130
    except (KeyError, OSError, RuntimeError, TypeError, ValueError) as exc:
        store.update_run(args.run_id, RunStatus.FAILED, error=str(exc))
        record = store.get_run(args.run_id)
        if record:
            write_terminal_files(store, args.run_id, record.run_directory)
        print(f"Mokume Studio worker failed: {exc}", file=sys.stderr)
        return 1


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
