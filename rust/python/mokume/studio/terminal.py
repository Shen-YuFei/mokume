"""Atomic terminal record publication for Mokume Studio runs."""

from __future__ import annotations

import json
import uuid
from dataclasses import dataclass
from pathlib import Path

from mokume.studio.models import (
    JobOperation,
    RunRecord,
    RunStatus,
    TERMINAL_RUN_STATUSES,
)
from mokume.studio.state import StateStore


@dataclass(frozen=True)
class TerminalResult:
    """Terminal status and optional outcome metadata."""

    status: RunStatus
    error: str | None = None
    provenance: dict | None = None


def write_terminal_files(
    store: StateStore,
    run_id: str,
    run_directory: str | Path,
    provenance: dict | None = None,
) -> None:
    """Refresh terminal files for an already published run."""
    record = store.get_run(run_id)
    if record is None or record.status not in TERMINAL_RUN_STATUSES:
        return
    directory = Path(run_directory)
    if not directory.is_dir():
        return
    parameters = _read_parameters(directory)
    terminal_provenance = provenance
    if terminal_provenance is None:
        terminal_provenance = _read_json(directory / "provenance.json")
    if terminal_provenance is None:
        terminal_provenance = _fallback_provenance(
            record,
            parameters,
            _artifact_payloads(store, run_id),
        )
    _write_terminal_snapshot(
        directory,
        record,
        store.events_after(run_id),
        terminal_provenance,
    )


def finalize_run(
    store: StateStore,
    run_id: str,
    run_directory: str | Path,
    result: TerminalResult,
) -> RunRecord:
    """Publish terminal files before making the terminal state observable."""
    directory = Path(run_directory)
    parameters = _read_parameters(directory)
    artifacts = _artifact_payloads(store, run_id)

    def publish(record: RunRecord, events: list[dict]) -> None:
        provenance = result.provenance
        if provenance is None:
            provenance = _fallback_provenance(record, parameters, artifacts)
        _write_terminal_snapshot(directory, record, events, provenance)

    return store.finalize_run(
        run_id,
        result.status,
        publish,
        error=result.error,
    )


def _read_parameters(run_directory: Path) -> dict:
    return _read_json(run_directory / "parameters.json") or {}


def _read_json(path: Path) -> dict | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return payload if isinstance(payload, dict) else None


def _artifact_payloads(store: StateStore, run_id: str) -> list[dict]:
    return [
        artifact.model_dump(mode="json") for artifact in store.list_artifacts(run_id)
    ]


def _fallback_provenance(
    record: RunRecord,
    parameters: dict,
    artifacts: list[dict],
) -> dict:
    return {
        "run_id": record.id,
        "command": record.argv,
        "operation": parameters.get("operation", JobOperation.NATIVE.value),
        "contract_version": 1,
        "knowledge_fingerprint": parameters.get("payload", {}).get(
            "knowledge_fingerprint"
        ),
        "approved_hash": record.approved_hash,
        "parameters": parameters,
        "threads": parameters.get("threads", 24),
        "inputs": parameters.get("input_snapshots", []),
        "artifacts": artifacts,
        "plan_source": parameters.get("payload", {}).get("plan_source"),
        "started_at": record.started_at,
        "finished_at": record.finished_at,
        "status": record.status.value,
        "error": record.error,
    }


def _write_terminal_snapshot(
    directory: Path,
    record: RunRecord,
    events: list[dict],
    provenance: dict,
) -> None:
    _write_json_atomic(directory / "provenance.json", provenance)
    _write_json_atomic(directory / "run.json", record.model_dump(mode="json"))
    temporary = directory / f".events.jsonl.{uuid.uuid4().hex}.tmp"
    with temporary.open("w", encoding="utf-8") as stream:
        for event in events:
            stream.write(json.dumps(event, sort_keys=True) + "\n")
    temporary.replace(directory / "events.jsonl")


def _write_json_atomic(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)
