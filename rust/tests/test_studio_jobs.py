"""Worker lifecycle contracts for Mokume Studio."""

from __future__ import annotations

import json
import os
import shutil
import time
from pathlib import Path

from mokume.studio.jobs import JobManager
from mokume.studio.models import RunRequest, RunStatus, TERMINAL_RUN_STATUSES
from mokume.studio.state import StateStore


def _project_with_peptides(tmp_path: Path) -> tuple[StateStore, Path]:
    """Create a guarded project with a real peptide table."""
    project_root = tmp_path / "project"
    project_root.mkdir()
    fixture = Path(__file__).parent / "example" / "peptides_small.csv"
    shutil.copy2(fixture, project_root / "peptides small.csv")
    store = StateStore(tmp_path / "state")
    store.open_project(str(project_root.resolve()))
    return store, project_root


def _request() -> RunRequest:
    """Build a real sum-aggregation request with paths containing spaces."""
    return RunRequest(
        argv=[
            "quantify",
            "peptides2protein",
            "--peptides",
            "peptides small.csv",
            "--quant-method",
            "sum",
            "--output",
            "protein result.tsv",
        ]
    )


def _wait_for_terminal(store: StateStore, run_id: str, timeout: float = 30):
    """Wait for a worker to reach a persisted terminal state."""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        record = store.get_run(run_id)
        if record and record.status in TERMINAL_RUN_STATUSES:
            return record
        time.sleep(0.05)
    raise TimeoutError(f"Studio run did not finish: {run_id}")


def _sleep_until_cancelled(ready) -> None:
    """Enter a process group and wait until the cancellation test stops it."""
    if os.name != "nt":
        os.setsid()
    ready.set()
    time.sleep(60)


def _assert_stage_events(store: StateStore, record) -> None:
    """Verify live and persisted worker stage events."""
    run_directory = Path(record.run_directory)
    stages = [
        event["payload"]
        for event in store.events_after(record.id)
        if event["type"] == "stage"
    ]
    assert [(stage["stage"], stage["status"]) for stage in stages] == [
        ("inputs", "running"),
        ("inputs", "succeeded"),
        ("workflow", "running"),
        ("workflow", "succeeded"),
        ("artifacts", "running"),
        ("artifacts", "succeeded"),
        ("provenance", "running"),
        ("provenance", "succeeded"),
    ]
    persisted_events = [
        json.loads(line)
        for line in (run_directory / "events.jsonl")
        .read_text(encoding="utf-8")
        .splitlines()
    ]
    assert any(
        event["type"] == "stage"
        and event["payload"]
        == {
            "elapsed_seconds": stages[-1]["elapsed_seconds"],
            "stage": "provenance",
            "status": "succeeded",
        }
        for event in persisted_events
    )


def test_real_worker_writes_output_artifact_and_provenance(tmp_path):
    """A real worker records output, input identity, provenance, and events."""
    store, project_root = _project_with_peptides(tmp_path)
    manager = JobManager(store)
    project = store.active_project()

    try:
        submitted = manager.submit(_request(), project)
        record = _wait_for_terminal(store, submitted.id)
    finally:
        manager.shutdown()

    assert record.status is RunStatus.SUCCEEDED
    assert (project_root / "protein result.tsv").is_file()
    provenance_path = Path(record.run_directory) / "provenance.json"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    assert provenance["approved_hash"] == record.approved_hash
    assert provenance["threads"] == 24
    assert len(provenance["inputs"]) == 1
    assert provenance["inputs"][0]["sha256"]
    assert len(store.list_artifacts(record.id)) == 1
    run_directory = Path(record.run_directory)
    assert json.loads((run_directory / "run.json").read_text())["status"] == "succeeded"
    assert (run_directory / "events.jsonl").read_text(encoding="utf-8").strip()
    _assert_stage_events(store, record)


class _SleepingJobManager(JobManager):
    """Replace the worker with a process that only waits for cancellation."""

    def _spawn(self, _spec):
        """Start a cancellable process group without invoking Mokume."""
        ready = self._context.Event()
        process = self._context.Process(target=_sleep_until_cancelled, args=(ready,))
        process.start()
        assert ready.wait(timeout=5)
        return process


def test_cancel_waits_for_worker_process_exit(tmp_path):
    """Cancellation waits for process exit and writes terminal records."""
    store, _project_root = _project_with_peptides(tmp_path)
    manager = _SleepingJobManager(store)
    project = store.active_project()
    submitted = manager.submit(_request(), project)
    process = getattr(manager, "_processes")[submitted.id]

    record = manager.cancel(submitted.id)

    assert not process.is_alive()
    assert record.status is RunStatus.CANCELLED
    provenance = json.loads(
        (Path(record.run_directory) / "provenance.json").read_text(encoding="utf-8")
    )
    assert provenance["status"] == "cancelled"
    assert (Path(record.run_directory) / "events.jsonl").is_file()


class _MutatingJobManager(JobManager):
    """Mutate a snapshotted input immediately before the worker starts."""

    def _spawn(self, spec):
        """Change the input after approval and then start the real worker."""
        input_path = Path(spec.parameters["input_snapshots"][0]["path"])
        input_path.write_text("changed after approval\n", encoding="utf-8")
        return super()._spawn(spec)


def test_worker_rejects_input_changed_after_approval(tmp_path):
    """The worker refuses an input whose metadata no longer matches approval."""
    store, project_root = _project_with_peptides(tmp_path)
    manager = _MutatingJobManager(store)
    project = store.active_project()

    try:
        submitted = manager.submit(_request(), project)
        record = _wait_for_terminal(store, submitted.id)
    finally:
        manager.shutdown()

    assert record.status is RunStatus.FAILED
    assert "changed after the run was approved" in (record.error or "")
    assert not (project_root / "protein result.tsv").exists()
    run_directory = Path(record.run_directory)
    assert json.loads((run_directory / "run.json").read_text())["status"] == "failed"
    assert (run_directory / "events.jsonl").is_file()
