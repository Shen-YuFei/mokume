"""Single-worker process manager for reproducible Mokume Studio runs."""

from __future__ import annotations

import hashlib
import importlib
import json
import multiprocessing
import os
import signal
import threading
import time
import uuid
from collections.abc import Callable
from multiprocessing.process import BaseProcess
from pathlib import Path
from typing import Any

from mokume.studio.catalog import command_paths, validate_and_canonicalize
from mokume.studio.models import (
    JobOperation,
    JobSpec,
    ProjectRecord,
    RunRecord,
    RunRequest,
    RunStatus,
    ScientificJobRequest,
    TERMINAL_RUN_STATUSES,
    utc_now,
)
from mokume.studio.paths import PathAccessError, ProjectPaths
from mokume.studio.state import StateStore


class JobConflictError(RuntimeError):
    """Raised when the single heavy-worker slot is already occupied."""


def canonical_hash(payload: dict) -> str:
    """Hash stable JSON used by approval and worker verification."""
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def validate_spec_integrity(spec: JobSpec) -> None:
    """Verify the immutable fields shared by manager and worker boundaries."""
    if canonical_hash(spec.parameters) != spec.approved_hash:
        raise RuntimeError("job parameters no longer match the approved hash")
    if spec.parameters.get("operation") != spec.operation.value:
        raise RuntimeError("job operation does not match immutable parameters")
    if spec.parameters.get("payload") != spec.payload:
        raise RuntimeError("job payload does not match immutable parameters")
    if spec.operation is JobOperation.NATIVE:
        if spec.parameters.get("argv") != spec.argv:
            raise RuntimeError("job argv does not match immutable parameters")
    elif spec.argv:
        raise RuntimeError("scientific jobs cannot carry native argv")


def path_snapshot(path: Path, guard: ProjectPaths) -> dict:
    """Capture input metadata without hashing file contents in the web process."""
    stat = path.stat()
    snapshot = {
        "path": str(path),
        "kind": "directory" if path.is_dir() else "file",
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }
    if path.is_dir():
        entries = []
        for candidate in sorted(path.rglob("*")):
            resolved = guard.resolve_existing(candidate)
            if resolved.is_file():
                item_stat = resolved.stat()
                entries.append(
                    {
                        "path": str(candidate.relative_to(path)),
                        "resolved_path": str(resolved),
                        "size": item_stat.st_size,
                        "mtime_ns": item_stat.st_mtime_ns,
                    }
                )
        snapshot["entries"] = entries
    return snapshot


def write_json_atomic(path: Path, payload: dict) -> None:
    """Replace one JSON record without exposing a partially written file."""
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def write_terminal_files(
    store: StateStore,
    run_id: str,
    run_directory: str | Path,
    provenance: dict | None = None,
) -> None:
    """Persist terminal run, provenance, and event records for every outcome."""
    record = store.get_run(run_id)
    if record is None or record.status not in TERMINAL_RUN_STATUSES:
        return
    directory = Path(run_directory)
    if not directory.is_dir():
        return
    parameters = _read_parameters(directory)
    provenance_path = directory / "provenance.json"
    if provenance is not None:
        write_json_atomic(provenance_path, provenance)
    elif not provenance_path.exists():
        write_json_atomic(
            provenance_path,
            _fallback_provenance(store, record, parameters),
        )
    write_json_atomic(directory / "run.json", record.model_dump(mode="json"))
    _write_events_atomic(directory / "events.jsonl", store.events_after(run_id))


def _read_parameters(run_directory: Path) -> dict:
    try:
        return json.loads(
            (run_directory / "parameters.json").read_text(encoding="utf-8")
        )
    except (OSError, ValueError):
        return {}


def _fallback_provenance(
    store: StateStore, record: RunRecord, parameters: dict
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
        "artifacts": [
            artifact.model_dump(mode="json")
            for artifact in store.list_artifacts(record.id)
        ],
        "plan_source": parameters.get("payload", {}).get("plan_source"),
        "started_at": record.started_at,
        "finished_at": record.finished_at,
        "status": record.status.value,
        "error": record.error,
    }


def _write_events_atomic(path: Path, events: list[dict]) -> None:
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        for event in events:
            stream.write(json.dumps(event, sort_keys=True) + "\n")
    temporary.replace(path)


def _run_worker_process(
    state_directory: str,
    run_id: str,
    run_directory: str,
    project_root: str,
    ready: Any,
) -> None:
    """Enter an isolated process group and execute the private worker entry point."""
    if os.name != "nt":
        os.setsid()
    ready.set()
    _redirect_process_streams(Path(run_directory))
    os.chdir(project_root)
    worker = importlib.import_module("mokume.studio.worker")
    exit_code = worker.main(["--state-directory", state_directory, "--run-id", run_id])
    raise SystemExit(exit_code)


def _redirect_process_streams(run_directory: Path) -> None:
    null_fd = os.open(os.devnull, os.O_RDONLY)
    stdout_fd = os.open(
        run_directory / "stdout.log", os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600
    )
    stderr_fd = os.open(
        run_directory / "stderr.log", os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600
    )
    try:
        os.dup2(null_fd, 0)
        os.dup2(stdout_fd, 1)
        os.dup2(stderr_fd, 2)
    finally:
        os.close(null_fd)
        os.close(stdout_fd)
        os.close(stderr_fd)


class JobManager:
    """Spawn, monitor, and truly cancel at most one heavy run."""

    def __init__(self, store: StateStore) -> None:
        self.store = store
        self._lock = threading.Lock()
        self._context = multiprocessing.get_context("spawn")
        self._processes: dict[str, BaseProcess] = {}
        self._recover_terminal_files()

    def submit(self, request: RunRequest, project: ProjectRecord) -> RunRecord:
        """Validate, persist, and spawn one immutable run specification."""
        spec = self._prepare_spec(request, project.root)
        return self.submit_spec(spec, project, " ".join(self._command_path(spec.argv)))

    def prepare_scientific_spec(
        self,
        request: ScientificJobRequest,
        project: ProjectRecord,
    ) -> JobSpec:
        """Build an unmaterialized, hash-bound scientific worker spec."""
        if request.operation is JobOperation.NATIVE:
            raise ValueError("scientific specs require a scientific operation")
        guard = ProjectPaths(project.root)
        identifier = str(uuid.uuid4())
        base_directory = guard.resolve_output(
            request.output_directory, allow_existing=True
        )
        if base_directory.exists() and not base_directory.is_dir():
            raise PathAccessError(f"run output is not a directory: {base_directory}")
        run_directory = base_directory / identifier
        if run_directory.exists():
            raise PathAccessError(f"run directory already exists: {run_directory}")
        snapshots = [
            path_snapshot(guard.resolve_existing(path), guard)
            for path in request.input_paths
        ]
        parameters = {
            "operation": request.operation.value,
            "payload": request.payload,
            "project_root": str(guard.root),
            "run_directory": str(run_directory),
            "threads": 24,
            "input_snapshots": snapshots,
        }
        return JobSpec(
            run_id=identifier,
            project_root=str(guard.root),
            run_directory=str(run_directory),
            operation=request.operation,
            argv=[],
            payload=request.payload,
            parameters=parameters,
            approved_hash=canonical_hash(parameters),
            created_at=utc_now(),
        )

    def submit_spec(
        self,
        spec: JobSpec,
        project: ProjectRecord,
        command: str,
        *,
        authorize: Callable[[], None] | None = None,
    ) -> RunRecord:
        """Revalidate and spawn a trusted server-created specification."""
        with self._lock:
            self._discard_finished_locked()
            if self._processes:
                raise JobConflictError("another Mokume run is already active")
            self._validate_spec(spec, project)
            if authorize is not None:
                authorize()
            self._materialize_spec(spec)
            return self._start_locked(spec, project, command)

    def _prepare_spec(self, request: RunRequest, project_root: str) -> JobSpec:
        canonical = validate_and_canonicalize(request.argv, project_root)
        guard = ProjectPaths(project_root)
        run_id = str(uuid.uuid4())
        base_directory = guard.resolve_output(
            request.output_directory, allow_existing=True
        )
        if base_directory.exists() and not base_directory.is_dir():
            raise PathAccessError(f"run output is not a directory: {base_directory}")
        run_directory = base_directory / run_id
        if run_directory.exists():  # practically impossible, kept deterministic
            raise PathAccessError(f"run directory already exists: {run_directory}")
        inputs, _outputs = command_paths(canonical)
        input_snapshots = [
            path_snapshot(guard.resolve_existing(path), guard) for path in inputs
        ]
        parameters = {
            "operation": JobOperation.NATIVE.value,
            "payload": {},
            "argv": canonical,
            "project_root": str(guard.root),
            "run_directory": str(run_directory),
            "threads": 24,
            "input_snapshots": input_snapshots,
        }
        return JobSpec(
            run_id=run_id,
            project_root=str(guard.root),
            run_directory=str(run_directory),
            operation=JobOperation.NATIVE,
            argv=canonical,
            payload={},
            parameters=parameters,
            approved_hash=canonical_hash(parameters),
            created_at=utc_now(),
        )

    def _validate_spec(self, spec: JobSpec, project: ProjectRecord) -> None:
        if spec.project_root != project.root:
            raise RuntimeError("job project no longer matches the active project")
        validate_spec_integrity(spec)
        guard = ProjectPaths(spec.project_root)
        snapshots = [
            path_snapshot(guard.resolve_existing(item["path"]), guard)
            for item in spec.parameters.get("input_snapshots", [])
        ]
        if snapshots != spec.parameters.get("input_snapshots"):
            raise RuntimeError("input files changed after the run was approved")

    def _materialize_spec(self, spec: JobSpec) -> None:
        guard = ProjectPaths(spec.project_root)
        directory = guard.resolve_output(spec.run_directory)
        directory.mkdir(parents=True)
        if spec.operation is JobOperation.NATIVE:
            self._prepare_command_outputs(spec.argv, guard)

    def _start_locked(
        self,
        spec: JobSpec,
        project: ProjectRecord,
        command: str,
    ) -> RunRecord:
        self._write_spec(spec)
        record = self.store.create_run(spec, project.id, command)
        self.store.update_run(spec.run_id, RunStatus.STARTING)
        try:
            process = self._spawn(spec)
            worker_pid = process.pid
            if worker_pid is None:
                raise RuntimeError("worker process did not expose a process ID")
        except (OSError, RuntimeError) as exc:
            self.store.update_run(spec.run_id, RunStatus.FAILED, error=str(exc))
            write_terminal_files(self.store, spec.run_id, spec.run_directory)
            raise
        self._processes[spec.run_id] = process
        self.store.set_worker_pid(spec.run_id, worker_pid)
        threading.Thread(
            target=self._monitor,
            args=(spec, process),
            daemon=True,
        ).start()
        return record.model_copy(
            update={"status": RunStatus.STARTING, "worker_pid": worker_pid}
        )

    def cancel(self, run_id: str) -> RunRecord:
        """Interrupt, terminate, then kill the worker group if necessary."""
        with self._lock:
            process = self._processes.get(run_id)
        record = self.store.get_run(run_id)
        if record is None:
            raise KeyError(run_id)
        if record.status in TERMINAL_RUN_STATUSES:
            return record
        if process is None or not process.is_alive():
            return record
        self.store.update_run(run_id, RunStatus.CANCELLING)
        try:
            self._interrupt(process)
        except ProcessLookupError:
            pass
        if not self._wait(process, 5):
            process.terminate()
        if not self._wait(process, 3):
            process.kill()
            process.join(timeout=3)
        self.store.update_run(run_id, RunStatus.CANCELLED)
        terminal = self.store.get_run(run_id) or record
        write_terminal_files(self.store, run_id, terminal.run_directory)
        return terminal

    def shutdown(self) -> None:
        """Cancel every active process during a graceful Studio shutdown."""
        with self._lock:
            run_ids = list(self._processes)
        for run_id in run_ids:
            try:
                self.cancel(run_id)
            except (KeyError, OSError, RuntimeError):
                self.store.update_run(
                    run_id,
                    RunStatus.INTERRUPTED,
                    error="Studio stopped during cancellation",
                )
                record = self.store.get_run(run_id)
                if record:
                    write_terminal_files(self.store, run_id, record.run_directory)

    def _write_spec(self, spec: JobSpec) -> None:
        path = self.store.spec_directory / f"{spec.run_id}.json"
        with path.open("x", encoding="utf-8") as stream:
            stream.write(spec.model_dump_json(indent=2))
            stream.write("\n")
        try:
            path.chmod(0o600)
        except OSError:
            pass
        run_directory = Path(spec.run_directory)
        (run_directory / "parameters.json").write_text(
            json.dumps(spec.parameters, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    def _spawn(self, spec: JobSpec) -> BaseProcess:
        ready = self._context.Event()
        process = self._context.Process(
            target=_run_worker_process,
            args=(
                str(self.store.directory),
                spec.run_id,
                spec.run_directory,
                spec.project_root,
                ready,
            ),
            name=f"mokume-studio-{spec.run_id}",
        )
        process.start()
        if not ready.wait(timeout=5):
            if process.is_alive():
                process.terminate()
            process.join(timeout=3)
            raise RuntimeError("worker process did not finish startup")
        return process

    def _monitor(self, spec: JobSpec, process: BaseProcess) -> None:
        run_id = spec.run_id
        run_directory = Path(spec.run_directory)
        stdout = run_directory / "stdout.log"
        stderr = run_directory / "stderr.log"
        offsets = {stdout: 0, stderr: 0}
        while process.is_alive():
            self._append_log_events(run_id, offsets)
            time.sleep(0.2)
        self._append_log_events(run_id, offsets)
        record = self.store.get_run(run_id)
        if record and record.status not in {
            RunStatus.CANCELLED,
            RunStatus.SUCCEEDED,
            RunStatus.FAILED,
        }:
            status = (
                RunStatus.CANCELLED
                if record.status is RunStatus.CANCELLING
                else RunStatus.FAILED
            )
            error = (
                None
                if status is RunStatus.CANCELLED
                else f"worker exited with code {process.exitcode}"
            )
            self.store.update_run(run_id, status, error=error)
        write_terminal_files(self.store, run_id, run_directory)
        with self._lock:
            self._processes.pop(run_id, None)

    def _append_log_events(self, run_id: str, offsets: dict[Path, int]) -> None:
        for path, offset in offsets.items():
            if not path.exists():
                continue
            with path.open("r", encoding="utf-8", errors="replace") as stream:
                stream.seek(offset)
                for line in stream:
                    self.store.add_event(
                        run_id,
                        "log",
                        {"stream": path.stem, "line": line.rstrip()[:8192]},
                    )
                offsets[path] = stream.tell()

    def _discard_finished_locked(self) -> None:
        finished = [
            run_id
            for run_id, process in self._processes.items()
            if not process.is_alive()
        ]
        for run_id in finished:
            self._processes.pop(run_id, None)

    def _recover_terminal_files(self) -> None:
        for record in self.store.list_runs():
            if record.status in TERMINAL_RUN_STATUSES:
                write_terminal_files(self.store, record.id, record.run_directory)

    @staticmethod
    def _prepare_command_outputs(argv: list[str], guard: ProjectPaths) -> None:
        _, outputs = command_paths(argv)
        for output in outputs:
            guarded = guard.resolve_output(output)
            guarded.parent.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def _command_path(argv: list[str]) -> list[str]:
        return argv[:2] if argv[:1] == ["quantify"] else argv[:1]

    @staticmethod
    def _interrupt(process: BaseProcess) -> None:
        if os.name == "nt":
            process.terminate()
        else:
            os.killpg(process.pid, signal.SIGINT)

    @staticmethod
    def _wait(process: BaseProcess, seconds: float) -> bool:
        process.join(timeout=seconds)
        return not process.is_alive()
