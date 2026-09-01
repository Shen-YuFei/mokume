"""SQLite-backed control state for local Mokume Studio sessions."""

from __future__ import annotations

import json
import sqlite3
import uuid
from collections.abc import Callable
from pathlib import Path

from platformdirs import user_state_dir

from mokume.studio.models import (
    ArtifactRecord,
    JobSpec,
    ProjectRecord,
    RunRecord,
    RunStatus,
    TERMINAL_RUN_STATUSES,
    utc_now,
)


SCHEMA = """
PRAGMA journal_mode = WAL;
PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS projects (
    id TEXT PRIMARY KEY,
    root TEXT NOT NULL UNIQUE,
    opened_at TEXT NOT NULL,
    active INTEGER NOT NULL DEFAULT 0
);

CREATE TABLE IF NOT EXISTS runs (
    id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL REFERENCES projects(id),
    status TEXT NOT NULL,
    command TEXT NOT NULL,
    argv_json TEXT NOT NULL,
    approved_hash TEXT NOT NULL,
    run_directory TEXT NOT NULL,
    worker_pid INTEGER,
    created_at TEXT NOT NULL,
    started_at TEXT,
    finished_at TEXT,
    error TEXT
);

CREATE TABLE IF NOT EXISTS run_events (
    sequence INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id TEXT NOT NULL REFERENCES runs(id),
    created_at TEXT NOT NULL,
    event_type TEXT NOT NULL,
    payload_json TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS artifacts (
    id TEXT PRIMARY KEY,
    run_id TEXT NOT NULL REFERENCES runs(id),
    path TEXT NOT NULL,
    media_type TEXT NOT NULL,
    size INTEGER NOT NULL,
    sha256 TEXT NOT NULL,
    UNIQUE(run_id, path)
);

CREATE TABLE IF NOT EXISTS approvals (
    id TEXT PRIMARY KEY,
    run_id TEXT,
    kind TEXT NOT NULL,
    payload_json TEXT NOT NULL,
    payload_hash TEXT NOT NULL,
    status TEXT NOT NULL,
    created_at TEXT NOT NULL,
    decided_at TEXT
);
"""


class StateStore:
    """Small synchronous store shared by the server and worker processes."""

    def __init__(self, directory: str | Path | None = None) -> None:
        root = Path(directory or user_state_dir("mokume", appauthor=False))
        root.mkdir(parents=True, exist_ok=True)
        self.directory = root.resolve()
        self.database = self.directory / "studio.sqlite3"
        self.spec_directory = self.directory / "jobs"
        self.spec_directory.mkdir(exist_ok=True)
        with self._connect() as connection:
            connection.executescript(SCHEMA)

    def _connect(self) -> sqlite3.Connection:
        connection = sqlite3.connect(self.database, timeout=30)
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA foreign_keys = ON")
        return connection

    def open_project(self, root: str) -> ProjectRecord:
        """Register and activate one canonical project root."""
        project_id = str(uuid.uuid5(uuid.NAMESPACE_URL, f"mokume-project:{root}"))
        opened_at = utc_now()
        with self._connect() as connection:
            connection.execute("UPDATE projects SET active = 0 WHERE active = 1")
            connection.execute(
                """
                INSERT INTO projects(id, root, opened_at, active)
                VALUES (?, ?, ?, 1)
                ON CONFLICT(root) DO UPDATE SET opened_at=excluded.opened_at, active=1
                """,
                (project_id, root, opened_at),
            )
        return ProjectRecord(id=project_id, root=root, opened_at=opened_at)

    def active_project(self) -> ProjectRecord | None:
        """Return the currently selected project, if any."""
        with self._connect() as connection:
            row = connection.execute(
                "SELECT id, root, opened_at FROM projects WHERE active = 1"
            ).fetchone()
        return ProjectRecord(**dict(row)) if row else None

    def get_project(self, project_id: str) -> ProjectRecord | None:
        """Return a previously registered project by opaque ID."""
        with self._connect() as connection:
            row = connection.execute(
                "SELECT id, root, opened_at FROM projects WHERE id=?", (project_id,)
            ).fetchone()
        return ProjectRecord(**dict(row)) if row else None

    def close_project(self) -> None:
        """Clear the active-project pointer without deleting run history."""
        with self._connect() as connection:
            connection.execute("UPDATE projects SET active = 0 WHERE active = 1")

    def create_run(
        self,
        spec: JobSpec,
        project_id: str,
        command: str,
    ) -> RunRecord:
        """Persist a queued run before any child process is created."""
        created_at = utc_now()
        with self._connect() as connection:
            connection.execute(
                """
                INSERT INTO runs(
                    id, project_id, status, command, argv_json, approved_hash,
                    run_directory, created_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    spec.run_id,
                    project_id,
                    RunStatus.QUEUED.value,
                    command,
                    json.dumps(spec.argv),
                    spec.approved_hash,
                    spec.run_directory,
                    created_at,
                ),
            )
        self.add_event(spec.run_id, "status", {"status": RunStatus.QUEUED.value})
        record = self.get_run(spec.run_id)
        if record is None:  # pragma: no cover - protects against storage corruption
            raise RuntimeError(f"run was not persisted: {spec.run_id}")
        return record

    def update_run(
        self,
        run_id: str,
        status: RunStatus,
        *,
        worker_pid: int | None = None,
        error: str | None = None,
    ) -> None:
        """Update lifecycle timestamps together with the requested status."""
        started_at = utc_now() if status is RunStatus.RUNNING else None
        finished_at = (
            utc_now()
            if status
            in {
                RunStatus.CANCELLED,
                RunStatus.SUCCEEDED,
                RunStatus.FAILED,
                RunStatus.INTERRUPTED,
            }
            else None
        )
        with self._connect() as connection:
            connection.execute(
                """
                UPDATE runs SET status=?, worker_pid=COALESCE(?, worker_pid),
                    started_at=COALESCE(?, started_at),
                    finished_at=COALESCE(?, finished_at), error=COALESCE(?, error)
                WHERE id=?
                """,
                (
                    status.value,
                    worker_pid,
                    started_at,
                    finished_at,
                    error,
                    run_id,
                ),
            )
        payload = {"status": status.value}
        if error:
            payload["error"] = error
        self.add_event(run_id, "status", payload)

    def finalize_run(
        self,
        run_id: str,
        status: RunStatus,
        publish: Callable[[RunRecord, list[dict]], None],
        *,
        error: str | None = None,
    ) -> RunRecord:
        """Publish terminal files before exposing the terminal database state."""
        if status not in TERMINAL_RUN_STATUSES:
            raise ValueError(f"run status is not terminal: {status.value}")
        payload = {"status": status.value, **({"error": error} if error else {})}
        with self._connect() as connection:
            connection.execute("BEGIN IMMEDIATE")
            current = connection.execute(
                "SELECT * FROM runs WHERE id = ?", (run_id,)
            ).fetchone()
            if current is None:
                raise KeyError(run_id)
            current_record = self._decode_run(current)
            if current_record.status in TERMINAL_RUN_STATUSES:
                return current_record
            connection.execute(
                """
                UPDATE runs SET status=?, finished_at=?, error=COALESCE(?, error)
                WHERE id=?
                """,
                (status.value, utc_now(), error, run_id),
            )
            connection.execute(
                """
                INSERT INTO run_events(run_id, created_at, event_type, payload_json)
                VALUES (?, ?, ?, ?)
                """,
                (run_id, utc_now(), "status", json.dumps(payload, sort_keys=True)),
            )
            row = connection.execute(
                "SELECT * FROM runs WHERE id = ?", (run_id,)
            ).fetchone()
            event_rows = connection.execute(
                """
                SELECT sequence, created_at, event_type, payload_json
                FROM run_events WHERE run_id=? ORDER BY sequence
                """,
                (run_id,),
            ).fetchall()
            record = self._decode_run(row)
            publish(record, self._decode_events(event_rows))
        return record

    def set_worker_pid(self, run_id: str, worker_pid: int) -> None:
        """Attach the spawned process identity without changing lifecycle state."""
        with self._connect() as connection:
            connection.execute(
                "UPDATE runs SET worker_pid=? WHERE id=?", (worker_pid, run_id)
            )

    def get_run(self, run_id: str) -> RunRecord | None:
        """Return one run by opaque ID."""
        with self._connect() as connection:
            row = connection.execute(
                "SELECT * FROM runs WHERE id = ?", (run_id,)
            ).fetchone()
        return self._decode_run(row) if row else None

    def list_runs(
        self,
        limit: int = 100,
        *,
        project_id: str | None = None,
    ) -> list[RunRecord]:
        """Return newest runs, optionally restricted to one project."""
        query = "SELECT * FROM runs"
        parameters: tuple[object, ...]
        if project_id is None:
            parameters = (limit,)
        else:
            query += " WHERE project_id=?"
            parameters = (project_id, limit)
        query += " ORDER BY created_at DESC LIMIT ?"
        with self._connect() as connection:
            rows = connection.execute(query, parameters).fetchall()
        return [self._decode_run(row) for row in rows]

    def has_active_run(self) -> bool:
        """Return whether a run is still able to read the active project."""
        statuses = (
            RunStatus.QUEUED.value,
            RunStatus.STARTING.value,
            RunStatus.RUNNING.value,
            RunStatus.CANCELLING.value,
        )
        with self._connect() as connection:
            row = connection.execute(
                "SELECT 1 FROM runs WHERE status IN (?, ?, ?, ?) LIMIT 1", statuses
            ).fetchone()
        return row is not None

    def add_event(self, run_id: str, event_type: str, payload: dict) -> int:
        """Append one resumable run event and return its sequence number."""
        with self._connect() as connection:
            cursor = connection.execute(
                """
                INSERT INTO run_events(run_id, created_at, event_type, payload_json)
                VALUES (?, ?, ?, ?)
                """,
                (run_id, utc_now(), event_type, json.dumps(payload, sort_keys=True)),
            )
            return int(cursor.lastrowid)

    def events_after(self, run_id: str, sequence: int = 0) -> list[dict]:
        """Read run events after a reconnect cursor."""
        with self._connect() as connection:
            rows = connection.execute(
                """
                SELECT sequence, created_at, event_type, payload_json
                FROM run_events WHERE run_id=? AND sequence>?
                ORDER BY sequence
                """,
                (run_id, sequence),
            ).fetchall()
        return self._decode_events(rows)

    @staticmethod
    def _decode_events(rows: list[sqlite3.Row]) -> list[dict]:
        return [
            {
                "sequence": row["sequence"],
                "created_at": row["created_at"],
                "type": row["event_type"],
                "payload": json.loads(row["payload_json"]),
            }
            for row in rows
        ]

    def interrupt_incomplete_runs(self) -> int:
        """Mark runs left active by an earlier server as interrupted."""
        statuses = (
            RunStatus.QUEUED.value,
            RunStatus.STARTING.value,
            RunStatus.RUNNING.value,
            RunStatus.CANCELLING.value,
        )
        with self._connect() as connection:
            rows = connection.execute(
                "SELECT id FROM runs WHERE status IN (?, ?, ?, ?)", statuses
            ).fetchall()
            cursor = connection.execute(
                """UPDATE runs SET status=?, finished_at=?, error=?
                WHERE status IN (?, ?, ?, ?)""",
                (
                    RunStatus.INTERRUPTED.value,
                    utc_now(),
                    "Studio stopped before the run completed",
                    *statuses,
                ),
            )
        for row in rows:
            self.add_event(
                row["id"],
                "status",
                {
                    "status": RunStatus.INTERRUPTED.value,
                    "error": "Studio stopped before the run completed",
                },
            )
        return cursor.rowcount

    def register_artifact(self, artifact: ArtifactRecord) -> None:
        """Persist a worker-verified output artifact."""
        with self._connect() as connection:
            connection.execute(
                """
                INSERT INTO artifacts(id, run_id, path, media_type, size, sha256)
                VALUES (?, ?, ?, ?, ?, ?)
                """,
                (
                    artifact.id,
                    artifact.run_id,
                    artifact.path,
                    artifact.media_type,
                    artifact.size,
                    artifact.sha256,
                ),
            )

    def get_artifact(self, artifact_id: str) -> ArtifactRecord | None:
        """Resolve an opaque artifact ID without accepting a filesystem path."""
        with self._connect() as connection:
            row = connection.execute(
                "SELECT * FROM artifacts WHERE id=?", (artifact_id,)
            ).fetchone()
        return ArtifactRecord(**dict(row)) if row else None

    def list_artifacts(self, run_id: str | None = None) -> list[ArtifactRecord]:
        """List registered outputs, optionally restricted to one run."""
        query = "SELECT * FROM artifacts"
        parameters: tuple[str, ...] = ()
        if run_id is not None:
            query += " WHERE run_id=?"
            parameters = (run_id,)
        query += " ORDER BY rowid"
        with self._connect() as connection:
            rows = connection.execute(query, parameters).fetchall()
        return [ArtifactRecord(**dict(row)) for row in rows]

    @staticmethod
    def _decode_run(row: sqlite3.Row) -> RunRecord:
        payload = dict(row)
        payload["argv"] = json.loads(payload.pop("argv_json"))
        return RunRecord(**payload)
