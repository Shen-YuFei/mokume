"""Deterministic control records for Studio scientific workflows."""

from __future__ import annotations

import hashlib
import json
import sqlite3
import uuid
from datetime import datetime, timedelta
from enum import Enum
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator

from mokume.studio.models import utc_now
from mokume.studio.state import StateStore


APPROVAL_LIFETIME = timedelta(minutes=30)
_MISSING = object()


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    )


def _payload_hash(payload: dict[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _expiry(created_at: str) -> str:
    created = datetime.fromisoformat(created_at)
    return (created + APPROVAL_LIFETIME).isoformat()


def _is_expired(expires_at: str, now: str) -> bool:
    return datetime.fromisoformat(expires_at) <= datetime.fromisoformat(now)


def _contains_api_key(value: Any) -> bool:
    if isinstance(value, dict):
        for key, nested in value.items():
            normalized = str(key).lower().replace("-", "_")
            if normalized in {"api_key", "apikey"} or _contains_api_key(nested):
                return True
    if isinstance(value, (list, tuple)):
        return any(_contains_api_key(item) for item in value)
    return False


class DatasetInspectionRequest(BaseModel):
    """Validated, serializable input for deterministic dataset inspection."""

    model_config = ConfigDict(extra="forbid")

    protein_matrix: str = Field(min_length=1)
    sdrf: str = Field(min_length=1)
    contrast: tuple[str, str]
    input_scale: Literal["linear", "log2"]
    peptide_counts: str | None = None
    data_type: Literal["LFQ", "DIA", "TMT"] | None = None
    quantification: str | None = None
    upstream_engine: str | None = None
    factor_column: str | None = None
    output_directory: str = "results/mokume"

    @field_validator(
        "protein_matrix",
        "sdrf",
        "peptide_counts",
        "quantification",
        "upstream_engine",
        "factor_column",
        "output_directory",
    )
    @classmethod
    def require_nonempty_strings(cls, value: str | None) -> str | None:
        """Reject blank path or metadata values without changing their spelling."""
        if value is not None and not value.strip():
            raise ValueError("value must not be blank")
        return value

    @field_validator("contrast")
    @classmethod
    def require_distinct_contrast(cls, value: tuple[str, str]) -> tuple[str, str]:
        """Require two non-empty, distinct contrast values."""
        normalized = tuple(item.strip() for item in value)
        if any(not item for item in normalized):
            raise ValueError("contrast values must not be blank")
        if normalized[0] == normalized[1]:
            raise ValueError("contrast values must be distinct")
        return normalized


class DatasetStatus(str, Enum):
    """Persisted lifecycle for one deterministic dataset inspection."""

    QUEUED = "queued"
    RUNNING = "running"
    READY = "ready"
    FAILED = "failed"


class DatasetRecord(BaseModel):
    """Control-plane reference to an inspection request and its result."""

    id: str
    project_id: str
    request: DatasetInspectionRequest
    status: DatasetStatus
    result: dict[str, Any] | None = None
    error: str | None = None
    created_at: str
    updated_at: str


class ApprovalStatus(str, Enum):
    """One-time authorization states for a scientific compute request."""

    PENDING = "pending"
    APPROVED = "approved"
    REJECTED = "rejected"
    CONSUMED = "consumed"
    EXPIRED = "expired"


class ApprovalRecord(BaseModel):
    """Hash-bound approval that can authorize at most one compute run."""

    id: str
    run_id: str | None = None
    kind: str
    payload: dict[str, Any]
    payload_hash: str
    status: ApprovalStatus
    created_at: str
    expires_at: str
    decided_at: str | None = None


class AgentThreadRecord(BaseModel):
    """Conversation state scoped to exactly one project and agent mode."""

    id: str
    project_id: str
    mode: str
    messages_json: str
    created_at: str
    updated_at: str

    @property
    def messages(self) -> list[dict[str, Any]]:
        """Return decoded messages while retaining canonical JSON in storage."""
        return json.loads(self.messages_json)


SCIENCE_SCHEMA = """
CREATE TABLE IF NOT EXISTS datasets (
    id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL REFERENCES projects(id),
    request_json TEXT NOT NULL,
    status TEXT NOT NULL,
    result_json TEXT,
    error TEXT,
    created_at TEXT NOT NULL,
    updated_at TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS datasets_project_updated
ON datasets(project_id, updated_at DESC);

CREATE TABLE IF NOT EXISTS agent_threads (
    id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL REFERENCES projects(id),
    mode TEXT NOT NULL,
    messages_json TEXT NOT NULL,
    created_at TEXT NOT NULL,
    updated_at TEXT NOT NULL
);
"""


class ScienceStore:
    """SQLite records layered on the existing Studio state database."""

    def __init__(self, state: StateStore) -> None:
        self.database = Path(state.database)
        with self._connect() as connection:
            connection.execute("PRAGMA journal_mode = WAL")
            connection.executescript(SCIENCE_SCHEMA)
            self._migrate_approvals(connection)

    def _connect(self) -> sqlite3.Connection:
        connection = sqlite3.connect(self.database, timeout=30)
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA foreign_keys = ON")
        return connection

    @staticmethod
    def _migrate_approvals(connection: sqlite3.Connection) -> None:
        columns = {
            row["name"]
            for row in connection.execute("PRAGMA table_info(approvals)").fetchall()
        }
        if not columns:
            raise RuntimeError("StateStore approvals table is missing")
        if "expires_at" not in columns:
            connection.execute("ALTER TABLE approvals ADD COLUMN expires_at TEXT")
        rows = connection.execute(
            "SELECT id, created_at FROM approvals WHERE expires_at IS NULL"
        ).fetchall()
        for row in rows:
            connection.execute(
                "UPDATE approvals SET expires_at=? WHERE id=?",
                (_expiry(row["created_at"]), row["id"]),
            )

    def create_dataset(
        self,
        project_id: str,
        request: DatasetInspectionRequest,
    ) -> DatasetRecord:
        """Persist one queued inspection request."""
        dataset_id = str(uuid.uuid4())
        timestamp = utc_now()
        with self._connect() as connection:
            connection.execute(
                """INSERT INTO datasets(
                    id, project_id, request_json, status, created_at, updated_at
                ) VALUES (?, ?, ?, ?, ?, ?)""",
                (
                    dataset_id,
                    project_id,
                    _canonical_json(request.model_dump(mode="json")),
                    DatasetStatus.QUEUED.value,
                    timestamp,
                    timestamp,
                ),
            )
        record = self.get_dataset(dataset_id)
        if record is None:  # pragma: no cover - protects against storage corruption
            raise RuntimeError(f"dataset was not persisted: {dataset_id}")
        return record

    def get_dataset(self, dataset_id: str) -> DatasetRecord | None:
        """Return one stored inspection record."""
        with self._connect() as connection:
            row = connection.execute(
                "SELECT * FROM datasets WHERE id=?", (dataset_id,)
            ).fetchone()
        return self._decode_dataset(row) if row else None

    def latest_dataset(self, project_id: str) -> DatasetRecord | None:
        """Return the most recently updated dataset for one project."""
        with self._connect() as connection:
            row = connection.execute(
                """SELECT * FROM datasets WHERE project_id=?
                ORDER BY updated_at DESC, rowid DESC LIMIT 1""",
                (project_id,),
            ).fetchone()
        return self._decode_dataset(row) if row else None

    def interrupt_incomplete_datasets(self) -> int:
        """Mark inspections left active by a server restart as failed."""
        with self._connect() as connection:
            cursor = connection.execute(
                """UPDATE datasets SET status=?, error=?, updated_at=?
                WHERE status IN (?, ?)""",
                (
                    DatasetStatus.FAILED.value,
                    "Studio restarted during dataset inspection",
                    utc_now(),
                    DatasetStatus.QUEUED.value,
                    DatasetStatus.RUNNING.value,
                ),
            )
        return cursor.rowcount

    def update_dataset(
        self,
        dataset_id: str,
        status: DatasetStatus,
        *,
        result: dict[str, Any] | None | object = _MISSING,
        error: str | None | object = _MISSING,
    ) -> DatasetRecord:
        """Update inspection state while preserving omitted result fields."""
        update_result = result is not _MISSING
        result_json = (
            None if result is None or result is _MISSING else _canonical_json(result)
        )
        update_error = error is not _MISSING
        with self._connect() as connection:
            cursor = connection.execute(
                """UPDATE datasets SET status=?, updated_at=?,
                    result_json=CASE WHEN ? THEN ? ELSE result_json END,
                    error=CASE WHEN ? THEN ? ELSE error END
                WHERE id=?""",
                (
                    status.value,
                    utc_now(),
                    update_result,
                    result_json,
                    update_error,
                    None if error is _MISSING else error,
                    dataset_id,
                ),
            )
        if cursor.rowcount != 1:
            raise KeyError(f"unknown dataset: {dataset_id}")
        record = self.get_dataset(dataset_id)
        if record is None:  # pragma: no cover - protects against storage corruption
            raise RuntimeError(f"dataset disappeared: {dataset_id}")
        return record

    def create_approval(
        self,
        kind: str,
        payload: dict[str, Any],
    ) -> ApprovalRecord:
        """Create a pending approval with a fixed thirty-minute lifetime."""
        if not kind.strip():
            raise ValueError("approval kind must not be blank")
        if _contains_api_key(payload):
            raise ValueError("API keys must not be stored in Studio state")
        approval_id = str(uuid.uuid4())
        created_at = utc_now()
        with self._connect() as connection:
            connection.execute(
                """INSERT INTO approvals(
                    id, run_id, kind, payload_json, payload_hash, status,
                    created_at, decided_at, expires_at
                ) VALUES (?, NULL, ?, ?, ?, ?, ?, NULL, ?)""",
                (
                    approval_id,
                    kind,
                    _canonical_json(payload),
                    _payload_hash(payload),
                    ApprovalStatus.PENDING.value,
                    created_at,
                    _expiry(created_at),
                ),
            )
        record = self.get_approval(approval_id)
        if record is None:  # pragma: no cover - protects against storage corruption
            raise RuntimeError(f"approval was not persisted: {approval_id}")
        return record

    def get_approval(self, approval_id: str) -> ApprovalRecord | None:
        """Return one approval without changing its lifecycle state."""
        with self._connect() as connection:
            row = connection.execute(
                "SELECT * FROM approvals WHERE id=?", (approval_id,)
            ).fetchone()
        return self._decode_approval(row) if row else None

    def decide_approval(
        self,
        approval_id: str,
        *,
        payload_hash: str,
        approved: bool,
    ) -> ApprovalRecord:
        """Atomically approve or reject an unchanged pending payload."""
        target = ApprovalStatus.APPROVED if approved else ApprovalStatus.REJECTED
        return self._transition_approval(
            approval_id,
            payload_hash,
            expected=ApprovalStatus.PENDING,
            target=target,
        )

    def consume_approval(
        self,
        approval_id: str,
        *,
        payload_hash: str,
    ) -> ApprovalRecord:
        """Atomically consume an approved payload exactly once."""
        return self._transition_approval(
            approval_id,
            payload_hash,
            expected=ApprovalStatus.APPROVED,
            target=ApprovalStatus.CONSUMED,
        )

    def _transition_approval(
        self,
        approval_id: str,
        payload_hash: str,
        *,
        expected: ApprovalStatus,
        target: ApprovalStatus,
    ) -> ApprovalRecord:
        now = utc_now()
        connection = self._connect()
        try:
            connection.execute("BEGIN IMMEDIATE")
            row = connection.execute(
                "SELECT * FROM approvals WHERE id=?", (approval_id,)
            ).fetchone()
            error = self._approval_transition_error(
                connection, row, payload_hash, expected, now
            )
            if error is None:
                cursor = connection.execute(
                    """UPDATE approvals SET status=?,
                        decided_at=COALESCE(decided_at, ?)
                    WHERE id=? AND status=? AND payload_hash=?""",
                    (
                        target.value,
                        now,
                        approval_id,
                        expected.value,
                        payload_hash,
                    ),
                )
                if cursor.rowcount != 1:  # pragma: no cover - write lock prevents it
                    connection.rollback()
                    raise ValueError("approval changed during transition")
                connection.commit()
            else:
                if error == "approval expired":
                    connection.commit()
                else:
                    connection.rollback()
                raise ValueError(error)
        finally:
            connection.close()
        record = self.get_approval(approval_id)
        if record is None:  # pragma: no cover - protects against storage corruption
            raise RuntimeError(f"approval disappeared: {approval_id}")
        return record

    @staticmethod
    def _approval_transition_error(
        connection: sqlite3.Connection,
        row: sqlite3.Row | None,
        payload_hash: str,
        expected: ApprovalStatus,
        now: str,
    ) -> str | None:
        if row is None:
            return "unknown approval"
        if row["payload_hash"] != payload_hash:
            return "approval payload hash mismatch"
        if row["status"] != expected.value:
            return f"approval is already {row['status']}"
        if _is_expired(row["expires_at"], now):
            connection.execute(
                """UPDATE approvals SET status=?,
                    decided_at=COALESCE(decided_at, ?) WHERE id=?""",
                (ApprovalStatus.EXPIRED.value, now, row["id"]),
            )
            return "approval expired"
        return None

    def link_run(
        self,
        approval_id: str,
        *,
        payload_hash: str,
        run_id: str,
    ) -> ApprovalRecord:
        """Attach a consumed approval to exactly one worker run."""
        with self._connect() as connection:
            cursor = connection.execute(
                """UPDATE approvals SET run_id=?
                WHERE id=? AND payload_hash=? AND status=? AND run_id IS NULL""",
                (
                    run_id,
                    approval_id,
                    payload_hash,
                    ApprovalStatus.CONSUMED.value,
                ),
            )
        if cursor.rowcount != 1:
            raise ValueError("approval cannot be linked to this run")
        record = self.get_approval(approval_id)
        if record is None:  # pragma: no cover - protects against storage corruption
            raise RuntimeError(f"approval disappeared: {approval_id}")
        return record

    def get_thread(
        self,
        thread_id: str,
        *,
        project_id: str,
        mode: str,
    ) -> AgentThreadRecord | None:
        """Read a thread only through its original project and mode scope."""
        with self._connect() as connection:
            row = connection.execute(
                "SELECT * FROM agent_threads WHERE id=?", (thread_id,)
            ).fetchone()
        if row and (row["project_id"] != project_id or row["mode"] != mode):
            raise ValueError("agent thread belongs to a different project or mode")
        return AgentThreadRecord(**dict(row)) if row else None

    def save_thread(
        self,
        thread_id: str,
        *,
        project_id: str,
        mode: str,
        messages: list[dict[str, Any]],
    ) -> AgentThreadRecord:
        """Create or update a canonically encoded, scope-bound conversation."""
        if not mode.strip():
            raise ValueError("agent mode must not be blank")
        if _contains_api_key(messages):
            raise ValueError("API keys must not be stored in Studio state")
        connection = self._connect()
        try:
            connection.execute("BEGIN IMMEDIATE")
            row = connection.execute(
                "SELECT project_id, mode FROM agent_threads WHERE id=?", (thread_id,)
            ).fetchone()
            if row and (row["project_id"] != project_id or row["mode"] != mode):
                connection.rollback()
                raise ValueError("agent thread belongs to a different project or mode")
            self._upsert_thread(connection, thread_id, project_id, mode, messages)
            connection.commit()
        finally:
            connection.close()
        record = self.get_thread(thread_id, project_id=project_id, mode=mode)
        if record is None:  # pragma: no cover - protects against storage corruption
            raise RuntimeError(f"agent thread was not persisted: {thread_id}")
        return record

    @staticmethod
    def _upsert_thread(
        connection: sqlite3.Connection,
        thread_id: str,
        project_id: str,
        mode: str,
        messages: list[dict[str, Any]],
    ) -> None:
        timestamp = utc_now()
        connection.execute(
            """INSERT INTO agent_threads(
                id, project_id, mode, messages_json, created_at, updated_at
            ) VALUES (?, ?, ?, ?, ?, ?)
            ON CONFLICT(id) DO UPDATE SET
                messages_json=excluded.messages_json,
                updated_at=excluded.updated_at""",
            (
                thread_id,
                project_id,
                mode,
                _canonical_json(messages),
                timestamp,
                timestamp,
            ),
        )

    @staticmethod
    def _decode_dataset(row: sqlite3.Row) -> DatasetRecord:
        payload = dict(row)
        payload["request"] = json.loads(payload.pop("request_json"))
        result_json = payload.pop("result_json")
        payload["result"] = json.loads(result_json) if result_json else None
        return DatasetRecord(**payload)

    @staticmethod
    def _decode_approval(row: sqlite3.Row) -> ApprovalRecord:
        payload = dict(row)
        payload["payload"] = json.loads(payload.pop("payload_json"))
        return ApprovalRecord(**payload)
