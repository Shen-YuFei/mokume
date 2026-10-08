"""Typed control-plane records used by Mokume Studio."""

from __future__ import annotations

from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any

from pydantic import BaseModel, ConfigDict, Field, field_validator


def utc_now() -> str:
    """Return an ISO-8601 timestamp in UTC."""
    return datetime.now(timezone.utc).isoformat()


class RunStatus(str, Enum):
    """Persisted lifecycle states for one compute worker."""

    QUEUED = "queued"
    STARTING = "starting"
    RUNNING = "running"
    CANCELLING = "cancelling"
    CANCELLED = "cancelled"
    SUCCEEDED = "succeeded"
    FAILED = "failed"
    INTERRUPTED = "interrupted"


class JobOperation(str, Enum):
    """Worker operations accepted from trusted server controllers."""

    NATIVE = "native"
    INSPECT_DATASET = "inspect_dataset"
    EVALUATE_RECOMMENDATION = "evaluate_recommendation"


TERMINAL_RUN_STATUSES = {
    RunStatus.CANCELLED,
    RunStatus.SUCCEEDED,
    RunStatus.FAILED,
    RunStatus.INTERRUPTED,
}


class FolderEntry(BaseModel):
    """One directory returned by the authenticated folder browser."""

    name: str
    path: str


class ProjectRecord(BaseModel):
    """The currently selected local project root."""

    id: str
    root: str
    opened_at: str


class OpenProjectRequest(BaseModel):
    """Request to select a readable local directory as the project root."""

    path: str = Field(min_length=1)


class RunRequest(BaseModel):
    """Canonical user-approved native command request."""

    argv: list[str] = Field(min_length=1)
    output_directory: str = "results/mokume"

    @field_validator("argv")
    @classmethod
    def reject_empty_arguments(cls, value: list[str]) -> list[str]:
        """Reject arguments that cannot safely cross process boundaries."""
        if any(not item or "\x00" in item for item in value):
            raise ValueError("argv entries must be non-empty and contain no NUL bytes")
        return value


class JobSpec(BaseModel):
    """Immutable specification reloaded and verified by a worker process."""

    model_config = ConfigDict(frozen=True)

    run_id: str
    project_root: str
    run_directory: str
    operation: JobOperation = JobOperation.NATIVE
    argv: list[str]
    payload: dict[str, Any] = Field(default_factory=dict)
    parameters: dict[str, Any]
    approved_hash: str
    created_at: str
    threads: int = 24

    @property
    def path(self) -> Path:
        """Return the immutable run directory as a path."""
        return Path(self.run_directory)


class ScientificJobRequest(BaseModel):
    """Trusted controller request for one isolated scientific worker."""

    operation: JobOperation
    payload: dict[str, Any]
    input_paths: list[str]
    output_directory: str = "results/mokume"


class RunRecord(BaseModel):
    """Persisted summary of one Studio run."""

    id: str
    project_id: str
    status: RunStatus
    command: str
    argv: list[str]
    approved_hash: str
    run_directory: str
    worker_pid: int | None = None
    created_at: str
    started_at: str | None = None
    finished_at: str | None = None
    error: str | None = None


class ArtifactRecord(BaseModel):
    """A worker output addressable through an opaque artifact ID."""

    id: str
    run_id: str
    path: str
    media_type: str
    size: int
    sha256: str


class ValidationRequest(BaseModel):
    """A parse-only command validation request."""

    argv: list[str] = Field(min_length=1)


class WorkflowTemplateDocument(BaseModel):
    """A portable Studio workflow configuration."""

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    schema_version: int = Field(alias="$schemaVersion")
    workflow: list[str] = Field(min_length=1)
    parameters: dict[str, Any]

    @field_validator("workflow")
    @classmethod
    def reject_invalid_workflow_parts(cls, value: list[str]) -> list[str]:
        """Keep workflow command parts safe and unambiguous."""
        if any(not part.strip() or "\x00" in part for part in value):
            raise ValueError(
                "workflow parts must be non-empty and contain no NUL bytes"
            )
        return value


class WorkflowTemplateWriteRequest(BaseModel):
    """Write one workflow template beneath the active project root."""

    model_config = ConfigDict(extra="forbid")

    path: str = Field(min_length=1)
    template: WorkflowTemplateDocument
    overwrite: bool = False

    @field_validator("path")
    @classmethod
    def reject_invalid_template_path(cls, value: str) -> str:
        """Reject paths that cannot be passed safely to the filesystem."""
        if "\x00" in value:
            raise ValueError("template path must not contain NUL bytes")
        return value


class ApprovalDecision(BaseModel):
    """Server-validated response to a pending AI or compute approval."""

    approved: bool
    payload_hash: str = Field(min_length=64, max_length=64, pattern=r"^[0-9a-f]{64}$")
