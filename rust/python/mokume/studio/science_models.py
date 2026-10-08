"""Validated records persisted by the Studio scientific control store."""

from __future__ import annotations

import json
from enum import Enum
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator


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


def _append_conversation_text(
    conversation: list[dict[str, str]], role: str, content: str
) -> None:
    text = content.strip()
    if not text:
        return
    if conversation and conversation[-1]["role"] == role:
        conversation[-1]["text"] = f"{conversation[-1]['text']}\n{text}"
    else:
        conversation.append({"role": role, "text": text})


class AgentThreadRecord(BaseModel):
    """Conversation state scoped to exactly one project and agent mode."""

    id: str
    project_id: str
    mode: str
    messages_json: str
    custom_title: str | None = None
    created_at: str
    updated_at: str

    @property
    def messages(self) -> list[dict[str, Any]]:
        """Return decoded messages while retaining canonical JSON in storage."""
        return json.loads(self.messages_json)

    @property
    def conversation(self) -> list[dict[str, str]]:
        """Return visible user, assistant, and model reasoning text for history."""
        conversation = []
        for message in self.messages:
            role = message.get("role")
            content = message.get("content")
            if role in {"user", "assistant"} and isinstance(content, str):
                _append_conversation_text(conversation, role, content)
            for part in message.get("parts", []):
                if not isinstance(part, dict):
                    continue
                part_kind = part.get("part_kind")
                part_role = (
                    "user"
                    if part_kind == "user-prompt"
                    else "assistant"
                    if part_kind == "text"
                    else "reasoning"
                    if part_kind == "thinking"
                    else None
                )
                if part_role and isinstance(part.get("content"), str):
                    _append_conversation_text(conversation, part_role, part["content"])
        return conversation

    @property
    def title(self) -> str:
        """Use the first user message as a compact history title."""
        if self.custom_title:
            return self.custom_title
        first = next(
            (item["text"] for item in self.conversation if item["role"] == "user"),
            "New conversation",
        )
        compact = " ".join(first.split())
        return f"{compact[:77]}…" if len(compact) > 78 else compact
