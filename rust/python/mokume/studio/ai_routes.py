"""Authenticated HTTP and AG-UI routes for optional Studio AI."""

from __future__ import annotations

import json
from collections.abc import Callable
from typing import Any, Protocol, cast

from fastapi import Depends, FastAPI, HTTPException, Query, Request, status
from pydantic_ai.messages import ModelMessagesTypeAdapter
from pydantic_ai.ui.ag_ui import AGUIAdapter

from mokume.studio.agent import AgentMode, StudioAgentDeps, agent_for
from mokume.studio.auth import Session
from mokume.studio.models import ApprovalDecision
from mokume.studio.providers import ProviderConfig, ProviderRegistry
from mokume.studio.science import DatasetInspectionRequest, ScienceStore
from mokume.studio.scientific import ScientificController
from mokume.studio.state import StateStore


Dependency = Callable[..., Any]


class AIRuntime(Protocol):
    """Runtime services required without importing the app module."""

    @property
    def providers(self) -> ProviderRegistry:
        """Return process-local provider credentials."""

    @property
    def science(self) -> ScienceStore:
        """Return persisted scientific control records."""

    @property
    def scientific(self) -> ScientificController:
        """Return the deterministic scientific controller."""

    @property
    def store(self) -> StateStore:
        """Return the Studio run-state store."""


def install_ai_routes(
    app: FastAPI,
    runtime: AIRuntime,
    require_session: Dependency,
    require_mutation: Dependency,
) -> None:
    """Install provider, inspection, approval, and AG-UI endpoints."""
    _install_provider_routes(app, runtime, require_session, require_mutation)
    _install_dataset_routes(app, runtime, require_session, require_mutation)
    _install_approval_routes(app, runtime, require_session, require_mutation)
    _install_agent_routes(app, runtime, require_session, require_mutation)


def _install_provider_routes(
    app: FastAPI,
    runtime: AIRuntime,
    require_session: Dependency,
    require_mutation: Dependency,
) -> None:
    @app.get("/api/ai/config")
    async def provider_config(session: Session = Depends(require_session)):
        return runtime.providers.summary(session.id)

    @app.post("/api/ai/config")
    async def configure_provider(
        payload: ProviderConfig,
        session: Session = Depends(require_mutation),
    ):
        return runtime.providers.save(session.id, payload)

    @app.delete("/api/ai/config", status_code=status.HTTP_204_NO_CONTENT)
    async def clear_provider(
        session: Session = Depends(require_mutation),
    ) -> None:
        runtime.providers.clear(session.id)


def _install_dataset_routes(
    app: FastAPI,
    runtime: AIRuntime,
    require_session: Dependency,
    require_mutation: Dependency,
) -> None:
    @app.post("/api/datasets/inspect", status_code=status.HTTP_202_ACCEPTED)
    async def inspect_dataset(
        payload: DatasetInspectionRequest,
        _session: Session = Depends(require_mutation),
    ) -> dict:
        project = _active_project(runtime)
        try:
            dataset, run = runtime.scientific.inspect(payload, project)
        except (OSError, RuntimeError, TypeError, ValueError) as exc:
            raise HTTPException(
                status.HTTP_422_UNPROCESSABLE_CONTENT, str(exc)
            ) from exc
        return {"dataset": dataset, "run": run}

    @app.get("/api/datasets/latest")
    async def latest_dataset(_session: Session = Depends(require_session)):
        project = _active_project(runtime)
        return runtime.science.latest_dataset(project.id)

    @app.get("/api/datasets/{dataset_id}")
    async def dataset(
        dataset_id: str,
        _session: Session = Depends(require_session),
    ):
        project = _active_project(runtime)
        record = runtime.science.get_dataset(dataset_id)
        if record is None or record.project_id != project.id:
            raise HTTPException(status.HTTP_404_NOT_FOUND, "Dataset not found")
        return record


def _install_approval_routes(
    app: FastAPI,
    runtime: AIRuntime,
    require_session: Dependency,
    require_mutation: Dependency,
) -> None:
    @app.get("/api/approvals/{approval_id}")
    async def approval(
        approval_id: str,
        _session: Session = Depends(require_session),
    ):
        project = _active_project(runtime)
        record = runtime.science.get_approval(approval_id)
        if record is None or record.payload.get("project_id") != project.id:
            raise HTTPException(status.HTTP_404_NOT_FOUND, "Approval not found")
        return record

    @app.post("/api/approvals/{approval_id}")
    async def decide_approval(
        approval_id: str,
        decision: ApprovalDecision,
        _session: Session = Depends(require_mutation),
    ):
        project = _active_project(runtime)
        record = runtime.science.get_approval(approval_id)
        if record is None or record.payload.get("project_id") != project.id:
            raise HTTPException(status.HTTP_404_NOT_FOUND, "Approval not found")
        try:
            return runtime.science.decide_approval(
                approval_id,
                payload_hash=decision.payload_hash,
                approved=decision.approved,
            )
        except ValueError as exc:
            raise HTTPException(status.HTTP_409_CONFLICT, str(exc)) from exc


def _install_agent_routes(
    app: FastAPI,
    runtime: AIRuntime,
    require_session: Dependency,
    require_mutation: Dependency,
) -> None:
    @app.post("/api/agent/run")
    async def run_agent(
        request: Request,
        session: Session = Depends(require_mutation),
    ):
        body = await _validated_agent_body(request)
        project = _active_project(runtime)
        mode, dataset_id = _agent_options(body)
        summary = runtime.providers.summary(session.id)
        if summary is None:
            raise HTTPException(
                status.HTTP_409_CONFLICT, "Configure an AI provider first"
            )
        try:
            model = runtime.providers.model_for(session.id)
            thread_id = cast(str, body["threadId"])
            thread = runtime.science.get_thread(
                thread_id,
                project_id=project.id,
                mode=mode,
            )
        except (LookupError, RuntimeError, ValueError) as exc:
            raise HTTPException(status.HTTP_409_CONFLICT, str(exc)) from exc
        history = (
            ModelMessagesTypeAdapter.validate_json(thread.messages_json)
            if thread
            else None
        )
        deps = StudioAgentDeps(
            controller=runtime.scientific,
            project=project,
            provider=summary.provider,
            model=summary.model,
            dataset_id=dataset_id,
        )

        async def save_messages(result) -> None:
            messages = json.loads(result.all_messages_json())
            runtime.science.save_thread(
                thread_id,
                project_id=project.id,
                mode=mode,
                messages=messages,
            )

        return await AGUIAdapter.dispatch_request(
            request,
            agent=agent_for(mode),
            model=model,
            deps=deps,
            message_history=history,
            conversation_id=thread_id,
            on_complete=save_messages,
            manage_system_prompt="server",
            allowed_file_url_schemes=frozenset(),
            allow_uploaded_files=False,
        )

    @app.get("/api/agent/threads/{thread_id}")
    async def agent_thread(
        thread_id: str,
        mode: AgentMode = Query(),
        _session: Session = Depends(require_session),
    ):
        project = _active_project(runtime)
        try:
            return runtime.science.get_thread(
                thread_id,
                project_id=project.id,
                mode=mode,
            )
        except ValueError as exc:
            raise HTTPException(status.HTTP_404_NOT_FOUND, str(exc)) from exc


async def _validated_agent_body(request: Request) -> dict[str, Any]:
    try:
        body = await request.json()
    except (UnicodeDecodeError, ValueError) as exc:
        raise HTTPException(status.HTTP_400_BAD_REQUEST, "Invalid AG-UI JSON") from exc
    if not isinstance(body, dict):
        raise HTTPException(
            status.HTTP_422_UNPROCESSABLE_CONTENT, "AG-UI body required"
        )
    for name in ("threadId", "runId"):
        value = body.get(name)
        if not isinstance(value, str) or not value.strip() or len(value) > 200:
            raise HTTPException(
                status.HTTP_422_UNPROCESSABLE_CONTENT, f"{name} must be a short string"
            )
    if body.get("tools") or body.get("context") or body.get("state"):
        raise HTTPException(
            status.HTTP_422_UNPROCESSABLE_CONTENT,
            "Client tools, context, and state are not accepted",
        )
    _validate_agent_messages(body.get("messages"), bool(body.get("resume")))
    return body


def _validate_agent_messages(messages: Any, has_resume: bool) -> None:
    if not isinstance(messages, list) or len(messages) > 1:
        raise HTTPException(
            status.HTTP_422_UNPROCESSABLE_CONTENT,
            "Send exactly the new user message, not client-side history",
        )
    if has_resume and messages:
        raise HTTPException(
            status.HTTP_422_UNPROCESSABLE_CONTENT,
            "Approval resume cannot include a new message",
        )
    if not has_resume and len(messages) != 1:
        raise HTTPException(
            status.HTTP_422_UNPROCESSABLE_CONTENT, "One user message is required"
        )
    if messages:
        message = messages[0]
        if not isinstance(message, dict):
            raise HTTPException(
                status.HTTP_422_UNPROCESSABLE_CONTENT,
                "Only a plain-text user message is accepted",
            )
        content = message.get("content")
        if message.get("role") != "user" or not isinstance(content, str):
            raise HTTPException(
                status.HTTP_422_UNPROCESSABLE_CONTENT,
                "Only a plain-text user message is accepted",
            )
        if not content.strip() or len(content) > 8000:
            raise HTTPException(
                status.HTTP_422_UNPROCESSABLE_CONTENT,
                "User message must contain 1 to 8000 characters",
            )


def _agent_options(body: dict[str, Any]) -> tuple[AgentMode, str | None]:
    forwarded = body.get("forwardedProps") or {}
    if not isinstance(forwarded, dict) or set(forwarded) - {"mode", "datasetId"}:
        raise HTTPException(
            status.HTTP_422_UNPROCESSABLE_CONTENT, "Invalid agent options"
        )
    mode = forwarded.get("mode", "ask")
    if mode not in {"ask", "plan"}:
        raise HTTPException(status.HTTP_422_UNPROCESSABLE_CONTENT, "Invalid agent mode")
    dataset_id = forwarded.get("datasetId")
    if dataset_id is not None and (
        not isinstance(dataset_id, str) or not dataset_id.strip()
    ):
        raise HTTPException(status.HTTP_422_UNPROCESSABLE_CONTENT, "Invalid dataset ID")
    return cast(AgentMode, mode), dataset_id


def _active_project(runtime: AIRuntime):
    project = runtime.store.active_project()
    if project is None:
        raise HTTPException(status.HTTP_409_CONFLICT, "Open a project folder first")
    return project
