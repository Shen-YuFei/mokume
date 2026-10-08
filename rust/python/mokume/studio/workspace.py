"""Workspace-scoped file browsing and workflow template persistence."""

from __future__ import annotations

import json
import uuid
from collections.abc import Callable
from pathlib import Path
from typing import Any

from fastapi import Depends, FastAPI, HTTPException, Query, status

from mokume.studio.models import (
    WorkflowTemplateDocument,
    WorkflowTemplateWriteRequest,
)
from mokume.studio.paths import PathAccessError, ProjectPaths, directory_children
from mokume.studio.state import StateStore


Dependency = Callable[..., Any]


def active_project_paths(store: StateStore) -> ProjectPaths:
    """Return the active workspace guard or an HTTP conflict."""
    project = store.active_project()
    if project is None:
        raise HTTPException(status.HTTP_409_CONFLICT, "Open a project folder first")
    try:
        return ProjectPaths(project.root)
    except PathAccessError as exc:
        raise HTTPException(status.HTTP_409_CONFLICT, str(exc)) from exc


def list_project_entries(guard: ProjectPaths, directory: Path) -> list[dict]:
    """Return folders before files while omitting unreadable or escaping links."""
    entries = []
    for child in directory_children(directory):
        try:
            resolved = guard.resolve_existing(child)
        except PathAccessError:
            continue
        entries.append(
            {
                "name": child.name,
                "path": str(resolved.relative_to(guard.root)),
                "kind": "directory" if resolved.is_dir() else "file",
                "size": resolved.stat().st_size if resolved.is_file() else None,
            }
        )
    entries.sort(
        key=lambda entry: (entry["kind"] != "directory", entry["name"].casefold())
    )
    return entries


def install_workflow_template_routes(
    app: FastAPI,
    store: StateStore,
    require_session: Dependency,
    require_mutation: Dependency,
) -> None:
    """Install workspace-only JSON workflow template endpoints."""

    @app.get("/api/workflow-template")
    async def read_workflow_template(
        path: str = Query(min_length=1),
        _session=Depends(require_session),
    ) -> dict:
        guard = active_project_paths(store)
        try:
            source = guard.resolve_existing(path)
            if not source.is_file() or source.suffix.casefold() != ".json":
                raise PathAccessError("workflow template must be a JSON file")
            document = WorkflowTemplateDocument.model_validate_json(
                source.read_text(encoding="utf-8")
            )
        except (OSError, ValueError) as exc:
            raise HTTPException(
                status.HTTP_422_UNPROCESSABLE_CONTENT, str(exc)
            ) from exc
        return {
            "path": guard.relative(source),
            "template": document.model_dump(mode="json", by_alias=True),
        }

    @app.put("/api/workflow-template")
    async def write_workflow_template(
        payload: WorkflowTemplateWriteRequest,
        _session=Depends(require_mutation),
    ) -> dict:
        guard = active_project_paths(store)
        try:
            if Path(payload.path).suffix.casefold() != ".json":
                raise PathAccessError("workflow template must be a JSON file")
            destination = guard.resolve_output(
                payload.path, allow_existing=payload.overwrite
            )
            if destination.exists() and not destination.is_file():
                raise PathAccessError("workflow template path is not a file")
            _write_workspace_json(
                destination,
                payload.template.model_dump(mode="json", by_alias=True),
            )
        except (OSError, ValueError) as exc:
            raise HTTPException(
                status.HTTP_422_UNPROCESSABLE_CONTENT, str(exc)
            ) from exc
        return {"path": guard.relative(destination)}


def _write_workspace_json(path: Path, payload: dict) -> None:
    """Atomically write JSON without leaving temporary files behind."""
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        temporary.write_text(
            json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)
