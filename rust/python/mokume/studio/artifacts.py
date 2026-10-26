"""Project-scoped artifact downloads and sandboxed report previews."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

from fastapi import Depends, FastAPI, HTTPException, Query, status
from fastapi.responses import FileResponse

from mokume.studio.auth import Session
from mokume.studio.models import ArtifactRecord
from mokume.studio.paths import PathAccessError, ProjectPaths


Dependency = Callable[..., Any]
ARTIFACT_PREVIEW_POLICY = (
    "sandbox allow-scripts; default-src 'none'; script-src 'unsafe-inline'; "
    "style-src 'unsafe-inline'; img-src data: blob:; font-src data:; "
    "frame-ancestors 'self'"
)


def _artifact_file(runtime: Any, artifact_id: str) -> tuple[ArtifactRecord, Path]:
    """Resolve an artifact only within the authenticated active project."""
    record = runtime.store.get_artifact(artifact_id)
    if record is None:
        raise HTTPException(status.HTTP_404_NOT_FOUND, "Artifact not found")
    run_record = runtime.store.get_run(record.run_id)
    project = runtime.store.active_project()
    if project is None or run_record is None or run_record.project_id != project.id:
        raise HTTPException(status.HTTP_404_NOT_FOUND, "Artifact project not found")
    try:
        return record, ProjectPaths(project.root).resolve_existing(record.path)
    except PathAccessError as exc:
        raise HTTPException(status.HTTP_404_NOT_FOUND, str(exc)) from exc


def install_artifact_routes(
    app: FastAPI, runtime: Any, require_session: Dependency
) -> None:
    """Install project-scoped artifact downloads and isolated HTML previews."""

    @app.get("/api/artifacts")
    async def artifacts(
        _session: Session = Depends(require_session),
        run_id: str | None = Query(default=None),
    ) -> dict:
        project = runtime.store.active_project()
        if project is None:
            return {"artifacts": []}
        if run_id is not None:
            record = runtime.store.get_run(run_id)
            if record is None or record.project_id != project.id:
                raise HTTPException(status.HTTP_404_NOT_FOUND, "Run not found")
            return {"artifacts": runtime.store.list_artifacts(run_id)}
        return {"artifacts": runtime.store.list_project_artifacts(project.id)}

    @app.get("/api/artifacts/{artifact_id}")
    async def artifact(artifact_id: str, _session: Session = Depends(require_session)):
        record, path = _artifact_file(runtime, artifact_id)
        headers = (
            {"Content-Security-Policy": "sandbox"}
            if record.media_type == "image/svg+xml"
            else {}
        )
        return FileResponse(
            path,
            media_type=record.media_type,
            filename=path.name,
            content_disposition_type=(
                "attachment" if record.media_type == "text/html" else "inline"
            ),
            headers=headers,
        )

    @app.get("/api/artifacts/{artifact_id}/preview")
    async def artifact_preview(
        artifact_id: str, _session: Session = Depends(require_session)
    ):
        record, path = _artifact_file(runtime, artifact_id)
        if record.media_type != "text/html":
            raise HTTPException(
                status.HTTP_415_UNSUPPORTED_MEDIA_TYPE,
                "Only HTML artifacts can be previewed",
            )
        return FileResponse(
            path,
            media_type="text/html",
            filename=path.name,
            content_disposition_type="inline",
            headers={
                "Content-Security-Policy": ARTIFACT_PREVIEW_POLICY,
                "X-Frame-Options": "SAMEORIGIN",
            },
        )
