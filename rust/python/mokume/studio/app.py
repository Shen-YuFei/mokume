"""FastAPI application factory for the loopback-only Mokume Studio UI."""

from __future__ import annotations

import asyncio
import importlib
import json
from collections.abc import Callable
from contextlib import asynccontextmanager
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from typing import Any

import psutil
from fastapi import Depends, FastAPI, Header, HTTPException, Query, Request, status
from fastapi.responses import (
    FileResponse,
    HTMLResponse,
    RedirectResponse,
    StreamingResponse,
)
from fastapi.templating import Jinja2Templates

from mokume.studio.ai_routes import install_ai_routes
from mokume.studio.auth import (
    CSRF_HEADER,
    SESSION_COOKIE,
    Session,
    SessionManager,
    csrf_is_allowed,
    loopback_origin_from_host,
    origin_is_allowed,
)
from mokume.studio.catalog import (
    CommandValidationError,
    command_schema,
    validate_and_canonicalize,
)
from mokume.studio.jobs import JobConflictError, JobManager
from mokume.studio.models import (
    TERMINAL_RUN_STATUSES,
    OpenProjectRequest,
    RunRequest,
    ValidationRequest,
)
from mokume.studio.paths import (
    PathAccessError,
    ProjectPaths,
    list_directories,
    readable_directory,
)
from mokume.studio.state import StateStore
from mokume.studio.providers import ProviderRegistry, default_provider_config_root
from mokume.studio.science import ScienceStore
from mokume.studio.scientific import ScientificController
from mokume.studio.workspace import (
    active_project_paths,
    install_workflow_template_routes,
    list_project_entries,
)


CONTENT_SECURITY_POLICY = "; ".join(
    (
        "default-src 'none'",
        "script-src 'self'",
        "style-src 'self'",
        "img-src 'self' data:",
        "connect-src 'self'",
        "font-src 'self'",
        "frame-src 'self'",
        "base-uri 'none'",
        "form-action 'self'",
        "frame-ancestors 'none'",
    )
)

Dependency = Callable[..., Any]


@dataclass(frozen=True)
class StudioAI:
    """AI and scientific services that are optional at the UI boundary."""

    providers: ProviderRegistry
    science: ScienceStore
    scientific: ScientificController


@dataclass(frozen=True)
class StudioRuntime:
    """Process-local services shared by the HTTP route groups."""

    sessions: SessionManager
    store: StateStore
    jobs: JobManager
    ai: StudioAI
    shutdown_callback: Callable[[], None] | None

    @property
    def providers(self) -> ProviderRegistry:
        """Expose session and persistent provider state to the AI route contract."""
        return self.ai.providers

    @property
    def science(self) -> ScienceStore:
        """Expose scientific records to the AI route contract."""
        return self.ai.science

    @property
    def scientific(self) -> ScientificController:
        """Expose the scientific controller to the AI route contract."""
        return self.ai.scientific


def create_app(
    *,
    startup_token: str,
    state_directory: str | Path | None = None,
    shutdown_callback: Callable[[], None] | None = None,
) -> FastAPI:
    """Build one single-user Studio app bound by the CLI to a loopback port."""
    store = StateStore(state_directory)
    store.interrupt_incomplete_runs()
    active_project = store.active_project()
    if active_project is not None:
        try:
            readable_directory(active_project.root)
        except PathAccessError:
            store.close_project()
    jobs = JobManager(store)
    science = ScienceStore(store)
    science.interrupt_incomplete_datasets()
    runtime = StudioRuntime(
        sessions=SessionManager(startup_token),
        store=store,
        jobs=jobs,
        ai=StudioAI(
            providers=ProviderRegistry(default_provider_config_root()),
            science=science,
            scientific=ScientificController(science, jobs),
        ),
        shutdown_callback=shutdown_callback,
    )

    @asynccontextmanager
    async def lifespan(_app: FastAPI):
        yield
        runtime.jobs.shutdown()

    app = FastAPI(
        title="Mokume Studio",
        docs_url=None,
        redoc_url=None,
        openapi_url=None,
        lifespan=lifespan,
    )
    app.state.runtime = runtime
    require_session, require_mutation = _security_dependencies(runtime)
    _install_security_middleware(app)
    _install_page_routes(app, runtime, require_session)
    _install_project_routes(app, runtime, require_session, require_mutation)
    _install_command_routes(app, runtime, require_session, require_mutation)
    _install_run_routes(app, runtime, require_session, require_mutation)
    _install_control_routes(app, runtime, require_session, require_mutation)
    install_ai_routes(app, runtime, require_session, require_mutation)
    return app


def _security_dependencies(runtime: StudioRuntime) -> tuple[Dependency, Dependency]:
    async def require_session(request: Request) -> Session:
        session = runtime.sessions.get(request.cookies.get(SESSION_COOKIE))
        if session is None:
            raise HTTPException(status.HTTP_401_UNAUTHORIZED, "Studio session required")
        return session

    async def require_mutation(
        request: Request,
        session: Session = Depends(require_session),
        csrf_token: str | None = Header(default=None, alias=CSRF_HEADER),
    ) -> Session:
        if not origin_is_allowed(
            request.headers.get("origin"), request.state.studio_origin
        ):
            raise HTTPException(status.HTTP_403_FORBIDDEN, "Origin rejected")
        if not csrf_is_allowed(csrf_token, session):
            raise HTTPException(status.HTTP_403_FORBIDDEN, "CSRF token rejected")
        return session

    return require_session, require_mutation


def _install_security_middleware(app: FastAPI) -> None:
    @app.middleware("http")
    async def secure_loopback_requests(request: Request, call_next):
        hosts = request.headers.getlist("host")
        origin = loopback_origin_from_host(hosts[0]) if len(hosts) == 1 else None
        if origin is None:
            return HTMLResponse("Invalid Host", status_code=status.HTTP_400_BAD_REQUEST)
        request.state.studio_origin = origin
        response = await call_next(request)
        response.headers.setdefault("Content-Security-Policy", CONTENT_SECURITY_POLICY)
        response.headers["Referrer-Policy"] = "no-referrer"
        response.headers["X-Content-Type-Options"] = "nosniff"
        response.headers["X-Frame-Options"] = "DENY"
        response.headers["Cache-Control"] = "no-store"
        return response


def _install_page_routes(
    app: FastAPI, runtime: StudioRuntime, require_session: Dependency
) -> None:
    template_root = files("mokume.studio").joinpath("templates")
    static_root = files("mokume.studio").joinpath("static")
    templates = Jinja2Templates(directory=str(template_root))

    @app.get("/api/health")
    async def health() -> dict:
        package = importlib.import_module("mokume")
        return {"status": "ready", "version": package.__version__}

    @app.get("/", response_class=HTMLResponse)
    async def index(request: Request, token: str | None = Query(default=None)):
        if token is not None:
            session = runtime.sessions.exchange_startup_token(token)
            if session is None:
                raise HTTPException(
                    status.HTTP_401_UNAUTHORIZED, "Invalid startup token"
                )
            response = RedirectResponse("/", status_code=status.HTTP_303_SEE_OTHER)
            response.set_cookie(
                SESSION_COOKIE,
                session.id,
                httponly=True,
                samesite="strict",
                secure=False,
                path="/",
            )
            return response
        await require_session(request)
        return templates.TemplateResponse(request=request, name="index.html")

    @app.get("/static/{asset}")
    async def static_asset(asset: str, _session: Session = Depends(require_session)):
        if asset not in {
            "mokume-favicon.png",
            "mokume-logo.png",
            "mokume-mark.png",
            "studio.css",
            "studio.js",
        }:
            raise HTTPException(status.HTTP_404_NOT_FOUND)
        return FileResponse(str(static_root.joinpath(asset)))

    @app.get("/api/session")
    async def session_info(
        request: Request,
        session: Session = Depends(require_session),
    ) -> dict:
        package = importlib.import_module("mokume")
        provider = runtime.providers.summary(session.id)
        return {
            "csrf_token": session.csrf_token,
            "version": package.__version__,
            "origin": request.state.studio_origin,
            "ai_configured": provider is not None,
            "ai_provider": provider,
        }


def _install_project_routes(
    app: FastAPI,
    runtime: StudioRuntime,
    require_session: Dependency,
    require_mutation: Dependency,
) -> None:
    @app.get("/api/project")
    async def project(_session: Session = Depends(require_session)):
        return runtime.store.active_project()

    @app.get("/api/folders")
    async def folders(
        _session: Session = Depends(require_session),
        path: str | None = Query(default=None),
    ) -> dict:
        try:
            directory, entries = list_directories(path)
        except PathAccessError as exc:
            raise HTTPException(status.HTTP_400_BAD_REQUEST, str(exc)) from exc
        parent = directory.parent if directory.parent != directory else None
        return {
            "path": str(directory),
            "parent": str(parent) if parent else None,
            "directories": [entry.model_dump() for entry in entries],
        }

    @app.post("/api/projects/open")
    async def open_project(
        payload: OpenProjectRequest,
        _session: Session = Depends(require_mutation),
    ):
        _require_idle(runtime.store)
        try:
            root = readable_directory(payload.path)
        except PathAccessError as exc:
            raise HTTPException(status.HTTP_400_BAD_REQUEST, str(exc)) from exc
        return runtime.store.open_project(str(root))

    @app.post("/api/projects/close", status_code=status.HTTP_204_NO_CONTENT)
    async def close_project(
        _session: Session = Depends(require_mutation),
    ) -> None:
        _require_idle(runtime.store)
        runtime.store.close_project()

    @app.get("/api/files")
    async def project_files(
        _session: Session = Depends(require_session),
        path: str = Query(default="."),
    ) -> dict:
        guard = active_project_paths(runtime.store)
        try:
            directory = guard.resolve_existing(path)
            if not directory.is_dir():
                raise PathAccessError(f"not a directory: {path}")
            entries = list_project_entries(guard, directory)
        except PathAccessError as exc:
            raise HTTPException(status.HTTP_400_BAD_REQUEST, str(exc)) from exc
        parent = directory.parent if directory != guard.root else None
        return {
            "path": guard.relative(directory),
            "parent": guard.relative(parent) if parent else None,
            "entries": entries,
        }


def _install_command_routes(
    app: FastAPI,
    runtime: StudioRuntime,
    require_session: Dependency,
    require_mutation: Dependency,
) -> None:
    install_workflow_template_routes(
        app, runtime.store, require_session, require_mutation
    )

    @app.get("/api/commands")
    async def commands(_session: Session = Depends(require_session)):
        return {"commands": command_schema()}

    @app.post("/api/commands/validate")
    async def validate_command(
        payload: ValidationRequest,
        _session: Session = Depends(require_mutation),
    ) -> dict:
        project = runtime.store.active_project()
        if project is None:
            raise HTTPException(status.HTTP_409_CONFLICT, "Open a project folder first")
        try:
            argv = validate_and_canonicalize(payload.argv, project.root)
        except (
            CommandValidationError,
            PathAccessError,
            RuntimeError,
            ValueError,
        ) as exc:
            raise HTTPException(status.HTTP_422_UNPROCESSABLE_ENTITY, str(exc)) from exc
        return {"valid": True, "argv": argv}


def _install_run_routes(
    app: FastAPI,
    runtime: StudioRuntime,
    require_session: Dependency,
    require_mutation: Dependency,
) -> None:
    _install_run_mutation_routes(app, runtime, require_mutation)
    _install_run_read_routes(app, runtime, require_session)


def _install_run_mutation_routes(
    app: FastAPI,
    runtime: StudioRuntime,
    require_mutation: Dependency,
) -> None:
    @app.post("/api/runs", status_code=status.HTTP_202_ACCEPTED)
    async def create_run(
        payload: RunRequest,
        _session: Session = Depends(require_mutation),
    ):
        project = runtime.store.active_project()
        if project is None:
            raise HTTPException(status.HTTP_409_CONFLICT, "Open a project folder first")
        try:
            return await asyncio.to_thread(runtime.jobs.submit, payload, project)
        except JobConflictError as exc:
            raise HTTPException(status.HTTP_409_CONFLICT, str(exc)) from exc
        except (
            CommandValidationError,
            PathAccessError,
            RuntimeError,
            ValueError,
        ) as exc:
            raise HTTPException(status.HTTP_422_UNPROCESSABLE_ENTITY, str(exc)) from exc

    @app.post("/api/runs/{run_id}/cancel")
    async def cancel_run(
        run_id: str,
        _session: Session = Depends(require_mutation),
    ):
        try:
            return await asyncio.to_thread(runtime.jobs.cancel, run_id)
        except KeyError as exc:
            raise HTTPException(status.HTTP_404_NOT_FOUND, "Run not found") from exc


def _install_run_read_routes(
    app: FastAPI,
    runtime: StudioRuntime,
    require_session: Dependency,
) -> None:

    @app.get("/api/runs")
    async def runs(_session: Session = Depends(require_session)) -> dict:
        project = runtime.store.active_project()
        return {
            "runs": []
            if project is None
            else runtime.store.list_runs(project_id=project.id)
        }

    @app.get("/api/runs/{run_id}")
    async def run(run_id: str, _session: Session = Depends(require_session)):
        record = runtime.store.get_run(run_id)
        if record is None:
            raise HTTPException(status.HTTP_404_NOT_FOUND, "Run not found")
        return record

    @app.get("/api/runs/{run_id}/events")
    async def run_events(
        run_id: str,
        _session: Session = Depends(require_session),
        last_event_id: str | None = Header(default=None, alias="Last-Event-ID"),
    ):
        if runtime.store.get_run(run_id) is None:
            raise HTTPException(status.HTTP_404_NOT_FOUND, "Run not found")
        try:
            cursor = int(last_event_id or 0)
        except ValueError as exc:
            raise HTTPException(
                status.HTTP_400_BAD_REQUEST, "Invalid Last-Event-ID"
            ) from exc
        return StreamingResponse(
            _event_stream(runtime.store, run_id, cursor), media_type="text/event-stream"
        )


def _install_control_routes(
    app: FastAPI,
    runtime: StudioRuntime,
    require_session: Dependency,
    require_mutation: Dependency,
) -> None:
    @app.get("/api/artifacts")
    async def artifacts(
        _session: Session = Depends(require_session),
        run_id: str | None = Query(default=None),
    ) -> dict:
        return {"artifacts": runtime.store.list_artifacts(run_id)}

    @app.get("/api/artifacts/{artifact_id}")
    async def artifact(
        artifact_id: str,
        _session: Session = Depends(require_session),
    ):
        record = runtime.store.get_artifact(artifact_id)
        if record is None:
            raise HTTPException(status.HTTP_404_NOT_FOUND, "Artifact not found")
        run_record = runtime.store.get_run(record.run_id)
        project = (
            runtime.store.get_project(run_record.project_id) if run_record else None
        )
        if project is None:
            raise HTTPException(status.HTTP_404_NOT_FOUND, "Artifact project not found")
        try:
            path = ProjectPaths(project.root).resolve_existing(record.path)
        except PathAccessError as exc:
            raise HTTPException(status.HTTP_404_NOT_FOUND, str(exc)) from exc
        disposition = "attachment" if record.media_type == "text/html" else "inline"
        return FileResponse(
            path,
            media_type=record.media_type,
            headers={"Content-Disposition": f'{disposition}; filename="{path.name}"'},
        )

    @app.get("/api/system")
    async def system_status(
        request: Request,
        _session: Session = Depends(require_session),
    ) -> dict:
        memory = psutil.virtual_memory()
        return {
            "origin": request.state.studio_origin,
            "state_directory": str(runtime.store.directory),
            "project": runtime.store.active_project(),
            "runs": len(runtime.store.list_runs()),
            "threads": 24,
            "memory": {
                "total_bytes": memory.total,
                "available_bytes": memory.available,
                "used_bytes": memory.total - memory.available,
                "percent": round(memory.percent, 1),
            },
        }

    @app.post("/api/studio/exit", status_code=status.HTTP_202_ACCEPTED)
    async def exit_studio(
        _session: Session = Depends(require_mutation),
    ) -> dict:
        _require_idle(runtime.store)
        if runtime.shutdown_callback is None:
            raise HTTPException(
                status.HTTP_501_NOT_IMPLEMENTED, "Shutdown is unavailable"
            )
        runtime.shutdown_callback()
        return {"status": "stopping"}


def _require_idle(store: StateStore) -> None:
    if store.has_active_run():
        raise HTTPException(
            status.HTTP_409_CONFLICT,
            "cannot change project while a run is active",
        )


async def _event_stream(store: StateStore, run_id: str, cursor: int):
    """Yield resumable server-sent events until the run reaches a terminal state."""
    while True:
        events = store.events_after(run_id, cursor)
        for event in events:
            cursor = event["sequence"]
            payload = json.dumps(event["payload"], separators=(",", ":"))
            yield f"id: {cursor}\nevent: {event['type']}\ndata: {payload}\n\n"
        record = store.get_run(run_id)
        if record is None or (record.status in TERMINAL_RUN_STATUSES and not events):
            break
        await asyncio.sleep(0.4)
