"""Security and local-project contracts for Mokume Studio."""

from __future__ import annotations

from importlib.resources import files

import httpx
import pytest

from mokume.studio.app import create_app
from mokume.studio.models import JobSpec, utc_now
from mokume.studio.paths import PathAccessError, ProjectPaths
from mokume.studio.state import StateStore


PORT = 18765
ORIGIN = f"http://127.0.0.1:{PORT}"
TOKEN = "test-startup-token"
pytestmark = pytest.mark.anyio


@pytest.fixture
def anyio_backend():
    """Run the async API tests on the asyncio backend."""
    return "asyncio"


@pytest.fixture(name="studio_client")
async def build_studio_client(tmp_path):
    """Create an in-memory authenticated-capable Studio HTTP client."""
    app = create_app(
        port=PORT,
        startup_token=TOKEN,
        state_directory=tmp_path / "state",
        testing=True,
    )
    transport = httpx.ASGITransport(app=app)
    async with httpx.AsyncClient(
        transport=transport, base_url="http://testserver"
    ) as test_client:
        test_client.app = app
        yield test_client


async def authenticate(studio_client: httpx.AsyncClient) -> str:
    """Exchange the one-time token and return the session CSRF token."""
    response = await studio_client.get(f"/?token={TOKEN}", follow_redirects=False)
    assert response.status_code == 303
    session = await studio_client.get("/api/session")
    assert session.status_code == 200
    return session.json()["csrf_token"]


def mutation_headers(csrf: str) -> dict[str, str]:
    """Build headers required by a state-changing Studio request."""
    return {"Origin": ORIGIN, "X-CSRF-Token": csrf}


async def test_health_is_ready_without_exposing_control_state(studio_client):
    """The public health check reveals no local control-plane paths."""
    response = await studio_client.get("/api/health")
    assert response.status_code == 200
    assert response.json()["status"] == "ready"
    assert "state_directory" not in response.json()


async def test_startup_token_is_one_time_and_becomes_http_only_cookie(studio_client):
    """The startup token is consumed once and becomes an HTTP-only cookie."""
    response = await studio_client.get(f"/?token={TOKEN}", follow_redirects=False)
    assert response.status_code == 303
    cookie = response.headers["set-cookie"]
    assert "HttpOnly" in cookie
    assert "SameSite=strict" in cookie
    assert (await studio_client.get("/")).status_code == 200

    transport = httpx.ASGITransport(app=studio_client.app)
    async with httpx.AsyncClient(
        transport=transport, base_url="http://testserver"
    ) as second_client:
        assert (await second_client.get(f"/?token={TOKEN}")).status_code == 401


async def test_control_api_requires_session(studio_client):
    """Control-plane endpoints reject requests without a Studio session."""
    assert (await studio_client.get("/api/project")).status_code == 401
    assert (await studio_client.get("/api/commands")).status_code == 401


async def test_mutation_requires_exact_origin_and_csrf(studio_client, tmp_path):
    """Mutations require both the exact loopback origin and the CSRF token."""
    csrf = await authenticate(studio_client)
    payload = {"path": str(tmp_path)}

    assert (
        await studio_client.post("/api/projects/open", json=payload)
    ).status_code == 403
    assert (
        await studio_client.post(
            "/api/projects/open",
            json=payload,
            headers={"Origin": ORIGIN, "X-CSRF-Token": "wrong"},
        )
    ).status_code == 403
    response = await studio_client.post(
        "/api/projects/open", json=payload, headers=mutation_headers(csrf)
    )
    assert response.status_code == 200
    assert response.json()["root"] == str(tmp_path.resolve())


async def test_invalid_host_is_rejected(studio_client):
    """Host-header attacks are rejected before route dispatch."""
    response = await studio_client.get(
        "/api/health", headers={"Host": "attacker.invalid"}
    )
    assert response.status_code == 400


async def test_project_guard_rejects_parent_and_symlink_escape(tmp_path):
    """Project paths reject both parent traversal and escaping symlinks."""
    project = tmp_path / "project"
    outside = tmp_path / "outside"
    project.mkdir()
    outside.mkdir()
    (outside / "secret.txt").write_text("secret", encoding="utf-8")
    (project / "escape").symlink_to(outside, target_is_directory=True)
    guard = ProjectPaths(project)

    with pytest.raises(PathAccessError, match="escapes project root"):
        guard.resolve_existing("../outside/secret.txt")
    with pytest.raises(PathAccessError, match="escapes project root"):
        guard.resolve_existing("escape/secret.txt")


async def test_state_restart_marks_active_runs_interrupted(tmp_path):
    """Restart recovery converts nonterminal runs into interrupted records."""
    store = StateStore(tmp_path / "state")
    project = store.open_project(str(tmp_path.resolve()))
    run_directory = tmp_path / "run-1"
    store.create_run(
        JobSpec(
            run_id="run-1",
            project_root=str(tmp_path.resolve()),
            run_directory=str(run_directory),
            argv=["correct-batches"],
            parameters={},
            approved_hash="hash",
            created_at=utc_now(),
        ),
        project.id,
        "correct-batches",
    )

    assert store.interrupt_incomplete_runs() == 1
    assert store.get_run("run-1").status.value == "interrupted"


async def test_wheel_resources_and_menu_bar_are_present(studio_client):
    """Packaged Studio resources expose the intentional workbench menu bar."""
    await authenticate(studio_client)
    package = files("mokume.studio")
    assert package.joinpath("templates/index.html").is_file()
    assert package.joinpath("static/studio.css").is_file()
    assert package.joinpath("static/studio.js").is_file()

    page = (await studio_client.get("/")).text
    for menu in ("File", "Run", "View", "Help"):
        assert f">{menu}<" in page
    assert ">Edit<" not in page
