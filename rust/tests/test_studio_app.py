"""Security and local-project contracts for Mokume Studio."""

from __future__ import annotations

from importlib.resources import files
from types import SimpleNamespace

import httpx
import pytest
import studio_test_support
from studio_test_support import assert_theme_and_layout_scripts, make_studio_app

from mokume.studio.models import JobSpec, RunStatus, utc_now
from mokume.studio.paths import PathAccessError, ProjectPaths
from mokume.studio.state import StateStore


PORT = 18765
ORIGIN = f"http://127.0.0.1:{PORT}"
TOKEN = "test-startup-token"
pytestmark = pytest.mark.anyio


@pytest.fixture(name="studio_client")
async def build_studio_client(tmp_path):
    """Create an in-memory authenticated-capable Studio HTTP client."""
    app = make_studio_app(TOKEN, tmp_path / "state")
    transport = httpx.ASGITransport(app=app)
    async with httpx.AsyncClient(transport=transport, base_url=ORIGIN) as test_client:
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
    async with httpx.AsyncClient(transport=transport, base_url=ORIGIN) as second_client:
        assert (await second_client.get(f"/?token={TOKEN}")).status_code == 401


async def test_control_api_requires_session(studio_client):
    """Control-plane endpoints reject requests without a Studio session."""
    assert (await studio_client.get("/api/project")).status_code == 401
    assert (await studio_client.get("/api/commands")).status_code == 401


async def test_static_assets_are_not_cached_between_local_updates(studio_client):
    """A refresh must load the current frontend instead of stale JavaScript."""
    await studio_test_support.assert_static_assets_not_cached(studio_client.app)


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


@pytest.mark.parametrize(
    "host",
    [
        "attacker.invalid:57828",
        "192.168.1.10:57828",
        "localhost",
        "localhost:0",
        "localhost:65536",
        "localhost:not-a-port",
        "[::1%25lo]:57828",
    ],
)
async def test_invalid_host_is_rejected(studio_client, host):
    """Host-header attacks are rejected before route dispatch."""
    response = await studio_client.get("/api/health", headers={"Host": host})
    assert response.status_code == 400


async def test_forwarded_loopback_host_uses_its_own_origin(tmp_path):
    """SSH and IDE loopback forwarding may replace the browser-side port."""
    forwarded_origin = "http://localhost:57828"
    app = make_studio_app(TOKEN, tmp_path / "forwarded-state")
    transport = httpx.ASGITransport(app=app)
    async with httpx.AsyncClient(
        transport=transport,
        base_url=forwarded_origin,
    ) as client:
        csrf = await authenticate(client)
        session = await client.get("/api/session")
        assert session.json()["origin"] == forwarded_origin
        rejected = await client.post(
            "/api/projects/open",
            json={"path": str(tmp_path)},
            headers={"Origin": ORIGIN, "X-CSRF-Token": csrf},
        )
        accepted = await client.post(
            "/api/projects/open",
            json={"path": str(tmp_path)},
            headers={"Origin": forwarded_origin, "X-CSRF-Token": csrf},
        )
        system = await client.get("/api/system")

    assert rejected.status_code == 403
    assert accepted.status_code == 200
    assert system.json()["origin"] == forwarded_origin
    async with httpx.AsyncClient(
        transport=httpx.ASGITransport(app=app),
        base_url=ORIGIN,
    ) as ipv6_client:
        ipv6 = await ipv6_client.get("/api/health", headers={"Host": "[::1]:57828"})
        assert ipv6.status_code == 200


async def test_system_status_reports_real_memory_usage(studio_client, monkeypatch):
    """Expose consistent host-memory totals for the header indicator."""
    monkeypatch.setattr(
        "mokume.studio.app.psutil.virtual_memory",
        lambda: SimpleNamespace(
            total=16 * 1024**3,
            available=6 * 1024**3,
            percent=62.5,
        ),
    )
    await authenticate(studio_client)

    response = await studio_client.get("/api/system")

    assert response.status_code == 200
    assert response.json()["memory"] == {
        "total_bytes": 16 * 1024**3,
        "available_bytes": 6 * 1024**3,
        "used_bytes": 10 * 1024**3,
        "percent": 62.5,
    }


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


async def test_project_files_list_nested_directories(studio_client, tmp_path):
    """The file tree can lazily request children beneath the project root."""
    project = tmp_path / "project"
    nested = project / "data" / "raw"
    nested.mkdir(parents=True)
    (nested / "sample.tsv").write_text("value\n", encoding="utf-8")
    (project / "z-folder").mkdir()
    (project / "a-file.txt").write_text("value\n", encoding="utf-8")
    csrf = await authenticate(studio_client)
    opened = await studio_client.post(
        "/api/projects/open",
        json={"path": str(project)},
        headers=mutation_headers(csrf),
    )
    assert opened.status_code == 200

    root = await studio_client.get("/api/files")
    data = await studio_client.get("/api/files", params={"path": "data"})
    raw = await studio_client.get("/api/files", params={"path": "data/raw"})

    assert [entry["path"] for entry in root.json()["entries"]] == [
        "data",
        "z-folder",
        "a-file.txt",
    ]
    assert data.json()["entries"][0]["path"] == "data/raw"
    assert raw.json()["entries"][0]["path"] == "data/raw/sample.tsv"


async def test_workflow_templates_round_trip_inside_the_project(
    studio_client, tmp_path
):
    """Template reads and writes stay on the Studio host inside the workspace."""
    project = tmp_path / "project"
    templates = project / "templates"
    templates.mkdir(parents=True)
    csrf = await authenticate(studio_client)
    headers = mutation_headers(csrf)
    opened = await studio_client.post(
        "/api/projects/open", json={"path": str(project)}, headers=headers
    )
    assert opened.status_code == 200
    document = {
        "$schemaVersion": 1,
        "workflow": ["quantify", "features2proteins"],
        "parameters": {"parquet": "input.parquet"},
    }

    saved = await studio_client.put(
        "/api/workflow-template",
        json={"path": "templates/protein.json", "template": document},
        headers=headers,
    )
    loaded = await studio_client.get(
        "/api/workflow-template", params={"path": "templates/protein.json"}
    )

    assert saved.status_code == 200
    assert saved.json() == {"path": "templates/protein.json"}
    assert loaded.status_code == 200
    assert loaded.json() == {
        "path": "templates/protein.json",
        "template": document,
    }


async def test_workflow_template_write_rejects_escape_and_unapproved_overwrite(
    studio_client, tmp_path
):
    """Template export cannot escape the workspace or silently replace a file."""
    project = tmp_path / "project"
    project.mkdir()
    existing = project / "workflow.json"
    existing.write_text("{}\n", encoding="utf-8")
    csrf = await authenticate(studio_client)
    headers = mutation_headers(csrf)
    await studio_client.post(
        "/api/projects/open", json={"path": str(project)}, headers=headers
    )
    document = {
        "$schemaVersion": 1,
        "workflow": ["correct-batches"],
        "parameters": {},
    }

    overwrite = await studio_client.put(
        "/api/workflow-template",
        json={"path": "workflow.json", "template": document},
        headers=headers,
    )
    approved_overwrite = await studio_client.put(
        "/api/workflow-template",
        json={
            "path": "workflow.json",
            "template": document,
            "overwrite": True,
        },
        headers=headers,
    )
    escape = await studio_client.put(
        "/api/workflow-template",
        json={"path": "../outside.json", "template": document},
        headers=headers,
    )

    assert overwrite.status_code == 422
    assert approved_overwrite.status_code == 200
    assert '"correct-batches"' in existing.read_text(encoding="utf-8")
    assert escape.status_code == 422
    assert not (tmp_path / "outside.json").exists()


async def test_active_run_blocks_project_switches(studio_client, tmp_path):
    """A queued worker keeps both project mutations locked."""
    csrf = await authenticate(studio_client)
    headers = mutation_headers(csrf)
    opened = await studio_client.post(
        "/api/projects/open", json={"path": str(tmp_path)}, headers=headers
    )
    project = opened.json()
    store = studio_client.app.state.runtime.store
    store.create_run(
        JobSpec(
            run_id="active-run",
            project_root=project["root"],
            run_directory=str(tmp_path / "active-run"),
            argv=["correct-batches"],
            parameters={},
            approved_hash="hash",
            created_at=utc_now(),
        ),
        project["id"],
        "correct-batches",
    )

    reopened = await studio_client.post(
        "/api/projects/open", json={"path": str(tmp_path)}, headers=headers
    )
    closed = await studio_client.post("/api/projects/close", headers=headers)
    exited = await studio_client.post("/api/studio/exit", headers=headers)

    assert reopened.status_code == 409
    assert closed.status_code == 409
    assert exited.status_code == 409
    assert store.active_project().id == project["id"]


async def test_run_history_is_scoped_to_the_active_project(studio_client, tmp_path):
    """Switching projects must not expose the previous project's run history."""
    csrf = await authenticate(studio_client)
    headers = mutation_headers(csrf)
    store = studio_client.app.state.runtime.store
    project_a_root = tmp_path / "project-a"
    project_b_root = tmp_path / "project-b"
    project_a_root.mkdir()
    project_b_root.mkdir()

    opened_a = await studio_client.post(
        "/api/projects/open",
        json={"path": str(project_a_root)},
        headers=headers,
    )
    project_a = opened_a.json()
    store.create_run(
        JobSpec(
            run_id="project-a-run",
            project_root=project_a["root"],
            run_directory=str(project_a_root / "run"),
            argv=["quantify", "peptides2protein"],
            parameters={},
            approved_hash="hash-a",
            created_at=utc_now(),
        ),
        project_a["id"],
        "quantify peptides2protein",
    )
    store.update_run("project-a-run", RunStatus.SUCCEEDED)

    runs_a = await studio_client.get("/api/runs")
    assert [run["id"] for run in runs_a.json()["runs"]] == ["project-a-run"]

    opened_b = await studio_client.post(
        "/api/projects/open",
        json={"path": str(project_b_root)},
        headers=headers,
    )
    assert opened_b.status_code == 200
    runs_b = await studio_client.get("/api/runs")
    assert runs_b.json() == {"runs": []}
    assert (await studio_client.get("/api/artifacts")).json() == {"artifacts": []}


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


async def test_state_restart_clears_an_unavailable_active_project(tmp_path):
    """A deleted temporary workspace must not block choosing a new folder."""
    state_directory = tmp_path / "state"
    project_root = tmp_path / "deleted-project"
    project_root.mkdir()
    StateStore(state_directory).open_project(str(project_root.resolve()))
    project_root.rmdir()

    app = make_studio_app(TOKEN, state_directory)
    transport = httpx.ASGITransport(app=app)
    async with httpx.AsyncClient(transport=transport, base_url=ORIGIN) as client:
        await authenticate(client)
        response = await client.get("/api/project")

    assert response.status_code == 200
    assert response.json() is None


def _assert_packaged_assets(package) -> None:
    assert package.joinpath("templates/index.html").is_file()
    assert package.joinpath("static/studio.css").is_file()
    assert package.joinpath("static/studio.js").is_file()
    assert package.joinpath("static/mokume-favicon.png").is_file()
    assert package.joinpath("static/mokume-logo.png").is_file()
    assert package.joinpath("static/mokume-mark.png").is_file()
    assert package.joinpath("static/LUCIDE_LICENSE.txt").is_file()
    assert package.joinpath("static/VSCODE_ICONS_LICENSE.txt").is_file()


def _assert_branding_and_menus(page: str) -> None:
    assert (
        '<link rel="icon" type="image/png" sizes="256x256" '
        'href="/static/mokume-favicon.png?v=transparent">'
    ) in page
    assert page.count('src="/static/mokume-mark.png"') == 1
    assert 'class="brand-logo"' not in page
    assert 'src="/static/mokume-logo.png"' not in page
    assert "<span>Mokume Studio</span>" in page
    assert '<span class="brand-mark">M</span>' not in page
    assert '<div class="welcome-mark">M</div>' not in page
    assert page.count('id="system-memory" class="system-memory unavailable"') == 1
    assert 'role="meter"' in page
    assert (
        page.index('id="project-chip"')
        < page.index('id="system-memory"')
        < page.index("</header>")
    )
    for menu in ("File", "Analysis", "View", "Help"):
        assert f">{menu}<" in page
    assert ">Edit<" not in page
    assert 'role="menuitemradio" data-language="zh-CN"' in page
    assert 'data-submenu="appearance-submenu"' in page
    assert 'data-submenu="language-submenu"' in page
    for icon in (
        "folder",
        "vscode-parquet",
        "vscode-config",
        "sdrf",
        "dna",
        "mass-spectrum",
        "feature-map",
        "msstats",
        "workflow",
    ):
        assert f'<symbol id="file-icon-{icon}"' in page

    assert 'data-action="import-workflow-template"' in page
    assert 'data-action="export-workflow-template"' in page
    assert 'id="workflow-template-dialog"' in page
    assert 'id="workflow-template-file" type="file"' not in page
    assert 'id="workflow-template-current"' in page
    assert 'id="workflow-template-name" type="text"' in page
    assert 'data-action="queue-command"' in page
    assert 'data-bottom-tab="queue"' in page
    assert 'data-bottom-tab="qc"' in page
    assert 'id="command-review-dialog"' in page
    assert 'id="run-details-dialog"' in page
    assert 'id="run-compare-dialog"' in page


def _assert_assistant_controls(page: str) -> None:
    assert 'class="assistant-chat-header"' in page
    assert 'class="assistant-title-chip"' not in page
    assert 'data-action="new-assistant-chat"' in page
    assert 'data-action="assistant-history"' in page
    assert (
        page.index('data-action="new-assistant-chat"')
        < page.index('data-action="assistant-history"')
        < page.index('id="configure-provider"')
    )
    assert 'data-action="inspect-dataset"' not in page
    assert 'id="dataset-dialog"' not in page
    assert 'id="dataset-form"' not in page
    assert 'id="conversation-dialog"' in page
    assert 'data-i18n-title="Collapse sidebar"' in page
    assert 'data-action="refresh-files"' in page
    assert 'data-action="collapse-folders"' in page
    assert 'data-i18n-title="Collapse bottom panel"' in page
    assert 'class="panel-restore-button sidebar-restore-button"' in page
    assert 'class="panel-restore-button assistant-restore-button"' in page
    assert 'id="assistant-mode-menu" class="assistant-mode-menu"' in page
    assert '<select id="assistant-mode"' not in page
    assert 'data-assistant-mode="ask"' in page
    assert 'data-assistant-mode="agent"' in page
    assert 'data-assistant-mode="plan"' not in page
    assert 'class="assistant-notice"' not in page
    assert 'id="assistant-model-name"' in page
    assert 'id="assistant-icon-ask"' in page
    assert 'id="assistant-icon-agent"' in page
    assert 'class="primary assistant-send-button"' in page
    assert ">Send</button>" not in page
    assert "Enter to send · Shift+Enter for a new line" not in page
    assert 'id="provider-persist" type="checkbox"' in page
    assert 'data-i18n="Save">Save</button>' in page
    assert "Save for Session" not in page
    thinking_options = ["default", "low", "medium", "high", "xhigh", "max", "custom"]
    thinking_positions = [
        page.index(f">{option}</option>") for option in thinking_options
    ]
    assert thinking_positions == sorted(thinking_positions)
    assert '<option value="off"' not in page
    assert '<option value="minimal"' not in page
    assert 'id="provider-thinking-custom"' in page
    assert 'placeholder="custom level"' in page
    for panel in ("sidebar", "assistant", "bottom"):
        assert f'data-resize-panel="{panel}"' in page
    for appearance in ("system", "light", "dark"):
        assert f'data-appearance="{appearance}"' in page


def _assert_workspace_scripts(script_text: str) -> None:
    assert 'row.addEventListener("dblclick"' not in script_text
    assert (
        'row.addEventListener("click", () => loadFolders(directory.path)' in script_text
    )
    assert "const projectRoot = state.project?.root || null;" in script_text
    assert "if (!projectRoot) throw error;" in script_text
    assert "await loadFolders();" in script_text
    assert "async function toggleDirectory" in script_text
    assert "function collapseFileTree" in script_text
    assert "async function openConversationHistory" in script_text
    assert "async function openStoredConversation" in script_text
    assert "async function renameStoredConversation" in script_text
    assert "async function deleteStoredConversation" in script_text
    assert 'className = "conversation-rename"' in script_text
    assert 'className = "conversation-delete"' in script_text
    assert 'className = "conversation-workspace"' in script_text
    assert '"Workspace: {path}", { path: workspace.root }' in script_text
    assert "function openDatasetDialog" not in script_text
    assert "function inspectDataset" not in script_text
    assert 'api("/api/datasets/inspect"' not in script_text
    assert "projectId: state.projectId" in script_text
    assert "const expandedPaths = expandedDirectoryPaths(tree);" in script_text
    assert "await restoreExpandedDirectories(tree, expandedPaths);" in script_text
    assert "row.dataset.path = entry.path;" in script_text
    assert "`/api/files?path=${encodeURIComponent(path)}`" in script_text
    assert 'openWorkflowTemplateDialog("import")' in script_text
    assert 'openWorkflowTemplateDialog("export")' in script_text
    assert "formState: agentWorkflowFormState()" in script_text
    assert 'event.type === "TOOL_CALL_RESULT"' in script_text
    assert "applyWorkflowParameterPatch(result)" in script_text
    assert 'group.className = "command-family"' in script_text
    assert "button.textContent = command.display_name" in script_text
    assert "function updateDifferentialExpressionPlotParameters" in script_text
    assert 'sampleNormalization === "condition-median"' in script_text
    numeric_block = script_text.split("const NUMERIC_VALUE_FLAGS", maxsplit=1)[1].split(
        "function commandArgv", maxsplit=1
    )[0]
    assert '"impute-shift"' in numeric_block
    assert '"irs-sdrf-value"' not in numeric_block
    assert "URL.createObjectURL" not in script_text


async def test_wheel_resources_and_menu_bar_are_present(studio_client):
    """Packaged Studio resources expose the intentional workbench menu bar."""
    await authenticate(studio_client)
    package = files("mokume.studio")
    page = (await studio_client.get("/")).text
    stylesheet_text = package.joinpath("static/studio.css").read_text(encoding="utf-8")
    script_text = package.joinpath("static/studio.js").read_text(encoding="utf-8")
    _assert_packaged_assets(package)
    _assert_branding_and_menus(page)
    _assert_assistant_controls(page)
    assert_theme_and_layout_scripts(script_text, stylesheet_text)
    _assert_workspace_scripts(script_text)
