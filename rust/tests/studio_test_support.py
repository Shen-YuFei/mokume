"""Shared construction helpers for optional Studio HTTP tests."""

import secrets
from pathlib import Path
from unittest.mock import patch

import httpx
from fastapi import FastAPI
from pydantic_ai.models.test import TestModel

from mokume.studio.app import create_app
from mokume.studio.providers import ProviderExecution, ProviderRegistry


AI_PORT = 18766
AI_ORIGIN = f"http://127.0.0.1:{AI_PORT}"
AI_TOKEN = "studio-ai-startup-token"
AI_SECRET = secrets.token_urlsafe(32)


def make_studio_app(token: str, state_directory: Path) -> FastAPI:
    """Create a testing-mode Studio app with one isolated state directory."""
    with patch(
        "mokume.studio.app.default_provider_config_root",
        return_value=state_directory,
    ):
        return create_app(
            startup_token=token,
            state_directory=state_directory,
        )


async def authenticated_project(
    client: httpx.AsyncClient, project: Path
) -> dict[str, str]:
    """Authenticate one client and open its test workspace."""
    response = await client.get(f"/?token={AI_TOKEN}", follow_redirects=False)
    assert response.status_code == 303
    session = (await client.get("/api/session")).json()
    headers = {"Origin": AI_ORIGIN, "X-CSRF-Token": session["csrf_token"]}
    opened = await client.post(
        "/api/projects/open",
        json={"path": str(project)},
        headers=headers,
    )
    assert opened.status_code == 200
    return headers


async def configure_provider(
    client: httpx.AsyncClient, headers: dict[str, str]
) -> httpx.Response:
    """Install the shared offline provider configuration for a test client."""
    response = await client.post(
        "/api/ai/config",
        json={
            "provider": "openai-responses",
            "model": "offline",
            "api_key": AI_SECRET,
        },
        headers=headers,
    )
    assert response.status_code == 200, response.text
    return response


def agent_body(**updates):
    """Return the minimum trusted AG-UI request envelope."""
    body = {
        "threadId": "ask-thread",
        "runId": "ask-run",
        "state": None,
        "messages": [{"id": "message-1", "role": "user", "content": "Help"}],
        "tools": [],
        "context": [],
        "forwardedProps": {"mode": "ask", "datasetId": None},
    }
    body.update(updates)
    return body


async def offline_ask(client, project_path: Path, monkeypatch):
    """Run one server-tool Ask exchange without a network model."""
    headers = await authenticated_project(client, project_path)
    await configure_provider(client, headers)
    project_id = (await client.get("/api/project")).json()["id"]
    offline = TestModel(call_tools=["get_analysis_context"])
    monkeypatch.setattr(
        ProviderRegistry,
        "execution_for",
        lambda _self, _session: ProviderExecution(offline, None, None),
    )
    response = await client.post(
        "/api/agent/run",
        json=agent_body(
            forwardedProps={
                "mode": "ask",
                "datasetId": None,
                "projectId": project_id,
            }
        ),
        headers={**headers, "Accept": "text/event-stream"},
    )
    return response, headers, project_id


async def assert_static_assets_not_cached(app) -> None:
    """Check that local refreshes load the current frontend assets."""
    route = next(
        route
        for route in app.routes
        if getattr(route, "path", None) == "/static/{asset}"
    )
    response = await route.endpoint("studio.js", _session=object())
    assert response.headers["cache-control"] == "no-store"


def assert_theme_and_layout_scripts(script_text: str, stylesheet_text: str) -> None:
    """Check the frontend behavior and layout hooks used by Studio."""
    _assert_theme_scripts(script_text, stylesheet_text)
    _assert_assistant_scripts(script_text, stylesheet_text)
    _assert_workbench_scripts(script_text, stylesheet_text)


def _assert_theme_scripts(script_text: str, stylesheet_text: str) -> None:
    """Check persisted appearance and conditional-field hooks."""
    assert 'const LANGUAGE_STORAGE_KEY = "mokume:language"' in script_text
    assert 'const APPEARANCE_STORAGE_KEY = "mokume:appearance"' in script_text
    assert "document.documentElement.lang = state.language" in script_text
    assert "document.documentElement.dataset.theme" in script_text
    assert ':root[data-theme="light"]' in stylesheet_text
    assert ".form-field[hidden] { display: none; }" in stylesheet_text


def _assert_assistant_scripts(script_text: str, stylesheet_text: str) -> None:
    """Check assistant streaming, Markdown, and composer layout hooks."""
    assert "async function renderAssistantMarkdown" in script_text
    assert "function appendAssistantThinking" in script_text
    assert "clearAssistantThinking(stream)" in script_text
    assert ".assistant-thinking-dot" in stylesheet_text
    assert "function appendAssistantReasoning" in script_text
    assert "function adjacentAssistantReasoning" in script_text
    assert "adjacentAssistantReasoning()" in script_text
    assert 'event.type === "REASONING_MESSAGE_CONTENT"' in script_text
    assert "reasoning: null" in script_text
    assert "reasoningMessageIds: new Set()" in script_text
    assert "stream.reasoningMessageIds.has(event.messageId)" in script_text
    assert ".assistant-reasoning-details" in stylesheet_text
    assert ".assistant-reasoning-details summary::before" not in stylesheet_text
    assert "function appendAssistantError" in script_text
    assert "error.assistantDisplayed = true" in script_text
    assert ".assistant-error .message-body" in stylesheet_text
    assert 'event.type === "TEXT_MESSAGE_END"' in script_text
    assert ".markdown-body table" in stylesheet_text
    assert (
        ".assistant-messages { flex: 1; min-height: 0; overflow-x: hidden;"
        in stylesheet_text
    )
    assert (
        ".markdown-body li > strong:first-child { white-space: normal; }"
        in stylesheet_text
    )
    assert ".markdown-body :not(pre) > code" in stylesheet_text
    assert ".markdown-body pre { max-width: 100%;" in stylesheet_text
    assert "table-layout: fixed;" in stylesheet_text
    assert "async function refreshSystemMemory" in script_text
    assert "SYSTEM_MEMORY_REFRESH_MS = 5000" in script_text
    assert ".system-memory-value" in stylesheet_text
    assert "padding-bottom: 2px;" in stylesheet_text
    assert "function providerThinkingLevel" in script_text
    assert "function restoreProviderThinkingLevel" in script_text
    assert 'value === "custom"' in script_text
    assert ".composer-context { display: flex; flex: 1 1 0;" in stylesheet_text
    assert ".composer-actions .assistant-model-button" in stylesheet_text
    assert "flex: 1 1 0; min-width: 0; max-width: 180px;" in stylesheet_text


def _assert_workbench_scripts(script_text: str, stylesheet_text: str) -> None:
    """Check menu, panel, review, queue, and QC hooks."""
    assert "function openSubmenu" in script_text
    assert (
        'trigger.addEventListener("mouseenter", () => openSubmenu(trigger))'
        in script_text
    )
    assert 'menu.addEventListener("mouseenter", cancelSubmenuClose)' in script_text
    assert "const PANEL_SIZE_LIMITS" in script_text
    assert "function bindPanelResizers" in script_text
    assert 'handle.addEventListener("pointerdown"' in script_text
    assert "max-width: min(420px, calc(100vw - 36px));" in stylesheet_text
    assert "overflow-wrap: anywhere; word-break: break-word;" in stylesheet_text
    assert "async function openCommandReview" in script_text
    assert "async function openRunDetails" in script_text
    assert "function executionStageSection" in script_text
    assert "async function openRunComparison" in script_text
    assert "async function processBatchQueue" in script_text
    assert "function renderQcPanel" in script_text
    assert ".stage-stepper" in stylesheet_text
    assert ".qc-metrics" in stylesheet_text
