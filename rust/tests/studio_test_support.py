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
