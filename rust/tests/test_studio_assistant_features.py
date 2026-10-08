"""Assistant form, knowledge, reasoning, and Markdown transport contracts."""

from __future__ import annotations

import json

import pytest
from pydantic_ai.messages import ToolReturnPart
from pydantic_ai.models.function import (
    AgentInfo,
    DeltaThinkingPart,
    DeltaToolCall,
    FunctionModel,
)
from studio_test_support import (
    agent_body,
    authenticated_project,
    configure_provider,
)

from mokume.studio.providers import ProviderExecution, ProviderRegistry


pytestmark = pytest.mark.anyio


async def _parameter_fill_stream(messages, _info: AgentInfo):
    """Request one form patch, then finish after the server tool returns."""
    if any(
        isinstance(part, ToolReturnPart)
        for message in messages
        for part in message.parts
    ):
        yield "The form is updated."
        return
    yield {
        0: DeltaToolCall(
            name="fill_workflow_parameters",
            json_args='{"parameters":{"memory":"8GB"}}',
            tool_call_id="fill-memory",
        )
    }


async def _knowledge_search_stream(messages, _info: AgentInfo):
    """Request one bounded knowledge search, then finish after its result."""
    if any(
        isinstance(part, ToolReturnPart)
        for message in messages
        for part in message.parts
    ):
        yield "The evidence is loaded."
        return
    yield {
        0: DeltaToolCall(
            name="search_knowledge",
            json_args=(
                '{"query":"DirectLFQ LFQ benchmark","data_type":"LFQ",'
                '"method":"directlfq"}'
            ),
            tool_call_id="search-knowledge",
        )
    }


async def _reasoning_stream(_messages, _info: AgentInfo):
    """Return provider-visible reasoning followed by the final answer."""
    yield {0: DeltaThinkingPart(content="Checked the current workspace context.")}
    yield "The analysis is ready."


def _use_stream_model(monkeypatch, stream_function) -> None:
    monkeypatch.setattr(
        ProviderRegistry,
        "execution_for",
        lambda _self, _session: ProviderExecution(
            FunctionModel(stream_function=stream_function), None, None
        ),
    )


def _stream_events(response) -> list[dict]:
    return [
        json.loads(line.removeprefix("data: "))
        for line in response.text.splitlines()
        if line.startswith("data: ")
    ]


def _tool_result(events: list[dict], name: str) -> dict:
    start = next(
        event
        for event in events
        if event["type"] == "TOOL_CALL_START" and event["toolCallName"] == name
    )
    return next(
        event
        for event in events
        if event["type"] == "TOOL_CALL_RESULT"
        and event["toolCallId"] == start["toolCallId"]
    )


def _forwarded(mode: str, project_id: str, form_state=None) -> dict:
    options = {"mode": mode, "datasetId": None, "projectId": project_id}
    if form_state is not None:
        options["formState"] = form_state
    return options


@pytest.mark.parametrize(
    "parameters",
    ({"unknown": "value"}, {"parquet": "../outside.parquet"}),
)
async def test_agent_rejects_untrusted_workflow_form_state(
    ai_client, tmp_path, parameters
):
    """Unknown parameters and escaping paths cannot enter the trusted form state."""
    client, _app = ai_client
    headers = await authenticated_project(client, tmp_path)
    await configure_provider(client, headers)
    project_id = (await client.get("/api/project")).json()["id"]
    form_state = {
        "workflow": ["quantify", "features2proteins"],
        "parameters": parameters,
    }

    response = await client.post(
        "/api/agent/run",
        json=agent_body(forwardedProps=_forwarded("agent", project_id, form_state)),
        headers=headers,
    )

    assert response.status_code == 422


async def test_agent_streams_a_validated_current_workflow_parameter_patch(
    ai_client, tmp_path, monkeypatch
):
    """Agent can fill the current form without receiving a client-owned tool."""
    client, _app = ai_client
    headers = await authenticated_project(client, tmp_path)
    await configure_provider(client, headers)
    project_id = (await client.get("/api/project")).json()["id"]
    _use_stream_model(monkeypatch, _parameter_fill_stream)
    form_state = {
        "workflow": ["quantify", "features2proteins"],
        "parameters": {"threads": "24"},
    }

    response = await client.post(
        "/api/agent/run",
        json=agent_body(
            threadId="agent-fill-thread",
            runId="agent-fill-run",
            forwardedProps=_forwarded("agent", project_id, form_state),
        ),
        headers={**headers, "Accept": "text/event-stream"},
    )

    assert response.status_code == 200, response.text
    events = _stream_events(response)
    result = _tool_result(events, "fill_workflow_parameters")
    assert json.loads(result["content"]) == {
        "type": "workflow_parameter_patch",
        "workflow": ["quantify", "features2proteins"],
        "parameters": {"memory": "8GB"},
    }
    assert any(event["type"] == "RUN_FINISHED" for event in events)


async def test_ask_streams_bounded_server_knowledge(ai_client, tmp_path, monkeypatch):
    """Ask can search server-owned evidence without client-provided tools."""
    client, _app = ai_client
    headers = await authenticated_project(client, tmp_path)
    await configure_provider(client, headers)
    project_id = (await client.get("/api/project")).json()["id"]
    _use_stream_model(monkeypatch, _knowledge_search_stream)

    response = await client.post(
        "/api/agent/run",
        json=agent_body(
            threadId="ask-knowledge-thread",
            runId="ask-knowledge-run",
            forwardedProps=_forwarded("ask", project_id),
        ),
        headers={**headers, "Accept": "text/event-stream"},
    )

    assert response.status_code == 200, response.text
    events = _stream_events(response)
    payload = json.loads(_tool_result(events, "search_knowledge")["content"])
    assert payload["scope"] == "explanation_only"
    assert payload["execution_authority"] is False
    assert payload["count"] == 1
    assert payload["results"][0]["id"] == "grid-lfq-preset"
    assert any(event["type"] == "RUN_FINISHED" for event in events)


async def test_ask_streams_and_persists_provider_visible_reasoning(
    ai_client, tmp_path, monkeypatch
):
    """Expose only an explicit ThinkingPart through AG-UI and saved history."""
    client, _app = ai_client
    headers = await authenticated_project(client, tmp_path)
    await configure_provider(client, headers)
    project_id = (await client.get("/api/project")).json()["id"]
    _use_stream_model(monkeypatch, _reasoning_stream)

    response = await client.post(
        "/api/agent/run",
        json=agent_body(
            threadId="ask-reasoning-thread",
            runId="ask-reasoning-run",
            forwardedProps=_forwarded("ask", project_id),
        ),
        headers={**headers, "Accept": "text/event-stream"},
    )

    assert response.status_code == 200, response.text
    events = _stream_events(response)
    assert any(event["type"] == "REASONING_MESSAGE_START" for event in events)
    assert any(
        event["type"] == "REASONING_MESSAGE_CONTENT"
        and event["delta"] == "Checked the current workspace context."
        for event in events
    )
    assert any(event["type"] == "REASONING_MESSAGE_END" for event in events)
    thread = await client.get("/api/agent/threads/ask-reasoning-thread?mode=ask")
    assert thread.status_code == 200
    assert thread.json()["conversation"] == [
        {"role": "user", "text": "Help"},
        {"role": "reasoning", "text": "Checked the current workspace context."},
        {"role": "assistant", "text": "The analysis is ready."},
    ]


async def test_assistant_markdown_is_formatted_without_allowing_raw_html(
    ai_client, tmp_path
):
    """Render common assistant Markdown while keeping model HTML inert."""
    client, _app = ai_client
    headers = await authenticated_project(client, tmp_path)
    response = await client.post(
        "/api/agent/markdown",
        json={
            "text": (
                "**Workspace**: `PXD`\n\n- one\n- two\n\n"
                "| A | B |\n| --- | --- |\n| 1 | 2 |\n\n"
                "<script>alert(1)</script>\n\n"
                "[unsafe](javascript:alert(1))"
            )
        },
        headers=headers,
    )

    assert response.status_code == 200
    html = response.json()["html"]
    assert "<strong>Workspace</strong>" in html
    assert "<code>PXD</code>" in html
    assert "<ul>" in html
    assert "<table>" in html
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in html
    assert "<script>" not in html
    assert 'href="javascript:' not in html
