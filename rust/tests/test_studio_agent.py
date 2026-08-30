"""Offline verification of Mokume Studio's fixed AI tool surfaces."""

from __future__ import annotations

from unittest.mock import Mock

import pytest
from pydantic_ai import DeferredToolRequests
from pydantic_ai.messages import ModelResponse, TextPart, ToolCallPart
from pydantic_ai.models.function import FunctionModel

from mokume.studio.agent import ASK_AGENT, PLAN_AGENT, StudioAgentDeps
from mokume.studio.models import ProjectRecord


pytestmark = pytest.mark.anyio


def agent_deps() -> StudioAgentDeps:
    """Build server-owned dependencies without loading data or a provider."""
    controller = Mock()
    controller.start_approved.side_effect = AssertionError(
        "approval-gated compute must not execute"
    )
    return StudioAgentDeps(
        controller=controller,
        project=ProjectRecord(id="project", root="/tmp/project", opened_at="now"),
        provider="offline",
        model="function-model",
    )


async def test_ask_agent_exposes_only_read_only_context():
    """Ask mode cannot see planning or compute tools."""
    exposed: list[str] = []

    async def answer(_messages, info):
        exposed.extend(tool.name for tool in info.function_tools)
        return ModelResponse(parts=[TextPart(content="read-only")])

    result = await ASK_AGENT.run(
        "Explain this dataset",
        model=FunctionModel(answer),
        deps=agent_deps(),
    )

    assert result.output == "read-only"
    assert exposed == ["get_analysis_context"]


async def test_plan_agent_exposes_fixed_tools_and_defers_compute_approval():
    """Plan mode exposes no arbitrary executor and defers its compute tool."""
    exposed: list[str] = []

    async def request_compute(_messages, info):
        exposed.extend(tool.name for tool in info.function_tools)
        return ModelResponse(
            parts=[
                ToolCallPart(
                    tool_name="run_approved_evaluation",
                    args={"approval_id": "approval", "payload_hash": "hash"},
                    tool_call_id="call-1",
                )
            ]
        )

    result = await PLAN_AGENT.run(
        "Run the approved plan",
        model=FunctionModel(request_compute),
        deps=agent_deps(),
    )

    assert set(exposed) == {
        "get_analysis_context",
        "prepare_evaluation",
        "run_approved_evaluation",
    }
    assert isinstance(result.output, DeferredToolRequests)
    assert result.output.calls == []
    assert [call.tool_name for call in result.output.approvals] == [
        "run_approved_evaluation"
    ]
