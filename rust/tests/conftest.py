"""Shared pytest configuration for the Rust-backed wheel."""

import os
import httpx
import pytest

from studio_test_support import AI_ORIGIN, AI_TOKEN, make_studio_app

# Headless matplotlib backend so the periphery plotting tests run without a
# display (must be set before matplotlib is first imported).
os.environ.setdefault("MPLBACKEND", "Agg")


@pytest.fixture
def anyio_backend():
    """Keep asynchronous HTTP and agent tests on the asyncio backend."""
    return "asyncio"


@pytest.fixture(name="ai_client")
async def build_ai_client(tmp_path):
    """Create one authenticated-capable Studio app and ASGI client."""
    app = make_studio_app(AI_TOKEN, tmp_path / "state")
    transport = httpx.ASGITransport(app=app)
    async with httpx.AsyncClient(transport=transport, base_url=AI_ORIGIN) as client:
        yield client, app
