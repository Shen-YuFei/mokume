"""Shared pytest configuration for the Rust-backed wheel."""

import os

import pytest

# Headless matplotlib backend so the periphery plotting tests run without a
# display (must be set before matplotlib is first imported).
os.environ.setdefault("MPLBACKEND", "Agg")


@pytest.fixture
def anyio_backend():
    """Keep asynchronous HTTP and agent tests on the asyncio backend."""
    return "asyncio"
