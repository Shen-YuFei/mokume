"""Shared construction helpers for optional Studio HTTP tests."""

from pathlib import Path

from fastapi import FastAPI

from mokume.studio.app import create_app


def make_studio_app(port: int, token: str, state_directory: Path) -> FastAPI:
    """Create a testing-mode Studio app with one isolated state directory."""
    return create_app(
        port=port,
        startup_token=token,
        state_directory=state_directory,
        testing=True,
    )
