"""Shared construction helpers for optional Studio HTTP tests."""

from pathlib import Path
from unittest.mock import patch

from fastapi import FastAPI

from mokume.studio.app import create_app


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
