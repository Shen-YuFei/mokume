"""Knowledge-location contracts shared by Studio and Plugin runtimes."""

from __future__ import annotations

import importlib
from pathlib import Path
import shutil
from types import SimpleNamespace

import pytest

from mokume.agentic.knowledge import load_knowledge_graph
from mokume.agentic.service import RecommendationService

REPOSITORY = Path(__file__).resolve().parents[2]
KNOWLEDGE_BUNDLE = (
    REPOSITORY / "rust" / "python" / "mokume" / "agentic" / "knowledge_bundle"
)
BUNDLED_KNOWLEDGE = KNOWLEDGE_BUNDLE / "knowledge.yaml"


def test_knowledge_resolution_prefers_explicit_then_environment_then_bundle(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Development overrides must not replace the installed default."""
    monkeypatch.delenv("MOKUME_AGENTIC_KNOWLEDGE", raising=False)
    bundled = load_knowledge_graph()
    assert bundled.fingerprint == load_knowledge_graph(BUNDLED_KNOWLEDGE).fingerprint
    assert RecommendationService() is not None

    copied = tmp_path / "knowledge"
    shutil.copytree(KNOWLEDGE_BUNDLE, copied)
    override = copied / "knowledge.yaml"
    override.write_text(
        override.read_text(encoding="utf-8").replace(
            "title: Evaluating differential expression methods for proteomics data",
            "title: Environment-selected knowledge snapshot",
            1,
        ),
        encoding="utf-8",
    )
    monkeypatch.setenv("MOKUME_AGENTIC_KNOWLEDGE", str(override))

    environment = load_knowledge_graph()
    assert environment.fingerprint != bundled.fingerprint
    assert load_knowledge_graph(BUNDLED_KNOWLEDGE).fingerprint == bundled.fingerprint


def test_mcp_main_uses_bundle_and_accepts_explicit_override(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The hidden CLI starts without plugin-root knowledge arguments."""
    module = importlib.import_module("mokume.agentic.mcp_server")
    observed: list[tuple[str, str | None]] = []

    def run(*, transport: str) -> None:
        observed.append(("transport", transport))

    def create_server(knowledge: str | None = None) -> SimpleNamespace:
        observed.append(("knowledge", knowledge))
        return SimpleNamespace(run=run)

    monkeypatch.setattr(module, "create_server", create_server)

    assert module.main([]) == 0
    assert observed == [("knowledge", None), ("transport", "stdio")]

    observed.clear()
    assert module.main(["--knowledge", str(BUNDLED_KNOWLEDGE)]) == 0
    assert observed == [
        ("knowledge", str(BUNDLED_KNOWLEDGE.resolve())),
        ("transport", "stdio"),
    ]
