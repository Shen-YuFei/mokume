"""Verify that the Python examples in the docs follow the actual code.

This guards against documentation drift: every ``python`` fenced block in the
listed files is parsed, every imported symbol is resolved against the installed
``mokume`` package, and every module-level attribute access (e.g.
``mokume.features2proteins``) is checked for existence. A renamed or removed
symbol, or a call to an API the pure-Python package does not expose (such as the
Rust-wheel-only ``mokume.features2proteins`` binding used by mistake in a
pure-Python example), fails this test instead of shipping a broken example.

It is a *static* check — it never executes the blocks — so illustrative snippets
that use placeholder variables or file paths (``peptides.csv``, ``protein_df``)
are fine. To extend coverage, add a file to ``_DOC_FILES``; only files whose
``python`` examples target the pure-Python ``mokume`` package belong here (the
Rust-wheel API pages document ``mokume._mokume`` and need the wheel installed).
"""

from __future__ import annotations

import ast
import importlib
import re
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]

_DOC_FILES = [
    _REPO_ROOT / "README.md",
    _REPO_ROOT / "docs" / "concepts" / "differential-expression.md",
]

_PYTHON_BLOCK = re.compile(r"^```python\s*$(.*?)^```\s*$", re.DOTALL | re.MULTILINE)


def _iter_doc_blocks():
    params = []
    for path in _DOC_FILES:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        for index, match in enumerate(_PYTHON_BLOCK.finditer(text)):
            params.append(pytest.param(match.group(1), id=f"{path.name}#{index}"))
    return params


def _resolve_from_import(module_name: str, symbol: str) -> None:
    module = importlib.import_module(module_name)
    if not hasattr(module, symbol):
        # ``from package import submodule`` is also valid.
        importlib.import_module(f"{module_name}.{symbol}")


@pytest.mark.parametrize("block", _iter_doc_blocks())
def test_doc_python_example_symbols_exist(block: str) -> None:
    """Resolve every import and module attribute used in one doc code block."""
    tree = ast.parse(block)  # SyntaxError => the example is malformed

    # 1. Every imported module / symbol must resolve. Record module aliases
    #    (``import pandas as pd`` -> ``pd``) for the attribute check below.
    module_aliases: dict[str, str] = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                importlib.import_module(alias.name)
                if alias.asname:
                    module_aliases[alias.asname] = alias.name
                elif "." not in alias.name:
                    module_aliases[alias.name] = alias.name
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            for alias in node.names:
                _resolve_from_import(node.module, alias.name)

    # 2. ``module.attr`` accesses on an imported module must exist on it — this
    #    is what catches a wheel-only API used in a pure-Python example.
    for node in ast.walk(tree):
        if (
            isinstance(node, ast.Attribute)
            and isinstance(node.value, ast.Name)
            and node.value.id in module_aliases
        ):
            module = importlib.import_module(module_aliases[node.value.id])
            assert hasattr(module, node.attr), (
                f"{module_aliases[node.value.id]}.{node.attr} does not exist; the "
                "documentation example references an API the package does not expose"
            )
