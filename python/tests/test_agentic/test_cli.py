"""CLI dependency-boundary tests for agentic analysis."""

# The loader is private because it exists only to enforce this command's import boundary.
# pylint: disable=protected-access

import builtins
import os
from pathlib import Path
import subprocess
import sys

import click
import pytest

from mokume.commands import agentic


_PYTHON_ROOT = Path(__file__).resolve().parents[2]
_CLI_WITHOUT_YAML = """
import builtins
import sys

original_import = builtins.__import__


def import_without_yaml(name, *args, **kwargs):
    if name == "yaml" or name.startswith("yaml."):
        raise ModuleNotFoundError("No module named 'yaml'", name="yaml")
    return original_import(name, *args, **kwargs)


builtins.__import__ = import_without_yaml
from mokume.mokume_cli import main

main()
"""


def _run_cli_without_yaml(
    *args: str, cwd: Path | None = None
) -> subprocess.CompletedProcess[str]:
    """Run one CLI process whose import system cannot load PyYAML."""
    return subprocess.run(
        [sys.executable, "-c", _CLI_WITHOUT_YAML, *args],
        check=False,
        capture_output=True,
        text=True,
        cwd=cwd,
        env={**os.environ, "PYTHONPATH": str(_PYTHON_ROOT)},
    )


@pytest.mark.parametrize(
    ("args", "expected"),
    [
        (("--help",), "agentic"),
        (("--version",), "mokume 0.1.0"),
        (("features2proteins", "--help"), "--quant-method"),
        (("agentic", "--help"), "optimize"),
        (("agentic", "optimize", "--help"), "--protein-matrix"),
    ],
)
def test_cli_help_does_not_require_agentic_extra(
    args: tuple[str, ...], expected: str
) -> None:
    """CLI discovery and help remain available without PyYAML."""
    result = _run_cli_without_yaml(*args)

    assert result.returncode == 0, result.stderr
    assert expected in result.stdout
    assert "Traceback" not in result.stderr


def test_agentic_command_reports_missing_extra_before_reading_inputs(tmp_path) -> None:
    """Selecting agentic gives an installation hint instead of a traceback."""
    protein_matrix = tmp_path / "proteins.csv"
    sdrf = tmp_path / "input.sdrf.tsv"
    output_dir = tmp_path / "optimization"
    protein_matrix.write_text("not,parsed\n", encoding="utf-8")
    sdrf.write_text("not\tparsed\n", encoding="utf-8")

    result = _run_cli_without_yaml(
        "agentic",
        "optimize",
        "--protein-matrix",
        str(protein_matrix),
        "--sdrf",
        str(sdrf),
        "--contrasts",
        "A vs B",
        "--no-llm",
        "--save-env",
        "--llm-base-url",
        "https://example.invalid/v1",
        "--output-dir",
        str(output_dir),
        cwd=tmp_path,
    )

    output = result.stdout + result.stderr
    assert result.returncode == 1
    assert "pip install mokume[agentic]" in output
    assert "Traceback" not in output
    assert "ParserError" not in output
    assert not (tmp_path / ".env").exists()
    assert not output_dir.exists()


@pytest.mark.parametrize("missing_name", ["mokume.internal", "optional_backend"])
def test_optimizer_loader_does_not_hide_other_import_failures(
    monkeypatch, missing_name: str
) -> None:
    """Only the dependency owned by the agentic extra gets rewritten."""
    original_import = builtins.__import__

    def fail_optimizer_import(name, *args, **kwargs):
        if name == "mokume.agentic.optimizer":
            raise ModuleNotFoundError("dependency import failed", name=missing_name)
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fail_optimizer_import)

    with pytest.raises(ModuleNotFoundError) as error:
        agentic._load_optimizer()

    assert error.value.name == missing_name


def test_optimizer_loader_explains_missing_yaml(monkeypatch) -> None:
    """A missing PyYAML dependency points at the complete agentic extra."""
    original_import = builtins.__import__

    def fail_optimizer_import(name, *args, **kwargs):
        if name == "mokume.agentic.optimizer":
            raise ModuleNotFoundError("No module named 'yaml'", name="yaml")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fail_optimizer_import)

    with pytest.raises(click.ClickException, match=r"mokume\[agentic\]"):
        agentic._load_optimizer()
