"""CLI and packaging contract tests for the optional Studio runtime."""

import importlib
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from mokume.studio import cli

try:
    import tomllib
except ModuleNotFoundError:
    tomllib = importlib.import_module("tomli")


def test_root_help_advertises_studio(monkeypatch, capsys):
    """Root help exposes the optional Studio command."""
    entrypoint = importlib.import_module("mokume.__main__")
    monkeypatch.setattr(entrypoint.sys, "argv", ["mokume", "--help"])

    with pytest.raises(SystemExit) as exc_info:
        entrypoint.main()

    assert exc_info.value.code == 0
    assert "studio" in capsys.readouterr().out


def test_studio_command_is_dispatched_without_rust(monkeypatch):
    """The Python entry point routes Studio without loading the Rust parser."""
    entrypoint = importlib.import_module("mokume.__main__")
    observed = {}
    module = SimpleNamespace(main=lambda argv: observed.update(argv=argv) or 0)
    real_import = entrypoint.importlib.import_module

    def import_module(name):
        if name == "mokume.commands.studio":
            return module
        return real_import(name)

    monkeypatch.setattr(entrypoint.importlib, "import_module", import_module)
    monkeypatch.setattr(
        entrypoint.sys,
        "argv",
        ["mokume", "studio", "--port", "9000", "--no-browser"],
    )

    with pytest.raises(SystemExit) as exc_info:
        entrypoint.main()

    assert exc_info.value.code == 0
    assert observed["argv"] == ["--port", "9000", "--no-browser"]


@pytest.mark.parametrize("value", ["0", "65536", "not-a-port"])
def test_invalid_port_is_rejected(value):
    """Ports outside the TCP range and non-integers are rejected."""
    with pytest.raises(SystemExit):
        cli.build_parser().parse_args(["--port", value])


def test_help_does_not_load_optional_runtime(monkeypatch):
    """Studio help remains available without optional web dependencies."""
    monkeypatch.setattr(
        cli,
        "_load_runtime",
        lambda: pytest.fail("Studio runtime loaded while rendering help"),
    )

    with pytest.raises(SystemExit) as exc_info:
        cli.main(["--help"])

    assert exc_info.value.code == 0


def test_missing_optional_runtime_has_install_hint(monkeypatch):
    """Missing optional dependencies produce an actionable install hint."""
    monkeypatch.setitem(sys.modules, "uvicorn", None)

    with pytest.raises(RuntimeError, match=r"mokume\[studio\]"):
        getattr(cli, "_load_runtime")()


def test_explicit_port_does_not_fall_back(monkeypatch):
    """An explicitly occupied port fails instead of silently changing ports."""
    attempted = []

    def occupied(port):
        attempted.append(port)
        raise OSError("occupied")

    monkeypatch.setattr(cli, "_bind_port", occupied)
    with pytest.raises(RuntimeError, match=rf"{cli.HOST}:9000"):
        getattr(cli, "_select_socket")(9000)

    assert attempted == [9000]


def test_default_port_scan_is_bounded(monkeypatch):
    """Automatic port selection only scans the documented bounded range."""
    attempted = []

    def occupied(port):
        attempted.append(port)
        raise OSError("occupied")

    monkeypatch.setattr(cli, "_bind_port", occupied)
    with pytest.raises(RuntimeError, match="no free Studio port"):
        getattr(cli, "_select_socket")(None)

    assert attempted == list(range(cli.DEFAULT_PORT, cli.DEFAULT_PORT + 50))


def test_clean_keyboard_interrupt_does_not_escape_server(monkeypatch):
    """A graceful Uvicorn interrupt does not leak a traceback to the shell."""

    def run(*, sockets):
        assert sockets == ["listener"]
        raise KeyboardInterrupt

    def server_factory(_config):
        return SimpleNamespace(started=True, should_exit=False, run=run)

    uvicorn = SimpleNamespace(
        Config=lambda *args, **kwargs: (args, kwargs),
        Server=server_factory,
    )
    monkeypatch.setattr(cli.secrets, "token_urlsafe", lambda _size: "token")

    getattr(cli, "_serve")(
        uvicorn,
        lambda **_kwargs: "app",
        "listener",
        9000,
        open_browser=False,
    )


def test_public_extras_match_runtime_capabilities():
    """The wheel exposes only supported feature extras and a complete union."""
    metadata = tomllib.loads(
        (Path(__file__).parents[1] / "pyproject.toml").read_text(encoding="utf-8")
    )
    extras = metadata["project"]["optional-dependencies"]
    assert set(extras) == {"analysis", "tissuemap", "plugin", "studio", "all"}
    assert set(extras["analysis"]) == {
        "numpy",
        "pandas",
        "scipy",
        "matplotlib",
        "seaborn",
        "scikit-learn",
        "plotly",
    }
    expected_studio = {
        "numpy",
        "pandas",
        "PyYAML>=6.0",
        "scipy",
        "scikit-learn",
        "statsmodels>=0.13",
        "fastapi",
        "uvicorn",
        "jinja2",
        "platformdirs",
        "pydantic-ai-slim[ag-ui,openai,anthropic,google]>=2.36,<3",
        "ag-ui-protocol>=0.1.19,<1",
    }

    assert expected_studio <= set(extras["studio"])
    assert set(extras["all"]) == set().union(
        *(set(extras[name]) for name in ("analysis", "tissuemap", "plugin", "studio"))
    )
