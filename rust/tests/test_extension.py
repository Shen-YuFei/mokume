"""Binding-level tests: the compiled extension loads and the CLI entry behaves."""

import pytest

import mokume
from mokume._mokume import run_cli


def test_version_is_nonempty():
    assert mokume.version()
    assert mokume.__version__ == mokume.version()


def test_compute_wrapper_rejects_bad_method(tmp_path):
    # A bad --method value is rejected by clap inside the extension and surfaced
    # as a normal RuntimeError -- no interpreter teardown.
    with pytest.raises(RuntimeError):
        mokume.peptides2protein(
            method="definitely-not-a-method",
            peptides="/nonexistent/peptides.csv",
            output=str(tmp_path / "out.tsv"),
        )


def test_run_cli_help_exits_zero():
    assert run_cli(["--help"]) == 0


def test_run_cli_subcommand_help_exits_zero():
    assert run_cli(["features2proteins", "--help"]) == 0


def test_run_cli_unknown_subcommand_is_nonzero():
    assert run_cli(["definitely-not-a-subcommand"]) != 0
