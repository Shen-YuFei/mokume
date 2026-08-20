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


def test_matrix_level_compute_api():
    matrix = [
        [2.0, 4.0, 8.0, 32.0, 64.0, 128.0],
        [8.0, 16.0, 32.0, 8.0, 16.0, 32.0],
        [32.0, None, 128.0, 4.0, 8.0, 16.0],
    ]
    normalized = mokume.normalize_matrix(
        matrix,
        "median",
        ["a1", "a2", "a3", "b1", "b2", "b3"],
        2,
    )
    assert len(normalized) == len(matrix)
    assert all(len(row) == 6 for row in normalized)

    imputed = mokume.impute_matrix(normalized, "mean", threads=2)
    assert imputed[2][1] is not None

    results = mokume.differential_expression(
        ["P1", "P2", "P3"],
        imputed,
        3,
        3,
        "limma",
        condition_a="case",
        condition_b="control",
        threads=2,
    )
    assert results
    assert {"ProteinName", "log2FC", "pvalue", "adj_pvalue", "significance"} <= set(
        results[0]
    )
