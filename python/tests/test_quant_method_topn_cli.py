"""``--quant-method top<N>`` is the only spelling of TopN on the CLI.

N used to travel in a companion option (``--topn`` / ``--topn_n``), which meant
one N had two spellings and, in the Rust kernel, ``top3 --topn 5`` silently ran
Top5 -- the digits in the method name were decorative. N now lives in the method
name alone. These tests pin that contract on the Python side; the Rust CLI is
held to the same one by ``mokume-cli``'s own suite.
"""

from pathlib import Path

import mokume.pipeline as pipeline
from click.testing import CliRunner

from mokume.mokume_cli import cli


def _make_input_files(tmp_path: Path) -> tuple[str, str]:
    parquet = tmp_path / "input.parquet"
    sdrf = tmp_path / "input.sdrf.tsv"
    parquet.write_text("placeholder", encoding="utf-8")
    sdrf.write_text("source name\nSample1\n", encoding="utf-8")
    return str(parquet), str(sdrf)


def _run(tmp_path, monkeypatch, method: str, extra: list[str] | None = None):
    """Invoke features2proteins with ``method``, capturing the pipeline kwargs."""
    parquet, sdrf = _make_input_files(tmp_path)
    captured: dict = {}

    def fake_run_pipeline(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(pipeline, "features_to_proteins", fake_run_pipeline)
    method_options = []
    if method.lower() == "pibaq":
        fasta = tmp_path / "proteome.fasta"
        fasta.write_text(">P1\nAAAAAAAK\n", encoding="utf-8")
        method_options = ["--fasta", str(fasta)]
    result = CliRunner().invoke(
        cli,
        [
            "features2proteins",
            "-p",
            parquet,
            "-o",
            "out.csv",
            "-s",
            sdrf,
            "--quant-method",
            method,
            *method_options,
            *(extra or []),
        ],
    )
    return result, captured


def test_top_n_carries_n_in_the_method_name(tmp_path, monkeypatch):
    for method, expected in [
        ("top1", "top1"),
        ("top3", "top3"),
        ("top5", "top5"),
        ("top10", "top10"),
        ("top100", "top100"),
    ]:
        result, captured = _run(tmp_path, monkeypatch, method)
        assert result.exit_code == 0, result.output
        assert captured["quant_method"] == expected


def test_bare_topn_normalizes_to_top3(tmp_path, monkeypatch):
    """``topn`` keeps the placeholder letter and means the canonical Top3."""
    result, captured = _run(tmp_path, monkeypatch, "topn")
    assert result.exit_code == 0, result.output
    assert captured["quant_method"] == "top3"


def test_method_name_is_case_insensitive(tmp_path, monkeypatch):
    result, captured = _run(tmp_path, monkeypatch, "TOP5")
    assert result.exit_code == 0, result.output
    assert captured["quant_method"] == "top5"


def test_top_without_a_numeral_is_rejected(tmp_path, monkeypatch):
    """``topa`` must fail rather than silently fall back to Top3."""
    for method in ["topa", "topx", "top", "topthree"]:
        result, _ = _run(tmp_path, monkeypatch, method)
        assert result.exit_code != 0, f"{method} should be rejected"
        assert "not a valid quantification method" in result.output


def test_top0_is_rejected(tmp_path, monkeypatch):
    """N = 0 parses as a numeral but averages nothing, so it is an error."""
    result, _ = _run(tmp_path, monkeypatch, "top0")
    assert result.exit_code != 0
    assert "must be an integer >= 1" in result.output


def test_topn_companion_option_is_gone(tmp_path, monkeypatch):
    """The old ``--topn`` flag must not be silently accepted and ignored."""
    result, _ = _run(tmp_path, monkeypatch, "top3", extra=["--topn", "5"])
    assert result.exit_code != 0
    assert "No such option" in result.output


def test_fixed_methods_still_work(tmp_path, monkeypatch):
    """The TopN param type must not disturb the non-TopN methods."""
    for method in [
        "pibaq",
        "maxlfq",
        "directlfq",
        "sum",
        "median",
        "spectral_count",
    ]:
        result, captured = _run(tmp_path, monkeypatch, method)
        assert result.exit_code == 0, result.output
        assert captured["quant_method"] == method


def test_removed_ibaq_method_name_is_rejected(tmp_path, monkeypatch):
    result, captured = _run(tmp_path, monkeypatch, "ibaq")
    assert result.exit_code != 0
    assert "not a valid quantification method" in result.output
    assert captured == {}


def test_unknown_method_lists_the_alternatives(tmp_path, monkeypatch):
    result, _ = _run(tmp_path, monkeypatch, "banana")
    assert result.exit_code != 0
    assert "not a valid quantification method" in result.output
    assert "top<N>" in result.output
