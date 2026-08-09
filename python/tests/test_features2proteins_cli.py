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


def test_features2proteins_passes_directlfq_and_batch_options(monkeypatch, tmp_path):
    parquet, sdrf = _make_input_files(tmp_path)
    captured = {}

    def fake_run_pipeline(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(pipeline, "features_to_proteins", fake_run_pipeline)

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
            "directlfq",
            "--directlfq-min-nonan",
            "3",
            "--batch-correction",
            "--batch-method",
            "column",
            "--batch-column",
            "characteristics[batch]",
            "--batch-covariates",
            "characteristics[sex],characteristics[organism part]",
            "--batch-nonparametric",
            "--batch-mean-only",
            "--batch-ref",
            "2",
            "--de-fdr-method",
            "ihw",
        ],
    )

    assert result.exit_code == 0
    assert captured["directlfq_min_nonan"] == 3
    assert captured["batch_correction"] is True
    assert captured["batch_method"] == "column"
    assert captured["batch_column"] == "characteristics[batch]"
    assert captured["batch_covariates"] == [
        "characteristics[sex]",
        "characteristics[organism part]",
    ]
    assert captured["batch_parametric"] is False
    assert captured["batch_mean_only"] is True
    assert captured["batch_ref"] == 2
    assert captured["de_fdr_method"] == "ihw"


def test_features2proteins_requires_batch_column_for_column_method(tmp_path):
    parquet, sdrf = _make_input_files(tmp_path)

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
            "--batch-correction",
            "--batch-method",
            "column",
        ],
    )

    assert result.exit_code != 0
    assert "requires --batch-column option" in result.output


def test_features2proteins_preserves_empty_ensemble_members(monkeypatch, tmp_path):
    parquet, _ = _make_input_files(tmp_path)
    captured = {}

    def fake_run_pipeline(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(pipeline, "features_to_proteins", fake_run_pipeline)

    result = CliRunner().invoke(
        cli,
        [
            "features2proteins",
            "-p",
            parquet,
            "-o",
            "out.csv",
            "--de-ensemble-methods",
            "",
        ],
    )

    assert result.exit_code == 0
    assert captured["de_ensemble_methods"] == [""]


def test_features2proteins_accepts_msstats_with_sdrf(monkeypatch, tmp_path):
    """The CLI accepts MSstats only as the selected feature input."""
    _, sdrf = _make_input_files(tmp_path)
    msstats = tmp_path / "input_msstats_in.csv"
    msstats.write_text("ProteinName,PeptideSequence,Intensity\n", encoding="utf-8")
    captured = {}

    def fake_run_pipeline(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(pipeline, "features_to_proteins", fake_run_pipeline)
    result = CliRunner().invoke(
        cli,
        [
            "features2proteins",
            "--msstats",
            str(msstats),
            "--sdrf",
            sdrf,
            "--output",
            "out.csv",
        ],
    )

    assert result.exit_code == 0
    assert captured["parquet"] is None
    assert captured["msstats"] == str(msstats)
