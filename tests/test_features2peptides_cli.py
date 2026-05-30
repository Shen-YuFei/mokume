"""CLI tests for features2peptides."""

from pathlib import Path

from click.testing import CliRunner

from mokume.mokume_cli import cli


def test_features2peptides_passes_keep_shared_peptides(monkeypatch, tmp_path: Path):
    """The CLI forwards --keep-shared-peptides to peptide normalization."""
    parquet = tmp_path / "input.parquet"
    parquet.write_text("placeholder", encoding="utf-8")
    captured = {}

    def fake_peptide_normalization(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(
        "mokume.commands.features2peptides.peptide_normalization",
        fake_peptide_normalization,
    )

    result = CliRunner().invoke(
        cli,
        [
            "features2peptides",
            "-p",
            str(parquet),
            "-o",
            str(tmp_path / "peptides.csv"),
            "--skip_normalization",
            "--keep-shared-peptides",
        ],
    )

    assert result.exit_code == 0
    assert captured["keep_shared_peptides"] is True
