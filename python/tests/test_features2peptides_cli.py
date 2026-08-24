"""CLI tests for features2peptides."""

import json
from pathlib import Path

from click.testing import CliRunner

from mokume.mokume_cli import cli


def _patch_peptide_normalization(monkeypatch):
    """Patch peptide_normalization and return captured keyword arguments."""
    captured = {}

    def fake_peptide_normalization(**kwargs):
        """Capture CLI arguments passed to peptide_normalization."""
        captured.update(kwargs)

    monkeypatch.setattr(
        "mokume.commands.features2peptides.peptide_normalization",
        fake_peptide_normalization,
    )
    return captured


def _invoke_features2peptides(tmp_path: Path, parquet: Path):
    """Run features2peptides with --keep-shared-peptides enabled."""
    return CliRunner().invoke(
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


def test_features2peptides_passes_keep_shared_peptides(monkeypatch, tmp_path: Path):
    """The CLI forwards --keep-shared-peptides to peptide normalization."""
    parquet = tmp_path / "input.parquet"
    parquet.write_text("placeholder", encoding="utf-8")
    captured = _patch_peptide_normalization(monkeypatch)

    result = _invoke_features2peptides(tmp_path, parquet)

    assert result.exit_code == 0
    assert captured["keep_shared_peptides"] is True


def test_features2peptides_rejects_unscoped_irs_option(tmp_path: Path):
    parquet = tmp_path / "input.parquet"
    parquet.write_text("placeholder", encoding="utf-8")

    result = CliRunner().invoke(
        cli,
        [
            "features2peptides",
            "-p",
            str(parquet),
            "-o",
            str(tmp_path / "peptides.csv"),
            "--irs_scope",
            "two_stage",
        ],
    )

    assert result.exit_code != 0
    assert "require --irs_channel or --irs_autodetect_regex" in result.output


def test_filter_config_controls_duplicate_length_and_unique_thresholds(
    monkeypatch,
    tmp_path: Path,
):
    """Filter config values replace defaults shared with legacy CLI options."""
    parquet = tmp_path / "input.parquet"
    parquet.write_text("placeholder", encoding="utf-8")
    config = tmp_path / "filters.json"
    config.write_text(
        json.dumps(
            {
                "peptide": {"min_peptide_length": 11},
                "protein": {"min_unique_peptides": 4},
            }
        ),
        encoding="utf-8",
    )
    captured = _patch_peptide_normalization(monkeypatch)

    result = CliRunner().invoke(
        cli,
        [
            "features2peptides",
            "-p",
            str(parquet),
            "-o",
            str(tmp_path / "peptides.csv"),
            "--filter-config",
            str(config),
        ],
    )

    assert result.exit_code == 0
    assert captured["min_aa"] == 11
    assert captured["min_unique"] == 4


def test_explicit_legacy_thresholds_override_filter_config(monkeypatch, tmp_path: Path):
    """Explicit legacy CLI thresholds remain authoritative over config values."""
    parquet = tmp_path / "input.parquet"
    parquet.write_text("placeholder", encoding="utf-8")
    config = tmp_path / "filters.json"
    config.write_text(
        json.dumps(
            {
                "peptide": {"min_peptide_length": 11},
                "protein": {"min_unique_peptides": 4},
            }
        ),
        encoding="utf-8",
    )
    captured = _patch_peptide_normalization(monkeypatch)

    result = CliRunner().invoke(
        cli,
        [
            "features2peptides",
            "-p",
            str(parquet),
            "-o",
            str(tmp_path / "peptides.csv"),
            "--filter-config",
            str(config),
            "--min_aa",
            "9",
            "--min_unique",
            "3",
        ],
    )

    assert result.exit_code == 0
    assert captured["min_aa"] == 9
    assert captured["min_unique"] == 3


def test_duplicate_min_unique_cli_options_are_rejected(tmp_path: Path):
    """Two CLI spellings must not silently override each other."""
    parquet = tmp_path / "input.parquet"
    parquet.write_text("placeholder", encoding="utf-8")

    result = CliRunner().invoke(
        cli,
        [
            "features2peptides",
            "-p",
            str(parquet),
            "-o",
            str(tmp_path / "peptides.csv"),
            "--min_unique",
            "3",
            "--filter-min-unique-peptides",
            "4",
        ],
    )

    assert result.exit_code != 0
    assert "Choose either --min_unique or --filter-min-unique-peptides" in result.output
