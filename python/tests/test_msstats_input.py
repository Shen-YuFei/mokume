"""Regression coverage for native SDRF and MSstats input."""

from pathlib import Path

import pandas as pd
import pytest

from mokume.io.feature import Feature
from mokume.pipeline.config import InputConfig, PipelineConfig, QuantificationConfig
from mokume.pipeline.features_to_proteins import features_to_proteins
from mokume.pipeline.runner import run_pipeline


def _write_sdrf(path: Path, rows: list[dict]) -> None:
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_lfq_msstats_runs_end_to_end(tmp_path):
    """LFQ MSstats rows resolve through SDRF and reach protein output."""
    sdrf = tmp_path / "input.sdrf.tsv"
    msstats = tmp_path / "input_msstats_in.csv"
    output = tmp_path / "proteins.csv"
    _write_sdrf(
        sdrf,
        [
            {
                "source name": "S1",
                "comment[data file]": "run1.raw",
                "comment[label]": "AC=MS:1002038;NT=label free sample",
                "comment[technical replicate]": "1",
                "comment[fraction identifier]": "1",
                "characteristics[biological replicate]": "1",
                "factor value[condition]": "A",
            },
            {
                "source name": "S2",
                "comment[data file]": "run2.raw",
                "comment[label]": "AC=MS:1002038;NT=label free sample",
                "comment[technical replicate]": "1",
                "comment[fraction identifier]": "1",
                "characteristics[biological replicate]": "1",
                "factor value[condition]": "B",
            },
        ],
    )
    pd.DataFrame(
        [
            {
                "ProteinName": "P1",
                "PeptideSequence": "PEPTIDEK",
                "PrecursorCharge": 2,
                "Run": 1,
                "Intensity": 100.0,
                "Reference": r"C:\data\run1.mzML",
            },
            {
                "ProteinName": "P1",
                "PeptideSequence": "PEPTIDEK",
                "PrecursorCharge": 2,
                "Run": 2,
                "Intensity": 200.0,
                "Reference": "run2.mzML",
            },
        ]
    ).to_csv(msstats, index=False)

    proteins = features_to_proteins(
        parquet=None,
        msstats=str(msstats),
        sdrf=str(sdrf),
        output=str(output),
        quant_method="sum",
        min_unique_peptides=1,
        run_normalization="none",
        sample_normalization="none",
    )

    row = proteins.set_index("ProteinName").loc["P1"]
    assert row["S1"] == 100.0
    assert row["S2"] == 200.0


def test_tmt10_channel_positions_use_sdrf_labels(tmp_path):
    """TMT channel positions map to canonical SDRF labels."""
    sdrf = tmp_path / "input.sdrf.tsv"
    msstats = tmp_path / "input_msstats_in.csv"
    labels = [
        "TMT126",
        "TMT127N",
        "TMT127C",
        "TMT128N",
        "TMT128C",
        "TMT129N",
        "TMT129C",
        "TMT130N",
        "TMT130C",
        "TMT131",
    ]
    _write_sdrf(
        sdrf,
        [
            {
                "source name": f"S{position}",
                "comment[data file]": "plex.raw",
                "comment[label]": label,
                "comment[technical replicate]": "1",
                "comment[fraction identifier]": "1",
                "characteristics[biological replicate]": str(position),
                "factor value[condition]": "A",
            }
            for position, label in enumerate(labels, start=1)
        ],
    )
    pd.DataFrame(
        [
            {
                "ProteinName": "P1",
                "PeptideSequence": ".(TMT6plex)PEPTIDEK",
                "Charge": 2,
                "Channel": channel,
                "Run": "1_1_1",
                "Intensity": intensity,
                "Reference": "plex.mzML_controllerType=0 controllerNumber=1 scan=1",
            }
            for channel, intensity in [(1, 100.0), (10, 1000.0)]
        ]
    ).to_csv(msstats, index=False)

    feature = Feature.from_msstats(str(msstats), str(sdrf))
    try:
        feature.enrich_with_sdrf(str(sdrf))
        rows = feature.parquet_db.execute(
            "SELECT channel, sample_accession, sequence, intensity "
            "FROM parquet_db ORDER BY channel"
        ).fetchall()
    finally:
        feature.parquet_db.close()

    assert rows == [
        ("TMT126", "S1", "PEPTIDEK", 100.0),
        ("TMT131", "S10", "PEPTIDEK", 1000.0),
    ]


def test_ratio_rejects_feature_level_msstats_input():
    """Ratio quantification rejects MSstats without PSM evidence."""
    config = PipelineConfig(
        input=InputConfig(msstats="input.csv", sdrf="input.sdrf.tsv"),
        quantification=QuantificationConfig(method="ratio"),
    )

    with pytest.raises(ValueError, match="PSM-level QPX input"):
        run_pipeline(config)
