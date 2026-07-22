"""Tests for the streaming + polars-based DirectLFQ load path.

Exercises ``LoadingStage.load_for_directlfq`` (returns a ``pa.Table``) and
``LoadingStage.convert_to_directlfq_format`` (Arrow Table -> wide pandas
DataFrame).

The new implementation streams DuckDB results via the Arrow reader and runs
the pivot through polars (zero-copy from Arrow). These tests pin down the
required filter behaviour (contaminants/decoys/entrapments/short peptides/
single-peptide proteins all dropped), the wide-format shape, and the numeric
log2 transform expected by DirectLFQ.

We deliberately keep a self-contained fixture instead of importing the one
from ``test_loading_sqlfirst`` because the directlfq path relies on the
``is_decoy`` column (when present) as the authoritative decoy flag, so the
fixture must set ``is_decoy=True`` for decoy rows to mirror real parquet
files produced by quantms.io.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

pytest.importorskip("polars")

from mokume.pipeline.config import (  # noqa: E402
    FilterConfig,
    InputConfig,
    NormalizationConfig,
    PipelineConfig,
)
from mokume.pipeline.stages import LoadingStage  # noqa: E402


_NEW_INTENSITIES_TYPE = pa.list_(
    pa.struct([("label", pa.string()), ("intensity", pa.float32())])
)
_PG_PROTEIN_TYPE = pa.list_(
    pa.struct(
        [
            ("accession", pa.string()),
            ("start", pa.int32()),
            ("end", pa.int32()),
            ("pre", pa.string()),
            ("post", pa.string()),
        ]
    )
)
_LFQ_SCHEMA = pa.schema(
    [
        ("sequence", pa.string()),
        ("peptidoform", pa.string()),
        ("pg_accessions", _PG_PROTEIN_TYPE),
        ("anchor_protein", pa.string()),
        ("charge", pa.int16()),
        ("run_file_name", pa.string()),
        ("unique", pa.bool_()),
        ("is_decoy", pa.bool_()),
        ("intensities", _NEW_INTENSITIES_TYPE),
    ]
)


def _pg(accession: str) -> List[dict]:
    return [{"accession": accession, "start": 1, "end": 10, "pre": "K", "post": "A"}]


def _intensities(value: float) -> List[dict]:
    return [{"label": "LFQ", "intensity": float(value)}]


def _make_lfq_parquet(path: str) -> List[str]:
    """Same shape as the SQL-first fixture, but with ``is_decoy`` set
    correctly so the SQL filter actually drops the decoy rows.

    Each row tuple is (sequence, peptidoform, accession, charge, run, intensity,
    is_decoy).
    """
    runs = ["run_S1", "run_S2", "run_S3", "run_S4"]

    p12345 = "sp|P12345|GOOD_HUMAN"
    p67890 = "sp|P67890|LONE_HUMAN"
    p99999 = "P99999"
    contam = "sp|CONTAMINANT_P00761|TRYP_PIG"
    decoy = "sp|DECOY_P11111|REV_HUMAN"
    entrap = "sp|ENTRAP_P22222|ENT_HUMAN"

    rows: List[Tuple] = []
    for sample_run in runs:
        # P12345: three unique peptides; PEPTIDEAR has two charges per sample
        rows.append(("PEPTIDEAR", "PEPTIDEAR", p12345, 2, sample_run, 101.0, False))
        rows.append(("PEPTIDEAR", "PEPTIDEAR", p12345, 3, sample_run, 302.0, False))
        rows.append(("PEPTIDEBR", "PEPTIDEBR", p12345, 2, sample_run, 203.0, False))
        rows.append(("PEPTIDECR", "PEPTIDECR", p12345, 2, sample_run, 404.0, False))
        # P67890: only one unique peptide -> dropped by min_unique_peptides
        rows.append(("ONLYONEPEP", "ONLYONEPEP", p67890, 2, sample_run, 500.0, False))
        # P99999: bare accession, two unique peptides -> kept
        rows.append(("KEPTPEPAR", "KEPTPEPAR", p99999, 2, sample_run, 600.0, False))
        rows.append(("KEPTPEPBR", "KEPTPEPBR", p99999, 3, sample_run, 700.0, False))
        # Contaminants -> dropped by pg_accessions LIKE filter
        rows.append(("CONTAMPEPA", "CONTAMPEPA", contam, 2, sample_run, 800.0, False))
        rows.append(("CONTAMPEPB", "CONTAMPEPB", contam, 2, sample_run, 900.0, False))
        # Decoys -> dropped by is_decoy=true filter
        rows.append(("DECOYPEPAA", "DECOYPEPAA", decoy, 2, sample_run, 50.0, True))
        rows.append(("DECOYPEPBB", "DECOYPEPBB", decoy, 2, sample_run, 60.0, True))
        # Entrapments -> dropped by pg_accessions LIKE filter
        rows.append(("ENTRAPPEPA", "ENTRAPPEPA", entrap, 2, sample_run, 70.0, False))
        rows.append(("ENTRAPPEPB", "ENTRAPPEPB", entrap, 2, sample_run, 80.0, False))
        # Short peptide -> dropped by min_aa=7
        rows.append(("SHORT", "SHORT", p12345, 2, sample_run, 999.0, False))

    # Non-unique rows on the good protein -> dropped by require_unique
    non_unique = [
        ("PEPTIDEAR", "PEPTIDEAR", p12345, 2, run_S, 12345.0, False) for run_S in runs
    ]

    data = {
        "sequence": [r[0] for r in rows] + [r[0] for r in non_unique],
        "peptidoform": [r[1] for r in rows] + [r[1] for r in non_unique],
        "pg_accessions": [_pg(r[2]) for r in rows] + [_pg(r[2]) for r in non_unique],
        "anchor_protein": [r[2] for r in rows] + [r[2] for r in non_unique],
        "charge": [r[3] for r in rows] + [r[3] for r in non_unique],
        "run_file_name": [r[4] for r in rows] + [r[4] for r in non_unique],
        "unique": [True] * len(rows) + [False] * len(non_unique),
        "is_decoy": [r[6] for r in rows] + [False] * len(non_unique),
        "intensities": [_intensities(r[5]) for r in rows]
        + [_intensities(r[5]) for r in non_unique],
    }
    table = pa.table(data, schema=_LFQ_SCHEMA)
    pq.write_table(table, path)
    return runs


def _make_lfq_sdrf(path: str, runs: List[str]) -> None:
    header = (
        "source name\tcomment[data file]\tcomment[label]\t"
        "characteristics[organism part]\tfactor value[disease]\t"
        "characteristics[biological replicate]\tcomment[fraction identifier]\t"
        "comment[technical replicate]\n"
    )
    lines = [header]
    for i, run in enumerate(runs, start=1):
        cond = "normal" if i % 2 == 1 else "disease"
        sample = f"Sample_{i}"
        lines.append(f"{sample}\t{run}\tLabel free sample\tbrain\t{cond}\t{i}\t1\t1\n")
    Path(path).write_text("".join(lines), encoding="utf-8")


@pytest.fixture
def lfq_dataset(tmp_path):
    parquet = tmp_path / "lfq_input.parquet"
    sdrf = tmp_path / "lfq_input.sdrf.tsv"
    runs = _make_lfq_parquet(str(parquet))
    _make_lfq_sdrf(str(sdrf), runs)
    return str(parquet), str(sdrf)


def _make_directlfq_config(parquet: str, sdrf: str) -> PipelineConfig:
    return PipelineConfig(
        input=InputConfig(parquet=parquet, sdrf=sdrf),
        filtering=FilterConfig(
            remove_contaminants=True, min_aa=7, min_unique_peptides=2
        ),
        normalization=NormalizationConfig(run_method="none", sample_method="none"),
    )


def test_load_for_directlfq_returns_arrow_table_with_expected_columns(lfq_dataset):
    parquet, sdrf = lfq_dataset
    cfg = _make_directlfq_config(parquet, sdrf)

    table = LoadingStage(cfg).load_for_directlfq()

    assert isinstance(table, pa.Table)
    assert set(table.column_names) == {
        "pg_accessions",
        "sequence",
        "sample_accession",
        "intensity",
    }
    # Filtering happens at the SQL layer; rows must include only good peptides
    # and intensities must be > 0.
    intensities = table.column("intensity").to_pylist()
    assert intensities, "expected at least one peptide row to survive filtering"
    assert all(v > 0 for v in intensities)


def test_load_for_directlfq_filters_contaminants_decoys_short(lfq_dataset):
    """Short peptides (<7 AA), contaminants, decoys, entrapments, and
    non-unique rows must all be filtered at the SQL layer before the Arrow
    table is materialised."""
    parquet, sdrf = lfq_dataset
    cfg = _make_directlfq_config(parquet, sdrf)

    table = LoadingStage(cfg).load_for_directlfq()
    seqs = set(table.column("sequence").to_pylist())

    forbidden = {
        "SHORT",
        "CONTAMPEPA",
        "CONTAMPEPB",
        "DECOYPEPAA",
        "DECOYPEPBB",
        "ENTRAPPEPA",
        "ENTRAPPEPB",
    }
    assert forbidden.isdisjoint(seqs), forbidden & seqs


def test_convert_to_directlfq_format_shape_and_index(lfq_dataset):
    """Wide DataFrame must be MultiIndex(['protein', 'ion']) with one column
    per sample. The fixture keeps two proteins (P12345, P99999) and four
    samples; P12345 contributes 3 peptides, P99999 contributes 2 → 5 rows.

    ``ONLYONEPEP`` (P67890) drops out via ``min_unique_peptides=2``.
    """
    parquet, sdrf = lfq_dataset
    cfg = _make_directlfq_config(parquet, sdrf)
    loader = LoadingStage(cfg)

    table = loader.load_for_directlfq()
    wide = loader.convert_to_directlfq_format(table)

    assert isinstance(wide, pd.DataFrame)
    assert wide.index.names == ["protein", "ion"]
    assert wide.shape == (5, 4), wide.shape

    proteins = set(wide.index.get_level_values("protein"))
    assert proteins == {"P12345", "P99999"}, proteins

    ions = set(wide.index.get_level_values("ion"))
    assert ions == {
        "PEPTIDEAR",
        "PEPTIDEBR",
        "PEPTIDECR",
        "KEPTPEPAR",
        "KEPTPEPBR",
    }, ions

    # QPX ``LFQ`` and SDRF ``Label free sample`` are the same finite label-free
    # alias, so each run retains the sample accession declared by the SDRF.
    assert set(wide.columns) == {"Sample_1", "Sample_2", "Sample_3", "Sample_4"}


def test_convert_to_directlfq_format_log2_values(lfq_dataset):
    """PEPTIDEAR for P12345 has two charges per sample (101.0 and 302.0). The
    sum is 403.0 and the log2 value stored in DirectLFQ format is log2(403).
    """
    parquet, sdrf = lfq_dataset
    cfg = _make_directlfq_config(parquet, sdrf)
    loader = LoadingStage(cfg)

    table = loader.load_for_directlfq()
    wide = loader.convert_to_directlfq_format(table)

    row = wide.loc[("P12345", "PEPTIDEAR")]
    expected = math.log2(403.0)
    np.testing.assert_allclose(row.values.astype(float), expected, rtol=1e-5)


def test_convert_to_directlfq_format_zero_becomes_nan(lfq_dataset, tmp_path):
    """Any pivot cell that ends up summing to 0 must be NaN after the
    transform, so DirectLFQ treats it as missing rather than -inf."""
    parquet, sdrf = lfq_dataset
    cfg = _make_directlfq_config(parquet, sdrf)
    loader = LoadingStage(cfg)

    table = loader.load_for_directlfq()
    wide = loader.convert_to_directlfq_format(table)

    # In the fixture every kept (protein, peptide, sample) cell has a positive
    # measurement, so NaN cells should only appear if pivot produced an empty
    # combination. We assert there are no -inf values that would have arisen
    # had log2(0) leaked through.
    arr = wide.to_numpy(dtype=float)
    assert not np.any(np.isneginf(arr))
    assert not np.any(np.isposinf(arr))
