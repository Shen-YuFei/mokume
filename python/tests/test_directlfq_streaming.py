from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mokume.pipeline.directlfq_streaming import estimate_protein_intensities_streamed
from mokume.quantification.directlfq import DirectLFQQuantification
from mokume.quantification.maxlfq import MaxLFQQuantification

lfq_estimation = pytest.importorskip("directlfq.protein_intensity_estimation")
lfq_config = pytest.importorskip("directlfq.config")
lfq_manager = pytest.importorskip("directlfq.lfq_manager")


def _make_normed_directlfq_frame() -> pd.DataFrame:
    index = pd.MultiIndex.from_tuples(
        [
            ("P1", "PEPA"),
            ("P1", "PEPB"),
            ("P1", "PEPC"),
            ("P2", "PEPD"),
            ("P2", "PEPE"),
        ],
        names=["protein", "ion"],
    )
    intensities = np.array(
        [
            [100.0, 120.0, 110.0, 90.0],
            [200.0, 210.0, 205.0, np.nan],
            [50.0, 60.0, np.nan, 45.0],
            [400.0, 390.0, 410.0, 405.0],
            [80.0, np.nan, 75.0, 85.0],
        ],
        dtype=float,
    )
    return pd.DataFrame(
        np.log2(intensities),
        index=index,
        columns=["S1", "S2", "S3", "S4"],
    )


def _make_long_directlfq_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            ("P1", "PEPA", "S1", 100.0),
            ("P1", "PEPA", "S2", 120.0),
            ("P1", "PEPA", "S3", 110.0),
            ("P1", "PEPB", "S1", 200.0),
            ("P1", "PEPB", "S2", 210.0),
            ("P1", "PEPB", "S3", 205.0),
            ("P1", "PEPC", "S1", 50.0),
            ("P1", "PEPC", "S2", 60.0),
            ("P2", "PEPD", "S1", 400.0),
            ("P2", "PEPD", "S2", 390.0),
            ("P2", "PEPD", "S3", 410.0),
            ("P2", "PEPE", "S1", 80.0),
            ("P2", "PEPE", "S3", 75.0),
        ],
        columns=["ProteinName", "PeptideSequence", "SampleID", "Intensity"],
    )


def _sort_long_result(df: pd.DataFrame) -> pd.DataFrame:
    return df.sort_values(["ProteinName", "SampleID"]).reset_index(drop=True)


def test_maxlfq_directlfq_wrapper_matches_upstream_run_lfq():
    peptide_df = _make_long_directlfq_frame()
    maxlfq = MaxLFQQuantification(min_peptides=1, threads=2)
    assert maxlfq.using_directlfq

    streamed_df = maxlfq.quantify(
        peptide_df,
        protein_column="ProteinName",
        peptide_column="PeptideSequence",
        intensity_column="Intensity",
        sample_column="SampleID",
    )

    directlfq = DirectLFQQuantification(min_nonan=1, num_cores=1)
    with tempfile.TemporaryDirectory() as temp_dir:
        input_path = directlfq._prepare_input_file(
            peptide_df,
            protein_column="ProteinName",
            peptide_column="PeptideSequence",
            intensity_column="Intensity",
            sample_column="SampleID",
            temp_dir=temp_dir,
        )
        lfq_manager.run_lfq(
            input_file=input_path,
            min_nonan=1,
            num_cores=1,
            compile_normalized_ion_table=False,
        )
        output_paths = list(Path(temp_dir).glob("*protein_intensities*.tsv"))
        assert len(output_paths) == 1
        upstream_wide = pd.read_csv(output_paths[0], sep="\t")
        upstream_df = directlfq._parse_wide_output(
            upstream_wide,
            protein_column="ProteinName",
            sample_column="SampleID",
        )

    pd.testing.assert_frame_equal(
        _sort_long_result(streamed_df),
        _sort_long_result(upstream_df),
        check_exact=False,
        rtol=1e-10,
        atol=1e-10,
    )


def test_streaming_directlfq_matches_upstream_sequential_estimation():
    lfq_config.set_global_protein_and_ion_id(protein_id="protein", quant_id="ion")
    lfq_config.set_compile_normalized_ion_table(False)
    normed_df = _make_normed_directlfq_frame()

    upstream_df, upstream_ions = lfq_estimation.estimate_protein_intensities(
        normed_df.copy(deep=True),
        min_nonan=1,
        num_samples_quadratic=2,
        num_cores=1,
    )
    streamed_df = estimate_protein_intensities_streamed(
        normed_df.copy(deep=True),
        min_nonan=1,
        num_samples_quadratic=2,
    )

    assert upstream_ions is None
    pd.testing.assert_frame_equal(
        streamed_df.sort_values("protein").reset_index(drop=True),
        upstream_df.sort_values("protein").reset_index(drop=True),
        check_exact=False,
        rtol=1e-10,
        atol=1e-10,
    )


def test_streaming_directlfq_sorts_noncontiguous_protein_groups():
    lfq_config.set_global_protein_and_ion_id(protein_id="protein", quant_id="ion")
    lfq_config.set_compile_normalized_ion_table(False)
    normed_df = _make_normed_directlfq_frame().iloc[[0, 3, 1, 4, 2]]
    sorted_df = normed_df.sort_index()

    upstream_df, _ = lfq_estimation.estimate_protein_intensities(
        sorted_df.copy(deep=True),
        min_nonan=1,
        num_samples_quadratic=2,
        num_cores=1,
    )
    streamed_df = estimate_protein_intensities_streamed(
        normed_df.copy(deep=True),
        min_nonan=1,
        num_samples_quadratic=2,
    )

    pd.testing.assert_frame_equal(
        streamed_df.sort_values("protein").reset_index(drop=True),
        upstream_df.sort_values("protein").reset_index(drop=True),
        check_exact=False,
        rtol=1e-10,
        atol=1e-10,
    )


def test_resolve_directlfq_num_cores_mirrors_joblib_sentinels():
    import os

    from mokume.quantification.maxlfq import _resolve_directlfq_num_cores

    n = os.cpu_count() or 1
    assert _resolve_directlfq_num_cores(4) == 4
    assert _resolve_directlfq_num_cores(1) == 1
    assert _resolve_directlfq_num_cores(-1) == n  # all cores, like joblib n_jobs=-1
    assert _resolve_directlfq_num_cores(-2) == max(1, n - 1)  # all but one
    assert _resolve_directlfq_num_cores(0) is None


def test_maxlfq_threads_minus_one_reaches_directlfq_as_all_cores(monkeypatch):
    """threads=-1 ('all cores') must reach DirectLFQ as a parallel core count,
    not silently collapse to sequential (num_cores=None)."""
    import os

    import mokume.quantification.directlfq as directlfq_mod
    from mokume.quantification.maxlfq import (
        MaxLFQQuantification,
        _resolve_directlfq_num_cores,
    )

    captured = {}

    class FakeDirectLFQ:
        def __init__(self, min_nonan, num_cores=None):
            captured["num_cores"] = num_cores

        def quantify(self, peptide_df, **kwargs):
            return pd.DataFrame({"ProteinName": [], "SampleID": [], "Intensity": []})

    maxlfq = MaxLFQQuantification(min_peptides=1, threads=-1)
    if not maxlfq.using_directlfq:
        pytest.skip("directlfq not available")

    monkeypatch.setattr(directlfq_mod, "DirectLFQQuantification", FakeDirectLFQ)
    maxlfq._quantify_with_directlfq(
        pd.DataFrame(
            {
                "ProteinName": ["P1"],
                "PeptideSequence": ["PEPA"],
                "SampleID": ["S1"],
                "Intensity": [10.0],
            }
        ),
        protein_column="ProteinName",
        peptide_column="PeptideSequence",
        intensity_column="Intensity",
        sample_column="SampleID",
    )
    assert captured["num_cores"] == _resolve_directlfq_num_cores(-1)
    assert captured["num_cores"] == (os.cpu_count() or 1)


def test_streaming_directlfq_preserves_empty_output_shape():
    lfq_config.set_global_protein_and_ion_id(protein_id="protein", quant_id="ion")
    normed_df = _make_normed_directlfq_frame()

    streamed_df = estimate_protein_intensities_streamed(
        normed_df,
        min_nonan=10,
        num_samples_quadratic=2,
    )

    assert list(streamed_df.columns) == ["protein", "S1", "S2", "S3", "S4"]
    assert streamed_df.empty
