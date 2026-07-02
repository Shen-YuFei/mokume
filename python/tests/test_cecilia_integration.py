"""
Integration test for the plugin-architecture pipeline (``run_dataset``).

Ported from main's ``tests/test_cecilia_integration.py`` but made runnable
without main's private cecilia NASH TMT dataset: it points at the small
dev fixture ``tests/example/feature_wide.parquet`` (a real PXD020192 LFQ
feature slice with its SDRF ``tests/example/PXD020192.sdrf.tsv``) so the
test is self-contained and green in CI.

It exercises the same contract main's test does:

  * ``PipelineConfig`` -> ``QuantificationPipeline(config).run_dataset()``
  * multiple quantification methods (median, maxlfq, top3 via the TopN
    pattern; ratio is skipped because no TMT/PSM fixture with reference
    samples ships with the dev branch)
  * median + IRS post-processing (guarded: skipped if the fixture has no
    reference samples to scale against)

Each test validates:
  1. The pipeline completes without error.
  2. ``dataset.proteins`` is populated and non-empty.
  3. ``dataset.uns`` provenance records at least two steps
     (``loading`` and ``quantification``).
"""

import os
import tempfile

import pandas as pd
import pytest

# Paths ---------------------------------------------------------------
EXAMPLE_DIR = os.path.join(os.path.dirname(__file__), "example")
PARQUET = os.path.join(EXAMPLE_DIR, "feature_wide.parquet")
SDRF = os.path.join(EXAMPLE_DIR, "PXD020192.sdrf.tsv")

# Skip the whole module only if the dev fixture itself is missing.
pytestmark = pytest.mark.skipif(
    not (os.path.exists(PARQUET) and os.path.exists(SDRF)),
    reason="Dev feature_wide fixture not available",
)

# Imports (after path setup) -------------------------------------------
from mokume.pipeline.config import (  # noqa: E402
    PipelineConfig,
    InputConfig,
    FilterConfig,
    NormalizationConfig,
    QuantificationConfig,
    IRSConfig,
    OutputConfig,
)
from mokume.pipeline.features_to_proteins import QuantificationPipeline  # noqa: E402
from mokume.core.dataset import QpxDataset  # noqa: E402


def _make_config(
    quant_method: str = "median",
    irs: bool = False,
    irs_remove_ref: bool = False,
    coverage_threshold: float = None,
    export_peptides: str = None,
) -> PipelineConfig:
    """Build a PipelineConfig pointing at the dev fixture.

    ``min_unique_peptides=1`` is used because the fixture is a tiny 500-row
    slice; the default of 2 leaves too few proteins to be a useful signal.
    """
    return PipelineConfig(
        input=InputConfig(parquet=PARQUET, sdrf=SDRF),
        filtering=FilterConfig(
            min_aa=7,
            min_unique_peptides=1,
            remove_contaminants=True,
        ),
        normalization=NormalizationConfig(
            run_method="median",
            sample_method="globalMedian",
        ),
        quantification=QuantificationConfig(
            method=quant_method,
            coverage_threshold=coverage_threshold,
        ),
        irs=IRSConfig(
            enabled=irs,
            remove_reference=irs_remove_ref,
        ),
        output=OutputConfig(
            export_peptides=export_peptides,
        ),
    )


def _assert_valid_result(dataset: QpxDataset, min_proteins: int = 1):
    """Common assertions for all pipeline results."""
    # Proteins must be populated and non-empty.
    assert dataset.proteins is not None, "proteins level is None"

    protein_df = dataset.proteins
    if isinstance(protein_df, pd.DataFrame):
        assert len(protein_df) >= min_proteins, (
            f"Expected >= {min_proteins} proteins, got {len(protein_df)}"
        )

    # Provenance must record at least loading + quantification.
    assert "provenance" in dataset.uns, "No provenance in uns"
    steps = dataset.uns["provenance"]["steps"]
    assert len(steps) >= 2, f"Expected >= 2 provenance steps, got {len(steps)}"

    step_names = [s["name"] for s in steps]
    assert "loading" in step_names
    assert "quantification" in step_names


# ======================================================================
# Tests
# ======================================================================


class TestMedianQuantification:
    """Standard flow: median summarization."""

    def test_median_pipeline(self):
        config = _make_config(quant_method="median")
        pipeline = QuantificationPipeline(config)
        dataset = pipeline.run_dataset()
        _assert_valid_result(dataset, min_proteins=10)

        # Quantification provenance carries the method name.
        quant_steps = [
            s
            for s in dataset.uns["provenance"]["steps"]
            if s["name"] == "quantification"
        ]
        assert quant_steps and quant_steps[0].get("method", "").lower() == "median"


class TestMaxLFQQuantification:
    """Standard flow: MaxLFQ (routes through the DirectLFQ-streaming path)."""

    def test_maxlfq_pipeline(self):
        config = _make_config(quant_method="maxlfq")
        pipeline = QuantificationPipeline(config)
        dataset = pipeline.run_dataset()
        # MaxLFQ's min-nonan filtering is aggressive on this tiny fixture, so
        # only require a non-empty result rather than a large protein count.
        _assert_valid_result(dataset, min_proteins=1)


class TestTopNQuantification:
    """Standard flow: TopN pattern (top3 resolved via the registry)."""

    def test_top3_pipeline(self):
        config = _make_config(quant_method="top3")
        pipeline = QuantificationPipeline(config)
        dataset = pipeline.run_dataset()
        _assert_valid_result(dataset, min_proteins=10)


class TestRatioQuantification:
    """Ratio flow: PS protocol log2 ratios.

    Skipped on the dev branch: the ratio flow requires a TMT PSM parquet
    plus an SDRF declaring reference/pooled samples, and no such fixture
    ships here (``feature_wide.parquet`` is LFQ and PXD020192 has no pooled
    samples). Kept as a skipped placeholder so the contract is documented.
    """

    @pytest.mark.skip(reason="No TMT/PSM ratio fixture with reference samples")
    def test_ratio_pipeline(self):  # pragma: no cover - documented placeholder
        config = _make_config(quant_method="ratio", coverage_threshold=0.65)
        pipeline = QuantificationPipeline(config)
        dataset = pipeline.run_dataset()
        _assert_valid_result(dataset, min_proteins=1)
        assert "ratio_config" in dataset.uns
        assert dataset.uns["ratio_config"]["reference_samples"]


class TestIRSNormalization:
    """Standard flow + IRS post-processing.

    The dev fixture (PXD020192) has no pooled/reference samples, so IRS is
    expected to no-op. The test asserts the pipeline still completes and,
    when IRS *did* run, that it is reflected in provenance.
    """

    def test_median_irs_pipeline(self):
        config = _make_config(
            quant_method="median",
            irs=True,
            irs_remove_ref=True,
        )
        pipeline = QuantificationPipeline(config)
        dataset = pipeline.run_dataset()
        _assert_valid_result(dataset, min_proteins=10)

        step_names = [s["name"] for s in dataset.uns["provenance"]["steps"]]
        if "irs_normalization" in step_names:
            # When IRS ran, provenance must reflect it (>= 3 steps).
            assert len(dataset.uns["provenance"]["steps"]) >= 3


class TestQpxDatasetSaveLoad:
    """Save/load roundtrip with real pipeline output."""

    def test_save_load_roundtrip(self):
        config = _make_config(quant_method="median")
        pipeline = QuantificationPipeline(config)
        dataset = pipeline.run_dataset()

        with tempfile.TemporaryDirectory() as tmpdir:
            dataset.save(tmpdir)
            loaded = QpxDataset.load(tmpdir)

            assert loaded.proteins is not None
            assert len(loaded.proteins) == len(dataset.proteins)
            assert loaded.uns.get("provenance") is not None


class TestPeptideExport:
    """Peptide export path populates a CSV and the peptides level."""

    def test_export_peptides(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            pep_path = os.path.join(tmpdir, "peptides.csv")
            config = _make_config(
                quant_method="median",
                export_peptides=pep_path,
            )
            pipeline = QuantificationPipeline(config)
            dataset = pipeline.run_dataset()
            _assert_valid_result(dataset, min_proteins=10)

            assert os.path.exists(pep_path), "Peptide export file not created"
            pep_df = pd.read_csv(pep_path)
            assert len(pep_df) > 0


class TestSchemaValidation:
    """Schema validation is wired in for the peptides level."""

    def test_peptides_schema_valid(self):
        config = _make_config(quant_method="median")
        pipeline = QuantificationPipeline(config)
        dataset = pipeline.run_dataset()

        if dataset.peptides is not None:
            errors = dataset.validate_level("peptides")
            assert errors == [], f"Schema errors: {errors}"


class TestProvenanceProvenanceDepth:
    """The core P3 contract: >= 2 provenance steps recorded."""

    def test_at_least_two_steps(self):
        config = _make_config(quant_method="median")
        pipeline = QuantificationPipeline(config)
        dataset = pipeline.run_dataset()
        assert len(dataset.uns["provenance"]["steps"]) >= 2
