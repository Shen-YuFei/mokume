"""
End-to-end tests for the OOP entry point ``run_pipeline``.

``mokume.pipeline.runner.run_pipeline`` is the documented metadata-driven
entry point: it resolves the quantification method from the plugin registry,
dispatches by ``input_level`` to one of the ``flows/`` modules, and runs
common post-processing. The production CLI uses ``features_to_proteins``
(i.e. ``QuantificationPipeline.run_dataset``) instead, so before these tests
no coverage exercised ``run_pipeline`` against real data.

The correctness oracle here is that, for the pure-Python backend,
``run_pipeline(config)`` must produce a protein matrix identical (within
1e-9) to ``QuantificationPipeline(config).run_dataset().get_level('proteins')``
for the same config. A divergence would mean a flow loads or normalizes
differently from the reference (e.g. double-IRS or double-normalization).

Both entry points are driven off the small dev fixture
``tests/example/feature_wide.parquet`` (a real PXD020192 LFQ feature slice
with its SDRF), so the test is self-contained.
"""

import os

import numpy as np
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
)
from mokume.pipeline.runner import run_pipeline  # noqa: E402
from mokume.pipeline.features_to_proteins import (  # noqa: E402
    QuantificationPipeline,
)


def _make_config(quant_method: str) -> PipelineConfig:
    """Build a PipelineConfig pointing at the dev fixture.

    ``min_unique_peptides=1`` is used because the fixture is a tiny slice;
    the default of 2 leaves too few proteins to be a useful signal.
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
        quantification=QuantificationConfig(method=quant_method),
    )


def _protein_matrix(dataset):
    """Return the protein matrix sorted by protein id and sample columns.

    Normalizes the first (protein-id) column name to ``"Protein"`` so the two
    entry points can be compared even if they label that column differently.
    """
    df = dataset.get_level("proteins")
    assert df is not None, "proteins level is None"
    df = df.copy()
    protein_col = df.columns[0]
    sample_cols = sorted(c for c in df.columns if c != protein_col)
    df = df.rename(columns={protein_col: "Protein"})
    df = df.sort_values("Protein").reset_index(drop=True)
    return df[["Protein", *sample_cols]], sample_cols


def _max_abs_diff(quant_method: str) -> float:
    """Run both entry points and return the max abs diff of protein matrices."""
    ds_run = run_pipeline(_make_config(quant_method))
    ds_ref = QuantificationPipeline(_make_config(quant_method)).run_dataset()

    m_run, cols_run = _protein_matrix(ds_run)
    m_ref, cols_ref = _protein_matrix(ds_ref)

    # Same proteins, same samples, same order.
    assert list(m_run["Protein"]) == list(m_ref["Protein"]), (
        f"{quant_method}: protein index differs between run_pipeline and run_dataset"
    )
    assert cols_run == cols_ref, (
        f"{quant_method}: sample columns differ "
        f"(run_pipeline={cols_run}, run_dataset={cols_ref})"
    )

    a_run = m_run[cols_run].to_numpy(dtype=float)
    a_ref = m_ref[cols_ref].to_numpy(dtype=float)

    # NaN in the same cell on both sides is agreement, not a diff.
    both_nan = np.isnan(a_run) & np.isnan(a_ref)
    diff = np.abs(np.where(both_nan, 0.0, a_run - a_ref))
    return float(np.nanmax(diff)) if diff.size else 0.0


class TestRunPipelineMatchesRunDataset:
    """run_pipeline must equal run_dataset within 1e-9 (pure-Python backend)."""

    @pytest.mark.parametrize("quant_method", ["median", "sum"])
    def test_standard_flow_matches(self, quant_method):
        assert _max_abs_diff(quant_method) < 1e-9

    def test_directlfq_flow_matches(self):
        pytest.importorskip("directlfq")
        pytest.importorskip("polars")
        assert _max_abs_diff("directlfq") < 1e-9


class TestRunPipelineProvenance:
    """run_pipeline populates proteins and records provenance steps."""

    def test_median_populates_proteins_and_provenance(self):
        dataset = run_pipeline(_make_config("median"))
        assert dataset.proteins is not None
        assert len(dataset.proteins) > 0

        steps = dataset.uns["provenance"]["steps"]
        step_names = [s["name"] for s in steps]
        assert "loading" in step_names
        assert "quantification" in step_names
