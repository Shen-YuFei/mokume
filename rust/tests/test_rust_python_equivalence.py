"""Rust <-> pure-Python equivalence for ``features2proteins``.

The compiled Rust kernel (``mokume._mokume``) and the pure-Python pipeline are
two implementations of the same protein-quantification contract. They share the
import name ``mokume`` and are mutually exclusive in a single interpreter, and CI
installs them in separate jobs with different dependency sets (the wheel job has
``_mokume`` but not ``duckdb``; the pure-Python job has ``duckdb`` but not
``_mokume``). Running both live in one process is therefore impossible.

Instead we compare the Rust kernel's output against a *frozen* pure-Python
reference checked in beside the fixture. The references were produced by
``QuantificationPipeline(...).run()`` (float64); the Rust kernel computes in
``f32``, so values agree to ~1e-7 relative and we assert ``rtol=1e-6``.

Regenerate the goldens after an intentional pure-Python change::

    cd python && python -c "\
from mokume.pipeline.config import PipelineConfig, InputConfig, QuantificationConfig; \
from mokume.pipeline.features_to_proteins import QuantificationPipeline; \
QuantificationPipeline(PipelineConfig(\
input=InputConfig(parquet='../rust/tests/example/feature_wide.parquet'), \
quantification=QuantificationConfig(method='sum'))).run().to_csv(\
'../rust/tests/example/feature_wide_sum_python.csv', index=False)"

The TMM goldens (``feature_wide_{sum,median}_tmm_python.csv`` on the sparse
fixture, ``feature_wide_dense_{sum,median}_tmm_python.csv`` on the dense one) are
produced the same way but with
``normalization=NormalizationConfig(sample_method='tmm')`` and, for the dense
case, ``input=InputConfig(parquet='.../feature_wide_dense.parquet')``. The dense
fixture packs 20 shared canonical peptides across every one of the 10 samples so
each sample clears TMM's >=10-valid-feature guard, making TMM actually shift the
output; the sparse ``feature_wide`` fixture has too few overlapping features, so
TMM legitimately falls back to factor 1.0 in BOTH implementations and equals the
un-normalized result.
"""

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import pytest

mokume = pytest.importorskip("mokume")

REPO = Path(__file__).resolve().parents[2]
FEATURE_PARQUET = REPO / "rust" / "tests" / "example" / "feature_wide.parquet"
DENSE_PARQUET = REPO / "rust" / "tests" / "example" / "feature_wide_dense.parquet"
GOLDEN_DIR = REPO / "rust" / "tests" / "example"


def _rust_features2proteins(
    method: str,
    out_path: Path,
    parquet: Path = FEATURE_PARQUET,
    sample_normalization: Optional[str] = None,
) -> pd.DataFrame:
    """Quantify via the compiled Rust kernel (in-process)."""
    if not hasattr(mokume, "_mokume"):
        pytest.skip("mokume._mokume (Rust kernel) not available in this interpreter")
    args = [
        "features2proteins",
        "--parquet",
        str(parquet),
        "--quant-method",
        method,
        "--output",
        str(out_path),
        "--output-format",
        "python-compatible",
    ]
    if sample_normalization is not None:
        args += ["--sample-normalization", sample_normalization]
    try:
        mokume.run(args)
    except RuntimeError:
        # The clap/PyO3 FFI wrapper can surface a benign post-run signal; the
        # CSV is still written. Re-raise only if it genuinely produced nothing.
        if not out_path.exists():
            raise
    return pd.read_csv(out_path)


def _canonical(df: pd.DataFrame) -> pd.DataFrame:
    """Protein x sample matrix indexed by protein, columns sorted, NaN -> 0."""
    df = df.rename(columns={df.columns[0]: "ProteinName"})
    df = df.set_index("ProteinName").sort_index()
    df = df.reindex(sorted(df.columns), axis=1)
    return df.fillna(0.0)


@pytest.mark.parametrize("method", ["sum", "median"])
def test_rust_matches_python_golden(method, tmp_path):
    """Rust output matches the frozen pure-Python reference within f32 tolerance."""
    golden_path = GOLDEN_DIR / f"feature_wide_{method}_python.csv"
    rust_df = _canonical(
        _rust_features2proteins(method, tmp_path / f"rust_{method}.csv")
    )
    py_df = _canonical(pd.read_csv(golden_path))

    assert list(rust_df.index) == list(py_df.index), "protein sets differ"
    assert list(rust_df.columns) == list(py_df.columns), "sample columns differ"

    rust_vals = rust_df.to_numpy(dtype=float)
    py_vals = py_df.to_numpy(dtype=float)
    assert np.allclose(rust_vals, py_vals, rtol=1e-6, atol=1e-2), (
        f"{method}: Rust output diverges from the pure-Python golden beyond f32 "
        f"tolerance (max abs diff {np.max(np.abs(rust_vals - py_vals)):.4g})"
    )


@pytest.mark.parametrize("method", ["sum", "median"])
@pytest.mark.parametrize(
    "fixture, parquet, golden_prefix",
    [
        ("sparse", FEATURE_PARQUET, "feature_wide"),
        ("dense", DENSE_PARQUET, "feature_wide_dense"),
    ],
)
def test_rust_tmm_matches_python_golden(
    method, fixture, parquet, golden_prefix, tmp_path
):
    """Rust ``--sample-normalization tmm`` matches the pure-Python TMM golden.

    Covers both the sparse fixture (where TMM falls back to factor 1.0 in both
    kernels) and the dense fixture (where TMM genuinely rescales), so the test
    pins the guard behaviour AND the active-normalization math.
    """
    golden_path = GOLDEN_DIR / f"{golden_prefix}_{method}_tmm_python.csv"
    rust_df = _canonical(
        _rust_features2proteins(
            method,
            tmp_path / f"rust_{fixture}_{method}_tmm.csv",
            parquet=parquet,
            sample_normalization="tmm",
        )
    )
    py_df = _canonical(pd.read_csv(golden_path))

    assert list(rust_df.index) == list(py_df.index), "protein sets differ"
    assert list(rust_df.columns) == list(py_df.columns), "sample columns differ"

    rust_vals = rust_df.to_numpy(dtype=float)
    py_vals = py_df.to_numpy(dtype=float)
    assert np.allclose(rust_vals, py_vals, rtol=1e-6, atol=1e-2), (
        f"tmm/{fixture}/{method}: Rust output diverges from the pure-Python TMM "
        f"golden beyond f32 tolerance "
        f"(max abs diff {np.max(np.abs(rust_vals - py_vals)):.4g})"
    )


@pytest.mark.parametrize("method", ["sum", "median"])
def test_rust_tmm_is_not_a_noop_on_dense(method, tmp_path):
    """On the dense fixture, Rust TMM must differ from ``none`` (regression guard).

    Before the Rust TMM implementation, ``--sample-normalization tmm`` was
    silently equivalent to ``none``; this asserts the two now diverge on a
    fixture whose samples clear TMM's >=10-valid-feature guard.
    """
    none_df = _canonical(
        _rust_features2proteins(
            method,
            tmp_path / f"rust_dense_{method}_none.csv",
            parquet=DENSE_PARQUET,
            sample_normalization="none",
        )
    )
    tmm_df = _canonical(
        _rust_features2proteins(
            method,
            tmp_path / f"rust_dense_{method}_tmm.csv",
            parquet=DENSE_PARQUET,
            sample_normalization="tmm",
        )
    )

    assert list(none_df.index) == list(tmm_df.index)
    assert list(none_df.columns) == list(tmm_df.columns)
    max_diff = np.max(np.abs(none_df.to_numpy(float) - tmm_df.to_numpy(float)))
    assert max_diff > 1.0, (
        f"tmm/{method}: TMM produced output identical to none on the dense "
        f"fixture (max abs diff {max_diff:.4g}) -- TMM is a no-op"
    )
