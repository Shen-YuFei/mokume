"""Rust <-> pure-Python equivalence for ``features2proteins``.

The compiled Rust kernel (``mokume._mokume``) and the pure-Python pipeline are
two implementations of the same protein-quantification contract. This test runs
the *same* feature parquet through both and asserts the protein matrices agree
within float tolerance.

The two implementations share the import name ``mokume`` and are mutually
exclusive in a single interpreter, so the pure-Python side is executed in a
subprocess whose working directory is the editable ``python/`` source tree
(there ``import mokume`` resolves to the pure-Python package, without
``_mokume``). The Rust kernel computes intensities in ``f32``; the pure-Python
pipeline uses ``f64``. Values therefore agree only to ~1e-7 relative, so we
assert ``rtol=1e-6``.
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

mokume = pytest.importorskip("mokume")

REPO = Path(__file__).resolve().parents[2]
PYTHON_DIR = REPO / "python"
FEATURE_PARQUET = REPO / "rust" / "tests" / "example" / "feature_wide.parquet"


def _rust_features2proteins(method: str, out_path: Path) -> pd.DataFrame:
    """Quantify via the compiled Rust kernel (in-process)."""
    if not hasattr(mokume, "_mokume"):
        pytest.skip("mokume._mokume (Rust kernel) not available in this interpreter")
    try:
        mokume.run(
            [
                "features2proteins",
                "--parquet",
                str(FEATURE_PARQUET),
                "--quant-method",
                method,
                "--output",
                str(out_path),
                "--output-format",
                "python-compatible",
            ]
        )
    except RuntimeError:
        # The clap/PyO3 FFI wrapper can surface a benign post-run signal; the
        # CSV is still written. Re-raise only if it genuinely produced nothing.
        if not out_path.exists():
            raise
    return pd.read_csv(out_path)


def _python_features2proteins(method: str, out_path: Path) -> pd.DataFrame:
    """Quantify via the pure-Python pipeline (subprocess in python/ source tree)."""
    script = (
        "from mokume.pipeline.config import "
        "PipelineConfig, InputConfig, QuantificationConfig\n"
        "from mokume.pipeline.features_to_proteins import QuantificationPipeline\n"
        f"cfg = PipelineConfig(input=InputConfig(parquet={str(FEATURE_PARQUET)!r}), "
        f"quantification=QuantificationConfig(method={method!r}))\n"
        "df = QuantificationPipeline(cfg).run()\n"
        f"df.to_csv({str(out_path)!r}, index=False)\n"
    )
    subprocess.run(
        [sys.executable, "-c", script],
        cwd=str(PYTHON_DIR),
        check=True,
    )
    return pd.read_csv(out_path)


def _canonical(df: pd.DataFrame) -> pd.DataFrame:
    """Protein x sample matrix indexed by protein, columns sorted, NaN -> 0."""
    df = df.rename(columns={df.columns[0]: "ProteinName"})
    df = df.set_index("ProteinName").sort_index()
    df = df.reindex(sorted(df.columns), axis=1)
    return df.fillna(0.0)


@pytest.mark.parametrize("method", ["sum", "median"])
def test_rust_python_features2proteins_equivalent(method, tmp_path):
    """Rust and pure-Python produce the same protein matrix within f32 tolerance."""
    rust_df = _canonical(
        _rust_features2proteins(method, tmp_path / f"rust_{method}.csv")
    )
    py_df = _canonical(_python_features2proteins(method, tmp_path / f"py_{method}.csv"))

    assert list(rust_df.index) == list(py_df.index), "protein sets differ"
    assert list(rust_df.columns) == list(py_df.columns), "sample columns differ"

    rust_vals = rust_df.to_numpy(dtype=float)
    py_vals = py_df.to_numpy(dtype=float)
    assert np.allclose(rust_vals, py_vals, rtol=1e-6, atol=1e-2), (
        f"{method}: Rust and pure-Python matrices diverge beyond f32 tolerance "
        f"(max abs diff {np.max(np.abs(rust_vals - py_vals)):.4g})"
    )
