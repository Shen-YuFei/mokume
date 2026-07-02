"""Guard-path tests for :mod:`mokume.pipeline.flows.rust`.

The module must import cleanly even when the compiled ``mokume._mokume``
extension is absent (the pure-Python environment). When the kernel is missing,
calling the flow must raise a clear ``RuntimeError`` naming ``mokume-rs`` /
``_mokume``. When the wheel IS present, the raise-assertion is skipped so this
test never fails in the compiled environment.
"""

from __future__ import annotations

import pytest

import mokume.pipeline.flows.rust as rust_flow
from mokume.pipeline.config import (
    BatchCorrectionConfig,
    InputConfig,
    PipelineConfig,
    RuntimeConfig,
)


def _rust_config() -> PipelineConfig:
    return PipelineConfig(
        input=InputConfig(parquet="in.parquet"),
        runtime=RuntimeConfig(backend="rust"),
    )


def test_module_imports_without_extension() -> None:
    # Importing the flow must never raise, kernel present or not.
    assert hasattr(rust_flow, "run")


def test_runtime_config_accepts_rust_backend() -> None:
    cfg = _rust_config()
    assert cfg.runtime.backend == "rust"


def test_runtime_config_rejects_unknown_backend() -> None:
    with pytest.raises(ValueError, match="backend"):
        RuntimeConfig(backend="julia")


def test_run_raises_when_kernel_absent() -> None:
    if rust_flow._kernel_run is not None:
        pytest.skip("mokume._mokume is installed; guard-raise path not exercised")

    with pytest.raises(RuntimeError) as excinfo:
        rust_flow.run(None, _rust_config())

    message = str(excinfo.value)
    assert "mokume-rs" in message
    assert "_mokume" in message


def test_run_pipeline_raises_when_kernel_absent() -> None:
    if rust_flow._kernel_run is not None:
        pytest.skip("mokume._mokume is installed; guard-raise path not exercised")

    from mokume.pipeline.runner import run_pipeline

    with pytest.raises(RuntimeError, match="mokume-rs"):
        run_pipeline(_rust_config())


@pytest.mark.parametrize("backend,expected_batch_calls", [("rust", 0), ("python", 1)])
def test_postprocess_batch_correction_guarded_by_backend(
    backend, expected_batch_calls, monkeypatch
) -> None:
    """ComBat runs in ``_postprocess`` only for the python backend.

    The Rust kernel already applies ComBat (cli_args emits ``--batch-correction``),
    so re-running it in ``_postprocess`` would double-correct the matrix.
    """
    import pandas as pd

    from mokume.core.dataset import QpxDataset
    from mokume.pipeline import runner
    from mokume.pipeline.stages import PostprocessingStage

    calls: list[int] = []
    monkeypatch.setattr(
        PostprocessingStage,
        "apply_batch_correction",
        lambda self, df, dataset=None: (calls.append(1), df)[1],
    )

    dataset = QpxDataset()
    dataset.proteins = pd.DataFrame(
        {"ProteinName": ["P1", "P2"], "S1": [1.0, 2.0], "S2": [3.0, 4.0]}
    )
    config = PipelineConfig(
        input=InputConfig(parquet="in.parquet"),
        runtime=RuntimeConfig(backend=backend),
        batch=BatchCorrectionConfig(enabled=True, method="run"),
    )

    runner._postprocess(dataset, config)
    assert len(calls) == expected_batch_calls
