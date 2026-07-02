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
from mokume.pipeline.config import InputConfig, PipelineConfig, RuntimeConfig


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
