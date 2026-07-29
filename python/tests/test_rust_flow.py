"""Guard-path tests for :mod:`mokume.pipeline.flows.rust`.

The module must import cleanly even when the compiled ``mokume._mokume``
extension is absent (the pure-Python environment). When the kernel is missing,
calling the flow must raise a clear ``RuntimeError`` naming ``mokume-rs`` /
``_mokume``. When the wheel IS present, the raise-assertion is skipped so this
test never fails in the compiled environment.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

import mokume.pipeline.flows.rust as rust_flow
from mokume.pipeline.config import (
    BatchCorrectionConfig,
    ImputationConfig,
    InputConfig,
    PipelineConfig,
    QuantificationConfig,
    RuntimeConfig,
)


def _rust_config(**overrides) -> PipelineConfig:
    overrides.setdefault("input", InputConfig(parquet="in.parquet"))
    overrides.setdefault("runtime", RuntimeConfig(backend="rust"))
    return PipelineConfig(**overrides)


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


@pytest.mark.parametrize(
    "config",
    [
        _rust_config(
            runtime=RuntimeConfig(backend="rust", duckdb_memory="8GB"),
        ),
        _rust_config(
            runtime=RuntimeConfig(backend="rust", duckdb_threads=2),
        ),
        _rust_config(
            quantification=QuantificationConfig(ion_alignment="hierarchical"),
        ),
    ],
)
def test_invalid_config_fails_before_kernel_or_temp_output(config, monkeypatch) -> None:
    def unexpected_call(*args, **kwargs):
        raise AssertionError("validation must happen before side effects")

    monkeypatch.setattr(rust_flow, "_kernel_run", unexpected_call)
    monkeypatch.setattr(rust_flow.tempfile, "TemporaryDirectory", unexpected_call)

    with pytest.raises(ValueError):
        rust_flow.run(None, config)


def test_run_reads_kernel_output_and_cleans_temp_directory(monkeypatch) -> None:
    output_paths: list[Path] = []
    observed_argv: list[str] = []

    def fake_kernel(argv) -> None:
        observed_argv.extend(argv)
        output = Path(argv[argv.index("--output") + 1])
        output_paths.append(output)
        pd.DataFrame(
            {
                "ProteinName": ["P1", "P2"],
                "sample-a": [100.0, 25.0],
                "sample-b": [50.0, 12.5],
            }
        ).to_csv(output, index=False)

    monkeypatch.setattr(rust_flow, "_kernel_run", fake_kernel)
    dataset = rust_flow.run(None, _rust_config())

    assert observed_argv[0] == "features2proteins"
    assert dataset.proteins.to_dict(orient="list") == {
        "ProteinName": ["P1", "P2"],
        "sample-a": [100.0, 25.0],
        "sample-b": [50.0, 12.5],
    }
    assert output_paths and not output_paths[0].parent.exists()


@pytest.mark.parametrize("create_empty", [False, True])
def test_run_rejects_missing_or_empty_kernel_output(create_empty, monkeypatch) -> None:
    output_paths: list[Path] = []

    def incomplete_kernel(argv) -> None:
        output = Path(argv[argv.index("--output") + 1])
        output_paths.append(output)
        if create_empty:
            output.touch()

    monkeypatch.setattr(rust_flow, "_kernel_run", incomplete_kernel)

    with pytest.raises(
        RuntimeError, match="without creating a non-empty protein output"
    ):
        rust_flow.run(None, _rust_config())

    assert output_paths and not output_paths[0].parent.exists()


def test_kernel_failure_removes_partial_output(monkeypatch) -> None:
    output_paths: list[Path] = []

    def failing_kernel(argv) -> None:
        output = Path(argv[argv.index("--output") + 1])
        output_paths.append(output)
        output.write_text("partial", encoding="utf-8")
        raise RuntimeError("kernel failed")

    monkeypatch.setattr(rust_flow, "_kernel_run", failing_kernel)

    with pytest.raises(RuntimeError, match="kernel failed"):
        rust_flow.run(None, _rust_config())

    assert output_paths and not output_paths[0].parent.exists()


@pytest.mark.parametrize("backend", ["rust", "python"])
def test_postprocess_applies_batch_for_both_backends(backend, monkeypatch) -> None:
    """ComBat runs in ``_postprocess`` for BOTH backends.

    Under the Python-owned post-processing model the Rust kernel only loads,
    filters, normalizes and quantifies; batch correction is applied in Python
    for either backend, so it must run exactly once regardless of ``backend``.
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
    assert len(calls) == 1


@pytest.mark.parametrize("backend", ["rust", "python"])
def test_postprocess_applies_imputation_when_enabled(backend, monkeypatch) -> None:
    """Imputation runs in ``_postprocess`` for BOTH backends when enabled.

    The Rust kernel does not impute, so ``ImputationStage.impute`` must be
    invoked from Python whenever ``config.imputation.enabled`` is set.
    """
    import pandas as pd

    from mokume.core.dataset import QpxDataset
    from mokume.pipeline import runner
    from mokume.pipeline.stages import ImputationStage

    calls: list[int] = []
    monkeypatch.setattr(
        ImputationStage,
        "impute",
        lambda self, df: (calls.append(1), df)[1],
    )

    dataset = QpxDataset()
    dataset.proteins = pd.DataFrame(
        {"ProteinName": ["P1", "P2"], "S1": [1.0, 2.0], "S2": [3.0, 4.0]}
    )
    config = PipelineConfig(
        input=InputConfig(parquet="in.parquet"),
        runtime=RuntimeConfig(backend=backend),
        imputation=ImputationConfig(enabled=True, method="minprob"),
    )

    runner._postprocess(dataset, config)
    assert len(calls) == 1
