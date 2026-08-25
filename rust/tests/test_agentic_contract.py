from __future__ import annotations

import pytest
import pandas as pd

from mokume.agentic.contract import method_contract, validate_config_values
from mokume.agentic.runner import ExperimentContext, _run_rust_de
from mokume.agentic.state import CandidateConfig


@pytest.mark.parametrize("de_method", ["rots", "limrots"])
@pytest.mark.parametrize("fdr_method", ["ihw", "bky", "storey"])
def test_standalone_rots_family_rejects_ineffective_fdr(
    de_method: str,
    fdr_method: str,
) -> None:
    config = CandidateConfig(
        name="ineffective-fdr",
        de_method=de_method,
        fdr_method=fdr_method,
    )

    with pytest.raises(ValueError, match="requires fdr_method='bh'"):
        validate_config_values(config.to_dict())


@pytest.mark.parametrize("de_method", ["rots", "limrots"])
def test_standalone_rots_family_accepts_native_fdr(de_method: str) -> None:
    validate_config_values(
        CandidateConfig(name="native-fdr", de_method=de_method).to_dict()
    )


def test_method_contract_exposes_conditional_fdr_catalog() -> None:
    contract = method_contract()

    assert contract["fdr_method_by_de_method"]["rots"] == ["bh"]
    assert contract["fdr_method_by_de_method"]["limrots"] == ["bh"]
    assert contract["fdr_method_by_de_method"]["ensemble"] == [
        "bh",
        "ihw",
        "bky",
        "storey",
    ]


def test_nonensemble_rejects_fake_ensemble_k_axis() -> None:
    """Single-method candidates must not expose an inactive ensemble axis."""
    config = CandidateConfig(name="single-method", ensemble_k=2)
    with pytest.raises(ValueError, match="must be null"):
        validate_config_values(config.to_dict())


def test_ensemble_requires_bounded_ensemble_k() -> None:
    """Ensemble candidates must accept a valid bounded voting threshold."""
    config = CandidateConfig(
        name="ensemble",
        de_method="ensemble",
        ensemble="limma,deqms,proda",
        ensemble_k=2,
    )
    validate_config_values(config.to_dict())


def test_runner_maps_auto_gate_and_omits_rots_fdr(monkeypatch) -> None:
    captured: dict[str, object] = {}

    def fake_de(*args, **kwargs):
        del args
        captured.update(kwargs)
        return []

    monkeypatch.setattr("mokume.agentic.runner.differential_expression", fake_de)
    frame = pd.DataFrame(
        {
            "ProteinName": ["P1"],
            "S1": [1.0],
            "S2": [2.0],
            "S3": [3.0],
            "S4": [4.0],
        }
    )
    context = ExperimentContext(
        sample_to_condition={"S1": "A", "S2": "A", "S3": "B", "S4": "B"},
        contrast=("A", "B"),
    )
    config = CandidateConfig(
        name="limrots-auto",
        de_method="limrots",
        log2fc_threshold="auto",
    )

    _run_rust_de(frame, config, context, None)

    assert captured["effect_size_gate"] == "mixture"
    assert "fdr_method" not in captured
