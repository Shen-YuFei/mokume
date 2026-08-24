from __future__ import annotations

import pytest

from mokume.agentic.contract import method_contract, validate_config_values
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
