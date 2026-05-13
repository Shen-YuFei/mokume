"""Regression tests for the LLM tool schema in mokume.agentic.llm_client."""

from mokume.agentic.llm_client import PROPOSAL_TOOL, _CONFIG_ITEM_SCHEMA


def _config_props() -> dict:
    return _CONFIG_ITEM_SCHEMA["properties"]


def test_de_method_includes_ensemble():
    """ensemble must be a permitted DE method choice."""
    assert "ensemble" in _config_props()["de_method"]["enum"]


def test_imputation_full_catalogue():
    """All runner-supported imputers are advertised to the LLM."""
    expected = {
        "none",
        "minprob",
        "mindet",
        "knn",
        "missforest",
        "seqknn",
        "qrilc",
        "mle",
        "mice",
        "nbavg",
        "gms",
        "bpca",
        "impseq",
        "impseqrob",
    }
    assert expected.issubset(set(_config_props()["imputation"]["enum"]))


def test_normalization_includes_mbqn_loess():
    """MBQN and LOESS must be exposed alongside the legacy choices."""
    norms = set(_config_props()["normalization"]["enum"])
    assert {"mbqn", "loess"}.issubset(norms)


def test_ensemble_fields_present_and_required():
    """ensemble + ensemble_k are declared and marked required (strict mode)."""
    props = _config_props()
    assert "ensemble" in props
    assert "ensemble_k" in props
    required = set(_CONFIG_ITEM_SCHEMA["required"])
    assert {"ensemble", "ensemble_k"}.issubset(required)


def test_strict_tool_mode_preserved():
    """Schema stays OpenAI strict-compatible (additionalProperties=False)."""
    assert _CONFIG_ITEM_SCHEMA["additionalProperties"] is False
    assert PROPOSAL_TOOL["function"]["strict"] is True
    assert PROPOSAL_TOOL["function"]["parameters"]["additionalProperties"] is False


def test_ensemble_k_bounds():
    """ensemble_k must be a bounded integer."""
    spec = _config_props()["ensemble_k"]
    assert spec["type"] == "integer"
    assert spec["minimum"] == 1
    assert spec["maximum"] == 5
