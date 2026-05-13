"""Tests for mokume.agentic.proposer."""

from mokume.agentic.config import AgenticConfig
from mokume.agentic.profiler import profile_data
from mokume.agentic.proposer import propose_configs, _items_to_configs


def test_propose_configs_no_llm(synthetic_protein_df, sample_to_condition):
    """propose_configs with use_llm=False falls back to rules."""
    profile = profile_data(synthetic_protein_df, sample_to_condition)
    config = AgenticConfig(use_llm=False)
    configs = propose_configs(profile, config)
    assert len(configs) > 0
    assert all(hasattr(c, "de_method") for c in configs)


def test_propose_configs_llm_fallback(synthetic_protein_df, sample_to_condition):
    """propose_configs with missing LLM package falls back gracefully."""
    profile = profile_data(synthetic_protein_df, sample_to_condition)
    config = AgenticConfig(use_llm=True)
    # This should fall back to rule-based (no API key / no package)
    configs = propose_configs(profile, config)
    assert len(configs) > 0


def test_items_to_configs():
    """Convert parsed JSON items to CandidateConfig list."""
    items = [
        {
            "name": "cfg1",
            "de_method": "deqms",
            "fdr_method": "ihw",
            "normalization": "median",
            "imputation": "none",
            "log2fc_threshold": 0.5,
            "reasoning": "test",
            "expected_outcome": "good",
        }
    ]
    configs = _items_to_configs(items)
    assert len(configs) == 1
    assert configs[0].de_method == "deqms"
    assert configs[0].normalization == "median"


def test_items_to_configs_ensemble_passthrough():
    """Ensemble + ensemble_k fields are forwarded into CandidateConfig."""
    items = [
        {
            "name": "ens",
            "de_method": "ensemble",
            "fdr_method": "bh",
            "normalization": "rlr",
            "imputation": "seqknn",
            "ensemble": "limma,deqms,proda",
            "ensemble_k": 3,
            "log2fc_threshold": 1.0,
            "reasoning": "consensus",
            "expected_outcome": "robust",
        }
    ]
    configs = _items_to_configs(items)
    assert configs[0].de_method == "ensemble"
    assert configs[0].ensemble == "limma,deqms,proda"
    assert configs[0].ensemble_k == 3


def test_items_to_configs_ensemble_defaults():
    """Missing ensemble fields default to 'none' / 2 (backward compat)."""
    items = [
        {
            "name": "no_ens",
            "de_method": "limma",
            "fdr_method": "bh",
            "normalization": "none",
            "imputation": "none",
            "log2fc_threshold": 0.5,
            "reasoning": "baseline",
            "expected_outcome": "stable",
        }
    ]
    configs = _items_to_configs(items)
    assert configs[0].ensemble == "none"
    assert configs[0].ensemble_k == 2


def test_items_to_configs_new_imputation_methods():
    """New imputation enum values round-trip through the proposer."""
    items = [
        {
            "name": f"cfg_{m}",
            "de_method": "limma",
            "fdr_method": "bh",
            "normalization": "none",
            "imputation": m,
            "log2fc_threshold": 0.5,
            "reasoning": "test",
            "expected_outcome": "ok",
        }
        for m in ("qrilc", "mle", "mice", "nbavg", "gms")
    ]
    configs = _items_to_configs(items)
    assert [c.imputation for c in configs] == ["qrilc", "mle", "mice", "nbavg", "gms"]
