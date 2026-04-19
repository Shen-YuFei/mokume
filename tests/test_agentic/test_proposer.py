"""Tests for mokume.agentic.proposer."""

from mokume.agentic.config import AgenticConfig
from mokume.agentic.profiler import profile_data
from mokume.agentic.proposer import propose_configs, _parse_configs_from_json


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


def test_parse_configs_from_json():
    """Parse well-formed LLM JSON response."""
    raw = '''Here are my suggestions:
    [
        {"name": "cfg1", "de_method": "deqms", "fdr_method": "ihw",
         "imputation": "none", "log2fc_threshold": 0.5,
         "reasoning": "test", "expected_outcome": "good"}
    ]
    '''
    configs = _parse_configs_from_json(raw)
    assert len(configs) == 1
    assert configs[0].de_method == "deqms"
