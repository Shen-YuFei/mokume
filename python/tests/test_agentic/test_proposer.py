"""Tests for mokume.agentic.proposer."""

from mokume.agentic.config import AgenticConfig
from mokume.agentic.profiler import profile_data
from mokume.agentic.proposer import (
    _dedup_configs,
    _format_seen_block,
    _items_to_configs,
    propose_configs,
)
from mokume.agentic.state import CandidateConfig


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
        for m in ("qrilc", "minprob", "mindet", "seqknn", "impseq")
    ]
    configs = _items_to_configs(items)
    assert [c.imputation for c in configs] == [
        "qrilc",
        "minprob",
        "mindet",
        "seqknn",
        "impseq",
    ]


def test_candidate_signature_stable_and_distinct():
    """signature() must be stable and discriminate config differences."""
    a = CandidateConfig(name="a", de_method="limma", normalization="median")
    a_renamed = CandidateConfig(name="A2", de_method="limma", normalization="median")
    b = CandidateConfig(name="b", de_method="limma", normalization="none")
    assert a.signature() == a_renamed.signature()  # name is not part of sig
    assert a.signature() != b.signature()


def test_dedup_configs_drops_seen():
    """_dedup_configs removes configs whose signature is already known."""
    c1 = CandidateConfig(name="c1", de_method="limma")
    c2 = CandidateConfig(name="c2", de_method="rots")
    deduped = _dedup_configs([c1, c2], {c1.signature()})
    assert [c.name for c in deduped] == ["c2"]


def test_dedup_configs_passthrough_when_seen_empty():
    """Empty seen set must not filter anything."""
    c1 = CandidateConfig(name="c1", de_method="limma")
    assert _dedup_configs([c1], set()) == [c1]
    assert _dedup_configs([c1], None) == [c1]


def test_format_seen_block_first_round():
    """No seen signatures yet -> friendly placeholder."""
    assert "first round" in _format_seen_block(None)
    assert "first round" in _format_seen_block(set())


def test_format_seen_block_sorted():
    """Seen signatures are listed deterministically."""
    sigs = {"limma|bh|none|none|fc0.5|none|k2", "rots|bh|none|none|fc0.5|none|k2"}
    out = _format_seen_block(sigs)
    assert "limma" in out and "rots" in out
    # Sorted output: lima... before rots... (alphabetical)
    assert out.index("limma") < out.index("rots")


def test_propose_configs_no_llm_dedup(synthetic_protein_df, sample_to_condition):
    """Rule-based fallback also honours the seen-signatures filter."""
    profile = profile_data(synthetic_protein_df, sample_to_condition)
    config = AgenticConfig(use_llm=False)
    first = propose_configs(profile, config)
    assert first, "first round should produce candidates"
    seen = {c.signature() for c in first}
    second = propose_configs(profile, config, seen=seen)
    # Either zero (full duplicates) or all-new signatures
    assert all(c.signature() not in seen for c in second)
