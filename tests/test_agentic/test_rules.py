"""Tests for mokume.agentic.rules."""

from mokume.agentic.profiler import profile_data
from mokume.agentic.rules import load_benchmarks, load_prompts, rule_propose


def test_rule_propose(synthetic_protein_df, sample_to_condition):
    """Rule engine generates candidates for LFQ data."""
    profile = profile_data(synthetic_protein_df, sample_to_condition)
    configs = rule_propose(profile)
    assert len(configs) > 0
    assert all(
        c.de_method in ("deqms", "limma", "rots", "proda", "msstats") for c in configs
    )
    assert all(c.fdr_method in ("bh", "ihw") for c in configs)


def test_load_prompts():
    """Prompt templates load correctly."""
    prompts = load_prompts()
    assert "proposal_system" in prompts
    assert "proposal_user" in prompts
    assert "reflection_system" in prompts
    assert "reflection_user" in prompts
    assert "report_system" in prompts
    assert "report_user" in prompts


def test_load_benchmarks():
    """Benchmark data loads correctly."""
    bm = load_benchmarks()
    assert "PXD001819" in bm
