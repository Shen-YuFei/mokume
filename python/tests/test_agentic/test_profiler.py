"""Tests for mokume.agentic.profiler."""

from mokume.agentic.profiler import profile_data


def test_profile_basic(synthetic_protein_df, sample_to_condition):
    """Profile returns correct shape and type info."""
    profile = profile_data(synthetic_protein_df, sample_to_condition)
    assert profile.n_proteins == 30
    assert profile.n_samples == 6
    assert profile.n_conditions == 2
    assert 0 < profile.missing_rate < 0.20
    assert profile.is_log_transformed is True
    assert profile.data_type == "LFQ"


def test_profile_to_dict(synthetic_protein_df, sample_to_condition):
    """to_dict() returns JSON-serializable dict."""
    profile = profile_data(synthetic_protein_df, sample_to_condition)
    d = profile.to_dict()
    assert isinstance(d, dict)
    assert "n_proteins" in d
    assert isinstance(d["intensity_range"], list)
