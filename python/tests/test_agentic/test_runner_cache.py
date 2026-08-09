"""Tests for the agentic preprocessing cache."""

import numpy as np
import pandas as pd
import pytest

from mokume.agentic.runner import PreprocessCache, run_experiment
from mokume.agentic.state import CandidateConfig


def _toy_matrix() -> pd.DataFrame:
    """Tiny protein matrix in log2 space."""
    np.random.seed(0)
    data = np.random.normal(20, 2, (8, 4))
    df = pd.DataFrame(data, columns=[f"S{i}" for i in range(1, 5)])
    df.insert(0, "protein", [f"P{i}" for i in range(8)])
    return df


def test_cache_hits_on_repeated_pair():
    """Same (norm, imp) pair the second time returns a cached copy."""
    cache = PreprocessCache()
    df = _toy_matrix()
    a = cache.get_or_compute("none", "none", df)
    b = cache.get_or_compute("none", "none", df)
    assert cache.misses == 1
    assert cache.hits == 1
    pd.testing.assert_frame_equal(a, b)


def test_cache_misses_on_different_pair():
    """Different (norm, imp) pairs are stored independently."""
    cache = PreprocessCache()
    df = _toy_matrix()
    cache.get_or_compute("none", "none", df)
    cache.get_or_compute("median", "none", df)
    cache.get_or_compute("none", "mindet", df)
    assert cache.misses == 3
    assert cache.hits == 0
    assert cache.stats()["unique_combos"] == 3


def test_cache_returns_copy_not_reference():
    """Mutating the returned frame must not poison the cache."""
    cache = PreprocessCache()
    df = _toy_matrix()
    first = cache.get_or_compute("none", "none", df)
    first.iloc[0, 1] = -99999
    second = cache.get_or_compute("none", "none", df)
    assert second.iloc[0, 1] != -99999


def test_invalid_ensemble_fails_before_preprocessing(monkeypatch):
    def unexpected_preprocessing(*args, **kwargs):
        raise AssertionError("preprocessing must not run")

    monkeypatch.setattr(PreprocessCache, "get_or_compute", unexpected_preprocessing)
    config = CandidateConfig(
        name="invalid_ensemble",
        de_method="ensemble",
        ensemble="limma,LIMMA",
        ensemble_k=1,
    )

    with pytest.raises(ValueError, match="duplicate ensemble"):
        run_experiment(
            config,
            _toy_matrix(),
            {f"S{i}": "A" if i < 3 else "B" for i in range(1, 5)},
            ("A", "B"),
            cache=PreprocessCache(),
        )


def test_ensemble_members_require_ensemble_de_method(monkeypatch):
    def unexpected_preprocessing(*args, **kwargs):
        raise AssertionError("preprocessing must not run")

    monkeypatch.setattr(PreprocessCache, "get_or_compute", unexpected_preprocessing)
    config = CandidateConfig(
        name="inconsistent_ensemble",
        de_method="limma",
        ensemble="limma,deqms",
        ensemble_k=1,
    )

    with pytest.raises(ValueError, match="only apply"):
        run_experiment(
            config,
            _toy_matrix(),
            {f"S{i}": "A" if i < 3 else "B" for i in range(1, 5)},
            ("A", "B"),
            cache=PreprocessCache(),
        )
