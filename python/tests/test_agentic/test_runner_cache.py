"""Tests for the agentic preprocessing cache."""

import numpy as np
import pandas as pd

from mokume.agentic.runner import PreprocessCache


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
