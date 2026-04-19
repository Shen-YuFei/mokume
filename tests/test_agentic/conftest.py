"""Shared fixtures for agentic tests."""

import numpy as np
import pandas as pd
import pytest


@pytest.fixture()
def synthetic_protein_df():
    """Create a synthetic protein matrix (30 proteins x 6 samples)."""
    np.random.seed(42)
    n_proteins = 30
    n_samples = 6
    proteins = [f"P{i:05d}" for i in range(n_proteins)]
    samples = [f"S{i}" for i in range(1, n_samples + 1)]

    # Log2 intensities ~ N(20, 2)
    data = np.random.normal(20, 2, (n_proteins, n_samples))

    # Spike first 5 proteins up in condition A (S1-S3)
    data[:5, :3] += 3.0  # ~8x fold change

    # Introduce 10% missing
    mask = np.random.random((n_proteins, n_samples)) < 0.10
    data[mask] = np.nan

    df = pd.DataFrame(data, index=proteins, columns=samples)
    df.insert(0, "protein", proteins)
    return df


@pytest.fixture()
def sample_to_condition():
    """Map samples to 2 conditions."""
    return {
        "S1": "A",
        "S2": "A",
        "S3": "A",
        "S4": "B",
        "S5": "B",
        "S6": "B",
    }


@pytest.fixture()
def ground_truth_proteins():
    """First 5 proteins are spiked-in."""
    return {f"P{i:05d}" for i in range(5)}
