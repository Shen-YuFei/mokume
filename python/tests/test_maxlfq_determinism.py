"""MaxLFQ must not depend on the order its input rows happen to arrive in.

Regression tests for bigbio/mokume#58: the reference peptide was chosen with
``np.argmax``, which breaks ties by row position, so an unordered upstream read
made protein quantities change between identical runs (observed: >2 log2 units on
one protein, and the differential-expression call count moving between runs).
"""

import numpy as np
import pandas as pd
import pytest

from mokume.quantification.maxlfq import (
    _maxlfq_solve_protein,
    _process_protein,
    _select_reference_peptide,
)


def _matrix():
    # 4 peptides x 4 samples; several share the same number of valid values so the
    # tie-break is exercised, and one has a missing value.
    return np.array(
        [
            [100.0, 200.0, 400.0, 800.0],
            [ 50.0, 100.0, 200.0, 400.0],
            [ 10.0,  20.0,  40.0, np.nan],
            [400.0, 800.0, 1600.0, 3200.0],
        ]
    )


@pytest.mark.parametrize("seed", range(8))
def test_solve_is_invariant_to_row_permutation(seed):
    """Permuting peptide rows must not change the protein quantities."""
    m = _matrix()
    base = _maxlfq_solve_protein(m)
    rng = np.random.default_rng(seed)
    perm = rng.permutation(m.shape[0])
    permuted = _maxlfq_solve_protein(m[perm, :])
    np.testing.assert_allclose(permuted, base, rtol=0, atol=0)


def test_reference_selection_is_invariant_to_row_permutation():
    """The chosen reference must be the same PEPTIDE regardless of its position."""
    m = _matrix()
    with np.errstate(divide="ignore", invalid="ignore"):
        log_m = np.log2(m)
    counts = np.sum(~np.isnan(log_m), axis=1)
    chosen = log_m[_select_reference_peptide(log_m, counts)]

    for seed in range(8):
        perm = np.random.default_rng(seed).permutation(m.shape[0])
        pm = log_m[perm, :]
        pc = np.sum(~np.isnan(pm), axis=1)
        np.testing.assert_allclose(pm[_select_reference_peptide(pm, pc)], chosen)


def test_reference_prefers_most_measured_then_most_intense():
    log_m = np.array(
        [
            [1.0, 1.0, np.nan],   # fewer measurements - never chosen
            [2.0, 2.0, 2.0],      # full, total 6
            [3.0, 3.0, 3.0],      # full, total 9 -> the reference
        ]
    )
    counts = np.sum(~np.isnan(log_m), axis=1)
    assert _select_reference_peptide(log_m, counts) == 2


def test_identical_rows_are_interchangeable():
    """Fully tied rows are numerically identical, so any choice is equivalent."""
    log_m = np.array([[2.0, 2.0], [2.0, 2.0]])
    counts = np.sum(~np.isnan(log_m), axis=1)
    idx = _select_reference_peptide(log_m, counts)
    assert idx in (0, 1)
    np.testing.assert_allclose(log_m[idx], [2.0, 2.0])


@pytest.mark.parametrize("seed", range(5))
def test_process_protein_is_invariant_to_dataframe_row_order(seed):
    """End-to-end: shuffling the input frame must not change the output.

    ``Series.unique()`` returns values in order of appearance, so without a canonical
    sort the pivot - and therefore the matrix handed to MaxLFQ - inherits the incoming
    row order.
    """
    rows = []
    for pep, scale in (("PEPTIDEA", 1.0), ("PEPTIDEB", 0.5), ("PEPTIDEC", 4.0)):
        for i, sample in enumerate(("S1", "S2", "S3")):
            rows.append({"protein": "P1", "peptide": pep, "sample": sample,
                         "intensity": 100.0 * scale * (i + 1)})
    df = pd.DataFrame(rows)
    samples = ["S1", "S2", "S3"]

    def run(frame):
        out = _process_protein("P1", frame, "peptide", "intensity", "sample", samples, 2)
        return {r["sample"]: r["intensity"] for r in out}

    base = run(df)
    shuffled = run(df.sample(frac=1.0, random_state=seed).reset_index(drop=True))
    assert set(base) == set(shuffled)
    for s in base:
        assert shuffled[s] == pytest.approx(base[s], rel=0, abs=0)
