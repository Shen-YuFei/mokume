"""Tests for mokume.analysis.effect_size_gate and the DE "auto" gate wiring.

The point of the estimator is that it must be free to return a gate *below* the
conventional fixed 0.5 (and below the old 0.1 floor) when the fold-change null
is narrow, as it is for compressed isobaric ratios. These tests pin that.
"""

import numpy as np
import pandas as pd
import pytest

from mokume.analysis.differential_expression import DifferentialExpression
from mokume.analysis.effect_size_gate import estimate_effect_size_gate

SAMPLE_TO_CONDITION = {
    "S1": "A",
    "S2": "A",
    "S3": "A",
    "S4": "B",
    "S5": "B",
    "S6": "B",
}


def _narrow_null_log2fc(n_null=900, n_signal=100, null_sd=0.01, seed=7):
    """Fold changes with a very tight null (sd=0.01) and a separable signal."""
    rng = np.random.default_rng(seed)
    return np.concatenate(
        [rng.normal(0.0, null_sd, n_null), rng.normal(0.3, 0.05, n_signal)]
    )


def _compressed_matrix(n_proteins=300, n_spike=40, seed=0):
    """Log2 matrix with an isobaric-like compressed effect (0.35) on a tight null."""
    rng = np.random.default_rng(seed)
    data = rng.normal(20.0, 0.05, (n_proteins, 6))
    data[:n_spike, :3] += 0.35
    df = pd.DataFrame(data, columns=list(SAMPLE_TO_CONDITION))
    df.insert(0, "protein", [f"P{i:04d}" for i in range(n_proteins)])
    return df


# --------------------------------------------------------------------------
# 4. The estimator goes below the removed 0.1 floor on a narrow null.
# --------------------------------------------------------------------------


def test_gate_goes_below_the_removed_floor_on_a_narrow_null():
    """A null of sd=0.01 yields a gate well under 0.1, proving no floor remains.

    A floor of 0.1 (or the fixed 0.5 fallback) would silently override the null
    width the data reports -- the exact miscalibration this module removes.
    """
    gate = estimate_effect_size_gate(_narrow_null_log2fc())
    assert gate == pytest.approx(0.0456, abs=1e-3)
    assert gate < 0.1, "gate was clamped at the removed floor"
    assert gate != 0.5, "gate silently fell back instead of being estimated"
    assert gate > 0


def test_gate_tracks_the_null_width():
    """A wider null must produce a larger gate: the gate follows the data."""
    narrow = estimate_effect_size_gate(_narrow_null_log2fc(null_sd=0.01))
    wide = estimate_effect_size_gate(_narrow_null_log2fc(null_sd=0.20))
    assert wide > narrow


@pytest.mark.parametrize("method", ["mixture", "null_quantile"])
def test_both_estimators_return_a_subfloor_gate_on_a_narrow_null(method):
    """Neither estimator imposes a lower clamp."""
    gate = estimate_effect_size_gate(_narrow_null_log2fc(), method=method)
    assert 0 < gate < 0.1


# --------------------------------------------------------------------------
# 4b. Degenerate inputs fall back (and only then).
# --------------------------------------------------------------------------


def test_gate_falls_back_when_too_few_points():
    """Fewer than 50 finite values: the null is not estimable at all."""
    assert estimate_effect_size_gate(np.array([0.1, 0.2, 0.3])) == 0.5


def test_gate_falls_back_just_under_the_50_point_threshold():
    """49 points fall back; 50 are enough to estimate."""
    x = _narrow_null_log2fc()
    assert estimate_effect_size_gate(x[:49]) == 0.5
    assert estimate_effect_size_gate(x[:50]) != 0.5


def test_gate_falls_back_on_constant_input():
    """A constant fold change has no separable signal component."""
    assert estimate_effect_size_gate(np.zeros(500)) == 0.5


def test_gate_falls_back_on_all_nan_input():
    """No finite values at all: fall back rather than emit NaN."""
    assert estimate_effect_size_gate(np.full(500, np.nan)) == 0.5


def test_gate_fallback_value_is_configurable():
    """The fallback is a caller-supplied default, not a hard-coded 0.5."""
    assert estimate_effect_size_gate(np.zeros(500), fallback=0.25) == 0.25
    assert estimate_effect_size_gate(np.array([1.0]), fallback=0.25) == 0.25


def test_gate_ignores_nonfinite_values():
    """NaN/inf are dropped, so padding a clean sample with them changes nothing."""
    x = _narrow_null_log2fc()
    padded = np.concatenate([x, [np.nan, np.inf, -np.inf]])
    assert estimate_effect_size_gate(padded) == pytest.approx(
        estimate_effect_size_gate(x)
    )


def test_gate_rejects_unknown_method():
    """An unknown gate method is a programming error, not a silent fallback."""
    with pytest.raises(ValueError, match="unknown gate method"):
        estimate_effect_size_gate(_narrow_null_log2fc(), method="nope")


def test_gate_warns_only_when_it_falls_back(caplog):
    """A data-derived gate logs no warning; a fallback always does."""
    with caplog.at_level("WARNING", logger="mokume.analysis.effect_size_gate"):
        estimate_effect_size_gate(_narrow_null_log2fc())
    assert caplog.text == ""

    with caplog.at_level("WARNING", logger="mokume.analysis.effect_size_gate"):
        estimate_effect_size_gate(np.zeros(500))
    assert "fallback" in caplog.text


# --------------------------------------------------------------------------
# 5. DifferentialExpression(log2fc_threshold="auto") end to end.
# --------------------------------------------------------------------------


def test_de_auto_gate_runs_end_to_end_and_is_data_derived():
    """ "auto" resolves to the estimated gate, not the hard-coded 0.5."""
    de = DifferentialExpression(method="limma", log2fc_threshold="auto", skip_log2=True)
    result = de.run(_compressed_matrix(), SAMPLE_TO_CONDITION, ("A", "B"))

    gate = estimate_effect_size_gate(result["log2FC"].to_numpy(dtype=float))
    assert gate == pytest.approx(0.165, abs=1e-2)
    assert gate < 0.5, "auto gate must adapt below the fixed default"
    assert "significance" in result.columns


def test_de_auto_gate_recovers_compressed_effects_the_fixed_gate_misses():
    """On compressed ratios the fixed 0.5 gate calls nothing; "auto" calls the spike.

    The 40 spiked proteins carry a true log2FC of 0.35 -- real, but below 0.5.
    This is the isobaric miscalibration the auto gate exists to fix.
    """
    matrix = _compressed_matrix()

    auto = DifferentialExpression(
        method="limma", log2fc_threshold="auto", skip_log2=True
    ).run(matrix, SAMPLE_TO_CONDITION, ("A", "B"))
    fixed = DifferentialExpression(
        method="limma", log2fc_threshold=0.5, skip_log2=True
    ).run(matrix, SAMPLE_TO_CONDITION, ("A", "B"))

    assert int((auto["significance"] == "UP").sum()) == 40
    assert int((fixed["significance"] == "UP").sum()) == 0
    # The spiked proteins are exactly the ones "auto" recovers.
    called = set(auto.loc[auto["significance"] == "UP", "ProteinName"])
    assert called == {f"P{i:04d}" for i in range(40)}


def test_de_fixed_gate_still_accepts_a_plain_float():
    """A numeric threshold is applied verbatim: no estimation is attempted.

    Every call must clear the requested 0.1 gate, and the 40 spiked proteins
    must all clear it. The call set is a superset of the spike: at a gate this
    loose one null protein (log2FC=0.133) also passes, which is a real false
    positive of the 0.1 operating point rather than a wiring defect -- so this
    asserts containment, not an exact count.
    """
    matrix = _compressed_matrix()
    result = DifferentialExpression(
        method="limma", log2fc_threshold=0.1, skip_log2=True
    ).run(matrix, SAMPLE_TO_CONDITION, ("A", "B"))

    up = result[result["significance"] == "UP"]
    assert (up["log2FC"] > 0.1).all()
    spiked = {f"P{i:04d}" for i in range(40)}
    assert spiked <= set(up["ProteinName"])
    # The gate is honoured as given, so it admits strictly more than the 0.165
    # the estimator would have chosen on this same matrix.
    auto = DifferentialExpression(
        method="limma", log2fc_threshold="auto", skip_log2=True
    ).run(matrix, SAMPLE_TO_CONDITION, ("A", "B"))
    assert len(up) > int((auto["significance"] == "UP").sum())
