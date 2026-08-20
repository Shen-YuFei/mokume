"""Tests for mokume.analysis.adaptive_fdr (Storey/BKY dispatch, pi0, reliability gate).

All p-value samples are drawn from seeded generators so every assertion below is
a fixed number, not a distributional guess.
"""

import numpy as np
import pytest

from mokume.analysis.adaptive_fdr import (
    _estimate_pi0_detail,
    adjust_pvalues,
    estimate_pi0,
    pi0_lower_bound,
    pi0_reliability,
    qvalues,
)

# --------------------------------------------------------------------------
# Fixtures: p-value populations with known shapes.
# --------------------------------------------------------------------------


def _two_group_pvalues(n_null=900, n_alt=300, seed=0):
    """Anti-conservative p-values: uniform nulls plus a spike near 0."""
    rng = np.random.default_rng(seed)
    return np.concatenate([rng.uniform(0, 1, n_null), rng.beta(0.3, 8.0, n_alt)])


def _mid_mass_pvalues(n=2000, seed=11):
    """Beta(2,2): mass in the middle, so the histogram has no spike near 0.

    pi0 lands well below 1 and n is large, so ``anticonservative`` is the only
    reliability criterion that fails -- this isolates that one branch.
    """
    return np.random.default_rng(seed).beta(2.0, 2.0, n)


def _conservative_pvalues(n=2000, seed=3):
    """Beta(3,1): right-skewed. Trips both at_boundary and anticonservative."""
    return np.random.default_rng(seed).beta(3.0, 1.0, n)


# --------------------------------------------------------------------------
# 1. bh / bky / storey each dispatch, and the adaptive pair is looser than BH.
# --------------------------------------------------------------------------


def test_adjust_pvalues_dispatches_all_three_methods():
    """Each of bh/bky/storey runs its own branch and reports itself as applied."""
    p = _two_group_pvalues()
    for method in ("bh", "bky", "storey"):
        adjusted, used = adjust_pvalues(p, method=method)
        assert used == method, f"{method} unexpectedly fell back to {used}"
        assert adjusted.shape == p.shape
        assert np.all(np.isfinite(adjusted))
        assert np.all((adjusted >= 0) & (adjusted <= 1))


def test_storey_and_bky_reject_more_than_bh_when_pi0_below_one():
    """With pi0 < 1 the adaptive procedures recover the budget BH wastes.

    Fixed counts for this seed: bh=137 < bky=144 < storey=155.
    """
    p = _two_group_pvalues()
    counts = {}
    for method in ("bh", "bky", "storey"):
        adjusted, used = adjust_pvalues(p, method=method)
        assert used == method
        counts[method] = int((adjusted < 0.05).sum())

    assert counts["bh"] == 137
    assert counts["bky"] == 144
    assert counts["storey"] == 155
    assert counts["bky"] > counts["bh"]
    assert counts["storey"] > counts["bh"]


def test_unknown_method_falls_back_to_bh():
    """A method this module does not implement (e.g. ihw) routes to BH."""
    p = _two_group_pvalues()
    adjusted, used = adjust_pvalues(p, method="ihw")
    assert used == "bh"
    bh_adjusted, _ = adjust_pvalues(p, method="bh")
    np.testing.assert_allclose(adjusted, bh_adjusted)


def test_storey_qvalues_known_answer():
    """q(p_(i)) = min over the tail of pi0 * m * p_(i) / i, NaN-aligned.

    m=2, pi0=0.5: q = [0.5*2*0.001/1, 0.5*2*0.5/2] = [0.001, 0.25].
    """
    q = qvalues(np.array([0.001, 0.5, np.nan]), pi0=0.5)
    assert q[0] == pytest.approx(0.001)
    assert q[1] == pytest.approx(0.25)
    assert np.isnan(q[2])


def test_qvalues_are_monotone_in_pvalue():
    """Storey q-values never decrease as the p-value grows."""
    p = _two_group_pvalues()
    q = qvalues(p, pi0=0.8)
    order = np.argsort(p)
    assert np.all(np.diff(q[order]) >= -1e-12)


def test_adjust_pvalues_aligns_nonfinite_input_to_nan():
    """Adjusted output keeps the input shape, with NaN where the input was not finite."""
    p = _two_group_pvalues()
    p_with_gaps = p.copy()
    p_with_gaps[[5, 50, 500]] = np.nan
    adjusted, used = adjust_pvalues(p_with_gaps, method="storey")
    assert used == "storey"
    assert adjusted.shape == p_with_gaps.shape
    assert np.all(np.isnan(adjusted[[5, 50, 500]]))
    assert np.isfinite(np.delete(adjusted, [5, 50, 500])).all()


# --------------------------------------------------------------------------
# 2. Reliability gate: each criterion triggers the BH fallback on its own.
# --------------------------------------------------------------------------


def test_reliability_too_few_triggers_alone():
    """n < 100 makes pi0 unstable: too_few fires while the others stay clean."""
    p = _two_group_pvalues(n_null=60, n_alt=20, seed=1)
    rel = pi0_reliability(p)
    assert rel["n"] == 80
    assert rel["too_few"] is True
    assert rel["at_boundary"] is False
    assert rel["anticonservative"] is True
    assert rel["reliable"] is False


def test_too_few_boundary_is_exactly_100():
    """99 hypotheses are 'too few', 100 are not."""
    rng = np.random.default_rng(2)
    p = rng.uniform(0, 1, 100)
    assert pi0_reliability(p)["too_few"] is False
    assert pi0_reliability(p[:99])["too_few"] is True


def test_reliability_at_boundary_triggers_alone():
    """pi0 ~ 1 is untrustworthy even on a clearly anti-conservative sample.

    pi0 is supplied so only the boundary criterion can fail.
    """
    p = _two_group_pvalues()
    rel = pi0_reliability(p, pi0=1.0)
    assert rel["at_boundary"] is True
    assert rel["too_few"] is False
    assert rel["anticonservative"] is True
    assert rel["reliable"] is False


def test_reliability_anticonservative_triggers_alone():
    """A histogram with no spike near 0 fails only the anticonservative check."""
    p = _mid_mass_pvalues()
    rel = pi0_reliability(p)
    assert rel["anticonservative"] is False
    assert rel["at_boundary"] is False
    assert rel["too_few"] is False
    assert rel["reliable"] is False


def test_reliable_when_all_criteria_pass():
    """A well-behaved two-group sample is trusted, so no fallback happens."""
    rel = pi0_reliability(_two_group_pvalues())
    assert rel["reliable"] is True
    assert rel["pi0"] == pytest.approx(0.7738, abs=1e-3)


@pytest.mark.parametrize("method", ["storey", "bky"])
def test_adjust_pvalues_falls_back_to_bh_when_too_few(method):
    """too_few alone forces method_used back to 'bh'."""
    p = _two_group_pvalues(n_null=60, n_alt=20, seed=1)
    adjusted, used = adjust_pvalues(p, method=method)
    assert used == "bh"
    bh_adjusted, _ = adjust_pvalues(p, method="bh")
    np.testing.assert_allclose(adjusted, bh_adjusted)


@pytest.mark.parametrize("method", ["storey", "bky"])
def test_adjust_pvalues_falls_back_to_bh_when_not_anticonservative(method):
    """A mid-mass histogram alone forces the fallback (n and pi0 are both fine)."""
    p = _mid_mass_pvalues()
    assert pi0_reliability(p)["anticonservative"] is False
    adjusted, used = adjust_pvalues(p, method=method)
    assert used == "bh"
    bh_adjusted, _ = adjust_pvalues(p, method="bh")
    np.testing.assert_allclose(adjusted, bh_adjusted)


@pytest.mark.parametrize("method", ["storey", "bky"])
def test_adjust_pvalues_falls_back_to_bh_at_pi0_boundary(method):
    """pi0 pinned at 1 forces the fallback.

    This right-skewed sample trips at_boundary *and* anticonservative together
    (a conservative histogram is what drives pi0 to 1), so it pins the boundary
    branch of the fallback rather than isolating it; the isolated boundary check
    is :func:`test_reliability_at_boundary_triggers_alone`.
    """
    p = _conservative_pvalues()
    rel = pi0_reliability(p)
    assert rel["at_boundary"] is True
    assert rel["pi0"] == pytest.approx(1.0)
    adjusted, used = adjust_pvalues(p, method=method)
    assert used == "bh"
    bh_adjusted, _ = adjust_pvalues(p, method="bh")
    np.testing.assert_allclose(adjusted, bh_adjusted)


def test_adjust_pvalues_all_nonfinite_returns_nan_without_fallback():
    """An all-NaN input yields all-NaN output and keeps the requested method name."""
    adjusted, used = adjust_pvalues(np.full(5, np.nan), method="storey")
    assert used == "storey"
    assert np.all(np.isnan(adjusted))


# --------------------------------------------------------------------------
# 3. pi0 conservative lower bound.
# --------------------------------------------------------------------------


def test_pi0_bound_never_lowers_the_estimate():
    """conservative_bound=True can only raise pi0, never lower it.

    Checked across 200 seeded two-group datasets: the bound is applied as
    max(pi0_hat, bound), so the bounded estimate dominates pointwise.
    """
    for seed in range(200):
        rng = np.random.default_rng(seed + 10_000)
        p = np.concatenate([rng.uniform(0, 1, 800), rng.beta(0.3, 8.0, 200)])
        bounded = estimate_pi0(p, conservative_bound=True)
        unbounded = estimate_pi0(p, conservative_bound=False)
        assert bounded >= unbounded - 1e-12, f"bound lowered pi0 at seed {seed}"


def test_pi0_lower_bound_is_the_curve_minimum():
    """The bound is min over the lambda grid of #{p >= lam} / ((1 - lam) m)."""
    p = _two_group_pvalues()
    lambdas = np.arange(0.05, 0.96, 0.05)
    expected = min(float(np.mean(p >= lam) / (1.0 - lam)) for lam in lambdas)
    assert pi0_lower_bound(p) == pytest.approx(expected)
    assert 0.0 < pi0_lower_bound(p) <= 1.0


def test_pi0_detail_reports_raw_bound_and_whether_it_applied():
    """_estimate_pi0_detail exposes pi0_raw / pi0_lower_bound / bound_applied."""
    p = _two_group_pvalues()
    detail = _estimate_pi0_detail(p)
    assert set(detail) == {"pi0", "pi0_raw", "pi0_lower_bound", "bound_applied"}
    assert detail["bound_applied"] is False  # bound does not bind on this seed
    assert detail["pi0"] == pytest.approx(detail["pi0_raw"])
    assert detail["pi0_lower_bound"] < detail["pi0_raw"]


def test_pi0_bound_binds_when_the_smoother_leaves_the_curve():
    """When the spline extrapolates below every observed point, the bound wins.

    This seed is one of the ~0.3% of two-group draws (6/2000 searched) where the
    smoother's lambda -> 1 extrapolation lands under the curve minimum. pi0 is
    then raised to the bound, which is exactly the documented behaviour.
    """
    rng = np.random.default_rng(793)
    n_null = int(rng.integers(100, 600))
    n_alt = int(rng.integers(100, 900))
    a = float(rng.uniform(0.05, 0.5))
    b = float(rng.uniform(5.0, 40.0))
    p = np.concatenate([rng.uniform(0, 1, n_null), rng.beta(a, b, n_alt)])

    detail = _estimate_pi0_detail(p)
    assert detail["bound_applied"] is True
    assert detail["pi0_raw"] == pytest.approx(0.2446, abs=1e-3)
    assert detail["pi0_lower_bound"] == pytest.approx(0.2493, abs=1e-3)
    assert detail["pi0"] == pytest.approx(detail["pi0_lower_bound"])
    assert detail["pi0"] > detail["pi0_raw"]
    # Turning the bound off returns the unbounded (more aggressive) estimate.
    assert estimate_pi0(p, conservative_bound=False) == pytest.approx(detail["pi0_raw"])


def test_pi0_reliability_exposes_bound_fields():
    """pi0_reliability surfaces the bound diagnostics alongside the gate flags."""
    rel = pi0_reliability(_two_group_pvalues())
    for key in ("pi0_raw", "pi0_lower_bound", "bound_applied"):
        assert key in rel
    assert rel["pi0"] == pytest.approx(rel["pi0_raw"])
    assert 0.0 < rel["pi0_lower_bound"] <= 1.0


def test_pi0_reliability_with_supplied_pi0_does_not_recompute_the_bound():
    """A caller-supplied pi0 is taken as given: the bound is reported as NaN."""
    rel = pi0_reliability(_two_group_pvalues(), pi0=0.6)
    assert rel["pi0"] == pytest.approx(0.6)
    assert rel["pi0_raw"] == pytest.approx(0.6)
    assert np.isnan(rel["pi0_lower_bound"])
    assert rel["bound_applied"] is False


# --------------------------------------------------------------------------
# pi0 estimator basics.
# --------------------------------------------------------------------------


def test_estimate_pi0_recovers_the_true_null_fraction():
    """With 900 nulls of 1200 the true pi0 is 0.75; the estimate lands near it."""
    p = _two_group_pvalues(n_null=900, n_alt=300)
    assert estimate_pi0(p) == pytest.approx(0.75, abs=0.05)


def test_estimate_pi0_is_one_for_pure_nulls():
    """An all-uniform sample carries no signal, so pi0 sits at the 1.0 boundary."""
    p = np.random.default_rng(3).uniform(0, 1, 3000)
    assert estimate_pi0(p) == pytest.approx(1.0, abs=0.05)


@pytest.mark.parametrize("method", ["smoother", "bootstrap", "fixed"])
def test_estimate_pi0_methods_agree_within_tolerance(method):
    """All three pi0 estimators land near the true 0.75 on the same sample."""
    p = _two_group_pvalues()
    assert estimate_pi0(p, method=method) == pytest.approx(0.75, abs=0.08)


def test_estimate_pi0_rejects_unknown_method():
    """An unknown pi0 method is a programming error, not a silent fallback."""
    with pytest.raises(ValueError, match="unknown pi0 method"):
        estimate_pi0(_two_group_pvalues(), method="nope")


def test_estimate_pi0_empty_input_is_one():
    """No hypotheses means nothing to reject: pi0 defaults to the safe 1.0."""
    assert estimate_pi0(np.array([])) == 1.0


def test_pvalues_outside_unit_interval_are_rejected():
    """p-values must lie in [0, 1]; anything else is a caller bug."""
    with pytest.raises(ValueError, match=r"p-values must lie in \[0, 1\]"):
        estimate_pi0(np.array([0.1, 1.5, 0.3]))
