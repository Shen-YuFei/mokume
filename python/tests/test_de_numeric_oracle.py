"""Cell-exact numeric oracle regression tests for the Python DE methods.

Motivation
----------
Every other DE test in this suite is shape/smoke only (does it return the
expected columns, does it flag *some* protein DE). None of them pins a numeric
cell to a reference value, so a *silent statistical regression* -- a flipped
sign, a mislabelled column -- passes CI unnoticed. A reviewer recently found
exactly two such bugs that this gap let through:

  * ROTS ``log2FC`` had its sign flipped (``mean(B) - mean(A)`` instead of the
    correct ``mean(A) - mean(B)``), and
  * limma's ``AveExpr`` column actually carried the ``log2FC`` value, not the
    per-protein row mean of the log2 matrix.

These tests pin the affected quantities to independently-derived reference
values so either regression fails CI.

Reference-value provenance
--------------------------
The fixture is the same 6-protein / A1..A3,B1..B3 raw-intensity matrix used by
the trusted Rust golden oracle
``rust/crates/mokume-pipeline/src/de.rs`` (constants ``RAW`` / ``EXPECTED`` /
``EXPECTED_AVE_EXPR`` / ``EXPECTED_ROTS_LOG2FC``), recreated here in Python so
the numbers line up and are auditable across languages.

Two provenance tiers are used, deliberately:

  1. **Independent derivation** (preferred). ``log2FC = mean(A) - mean(B)`` on
     log2 data, ``AveExpr = row mean of the log2 matrix``, and the ROTS sign
     are all standard and hand-computable. They are computed here directly from
     the fixture with numpy -- NOT by calling the function under test -- so the
     check is not circular. (These independently-derived numbers also match the
     Rust ``de.rs`` constants to floating-point noise, which is a nice
     cross-check but is not what the assertions rely on.)

  2. **Pinned to the Rust golden oracle** (for eBayes-moderated quantities that
     are not trivially hand-derivable). The moderated ``t_stat`` / ``pvalue``
     come from ``de.rs::EXPECTED`` (the cross-validated Rust golden values that
     yasset's review implicitly accepted -- he flagged the ``AveExpr`` *label*,
     not the moderated numbers). These are clearly marked below.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from mokume.analysis.limma import run_limma
from mokume.analysis.rots import run_rots

# ---------------------------------------------------------------------------
# Fixture: the 6-protein / 3-vs-3 RAW intensity matrix from the Rust de.rs
# oracle (rust/crates/mokume-pipeline/src/de.rs, constant `RAW`). Proteins are
# rows, samples are columns; A1..A3 are group A, B1..B3 are group B. Every cell
# is strictly positive, so log2 never produces a NaN and every protein is
# testable (>= 2 finite values per group).
# ---------------------------------------------------------------------------

PROTEINS = ["P1", "P2", "P3", "P4", "P5", "P6"]
SAMPLES_A = ["A1", "A2", "A3"]
SAMPLES_B = ["B1", "B2", "B3"]
SAMPLES = SAMPLES_A + SAMPLES_B

RAW = np.array(
    [
        [1000.0, 1149.0, 870.0, 4000.0, 4290.0, 3732.0],
        [32000.0, 34000.0, 30000.0, 8000.0, 9200.0, 6800.0],
        [256.0, 270.0, 244.0, 257.0, 248.0, 256.0],
        [1048576.0, 4194304.0, 262144.0, 2097152.0, 524288.0, 8388608.0],
        [32.0, 32.2, 31.8, 64.0, 64.4, 63.6],
        [512.0, 724.0, 362.0, 480.0, 660.0, 350.0],
    ]
)


@pytest.fixture(scope="module")
def log2_df() -> pd.DataFrame:
    """log2 protein x sample matrix, the input the DE methods expect."""
    return pd.DataFrame(np.log2(RAW), index=PROTEINS, columns=SAMPLES)


# ---------------------------------------------------------------------------
# Independently-derived reference values (tier 1). Computed straight from the
# fixture with numpy -- NOT from the function under test.
# ---------------------------------------------------------------------------

_LOG2 = np.log2(RAW)
# log2FC = mean(group A) - mean(group B) on the log2 matrix.
REF_LOG2FC = dict(zip(PROTEINS, _LOG2[:, :3].mean(axis=1) - _LOG2[:, 3:].mean(axis=1)))
# AveExpr = row mean of the log2 matrix over ALL samples (limma's Amean).
REF_AVE_EXPR = dict(zip(PROTEINS, _LOG2.mean(axis=1)))
# Per-group means, used to check the sign of log2FC against the data direction.
REF_MEAN_A = dict(zip(PROTEINS, _LOG2[:, :3].mean(axis=1)))
REF_MEAN_B = dict(zip(PROTEINS, _LOG2[:, 3:].mean(axis=1)))


# ---------------------------------------------------------------------------
# eBayes-MODERATED reference values (tier 2). Pinned to the Rust golden oracle
# `de.rs::EXPECTED` (protein, log2FC, pvalue, adj_pvalue, t_stat, significance)
# -- the cross-validated moderated numbers the review accepted. These are NOT
# hand-derivable, so they are pinned, not independently recomputed. The
# moderated t deliberately differs from an ordinary two-sample t (e.g. P5:
# moderated -58.99 vs ordinary -135.83), so this pins the eBayes shrinkage path,
# not a plain t-test.
# ---------------------------------------------------------------------------

REF_MOD_T = {
    "P1": -16.39236762936801,
    "P2": 15.653604075943404,
    "P5": -58.99007257281155,
}
REF_MOD_P = {
    "P1": 3.0352542907348282e-05,
    "P2": 3.737705510611045e-05,
    "P5": 8.847368379958774e-08,
}


# ---------------------------------------------------------------------------
# Sanity: the independent tier-1 values agree with the Rust de.rs constants.
# This is a cross-language audit of the fixture, not the load-bearing check --
# the assertions below stand on the numpy derivation regardless.
# ---------------------------------------------------------------------------


def test_fixture_matches_rust_oracle_constants():
    """Independent numpy log2FC/AveExpr == de.rs EXPECTED / EXPECTED_AVE_EXPR."""
    # From rust/crates/mokume-pipeline/src/de.rs (EXPECTED[*].1 and
    # EXPECTED_AVE_EXPR). Auditable cross-check that the Python fixture is the
    # same matrix as the Rust golden oracle.
    rust_log2fc = {
        "P1": -2.0004868432854863,
        "P2": 2.0090616097754097,
        "P3": 0.015910691677655464,
        "P4": -1.0,
        "P5": -1.0,
        "P6": 0.09175595082658106,
    }
    rust_ave_expr = {
        "P1": 10.96584974099084,
        "P2": 13.959371292060668,
        "P3": 7.994562299032412,
        "P4": 20.5,
        "P5": 5.499981214541417,
        "P6": 8.954019282642177,
    }
    for protein in PROTEINS:
        assert REF_LOG2FC[protein] == pytest.approx(rust_log2fc[protein], abs=1e-12)
        assert REF_AVE_EXPR[protein] == pytest.approx(rust_ave_expr[protein], abs=1e-12)


# ---------------------------------------------------------------------------
# limma: log2FC AND AveExpr -- the exact regression yasset caught.
# ---------------------------------------------------------------------------


def test_limma_log2fc_is_mean_a_minus_mean_b(log2_df):
    """limma log2FC == independently-derived mean(A) - mean(B) on log2 data."""
    res = run_limma(log2_df, SAMPLES_A, SAMPLES_B, "groupA", "groupB")
    res = res.set_index("ProteinName")
    for protein in PROTEINS:
        assert res.loc[protein, "log2FC"] == pytest.approx(
            REF_LOG2FC[protein], rel=1e-9, abs=1e-12
        )


def test_limma_avexpr_is_rowmean_not_log2fc(log2_df):
    """limma AveExpr == row mean of the log2 matrix, and is NOT the log2FC.

    Pins the exact bug the reviewer found: AveExpr had been the log2FC value.
    The two assertions are complementary --
      * ``AveExpr == row mean``  (positive statement of what it must be), and
      * ``AveExpr != log2FC``    (the specific regression that must not recur).
    For this fixture every protein has ``|AveExpr - log2FC|`` far larger than
    any tolerance (smallest gap is P4: 21.5), so the inequality is unambiguous.
    """
    res = run_limma(log2_df, SAMPLES_A, SAMPLES_B, "groupA", "groupB")
    res = res.set_index("ProteinName")
    for protein in PROTEINS:
        ave = res.loc[protein, "AveExpr"]
        log2fc = res.loc[protein, "log2FC"]
        # AveExpr is the row mean of the log2 matrix (limma's Amean).
        assert ave == pytest.approx(REF_AVE_EXPR[protein], rel=1e-9, abs=1e-12)
        # ...and therefore must NOT equal the log2FC (the caught regression).
        assert abs(ave - log2fc) > 1.0, (
            f"{protein}: AveExpr ({ave}) collapsed onto log2FC ({log2fc}) -- "
            "the AveExpr==log2FC regression is back"
        )


# ---------------------------------------------------------------------------
# limma: one eBayes-moderated quantity pinned to the Rust golden oracle.
# ---------------------------------------------------------------------------


def test_limma_moderated_t_and_p_match_rust_oracle(log2_df):
    """limma moderated t_stat / pvalue == de.rs::EXPECTED (cross-validated).

    The moderated t is NOT an ordinary two-sample t (eBayes shrinks the
    variance toward the prior), so this pins the empirical-Bayes moderation
    path itself. Provenance: rust/crates/mokume-pipeline/src/de.rs::EXPECTED.
    """
    res = run_limma(log2_df, SAMPLES_A, SAMPLES_B, "groupA", "groupB")
    res = res.set_index("ProteinName")
    for protein, ref_t in REF_MOD_T.items():
        assert res.loc[protein, "t_stat"] == pytest.approx(ref_t, rel=1e-6)
    for protein, ref_p in REF_MOD_P.items():
        assert res.loc[protein, "pvalue"] == pytest.approx(ref_p, rel=1e-6)


# ---------------------------------------------------------------------------
# ROTS: log2FC sign (A - B), the second regression yasset caught.
# ---------------------------------------------------------------------------


def test_rots_log2fc_sign_is_a_minus_b(log2_df):
    """ROTS log2FC == mean(A) - mean(B): positive for a protein high in A.

    ROTS' log2FC (rots.py: ``nanmean(A) - nanmean(B)``) is deterministic and
    RNG-independent, even though ROTS' permutation p-values are not. This pins
    both the numeric value (independently derived) and, crucially, the SIGN:

      * P2 is HIGH in A (mean_A > mean_B) -> log2FC must be POSITIVE, and
      * P1 is LOW in A  (mean_A < mean_B) -> log2FC must be NEGATIVE.

    The pre-fix bug computed ``mean(B) - mean(A)``, flipping every sign; that
    would make P2 negative and P1 positive and fail this test.
    """
    # Guard the fixture's own premise: P2 is genuinely high in A, P1 low in A.
    assert REF_MEAN_A["P2"] > REF_MEAN_B["P2"]
    assert REF_MEAN_A["P1"] < REF_MEAN_B["P1"]

    res = run_rots(
        log2_df, SAMPLES_A, SAMPLES_B, "groupA", "groupB", n_boot=30, seed=42
    )
    res = res.set_index("ProteinName")

    for protein in PROTEINS:
        got = res.loc[protein, "log2FC"]
        # Numeric value: independently-derived mean(A) - mean(B).
        assert got == pytest.approx(REF_LOG2FC[protein], rel=1e-9, abs=1e-12)
        # Sign must agree with the data direction (the flip the reviewer caught).
        assert np.sign(got) == np.sign(REF_MEAN_A[protein] - REF_MEAN_B[protein])

    # Explicit, human-readable sign checks on the two clearest proteins.
    assert res.loc["P2", "log2FC"] > 0, "P2 is high in A; log2FC must be > 0"
    assert res.loc["P1", "log2FC"] < 0, "P1 is low in A; log2FC must be < 0"


def test_rots_log2fc_matches_limma_log2fc_sign(log2_df):
    """ROTS and limma share the mean(A) - mean(B) convention (same sign).

    Both methods must report log2FC in the SAME direction; a sign flip in one
    but not the other (the bug) would break this cross-method agreement.
    """
    limma_res = run_limma(log2_df, SAMPLES_A, SAMPLES_B, "groupA", "groupB")
    limma_res = limma_res.set_index("ProteinName")["log2FC"]
    rots_res = run_rots(
        log2_df, SAMPLES_A, SAMPLES_B, "groupA", "groupB", n_boot=30, seed=42
    )
    rots_res = rots_res.set_index("ProteinName")["log2FC"]
    for protein in PROTEINS:
        assert rots_res.loc[protein] == pytest.approx(
            limma_res.loc[protein], rel=1e-9, abs=1e-12
        )
