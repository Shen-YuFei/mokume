"""Tests for the piBAQ (paralog-aware iBAQ) algorithm.

Covers:

* UniProt isoform canonicalisation (``canonicalize_isoform`` /
  ``canonicalize_isoforms_map``)
* Shared-peptide connected-component family discovery
* YAML override loading and merge semantics
* Proportional and equal shared-peptide allocation in ``compute_pibaq``
* End-to-end integration with a toy FASTA + peptide table, including the
  per-protein TPA mode
"""

from __future__ import annotations

import inspect
import textwrap
from pathlib import Path
from typing import Dict, Set

import pandas as pd
import pytest

from mokume.core.constants import (
    CONDITION,
    EVIDENCE_FAMILY_ONLY,
    EVIDENCE_HIGH,
    EVIDENCE_LEVEL,
    EVIDENCE_MEDIUM,
    FAMILY_ID,
    FAMILY_SIZE,
    PIBAQ,
    MOLECULARWEIGHT,
    NORM_INTENSITY,
    PEPTIDE_CANONICAL,
    PROTEIN_NAME,
    SAMPLE_ID,
    TPA,
)
from mokume.io.fasta import (
    canonicalize_isoform,
    canonicalize_isoforms_map,
    digest_fasta_full,
    invert_peptide_index,
)
from mokume.quantification.families import (
    Family,
    discover_families,
    families_by_member,
    load_families_yaml,
    merge_overrides,
)
from mokume.quantification.pibaq import (
    _membership_mask,
    compute_pibaq,
    peptides_to_protein,
)
from mokume.quantification import get_quantification_method


class _NonIterableContainer:
    """Container proving membership does not enumerate its key domain."""

    def __init__(self, values):
        self._values = set(values)

    def __contains__(self, value):
        return value in self._values

    def __iter__(self):
        raise AssertionError("membership filtering must not materialize the container")


def test_compute_pibaq_preserves_public_signature():
    """The static-analysis refactor must not narrow the public call contract."""
    parameters = inspect.signature(compute_pibaq).parameters
    assert tuple(parameters) == (
        "peptide_df",
        "accession_to_peptides",
        "peptide_to_accessions",
        "families",
        "mw_map",
        "min_anchors",
        "high_anchor_threshold",
        "extra_group_cols",
    )
    assert all(
        parameter.kind is inspect.Parameter.POSITIONAL_OR_KEYWORD
        for parameter in tuple(parameters.values())[:4]
    )
    assert all(
        parameter.kind is inspect.Parameter.KEYWORD_ONLY
        for parameter in tuple(parameters.values())[4:]
    )


def test_membership_mask_probes_observed_values_without_iterating_container():
    values = pd.Series(
        [
            "present",
            "absent",
            None,
            float("nan"),
            pd.NA,
            ["present"],
            {"present"},
            {"present": True},
            "present",
        ],
        dtype=object,
    )
    container = _NonIterableContainer({"present", "unobserved"})

    mask = _membership_mask(values, container)
    empty_mask = _membership_mask(pd.Series(dtype=object), container)

    assert mask.tolist() == [
        True,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        True,
    ]
    assert mask.dtype == bool
    assert empty_mask.empty
    assert empty_mask.dtype == bool


# ---------------------------------------------------------------------------
# canonicalize_isoform / canonicalize_isoforms_map
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("P05067", "P05067"),
        ("P05067-1", "P05067"),
        ("P05067-2", "P05067"),
        ("P70255-7", "P70255"),
        ("P02768ups", "P02768ups"),  # non-numeric suffix preserved
        ("Q14213-PRO_0000010047", "Q14213-PRO_0000010047"),  # not pure digits
        ("", ""),
    ],
)
def test_canonicalize_isoform_handles_uniprot_conventions(raw, expected):
    """``-N`` strip rule must apply only when the suffix is purely numeric."""
    assert canonicalize_isoform(raw) == expected


def test_canonicalize_isoforms_map_unions_peptide_sets():
    """All peptides across -2 / -3 isoforms collapse onto the canonical."""
    acc_to_peps = {
        "P70255": {"AAAR", "BBBK"},
        "P70255-2": {"AAAR", "CCCR"},
        "P70255-3": {"DDDR"},
        "Q12345": {"EEER"},  # isolated
    }
    canonical, alias = canonicalize_isoforms_map(acc_to_peps)
    assert canonical["P70255"] == {"AAAR", "BBBK", "CCCR", "DDDR"}
    assert canonical["Q12345"] == {"EEER"}
    assert alias["P70255-2"] == "P70255"
    assert alias["P70255-3"] == "P70255"
    assert alias["P70255"] == "P70255"
    assert alias["Q12345"] == "Q12345"


def test_invert_peptide_index_round_trip():
    """``invert_peptide_index`` is the inverse of ``acc_to_peps``."""
    acc_to_peps = {"A": {"p1", "p2"}, "B": {"p2", "p3"}}
    pep_to_accs = invert_peptide_index(acc_to_peps)
    assert pep_to_accs == {"p1": {"A"}, "p2": {"A", "B"}, "p3": {"B"}}


# ---------------------------------------------------------------------------
# discover_families (connected components)
# ---------------------------------------------------------------------------


def _toy_paralog_indices() -> tuple[Dict[str, Set[str]], Dict[str, Set[str]]]:
    """Two paralog families + one isolated protein, mirroring Hemna patterns.

    Family ACT (3 members, fully overlapping): all share peptides ``s1``..``s4``.
    Family TPM (2 members): share ``t1``, ``t2``; each has one unique anchor.
    Isolated protein ``LONE``: no peptide overlap with any other accession.
    """
    acc_to_peps = {
        "ACTC1": {"s1", "s2", "s3", "s4"},
        "ACTB": {"s1", "s2", "s3", "s4"},
        "ACTG1": {"s1", "s2", "s3", "s4"},
        "TPM1": {"t1", "t2", "u_tpm1"},
        "TPM2": {"t1", "t2", "u_tpm2"},
        "LONE": {"x1", "x2", "x3"},
    }
    return acc_to_peps, invert_peptide_index(acc_to_peps)


def test_discover_families_groups_paralogs_and_keeps_singletons():
    acc_to_peps, pep_to_accs = _toy_paralog_indices()
    families = discover_families(acc_to_peps, pep_to_accs, min_shared=2)

    by_member = families_by_member(families)
    assert by_member["ACTC1"].members == ("ACTB", "ACTC1", "ACTG1")
    assert by_member["TPM1"].members == ("TPM1", "TPM2")
    assert by_member["LONE"].members == ("LONE",)
    assert by_member["LONE"].size == 1


def test_discover_families_min_shared_threshold_isolates_low_overlap_pairs():
    """When ``min_shared`` is raised, weakly-overlapping pairs split apart."""
    acc_to_peps, pep_to_accs = _toy_paralog_indices()
    families = discover_families(acc_to_peps, pep_to_accs, min_shared=5)
    # All edges removed -> 6 singleton families
    assert len(families) == 6
    assert all(f.size == 1 for f in families)


def test_discover_families_rejects_invalid_min_shared():
    with pytest.raises(ValueError):
        discover_families({}, {}, min_shared=0)


# ---------------------------------------------------------------------------
# load_families_yaml + merge_overrides
# ---------------------------------------------------------------------------


def _write_yaml(path: Path, body: str) -> Path:
    path.write_text(textwrap.dedent(body))
    return path


def test_load_families_yaml_parses_well_formed_document(tmp_path: Path):
    yaml_path = _write_yaml(
        tmp_path / "families.yaml",
        """\
        families:
          - name: ACT
            members: [ACTC1, ACTB, ACTG1]
          - name: TPM
            members: [TPM1, TPM2]
        """,
    )
    families = load_families_yaml(yaml_path)
    assert {f.family_id for f in families} == {"ACT", "TPM"}
    assert families[0].members == ("ACTB", "ACTC1", "ACTG1")


@pytest.mark.parametrize(
    "body, err_substr",
    [
        ("families: not_a_list", "must be a list"),
        ("families:\n  - members: [A]", "is missing a 'name' string"),
        ("families:\n  - name: X\n    members: []", "has no members"),
        (
            "families:\n  - name: X\n    members: [A]\n  - name: X\n    members: [B]",
            "duplicate family name",
        ),
    ],
)
def test_load_families_yaml_rejects_malformed_documents(
    tmp_path: Path, body, err_substr
):
    yaml_path = _write_yaml(tmp_path / "bad.yaml", body)
    with pytest.raises(ValueError, match=err_substr):
        load_families_yaml(yaml_path)


def test_merge_overrides_drops_consumed_members_and_keeps_remainder():
    auto = [
        Family(family_id="ACT", members=("ACTB", "ACTC1", "ACTG1")),
        Family(family_id="LONE", members=("LONE",)),
    ]
    override = [Family(family_id="ACT_PINNED", members=("ACTC1", "ACTB"))]
    merged = merge_overrides(auto, override)
    by_id = {f.family_id: f for f in merged}
    assert by_id["ACT_PINNED"].members == ("ACTB", "ACTC1")
    # ACTG1 is the only auto member not pinned; it gets re-emitted as a smaller family.
    leftover = [
        f for f in merged if "ACTG1" in f.members and f.family_id != "ACT_PINNED"
    ]
    assert len(leftover) == 1
    assert leftover[0].members == ("ACTG1",)
    # Untouched singleton survives verbatim.
    assert by_id["LONE"].members == ("LONE",)


def test_merge_overrides_no_op_when_overrides_absent():
    auto = [Family(family_id="A", members=("a",))]
    assert merge_overrides(auto, None) == auto
    assert merge_overrides(auto, []) == auto


# ---------------------------------------------------------------------------
# compute_pibaq core
# ---------------------------------------------------------------------------


def _build_peptide_df(rows):
    """Build a minimal long-format peptide DataFrame for piBAQ tests."""
    return pd.DataFrame(
        rows,
        columns=[PROTEIN_NAME, PEPTIDE_CANONICAL, SAMPLE_ID, CONDITION, NORM_INTENSITY],
    )


def test_compute_pibaq_singleton_family_matches_baseline_formula():
    """For an isolated protein the per-protein branch reduces to
    ``sum(intensity) / n_proteotypic``, the legacy baseline formula."""
    acc_to_peps = {"LONE": {"p1", "p2", "p3"}}
    pep_to_accs = invert_peptide_index(acc_to_peps)
    families = discover_families(acc_to_peps, pep_to_accs)
    peptide_df = _build_peptide_df(
        [
            ("LONE", "p1", "S1", "C1", 100.0),
            ("LONE", "p2", "S1", "C1", 200.0),
            # p3 not observed -> 2 anchors observed, denominator stays 3 (FASTA)
        ]
    )

    result = compute_pibaq(
        peptide_df,
        acc_to_peps,
        pep_to_accs,
        families,
    )
    assert len(result) == 1
    row = result.iloc[0]
    assert row[PROTEIN_NAME] == "LONE"
    assert row[NORM_INTENSITY] == pytest.approx(300.0)
    assert row[PIBAQ] == pytest.approx(300.0 / 3)
    assert row[FAMILY_ID] == "LONE"
    assert row[FAMILY_SIZE] == 1
    assert row[EVIDENCE_LEVEL] == EVIDENCE_MEDIUM  # 2 anchors


def test_compute_pibaq_proportional_branch_splits_shared_peptide_by_anchor_ratio():
    """Two paralogs in the SAME family, each with their own unique anchor.

    Member A has 3 unique anchors observed, member B has 1; shared-peptide
    intensity is split 3:1 between them. The two paralogs share two distinct
    peptides (``sh1``, ``sh2``) so they cluster into the same family at the
    default ``min_shared=2``; the test then observes only one of those shared
    peptides at quantification time.
    """
    acc_to_peps = {
        "A": {"sh1", "sh2", "ua1", "ua2", "ua3"},
        "B": {"sh1", "sh2", "ub1"},
    }
    pep_to_accs = invert_peptide_index(acc_to_peps)
    families = discover_families(acc_to_peps, pep_to_accs)
    assert families[0].size == 2, "fixture invariant: A and B must co-cluster"
    peptide_df = _build_peptide_df(
        [
            ("A", "ua1", "S1", "C1", 100.0),
            ("A", "ua2", "S1", "C1", 100.0),
            ("A", "ua3", "S1", "C1", 100.0),
            ("B", "ub1", "S1", "C1", 100.0),
            # One of the two shared peptides observed via razor row on A.
            ("A", "sh1", "S1", "C1", 400.0),
        ]
    )

    result = compute_pibaq(peptide_df, acc_to_peps, pep_to_accs, families)
    by_acc = result.set_index(PROTEIN_NAME)
    # A gets 3 unique + 3/4 of shared = 300 + 300 = 600; denominator = 5 total
    # potential peptides (3 proteotypic ua* + 2 shared sh*).
    assert by_acc.loc["A", NORM_INTENSITY] == pytest.approx(300 + 400 * 0.75)
    assert by_acc.loc["A", PIBAQ] == pytest.approx((300 + 300) / 5)
    # B gets 1 unique + 1/4 of shared = 100 + 100 = 200; denominator = 3 total
    # potential peptides (1 proteotypic ub1 + 2 shared sh*).
    assert by_acc.loc["B", NORM_INTENSITY] == pytest.approx(100 + 400 * 0.25)
    assert by_acc.loc["B", PIBAQ] == pytest.approx((100 + 100) / 3)
    # Both rows share the same family metadata.
    assert set(result[FAMILY_ID]) == {by_acc.loc["A", FAMILY_ID]}
    assert all(result[FAMILY_SIZE] == 2)
    assert set(result[EVIDENCE_LEVEL]) == {EVIDENCE_MEDIUM}


def test_compute_pibaq_without_anchors_splits_shared_signal_equally():
    """An unresolved family conserves signal by splitting it among members."""
    acc_to_peps = {
        "ACTC1": {"s1", "s2", "s3"},
        "ACTB": {"s1", "s2", "s3"},
        "ACTG1": {"s1", "s2", "s3"},
    }
    pep_to_accs = invert_peptide_index(acc_to_peps)
    families = discover_families(acc_to_peps, pep_to_accs)
    peptide_df = _build_peptide_df(
        [
            ("ACTC1", "s1", "S1", "C1", 100.0),
            ("ACTC1", "s2", "S1", "C1", 200.0),
            ("ACTB", "s3", "S1", "C1", 300.0),
        ]
    )

    result = compute_pibaq(peptide_df, acc_to_peps, pep_to_accs, families)
    assert len(result) == 3
    # The 600 total intensity is allocated once, not replicated three times.
    assert result[NORM_INTENSITY].sum() == pytest.approx(600.0)
    by_acc = result.set_index(PROTEIN_NAME)
    for member in acc_to_peps:
        assert by_acc.loc[member, NORM_INTENSITY] == pytest.approx(200.0)
        assert by_acc.loc[member, PIBAQ] == pytest.approx(200.0 / 3.0)
    assert set(result[EVIDENCE_LEVEL]) == {EVIDENCE_FAMILY_ONLY}
    assert set(result[FAMILY_SIZE]) == {3}


def test_compute_pibaq_shared_allocation_is_per_sample_intensity_weighted():
    """Shared intensity is split *per sample* by each member's unique-peptide
    intensity in that sample -- not by a single global anchor-count ratio.

    Member ``A`` is observed (with unique peptides) only in ``S1``; in ``S2``
    only ``B`` is present. The shared peptide ``sh1`` is also observed alone
    in ``S3``. Allocation must give ``A`` exactly zero in ``S2`` and split the
    unanchored ``S3`` signal equally.
    """
    acc_to_peps = {
        "A": {"sh1", "sh2", "ua1", "ua2", "ua3"},
        "B": {"sh1", "sh2", "ub1"},
    }
    pep_to_accs = invert_peptide_index(acc_to_peps)
    families = discover_families(acc_to_peps, pep_to_accs)
    assert families[0].size == 2, "fixture invariant: A and B must co-cluster"
    peptide_df = _build_peptide_df(
        [
            # S1: A fully present (3 unique), B present (1 unique), shared via A.
            ("A", "ua1", "S1", "C1", 100.0),
            ("A", "ua2", "S1", "C1", 100.0),
            ("A", "ua3", "S1", "C1", 100.0),
            ("B", "ub1", "S1", "C1", 100.0),
            ("A", "sh1", "S1", "C1", 400.0),
            # S2: A absent (no unique peptides); only B + the shared peptide.
            ("B", "ub1", "S2", "C1", 100.0),
            ("B", "sh1", "S2", "C1", 400.0),
            # S3: neither member has unique signal.
            ("A", "sh1", "S3", "C1", 400.0),
        ]
    )

    result = compute_pibaq(peptide_df, acc_to_peps, pep_to_accs, families)
    by = result.set_index([PROTEIN_NAME, SAMPLE_ID])

    # S1 (both present): A's unique intensity (300) is 3x B's (100), so the
    # 400 shared splits 300/100 -> A=600, B=200.
    assert by.loc[("A", "S1"), NORM_INTENSITY] == pytest.approx(600.0)
    assert by.loc[("B", "S1"), NORM_INTENSITY] == pytest.approx(200.0)
    assert by.loc[("A", "S1"), PIBAQ] == pytest.approx(600.0 / 5)

    # S2 (A absent): a zero-weight member receives no phantom row or signal.
    assert ("A", "S2") not in by.index
    assert by.loc[("B", "S2"), NORM_INTENSITY] == pytest.approx(500.0)
    assert by.loc[("B", "S2"), PIBAQ] == pytest.approx(500.0 / 3)

    # S3 (no anchors): shared signal is split equally among mapped members.
    assert by.loc[("A", "S3"), NORM_INTENSITY] == pytest.approx(200.0)
    assert by.loc[("B", "S3"), NORM_INTENSITY] == pytest.approx(200.0)


def test_compute_pibaq_partial_family_stays_proportional():
    """A zero-anchor member receives no shared signal when peers are anchored."""
    acc_to_peps = {
        "P1": {"u1", "u2", "sh1", "sh2"},  # 2 unique anchors
        "P2": {"u3", "sh1", "sh2"},  # 1 unique anchor
        "P3": {"sh1", "sh2"},  # 0 unique anchors
    }
    pep_to_accs = invert_peptide_index(acc_to_peps)
    families = discover_families(acc_to_peps, pep_to_accs)
    assert families[0].size == 3, "fixture invariant: P1/P2/P3 co-cluster"
    peptide_df = _build_peptide_df(
        [
            ("P1", "u1", "S1", "C1", 100.0),
            ("P1", "u2", "S1", "C1", 100.0),
            ("P2", "u3", "S1", "C1", 100.0),
            ("P1", "sh1", "S1", "C1", 300.0),
        ]
    )

    result = compute_pibaq(peptide_df, acc_to_peps, pep_to_accs, families)
    by = result.set_index(PROTEIN_NAME)

    assert set(result[PROTEIN_NAME]) == {"P1", "P2"}
    # P1: (200 unique + 200 shared) / 4 owned peptides = 100.
    assert by.loc["P1", PIBAQ] == pytest.approx(100.0)
    # P2: (100 unique + 100 shared) / 3 owned peptides = 66.67.
    assert by.loc["P2", PIBAQ] == pytest.approx(200.0 / 3)
    assert "P3" not in by.index
    # At least one member is supported, so the family is member-resolving.
    assert set(result[EVIDENCE_LEVEL]) == {EVIDENCE_MEDIUM}

    family_only = compute_pibaq(
        peptide_df,
        acc_to_peps,
        pep_to_accs,
        families,
        min_anchors=3,
    )
    by = family_only.set_index(PROTEIN_NAME)
    # No member reaches the stricter threshold, so the 300 shared intensity
    # is split 100/100/100 while each member keeps its own unique signal.
    assert by.loc["P1", NORM_INTENSITY] == pytest.approx(300.0)
    assert by.loc["P2", NORM_INTENSITY] == pytest.approx(200.0)
    assert by.loc["P3", NORM_INTENSITY] == pytest.approx(100.0)
    assert family_only[NORM_INTENSITY].sum() == pytest.approx(600.0)
    assert set(family_only[EVIDENCE_LEVEL]) == {EVIDENCE_FAMILY_ONLY}


def test_compute_pibaq_evidence_high_at_threshold():
    """``high`` label requires every family member to clear the high threshold."""
    acc_to_peps = {
        "A": {"u1", "u2", "u3", "u4", "sh"},
        "B": {"v1", "v2", "v3", "v4", "sh"},
    }
    pep_to_accs = invert_peptide_index(acc_to_peps)
    families = discover_families(acc_to_peps, pep_to_accs)
    rows = []
    for pep in ["u1", "u2", "u3", "u4"]:
        rows.append(("A", pep, "S1", "C1", 10.0))
    for pep in ["v1", "v2", "v3", "v4"]:
        rows.append(("B", pep, "S1", "C1", 10.0))
    rows.append(("A", "sh", "S1", "C1", 50.0))
    peptide_df = _build_peptide_df(rows)

    result = compute_pibaq(peptide_df, acc_to_peps, pep_to_accs, families)
    assert set(result[EVIDENCE_LEVEL]) == {EVIDENCE_HIGH}


def test_compute_pibaq_tpa_uses_per_member_mw_in_proportional_branch():
    acc_to_peps = {
        "A": {"ua", "sh"},
        "B": {"ub", "sh"},
    }
    pep_to_accs = invert_peptide_index(acc_to_peps)
    families = discover_families(acc_to_peps, pep_to_accs)
    peptide_df = _build_peptide_df(
        [
            ("A", "ua", "S1", "C1", 100.0),
            ("B", "ub", "S1", "C1", 100.0),
            ("A", "sh", "S1", "C1", 200.0),
        ]
    )
    mw_map = {"A": 50.0, "B": 25.0}
    result = compute_pibaq(
        peptide_df, acc_to_peps, pep_to_accs, families, mw_map=mw_map
    )
    by_acc = result.set_index(PROTEIN_NAME)
    assert by_acc.loc["A", MOLECULARWEIGHT] == pytest.approx(50.0)
    assert by_acc.loc["B", MOLECULARWEIGHT] == pytest.approx(25.0)
    assert by_acc.loc["A", TPA] == pytest.approx(by_acc.loc["A", NORM_INTENSITY] / 50.0)


def test_compute_pibaq_tpa_uses_member_mw_without_anchors():
    acc_to_peps = {
        "ACTC1": {"s1", "s2"},
        "ACTB": {"s1", "s2"},
        "ACTG1": {"s1", "s2"},
    }
    pep_to_accs = invert_peptide_index(acc_to_peps)
    families = discover_families(acc_to_peps, pep_to_accs)
    peptide_df = _build_peptide_df(
        [
            ("ACTC1", "s1", "S1", "C1", 100.0),
            ("ACTB", "s2", "S1", "C1", 200.0),
        ]
    )
    mw_map = {"ACTC1": 10.0, "ACTB": 20.0, "ACTG1": 30.0}
    result = compute_pibaq(
        peptide_df, acc_to_peps, pep_to_accs, families, mw_map=mw_map
    )
    by_acc = result.set_index(PROTEIN_NAME)
    assert list(by_acc.loc[["ACTC1", "ACTB", "ACTG1"], MOLECULARWEIGHT]) == [
        10.0,
        20.0,
        30.0,
    ]
    assert list(by_acc.loc[["ACTC1", "ACTB", "ACTG1"], TPA]) == pytest.approx(
        [10.0, 5.0, 100.0 / 30.0]
    )


def test_compute_pibaq_tpa_finalizes_mixed_branch_weights_together():
    """Every family uses its members' own molecular weights."""
    acc_to_peps = {
        "A": {"ua"},
        "B": {"s1", "s2"},
        "C": {"s1", "s2"},
    }
    pep_to_accs = invert_peptide_index(acc_to_peps)
    families = discover_families(acc_to_peps, pep_to_accs)
    peptide_df = _build_peptide_df(
        [
            ("A", "ua", "S1", "C1", 100.0),
            ("B", "s1", "S1", "C1", 100.0),
            ("C", "s2", "S1", "C1", 200.0),
        ]
    )

    result = compute_pibaq(
        peptide_df,
        acc_to_peps,
        pep_to_accs,
        families,
        mw_map={"A": 10.0, "B": 20.0, "C": 30.0},
    ).set_index(PROTEIN_NAME)

    assert result.loc["A", MOLECULARWEIGHT] == pytest.approx(10.0)
    assert result.loc["A", TPA] == pytest.approx(10.0)
    assert list(result.loc[["B", "C"], MOLECULARWEIGHT]) == [20.0, 30.0]
    assert list(result.loc[["B", "C"], TPA]) == pytest.approx([7.5, 5.0])


def test_compute_pibaq_empty_inputs_return_typed_empty_frame():
    families = [Family(family_id="A", members=("a",))]
    empty = pd.DataFrame(
        columns=[PROTEIN_NAME, PEPTIDE_CANONICAL, SAMPLE_ID, CONDITION, NORM_INTENSITY]
    )
    result = compute_pibaq(empty, {"a": {"p1"}}, {"p1": {"a"}}, families)
    assert list(result.columns)[:4] == [
        PROTEIN_NAME,
        SAMPLE_ID,
        CONDITION,
        NORM_INTENSITY,
    ]
    assert FAMILY_ID in result.columns
    assert result.empty


# ---------------------------------------------------------------------------
# End-to-end: peptides_to_protein CLI driver
# ---------------------------------------------------------------------------


# Two paralogs sharing TWO tryptic peptides (so they cluster into one family
# at the default ``min_shared=2``) plus a unique-to-each tail; complemented
# by an independent control protein with three unique peptides. The unique
# tails contain no internal K/R so the digester emits one peptide per tail
# at min_aa=7.
PIBAQ_TOY_FASTA = textwrap.dedent(
    """\
    >sp|P10001|PARA_HUMAN Toy paralog A
    SHAREDPEPTIDEONEK
    SHAREDPEPTIDETWOK
    AAAAAAAAAAAAAAAK
    >sp|P10002|PARB_HUMAN Toy paralog B
    SHAREDPEPTIDEONEK
    SHAREDPEPTIDETWOK
    CCCCCCCCCCCCCCCK
    >sp|P10003|INDEP_HUMAN Independent
    DDDDDDDDDDDDDDK
    EEEEEEEEEEEEEEK
    FFFFFFFFFFFFFFK
    """
).strip()


@pytest.fixture()
def pibaq_fasta(tmp_path: Path) -> Path:
    path = tmp_path / "pibaq_toy.fasta"
    path.write_text(PIBAQ_TOY_FASTA + "\n")
    return path


def test_peptides_to_protein_pibaq_default_routes_to_compute_pibaq(
    pibaq_fasta: Path, tmp_path: Path
):
    """Run the public CLI driver and confirm the new metadata columns are emitted."""
    peptide_table = tmp_path / "peptides.csv"
    pd.DataFrame(
        [
            # P10001 (PARA): observes its unique tail and one shared peptide.
            ("sp|P10001|PARA_HUMAN", "SHAREDPEPTIDEONEK", "S1", "C1", 1000.0),
            ("sp|P10001|PARA_HUMAN", "AAAAAAAAAAAAAAAK", "S1", "C1", 500.0),
            # P10002 (PARB): observes its unique tail only.
            ("sp|P10002|PARB_HUMAN", "CCCCCCCCCCCCCCCK", "S1", "C1", 500.0),
            # Independent control: three unique peptides observed.
            ("sp|P10003|INDEP_HUMAN", "DDDDDDDDDDDDDDK", "S1", "C1", 1500.0),
            ("sp|P10003|INDEP_HUMAN", "EEEEEEEEEEEEEEK", "S1", "C1", 1500.0),
            ("sp|P10003|INDEP_HUMAN", "FFFFFFFFFFFFFFK", "S1", "C1", 1500.0),
        ],
        columns=[PROTEIN_NAME, PEPTIDE_CANONICAL, SAMPLE_ID, CONDITION, NORM_INTENSITY],
    ).to_csv(peptide_table, index=False)

    output_tsv = tmp_path / "out.tsv"
    peptides_to_protein(
        str(pibaq_fasta),
        str(peptide_table),
        "Trypsin",
        False,
        7,
        40,
        False,
        False,
        0,
        0.0,
        "",
        str(output_tsv),
        False,
        str(tmp_path / "qc.pdf"),
    )
    result = pd.read_csv(output_tsv, sep="\t")
    by_acc = result.set_index(PROTEIN_NAME)

    # Each detected protein appears with its canonical accession (sp|...| stripped).
    assert set(by_acc.index) == {"P10001", "P10002", "P10003"}

    # piBAQ metadata columns must be present and consistently populated.
    for col in (FAMILY_ID, EVIDENCE_LEVEL, FAMILY_SIZE):
        assert col in result.columns

    # The two paralogs share a family of size 2; the independent protein is a singleton.
    assert by_acc.loc["P10001", FAMILY_ID] == by_acc.loc["P10002", FAMILY_ID]
    assert by_acc.loc["P10001", FAMILY_SIZE] == 2
    assert by_acc.loc["P10003", FAMILY_SIZE] == 1


def test_digest_fasta_full_canonicalises_isoforms(tmp_path: Path):
    fasta_path = tmp_path / "isoform.fasta"
    fasta_path.write_text(
        textwrap.dedent(
            """\
            >sp|P10004|CANON_HUMAN Canonical
            DDDDDDDDDDDDDDK
            EEEEEEEEEEEEEEK
            >sp|P10004-2|CANON_HUMAN Isoform 2
            EEEEEEEEEEEEEEK
            GGGGGGGGGGGGGGK
            """
        ).strip()
        + "\n"
    )
    acc_to_peps, pep_to_accs, _ = digest_fasta_full(
        str(fasta_path),
        enzyme="Trypsin",
        min_aa=7,
        max_aa=40,
        canonicalize_isoforms=True,
        compute_mw=False,
    )
    # Both entries collapse onto the canonical accession.
    assert set(acc_to_peps) == {"P10004"}
    assert acc_to_peps["P10004"] >= {
        "DDDDDDDDDDDDDDK",
        "EEEEEEEEEEEEEEK",
        "GGGGGGGGGGGGGGK",
    }
    assert "P10004-2" not in acc_to_peps


def test_total_histone_intensities_matches_accession_protein_names():
    """Regression for issue #19.

    ``ProteinName`` holds UniProt accessions (from ``pg_accessions``), so histones
    must be matched by accession (``histone_proteins``), not only by entry name
    (``histone_entries``); otherwise no histone is found, the intensity falls back
    to 1.0, and ProteomicsRuler CopyNumber is inflated.
    """
    from mokume.core.constants import NORM_INTENSITY, PROTEIN_NAME
    from mokume.model.organism import OrganismDescription
    from mokume.quantification.pibaq import ConcentrationWeightByProteomicRuler

    organism = OrganismDescription(
        name="test",
        genome_size=1,
        histone_proteins=["P0C0S8", "P62805"],  # accessions (match ProteinName)
        histone_entries=["H4_HUMAN"],  # entry name (legacy match, no accession row)
    )
    ruler = ConcentrationWeightByProteomicRuler(
        organism, ploidy=2, concentration_per_cell=1.0
    )
    protein_intensities = pd.DataFrame(
        {
            PROTEIN_NAME: ["P0C0S8", "P62805", "Q99999"],
            NORM_INTENSITY: [10.0, 20.0, 5.0],
        }
    )
    # The two histone accessions sum to 30.0 (not the 1.0 fallback); the
    # non-histone accession is excluded.
    assert ruler.total_histone_intensities(protein_intensities) == 30.0


def test_public_factory_uses_only_pibaq_name():
    canonical = get_quantification_method("pibaq", fasta="unused.fasta")
    assert canonical.name == "piBAQ"

    with pytest.raises(ValueError, match="Unknown quantification method"):
        get_quantification_method("ibaq", fasta="unused.fasta")
