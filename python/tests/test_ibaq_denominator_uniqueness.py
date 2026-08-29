"""Tests for the cross-protein uniqueness contract of ``extract_fasta``.

These tests pin down the behaviour fixed in ``fix/ibaq-denominator``: the
iBAQ denominator must count only peptides that are unique across the full
FASTA database, mirroring the proteotypic numerator. Counting all theoretical
peptides per protein silently deflates large homologous families (myosin,
tubulin, actin, histone) by 3-20x.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from mokume.io.fasta import extract_fasta as extract_fasta_io
from mokume.quantification.pibaq import extract_fasta

# Two toy proteins that share three tryptic peptides and each carry one
# extra cleanly-tryptic peptide that is unique to that accession. The unique
# tails contain no internal K or R so trypsin produces exactly one peptide
# per tail at min_aa=7. With the historical (intra-protein-dedup) denominator
# both proteins would report 4 theoretical peptides; with the cross-protein
# unique denominator they each report 1.
TOY_HOMOLOGS = textwrap.dedent(
    """\
    >sp|P00001|ACTB_TEST Toy beta-actin homolog
    MAGINRSILEGDIVHFGNERFFCPETLFQPSFLGMESCGIHETTFNSIMK
    GDFEEMATAASSSLEK
    AAAAAAAAAAAAAAAK
    >sp|P00002|ACTG_TEST Toy gamma-actin homolog
    MAGINRSILEGDIVHFGNERFFCPETLFQPSFLGMESCGIHETTFNSIMK
    GDFEEMATAASSSLEK
    CCCCCCCCCCCCCCCK
    """
).strip()

# A completely independent protein with no homologs in the FASTA. Each
# tryptic peptide has no internal K/R so they survive the digest as single
# peptides and all count as unique.
TOY_INDEPENDENT = textwrap.dedent(
    """\
    >sp|P00003|INDEP_TEST Independent toy protein
    SSSSSSSSSSSSSK
    TTTTTTTTTTTTTK
    YYYYYYYYYYYYYK
    """
).strip()


@pytest.fixture()
def homolog_fasta(tmp_path: Path) -> Path:
    path = tmp_path / "toy_homologs.fasta"
    path.write_text(TOY_HOMOLOGS + "\n" + TOY_INDEPENDENT + "\n")
    return path


def test_homolog_count_reflects_only_proteotypic_peptides(homolog_fasta):
    """ACTB_TEST and ACTG_TEST share two peptides; only the trailing unique
    peptide on each protein should contribute to the denominator."""
    uniquepepcounts, mw_dict, found = extract_fasta(
        fasta=str(homolog_fasta),
        enzyme="Trypsin",
        proteins=["sp|P00001|ACTB_TEST", "sp|P00002|ACTG_TEST"],
        min_aa=7,
        max_aa=40,
        tpa=False,
    )

    assert found == {"sp|P00001|ACTB_TEST", "sp|P00002|ACTG_TEST"}
    assert uniquepepcounts["sp|P00001|ACTB_TEST"] == 1
    assert uniquepepcounts["sp|P00002|ACTG_TEST"] == 1
    assert mw_dict == {}


def test_independent_protein_counts_all_peptides(homolog_fasta):
    """A protein with no homolog should report every digested peptide as unique."""
    uniquepepcounts, _, _ = extract_fasta(
        fasta=str(homolog_fasta),
        enzyme="Trypsin",
        proteins=["sp|P00003|INDEP_TEST"],
        min_aa=7,
        max_aa=40,
        tpa=False,
    )
    assert uniquepepcounts["sp|P00003|INDEP_TEST"] == 3


def test_tpa_returns_molecular_weight(homolog_fasta):
    uniquepepcounts, mw_dict, _ = extract_fasta(
        fasta=str(homolog_fasta),
        enzyme="Trypsin",
        proteins=["sp|P00001|ACTB_TEST"],
        min_aa=7,
        max_aa=40,
        tpa=True,
    )
    assert "sp|P00001|ACTB_TEST" in mw_dict
    assert mw_dict["sp|P00001|ACTB_TEST"] > 0
    assert uniquepepcounts["sp|P00001|ACTB_TEST"] == 1


def test_unknown_proteins_raise(tmp_path):
    fasta = tmp_path / "tiny.fasta"
    fasta.write_text(TOY_INDEPENDENT + "\n")
    with pytest.raises(ValueError, match="None of the .* proteins were found"):
        extract_fasta(
            fasta=str(fasta),
            enzyme="Trypsin",
            proteins=["sp|NOTFOUND|GONE"],
            min_aa=7,
            max_aa=40,
            tpa=False,
        )


def test_keys_preserve_caller_protein_names(homolog_fasta):
    """The returned ``uniquepepcounts`` dict must use the caller's protein
    names verbatim so downstream maps keyed by the input list still work."""
    proteins_with_prefix = [
        "sp|P00001|ACTB_TEST",
        "sp|P00002|ACTG_TEST",
    ]
    uniquepepcounts, _, found = extract_fasta(
        fasta=str(homolog_fasta),
        enzyme="Trypsin",
        proteins=proteins_with_prefix,
        min_aa=7,
        max_aa=40,
        tpa=False,
    )
    assert set(uniquepepcounts.keys()) == set(proteins_with_prefix)
    assert found == set(proteins_with_prefix)


def test_io_fasta_extract_fasta_matches_quantification(homolog_fasta):
    """The two public ``extract_fasta`` entry points must agree on
    cross-protein unique counts (they previously both shared the same bug)."""
    proteins = [
        "sp|P00001|ACTB_TEST",
        "sp|P00002|ACTG_TEST",
        "sp|P00003|INDEP_TEST",
    ]
    counts_quant, _, _ = extract_fasta(
        fasta=str(homolog_fasta),
        enzyme="Trypsin",
        proteins=proteins,
        min_aa=7,
        max_aa=40,
        tpa=False,
    )
    counts_io, _, _ = extract_fasta_io(
        fasta=str(homolog_fasta),
        enzyme="Trypsin",
        proteins=proteins,
        min_aa=7,
        max_aa=40,
        tpa=False,
    )
    assert counts_quant == counts_io
    assert counts_io["sp|P00001|ACTB_TEST"] == 1
    assert counts_io["sp|P00002|ACTG_TEST"] == 1
    assert counts_io["sp|P00003|INDEP_TEST"] == 3


def test_unparseable_sequence_does_not_crash(tmp_path):
    """A FASTA entry whose sequence cannot be parsed by pyOpenMS (lowercase
    or punctuation, varies across versions) must not prevent clean siblings
    from being quantified."""
    fasta_text = (
        ">sp|P00010|GOOD_TEST Clean protein\n"
        "AAAAAAAAAAAAAAAK\n"
        "CCCCCCCCCCCCCCCK\n"
        ">sp|P00011|JUNK_TEST Unparseable junk content\n"
        "abc!@#xyz\n"
    )
    fasta_path = tmp_path / "with_junk.fasta"
    fasta_path.write_text(fasta_text)
    counts, _, found = extract_fasta(
        fasta=str(fasta_path),
        enzyme="Trypsin",
        proteins=["sp|P00010|GOOD_TEST", "sp|P00011|JUNK_TEST"],
        min_aa=7,
        max_aa=40,
        tpa=False,
    )
    assert "sp|P00010|GOOD_TEST" in found
    assert counts["sp|P00010|GOOD_TEST"] == 2
