"""Tests for the shared-aware piBAQ numerator contract.

``peptides_to_protein`` accepts arbitrary peptide tables, including tables
that still carry a ``unique`` flag from QPX. piBAQ must not trust that flag as
an upstream filter: it derives peptide sharing from the FASTA digest, collapses
razor-mirror shared rows once, and re-allocates shared intensity through the
same family-aware path used by ``features2proteins --quant-method pibaq``.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import pandas as pd
import pytest

from mokume.quantification.pibaq import peptides_to_protein

# Two toy homologs share one tryptic peptide; the unique tails carry no
# internal K/R so trypsin produces exactly one extra unique peptide each.
# A third unrelated protein has three independent unique peptides.
TOY_FASTA = textwrap.dedent(
    """\
    >sp|P00001|ACTB_TEST Toy homolog A
    AAAAAAAAAAAAAAAK
    GDFEEMATAASSSLEK
    >sp|P00002|ACTG_TEST Toy homolog B
    CCCCCCCCCCCCCCCK
    GDFEEMATAASSSLEK
    >sp|P00003|INDEP_TEST Independent
    SSSSSSSSSSSSSK
    TTTTTTTTTTTTTK
    YYYYYYYYYYYYYK
    """
).strip()


@pytest.fixture()
def toy_fasta(tmp_path: Path) -> Path:
    path = tmp_path / "toy.fasta"
    path.write_text(TOY_FASTA + "\n")
    return path


def _build_peptide_table(with_unique_col: bool) -> pd.DataFrame:
    rows = [
        # Protein, Peptide, Intensity, unique
        ("sp|P00001|ACTB_TEST", "AAAAAAAAAAAAAAAK", 1000.0, 1),
        ("sp|P00001|ACTB_TEST", "GDFEEMATAASSSLEK", 500.0, 0),
        ("sp|P00002|ACTG_TEST", "CCCCCCCCCCCCCCCK", 2000.0, 1),
        ("sp|P00002|ACTG_TEST", "GDFEEMATAASSSLEK", 500.0, 0),
        ("sp|P00003|INDEP_TEST", "SSSSSSSSSSSSSK", 300.0, 1),
        ("sp|P00003|INDEP_TEST", "TTTTTTTTTTTTTK", 300.0, 1),
        ("sp|P00003|INDEP_TEST", "YYYYYYYYYYYYYK", 300.0, 1),
    ]
    df = pd.DataFrame(
        rows, columns=["ProteinName", "PeptideCanonical", "NormIntensity", "unique"]
    )
    df["SampleID"] = "S1"
    df["BioReplicate"] = 1
    df["Condition"] = "C1"
    if not with_unique_col:
        df = df.drop(columns=["unique"])
    return df


def _run_peptides_to_protein(
    toy_fasta_path: Path, peptides_path: Path, output_path: Path
):
    peptides_to_protein(
        fasta=str(toy_fasta_path),
        peptides=str(peptides_path),
        enzyme="Trypsin",
        normalize=False,
        min_aa=7,
        max_aa=40,
        tpa=True,
        ruler=False,
        ploidy=0,
        cpc=0.0,
        organism="",
        output=str(output_path),
        verbose=False,
        qc_report=str(output_path.with_suffix(".pdf")),
    )
    return pd.read_csv(output_path, sep="\t")


def test_pibaq_path_re_allocates_shared_peptides_via_fasta_mapping(toy_fasta, tmp_path):
    """piBAQ bypasses the upstream ``unique`` filter and re-allocates
    shared peptide intensity from the FASTA digest.

    The toy FASTA has ACTB and ACTG sharing exactly one peptide
    (``GDFEEM...K``) and one unique tail each; at ``min_shared=2`` they
    fall into separate singleton families. Razor-mirror rows for the
    shared peptide collapse to a single observation (``max`` dedup) and
    razor-assign to the family with the highest unique-anchor count;
    when anchors tie the lex-first family_id wins, so ``P00001`` claims
    the 500 of shared signal regardless of the upstream razor flag.
    """
    csv_path = tmp_path / "peptides.csv"
    _build_peptide_table(with_unique_col=True).to_csv(csv_path, index=False)
    res = _run_peptides_to_protein(toy_fasta, csv_path, tmp_path / "out.tsv")

    by_protein = res.set_index("ProteinName")["NormIntensity"].to_dict()
    assert by_protein["P00001"] == pytest.approx(1500.0)  # 1000 unique + 500 shared
    assert by_protein["P00002"] == pytest.approx(2000.0)  # 2000 unique only
    assert by_protein["P00003"] == pytest.approx(900.0)

    # TPA = NormIntensity / MolecularWeight; the per-protein MW is fixed
    # by the FASTA, so the TPA ratio inherits the shared-aware numerator.
    tpa_by_protein = res.set_index("ProteinName")["TPA"].to_dict()
    assert all(v > 0 for v in tpa_by_protein.values())


def test_pibaq_path_handles_input_without_unique_column(toy_fasta, tmp_path):
    """With or without the ``unique`` column the piBAQ path produces the
    same numerator -- it derives proteotypic / shared status from the
    FASTA digest, not the upstream razor flag."""
    csv_path = tmp_path / "peptides_no_unique.csv"
    _build_peptide_table(with_unique_col=False).to_csv(csv_path, index=False)
    res = _run_peptides_to_protein(toy_fasta, csv_path, tmp_path / "out.tsv")

    by_protein = res.set_index("ProteinName")["NormIntensity"].to_dict()
    # Identical to the with-unique-col case: ``max`` deduplication of the
    # razor mirror rows for ``GDFEEM...K`` reduces the shared peptide to a
    # single 500 observation, razor-assigned to the lex-first family.
    assert by_protein["P00001"] == pytest.approx(1500.0)
    assert by_protein["P00002"] == pytest.approx(2000.0)
    assert by_protein["P00003"] == pytest.approx(900.0)


def test_unique_column_accepts_boolean_via_parquet(toy_fasta, tmp_path):
    """Parquet preserves boolean dtype across the round-trip; the QPX
    feature schema declares ``unique`` as boolean. The filter must accept
    ``True`` as proteotypic and ``False`` as shared."""
    df = _build_peptide_table(with_unique_col=True)
    df["unique"] = df["unique"].astype(bool)
    parquet_path = tmp_path / "peptides.parquet"
    df.to_parquet(parquet_path, index=False)

    res = _run_peptides_to_protein(toy_fasta, parquet_path, tmp_path / "out.tsv")
    by_protein = res.set_index("ProteinName")["NormIntensity"].to_dict()
    # Identical to the CSV path: piBAQ re-allocates the shared peptide
    # from the FASTA digest regardless of how the unique flag was encoded.
    assert by_protein["P00001"] == pytest.approx(1500.0)
    assert by_protein["P00002"] == pytest.approx(2000.0)
    assert by_protein["P00003"] == pytest.approx(900.0)
