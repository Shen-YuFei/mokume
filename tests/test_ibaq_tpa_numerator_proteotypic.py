"""Tests for the proteotypic-only numerator contract of `peptides_to_protein`.

The mokume pipeline (``features2proteins``) already strips ``unique != 1``
rows upstream and drops the ``unique`` column before reaching iBAQ. But the
``peptides2protein`` CLI accepts arbitrary peptide tables (including raw
QPX feature parquets), and the previous implementation summed
``NormIntensity`` over every row reaching it — double-counting shared
homologue signal in both the iBAQ and TPA numerators (Hemna's concern
about iBAQpy's TPA path).

These tests pin down a defensive filter inside ``peptides_to_protein``:
when the input carries a ``unique`` column, only proteotypic rows
(``unique`` in ``{1, "1", True}``) flow into the per-protein aggregation.
When the input does not carry the column (the pipeline path), behaviour is
unchanged.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import pandas as pd
import pytest

from mokume.quantification.ibaq import peptides_to_protein

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


def test_unique_column_filters_shared_peptides_from_numerator(toy_fasta, tmp_path):
    """When the peptide table carries `unique`, the iBAQ and TPA numerators
    must aggregate only proteotypic rows. Shared `GDFEEM...K` (500) must
    not contribute to either ACTB_TEST or ACTG_TEST."""
    csv_path = tmp_path / "peptides.csv"
    _build_peptide_table(with_unique_col=True).to_csv(csv_path, index=False)
    res = _run_peptides_to_protein(toy_fasta, csv_path, tmp_path / "out.tsv")

    by_protein = res.set_index("ProteinName")["NormIntensity"].to_dict()
    assert by_protein["sp|P00001|ACTB_TEST"] == pytest.approx(1000.0)
    assert by_protein["sp|P00002|ACTG_TEST"] == pytest.approx(2000.0)
    assert by_protein["sp|P00003|INDEP_TEST"] == pytest.approx(900.0)

    # TPA = NormIntensity / MolecularWeight; the per-protein MW is fixed
    # by the FASTA, so the TPA ratio must inherit the proteotypic-only
    # numerator. We assert relative ordering rather than absolute values
    # because the toy MW is deterministic but not worth hard-coding.
    tpa_by_protein = res.set_index("ProteinName")["TPA"].to_dict()
    assert tpa_by_protein["sp|P00002|ACTG_TEST"] > tpa_by_protein["sp|P00001|ACTB_TEST"]
    assert all(v > 0 for v in tpa_by_protein.values())


def test_missing_unique_column_is_backward_compatible(toy_fasta, tmp_path):
    """When the peptide table omits `unique` (the pipeline path: stages.py
    strips the column after filtering), behaviour is unchanged and every
    surviving row sums into the numerator. The caller is responsible for
    filtering upstream."""
    csv_path = tmp_path / "peptides_no_unique.csv"
    _build_peptide_table(with_unique_col=False).to_csv(csv_path, index=False)
    res = _run_peptides_to_protein(toy_fasta, csv_path, tmp_path / "out.tsv")

    by_protein = res.set_index("ProteinName")["NormIntensity"].to_dict()
    # No filter applied → shared peptide intensity (500) folds into both
    # ACTB_TEST and ACTG_TEST exactly like the historical behaviour.
    assert by_protein["sp|P00001|ACTB_TEST"] == pytest.approx(1500.0)
    assert by_protein["sp|P00002|ACTG_TEST"] == pytest.approx(2500.0)
    assert by_protein["sp|P00003|INDEP_TEST"] == pytest.approx(900.0)


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
    assert by_protein["sp|P00001|ACTB_TEST"] == pytest.approx(1000.0)
    assert by_protein["sp|P00002|ACTG_TEST"] == pytest.approx(2000.0)
    assert by_protein["sp|P00003|INDEP_TEST"] == pytest.approx(900.0)
