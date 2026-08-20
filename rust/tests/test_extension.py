"""Binding-level tests: the compiled extension loads and the CLI entry behaves."""

import csv
import sys
from pathlib import Path

import pandas as pd
import pytest
import pyopenms

import mokume
from mokume.__main__ import main as console_main
from mokume._mokume import run_cli
from mokume._pibaq_digest import build_pibaq_digest
from mokume.quantification.pibaq import (
    peptides_to_protein as python_peptides_to_protein,
)

ProteaseDB = getattr(pyopenms, "ProteaseDB")
FEATURE_PARQUET = Path(__file__).parent / "example" / "feature_wide.parquet"
RUNTIME_PROTEASES = [entry["name"] for entry in mokume.protease_catalog()]


def _pyopenms_text(value):
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def test_version_is_nonempty():
    assert mokume.version()
    assert mokume.__version__ == mokume.version()


def test_compute_wrapper_rejects_bad_method(tmp_path):
    # A bad --method value is rejected by clap inside the extension and surfaced
    # as a normal RuntimeError -- no interpreter teardown.
    with pytest.raises(RuntimeError):
        mokume.peptides2protein(
            method="definitely-not-a-method",
            peptides="/nonexistent/peptides.csv",
            output=str(tmp_path / "out.tsv"),
        )


def test_run_cli_help_exits_zero():
    assert run_cli(["--help"]) == 0


def test_run_cli_subcommand_help_exits_zero():
    assert run_cli(["features2proteins", "--help"]) == 0


def test_run_cli_unknown_subcommand_is_nonzero():
    assert run_cli(["definitely-not-a-subcommand"]) != 0


def test_runtime_protease_catalog_matches_installed_pyopenms():
    database = ProteaseDB()
    names = []
    database.getAllNames(names)
    expected = []
    for raw_name in names:
        enzyme = database.getEnzyme(raw_name)
        expected.append(
            {
                "name": _pyopenms_text(enzyme.getName()),
                "regex": _pyopenms_text(enzyme.getRegEx()),
                "synonyms": sorted(
                    _pyopenms_text(value) for value in enzyme.getSynonyms()
                ),
            }
        )
    catalog = mokume.protease_catalog()

    assert catalog == sorted(expected, key=lambda entry: entry["name"])


def test_runtime_provider_digests_every_installed_protease(tmp_path):
    fasta = tmp_path / "catalog.fasta"
    fasta.write_text(
        ">P1\nAKPRAKRDDEFLP\n>P2\nAKPRAKRDDEYLP\n",
        encoding="utf-8",
    )
    catalog_hashes = set()
    for entry in mokume.protease_catalog():
        mapping, provenance = build_pibaq_digest(
            (str(fasta), entry["name"], 1, 1000, 0)
        )
        assert mapping["P1"]
        assert mapping["P2"]
        assert provenance["enzyme"] == entry["name"]
        assert provenance["min_aa"] == 1
        assert provenance["max_aa"] == 1000
        assert provenance["missed_cleavages"] == 0
        assert len(provenance["catalog_hash"]) == 64
        catalog_hashes.add(provenance["catalog_hash"])
        for synonym in entry["synonyms"]:
            synonym_mapping, synonym_provenance = build_pibaq_digest(
                (str(fasta), synonym, 1, 1000, 0)
            )
            assert synonym_mapping == mapping
            assert synonym_provenance["enzyme"] == entry["name"]

    assert len(catalog_hashes) == 1


def test_runtime_digest_preserves_fasta_accession_contract(tmp_path):
    fasta = tmp_path / "accessions.fasta"
    fasta.write_text(
        ">sp|P12345-2|PROT_ISOFORM\nAAAA\n"
        ">tr|P12345|PROT_CANONICAL\nCCCC\n"
        ">P02768ups|ALBU_HUMAN_UPS\nDDDD\n"
        ">O13547\nEEEE\n",
        encoding="utf-8",
    )

    mapping, _ = build_pibaq_digest((str(fasta), "no cleavage", 1, 100, 0))

    assert mapping == {
        "P12345": {"AAAA", "CCCC"},
        "P02768ups": {"DDDD"},
        "O13547": {"EEEE"},
    }


@pytest.mark.parametrize(
    "enzyme",
    RUNTIME_PROTEASES,
)
def test_rust_pibaq_matches_python_for_runtime_catalog(tmp_path, enzyme):
    fasta = tmp_path / "proteome.fasta"
    fasta.write_text(
        ">sp|P10001|ALPHA\nAKPRAKRDDEFLPMSTY\n"
        ">tr|P10002|BETA\nAKPRAKRDDEYLPQSTY\n"
        ">P10003\nMNGHIVCKTELA\n",
        encoding="utf-8",
    )
    mapping, _ = build_pibaq_digest((str(fasta), enzyme, 1, 1000, 0))
    peptide_owners = {}
    for protein in sorted(mapping):
        for peptide in sorted(mapping[protein]):
            peptide_owners.setdefault(peptide, protein)

    observations = []
    for index, peptide in enumerate(sorted(peptide_owners), start=1):
        protein = peptide_owners[peptide]
        observations.extend(
            [
                (protein, peptide, "S1", "A", 100.0 + index),
                (protein, peptide, "S2", "B", 250.0 + 2 * index),
            ]
        )
    peptide_table = tmp_path / "peptides.csv"
    pd.DataFrame(
        observations,
        columns=[
            "ProteinName",
            "PeptideCanonical",
            "SampleID",
            "Condition",
            "NormIntensity",
        ],
    ).to_csv(peptide_table, index=False)
    rust_output = tmp_path / "rust.tsv"
    python_output = tmp_path / "python.tsv"

    python_peptides_to_protein(
        fasta=str(fasta),
        peptides=str(peptide_table),
        enzyme=enzyme,
        normalize=False,
        min_aa=1,
        max_aa=1000,
        tpa=False,
        ruler=False,
        ploidy=0,
        cpc=0.0,
        organism="",
        output=str(python_output),
        verbose=False,
        qc_report=str(tmp_path / "unused.pdf"),
        min_shared=2,
        min_anchors=1,
        high_anchor_threshold=3,
    )
    mokume.peptides2protein(
        method="pibaq",
        peptides=str(peptide_table),
        fasta=str(fasta),
        enzyme=enzyme,
        min_aa=1,
        max_aa=1000,
        min_shared=2,
        min_anchors=1,
        high_anchor_threshold=3,
        output=str(rust_output),
    )

    columns = [
        "ProteinName",
        "SampleID",
        "Condition",
        "NormIntensity",
        "PiBAQ",
        "FamilyId",
        "EvidenceLevel",
        "FamilySize",
    ]
    sort_by = ["ProteinName", "SampleID"]
    actual = (
        pd.read_csv(rust_output, sep="\t").sort_values(sort_by).reset_index(drop=True)
    )
    expected = (
        pd.read_csv(python_output, sep="\t").sort_values(sort_by).reset_index(drop=True)
    )
    pd.testing.assert_frame_equal(
        actual[columns],
        expected[columns],
        check_dtype=False,
        check_exact=False,
        rtol=1e-12,
        atol=1e-12,
    )


def test_features2proteins_pibaq_uses_runtime_digest(tmp_path):
    fasta = tmp_path / "proteome.fasta"
    fasta.write_text(
        ">sp|Q86U42|PABP2_HUMAN\nAAAAAAAAAAGAAGGR\n",
        encoding="utf-8",
    )
    output = tmp_path / "proteins.csv"

    mokume.features2proteins(
        parquet=str(FEATURE_PARQUET),
        output=str(output),
        quant_method="pibaq",
        fasta=str(fasta),
        pibaq_enzyme="no cleavage",
        min_aa=1,
        pibaq_max_aa=1000,
        min_unique=1,
        run_normalization="none",
        sample_normalization="none",
    )

    with output.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 1
    assert rows[0]["ProteinName"] == "Q86U42"


def test_console_pibaq_uses_runtime_digest(tmp_path, monkeypatch):
    fasta = tmp_path / "proteome.fasta"
    fasta.write_text(">P1\nAKPRAKRDDEFLP\n", encoding="utf-8")
    peptide_table = tmp_path / "peptides.csv"
    peptide_table.write_text(
        "ProteinName,PeptideCanonical,SampleID,Condition,NormIntensity\n"
        "P1,AKPRAKRDDEFLP,S1,A,100.0\n",
        encoding="utf-8",
    )
    output = tmp_path / "proteins.tsv"

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "mokume",
            "peptides2protein",
            "--method",
            "pibaq",
            "--peptides",
            str(peptide_table),
            "--fasta",
            str(fasta),
            "--enzyme",
            "no cleavage",
            "--min-aa",
            "1",
            "--max-aa",
            "1000",
            "--output",
            str(output),
        ],
    )

    with pytest.raises(SystemExit) as exit_info:
        console_main()

    assert exit_info.value.code == 0
    with output.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["ProteinName"] == "P1"
    assert float(rows[0]["PiBAQ"]) > 0


def test_matrix_level_compute_api():
    matrix = [
        [2.0, 4.0, 8.0, 32.0, 64.0, 128.0],
        [8.0, 16.0, 32.0, 8.0, 16.0, 32.0],
        [32.0, None, 128.0, 4.0, 8.0, 16.0],
    ]
    normalized = mokume.normalize_matrix(
        matrix,
        "median",
        ["a1", "a2", "a3", "b1", "b2", "b3"],
        2,
    )
    assert len(normalized) == len(matrix)
    assert all(len(row) == 6 for row in normalized)

    imputed = mokume.impute_matrix(normalized, "mean", threads=2)
    assert imputed[2][1] is not None

    results = mokume.differential_expression(
        ["P1", "P2", "P3"],
        imputed,
        3,
        3,
        "limma",
        condition_a="case",
        condition_b="control",
        threads=2,
    )
    assert results
    assert {"ProteinName", "log2FC", "pvalue", "adj_pvalue", "significance"} <= set(
        results[0]
    )
