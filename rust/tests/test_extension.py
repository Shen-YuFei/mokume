"""Binding-level tests: the compiled extension loads and the CLI entry behaves."""

import csv
import functools
import importlib.util
import inspect
import sys
from pathlib import Path

import pandas as pd
import pytest
import pyopenms

import mokume
from mokume.__main__ import main as console_main
from mokume._mokume import run_cli
from mokume._pibaq_digest import build_pibaq_digest
from mokume.model.organism import OrganismDescription
from mokume.quantification.families import discover_families
from mokume.quantification.pibaq import (
    ConcentrationWeightByProteomicRuler,
    compute_pibaq,
    normalize_pibaq,
    peptides_to_protein as compatibility_peptides_to_protein,
)

ProteaseDB = getattr(pyopenms, "ProteaseDB")
FEATURE_PARQUET = Path(__file__).parent / "example" / "feature_wide.parquet"
RUNTIME_PROTEASES = [entry["name"] for entry in mokume.protease_catalog()]


@functools.cache
def _reference_pibaq_module():
    """Load the canonical Python implementation under an isolated module name."""
    source = (
        Path(__file__).parents[2] / "python" / "mokume" / "quantification" / "pibaq.py"
    )
    allocation_source = source.with_name("_pibaq_allocation.py")
    allocation_name = "mokume.quantification._pibaq_allocation"
    name = "_mokume_canonical_pibaq_reference"
    spec = importlib.util.spec_from_file_location(name, source)
    allocation_spec = importlib.util.spec_from_file_location(
        allocation_name, allocation_source
    )
    if (
        spec is None
        or spec.loader is None
        or allocation_spec is None
        or allocation_spec.loader is None
    ):
        raise RuntimeError(f"Could not load canonical piBAQ reference from {source}")
    module = importlib.util.module_from_spec(spec)
    allocation_module = importlib.util.module_from_spec(allocation_spec)
    previous_allocation = sys.modules.get(allocation_name)
    sys.modules[name] = module
    sys.modules[allocation_name] = allocation_module
    try:
        allocation_spec.loader.exec_module(allocation_module)
        spec.loader.exec_module(module)
    finally:
        if previous_allocation is None:
            sys.modules.pop(allocation_name, None)
        else:
            sys.modules[allocation_name] = previous_allocation
    return module


def _pyopenms_text(value):
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def test_version_is_nonempty():
    assert mokume.version()
    assert mokume.__version__ == mokume.version()


def test_native_compute_pibaq_preserves_public_signature():
    parameters = list(inspect.signature(compute_pibaq).parameters.values())
    expected_names = (
        "peptide_df accession_to_peptides peptide_to_accessions families mw_map "
        "min_anchors high_anchor_threshold extra_group_cols"
    ).split()
    assert [parameter.name for parameter in parameters] == expected_names
    assert {parameter.kind for parameter in parameters[:4]} == {
        inspect.Parameter.POSITIONAL_OR_KEYWORD
    }
    assert {parameter.kind for parameter in parameters[4:]} == {
        inspect.Parameter.KEYWORD_ONLY
    }


def test_pibaq_compatibility_wrapper_dispatches_native(monkeypatch, tmp_path):
    """The legacy Python call shape maps losslessly onto the native command."""
    captured = {}
    monkeypatch.setattr(
        mokume, "peptides2protein", lambda **kwargs: captured.update(kwargs)
    )

    compatibility_peptides_to_protein(
        "proteome.fasta",
        "peptides.csv",
        "Trypsin",
        False,
        7,
        30,
        False,
        False,
        0,
        0.0,
        "",
        str(tmp_path / "proteins.tsv"),
        False,
        str(tmp_path / "qc.pdf"),
        families_yaml="families.yaml",
        min_shared=2,
        min_anchors=1,
        high_anchor_threshold=3,
    )

    assert captured == {
        "quant_method": "pibaq",
        "fasta": "proteome.fasta",
        "peptides": "peptides.csv",
        "enzyme": "Trypsin",
        "normalize": False,
        "min_aa": 7,
        "max_aa": 30,
        "tpa": False,
        "ruler": False,
        "ploidy": 0,
        "cpc": 0.0,
        "organism": "",
        "output": str(tmp_path / "proteins.tsv"),
        "min_shared": 2,
        "min_anchors": 1,
        "high_anchor_threshold": 3,
        "families": "families.yaml",
    }


def test_compute_wrapper_rejects_bad_method(tmp_path):
    # A bad --quant-method value is rejected by clap inside the extension and surfaced
    # as a normal RuntimeError -- no interpreter teardown.
    with pytest.raises(RuntimeError):
        mokume.peptides2protein(
            quant_method="definitely-not-a-method",
            peptides="/nonexistent/peptides.csv",
            output=str(tmp_path / "out.tsv"),
        )


def test_run_cli_help_exits_zero():
    assert run_cli(["--help"]) == 0


def test_run_cli_subcommand_help_exits_zero():
    assert run_cli(["quantify", "features2proteins", "--help"]) == 0


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


def _runtime_catalog_case(tmp_path, enzyme):
    """Build a two-sample table from one runtime pyOpenMS digest."""
    fasta = tmp_path / "proteome.fasta"
    fasta.write_text(
        ">sp|P10001|ALPHA\nAKPRAKRDDEFLPMSTY\n"
        ">tr|P10002|BETA\nAKPRAKRDDEYLPQSTY\n"
        ">P10003\nMNGHIVCKTELA\n",
        encoding="utf-8",
    )
    mapping, _ = build_pibaq_digest((str(fasta), enzyme, 1, 1000, 0))
    peptide_to_accessions = {}
    peptide_owners = {}
    for protein in sorted(mapping):
        for peptide in sorted(mapping[protein]):
            peptide_to_accessions.setdefault(peptide, set()).add(protein)
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
    return fasta, mapping, peptide_to_accessions, peptide_table


def _assert_pibaq_equal(actual, expected, extra_sort=()):
    """Compare public piBAQ columns after deterministic row sorting."""
    columns = list(expected.columns)
    sort_by = ["ProteinName", "SampleID", "Condition", *extra_sort]
    actual = actual.sort_values(sort_by).reset_index(drop=True)
    expected = expected.sort_values(sort_by).reset_index(drop=True)
    pd.testing.assert_frame_equal(
        actual[columns],
        expected[columns],
        check_dtype=False,
        check_exact=False,
        rtol=1e-12,
        atol=1e-12,
    )


@pytest.mark.parametrize(
    "enzyme",
    RUNTIME_PROTEASES,
)
def test_rust_pibaq_matches_python_for_runtime_catalog(tmp_path, enzyme):
    fasta, mapping, peptide_to_accessions, peptide_table = _runtime_catalog_case(
        tmp_path, enzyme
    )
    rust_output = tmp_path / "rust.tsv"
    families = discover_families(mapping, peptide_to_accessions, min_shared=2)
    peptide_df = pd.read_csv(peptide_table)
    options = {"min_anchors": 1, "high_anchor_threshold": 3}
    expected = _reference_pibaq_module().compute_pibaq(
        peptide_df, mapping, peptide_to_accessions, families, **options
    )
    compatibility = compute_pibaq(
        peptide_df, mapping, peptide_to_accessions, families, **options
    )
    mokume.peptides2protein(
        quant_method="pibaq",
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
    _assert_pibaq_equal(compatibility, expected)
    _assert_pibaq_equal(pd.read_csv(rust_output, sep="\t"), expected)


def _tpa_extra_group_case():
    """Build a two-fraction shared-peptide case with molecular weights."""
    mapping = {
        "A": {"ua", "shared"},
        "B": {"ub", "shared"},
    }
    peptide_to_accessions = {
        "ua": {"A"},
        "ub": {"B"},
        "shared": {"A", "B"},
    }
    families = discover_families(mapping, peptide_to_accessions, min_shared=1)
    peptide_df = pd.DataFrame(
        [
            ("A", "ua", "S1", "C1", "F1", 100.0),
            ("B", "ub", "S1", "C1", "F1", 300.0),
            ("A", "shared", "S1", "C1", "F1", 400.0),
            ("A", "ua", "S1", "C1", "F2", 200.0),
            ("B", "ub", "S1", "C1", "F2", 200.0),
            ("B", "shared", "S1", "C1", "F2", 600.0),
        ],
        columns=[
            "ProteinName",
            "PeptideCanonical",
            "SampleID",
            "Condition",
            "Fraction",
            "NormIntensity",
        ],
    )
    options = {
        "mw_map": {"A": 50.0, "B": 25.0},
        "extra_group_cols": ["Fraction"],
    }
    return peptide_df, mapping, peptide_to_accessions, families, options


def test_native_compute_pibaq_matches_reference_with_tpa_and_extra_groups():
    peptide_df, mapping, peptide_to_accessions, families, options = (
        _tpa_extra_group_case()
    )
    expected = _reference_pibaq_module().compute_pibaq(
        peptide_df,
        mapping,
        peptide_to_accessions,
        families,
        **options,
    )
    actual = compute_pibaq(
        peptide_df,
        mapping,
        peptide_to_accessions,
        families,
        **options,
    )
    _assert_pibaq_equal(actual, expected, extra_sort=("Fraction",))


def test_rust_python_pibaq_periphery_matches_canonical_reference():
    reference = _reference_pibaq_module()
    normalization_input = pd.DataFrame(
        {
            "SampleID": ["S1", "S1", "S2", "S2"],
            "Condition": ["C1", "C1", "C2", "C2"],
            "PiBAQ": [1.0, 3.0, 2.0, 2.0],
        }
    )
    expected_normalized = reference.normalize_pibaq(normalization_input.copy())
    actual_normalized = normalize_pibaq(normalization_input.copy())
    pd.testing.assert_frame_equal(actual_normalized, expected_normalized)

    organism = OrganismDescription(
        name="pibaq-parity",
        genome_size=1,
        histone_proteins=["H1"],
        histone_entries=["HISTONE_ONE"],
    )
    ruler_input = pd.DataFrame(
        {
            "ProteinName": ["H1", "P1", "HISTONE_ONE", "P2"],
            "Condition": ["C1", "C1", "C2", "C2"],
            "NormIntensity": [20.0, 5.0, 30.0, 10.0],
            "MolecularWeight": [10.0, 20.0, 15.0, 25.0],
        }
    )
    expected_ruler = reference.ConcentrationWeightByProteomicRuler(
        organism, 2, 1.0
    ).apply_by_condition(ruler_input)
    actual_ruler = ConcentrationWeightByProteomicRuler(
        organism, 2, 1.0
    ).apply_by_condition(ruler_input)
    pd.testing.assert_frame_equal(actual_ruler, expected_ruler)


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
            "quantify",
            "peptides2protein",
            "--quant-method",
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
