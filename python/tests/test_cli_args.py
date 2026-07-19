"""Unit tests for :func:`mokume.core.cli_args.build_features2proteins_argv`.

These exercise the single-source PipelineConfig -> ``features2proteins`` argv
mapping. They need no compiled extension: the function is pure.
"""

from __future__ import annotations

from dataclasses import fields

import pytest

from mokume.core.cli_args import (
    HYBRID_FIELD_OWNERS,
    HybridFieldOwner,
    build_features2proteins_argv,
)
from mokume.pipeline.config import (
    BatchCorrectionConfig,
    DEConfig,
    FilterConfig,
    ImputationConfig,
    IRSConfig,
    InputConfig,
    NormalizationConfig,
    OutputConfig,
    PipelineConfig,
    QuantificationConfig,
    RuntimeConfig,
)


def _flag_value(argv: list[str], flag: str) -> str:
    """Return the token following ``flag`` in ``argv``."""
    assert argv.count(flag) == 1
    idx = argv.index(flag)
    return argv[idx + 1]


def _minimal_config(**overrides) -> PipelineConfig:
    """Build a minimal PipelineConfig; only InputConfig.parquet is required."""
    overrides.setdefault("input", InputConfig(parquet="in.parquet"))
    return PipelineConfig(**overrides)


def test_every_config_leaf_has_an_explicit_hybrid_owner() -> None:
    config = _minimal_config()
    expected = {
        f"{section.name}.{item.name}"
        for section in fields(config)
        for item in fields(getattr(config, section.name))
    }

    assert set(HYBRID_FIELD_OWNERS) == expected
    assert all(
        isinstance(owner, HybridFieldOwner) for owner in HYBRID_FIELD_OWNERS.values()
    )


def test_argv0_is_subcommand() -> None:
    argv = build_features2proteins_argv(_minimal_config(), "out.csv")
    assert argv[0] == "features2proteins"


def test_required_flags_present_with_values() -> None:
    config = _minimal_config(
        quantification=QuantificationConfig(method="sum"),
    )
    argv = build_features2proteins_argv(config, "out.csv")

    assert _flag_value(argv, "--parquet") == "in.parquet"
    assert _flag_value(argv, "--output") == "out.csv"
    assert _flag_value(argv, "--output-format") == "python-compatible"
    assert _flag_value(argv, "--quant-method") == "sum"


def test_quant_method_is_lowercased() -> None:
    config = _minimal_config(quantification=QuantificationConfig(method="MaxLFQ"))
    argv = build_features2proteins_argv(config, "out.csv")
    assert _flag_value(argv, "--quant-method") == "maxlfq"


def test_sample_normalization_default_lowercased() -> None:
    # NormalizationConfig default sample_method is "globalMedian".
    argv = build_features2proteins_argv(_minimal_config(), "out.csv")
    assert _flag_value(argv, "--sample-normalization") == "globalmedian"
    assert _flag_value(argv, "--run-normalization") == "median"


def test_none_input_fields_are_omitted() -> None:
    # No sdrf / fasta provided -> those flags must not appear.
    argv = build_features2proteins_argv(_minimal_config(), "out.csv")
    assert "--sdrf" not in argv
    assert "--fasta" not in argv


def test_optional_input_fields_emitted_when_set() -> None:
    config = _minimal_config(
        input=InputConfig(parquet="in.parquet", sdrf="s.tsv", fasta_file="db.fasta")
    )
    argv = build_features2proteins_argv(config, "out.csv")
    assert _flag_value(argv, "--sdrf") == "s.tsv"
    assert _flag_value(argv, "--fasta") == "db.fasta"


def test_remove_contaminants_toggle() -> None:
    remove = build_features2proteins_argv(
        _minimal_config(filtering=FilterConfig(remove_contaminants=True)), "out.csv"
    )
    assert "--remove-contaminants" in remove
    assert "--keep-contaminants" not in remove

    keep = build_features2proteins_argv(
        _minimal_config(filtering=FilterConfig(remove_contaminants=False)), "out.csv"
    )
    assert "--keep-contaminants" in keep
    assert "--remove-contaminants" not in keep


def test_filter_thresholds_are_stringified() -> None:
    config = _minimal_config(filtering=FilterConfig(min_aa=9, min_unique_peptides=3))
    argv = build_features2proteins_argv(config, "out.csv")
    assert _flag_value(argv, "--min-aa") == "9"
    assert _flag_value(argv, "--min-unique") == "3"


def test_ibaq_families_omitted_when_none() -> None:
    argv = build_features2proteins_argv(_minimal_config(), "out.csv")
    assert "--ibaq-families" not in argv


def test_ibaq_families_emitted_when_set() -> None:
    config = _minimal_config(
        quantification=QuantificationConfig(ibaq_families_yaml="fam.yaml")
    )
    argv = build_features2proteins_argv(config, "out.csv")
    assert _flag_value(argv, "--ibaq-families") == "fam.yaml"


def test_directlfq_cores_omitted_when_none() -> None:
    argv = build_features2proteins_argv(_minimal_config(), "out.csv")
    assert "--directlfq-cores" not in argv


def test_directlfq_cores_emitted_when_set() -> None:
    config = _minimal_config(quantification=QuantificationConfig(directlfq_num_cores=4))
    argv = build_features2proteins_argv(config, "out.csv")
    assert _flag_value(argv, "--directlfq-cores") == "4"


def test_batch_correction_flags_never_emitted() -> None:
    # Batch correction is Python-owned post-processing, so the kernel argv must
    # never carry --batch-correction / --batch-method -- even when batch is
    # enabled. Emitting them would double-correct the matrix.
    off = build_features2proteins_argv(_minimal_config(), "out.csv")
    assert "--batch-correction" not in off
    assert "--batch-method" not in off

    config = _minimal_config(batch=BatchCorrectionConfig(enabled=True, method="run"))
    on = build_features2proteins_argv(config, "out.csv")
    assert "--batch-correction" not in on
    assert "--batch-method" not in on


def test_other_python_owned_postprocessing_is_not_emitted() -> None:
    config = _minimal_config(
        quantification=QuantificationConfig(coverage_threshold=0.75),
        irs=IRSConfig(enabled=True, reference_samples=["pool-a"]),
        imputation=ImputationConfig(enabled=True, method="minprob"),
        de=DEConfig(enabled=True, method="limma"),
        output=OutputConfig(
            plot_volcano=True,
            interactive_report=True,
            export_anndata=True,
        ),
    )
    argv = build_features2proteins_argv(config, "out.csv")

    for flag in [
        "--coverage-threshold",
        "--irs",
        "--irs-reference-samples",
        "--impute",
        "--de",
        "--plot-volcano",
        "--interactive-report",
        "--export-anndata",
    ]:
        assert flag not in argv


def test_top3_forwarded_as_top3() -> None:
    config = _minimal_config(quantification=QuantificationConfig(method="top3"))
    argv = build_features2proteins_argv(config, "out.csv")
    assert _flag_value(argv, "--quant-method") == "top3"
    assert "--topn" not in argv


def test_topn_other_than_three_translated_to_topn_flag() -> None:
    config = _minimal_config(quantification=QuantificationConfig(method="top5"))
    argv = build_features2proteins_argv(config, "out.csv")
    assert _flag_value(argv, "--quant-method") == "topn"
    assert _flag_value(argv, "--topn") == "5"


def test_normalization_run_method_lowercased() -> None:
    config = _minimal_config(
        normalization=NormalizationConfig(run_method="MEAN", sample_method="TMM")
    )
    argv = build_features2proteins_argv(config, "out.csv")
    assert _flag_value(argv, "--run-normalization") == "mean"
    assert _flag_value(argv, "--sample-normalization") == "tmm"


def test_normalization_proteins_file_is_forwarded() -> None:
    config = _minimal_config(
        normalization=NormalizationConfig(proteins_file="stable-proteins.txt")
    )
    argv = build_features2proteins_argv(config, "out.csv")
    assert _flag_value(argv, "--normalization-proteins") == "stable-proteins.txt"


def test_nondefault_ratio_fraction_merge_is_forwarded() -> None:
    config = _minimal_config(
        quantification=QuantificationConfig(method="ratio", ratio_fraction_merge="max")
    )
    argv = build_features2proteins_argv(config, "out.csv")
    assert _flag_value(argv, "--ratio-fraction-merge") == "max"


def test_default_ratio_fraction_merge_keeps_argv_stable() -> None:
    config = _minimal_config(quantification=QuantificationConfig(method="ratio"))
    argv = build_features2proteins_argv(config, "out.csv")
    assert "--ratio-fraction-merge" not in argv


def test_ratio_reference_context_is_forwarded_without_enabling_rust_irs() -> None:
    config = _minimal_config(
        quantification=QuantificationConfig(method="ratio"),
        irs=IRSConfig(
            reference_samples=["pool-a", "pool-b"],
            reference_regex="bridge|reference",
        ),
    )
    argv = build_features2proteins_argv(config, "out.csv")

    assert [
        argv[index + 1]
        for index, token in enumerate(argv)
        if token == "--irs-reference-sample"
    ] == ["pool-a", "pool-b"]
    assert _flag_value(argv, "--irs-reference-regex") == "bridge|reference"
    assert "--irs" not in argv


def test_ratio_reference_sample_with_comma_is_forwarded_losslessly() -> None:
    config = _minimal_config(
        quantification=QuantificationConfig(method="ratio"),
        irs=IRSConfig(reference_samples=["pool, batch A"]),
    )
    argv = build_features2proteins_argv(config, "out.csv")

    assert _flag_value(argv, "--irs-reference-sample") == "pool, batch A"
    assert "--irs-reference-samples" not in argv


def test_empty_ratio_reference_list_preserves_automatic_detection() -> None:
    config = _minimal_config(
        quantification=QuantificationConfig(method="ratio"),
        irs=IRSConfig(reference_samples=[]),
    )
    argv = build_features2proteins_argv(config, "out.csv")
    assert "--irs-reference-sample" not in argv


def test_quantification_exports_are_forwarded_to_compatible_methods() -> None:
    peptide_config = _minimal_config(
        quantification=QuantificationConfig(method="sum"),
        output=OutputConfig(export_peptides="peptides.csv"),
    )
    peptide_argv = build_features2proteins_argv(peptide_config, "out.csv")
    assert _flag_value(peptide_argv, "--export-peptides") == "peptides.csv"

    ion_config = _minimal_config(
        quantification=QuantificationConfig(method="directlfq"),
        output=OutputConfig(export_ions="ions.csv"),
    )
    ion_argv = build_features2proteins_argv(ion_config, "out.csv")
    assert _flag_value(ion_argv, "--export-ions") == "ions.csv"


@pytest.mark.parametrize(
    ("config", "message"),
    [
        (
            _minimal_config(runtime=RuntimeConfig(backend="rust", duckdb_memory="8GB")),
            "duckdb_memory",
        ),
        (
            _minimal_config(runtime=RuntimeConfig(backend="rust", duckdb_threads=2)),
            "duckdb_threads",
        ),
        (
            _minimal_config(
                quantification=QuantificationConfig(ion_alignment="hierarchical")
            ),
            "ion_alignment",
        ),
    ],
)
def test_unsupported_hybrid_settings_are_rejected(config, message) -> None:
    with pytest.raises(ValueError, match=message):
        build_features2proteins_argv(config, "out.csv")


def test_disabled_ion_alignment_variants_are_accepted_and_omitted() -> None:
    for ion_alignment in [None, "none", " NONE "]:
        config = _minimal_config(
            quantification=QuantificationConfig(ion_alignment=ion_alignment)
        )
        argv = build_features2proteins_argv(config, "out.csv")
        assert "--ion-alignment" not in argv


@pytest.mark.parametrize(
    ("method", "output", "flag"),
    [
        (
            "directlfq",
            OutputConfig(export_peptides="peptides.csv"),
            "--export-peptides",
        ),
        ("sum", OutputConfig(export_ions="ions.csv"), "--export-ions"),
    ],
)
def test_incompatible_exports_are_forwarded_for_rust_validation(
    method, output, flag
) -> None:
    config = _minimal_config(
        quantification=QuantificationConfig(method=method), output=output
    )
    argv = build_features2proteins_argv(config, "out.csv")
    assert _flag_value(argv, flag) in {"peptides.csv", "ions.csv"}


def test_all_tokens_are_strings() -> None:
    config = _minimal_config(
        input=InputConfig(parquet="in.parquet", sdrf="s.tsv", fasta_file="db.fasta"),
        quantification=QuantificationConfig(directlfq_num_cores=2),
        batch=BatchCorrectionConfig(enabled=True, method="column"),
    )
    argv = build_features2proteins_argv(config, "out.csv")
    assert all(isinstance(token, str) for token in argv)
