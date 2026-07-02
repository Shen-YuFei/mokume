"""Unit tests for :func:`mokume.core.cli_args.build_features2proteins_argv`.

These exercise the single-source PipelineConfig -> ``features2proteins`` argv
mapping. They need no compiled extension: the function is pure.
"""

from __future__ import annotations

from mokume.core.cli_args import build_features2proteins_argv
from mokume.pipeline.config import (
    BatchCorrectionConfig,
    FilterConfig,
    InputConfig,
    NormalizationConfig,
    PipelineConfig,
    QuantificationConfig,
)


def _flag_value(argv: list[str], flag: str) -> str:
    """Return the token following ``flag`` in ``argv``."""
    idx = argv.index(flag)
    return argv[idx + 1]


def _minimal_config(**overrides) -> PipelineConfig:
    """Build a minimal PipelineConfig; only InputConfig.parquet is required."""
    overrides.setdefault("input", InputConfig(parquet="in.parquet"))
    return PipelineConfig(**overrides)


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


def test_batch_correction_flags() -> None:
    off = build_features2proteins_argv(_minimal_config(), "out.csv")
    assert "--batch-correction" not in off
    assert "--batch-method" not in off

    config = _minimal_config(batch=BatchCorrectionConfig(enabled=True, method="run"))
    on = build_features2proteins_argv(config, "out.csv")
    assert "--batch-correction" in on
    assert _flag_value(on, "--batch-method") == "run"


def test_normalization_run_method_lowercased() -> None:
    config = _minimal_config(
        normalization=NormalizationConfig(run_method="MEAN", sample_method="TMM")
    )
    argv = build_features2proteins_argv(config, "out.csv")
    assert _flag_value(argv, "--run-normalization") == "mean"
    assert _flag_value(argv, "--sample-normalization") == "tmm"


def test_all_tokens_are_strings() -> None:
    config = _minimal_config(
        input=InputConfig(parquet="in.parquet", sdrf="s.tsv", fasta_file="db.fasta"),
        quantification=QuantificationConfig(directlfq_num_cores=2),
        batch=BatchCorrectionConfig(enabled=True, method="column"),
    )
    argv = build_features2proteins_argv(config, "out.csv")
    assert all(isinstance(token, str) for token in argv)
