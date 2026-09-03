"""Machine-readable native command catalog and safe argv canonicalization."""

from __future__ import annotations

import importlib
from collections.abc import Iterable
from pathlib import Path
from typing import Any

from mokume.studio.parameter_help import PARAMETER_HELP
from mokume.studio.periphery_commands import (
    PERIPHERY_COMMAND_PATHS,
    PERIPHERY_COMMAND_SPECS,
    periphery_output_paths,
    run_periphery_command,
    validate_periphery_command,
)
from mokume.studio.paths import ProjectPaths


NATIVE_COMMAND_PATHS = {
    ("quantify", "features2proteins"),
    ("quantify", "features2peptides"),
    ("quantify", "peptides2protein"),
    ("correct-batches",),
}

WORKFLOW_CATEGORIES = {
    ("quantify", "features2proteins"): "Quantification",
    ("quantify", "features2peptides"): "Quantification",
    ("quantify", "peptides2protein"): "Quantification",
    ("correct-batches",): "Batch correction",
    ("tissuemap",): "Tissue atlas",
    ("plot", "pca"): "Visualization",
    ("plot", "tsne"): "Visualization",
    ("plot", "de"): "Visualization",
    ("interactive-report",): "Reports",
}

WORKFLOW_DISPLAY_NAMES = {
    ("quantify", "features2proteins"): "Feature → Protein",
    ("quantify", "features2peptides"): "Feature → Peptide",
    ("quantify", "peptides2protein"): "Peptide → Protein",
    ("correct-batches",): "Batch Correction",
    ("tissuemap",): "TissueMap",
    ("plot", "pca"): "PCA",
    ("plot", "tsne"): "t-SNE",
    ("plot", "de"): "DE",
    ("interactive-report",): "Interactive Report",
}

STUDIO_HIDDEN_FLAGS = {
    ("quantify", "features2peptides"): {"generate_filter_config"},
}

# Native workflows keep Clap as their parameter contract; Python periphery
# workflows mirror their argparse contracts. This map controls presentation only.
COMMAND_PRESENTATION: dict[tuple[str, ...], tuple[dict[str, Any], ...]] = {
    ("quantify", "features2proteins"): (
        {
            "id": "input-output",
            "title": "Input & output",
            "flags": "parquet msstats psm output sdrf export_peptides export_ions".split(),
            "common": ("parquet", "msstats", "psm", "output", "sdrf"),
        },
        {
            "id": "quantification",
            "title": "Quantification",
            "flags": (
                "quant_method min_aa min_unique keep_contaminants fasta pibaq_enzyme "
                "pibaq_max_aa pibaq_min_shared pibaq_families_yaml "
                "pibaq_min_anchors directlfq_min_nonan "
                "directlfq_num_samples_quadratic coverage_threshold "
                "min_sample_correlation ratio_fraction_merge"
            ).split(),
            "common": ("quant_method",),
        },
        {
            "id": "normalization-correction",
            "title": "Normalization & correction",
            "flags": (
                "run_normalization sample_normalization normalization_proteins "
                "batch_correction batch_method batch_column batch_covariate "
                "batch_nonparametric batch_mean_only batch_ref irs "
                "irs_reference_sample irs_sdrf_column irs_sdrf_value "
                "irs_reference_regex irs_stat irs_remove_reference"
            ).split(),
            "common": (
                "run_normalization sample_normalization batch_correction irs"
            ).split(),
        },
        {
            "id": "imputation-de",
            "title": "Imputation & differential expression",
            "flags": (
                "impute_method impute_quantile impute_shift impute_scale "
                "impute_n_neighbors de_contrast de_contrast_file de_method "
                "de_ensemble_method de_ensemble_min_k de_log2fc_threshold "
                "de_effect_size_gate de_fdr_threshold de_fdr_method de_output"
            ).split(),
            "common": (
                "impute_method de_contrast de_contrast_file de_log2fc_threshold "
                "de_fdr_threshold de_output"
            ).split(),
        },
        {
            "id": "runtime",
            "title": "Runtime",
            "flags": ("memory", "threads"),
            "common": ("memory", "threads"),
        },
    ),
    ("quantify", "features2peptides"): (
        {
            "id": "input-output",
            "title": "Input & output",
            "flags": (
                "parquet",
                "sdrf",
                "output",
                "save_parquet",
            ),
            "common": ("parquet", "sdrf", "output"),
        },
        {
            "id": "filtering",
            "title": "Filtering",
            "flags": (
                "min_aa",
                "min_unique",
                "keep_shared_peptides",
                "remove_ids",
                "remove_decoy_contaminants",
                "remove_low_frequency_peptides",
                "filter_config",
                "filter_min_intensity",
                "filter_cv_threshold",
                "filter_charge_state",
                "filter_max_missed_cleavages",
                "filter_peptide_fdr",
                "filter_score",
                "filter_exclude_modification",
                "filter_protein_fdr",
                "filter_min_features",
                "filter_max_missing_rate",
            ),
            "common": (
                "min_aa",
                "min_unique",
                "remove_decoy_contaminants",
            ),
        },
        {
            "id": "normalization-aggregation",
            "title": "Normalization & aggregation",
            "flags": (
                "skip_normalization",
                "run_normalization",
                "sample_normalization",
                "log2",
                "aggregation_level",
            ),
            "common": (
                "run_normalization",
                "sample_normalization",
                "log2",
                "aggregation_level",
            ),
        },
        {
            "id": "irs",
            "title": "Internal reference scaling",
            "flags": (
                "irs_channel",
                "irs_autodetect_regex",
                "irs_stat",
                "irs_scope",
            ),
            "common": ("irs_channel", "irs_autodetect_regex"),
        },
    ),
    ("quantify", "peptides2protein"): (
        {
            "id": "input-output",
            "title": "Input & output",
            "flags": ("fasta", "peptides", "output"),
            "common": ("fasta", "peptides", "output"),
        },
        {
            "id": "quantification",
            "title": "Quantification",
            "flags": (
                "quant_method",
                "enzyme",
                "normalize",
                "min_aa",
                "max_aa",
                "directlfq_min_nonan",
                "families_yaml",
                "min_shared",
                "min_anchors",
                "high_anchor_threshold",
            ),
            "common": (
                "quant_method",
                "normalize",
            ),
        },
        {
            "id": "absolute-abundance",
            "title": "Absolute abundance",
            "flags": ("tpa", "ruler", "ploidy", "organism", "cpc"),
            "common": (),
        },
        {
            "id": "qc-runtime",
            "title": "QC & runtime",
            "flags": ("qc_report", "threads"),
            "common": ("qc_report", "threads"),
        },
    ),
    ("correct-batches",): (
        {
            "id": "input-output",
            "title": "Input & output",
            "flags": ("input", "pattern", "comment", "sep", "output"),
            "common": ("input", "pattern", "output"),
        },
        {
            "id": "columns-export",
            "title": "Column mapping & export",
            "flags": (
                "sample_id_column",
                "protein_id_column",
                "pibaq_raw_column",
                "pibaq_corrected_column",
                "export_anndata",
            ),
            "common": (),
        },
    ),
    ("tissuemap",): (
        {
            "id": "input-output",
            "title": "Input & output",
            "flags": ("input", "outdir", "config", "tmt_dataset"),
            "common": ("input", "outdir", "config"),
        },
        {
            "id": "embedding",
            "title": "Embedding",
            "flags": ("impute_method", "embedding_method"),
            "common": ("impute_method", "embedding_method"),
        },
        {
            "id": "runtime",
            "title": "Runtime",
            "flags": ("threads", "dpi"),
            "common": ("threads",),
        },
    ),
    ("plot", "pca"): (
        {
            "id": "input-output",
            "title": "Input & output",
            "flags": ("protein_matrix", "sdrf", "output"),
            "common": ("protein_matrix", "sdrf", "output"),
        },
    ),
    ("plot", "tsne"): (
        {
            "id": "input-output",
            "title": "Input & output",
            "flags": ("input", "pattern", "output"),
            "common": ("input", "pattern", "output"),
        },
    ),
    ("plot", "de"): (
        {
            "id": "input-output",
            "title": "Input & output",
            "flags": ("protein_matrix", "outdir", "sdrf"),
            "common": ("protein_matrix", "outdir", "sdrf"),
        },
        {
            "id": "plots-contrasts",
            "title": "Plots & contrasts",
            "flags": ("volcano", "heatmap", "contrast"),
            "common": ("volcano", "heatmap", "contrast"),
        },
        {
            "id": "thresholds-highlights",
            "title": "Thresholds & highlights",
            "flags": ("log2fc", "fdr", "highlight_protein"),
            "common": ("log2fc", "fdr"),
        },
    ),
    ("interactive-report",): (
        {
            "id": "input-output",
            "title": "Input & output",
            "flags": ("protein_matrix", "sdrf", "output"),
            "common": ("protein_matrix", "sdrf", "output"),
        },
        {
            "id": "contrasts",
            "title": "Contrasts",
            "flags": ("contrast",),
            "common": ("contrast",),
        },
        {
            "id": "thresholds-highlights",
            "title": "Thresholds & highlights",
            "flags": ("log2fc", "fdr", "highlight_protein"),
            "common": ("log2fc", "fdr"),
        },
    ),
}

OUTPUT_ARGUMENTS = {
    "de-output",
    "export-ions",
    "export-peptides",
    "generate-filter-config",
    "log-file",
    "outdir",
    "output",
    "qc-report",
}


class CommandValidationError(ValueError):
    """Raised when a Studio command does not match the installed CLI contract."""


def command_schema() -> list[dict[str, Any]]:
    """Return every compute and analysis workflow supported by Studio."""
    package = importlib.import_module("mokume")
    schema = getattr(package, "command_schema")()
    commands = schema.get("commands", schema) if isinstance(schema, dict) else schema
    return [
        _with_presentation(item)
        for item in (
            *(
                command
                for command in commands
                if tuple(command.get("path", ())) in NATIVE_COMMAND_PATHS
            ),
            *PERIPHERY_COMMAND_SPECS,
        )
    ]


def _with_presentation(command: dict[str, Any]) -> dict[str, Any]:
    """Attach Studio-only grouping without changing native flag metadata."""
    presented = dict(command)
    path = tuple(command["path"])
    presented["category"] = WORKFLOW_CATEGORIES[path]
    presented["display_name"] = WORKFLOW_DISPLAY_NAMES[path]
    descriptions = PARAMETER_HELP[path]
    hidden = STUDIO_HIDDEN_FLAGS.get(path, set())
    presented["flags"] = [
        {
            **flag,
            "help": flag.get("help") or descriptions.get(str(flag.get("id"))),
            **({"studio_hidden": True} if flag.get("id") in hidden else {}),
        }
        for flag in command.get("flags", ())
    ]
    groups = COMMAND_PRESENTATION[path]
    presented["presentation"] = {
        "groups": [
            {
                "id": group["id"],
                "title": group["title"],
                "flags": list(group["flags"]),
                "common": list(group["common"]),
            }
            for group in groups
        ]
    }
    return presented


def validate_and_canonicalize(
    argv: Iterable[str], project_root: str | Path
) -> list[str]:
    """Validate workflow argv and guard every path beneath the active project."""
    original = list(argv)
    spec = _command_spec(original)
    paths = ProjectPaths(project_root)
    canonical = _canonicalize_paths(original, spec, paths)
    path = tuple(spec["path"])
    if path in NATIVE_COMMAND_PATHS:
        package = importlib.import_module("mokume")
        getattr(package, "validate_args")(canonical)
    else:
        validate_periphery_command(canonical)
    _validate_output_paths(canonical, paths)
    return canonical


def execute_command(argv: list[str]) -> None:
    """Run a Studio-supported command through its owning implementation."""
    spec = _command_spec(argv)
    if tuple(spec["path"]) in PERIPHERY_COMMAND_PATHS:
        run_periphery_command(argv)
        return
    package = importlib.import_module("mokume")
    getattr(package, "run")(argv)


def command_paths(argv: Iterable[str]) -> tuple[list[Path], list[Path]]:
    """Return explicit input and output paths from already-canonical argv."""
    tokens = list(argv)
    spec = _command_spec(tokens)
    options = _option_index(spec.get("flags", ()))
    inputs: list[Path] = []
    outputs: list[Path] = []
    index = len(spec["path"])
    while index < len(tokens):
        option, inline_value = _split_option(tokens[index])
        argument = options.get(option)
        if argument is None:
            index += 1
            continue
        count = _value_count(argument)
        values = ([inline_value] if inline_value is not None else []) + tokens[
            index + 1 : index + 1 + count - (inline_value is not None)
        ]
        target = outputs if str(argument.get("long")) in OUTPUT_ARGUMENTS else inputs
        target.extend(_path_values(argument, values))
        index += 1 + count - (inline_value is not None)
    return inputs, periphery_output_paths(tokens, outputs)


def _validate_output_paths(argv: list[str], paths: ProjectPaths) -> None:
    outputs = [paths.resolve_output(path) for path in command_paths(argv)[1]]
    if len(outputs) != len(set(outputs)):
        raise CommandValidationError(
            "workflow resolves multiple outputs to the same path"
        )


def _path_values(argument: dict[str, Any], values: list[str]) -> list[Path]:
    names = argument.get("value_names") or []
    return [
        Path(value)
        for position, value in enumerate(values)
        if names and names[min(position, len(names) - 1)] in {"FILE", "DIR"}
    ]


def _command_spec(argv: list[str]) -> dict[str, Any]:
    for spec in command_schema():
        path = list(spec["path"])
        if argv[: len(path)] == path:
            return spec
    raise CommandValidationError("command is not available in Mokume Studio")


def _canonicalize_paths(
    argv: list[str], spec: dict[str, Any], paths: ProjectPaths
) -> list[str]:
    command_length = len(spec["path"])
    options = _option_index(spec.get("flags", ()))
    canonical = argv[:command_length]
    index = command_length
    while index < len(argv):
        token = argv[index]
        option, inline_value = _split_option(token)
        argument = options.get(option)
        if argument is None:
            canonical.append(token)
            index += 1
            continue
        canonical.append(option)
        value_count = _value_count(argument)
        values = ([inline_value] if inline_value is not None else []) + argv[
            index + 1 : index + 1 + value_count - (inline_value is not None)
        ]
        if len(values) != value_count:
            raise CommandValidationError(f"missing value for {option}")
        canonical.extend(_canonical_values(argument, option, values, paths))
        index += 1 + value_count - (inline_value is not None)
    return canonical


def _option_index(arguments: Iterable[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    options: dict[str, dict[str, Any]] = {}
    for argument in arguments:
        if argument.get("long"):
            options[f"--{argument['long']}"] = argument
        if argument.get("short"):
            options[f"-{argument['short']}"] = argument
    return options


def _split_option(token: str) -> tuple[str, str | None]:
    if token.startswith("--") and "=" in token:
        return tuple(token.split("=", 1))  # type: ignore[return-value]
    return token, None


def _value_count(argument: dict[str, Any]) -> int:
    arity = argument.get("value_arity") or {}
    maximum = arity.get("max")
    minimum = int(arity.get("min", 0))
    if maximum == 0 or not argument.get("value_names"):
        return 0
    if maximum is None or maximum != minimum:
        return max(1, int(minimum))
    return int(maximum)


def _canonical_values(
    argument: dict[str, Any],
    option: str,
    values: list[str],
    paths: ProjectPaths,
) -> list[str]:
    names = argument.get("value_names") or []
    canonical: list[str] = []
    for position, value in enumerate(values):
        value_name = names[min(position, len(names) - 1)] if names else ""
        if value_name not in {"FILE", "DIR"}:
            canonical.append(value)
            continue
        long_name = str(argument.get("long") or option.lstrip("-"))
        resolved = (
            paths.resolve_output(value)
            if long_name in OUTPUT_ARGUMENTS
            else paths.resolve_existing(value)
        )
        canonical.append(str(resolved))
    return canonical
