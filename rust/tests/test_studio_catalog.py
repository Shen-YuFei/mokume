"""Presentation guarantees for the Studio command catalog."""

from __future__ import annotations

from collections import Counter

import pytest

from mokume.studio.catalog import (
    WORKFLOW_DISPLAY_NAMES,
    command_paths,
    command_schema,
    execute_command,
    validate_and_canonicalize,
)


CRITICAL_COMMON_FLAGS = {
    ("quantify", "features2proteins"): {
        "parquet",
        "msstats",
        "psm",
        "output",
    },
    ("quantify", "features2peptides"): {"parquet", "output"},
    ("quantify", "peptides2protein"): {"peptides", "output"},
    ("correct-batches",): {"input", "output"},
    ("tissuemap",): {"input", "outdir"},
    ("plot", "pca"): {"protein_matrix", "sdrf", "output"},
    ("plot", "tsne"): {"input", "output"},
    ("plot", "de"): {"protein_matrix", "outdir", "contrast"},
    ("interactive-report",): {"protein_matrix", "sdrf", "output", "contrast"},
}

EXPECTED_WORKFLOWS = {
    ("quantify", "features2proteins"),
    ("quantify", "features2peptides"),
    ("quantify", "peptides2protein"),
    ("correct-batches",),
    ("tissuemap",),
    ("plot", "pca"),
    ("plot", "tsne"),
    ("plot", "de"),
    ("interactive-report",),
}

MINIMUM_PERIPHERY_COMMANDS = (
    ["tissuemap", "--input", "data", "--outdir", "atlas", "--threads", "24"],
    [
        "plot",
        "pca",
        "--protein-matrix",
        "proteins.csv",
        "--sdrf",
        "sdrf.tsv",
        "--output",
        "pca.png",
    ],
    ["plot", "tsne", "--input", "data", "--output", "tsne.pdf"],
    [
        "plot",
        "de",
        "--protein-matrix",
        "proteins.csv",
        "--outdir",
        "plots",
        "--volcano",
        "--contrast",
        "comparison",
        "a",
        "b",
        "de.csv",
    ],
    [
        "interactive-report",
        "--protein-matrix",
        "proteins.csv",
        "--sdrf",
        "sdrf.tsv",
        "--output",
        "report.html",
        "--contrast",
        "comparison",
        "a",
        "b",
        "de.csv",
    ],
)


def test_studio_exposes_all_analysis_workflows_by_category():
    """Studio includes every user-facing compute/analysis leaf command."""
    commands = command_schema()

    assert {tuple(command["path"]) for command in commands} == EXPECTED_WORKFLOWS
    assert all(command["category"] for command in commands)
    displayed = {
        tuple(command["path"]): command["display_name"] for command in commands
    }
    assert displayed == WORKFLOW_DISPLAY_NAMES
    assert displayed[("plot", "tsne")] == "t-SNE"


@pytest.mark.parametrize(
    "command", command_schema(), ids=lambda item: "-".join(item["path"])
)
def test_studio_groups_every_command_flag_once(command):
    """Grouping is complete while required flags always remain visible."""
    flags = {
        flag["id"]: flag
        for flag in command["flags"]
        if not flag.get("global") and not flag.get("studio_hidden")
    }
    groups = command["presentation"]["groups"]
    grouped = [flag_id for group in groups for flag_id in group["flags"]]
    common = {flag_id for group in groups for flag_id in group["common"]}

    assert all(set(group) == {"id", "title", "flags", "common"} for group in groups)
    assert Counter(grouped) == Counter(flags.keys())
    assert common <= flags.keys()
    assert {
        flag_id for flag_id, flag in flags.items() if flag.get("required")
    } <= common


def test_key_inputs_and_outputs_are_not_hidden_as_advanced():
    """Important files stay visible even where the native schema makes them optional."""
    commands = {tuple(command["path"]): command for command in command_schema()}

    for path, expected in CRITICAL_COMMON_FLAGS.items():
        groups = commands[path]["presentation"]["groups"]
        common = {flag_id for group in groups for flag_id in group["common"]}
        assert expected <= common


def test_every_studio_parameter_has_a_description():
    """Workflow forms never fall back to bare value placeholders."""
    for command in command_schema():
        missing = [
            flag["id"]
            for flag in command["flags"]
            if not flag.get("global")
            and not flag.get("studio_hidden")
            and not flag.get("help")
        ]
        assert missing == [], f"{' '.join(command['path'])}: {missing}"


def test_filter_config_generator_is_hidden_from_studio():
    """Template generation is a Studio action, not an analysis parameter."""
    command = next(
        command
        for command in command_schema()
        if command["path"] == ["quantify", "features2peptides"]
    )
    generator = next(
        flag for flag in command["flags"] if flag["id"] == "generate_filter_config"
    )
    grouped = {
        flag_id
        for group in command["presentation"]["groups"]
        for flag_id in group["flags"]
    }

    assert generator["studio_hidden"] is True
    assert generator["id"] not in grouped


def test_memory_is_a_common_runtime_parameter():
    """Memory stays visible beside threads instead of moving under Advanced."""
    commands = {tuple(command["path"]): command for command in command_schema()}
    groups = commands[("quantify", "features2proteins")]["presentation"]["groups"]
    runtime = next(group for group in groups if group["id"] == "runtime")

    assert "memory" in runtime["common"]


def test_features2proteins_keeps_only_core_workflow_choices_common():
    """Method-specific tuning stays under Advanced until its method is selected."""
    commands = {tuple(command["path"]): command for command in command_schema()}
    groups = {
        group["id"]: group
        for group in commands[("quantify", "features2proteins")]["presentation"][
            "groups"
        ]
    }

    assert groups["quantification"]["common"] == ["quant_method"]
    assert groups["normalization-correction"]["common"] == [
        "run_normalization",
        "sample_normalization",
        "batch_correction",
        "irs",
    ]
    assert groups["imputation-de"]["common"] == [
        "impute_method",
        "de_contrast",
        "de_contrast_file",
        "de_log2fc_threshold",
        "de_fdr_threshold",
        "de_output",
    ]
    assert groups["runtime"]["common"] == ["memory", "threads"]


def test_secondary_workflows_hide_specialized_defaults():
    """Secondary workflows expose core choices and keep tuning in Advanced."""
    commands = {tuple(command["path"]): command for command in command_schema()}
    peptide_groups = {
        group["id"]: group
        for group in commands[("quantify", "features2peptides")]["presentation"][
            "groups"
        ]
    }
    protein_groups = {
        group["id"]: group
        for group in commands[("quantify", "peptides2protein")]["presentation"][
            "groups"
        ]
    }

    assert (
        "skip_normalization"
        not in peptide_groups["normalization-aggregation"]["common"]
    )
    assert protein_groups["quantification"]["common"] == [
        "quant_method",
        "normalize",
    ]
    assert protein_groups["absolute-abundance"]["common"] == []
    assert protein_groups["qc-runtime"]["common"] == ["qc_report", "threads"]


def test_periphery_paths_are_canonical_and_report_outputs_are_expanded(tmp_path):
    """Python workflows retain workspace guards and record every report artifact."""
    project = tmp_path / "project"
    project.mkdir()
    for name in ("proteins.csv", "sdrf.tsv", "a.csv", "b.csv"):
        (project / name).write_text("fixture\n", encoding="utf-8")

    canonical = validate_and_canonicalize(
        [
            "interactive-report",
            "--protein-matrix",
            "proteins.csv",
            "--sdrf",
            "sdrf.tsv",
            "--output",
            "report.html",
            "--contrast",
            "a",
            "control",
            "case",
            "a.csv",
            "--contrast",
            "b",
            "control",
            "case",
            "b.csv",
        ],
        project,
    )
    inputs, outputs = command_paths(canonical)

    assert set(inputs) == {
        project / "proteins.csv",
        project / "sdrf.tsv",
        project / "a.csv",
        project / "b.csv",
    }
    assert outputs == [project / "report_a.html", project / "report_b.html"]


def test_each_periphery_workflow_accepts_its_minimum_parameters(tmp_path):
    """Every added workflow validates through its owning command parser."""
    project = tmp_path / "project"
    data = project / "data"
    data.mkdir(parents=True)
    for name in ("proteins.csv", "sdrf.tsv", "de.csv"):
        (project / name).write_text("fixture\n", encoding="utf-8")
    for argv in MINIMUM_PERIPHERY_COMMANDS:
        prefix_length = 2 if argv[0] == "plot" else 1
        canonical = validate_and_canonicalize(argv, project)
        assert canonical[:prefix_length] == argv[:prefix_length]


def test_periphery_rejects_unsafe_or_duplicate_contrast_output_keys(tmp_path):
    """Contrast-derived filenames cannot escape or overwrite one another."""
    project = tmp_path / "project"
    project.mkdir()
    for name in ("proteins.csv", "sdrf.tsv", "de.csv"):
        (project / name).write_text("fixture\n", encoding="utf-8")

    base = [
        "interactive-report",
        "--protein-matrix",
        "proteins.csv",
        "--sdrf",
        "sdrf.tsv",
        "--output",
        "report.html",
    ]
    with pytest.raises(ValueError, match="file-safe"):
        validate_and_canonicalize(
            [*base, "--contrast", "a/../../escape", "a", "b", "de.csv"],
            project,
        )
    with pytest.raises(ValueError, match="same path"):
        validate_and_canonicalize(
            [
                *base,
                "--contrast",
                "same",
                "a",
                "b",
                "de.csv",
                "--contrast",
                "same",
                "a",
                "b",
                "de.csv",
            ],
            project,
        )


def test_plot_de_requires_an_output_type(tmp_path):
    """Studio validation applies the command's cross-parameter contract."""
    project = tmp_path / "project"
    project.mkdir()
    for name in ("proteins.csv", "de.csv"):
        (project / name).write_text("fixture\n", encoding="utf-8")

    with pytest.raises(ValueError, match="select --volcano or --heatmap"):
        validate_and_canonicalize(
            [
                "plot",
                "de",
                "--protein-matrix",
                "proteins.csv",
                "--outdir",
                "plots",
                "--contrast",
                "comparison",
                "a",
                "b",
                "de.csv",
            ],
            project,
        )


def test_execute_command_routes_python_workflows(monkeypatch):
    """The Studio worker sends periphery workflows to their Python owner."""
    calls = []

    def fake_main(args, mode):
        calls.append((args, mode))
        return 0

    module = __import__("mokume.commands.de_plots", fromlist=["main"])
    monkeypatch.setattr(module, "main", fake_main)

    execute_command(["plot", "pca", "--protein-matrix", "proteins.csv"])

    assert calls == [(["--protein-matrix", "proteins.csv"], "pca")]
