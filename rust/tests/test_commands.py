"""The relocated periphery command modules import and expose ``main(argv)``."""

import importlib
import sys
from types import SimpleNamespace

import pandas as pd
import pytest

COMMAND_MODULES = [
    "visualize",
    "tissuemap",
    "de_plots",
    "interactive_report",
    "peptides2protein_qc",
    "peptides2protein_pibaq",
]


@pytest.mark.parametrize("name", COMMAND_MODULES)
def test_command_module_imports_with_main(name):
    """Each relocated periphery command must expose a callable entrypoint."""
    module = importlib.import_module(f"mokume.commands.{name}")
    assert callable(getattr(module, "main", None))


def test_pibaq_qc_wrapper_accepts_global_options_before_subcommand(monkeypatch):
    """Global options must not prevent post-command piBAQ QC rendering."""
    observed = {}
    package = SimpleNamespace(
        _run_cli=lambda _args: 0,
        peptides2protein_qc=lambda **kwargs: observed.update(kwargs),
    )
    entrypoint = importlib.import_module("mokume.__main__")
    argv = [
        "mokume",
        "--log-level",
        "warn",
        "peptides2protein",
        "--method",
        "pibaq",
        "--verbose",
        "--normalize",
        "--qc-report",
        "qc.pdf",
        "-o",
        "proteins.tsv",
    ]
    monkeypatch.setattr(entrypoint.sys, "argv", argv)
    monkeypatch.setattr(entrypoint.importlib, "import_module", lambda _name: package)

    with pytest.raises(SystemExit) as exc_info:
        entrypoint.main()

    assert exc_info.value.code == 0
    assert observed == {
        "protein_table": "proteins.tsv",
        "qc_report": "qc.pdf",
        "plot_column": "PiBAQPpb",
        "tpa": False,
        "ruler": False,
    }


def test_pibaq_qc_wrapper_accepts_equals_syntax_without_verbose():
    """An explicit QC path requests rendering without a separate verbose flag."""
    observed = {}
    package = SimpleNamespace(
        peptides2protein_qc=lambda **kwargs: observed.update(kwargs)
    )
    entrypoint = importlib.import_module("mokume.__main__")

    entrypoint._render_requested_pibaq_qc(
        [
            "peptides2protein",
            "--method=pibaq",
            "--qc-report=qc.pdf",
            "--output=proteins.tsv",
        ],
        package,
    )

    assert observed["qc_report"] == "qc.pdf"


def test_python_compute_wrapper_runs_requested_pibaq_qc(monkeypatch):
    """The in-process Python wrapper must not accept qc_report as a no-op."""
    observed = {}
    package = importlib.import_module("mokume")
    monkeypatch.setattr(package, "_prepare_pibaq_digest", lambda _args: None)
    monkeypatch.setattr(package, "_native_run", lambda _args: None)
    monkeypatch.setattr(
        package,
        "peptides2protein_qc",
        lambda **kwargs: observed.update(kwargs),
    )

    package.peptides2protein(
        method="pibaq",
        peptides="peptides.tsv",
        fasta="proteins.fasta",
        output="proteins.tsv",
        qc_report="qc.pdf",
    )

    assert observed["protein_table"] == "proteins.tsv"
    assert observed["qc_report"] == "qc.pdf"


def test_non_pibaq_verbose_does_not_render_pibaq_qc():
    """Verbose output for another method must not invoke the piBAQ renderer."""
    package = SimpleNamespace(
        peptides2protein_qc=lambda **_kwargs: pytest.fail("unexpected piBAQ QC")
    )
    entrypoint = importlib.import_module("mokume.__main__")

    entrypoint._render_requested_pibaq_qc(
        [
            "peptides2protein",
            "--method",
            "top3",
            "--verbose",
            "--output",
            "proteins.tsv",
        ],
        package,
    )


def test_pibaq_compat_command_qc_report_enables_rendering(monkeypatch):
    """The compatibility command treats an explicit QC path as a render request."""
    captured = {}
    module = importlib.import_module("mokume.commands.peptides2protein_pibaq")
    pibaq_module = SimpleNamespace(
        peptides_to_protein=lambda **kwargs: captured.update(kwargs)
    )
    monkeypatch.setattr(
        module.importlib,
        "import_module",
        lambda _name: pibaq_module,
    )

    result = module.main(
        [
            "--peptides",
            "peptides.tsv",
            "--fasta",
            "proteins.fasta",
            "--enzyme",
            "Trypsin",
            "--output",
            "proteins.tsv",
            "--qc-report",
            "qc.pdf",
        ]
    )

    assert result == 0
    assert captured["verbose"] is True
    assert captured["qc_report"] == "qc.pdf"


def test_tissuemap_config_values_are_not_overridden_by_parser_defaults(
    tmp_path, monkeypatch
):
    """Omitted CLI options must preserve values loaded from TissueMap YAML."""
    scan_dir = tmp_path / "datasets"
    scan_dir.mkdir()
    output_dir = tmp_path / "configured-output"
    config_path = tmp_path / "tissuemap.yaml"
    config_path.write_text(
        "\n".join(
            (
                "n_jobs: 3",
                "input:",
                f"  scan_dir: {scan_dir}",
                "output:",
                f"  output_dir: {output_dir}",
            )
        ),
        encoding="utf-8",
    )
    captured = {}

    def capture_pipeline(config):
        """Return a no-op pipeline while retaining its resolved configuration."""
        captured["config"] = config
        return SimpleNamespace(run=lambda: None)

    monkeypatch.setitem(
        sys.modules,
        "mokume.tissuemap.pipeline",
        SimpleNamespace(TissueMapPipeline=capture_pipeline),
    )
    module = importlib.import_module("mokume.commands.tissuemap")

    assert module.main(["--config", str(config_path)]) == 0

    config = captured["config"]
    assert config.n_jobs == 3
    assert config.input.scan_dir == scan_dir
    assert config.output.output_dir == output_dir


@pytest.mark.parametrize(
    "content,match",
    [
        ("unknown: true\n", "Unknown TissueMap config keys"),
        ("input:\n  unknown: true\n", "Unknown keys in 'input' section"),
    ],
)
def test_tissuemap_config_rejects_unknown_keys(tmp_path, content, match):
    """Unknown top-level and section keys must fail closed."""
    config_path = tmp_path / "tissuemap.yaml"
    config_path.write_text(content, encoding="utf-8")
    config_module = importlib.import_module("mokume.tissuemap.config")

    with pytest.raises(ValueError, match=match):
        config_module.load_config(config_path)


def test_visualize_reads_tab_delimited_protein_tables(tmp_path, monkeypatch):
    """Visualization must honor the delimiter of protein table inputs."""
    protein_path = tmp_path / "demo-proteins.tsv"
    protein_path.write_text(
        "ProteinName\tSampleID\tPiBAQLog\n"
        "P1\tS1\t1.0\n"
        "P2\tS1\t2.0\n"
        "P1\tS2\t1.1\n"
        "P2\tS2\t2.1\n",
        encoding="utf-8",
    )
    module = importlib.import_module("mokume.commands.visualize")
    plotting = importlib.import_module("mokume.plotting")

    monkeypatch.setattr(
        plotting,
        "compute_pca_with_plot",
        lambda frame, n_components: frame.iloc[:, :1],
    )
    monkeypatch.setattr(plotting, "plot_tsne", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        module,
        "_compute_tsne_compat",
        lambda frame: pd.DataFrame(
            {"tSNE1": [0.0] * len(frame), "tSNE2": [0.0] * len(frame)},
            index=frame.index,
        ),
    )

    assert module.main(["--folder", str(tmp_path)]) == 0


@pytest.mark.parametrize(
    "module_name,args,match",
    [
        (
            "de_plots",
            ["--protein-matrix", "proteins.csv", "--plot-dir", "plots"],
            "select --volcano, --heatmap, or --pca",
        ),
        (
            "de_plots",
            [
                "--protein-matrix",
                "proteins.csv",
                "--plot-dir",
                "plots",
                "--pca",
            ],
            "--heatmap/--pca require --sdrf",
        ),
        (
            "interactive_report",
            [
                "--protein-matrix",
                "proteins.csv",
                "--sdrf",
                "data.sdrf.tsv",
                "--report-output",
                "report.html",
                "--plot-dir",
                "plots",
            ],
            "choose --report-output or --plot-dir",
        ),
        (
            "peptides2protein_pibaq",
            [
                "--peptides",
                "peptides.tsv",
                "--fasta",
                "proteins.fasta",
                "--enzyme",
                "Trypsin",
                "--output",
                "proteins.tsv",
                "--ploidy",
                "2",
            ],
            "--ploidy/--organism/--cpc require --ruler",
        ),
    ],
)
def test_periphery_commands_reject_options_that_would_be_ignored(
    module_name, args, match
):
    """Periphery commands must reject options outside their active path."""
    module = importlib.import_module(f"mokume.commands.{module_name}")

    with pytest.raises(SystemExit, match=match):
        module.main(args)


def test_tissuemap_generate_config_rejects_run_options(tmp_path):
    """Config generation must reject execution-only overrides."""
    module = importlib.import_module("mokume.commands.tissuemap")

    code = module.main(
        [
            "--generate-config",
            str(tmp_path / "config.yaml"),
            "--n-jobs",
            "24",
        ]
    )

    assert code == 2
    assert not (tmp_path / "config.yaml").exists()
