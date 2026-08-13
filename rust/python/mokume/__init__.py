"""mokume — a proteomics quantification toolkit.

The compute-heavy pipeline (expression matrix, normalization, imputation,
differential expression, batch correction) lives in a Rust kernel, exposed here
through the compiled ``mokume._mokume`` extension that maturin builds from
``crates/mokume-py``. The Python periphery lives in :mod:`mokume.commands` as
first-class modules. Plotting and reporting consume kernel tables, TissueMap
derives downstream atlas outputs from QPX data, and explicit fallbacks compute
operations the kernel does not provide.

This is the PyO3/maturin layout used by projects such as polars and
pydantic-core (Python imports a compiled Rust extension). The compute commands
run *in-process* (no subprocess) through the same clap parsing + dispatch the
standalone CLI binary uses, so the flag handling stays single-sourced in Rust.

Each wrapper maps keyword arguments to CLI flags (``key=value`` -> ``--key
value`` with ``_`` rewritten to ``-``; ``key=True`` -> ``--key``; a list value
repeats the flag once per item; ``None`` / ``False`` are skipped). For full
control, call :func:`run` with an explicit argument list.
"""

import importlib
import importlib.metadata
import warnings

from mokume._mokume import run as _run
from mokume._mokume import version

# `mokume-rs` (this Rust kernel) and `mokume` (pure Python) both install the
# `mokume` import package, so pip silently overwrites files when both are
# present. Warn so the user keeps only one.
try:
    importlib.metadata.distribution("mokume")
    warnings.warn(
        "Both 'mokume-rs' (Rust kernel) and 'mokume' (pure-Python) are installed; "
        "they share the 'mokume' import name and overwrite each other's files. "
        "Keep only one: uninstall the other (`pip uninstall mokume`).",
        RuntimeWarning,
        stacklevel=2,
    )
except importlib.metadata.PackageNotFoundError:
    pass

__all__ = [
    "version",
    "run",
    "features2peptides",
    "features2proteins",
    "peptides2protein",
    "correct_batches",
    "tsne_visualization",
    "tissuemap",
    "peptides2protein_qc",
    "peptides2protein_ibaq",
    "de_plots",
    "interactive_report",
    "qc_report",
    "workflow_comparison",
    "impute",
]
__version__ = version()

_BOOTSTRAPPED = []


def _bootstrap():
    """Apply the periphery's warning filters + logging once, lazily.

    Kept out of import time so a pure ``import mokume`` for the compute path does
    not pull in the logging / warnings setup the plotting commands want.
    """
    if _BOOTSTRAPPED:
        return

    # Lazy: importing mokume.core would pull in the periphery's heavier modules.
    from mokume.core.logging_config import initialize_logging

    # Suppress the numpy matrix deprecation and the pyopenms OPENMS_DATA_PATH
    # false positive, matching the upstream package import behaviour.
    warnings.filterwarnings(
        "ignore",
        category=PendingDeprecationWarning,
        module="numpy.matrixlib.defmatrix",
    )
    warnings.filterwarnings("ignore", message=".*OPENMS_DATA_PATH.*")
    initialize_logging()
    _BOOTSTRAPPED.append(True)


def _flags(kwargs):
    """Translate ``key=value`` keyword arguments into a CLI flag list."""
    args = []
    for key, value in kwargs.items():
        if value is None or value is False:
            continue
        flag = "--" + key.replace("_", "-")
        if value is True:
            args.append(flag)
        elif isinstance(value, (list, tuple)):
            for item in value:
                args.append(flag)
                args.append(str(item))
        else:
            args.append(flag)
            args.append(str(value))
    return args


def _build_args(command, kwargs):
    """Prefix the subcommand name to the flags built from ``kwargs``."""
    return [command, *_flags(kwargs)]


def run(args):
    """Run a mokume subcommand in-process from an explicit argument list.

    Example::

        mokume.run(["features2proteins", "--parquet", "x.parquet",
                    "--output", "y.csv"])
    """
    _run(list(args))


def features2peptides(**kwargs):
    """Run ``features2peptides`` (feature parquet -> peptide-level output)."""
    _run(_build_args("features2peptides", kwargs))


def features2proteins(**kwargs):
    """Run ``features2proteins`` (feature parquet -> protein matrix, optional DE)."""
    _run(_build_args("features2proteins", kwargs))


def peptides2protein(**kwargs):
    """Run ``peptides2protein`` (peptide-level input -> protein quantities)."""
    _run(_build_args("peptides2protein", kwargs))


def correct_batches(**kwargs):
    """Run ``correct-batches`` (ComBat batch-effect correction on iBAQ output)."""
    _run(_build_args("correct-batches", kwargs))


def _run_command(module_name, argv):
    """Run a :mod:`mokume.commands` periphery command, raising on failure.

    The command modules signal failure with a non-zero return or ``SystemExit``;
    both are surfaced here as a ``RuntimeError`` so library callers get a normal
    exception instead of an interpreter exit.
    """
    _bootstrap()
    module = importlib.import_module(f"mokume.commands.{module_name}")
    try:
        code = module.main(list(argv))
    except SystemExit as exc:
        code = exc.code
        if isinstance(code, str):
            raise RuntimeError(code) from None
        if code not in (0, None):
            raise RuntimeError(f"{module_name} failed with exit code {code}") from None
        return
    if code not in (0, None):
        raise RuntimeError(f"{module_name} failed with exit code {code}")


def tsne_visualization(**kwargs):
    """Render the t-SNE plot for a folder of protein files (``plotting`` extra)."""
    _run_command("visualize", _flags(kwargs))


def tissuemap(**kwargs):
    """Run the per-dataset tissue proteome analysis (``tissuemap`` extra)."""
    _run_command("tissuemap", _flags(kwargs))


def peptides2protein_qc(**kwargs):
    """Render the iBAQ QC report from a protein table (``plotting`` extra)."""
    _run_command("peptides2protein_qc", _flags(kwargs))


def peptides2protein_ibaq(**kwargs):
    """Compute iBAQ in pure Python for non-Rust-ported enzymes (``ibaq`` extra)."""
    _run_command("peptides2protein_ibaq", _flags(kwargs))


def de_plots(args):
    """Render features2proteins DE plots from kernel-written CSVs.

    Takes an explicit argument list (the per-contrast ``--contrast KEY A B CSV``
    flag repeats, which keyword arguments cannot express). See
    ``python -m mokume.commands.de_plots --help``.
    """
    _run_command("de_plots", args)


def interactive_report(args):
    """Render the features2proteins interactive HTML report from kernel CSVs.

    Takes an explicit argument list (see ``de_plots`` for why) — run
    ``python -m mokume.commands.interactive_report --help`` for the flags.
    """
    _run_command("interactive_report", args)


def qc_report(
    protein_matrix,
    sdrf=None,
    output="qc_report.html",
    de_results=None,
    title="QC Report",
    is_log2=False,
):
    """Build a single-matrix QC HTML report (``analysis`` extra).

    Computes PCA / t-SNE / silhouette / variance decomposition / CV / missing-value
    and DE-quality metrics from a protein matrix CSV and writes an interactive HTML
    report. ``sdrf`` supplies the sample -> condition grouping; ``de_results`` is an
    optional DE result CSV. Returns the output path. For volcano gene-highlighting
    call :func:`mokume.reports.qc_report.generate_qc_report` directly.
    """
    _bootstrap()
    import pandas as pd

    from mokume.reports.qc_report import generate_qc_report

    protein_df = pd.read_csv(protein_matrix)
    sample_to_condition = {}
    if sdrf:
        from mokume.normalization.irs import detect_condition_from_sdrf

        sample_to_condition = detect_condition_from_sdrf(sdrf)
    de_df = pd.read_csv(de_results) if de_results else None
    return generate_qc_report(
        protein_df,
        sample_to_condition,
        de_results=de_df,
        output_html=str(output),
        title=title,
        is_log2=is_log2,
    )


def workflow_comparison(
    workflows,
    output="workflow_comparison.html",
    title="Workflow Comparison",
    marker_genes=None,
):
    """Build an HTML report comparing several quantification workflows (``analysis`` extra).

    ``workflows`` is a list of dicts, one per workflow, each with ``name`` and
    either ``protein_df`` (a DataFrame) or ``protein_matrix`` (a CSV path); plus
    optional ``de_results`` (DataFrame or CSV path), ``sample_to_condition`` (dict)
    or ``sdrf`` (path), and ``is_log2``. Returns the output path.
    """
    _bootstrap()
    import pandas as pd

    from mokume.reports.workflow_comparison import generate_comparison_report

    built = []
    for workflow in workflows:
        entry = {"name": workflow["name"], "is_log2": workflow.get("is_log2", False)}
        protein = workflow.get("protein_df")
        entry["protein_df"] = (
            protein if protein is not None else pd.read_csv(workflow["protein_matrix"])
        )
        de_results = workflow.get("de_results")
        if de_results is not None:
            entry["de_results"] = (
                de_results
                if hasattr(de_results, "columns")
                else pd.read_csv(de_results)
            )
        condition = workflow.get("sample_to_condition")
        if condition is None and workflow.get("sdrf"):
            from mokume.normalization.irs import detect_condition_from_sdrf

            condition = detect_condition_from_sdrf(workflow["sdrf"])
        entry["sample_to_condition"] = condition or {}
        built.append(entry)
    return generate_comparison_report(
        built, output_html=str(output), title=title, marker_genes=marker_genes
    )


def impute(matrix, method, output=None, **kwargs):
    """Impute missing values with the pure-Python imputers (``analysis`` extra).

    Reaches the imputers the Rust kernel does not reproduce (``missforest`` wraps
    a scikit-learn estimator), plus every other supported method. ``matrix`` is a
    wide protein matrix CSV or DataFrame; writes ``output`` if given and returns
    the imputed DataFrame.
    """
    _bootstrap()
    import pandas as pd

    from mokume.imputation import impute_missing_values

    frame = matrix if hasattr(matrix, "columns") else pd.read_csv(matrix, index_col=0)
    result = impute_missing_values(frame, method=method, **kwargs)
    if output:
        result.to_csv(output)
    return result
