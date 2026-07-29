"""Single-source mapping from :class:`PipelineConfig` to Rust-kernel argv.

The Rust kernel (``mokume._mokume.run``) reuses the ``features2proteins``
subcommand's argument parser. When the pipeline is asked to run with
``config.runtime.backend == "rust"`` the Python side must hand the kernel the
exact same argument vector a user would type on the command line. Keeping that
translation in one pure function means the flow module (and its tests) never
have to duplicate the flag spelling, and the mapping stays unit-testable
without the compiled extension present.

The flag names and their allowed values were verified against
``mokume features2proteins --help`` from the installed ``mokume-rs`` wheel.
"""

import re
from dataclasses import fields
from enum import Enum

from typing import List

from mokume.pipeline.config import IRSConfig, PipelineConfig

# Read the default off the dataclass rather than repeating the literal, so the
# "only forward a non-default regex" rule cannot drift from the actual default.
_DEFAULT_REFERENCE_REGEX = next(
    field.default for field in fields(IRSConfig) if field.name == "reference_regex"
)


class HybridFieldOwner(str, Enum):
    """Effective owner of a leaf :class:`PipelineConfig` field."""

    RUST_KERNEL = "rust-kernel"
    PYTHON_POSTPROCESS = "python-postprocess"
    SHARED_CONTEXT = "shared-context"
    HYBRID_CONTROL = "hybrid-control"
    UNSUPPORTED = "unsupported"


# This inventory is executable documentation. ``test_cli_args`` derives every
# leaf field from the dataclasses and fails when a field is added without an
# ownership decision.
HYBRID_FIELD_OWNERS = {
    # Inputs.
    "input.parquet": HybridFieldOwner.RUST_KERNEL,
    "input.sdrf": HybridFieldOwner.SHARED_CONTEXT,
    "input.fasta_file": HybridFieldOwner.RUST_KERNEL,
    "input.qpx_dir": HybridFieldOwner.PYTHON_POSTPROCESS,
    # Filtering and normalization performed by the Rust kernel.
    "filtering.min_aa": HybridFieldOwner.RUST_KERNEL,
    "filtering.min_unique_peptides": HybridFieldOwner.RUST_KERNEL,
    "filtering.remove_contaminants": HybridFieldOwner.RUST_KERNEL,
    "normalization.run_method": HybridFieldOwner.RUST_KERNEL,
    "normalization.sample_method": HybridFieldOwner.RUST_KERNEL,
    "normalization.proteins_file": HybridFieldOwner.RUST_KERNEL,
    # Quantification. The method is also used by Python postprocessing to
    # select scale-sensitive behavior such as ratio handling.
    "quantification.method": HybridFieldOwner.SHARED_CONTEXT,
    "quantification.ion_alignment": HybridFieldOwner.UNSUPPORTED,
    "quantification.coverage_threshold": HybridFieldOwner.PYTHON_POSTPROCESS,
    "quantification.ratio_fraction_merge": HybridFieldOwner.RUST_KERNEL,
    "quantification.directlfq_num_cores": HybridFieldOwner.RUST_KERNEL,
    "quantification.directlfq_min_nonan": HybridFieldOwner.RUST_KERNEL,
    "quantification.directlfq_num_samples_quadratic": HybridFieldOwner.RUST_KERNEL,
    "quantification.ibaq_enzyme": HybridFieldOwner.RUST_KERNEL,
    "quantification.ibaq_max_aa": HybridFieldOwner.RUST_KERNEL,
    "quantification.ibaq_min_shared": HybridFieldOwner.RUST_KERNEL,
    "quantification.ibaq_families_yaml": HybridFieldOwner.RUST_KERNEL,
    "quantification.ibaq_min_anchors": HybridFieldOwner.RUST_KERNEL,
    "quantification.ibaq_high_anchor_threshold": HybridFieldOwner.RUST_KERNEL,
    # Python-owned postprocessing.
    "irs.enabled": HybridFieldOwner.PYTHON_POSTPROCESS,
    # Ratio aggregation consumes reference selection in Rust even when the IRS
    # postprocessing stage itself remains disabled there.
    "irs.reference_samples": HybridFieldOwner.SHARED_CONTEXT,
    "irs.sdrf_column": HybridFieldOwner.PYTHON_POSTPROCESS,
    "irs.sdrf_values": HybridFieldOwner.PYTHON_POSTPROCESS,
    "irs.reference_regex": HybridFieldOwner.SHARED_CONTEXT,
    "irs.stat": HybridFieldOwner.PYTHON_POSTPROCESS,
    "irs.remove_reference": HybridFieldOwner.PYTHON_POSTPROCESS,
    "batch.enabled": HybridFieldOwner.PYTHON_POSTPROCESS,
    "batch.method": HybridFieldOwner.PYTHON_POSTPROCESS,
    "batch.column": HybridFieldOwner.PYTHON_POSTPROCESS,
    "batch.covariates": HybridFieldOwner.PYTHON_POSTPROCESS,
    "batch.parametric": HybridFieldOwner.PYTHON_POSTPROCESS,
    "batch.mean_only": HybridFieldOwner.PYTHON_POSTPROCESS,
    "batch.ref_batch": HybridFieldOwner.PYTHON_POSTPROCESS,
    "imputation.enabled": HybridFieldOwner.PYTHON_POSTPROCESS,
    "imputation.method": HybridFieldOwner.PYTHON_POSTPROCESS,
    "imputation.quantile": HybridFieldOwner.PYTHON_POSTPROCESS,
    "imputation.shift": HybridFieldOwner.PYTHON_POSTPROCESS,
    "imputation.scale": HybridFieldOwner.PYTHON_POSTPROCESS,
    "imputation.n_neighbors": HybridFieldOwner.PYTHON_POSTPROCESS,
    "de.enabled": HybridFieldOwner.PYTHON_POSTPROCESS,
    "de.contrasts": HybridFieldOwner.PYTHON_POSTPROCESS,
    "de.method": HybridFieldOwner.PYTHON_POSTPROCESS,
    "de.log2fc_threshold": HybridFieldOwner.PYTHON_POSTPROCESS,
    "de.fdr_threshold": HybridFieldOwner.PYTHON_POSTPROCESS,
    "de.fdr_method": HybridFieldOwner.PYTHON_POSTPROCESS,
    "de.output": HybridFieldOwner.PYTHON_POSTPROCESS,
    "de.ensemble_methods": HybridFieldOwner.PYTHON_POSTPROCESS,
    "de.ensemble_min_k": HybridFieldOwner.PYTHON_POSTPROCESS,
    # Quantification exports are produced by Rust. Presentation and AnnData
    # remain in the shared Python postprocessing path.
    "output.export_peptides": HybridFieldOwner.RUST_KERNEL,
    "output.export_ions": HybridFieldOwner.RUST_KERNEL,
    "output.plot_dir": HybridFieldOwner.PYTHON_POSTPROCESS,
    "output.plot_volcano": HybridFieldOwner.PYTHON_POSTPROCESS,
    "output.plot_heatmap": HybridFieldOwner.PYTHON_POSTPROCESS,
    "output.plot_pca": HybridFieldOwner.PYTHON_POSTPROCESS,
    "output.highlight_genes": HybridFieldOwner.PYTHON_POSTPROCESS,
    "output.interactive_report": HybridFieldOwner.PYTHON_POSTPROCESS,
    "output.report_output": HybridFieldOwner.PYTHON_POSTPROCESS,
    "output.export_anndata": HybridFieldOwner.PYTHON_POSTPROCESS,
    "output.anndata_view": HybridFieldOwner.PYTHON_POSTPROCESS,
    # The backend selects this adapter. Rust does not currently enforce the
    # documented DuckDB resource semantics in-process, so accepting either
    # resource field would still silently change the request.
    "runtime.duckdb_memory": HybridFieldOwner.UNSUPPORTED,
    "runtime.duckdb_threads": HybridFieldOwner.UNSUPPORTED,
    "runtime.backend": HybridFieldOwner.HYBRID_CONTROL,
}


def _filtering_argv(config: PipelineConfig) -> List[str]:
    """Peptide/protein filtering flags."""
    argv = [
        "--min-aa",
        str(config.filtering.min_aa),
        "--min-unique",
        str(config.filtering.min_unique_peptides),
    ]
    argv.append(
        "--remove-contaminants"
        if config.filtering.remove_contaminants
        else "--keep-contaminants"
    )
    return argv


def _normalization_argv(config: PipelineConfig) -> List[str]:
    """Run/sample normalization flags and the optional protein subset."""
    argv = [
        "--run-normalization",
        config.normalization.run_method.lower(),
        "--sample-normalization",
        config.normalization.sample_method.lower(),
    ]
    if config.normalization.proteins_file is not None:
        argv += ["--normalization-proteins", str(config.normalization.proteins_file)]
    return argv


def _ratio_reference_argv(config: PipelineConfig) -> List[str]:
    """Reference-sample context for ``ratio`` quantification.

    An explicitly empty list falls through to SDRF/regex detection rather than
    becoming an empty Rust override, matching Python ratio loading.
    """
    argv: List[str] = []
    samples = [str(sample) for sample in config.irs.reference_samples or []]
    if samples:
        # ``--irs-reference-sample`` (repeatable) is newer than the comma-joined
        # ``--irs-reference-samples``, and mokume / mokume-rs are installed
        # separately, so a current Python package can sit on an older wheel.
        # Prefer the long-standing plural form and only reach for the new flag
        # when a name actually contains a comma -- the case it exists for.
        if any("," in sample for sample in samples):
            for sample in samples:
                argv += ["--irs-reference-sample", sample]
        else:
            argv += ["--irs-reference-samples", ",".join(samples)]
    if config.irs.reference_regex != _DEFAULT_REFERENCE_REGEX:
        argv += ["--irs-reference-regex", config.irs.reference_regex]
    return argv


def validate_hybrid_config(config: PipelineConfig) -> None:
    """Reject settings the in-process Rust profile cannot honor."""
    if config.runtime.duckdb_memory is not None:
        raise ValueError(
            "The Rust backend cannot enforce runtime.duckdb_memory; "
            "use backend='python' or an external process memory limit."
        )
    if config.runtime.duckdb_threads is not None:
        raise ValueError(
            "The Rust backend cannot reliably enforce per-call "
            "runtime.duckdb_threads; use backend='python' until the Rust "
            "kernel uses an operation-scoped thread pool."
        )

    quant = config.quantification
    if (
        quant.ion_alignment is not None
        and quant.ion_alignment.strip().lower() != "none"
    ):
        raise ValueError(
            "The Rust backend does not support quantification.ion_alignment="
            f"{quant.ion_alignment!r}; only None or 'none' is supported."
        )


def build_features2proteins_argv(config: PipelineConfig, output_path: str) -> List[str]:
    """Translate a :class:`PipelineConfig` into a ``features2proteins`` argv.

    Parameters
    ----------
    config : PipelineConfig
        The pipeline configuration to translate.
    output_path : str
        Path the kernel should write its wide protein matrix to. Emitted via
        ``--output`` and forced to ``--output-format python-compatible`` so the
        CSV matches the pure-Python flow's ``ProteinName`` + sample-column
        layout.

    Returns
    -------
    list[str]
        Argument vector with ``"features2proteins"`` as ``argv[0]``. Only flags
        whose config value is meaningful (not ``None``) are emitted; unknown
        flags are never produced.

    Raises
    ------
    ValueError
        If the configuration asks for something the hybrid path cannot honor
        (see :func:`validate_hybrid_config`), or if a reference sample name
        cannot be expressed with the installed kernel's flags.
    """
    validate_hybrid_config(config)
    argv: List[str] = ["features2proteins"]

    # Inputs and outputs.
    argv += ["--parquet", config.input.parquet]
    if config.input.sdrf is not None:
        argv += ["--sdrf", config.input.sdrf]
    if config.input.fasta_file is not None:
        argv += ["--fasta", config.input.fasta_file]
    argv += ["--output", output_path]
    argv += ["--output-format", "python-compatible"]

    # Quantification method. The kernel only understands ``top3`` and the
    # generic ``topn`` (with ``--topn N``); a ``top{N}`` request for any other N
    # is translated to ``topn`` + ``--topn N``. Every other method forwards
    # verbatim (lower-cased).
    method = config.quantification.method.lower()
    topn_match = re.fullmatch(r"top(\d+)", method)
    if topn_match and topn_match.group(1) != "3":
        argv += ["--quant-method", "topn"]
        argv += ["--topn", topn_match.group(1)]
    else:
        argv += ["--quant-method", method]

    argv += _filtering_argv(config)
    argv += _normalization_argv(config)

    # piBAQ knobs.
    quant = config.quantification
    # ``mean`` is already the kernel default. Keep the default argv stable and
    # only forward a setting that changes ratio quantification semantics.
    if quant.ratio_fraction_merge.lower() != "mean":
        argv += ["--ratio-fraction-merge", quant.ratio_fraction_merge.lower()]
    if method == "ratio":
        argv += _ratio_reference_argv(config)
    argv += ["--ibaq-enzyme", quant.ibaq_enzyme]
    argv += ["--ibaq-max-aa", str(quant.ibaq_max_aa)]
    argv += ["--ibaq-min-shared", str(quant.ibaq_min_shared)]
    if quant.ibaq_families_yaml is not None:
        argv += ["--ibaq-families", quant.ibaq_families_yaml]
    argv += ["--ibaq-min-anchors", str(quant.ibaq_min_anchors)]
    argv += ["--ibaq-high-anchor-threshold", str(quant.ibaq_high_anchor_threshold)]

    # DirectLFQ knobs.
    if quant.directlfq_num_cores is not None:
        argv += ["--directlfq-cores", str(quant.directlfq_num_cores)]
    argv += ["--directlfq-min-nonan", str(quant.directlfq_min_nonan)]
    argv += [
        "--directlfq-num-samples-quadratic",
        str(quant.directlfq_num_samples_quadratic),
    ]

    # Quantification-stage exports belong to the kernel. Python-owned plots,
    # reports, and AnnData export stay in ``runner._postprocess``.
    if config.output.export_peptides is not None:
        argv += ["--export-peptides", str(config.output.export_peptides)]
    if config.output.export_ions is not None:
        argv += ["--export-ions", str(config.output.export_ions)]

    # Batch correction is intentionally NOT emitted: under the Python-owned
    # post-processing model the kernel only loads, filters, normalizes and
    # quantifies. ComBat runs in ``runner._postprocess`` for both backends, so
    # emitting ``--batch-correction`` here would double-correct the matrix.

    return argv
