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
from typing import List

from mokume.pipeline.config import PipelineConfig


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
        flags are never produced. The function is pure and side-effect free.
    """
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

    # Filtering.
    argv += ["--min-aa", str(config.filtering.min_aa)]
    argv += ["--min-unique", str(config.filtering.min_unique_peptides)]
    if config.filtering.remove_contaminants:
        argv += ["--remove-contaminants"]
    else:
        argv += ["--keep-contaminants"]

    # Normalization.
    argv += ["--run-normalization", config.normalization.run_method.lower()]
    argv += ["--sample-normalization", config.normalization.sample_method.lower()]

    # piBAQ knobs.
    quant = config.quantification
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

    # Batch correction is intentionally NOT emitted: under the Python-owned
    # post-processing model the kernel only loads, filters, normalizes and
    # quantifies. ComBat runs in ``runner._postprocess`` for both backends, so
    # emitting ``--batch-correction`` here would double-correct the matrix.

    return argv
