"""Rust-kernel quantification flow.

Routes ``run_pipeline`` to the compiled ``mokume-rs`` kernel when
``config.runtime.backend == "rust"``. The kernel runs the whole
features-to-proteins path in Rust and writes a wide protein x sample CSV; this
module reads that CSV back into a :class:`QpxDataset` whose ``proteins`` frame
has the SAME shape the pure-Python :mod:`mokume.pipeline.flows.standard` flow
produces (``ProteinName`` plus one column per sample), so the shared
``runner._postprocess`` stays backend-agnostic.

The extension is import-guarded: environments that ship only the pure-Python
package have no ``mokume._mokume`` and must still be able to
``import mokume.pipeline.flows.rust`` without error. The missing-kernel case is
surfaced only when :func:`run` is actually called.
"""

import os
import tempfile

import pandas as pd

from mokume.core.cli_args import build_features2proteins_argv
from mokume.core.dataset import QpxDataset
from mokume.core.logger import get_logger
from mokume.pipeline.config import PipelineConfig
from mokume.quantification.base import QuantificationMethod

try:
    from mokume._mokume import run as _kernel_run  # noqa: F401
except ImportError:
    _kernel_run = None

logger = get_logger("mokume.pipeline.flows.rust")


def run(method: QuantificationMethod, config: PipelineConfig) -> QpxDataset:
    """Execute the quantification flow through the Rust kernel.

    Parameters
    ----------
    method : QuantificationMethod
        The resolved quantification method (used for metadata only; the actual
        dispatch happens inside the kernel via ``--quant-method``).
    config : PipelineConfig
        Pipeline configuration.

    Returns
    -------
    QpxDataset
        Dataset with ``proteins`` populated as a wide ``ProteinName`` + sample
        matrix, matching the pure-Python standard flow.

    Raises
    ------
    RuntimeError
        If the ``mokume-rs`` extension (``mokume._mokume``) is not installed.
    """
    _ = method  # metadata only; kernel selects the method from the argv
    if _kernel_run is None:
        raise RuntimeError(
            "The Rust backend requires the mokume-rs wheel (mokume._mokume), "
            "which is not installed in this environment. Install mokume-rs or "
            'use backend="python".'
        )

    dataset = QpxDataset()

    # The kernel writes its wide protein matrix to a CSV; hand it a temp path.
    tmp_fd, tmp_out = tempfile.mkstemp(suffix=".csv")
    os.close(tmp_fd)
    try:
        argv = build_features2proteins_argv(config, tmp_out)
        logger.info("Running Rust kernel: features2proteins")
        logger.debug("Kernel argv: %s", argv)
        _kernel_run(argv)

        protein_df = pd.read_csv(tmp_out)
    finally:
        if os.path.exists(tmp_out):
            os.remove(tmp_out)

    dataset.proteins = protein_df
    dataset.record_step(
        "quantification",
        method=config.quantification.method,
        backend="rust",
        rows_out=len(protein_df),
    )
    logger.info("Rust quantification complete: %d proteins", len(protein_df))

    return dataset
