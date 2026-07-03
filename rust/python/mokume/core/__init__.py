"""
Core modules for the mokume package.

This module provides fundamental utilities including constants, logging,
and asynchronous file writing capabilities.
"""

from mokume.core.constants import (
    PROTEIN_NAME,
    PEPTIDE_SEQUENCE,
    PEPTIDE_CANONICAL,
    PEPTIDE_CHARGE,
    CHANNEL,
    CONDITION,
    BIOREPLICATE,
    TECHREPLICATE,
    RUN,
    FRACTION,
    INTENSITY,
    NORM_INTENSITY,
    SAMPLE_ID,
    IBAQ,
    IBAQ_NORMALIZED,
    IBAQ_LOG,
    IBAQ_PPB,
    IBAQ_BEC,
    TPA,
    MOLECULARWEIGHT,
    COPYNUMBER,
    CONCENTRATION_NM,
    WEIGHT_NG,
    MOLES_NMOL,
    PARQUET_COLUMNS,
    parquet_map,
)
from mokume.core.logger import (
    get_logger,
    configure_logging,
    log_execution_time,
    log_function_call,
)

__all__ = [
    # Constants
    "PROTEIN_NAME",
    "PEPTIDE_SEQUENCE",
    "PEPTIDE_CANONICAL",
    "PEPTIDE_CHARGE",
    "CHANNEL",
    "CONDITION",
    "BIOREPLICATE",
    "TECHREPLICATE",
    "RUN",
    "FRACTION",
    "INTENSITY",
    "NORM_INTENSITY",
    "SAMPLE_ID",
    "IBAQ",
    "IBAQ_NORMALIZED",
    "IBAQ_LOG",
    "IBAQ_PPB",
    "IBAQ_BEC",
    "TPA",
    "MOLECULARWEIGHT",
    "COPYNUMBER",
    "CONCENTRATION_NM",
    "WEIGHT_NG",
    "MOLES_NMOL",
    "PARQUET_COLUMNS",
    "parquet_map",
    # Logger
    "get_logger",
    "configure_logging",
    "log_execution_time",
    "log_function_call",
]


# write_queue pulls in pyarrow, which only the parquet-writing compute path needs
# (the ``ibaq`` extra). Import it lazily so importing mokume.core for logging /
# constants / the plotting periphery does not require pyarrow. WriteCSVTask /
# WriteParquetTask stay reachable as ``mokume.core.WriteCSVTask`` via this hook,
# but are kept out of ``__all__`` since they are not eagerly bound.
def __getattr__(name):
    if name in ("WriteCSVTask", "WriteParquetTask"):
        from mokume.core import write_queue

        return getattr(write_queue, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
