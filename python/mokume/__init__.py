"""
mokume - A comprehensive proteomics quantification library.

This package provides tools for processing and analyzing proteomics data using
multiple protein quantification methods including piBAQ (paralog-aware iBAQ),
Top3, TopN, and MaxLFQ.
"""

import importlib.metadata
import warnings

from mokume._lazy import module_api
from mokume.core.logging_config import initialize_logging

# Suppress numpy matrix deprecation warning
warnings.filterwarnings(
    "ignore", category=PendingDeprecationWarning, module="numpy.matrixlib.defmatrix"
)
# Suppress pyopenms false-positive OPENMS_DATA_PATH warning
warnings.filterwarnings("ignore", message=".*OPENMS_DATA_PATH.*")

# `mokume-py` (pure Python) and `mokume` (Rust kernel) both install the `mokume`
# import package, so pip silently overwrites files when both are present. Warn so
# the user keeps only one.
try:
    importlib.metadata.distribution("mokume")
    warnings.warn(
        "Both 'mokume-py' (pure-Python) and 'mokume' (Rust kernel) are installed; "
        "they share the 'mokume' import name and overwrite each other's files. "
        "Keep only one: uninstall the other (`pip uninstall mokume`).",
        RuntimeWarning,
        stacklevel=2,
    )
except importlib.metadata.PackageNotFoundError:
    pass

__version__ = "0.2.0"

# Initialize logging with default settings
# Users can override these settings by calling initialize_logging with their own settings
initialize_logging()

_LAZY_EXPORTS = {
    "is_directlfq_available": (
        "mokume.quantification.directlfq",
        "is_directlfq_available",
    ),
    "is_batch_correction_available": (
        "mokume.postprocessing.batch_correction",
        "is_batch_correction_available",
    ),
}


__getattr__, __dir__ = module_api(_LAZY_EXPORTS, globals(), __name__)


__all__ = [
    "__version__",
    "is_directlfq_available",
    "is_batch_correction_available",
]
