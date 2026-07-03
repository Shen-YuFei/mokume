"""
mokume - A comprehensive proteomics quantification library.

This package provides tools for processing and analyzing proteomics data using
multiple protein quantification methods including iBAQ (intensity-based absolute
quantification), Top3, TopN, and MaxLFQ.
"""

import importlib.metadata
import logging
import warnings

from mokume.core.logging_config import initialize_logging
from mokume.postprocessing.batch_correction import is_batch_correction_available
from mokume.quantification.directlfq import is_directlfq_available

# Suppress numpy matrix deprecation warning
warnings.filterwarnings(
    "ignore", category=PendingDeprecationWarning, module="numpy.matrixlib.defmatrix"
)
# Suppress pyopenms false-positive OPENMS_DATA_PATH warning
warnings.filterwarnings("ignore", message=".*OPENMS_DATA_PATH.*")

# `mokume` (pure Python) and `mokume-rs` (Rust kernel) both install the `mokume`
# import package, so pip silently overwrites files when both are present. Warn so
# the user keeps only one.
try:
    importlib.metadata.distribution("mokume-rs")
    warnings.warn(
        "Both 'mokume' (pure-Python) and 'mokume-rs' (Rust kernel) are installed; "
        "they share the 'mokume' import name and overwrite each other's files. "
        "Keep only one: uninstall the other (`pip uninstall mokume-rs`).",
        RuntimeWarning,
        stacklevel=2,
    )
except importlib.metadata.PackageNotFoundError:
    pass

__version__ = "0.1.0"

# Initialize logging with default settings
# Users can override these settings by calling initialize_logging with their own settings
initialize_logging()

# Auto-load all built-in plugins so `import mokume` (or the first registry
# access) populates every registration group. `_plugins` imports each plugin
# module in isolation and swallows optional-dependency ImportErrors, so this is
# guarded once more here against any unexpected import-time failure to keep
# `import mokume` robust.
try:
    from mokume import _plugins as _plugins  # noqa: F401  (import for side effects)
except Exception:  # pragma: no cover - defensive; plugins are optional at import
    logging.getLogger(__name__).debug(
        "Plugin auto-loading failed; plugins can still be registered on demand.",
        exc_info=True,
    )

__all__ = [
    "__version__",
    "is_directlfq_available",
    "is_batch_correction_available",
]
