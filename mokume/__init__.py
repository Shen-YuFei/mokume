"""
mokume - A comprehensive proteomics quantification library.

This package provides tools for processing and analyzing proteomics data using
multiple protein quantification methods including iBAQ (intensity-based absolute
quantification), Top3, TopN, and MaxLFQ.
"""

import warnings

# Suppress numpy matrix deprecation warning
warnings.filterwarnings(
    "ignore", category=PendingDeprecationWarning, module="numpy.matrixlib.defmatrix"
)
# Suppress pyopenms false-positive OPENMS_DATA_PATH warning
warnings.filterwarnings("ignore", message=".*OPENMS_DATA_PATH.*")

__version__ = "0.1.0"

# Import logging configuration
from mokume.core.logging_config import initialize_logging

# Initialize logging with default settings
# Users can override these settings by calling initialize_logging with their own settings
initialize_logging()

# Availability checks for optional dependencies
from mokume.quantification.directlfq import is_directlfq_available
from mokume.postprocessing.batch_correction import is_batch_correction_available

__all__ = [
    "__version__",
    "is_directlfq_available",
    "is_batch_correction_available",
]
