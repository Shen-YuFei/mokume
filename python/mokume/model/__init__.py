"""
Data models and enumerations for the mokume package.

This module provides enumerations and data classes for:
- Labeling types (TMT, iTRAQ, LFQ)
- Normalization methods
- Quantification methods
- Summarization strategies
- Organism metadata
"""

from mokume.model.labeling import (
    QuantificationCategory,
    IsobaricLabel,
    IsobaricLabelSpec,
    TMT6plex,
    TMT10plex,
    TMT11plex,
    TMT16plex,
    ITRAQ4plex,
    ITRAQ8plex,
)
from mokume.model.normalization import (
    FeatureNormalizationMethod,
    PeptideNormalizationMethod,
    # Preferred names (aliases)
    RunNormalizationMethod,
    SampleNormalizationMethod,
)
from mokume.model.organism import OrganismDescription
from mokume.model.quantification import QuantificationMethod
from mokume.model.summarization import SummarizationMethod
from mokume.model.batch_correction import (
    BatchDetectionMethod,
    BatchCorrectionConfig,
)

__all__ = [
    # Labeling
    "QuantificationCategory",
    "IsobaricLabel",
    "IsobaricLabelSpec",
    "TMT6plex",
    "TMT10plex",
    "TMT11plex",
    "TMT16plex",
    "ITRAQ4plex",
    "ITRAQ8plex",
    # Normalization
    "FeatureNormalizationMethod",
    "PeptideNormalizationMethod",
    # Normalization (preferred names)
    "RunNormalizationMethod",
    "SampleNormalizationMethod",
    # Organism
    "OrganismDescription",
    # Quantification
    "QuantificationMethod",
    # Summarization
    "SummarizationMethod",
    # Batch correction
    "BatchDetectionMethod",
    "BatchCorrectionConfig",
]
