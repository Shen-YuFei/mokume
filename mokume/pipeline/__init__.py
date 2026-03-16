"""
Unified pipelines for proteomics data processing.

This module provides high-level pipelines that combine multiple processing
steps into single, easy-to-use functions.
"""

from mokume.pipeline.config import (
    PipelineConfig,
    InputConfig,
    FilterConfig,
    NormalizationConfig,
    QuantificationConfig,
    IRSConfig,
    BatchCorrectionConfig,
    DEConfig,
    OutputConfig,
)
from mokume.pipeline.features_to_proteins import (
    QuantificationPipeline,
    features_to_proteins,
)
from mokume.pipeline.stages import (
    LoadingStage,
    NormalizationStage,
    QuantificationStage,
    PostprocessingStage,
)

__all__ = [
    "PipelineConfig",
    "InputConfig",
    "FilterConfig",
    "NormalizationConfig",
    "QuantificationConfig",
    "IRSConfig",
    "BatchCorrectionConfig",
    "DEConfig",
    "OutputConfig",
    "QuantificationPipeline",
    "features_to_proteins",
    "LoadingStage",
    "NormalizationStage",
    "QuantificationStage",
    "PostprocessingStage",
]
