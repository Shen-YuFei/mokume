"""
CLI commands for the mokume package.

This module provides Click commands for the mokume CLI.
"""

from mokume.commands.features2peptides import features2parquet
from mokume.commands.features2proteins import features2proteins
from mokume.commands.peptides2protein import peptides2protein
from mokume.commands.visualize import tsne_visualization
from mokume.commands.batch_correct import correct_batches

__all__ = [
    "features2parquet",
    "features2proteins",  # Unified pipeline (recommended)
    "peptides2protein",
    "tsne_visualization",
    "correct_batches",
]
