"""
Statistical analysis module for the mokume package.

Provides differential expression analysis and related statistical methods.
"""

from mokume.analysis.deqms import run_deqms, spectra_count_ebayes
from mokume.analysis.differential_expression import DifferentialExpression
from mokume.analysis.limrots import run_limrots
from mokume.analysis.proda import run_proda

__all__ = [
    "DifferentialExpression",
    "run_deqms",
    "run_limrots",
    "run_proda",
    "spectra_count_ebayes",
]
