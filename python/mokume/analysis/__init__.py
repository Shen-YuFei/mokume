"""
Statistical analysis module for the mokume package.

Provides differential expression analysis and related statistical methods.
"""

from mokume.analysis.deqms import run_deqms
from mokume.analysis.differential_expression import DifferentialExpression
from mokume.analysis.ensemble import run_ensemble, combine_de_results
from mokume.analysis.limma import run_limma
from mokume.analysis.limrots import run_limrots
from mokume.analysis.proda import run_proda
from mokume.analysis.rots import run_rots

__all__ = [
    "DifferentialExpression",
    "run_deqms",
    "run_limma",
    "run_limrots",
    "run_proda",
    "run_rots",
    "run_ensemble",
    "combine_de_results",
]
