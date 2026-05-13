"""State definitions for the agentic optimization workflow."""

from dataclasses import dataclass, field
from typing import Any


@dataclass
class CandidateConfig:
    """A candidate analysis configuration to test."""

    name: str
    de_method: str = "deqms"
    fdr_method: str = "bh"
    normalization: str = "none"
    imputation: str = "none"
    log2fc_threshold: float = 0.5
    ensemble: str = "none"  # "none" or comma-separated methods
    ensemble_k: int = 2  # min agreement for top-k consensus
    reasoning: str = ""
    expected_outcome: str = ""

    def to_dict(self) -> dict:
        """Serialize to dictionary."""
        return {
            "name": self.name,
            "de_method": self.de_method,
            "fdr_method": self.fdr_method,
            "normalization": self.normalization,
            "imputation": self.imputation,
            "log2fc_threshold": self.log2fc_threshold,
            "ensemble": self.ensemble,
            "ensemble_k": self.ensemble_k,
            "reasoning": self.reasoning,
            "expected_outcome": self.expected_outcome,
        }


@dataclass
class EvaluationResult:
    """Evaluation metrics for a single configuration run."""

    config_name: str
    config: dict
    # Ground truth metrics (None if no ground truth)
    tp: int | None = None
    fp: int | None = None
    fn: int | None = None
    auc: float | None = None
    sensitivity: float | None = None
    specificity: float | None = None
    # Universal metrics
    n_de_up: int = 0
    n_de_down: int = 0
    median_cv: float | None = None
    missing_rate: float = 0.0
    # Composite score
    score: float = 0.0


@dataclass
class ReflectionResult:
    """Output from the reflection step."""

    converged: bool
    next_configs: list[CandidateConfig] = field(default_factory=list)
    analysis: str = ""
    adjustments: list[str] = field(default_factory=list)


@dataclass
class RoundResult:
    """Results from a single optimization round."""

    round_num: int
    configs: list[CandidateConfig]
    results: list[EvaluationResult]
    best_config_name: str = ""
    reflection: ReflectionResult | None = None


@dataclass
class AuditEntry:
    """A single entry in the audit trail."""

    step: str
    round_num: int
    data: Any = None


@dataclass
class AgenticState:
    """Full state of the optimization workflow."""

    rounds: list[RoundResult] = field(default_factory=list)
    total_experiments: int = 0
    converged: bool = False
    best_config: CandidateConfig | None = None
    best_score: float = 0.0
    audit_trail: list[AuditEntry] = field(default_factory=list)
