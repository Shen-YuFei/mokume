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

    def signature(self) -> str:
        """Return a stable identifier for dedup (excludes free-text fields)."""
        return "|".join(
            [
                self.de_method,
                self.fdr_method,
                self.normalization,
                self.imputation,
                f"fc{self.log2fc_threshold:g}",
                self.ensemble,
                f"k{self.ensemble_k}",
            ]
        )

    @classmethod
    def from_dict(cls, data: dict) -> "CandidateConfig":
        """Rebuild a CandidateConfig from a to_dict() payload."""
        return cls(
            name=data.get("name", "restored"),
            de_method=data.get("de_method", "deqms"),
            fdr_method=data.get("fdr_method", "bh"),
            normalization=data.get("normalization", "none"),
            imputation=data.get("imputation", "none"),
            log2fc_threshold=float(data.get("log2fc_threshold", 0.5)),
            ensemble=data.get("ensemble", "none"),
            ensemble_k=int(data.get("ensemble_k", 2)),
            reasoning=data.get("reasoning", ""),
            expected_outcome=data.get("expected_outcome", ""),
        )


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

    def to_dict(self) -> dict:
        """Serialize to JSON-friendly dict."""
        return {
            "config_name": self.config_name,
            "config": self.config,
            "tp": self.tp,
            "fp": self.fp,
            "fn": self.fn,
            "auc": self.auc,
            "sensitivity": self.sensitivity,
            "specificity": self.specificity,
            "n_de_up": self.n_de_up,
            "n_de_down": self.n_de_down,
            "median_cv": self.median_cv,
            "missing_rate": self.missing_rate,
            "score": self.score,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "EvaluationResult":
        """Rebuild an EvaluationResult from a to_dict() payload."""
        return cls(
            config_name=data.get("config_name", ""),
            config=data.get("config", {}) or {},
            tp=data.get("tp"),
            fp=data.get("fp"),
            fn=data.get("fn"),
            auc=data.get("auc"),
            sensitivity=data.get("sensitivity"),
            specificity=data.get("specificity"),
            n_de_up=int(data.get("n_de_up", 0)),
            n_de_down=int(data.get("n_de_down", 0)),
            median_cv=data.get("median_cv"),
            missing_rate=float(data.get("missing_rate", 0.0)),
            score=float(data.get("score", 0.0)),
        )


@dataclass
class ReflectionResult:
    """Output from the reflection step."""

    converged: bool
    next_configs: list[CandidateConfig] = field(default_factory=list)
    analysis: str = ""
    adjustments: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Serialize to JSON-friendly dict."""
        return {
            "converged": self.converged,
            "next_configs": [c.to_dict() for c in self.next_configs],
            "analysis": self.analysis,
            "adjustments": list(self.adjustments),
        }

    @classmethod
    def from_dict(cls, data: dict) -> "ReflectionResult":
        """Rebuild a ReflectionResult from a to_dict() payload."""
        return cls(
            converged=bool(data.get("converged", False)),
            next_configs=[
                CandidateConfig.from_dict(c) for c in data.get("next_configs", [])
            ],
            analysis=data.get("analysis", ""),
            adjustments=list(data.get("adjustments", [])),
        )


@dataclass
class RoundResult:
    """Results from a single optimization round."""

    round_num: int
    configs: list[CandidateConfig]
    results: list[EvaluationResult]
    best_config_name: str = ""
    reflection: ReflectionResult | None = None

    def to_dict(self) -> dict:
        """Serialize to JSON-friendly dict."""
        return {
            "round_num": self.round_num,
            "configs": [c.to_dict() for c in self.configs],
            "results": [r.to_dict() for r in self.results],
            "best_config_name": self.best_config_name,
            "reflection": self.reflection.to_dict() if self.reflection else None,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "RoundResult":
        """Rebuild a RoundResult from a to_dict() payload."""
        reflection = data.get("reflection")
        return cls(
            round_num=int(data.get("round_num", 0)),
            configs=[CandidateConfig.from_dict(c) for c in data.get("configs", [])],
            results=[EvaluationResult.from_dict(r) for r in data.get("results", [])],
            best_config_name=data.get("best_config_name", ""),
            reflection=ReflectionResult.from_dict(reflection) if reflection else None,
        )


@dataclass
class AuditEntry:
    """A single entry in the audit trail."""

    step: str
    round_num: int
    data: Any = None

    def to_dict(self) -> dict:
        """Serialize to JSON-friendly dict."""
        return {"step": self.step, "round_num": self.round_num, "data": self.data}

    @classmethod
    def from_dict(cls, data: dict) -> "AuditEntry":
        """Rebuild an AuditEntry from a to_dict() payload."""
        return cls(
            step=data.get("step", ""),
            round_num=int(data.get("round_num", 0)),
            data=data.get("data"),
        )


@dataclass
class AgenticState:
    """Full state of the optimization workflow."""

    rounds: list[RoundResult] = field(default_factory=list)
    total_experiments: int = 0
    converged: bool = False
    best_config: CandidateConfig | None = None
    best_score: float = 0.0
    audit_trail: list[AuditEntry] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Serialize the full optimization state to a JSON-friendly dict."""
        return {
            "rounds": [r.to_dict() for r in self.rounds],
            "total_experiments": self.total_experiments,
            "converged": self.converged,
            "best_config": self.best_config.to_dict() if self.best_config else None,
            "best_score": self.best_score,
            "audit_trail": [a.to_dict() for a in self.audit_trail],
        }

    @classmethod
    def from_dict(cls, data: dict) -> "AgenticState":
        """Rebuild an AgenticState from a to_dict() payload."""
        best = data.get("best_config")
        return cls(
            rounds=[RoundResult.from_dict(r) for r in data.get("rounds", [])],
            total_experiments=int(data.get("total_experiments", 0)),
            converged=bool(data.get("converged", False)),
            best_config=CandidateConfig.from_dict(best) if best else None,
            best_score=float(data.get("best_score", 0.0)),
            audit_trail=[AuditEntry.from_dict(a) for a in data.get("audit_trail", [])],
        )
