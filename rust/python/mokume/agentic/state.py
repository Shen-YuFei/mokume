"""Typed candidate and evaluation records for the plugin service."""

from __future__ import annotations

from typing import Any, NamedTuple


class CandidateConfig(NamedTuple):
    """One bounded protein-matrix analysis configuration."""

    name: str
    de_method: str = "deqms"
    fdr_method: str = "bh"
    normalization: str = "none"
    imputation: str = "none"
    log2fc_threshold: float | str = 0.5
    fdr_threshold: float = 0.05
    ensemble: str = "none"
    ensemble_k: int | None = None
    reasoning: str = ""
    expected_outcome: str = ""
    evidence_refs: tuple[str, ...] | list[str] = ()
    confidence: str = "unrated"
    limitations: tuple[str, ...] | list[str] = ()
    generated_by: str = "host"

    def to_dict(self) -> dict[str, Any]:
        """Serialize the candidate for validation and audit output."""
        return {
            "name": self.name,
            "de_method": self.de_method,
            "fdr_method": self.fdr_method,
            "normalization": self.normalization,
            "imputation": self.imputation,
            "log2fc_threshold": self.log2fc_threshold,
            "fdr_threshold": self.fdr_threshold,
            "ensemble": self.ensemble,
            "ensemble_k": self.ensemble_k,
            "reasoning": self.reasoning,
            "expected_outcome": self.expected_outcome,
            "evidence_refs": list(self.evidence_refs),
            "confidence": self.confidence,
            "limitations": list(self.limitations),
            "generated_by": self.generated_by,
        }

    def signature(self) -> str:
        """Return a stable identity for executable method choices."""
        return "|".join(
            [
                self.de_method,
                self.fdr_method,
                self.normalization,
                self.imputation,
                f"fc{self.log2fc_threshold}",
                f"a{self.fdr_threshold:g}",
                self.ensemble,
                f"k{self.ensemble_k if self.ensemble_k is not None else 'none'}",
            ]
        )


class GroundTruthMetrics(NamedTuple):
    """Score A and directional metrics available only with known truth."""

    tp: int | None = None
    fp: int | None = None
    fn: int | None = None
    tn: int | None = None
    auc: float | None = None
    sensitivity: float | None = None
    specificity: float | None = None
    truth_direction_correct: int | None = None
    truth_direction_incorrect: int | None = None
    truth_direction_accuracy: float | None = None
    recall_at_emp_fdr: float | None = None
    pauc: float | None = None
    pauc001: float | None = None
    pauc005: float | None = None
    nmcc: float | None = None
    gmean: float | None = None
    recall_emp_fdr_curve: tuple[tuple[float, float], ...] | None = None

    def to_dict(self) -> dict[str, Any]:
        """Serialize the flattened public Score A fields."""
        return {
            "tp": self.tp,
            "fp": self.fp,
            "fn": self.fn,
            "tn": self.tn,
            "auc": self.auc,
            "sensitivity": self.sensitivity,
            "specificity": self.specificity,
            "truth_direction_correct": self.truth_direction_correct,
            "truth_direction_incorrect": self.truth_direction_incorrect,
            "truth_direction_accuracy": self.truth_direction_accuracy,
            "recall_at_emp_fdr": self.recall_at_emp_fdr,
            "pauc": self.pauc,
            "pauc001": self.pauc001,
            "pauc005": self.pauc005,
            "nmcc": self.nmcc,
            "gmean": self.gmean,
            "recall_emp_fdr_curve": (
                [[alpha, recall] for alpha, recall in self.recall_emp_fdr_curve]
                if self.recall_emp_fdr_curve is not None
                else None
            ),
        }


class FdrCalibrationMetrics(NamedTuple):
    """Side-by-side BH and adaptive-FDR diagnostics."""

    n_called_bh: int | None = None
    emp_fdr_bh: float | None = None
    recall_at_alpha_bh: float | None = None
    n_called_adaptive: int | None = None
    emp_fdr_adaptive: float | None = None
    recall_at_alpha_adaptive: float | None = None
    adaptive_method_used: str | None = None
    adaptive_pi0: float | None = None

    def to_dict(self) -> dict[str, Any]:
        """Serialize the flattened public FDR calibration fields."""
        return {
            "n_called_bh": self.n_called_bh,
            "emp_fdr_bh": self.emp_fdr_bh,
            "recall_at_alpha_bh": self.recall_at_alpha_bh,
            "n_called_adaptive": self.n_called_adaptive,
            "emp_fdr_adaptive": self.emp_fdr_adaptive,
            "recall_at_alpha_adaptive": self.recall_at_alpha_adaptive,
            "adaptive_method_used": self.adaptive_method_used,
            "adaptive_pi0": self.adaptive_pi0,
        }


class EvaluationResult(NamedTuple):
    """Metrics for one candidate, with Score A fields optional without truth."""

    config_name: str
    config: dict[str, Any]
    truth_metrics: GroundTruthMetrics = GroundTruthMetrics()
    fdr_calibration: FdrCalibrationMetrics = FdrCalibrationMetrics()
    n_de_up: int = 0
    n_de_down: int = 0
    median_cv: float | None = None
    missing_rate: float = 0.0
    score_a: float | None = None

    def with_score_a(self, score_a: float | None) -> EvaluationResult:
        """Return a copy carrying the computed Score A value."""
        return EvaluationResult(
            config_name=self.config_name,
            config=self.config,
            truth_metrics=self.truth_metrics,
            fdr_calibration=self.fdr_calibration,
            n_de_up=self.n_de_up,
            n_de_down=self.n_de_down,
            median_cv=self.median_cv,
            missing_rate=self.missing_rate,
            score_a=score_a,
        )

    def to_dict(self) -> dict[str, Any]:
        """Serialize the measurement to a JSON-compatible dictionary."""
        payload = {
            "config_name": self.config_name,
            "config": self.config,
            "n_de_up": self.n_de_up,
            "n_de_down": self.n_de_down,
            "median_cv": self.median_cv,
            "missing_rate": self.missing_rate,
            "score_a": self.score_a,
        }
        payload.update(self.truth_metrics.to_dict())
        payload.update(self.fdr_calibration.to_dict())
        return payload
