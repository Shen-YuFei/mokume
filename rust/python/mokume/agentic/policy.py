"""Deterministic policy checks over the agentic knowledge graph."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

from mokume.agentic.contract import is_supported_quantification
from mokume.agentic.knowledge import (
    EvidenceRecord,
    KnowledgeGraph,
    load_knowledge_graph,
)
from mokume.agentic.profiler import DataProfile


Severity = Literal["info", "warning", "error"]


@dataclass(frozen=True)
class Diagnostic:
    """A machine-readable policy finding."""

    code: str
    severity: Severity
    message: str
    evidence_refs: tuple[str, ...] = ()
    details: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Serialize a diagnostic for model context and audit output."""
        return {
            "code": self.code,
            "severity": self.severity,
            "message": self.message,
            "evidence_refs": list(self.evidence_refs),
            "details": dict(self.details),
        }


@dataclass(frozen=True)
class PolicyDecision:
    """Applicable evidence plus diagnostics produced without an LLM."""

    selected_evidence: tuple[EvidenceRecord, ...]
    diagnostics: tuple[Diagnostic, ...]

    @property
    def generation_allowed(self) -> bool:
        """Return whether policy permits a generated recommendation."""
        return not any(item.severity == "error" for item in self.diagnostics)


def evaluate_policy(
    profile: DataProfile,
    graph: KnowledgeGraph | None = None,
    limit: int = 4,
) -> PolicyDecision:
    """Apply compatibility and evidence-quality rules to a data profile."""
    graph = graph or load_knowledge_graph()
    diagnostics = _profile_diagnostics(profile)
    selected, compatibility_diagnostics = _select_compatible_evidence(
        profile,
        graph,
        limit,
    )
    diagnostics.extend(compatibility_diagnostics)

    if not selected:
        diagnostics.append(
            Diagnostic(
                code="NO_APPLICABLE_EVIDENCE",
                severity="error",
                message=(
                    f"No eligible benchmark evidence matches data type "
                    f"{profile.data_type!r} and the declared upstream configuration."
                ),
            )
        )
    else:
        diagnostics.extend(_evidence_diagnostics(profile, selected))
    return PolicyDecision(tuple(selected), tuple(diagnostics))


def _profile_diagnostics(profile: DataProfile) -> list[Diagnostic]:
    """Check profile-level preconditions before selecting evidence."""
    diagnostics: list[Diagnostic] = []
    if profile.n_conditions < 2:
        diagnostics.append(
            Diagnostic(
                code="INSUFFICIENT_CONDITIONS",
                severity="error",
                message="Differential expression requires at least two conditions.",
            )
        )
    if min(profile.samples_per_condition.values(), default=0) < 2:
        diagnostics.append(
            Diagnostic(
                code="INSUFFICIENT_REPLICATES",
                severity="error",
                message="Differential expression requires at least two samples per group.",
            )
        )
    if profile.data_type_source == "inferred":
        diagnostics.append(
            Diagnostic(
                code="PROFILE_DATA_TYPE_INFERRED",
                severity="warning",
                message=(
                    f"Data type {profile.data_type!r} was inferred from sample names; "
                    "declare it explicitly when using setting-specific evidence."
                ),
            )
        )
    if profile.quantification is not None and not is_supported_quantification(
        profile.quantification
    ):
        diagnostics.append(
            Diagnostic(
                code="UNKNOWN_QUANTIFICATION_METADATA",
                severity="warning",
                message=(
                    f"Quantification {profile.quantification!r} is not a canonical "
                    "Mokume method name; engine-independent evidence may still apply."
                ),
            )
        )
    return diagnostics


def _select_compatible_evidence(
    profile: DataProfile,
    graph: KnowledgeGraph,
    limit: int,
) -> tuple[list[EvidenceRecord], list[Diagnostic]]:
    """Select compatible records and explain every excluded near-match."""
    selected: list[EvidenceRecord] = []
    diagnostics: list[Diagnostic] = []
    for record in graph.matching(profile):
        mismatch = _compatibility_diagnostic(profile, record)
        if mismatch is not None:
            diagnostics.append(mismatch)
            continue
        selected.append(record)
        if len(selected) == limit:
            break
    return selected, diagnostics


def _compatibility_diagnostic(
    profile: DataProfile,
    record: EvidenceRecord,
) -> Diagnostic | None:
    """Return the reason a matching-type record cannot be used, if any."""
    evidence_engine = record.applicability.upstream_engine
    if (
        profile.upstream_engine is not None
        and evidence_engine is not None
        and profile.upstream_engine.casefold() != evidence_engine.casefold()
    ):
        return Diagnostic(
            code="UPSTREAM_ENGINE_MISMATCH",
            severity="info",
            message=(
                f"Evidence {record.id} used {evidence_engine}, but the current "
                f"profile declares {profile.upstream_engine}; it was excluded."
            ),
            evidence_refs=(record.id,),
        )
    evidence_quant = record.pipeline.quantification
    if (
        profile.quantification is not None
        and evidence_quant is not None
        and profile.quantification != evidence_quant
    ):
        return Diagnostic(
            code="UPSTREAM_QUANTIFICATION_MISMATCH",
            severity="warning",
            message=(
                f"Evidence {record.id} used {evidence_quant}, but the current "
                f"matrix declares {profile.quantification}; it was excluded."
            ),
            evidence_refs=(record.id,),
        )
    return None


def _evidence_diagnostics(
    profile: DataProfile,
    selected: list[EvidenceRecord],
) -> list[Diagnostic]:
    """Report quality limits and controls for the selected evidence."""
    selected_refs = tuple(item.id for item in selected)
    return [
        *_source_quality_diagnostics(profile, selected, selected_refs),
        *_benchmark_profile_diagnostics(profile, selected),
        *_profile_control_diagnostics(profile, selected_refs),
    ]


def _benchmark_profile_diagnostics(
    profile: DataProfile,
    selected: list[EvidenceRecord],
) -> list[Diagnostic]:
    """Flag inputs outside the measured coverage of a benchmark preset."""
    diagnostics: list[Diagnostic] = []
    observed = {
        "n_proteins": profile.n_proteins,
        "n_samples": profile.n_samples,
        "missing_rate": profile.missing_rate,
    }
    for record in selected:
        reference = record.reference_profile
        if reference is None:
            continue
        ranges = {
            "n_proteins": reference.n_proteins,
            "n_samples": reference.n_samples,
            "missing_rate": reference.missing_rate,
        }
        outside = {
            field: {
                "observed": value,
                "reference_min": ranges[field][0],
                "reference_max": ranges[field][1],
            }
            for field, value in observed.items()
            if not ranges[field][0] <= value <= ranges[field][1]
        }
        if outside:
            diagnostics.append(
                Diagnostic(
                    code="OUTSIDE_BENCHMARK_PROFILE",
                    severity="warning",
                    message=(
                        f"The current dataset is outside the measured profile of "
                        f"evidence {record.id} for: {', '.join(sorted(outside))}. "
                        "Treat this preset as low-confidence exploratory evidence."
                    ),
                    evidence_refs=(record.id,),
                    details=outside,
                )
            )
    return diagnostics


def _source_quality_diagnostics(
    profile: DataProfile,
    selected: list[EvidenceRecord],
    selected_refs: tuple[str, ...],
) -> list[Diagnostic]:
    """Report transfer and method limitations in selected sources."""
    diagnostics: list[Diagnostic] = []
    if profile.quantification is None and any(
        item.pipeline.quantification is not None for item in selected
    ):
        diagnostics.append(
            Diagnostic(
                code="UPSTREAM_QUANTIFICATION_UNKNOWN",
                severity="warning",
                message=(
                    "The current matrix does not declare its quantification method. "
                    "Full-pipeline presets are conditional and only their executable "
                    "protein-matrix slice may be tested."
                ),
                evidence_refs=selected_refs,
            )
        )
    diagnostics.extend(_transfer_diagnostics(selected))
    deqms = tuple(item.id for item in selected if item.pipeline.de_method == "deqms")
    if deqms and not profile.has_peptide_counts:
        diagnostics.append(
            Diagnostic(
                code="DEQMS_COUNT_FALLBACK",
                severity="warning",
                message=(
                    "DEqMS evidence is available, but peptide counts are absent; Mokume "
                    "will fall back to count=1 and lose count-dependent moderation."
                ),
                evidence_refs=deqms,
            )
        )
    return diagnostics


def _transfer_diagnostics(selected: list[EvidenceRecord]) -> list[Diagnostic]:
    """Report provisional presets and weak held-out transfer."""
    diagnostics: list[Diagnostic] = []
    provisional = tuple(
        item.id for item in selected if item.status == "candidate_for_review"
    )
    if provisional:
        diagnostics.append(
            Diagnostic(
                code="PROVISIONAL_PRESET",
                severity="warning",
                message="Selected Grid presets are review candidates, not published defaults.",
                evidence_refs=provisional,
            )
        )
    unstable = tuple(
        item.id
        for item in selected
        if float(item.metrics.get("mean_lodo_rank_spearman", 1.0)) < 0.2
    )
    if unstable:
        diagnostics.append(
            Diagnostic(
                code="LOW_TRANSFER_STABILITY",
                severity="warning",
                message="Held-out rank transfer is weak for part of the selected evidence.",
                evidence_refs=unstable,
            )
        )
    return diagnostics


def _profile_control_diagnostics(
    profile: DataProfile,
    selected_refs: tuple[str, ...],
) -> list[Diagnostic]:
    """Request controls justified by current-matrix characteristics."""
    diagnostics: list[Diagnostic] = []
    if profile.missing_rate < 0.15:
        diagnostics.append(
            Diagnostic(
                code="LOW_MISSINGNESS_CONTROL",
                severity="info",
                message=(
                    "Missingness is below 15%; include a no-imputation control "
                    "rather than assuming imputation improves the matrix."
                ),
                evidence_refs=selected_refs,
            )
        )
    if profile.data_type == "TMT":
        diagnostics.append(
            Diagnostic(
                code="TMT_AUTO_GATE_CONTROL",
                severity="info",
                message=(
                    "Include an auto effect-size-gate control for possible isobaric "
                    "ratio compression; do not replace the evidence's fixed gate."
                ),
                evidence_refs=selected_refs,
            )
        )
    return diagnostics
