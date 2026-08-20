"""Deterministic candidate generation from policy-selected evidence."""

from __future__ import annotations

from mokume.agentic.context import BoundContext, bind_context, required_limitations
from mokume.agentic.knowledge import EvidenceRecord
from mokume.agentic.profiler import DataProfile
from mokume.agentic.state import CandidateConfig
from mokume.core.logger import get_logger


logger = get_logger("mokume.agentic.rules")
MAX_RULE_CANDIDATES = 5


def rule_propose(
    profile: DataProfile,
    context: BoundContext | None = None,
) -> list[CandidateConfig]:
    """Generate a small, traceable search set from bound evidence."""
    context = context or bind_context(profile)
    if not context.scope.generation_allowed:
        errors = [
            item.message for item in context.diagnostics if item.severity == "error"
        ]
        raise ValueError("Agentic policy rejected generation: " + "; ".join(errors))

    diagnostic_codes = {item.code for item in context.diagnostics}
    configs: list[CandidateConfig] = []
    # One setting preset plus one independent published reference leaves room
    # for controls without letting deterministic fallback consume the whole
    # experiment budget with near-duplicate priors.
    for record in context.evidence[:2]:
        _append_unique(configs, _from_evidence(record, context))

    _append_unique(configs, _baseline_control(context, diagnostic_codes))
    _append_no_imputation_control(configs, diagnostic_codes)
    _append_deqms_control(configs, diagnostic_codes)

    result = configs[:MAX_RULE_CANDIDATES]
    logger.info("Policy generated %d traceable candidate configs", len(result))
    return result


def _baseline_control(
    context: BoundContext,
    diagnostic_codes: set[str],
) -> CandidateConfig:
    """Build the conservative non-imputed limma control."""
    baseline_refs = [item.id for item in context.evidence[:2]]
    baseline_limitations = list(required_limitations(context, baseline_refs))
    baseline_limitations.append(
        "This conservative control is a search baseline, not a benchmark winner."
    )
    return CandidateConfig(
        name="control_limma_bh_none_none",
        de_method="limma",
        fdr_method="bh",
        normalization="none",
        imputation="none",
        log2fc_threshold=(
            "auto" if "TMT_AUTO_GATE_CONTROL" in diagnostic_codes else 0.5
        ),
        reasoning="Conservative control for comparison with evidence-derived configs.",
        expected_outcome="Stable reference point for the bounded search.",
        evidence_refs=baseline_refs,
        confidence="low",
        limitations=baseline_limitations,
        generated_by="policy",
    )


def _append_no_imputation_control(
    configs: list[CandidateConfig],
    diagnostic_codes: set[str],
) -> None:
    """Add a one-axis no-imputation control when missingness is low."""
    if (
        configs
        and "LOW_MISSINGNESS_CONTROL" in diagnostic_codes
        and configs[0].imputation != "none"
    ):
        source = configs[0]
        _append_unique(
            configs,
            CandidateConfig(
                name=f"control_no_imputation_{source.de_method}",
                de_method=source.de_method,
                fdr_method=source.fdr_method,
                normalization=source.normalization,
                imputation="none",
                log2fc_threshold=source.log2fc_threshold,
                reasoning="Low-missingness no-imputation control for the leading prior.",
                expected_outcome="Tests whether imputation adds value on this matrix.",
                evidence_refs=list(source.evidence_refs),
                confidence="low",
                limitations=[
                    *source.limitations,
                    "This one-axis control was not itself selected as the preset.",
                ],
                generated_by="policy",
            ),
        )


def _append_deqms_control(
    configs: list[CandidateConfig],
    diagnostic_codes: set[str],
) -> None:
    """Add a count-independent control when DEqMS lacks peptide counts."""
    if (
        configs
        and "DEQMS_COUNT_FALLBACK" in diagnostic_codes
        and configs[0].de_method == "deqms"
    ):
        source = configs[0]
        _append_unique(
            configs,
            CandidateConfig(
                name=f"control_limma_{source.normalization}_{source.imputation}",
                de_method="limma",
                fdr_method=source.fdr_method,
                normalization=source.normalization,
                imputation=source.imputation,
                log2fc_threshold=source.log2fc_threshold,
                reasoning="Count-independent control because peptide counts are unavailable.",
                expected_outcome="Separates DEqMS count fallback from preprocessing effects.",
                evidence_refs=list(source.evidence_refs),
                confidence="low",
                limitations=[
                    *source.limitations,
                    "This one-axis control was not itself selected as the preset.",
                ],
                generated_by="policy",
            ),
        )


def _from_evidence(
    record: EvidenceRecord,
    context: BoundContext,
) -> CandidateConfig:
    """Project a full pipeline record onto executable matrix-level axes."""
    pipeline = record.pipeline
    limitations = list(required_limitations(context, [record.id]))
    if pipeline.quantification is not None:
        limitations.append(
            f"Evidence used quantification={pipeline.quantification}; this entry point "
            "does not change the matrix's upstream quantification."
        )
    return CandidateConfig(
        name=record.id.replace("-", "_"),
        de_method=pipeline.de_method,
        fdr_method=pipeline.fdr_method,
        normalization=pipeline.normalization,
        imputation=pipeline.imputation,
        log2fc_threshold=pipeline.log2fc_threshold,
        ensemble=pipeline.ensemble,
        ensemble_k=pipeline.ensemble_k,
        reasoning=f"Executable matrix-level slice of evidence {record.id}.",
        expected_outcome="Benchmark-informed candidate; current-data performance is unknown.",
        evidence_refs=[record.id],
        confidence=record.confidence,
        limitations=_deduplicate(limitations),
        generated_by="policy",
    )


def _append_unique(
    configs: list[CandidateConfig],
    candidate: CandidateConfig,
) -> None:
    """Append a candidate only when its executable signature is new."""
    if all(existing.signature() != candidate.signature() for existing in configs):
        configs.append(candidate)


def _deduplicate(items: list[str]) -> list[str]:
    """Deduplicate strings while preserving their first-seen order."""
    return list(dict.fromkeys(items))
