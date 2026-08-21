"""Bind typed profile, policy, evidence, and contract blocks for generation."""

from __future__ import annotations

from dataclasses import dataclass
import json
from typing import Any

from .contract import (
    CONFIDENCE_LEVELS,
    EXECUTABLE_AXES,
    GENERATED_BLOCK_FIELDS,
    GENERATED_CONFIG_FIELDS,
    method_contract,
    validate_config_values,
)
from .knowledge import (
    EvidenceRecord,
    KnowledgeGraph,
    load_knowledge_graph,
)
from .policy import Diagnostic, PolicyDecision, evaluate_policy
from .profiler import DataProfile


@dataclass(frozen=True)
class ContextBlock:
    """A named, read-only input block exposed to the generator."""

    id: str
    type: str
    content: Any

    def to_dict(self) -> dict[str, Any]:
        """Serialize a context block."""
        return {"id": self.id, "type": self.type, "content": self.content}


@dataclass(frozen=True)
class GenerationScope:
    """Read/write and evidence-reference boundary for the host model."""

    readable_blocks: tuple[str, ...]
    writable_block: str
    allowed_evidence_refs: tuple[str, ...]
    generation_allowed: bool

    def to_dict(self) -> dict[str, Any]:
        """Serialize the generation boundary."""
        return {
            "readable_blocks": list(self.readable_blocks),
            "writable_block": self.writable_block,
            "allowed_evidence_refs": list(self.allowed_evidence_refs),
            "generation_allowed": self.generation_allowed,
            "rules": [
                "Do not mutate or invent any readable context block.",
                "Every evidence reference must be in allowed_evidence_refs.",
                "Do not describe a benchmark candidate as optimal for the current dataset.",
                "Quantification is evidence context, not an executable axis at this entry point.",
            ],
        }


@dataclass(frozen=True)
class BoundContext:
    """Parametric blocks and the scope governing a generated block."""

    blocks: tuple[ContextBlock, ...]
    scope: GenerationScope
    policy: PolicyDecision
    graph: KnowledgeGraph

    def to_dict(self) -> dict[str, Any]:
        """Serialize the exact context returned to the host model."""
        return {
            "knowledge_fingerprint": self.graph.fingerprint,
            "blocks": [block.to_dict() for block in self.blocks],
            "generation_scope": self.scope.to_dict(),
        }

    def to_json(self) -> str:
        """Serialize the bound context deterministically for prompts and audit."""
        return json.dumps(self.to_dict(), indent=2, sort_keys=True)

    @property
    def evidence(self) -> tuple[EvidenceRecord, ...]:
        """Return evidence selected by deterministic policy."""
        return self.policy.selected_evidence

    @property
    def diagnostics(self) -> tuple[Diagnostic, ...]:
        """Return deterministic policy findings."""
        return self.policy.diagnostics


@dataclass(frozen=True)
class GeneratedRecommendationBlock:
    """Validated writable block returned by the proposal model."""

    configs: tuple[dict[str, Any], ...]
    evidence_refs: tuple[str, ...]
    confidence: str
    limitations: tuple[str, ...]
    abstain_reason: str | None


def bind_context(
    profile: DataProfile,
    graph: KnowledgeGraph | None = None,
) -> BoundContext:
    """Create the minimal read-only context permitted for recommendation."""
    graph = graph or load_knowledge_graph()
    decision = evaluate_policy(profile, graph)
    evidence = [
        record.to_context_dict(graph.sources[record.source_id])
        for record in decision.selected_evidence
    ]
    selected_ids = {record.id for record in decision.selected_evidence}
    evidence_graph = {
        "nodes": evidence,
        "edges": [
            edge.to_dict() for edge in graph.edges if edge.source in selected_ids
        ],
    }
    blocks = (
        ContextBlock("profile", "ProfileBlock", profile.to_dict()),
        ContextBlock("contract", "ContractBlock", method_contract()),
        ContextBlock(
            "diagnostics",
            "DiagnosticBlock",
            [item.to_dict() for item in decision.diagnostics],
        ),
        ContextBlock("evidence", "EvidenceBlock", evidence_graph),
    )
    scope = GenerationScope(
        readable_blocks=tuple(block.id for block in blocks),
        writable_block="GeneratedRecommendationBlock",
        allowed_evidence_refs=tuple(item.id for item in decision.selected_evidence),
        generation_allowed=decision.generation_allowed,
    )
    return BoundContext(blocks, scope, decision, graph)


def validate_generated_recommendation(
    payload: dict[str, Any],
    context: BoundContext,
) -> GeneratedRecommendationBlock:
    """Validate semantic constraints after tool or JSON-schema parsing."""
    if context.scope.writable_block != "GeneratedRecommendationBlock":
        raise ValueError("Context does not permit GeneratedRecommendationBlock writes")
    return _validate_recommendation_payload(payload, context)


def _validate_recommendation_payload(
    payload: dict[str, Any],
    context: BoundContext,
) -> GeneratedRecommendationBlock:
    """Validate a recommendation payload within its enclosing writable block."""
    if not isinstance(payload, dict):
        raise ValueError("GeneratedRecommendationBlock must be an object")
    _validate_block_fields(payload)
    configs = payload["configs"]
    evidence_refs = payload["evidence_refs"]
    confidence = payload["confidence"]
    limitations = payload["limitations"]
    abstain_reason = payload["abstain_reason"]
    _validate_block_values(
        configs,
        evidence_refs,
        confidence,
        limitations,
        abstain_reason,
    )
    block = GeneratedRecommendationBlock(
        configs=tuple(configs),
        evidence_refs=tuple(evidence_refs),
        confidence=confidence,
        limitations=tuple(limitations),
        abstain_reason=abstain_reason,
    )
    _validate_evidence_scope(block.evidence_refs, context)
    _validate_generation_decision(block, context)
    return block


def _validate_block_fields(payload: dict[str, Any]) -> None:
    """Require the exact writable-block contract."""
    if set(payload) != set(GENERATED_BLOCK_FIELDS):
        raise ValueError(
            "GeneratedRecommendationBlock fields must be exactly "
            f"{sorted(GENERATED_BLOCK_FIELDS)}"
        )


def _validate_block_values(
    configs: Any,
    evidence_refs: Any,
    confidence: Any,
    limitations: Any,
    abstain_reason: Any,
) -> None:
    """Validate top-level value types before semantic checks."""
    if not isinstance(configs, list) or len(configs) > 5:
        raise ValueError("configs must be a list with at most five entries")
    if not isinstance(evidence_refs, list) or not all(
        isinstance(item, str) for item in evidence_refs
    ):
        raise ValueError("evidence_refs must be a list of strings")
    if len(evidence_refs) != len(set(evidence_refs)):
        raise ValueError("evidence_refs must not contain duplicates")
    if confidence not in CONFIDENCE_LEVELS:
        raise ValueError(f"Unsupported confidence level: {confidence!r}")
    if not isinstance(limitations, list) or not all(
        isinstance(item, str) and item.strip() for item in limitations
    ):
        raise ValueError("limitations must be a list of non-empty strings")
    if abstain_reason is not None and (
        not isinstance(abstain_reason, str) or not abstain_reason.strip()
    ):
        raise ValueError("abstain_reason must be a non-empty string or null")


def _validate_evidence_scope(
    evidence_refs: tuple[str, ...],
    context: BoundContext,
) -> None:
    """Reject references outside the bound evidence graph."""
    unknown_refs = set(evidence_refs) - set(context.scope.allowed_evidence_refs)
    if unknown_refs:
        raise ValueError(
            f"Generated recommendation cited unknown evidence: {unknown_refs}"
        )


def _validate_generation_decision(
    block: GeneratedRecommendationBlock,
    context: BoundContext,
) -> None:
    """Enforce the policy-controlled recommendation or abstention branch."""
    if block.abstain_reason is not None:
        if block.configs or block.evidence_refs or block.confidence != "low":
            raise ValueError(
                "Abstention requires empty configs and evidence_refs with low confidence"
            )
        return
    if not context.scope.generation_allowed:
        raise ValueError("Policy errors require abstention")
    if block.configs:
        if not block.evidence_refs:
            raise ValueError(
                "Generated configs must cite at least one allowed evidence id"
            )
        if not block.limitations:
            raise ValueError("Generated configs must state at least one limitation")
        missing_limitations = [
            item
            for item in required_limitations(context, block.evidence_refs)
            if item not in block.limitations
        ]
        if missing_limitations:
            raise ValueError(
                "Generated recommendation omitted required limitations: "
                f"{missing_limitations}"
            )
        _validate_generated_configs(list(block.configs))
        _validate_confidence(block.confidence, list(block.evidence_refs), context)
    elif not block.abstain_reason:
        raise ValueError("An empty config list requires abstain_reason")


def _validate_generated_configs(configs: list[Any]) -> None:
    """Validate every executable candidate against the runtime contract."""
    names: set[str] = set()
    signatures: set[tuple[Any, ...]] = set()
    for item in configs:
        if not isinstance(item, dict):
            raise ValueError("Each generated config must be an object")
        if set(item) != set(GENERATED_CONFIG_FIELDS):
            raise ValueError(
                "Generated config fields must be exactly "
                f"{sorted(GENERATED_CONFIG_FIELDS)}"
            )
        for field in ("name", "reasoning", "expected_outcome"):
            value = item[field]
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{field} must be a non-empty string")
        if item["name"] in names:
            raise ValueError("Generated candidate names must be unique")
        names.add(item["name"])
        validate_config_values(item)
        signature = tuple(item[field] for field in EXECUTABLE_AXES)
        if signature in signatures:
            raise ValueError(
                "Generated candidates must have unique executable settings"
            )
        signatures.add(signature)


def required_limitations(
    context: BoundContext,
    evidence_refs: list[str] | tuple[str, ...],
) -> tuple[str, ...]:
    """Return limitations that a recommendation must reproduce verbatim."""
    refs = set(evidence_refs)
    limitations: list[str] = []
    for evidence_id in evidence_refs:
        limitations.extend(context.graph.evidence[evidence_id].limitations)
    for diagnostic in context.diagnostics:
        diagnostic_refs = set(diagnostic.evidence_refs)
        if not diagnostic_refs or diagnostic_refs & refs:
            limitations.append(diagnostic.message)
    return tuple(dict.fromkeys(limitations))


def _validate_confidence(
    confidence: str,
    evidence_refs: list[str],
    context: BoundContext,
) -> None:
    """Prevent generated confidence from exceeding cited evidence."""
    rank = {name: index for index, name in enumerate(CONFIDENCE_LEVELS)}
    cited = [context.graph.evidence[item].confidence for item in evidence_refs]
    ceiling = min(cited, key=rank.__getitem__)
    if any(
        diagnostic.code == "OUTSIDE_BENCHMARK_PROFILE"
        and set(diagnostic.evidence_refs) & set(evidence_refs)
        for diagnostic in context.diagnostics
    ):
        ceiling = "low"
    if rank[confidence] > rank[ceiling]:
        raise ValueError(
            f"Generated confidence {confidence!r} exceeds cited evidence "
            f"confidence ceiling {ceiling!r}"
        )
