"""Bounded text search over validated Mokume knowledge records."""

from __future__ import annotations

from collections.abc import Collection, Iterable
import re

from mokume.agentic.dataset_knowledge import DatasetEvidence
from mokume.agentic.knowledge import EvidenceRecord, KnowledgeGraph


SearchRecord = EvidenceRecord | DatasetEvidence
RankedRecord = tuple[int, int, str, SearchRecord]


def method_key(value: str) -> str:
    """Normalize a method name for exact method filtering."""
    return re.sub(r"[^a-z0-9]+", "", value.casefold())


def matching_knowledge_records(
    graph: KnowledgeGraph,
    query_tokens: set[str],
    data_type: str | None,
    selected_method: str | None,
    limit: int,
) -> list[SearchRecord]:
    """Rank method evidence and dataset summaries for one validated query."""
    requested = _requested_datasets(graph.datasets.values(), query_tokens)
    ranked = _rank_evidence(
        graph, query_tokens, data_type, selected_method, requested.values()
    )
    ranked.extend(
        _rank_datasets(graph, query_tokens, data_type, selected_method, requested)
    )
    ranked.sort(key=lambda item: item[:3])
    return [record for _, _, _, record in ranked[:limit]]


def _requested_datasets(
    records: Iterable[DatasetEvidence], query_tokens: set[str]
) -> dict[str, DatasetEvidence]:
    return {
        record.id: record
        for record in records
        if _dataset_identifiers(record) & query_tokens
    }


def _rank_evidence(
    graph: KnowledgeGraph,
    query_tokens: set[str],
    data_type: str | None,
    selected_method: str | None,
    requested: Collection[DatasetEvidence],
) -> list[RankedRecord]:
    ranked: list[RankedRecord] = []
    for record in graph.evidence.values():
        if requested and not _record_references_datasets(record, requested):
            continue
        if data_type and record.applicability.data_type.upper() != data_type:
            continue
        if selected_method and selected_method not in _record_method_keys(record):
            continue
        search_text = _record_search_text(graph, record)
        score = sum(token in search_text for token in query_tokens)
        if score:
            ranked.append((-score, record.priority, record.id, record))
    return ranked


def _rank_datasets(
    graph: KnowledgeGraph,
    query_tokens: set[str],
    data_type: str | None,
    selected_method: str | None,
    requested: dict[str, DatasetEvidence],
) -> list[RankedRecord]:
    ranked: list[RankedRecord] = []
    for record in graph.datasets.values():
        if requested and record.id not in requested:
            continue
        if data_type and record.data_type != data_type:
            continue
        identifier_match = bool(_dataset_identifiers(record) & query_tokens)
        if selected_method and (
            not identifier_match or selected_method not in _dataset_method_keys(record)
        ):
            continue
        search_text = _dataset_search_text(graph, record)
        score = sum(token in search_text for token in query_tokens)
        if score:
            ranked.append((-score - 10 * identifier_match, 30, record.id, record))
    return ranked


def _dataset_identifiers(record: DatasetEvidence) -> set[str]:
    return {record.dataset_id.casefold(), record.project_accession.casefold()}


def _record_method_keys(record: EvidenceRecord) -> set[str]:
    values = [
        record.pipeline.quantification,
        record.pipeline.normalization,
        record.pipeline.imputation,
        record.pipeline.de_method,
        record.pipeline.fdr_method,
        record.pipeline.ensemble,
        record.applicability.setting,
        record.applicability.upstream_engine,
    ]
    return {method_key(str(value)) for value in values if value is not None}


def _record_search_text(graph: KnowledgeGraph, record: EvidenceRecord) -> str:
    source = graph.sources[record.source_id]
    values = [
        record.id,
        record.source_id,
        record.kind,
        record.status,
        record.confidence,
        *record.applicability.to_dict().values(),
        *record.pipeline.to_dict().values(),
        *record.metrics.keys(),
        *record.metrics.values(),
        *record.limitations,
        source.id,
        source.kind,
        source.title,
        source.trust,
        source.status,
    ]
    if record.reference_profile is not None:
        values.extend(record.reference_profile.projects)
    return " ".join(str(value).casefold() for value in values if value is not None)


def _record_references_datasets(
    record: EvidenceRecord, datasets: Collection[DatasetEvidence]
) -> bool:
    if record.reference_profile is None:
        return False
    referenced = {project.casefold() for project in record.reference_profile.projects}
    return any(_dataset_identifiers(dataset) & referenced for dataset in datasets)


def _dataset_method_keys(record: DatasetEvidence) -> set[str]:
    return {method_key(value) for value in record.method_values()}


def _dataset_search_text(graph: KnowledgeGraph, record: DatasetEvidence) -> str:
    source = graph.sources[record.source_id]
    values = [
        record.id,
        record.dataset_id,
        record.project_accession,
        record.data_type,
        record.design_family,
        record.factor_column,
        record.condition_map,
        record.contrasts,
        record.ground_truth,
        record.benchmark,
        record.held_out_evaluation,
        source.id,
        source.kind,
        source.title,
        source.trust,
        source.status,
    ]
    return " ".join(str(value).casefold() for value in values if value is not None)
