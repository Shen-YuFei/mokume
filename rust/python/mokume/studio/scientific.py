"""Server-side scientific planning boundary for Mokume Studio."""

from __future__ import annotations

from pathlib import Path
import re
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator

from mokume.agentic.knowledge import EvidenceRecord, load_knowledge_graph
from mokume.studio.jobs import JobManager
from mokume.studio.models import (
    JobOperation,
    JobSpec,
    ProjectRecord,
    RunRecord,
    ScientificJobRequest,
)
from mokume.studio.paths import ProjectPaths
from mokume.studio.science import (
    ApprovalRecord,
    DatasetInspectionRequest,
    DatasetRecord,
    DatasetStatus,
    ScienceStore,
)


_KNOWLEDGE_DATA_TYPES = {"DIA", "LFQ", "TMT"}
_SEARCH_TOKEN = re.compile(r"[^\W_]+(?:[.+-][^\W_]+)*", re.UNICODE)


def workspace_identity(project: ProjectRecord) -> dict[str, str]:
    """Return the stable absolute identity shown for one active workspace."""
    return {
        "id": project.id,
        "name": Path(project.root).name or project.root,
        "root": project.root,
    }


class EvaluationPlanRequest(BaseModel):
    """Model-selected subset of deterministic policy candidates."""

    model_config = ConfigDict(extra="forbid")

    dataset_id: str = Field(min_length=1)
    config_names: list[str] = Field(min_length=1, max_length=5)
    output_directory: str = Field(min_length=1)
    ground_truth: str | None = None
    expected_direction: Literal["UP", "DOWN"] | None = None

    @field_validator("config_names")
    @classmethod
    def unique_config_names(cls, value: list[str]) -> list[str]:
        """Reject blank or repeated candidate identities."""
        if any(not item.strip() for item in value):
            raise ValueError("config names must not be blank")
        if len(set(value)) != len(value):
            raise ValueError("config names must be unique")
        return value


class ScientificController:
    """Bind Studio state, guarded paths, worker jobs, and Mokume policy."""

    def __init__(self, store: ScienceStore, jobs: JobManager) -> None:
        self.store = store
        self.jobs = jobs
        self.knowledge = load_knowledge_graph()
        self.knowledge_fingerprint = self.knowledge.fingerprint

    def inspect(
        self,
        request: DatasetInspectionRequest,
        project: ProjectRecord,
    ) -> tuple[DatasetRecord, RunRecord]:
        """Queue a read-only dataset inspection in the isolated worker."""
        canonical, input_paths = self._canonical_inspection(request, project)
        dataset = self.store.create_dataset(project.id, canonical)
        payload = {
            "dataset_id": dataset.id,
            "knowledge_fingerprint": self.knowledge_fingerprint,
            "request": self._inspection_service_request(canonical),
        }
        try:
            spec = self.jobs.prepare_scientific_spec(
                ScientificJobRequest(
                    operation=JobOperation.INSPECT_DATASET,
                    payload=payload,
                    input_paths=input_paths,
                    output_directory=canonical.output_directory,
                ),
                project,
            )
            run = self.jobs.submit_spec(spec, project, "inspect_dataset")
        except (OSError, RuntimeError, TypeError, ValueError) as exc:
            self.store.update_dataset(dataset.id, DatasetStatus.FAILED, error=str(exc))
            raise
        return dataset, run

    def context(
        self,
        project: ProjectRecord,
        dataset_id: str | None = None,
    ) -> dict[str, Any]:
        """Return model-safe structured context without matrix rows or sample IDs."""
        dataset = (
            self.store.get_dataset(dataset_id)
            if dataset_id
            else self.store.latest_dataset(project.id)
        )
        if dataset is not None and dataset.project_id != project.id:
            raise ValueError("dataset belongs to a different project")
        payload: dict[str, Any] = {
            "project_id": project.id,
            "workspace": {
                **workspace_identity(project),
                "access": "workspace_only",
            },
            "mokume_threads": 24,
            "knowledge_fingerprint": self.knowledge_fingerprint,
            "capabilities": [
                "inspect_dataset",
                "search_knowledge",
                "evaluate_recommendation",
            ],
            "disclosure": {"metadata": True, "raw_rows": False},
            "dataset": None,
        }
        if dataset is not None:
            payload["dataset"] = self._dataset_context(dataset, project)
        return payload

    def search_knowledge(
        self,
        query: str,
        *,
        data_type: str | None = None,
        method: str | None = None,
        limit: int = 5,
    ) -> dict[str, Any]:
        """Search the validated bundled evidence without reading arbitrary files."""
        (
            normalized_query,
            normalized_type,
            normalized_method,
            query_tokens,
            method_key,
        ) = self._knowledge_search_filters(query, data_type, method, limit)
        results = []
        for record in self._matching_knowledge_records(
            query_tokens, normalized_type, method_key, limit
        ):
            result = record.to_context_dict(self.knowledge.sources[record.source_id])
            result["eligible_as_prior"] = record.eligible_as_prior
            results.append(result)
        return {
            "scope": "explanation_only",
            "execution_authority": False,
            "knowledge_fingerprint": self.knowledge_fingerprint,
            "query": normalized_query,
            "filters": {
                "data_type": normalized_type,
                "method": normalized_method,
            },
            "count": len(results),
            "results": results,
        }

    def _knowledge_search_filters(
        self,
        query: str,
        data_type: str | None,
        method: str | None,
        limit: int,
    ) -> tuple[str, str | None, str | None, set[str], str | None]:
        normalized_query = self._search_value(query, "query", 200)
        normalized_type = (
            self._search_value(data_type, "data_type", 3).upper()
            if data_type is not None
            else None
        )
        if normalized_type is not None and normalized_type not in _KNOWLEDGE_DATA_TYPES:
            raise ValueError("data_type must be DIA, LFQ, or TMT")
        normalized_method = (
            self._search_value(method, "method", 64) if method is not None else None
        )
        if isinstance(limit, bool) or not isinstance(limit, int) or not 1 <= limit <= 5:
            raise ValueError("limit must be an integer from 1 to 5")
        query_tokens = {
            token.casefold() for token in _SEARCH_TOKEN.findall(normalized_query)
        }
        if not query_tokens:
            raise ValueError("query must contain searchable text")
        method_key = self._method_key(normalized_method) if normalized_method else None
        if normalized_method and not method_key:
            raise ValueError("method must contain searchable text")
        return (
            normalized_query,
            normalized_type,
            normalized_method,
            query_tokens,
            method_key,
        )

    def _matching_knowledge_records(
        self,
        query_tokens: set[str],
        data_type: str | None,
        method_key: str | None,
        limit: int,
    ) -> list[EvidenceRecord]:
        ranked: list[tuple[int, int, str, EvidenceRecord]] = []
        for record in self.knowledge.evidence.values():
            if data_type and record.applicability.data_type.upper() != data_type:
                continue
            if method_key and method_key not in self._record_method_keys(record):
                continue
            search_text = self._record_search_text(record)
            score = sum(token in search_text for token in query_tokens)
            if score:
                ranked.append((-score, record.priority, record.id, record))
        ranked.sort(key=lambda item: item[:3])
        return [record for _, _, _, record in ranked[:limit]]

    @staticmethod
    def _search_value(value: str, name: str, maximum: int) -> str:
        if not isinstance(value, str):
            raise ValueError(f"{name} must be text")
        normalized = value.strip()
        if not normalized:
            raise ValueError(f"{name} must not be blank")
        if "\x00" in normalized or len(normalized) > maximum:
            raise ValueError(f"{name} is invalid or too long")
        return normalized

    @staticmethod
    def _method_key(value: str) -> str:
        return re.sub(r"[^a-z0-9]+", "", value.casefold())

    def _record_method_keys(self, record: EvidenceRecord) -> set[str]:
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
        return {self._method_key(str(value)) for value in values if value is not None}

    def _record_search_text(self, record: EvidenceRecord) -> str:
        source = self.knowledge.sources[record.source_id]
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

    def prepare_evaluation(
        self,
        request: EvaluationPlanRequest,
        project: ProjectRecord,
        *,
        provider: str,
        model: str,
    ) -> ApprovalRecord:
        """Create a durable final-parameter approval from policy candidates only."""
        dataset = self._ready_dataset(request.dataset_id, project)
        recommendation = self._selected_recommendation(dataset, request.config_names)
        service_request, inputs = self._evaluation_service_request(
            dataset,
            request,
            project,
            recommendation,
        )
        payload = {
            "dataset_id": dataset.id,
            "knowledge_fingerprint": self.knowledge_fingerprint,
            "plan_source": {"provider": provider, "model": model},
            "request": service_request,
        }
        spec = self.jobs.prepare_scientific_spec(
            ScientificJobRequest(
                operation=JobOperation.EVALUATE_RECOMMENDATION,
                payload=payload,
                input_paths=inputs,
                output_directory=dataset.request.output_directory,
            ),
            project,
        )
        approval_payload = {
            "project_id": project.id,
            "dataset_id": dataset.id,
            "job_spec": spec.model_dump(mode="json"),
            "card": self._approval_card(spec, service_request, request.config_names),
        }
        return self.store.create_approval("evaluate_recommendation", approval_payload)

    def start_approved(
        self,
        approval_id: str,
        payload_hash: str,
        project: ProjectRecord,
    ) -> RunRecord:
        """Consume one unchanged approval and submit its immutable worker spec."""
        approval = self.store.get_approval(approval_id)
        if approval is None:
            raise ValueError("approval not found")
        if approval.payload_hash != payload_hash:
            raise ValueError("approval payload hash mismatch")
        if approval.kind != "evaluate_recommendation":
            raise ValueError("approval does not authorize an evaluation")
        if approval.payload.get("project_id") != project.id:
            raise ValueError("approval belongs to a different project")
        spec = JobSpec.model_validate(approval.payload.get("job_spec"))
        if spec.approved_hash != approval.payload["card"].get("approved_hash"):
            raise ValueError("approval card no longer matches the job specification")

        def consume() -> None:
            self.store.consume_approval(approval_id, payload_hash=payload_hash)

        run = self.jobs.submit_spec(
            spec,
            project,
            "evaluate_recommendation",
            authorize=consume,
        )
        self.store.link_run(
            approval_id,
            payload_hash=payload_hash,
            run_id=run.id,
        )
        return run

    @staticmethod
    def _canonical_inspection(
        request: DatasetInspectionRequest,
        project: ProjectRecord,
    ) -> tuple[DatasetInspectionRequest, list[str]]:
        guard = ProjectPaths(project.root)
        matrix = str(guard.resolve_existing(request.protein_matrix))
        sdrf = str(guard.resolve_existing(request.sdrf))
        peptide_counts = (
            str(guard.resolve_existing(request.peptide_counts))
            if request.peptide_counts
            else None
        )
        guard.resolve_output(request.output_directory, allow_existing=True)
        canonical = request.model_copy(
            update={
                "protein_matrix": matrix,
                "sdrf": sdrf,
                "peptide_counts": peptide_counts,
            }
        )
        inputs = [matrix, sdrf]
        if peptide_counts:
            inputs.append(peptide_counts)
        return canonical, inputs

    @staticmethod
    def _inspection_service_request(request: DatasetInspectionRequest) -> dict:
        metadata = {
            "data_type": request.data_type,
            "quantification": request.quantification,
            "upstream_engine": request.upstream_engine,
            "factor_column": request.factor_column,
        }
        return {
            "protein_matrix": request.protein_matrix,
            "sdrf": request.sdrf,
            "input_scale": request.input_scale,
            "contrast": list(request.contrast),
            "peptide_counts": request.peptide_counts,
            "metadata": metadata,
        }

    def _ready_dataset(
        self,
        dataset_id: str,
        project: ProjectRecord,
    ) -> DatasetRecord:
        dataset = self.store.get_dataset(dataset_id)
        if dataset is None:
            raise ValueError("dataset not found")
        if dataset.project_id != project.id:
            raise ValueError("dataset belongs to a different project")
        if dataset.status is not DatasetStatus.READY or dataset.result is None:
            raise ValueError("dataset inspection is not ready")
        if dataset.result.get("context", {}).get("knowledge_fingerprint") != (
            self.knowledge_fingerprint
        ):
            raise ValueError("dataset inspection used a different knowledge snapshot")
        return dataset

    @staticmethod
    def _selected_recommendation(
        dataset: DatasetRecord,
        config_names: list[str],
    ) -> dict:
        policy = dataset.result["policy_recommendation"]
        if policy.get("abstain_reason"):
            raise ValueError(f"Mokume policy abstained: {policy['abstain_reason']}")
        by_name = {item["name"]: item for item in policy.get("configs", [])}
        unknown = [name for name in config_names if name not in by_name]
        if unknown:
            raise ValueError(f"configs are outside Mokume policy: {unknown}")
        recommendation = dict(policy)
        recommendation["configs"] = [by_name[name] for name in config_names]
        return recommendation

    @staticmethod
    def _evaluation_service_request(
        dataset: DatasetRecord,
        plan: EvaluationPlanRequest,
        project: ProjectRecord,
        recommendation: dict,
    ) -> tuple[dict, list[str]]:
        guard = ProjectPaths(project.root)
        output = str(guard.resolve_output(plan.output_directory))
        ground_truth = (
            str(guard.resolve_existing(plan.ground_truth))
            if plan.ground_truth
            else None
        )
        if bool(ground_truth) != bool(plan.expected_direction):
            raise ValueError(
                "ground_truth and expected_direction must be supplied together"
            )
        request = dataset.request
        options = {
            "input_scale": request.input_scale,
            "peptide_counts": request.peptide_counts,
            "ground_truth": ground_truth,
            "expected_direction": plan.expected_direction,
            "fdr_threshold": 0.05,
            "threads": 24,
            "data_type": request.data_type,
            "quantification": request.quantification,
            "upstream_engine": request.upstream_engine,
            "factor_column": request.factor_column,
        }
        service_request = {
            "protein_matrix": request.protein_matrix,
            "sdrf": request.sdrf,
            "contrast": list(request.contrast),
            "recommendation": recommendation,
            "output_dir": output,
            "options": options,
        }
        inputs = [request.protein_matrix, request.sdrf]
        if request.peptide_counts:
            inputs.append(request.peptide_counts)
        if ground_truth:
            inputs.append(ground_truth)
        return service_request, inputs

    def _dataset_context(
        self,
        dataset: DatasetRecord,
        project: ProjectRecord,
    ) -> dict[str, Any]:
        guard = ProjectPaths(project.root)
        request = dataset.request
        result = dataset.result or {}
        return {
            "id": dataset.id,
            "status": dataset.status.value,
            "files": {
                "protein_matrix": guard.relative(request.protein_matrix),
                "sdrf": guard.relative(request.sdrf),
                "peptide_counts": (
                    guard.relative(request.peptide_counts)
                    if request.peptide_counts
                    else None
                ),
            },
            "contrast": list(request.contrast),
            "input_scale": request.input_scale,
            "profile": result.get("profile"),
            "context": result.get("context"),
            "policy_recommendation": result.get("policy_recommendation"),
            "ranking_contract": result.get("ranking_contract"),
        }

    @staticmethod
    def _approval_card(
        spec: JobSpec,
        service_request: dict,
        config_names: list[str],
    ) -> dict[str, Any]:
        options = service_request["options"]
        return {
            "tool": JobOperation.EVALUATE_RECOMMENDATION.value,
            "contract_version": 1,
            "approved_hash": spec.approved_hash,
            "inputs": spec.parameters["input_snapshots"],
            "contrast": service_request["contrast"],
            "configs": config_names,
            "output_directory": service_request["output_dir"],
            "threads": spec.threads,
            "ground_truth": options["ground_truth"],
            "expected_direction": options["expected_direction"],
            "overwrite": False,
            "model_disclosure": "profile metadata only; raw rows excluded",
        }
