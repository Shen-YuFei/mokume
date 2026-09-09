"""Dataset-level summaries extracted from validated benchmark artifacts."""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any, NamedTuple, TYPE_CHECKING

import yaml

if TYPE_CHECKING:
    from mokume.agentic.knowledge import SourceEnvelope


_DATA_TYPES = {"DIA", "LFQ", "TMT"}


class DatasetEvidence(NamedTuple):
    """One searchable dataset summary that is never a recommendation prior."""

    id: str
    source_id: str
    dataset_id: str
    project_accession: str
    data_type: str
    design_family: str
    factor_column: str
    condition_map: dict[str, str]
    contrasts: tuple[dict[str, Any], ...]
    ground_truth: dict[str, Any]
    dataset_url: str
    preset_eligible: bool
    benchmark: dict[str, Any]
    held_out_evaluation: dict[str, Any] | None

    def to_context_dict(self, source: SourceEnvelope) -> dict[str, Any]:
        """Serialize bounded dataset evidence for explanation-only search."""
        limitations = [
            "Dataset summaries are explanation-only and are not recommendation priors."
        ]
        if not self.preset_eligible:
            limitations.append(
                "This dataset is excluded from cross-dataset preset aggregation."
            )
        if self.held_out_evaluation is None:
            limitations.append("No leave-one-dataset-out evaluation is available.")
        return {
            "id": self.id,
            "kind": "dataset_benchmark_summary",
            "status": "complete" if self.benchmark["grid_complete"] else "incomplete",
            "eligible_as_prior": False,
            "preset_eligible": self.preset_eligible,
            "applicability": {"data_type": self.data_type},
            "dataset": {
                "id": self.dataset_id,
                "project_accession": self.project_accession,
                "url": self.dataset_url,
            },
            "design": {
                "family": self.design_family,
                "factor_column": self.factor_column,
                "condition_map": dict(self.condition_map),
                "contrasts": [dict(item) for item in self.contrasts],
            },
            "ground_truth": dict(self.ground_truth),
            "benchmark": dict(self.benchmark),
            "held_out_evaluation": (
                dict(self.held_out_evaluation)
                if self.held_out_evaluation is not None
                else None
            ),
            "limitations": limitations,
            "source": source.to_dict(),
        }

    def method_values(self) -> tuple[str, ...]:
        """Return methods selected in this dataset's held-out evaluation."""
        if self.held_out_evaluation is None:
            return ()
        selected = self.held_out_evaluation["selected_pipeline"]
        return tuple(str(value) for value in selected.values() if value is not None)


def load_dataset_evidence(
    source_id: str, source_root: Path
) -> dict[str, DatasetEvidence]:
    """Join one benchmark manifest with its project status and held-out results."""
    manifest = _yaml_mapping(source_root / "benchmark_manifest.yaml")
    projects = _mapping(manifest.get("projects"), "benchmark projects")
    search_space = _mapping(manifest.get("search_space"), "benchmark search_space")
    contract_id = _text(search_space.get("contract_id"), "grid contract id")
    statuses = _keyed_tsv(source_root / "project_status.tsv", "project")
    held_out = _keyed_tsv(source_root / "leave_one_dataset_out.tsv", "held_out_project")
    if set(projects) != set(statuses):
        raise ValueError("Benchmark manifest and project status datasets differ")
    mismatched_contracts = [
        dataset_id
        for dataset_id, row in statuses.items()
        if row.get("grid_contract_id") != contract_id
    ]
    if mismatched_contracts:
        raise ValueError(f"Grid contract mismatch for {sorted(mismatched_contracts)}")
    unknown_held_out = set(held_out) - set(projects)
    if unknown_held_out:
        raise ValueError(
            f"Unknown held-out benchmark datasets: {sorted(unknown_held_out)}"
        )

    records: dict[str, DatasetEvidence] = {}
    for dataset_id, project_item in projects.items():
        project = _mapping(project_item, f"benchmark project {dataset_id}")
        record = _dataset_record(
            source_id,
            _text(dataset_id, "dataset id"),
            project,
            statuses[dataset_id],
            held_out.get(dataset_id),
        )
        if record.id in records:
            raise ValueError(f"Duplicate dataset evidence id: {record.id}")
        records[record.id] = record
    return records


def _dataset_record(
    source_id: str,
    dataset_id: str,
    project: dict[str, Any],
    status: dict[str, str],
    held_out: dict[str, str] | None,
) -> DatasetEvidence:
    data_type = _text(project.get("setting"), f"{dataset_id} setting").upper()
    if data_type not in _DATA_TYPES or status.get("setting") != data_type:
        raise ValueError(f"Invalid or inconsistent data type for {dataset_id}")
    benchmark = _mapping(project.get("benchmark"), f"{dataset_id} benchmark")
    preset_eligible = _bool(benchmark.get("preset_eligible", True), "preset_eligible")
    return DatasetEvidence(
        id=f"dataset-{dataset_id}",
        source_id=source_id,
        dataset_id=dataset_id,
        project_accession=_text(project.get("project_accession"), "project accession"),
        data_type=data_type,
        design_family=_text(benchmark.get("design_family"), "design family"),
        factor_column=_text(project.get("factor_column"), "factor column"),
        condition_map=_string_mapping(project.get("condition_map", {})),
        contrasts=_contrasts(benchmark.get("contrasts"), dataset_id),
        ground_truth=_ground_truth(benchmark.get("ground_truth"), dataset_id),
        dataset_url=_text(benchmark.get("dataset_url"), "dataset URL"),
        preset_eligible=preset_eligible,
        benchmark=_benchmark_status(status),
        held_out_evaluation=_held_out_evaluation(held_out, data_type),
    )


def _contrasts(value: Any, dataset_id: str) -> tuple[dict[str, Any], ...]:
    if not isinstance(value, list) or not value:
        raise ValueError(f"Benchmark contrasts must be a non-empty list: {dataset_id}")
    contrasts = []
    for item in value:
        contrast = _mapping(item, f"{dataset_id} contrast")
        result = {
            "id": _text(contrast.get("id"), "contrast id"),
            "group1": _text(contrast.get("group1"), "contrast group1"),
            "group2": _text(contrast.get("group2"), "contrast group2"),
            "expected_direction": _text(
                contrast.get("expected_direction"), "expected direction"
            ),
            "expected_log2fc": _optional_number(
                contrast.get("expected_log2fc"), "expected log2fc"
            ),
        }
        if contrast.get("ratio_type") is not None:
            result["ratio_type"] = _text(contrast["ratio_type"], "ratio type")
        contrasts.append(result)
    return tuple(contrasts)


def _ground_truth(value: Any, dataset_id: str) -> dict[str, Any]:
    item = _mapping(value, f"{dataset_id} ground truth")
    source = _mapping(item.get("source"), f"{dataset_id} ground truth source")
    return {
        "entity_count": _integer(item.get("entity_count"), "ground truth count"),
        "source": {str(key): scalar for key, scalar in source.items()},
    }


def _benchmark_status(row: dict[str, str]) -> dict[str, Any]:
    return {
        "grid_contract_id": _text(row.get("grid_contract_id"), "grid contract id"),
        "declared_candidates": _integer(
            row.get("declared_candidates"), "declared candidates"
        ),
        "successful_candidates": _integer(
            row.get("successful_candidates"), "successful candidates"
        ),
        "failed_candidates": _integer(
            row.get("failed_candidates"), "failed candidates"
        ),
        "contrast_count": _integer(row.get("contrast_count"), "contrast count"),
        "grid_complete": _bool(row.get("grid_complete"), "grid complete"),
        "grid_contract_fingerprint": _text(
            row.get("grid_contract_fingerprint"), "grid contract fingerprint"
        ),
        "project_input_fingerprint": _text(
            row.get("project_input_fingerprint"), "project input fingerprint"
        ),
    }


def _held_out_evaluation(
    row: dict[str, str] | None, data_type: str
) -> dict[str, Any] | None:
    if row is None:
        return None
    if row.get("setting") != data_type:
        raise ValueError("Held-out evaluation data type does not match its dataset")
    return {
        "status": _text(row.get("status"), "held-out status"),
        "training_projects": [
            value for value in row.get("training_projects", "").split(",") if value
        ],
        "ranked_workflows": _integer(row.get("ranked_workflows"), "ranked workflows"),
        "held_out_rank_spearman": _number(
            row.get("held_out_rank_spearman"), "held-out rank Spearman"
        ),
        "selected_pipeline": {
            "id": _text(row.get("selected_config_id"), "selected config id"),
            "quantification": _text(row.get("quant"), "selected quantification"),
            "normalization": _text(row.get("normalization"), "selected normalization"),
            "imputation": _text(row.get("imputation"), "selected imputation"),
            "de_method": _text(row.get("de_method"), "selected DE method"),
            "fdr_method": _text(row.get("fdr_method"), "selected FDR method"),
            "log2fc_threshold": _number(row.get("log2fc"), "selected log2fc"),
        },
        "held_out_benchmark_rank": _number(
            row.get("held_out_benchmark_rank"), "held-out benchmark rank"
        ),
        "held_out_score_a": _number(row.get("held_out_score_A"), "held-out Score A"),
        "oracle_config_id": _text(row.get("oracle_config_id"), "oracle config id"),
        "oracle_score_a": _number(row.get("oracle_score_A"), "oracle Score A"),
        "score_a_regret": _number(row.get("score_A_regret"), "Score A regret"),
    }


def _keyed_tsv(path: Path, key: str) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if key not in (reader.fieldnames or []):
            raise ValueError(f"Benchmark table lacks key column {key!r}: {path.name}")
        rows: dict[str, dict[str, str]] = {}
        for row in reader:
            identity = _text(row.get(key), f"{path.name} {key}")
            if identity in rows:
                raise ValueError(f"Duplicate benchmark table key: {identity}")
            rows[identity] = row
    return rows


def _yaml_mapping(path: Path) -> dict[str, Any]:
    return _mapping(yaml.safe_load(path.read_text(encoding="utf-8")), path.name)


def _mapping(value: Any, name: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be an object")
    return value


def _string_mapping(value: Any) -> dict[str, str]:
    item = _mapping(value, "condition map")
    if not all(
        isinstance(key, str) and isinstance(entry, str) for key, entry in item.items()
    ):
        raise ValueError("condition map must contain text keys and values")
    return dict(item)


def _text(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be non-empty text")
    return value.strip()


def _integer(value: Any, name: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a non-negative integer")
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a non-negative integer") from exc
    if result < 0 or str(result) != str(value):
        raise ValueError(f"{name} must be a non-negative integer")
    return result


def _number(value: Any, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite number")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite number") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite number")
    return result


def _optional_number(value: Any, name: str) -> float | None:
    return None if value is None else _number(value, name)


def _bool(value: Any, name: str) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str) and value.casefold() in {"true", "false"}:
        return value.casefold() == "true"
    raise ValueError(f"{name} must be boolean")
