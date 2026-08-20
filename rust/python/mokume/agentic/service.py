"""Deterministic service boundary used by the Mokume Plugin MCP server."""

from __future__ import annotations

import csv
from dataclasses import dataclass
import json
from pathlib import Path
import re
from typing import Any, NamedTuple

import numpy as np
import pandas as pd

from mokume.agentic.context import (
    BoundContext,
    GeneratedRecommendationBlock,
    bind_context,
    validate_generated_recommendation,
)
from mokume.agentic.contract import GENERATED_CONFIG_FIELDS
from mokume.agentic.evaluator import (
    compute_score_ground_truth,
    evaluate,
    method_sensitivity,
)
from mokume.agentic.knowledge import KnowledgeGraph, load_knowledge_graph
from mokume.agentic.profiler import AcquisitionMetadata, DataProfile, profile_data
from mokume.agentic.rules import rule_propose
from mokume.agentic.runner import (
    ExperimentContext,
    PreprocessCache,
    run_experiment,
)
from mokume.agentic.state import CandidateConfig, EvaluationResult
from mokume.normalization.irs import detect_condition_from_sdrf

_METADATA_FIELDS = frozenset(
    {"data_type", "quantification", "upstream_engine", "factor_column"}
)
_OPTION_FIELDS = _METADATA_FIELDS | {
    "ground_truth",
    "expected_direction",
    "fdr_threshold",
    "input_scale",
    "threads",
}


@dataclass(frozen=True)
class DatasetMetadata:
    """Declared acquisition facts and the SDRF factor selection."""

    acquisition: AcquisitionMetadata
    factor_column: str | None = None


@dataclass(frozen=True)
class EvaluationOptions:
    """Validated optional settings for one recommendation round."""

    metadata: DatasetMetadata
    ground_truth: str | None = None
    expected_direction: str | None = None
    fdr_threshold: float = 0.05
    input_scale: str = "auto"
    threads: int = 24


@dataclass(frozen=True)
class EvaluationRequest:
    """Required and optional inputs for one evaluation tool call."""

    protein_matrix: str
    sdrf: str
    contrast: list[str]
    recommendation: dict[str, Any]
    output_dir: str
    options: dict[str, Any] | None = None


@dataclass(frozen=True)
class EvaluationRun:
    """Inputs shared by every candidate in one evaluation round."""

    experiment: ExperimentContext
    ground_truth: set[str] | None
    expected_direction: str | None
    output_dir: Path


class PreparedEvaluation(NamedTuple):
    """Validated request state ready for execution or abstention."""

    runtime: EvaluationOptions
    contrast: tuple[str, str]
    destination: Path
    frame: pd.DataFrame
    profile: DataProfile
    conditions: dict[str, str]
    context: BoundContext
    generated: GeneratedRecommendationBlock


class RecommendationService:
    """Bind one immutable knowledge graph to the two public MCP operations."""

    def __init__(self, knowledge: str) -> None:
        knowledge_path = str(Path(knowledge).resolve())
        _require_file(knowledge_path, "knowledge")
        self._graph = load_knowledge_graph(knowledge_path)

    def inspect_dataset(
        self,
        protein_matrix: str,
        sdrf: str,
        metadata: dict[str, str | None] | None = None,
    ) -> dict[str, Any]:
        """Profile a matrix and bind only compatible, traceable evidence."""
        _require_file(protein_matrix, "protein_matrix")
        _require_file(sdrf, "sdrf")
        frame = _read_matrix(protein_matrix)
        profile, _, context = self._profile(
            frame,
            sdrf,
            _parse_dataset_metadata(metadata),
        )
        return {
            "profile": profile.to_dict(),
            "context": context.to_dict(),
            "policy_recommendation": _policy_recommendation(profile, context),
            "ranking_contract": "unranked_without_ground_truth",
        }

    def evaluate_recommendation(
        self,
        request: EvaluationRequest,
    ) -> dict[str, Any]:
        """Validate and evaluate one host-generated recommendation block."""
        prepared = self._prepare_evaluation(request)
        if prepared.generated.abstain_reason is not None:
            return _write_evaluation(
                prepared.destination,
                _abstention_payload(self._graph, prepared),
            )

        linear_frame, resolved_scale = _as_linear(
            prepared.frame,
            prepared.runtime.input_scale,
            prepared.profile.is_log_transformed,
        )
        experiment = ExperimentContext(
            sample_to_condition=prepared.conditions,
            contrast=prepared.contrast,
            threads=prepared.runtime.threads,
        )
        run = EvaluationRun(
            experiment=experiment,
            ground_truth=_load_ground_truth(prepared.runtime.ground_truth),
            expected_direction=prepared.runtime.expected_direction,
            output_dir=prepared.destination,
        )
        candidates = [
            _candidate(item, prepared.generated, prepared.runtime.fdr_threshold)
            for item in prepared.generated.configs
        ]
        payload = _run_candidates(candidates, linear_frame, run)
        payload.update(
            _audit_fields(
                self._graph,
                prepared.generated,
                prepared.context,
                input_scale=resolved_scale,
            )
        )
        return _write_evaluation(prepared.destination, payload)

    def _prepare_evaluation(
        self,
        request: EvaluationRequest,
    ) -> PreparedEvaluation:
        """Validate and bind one evaluation request without writing output."""
        runtime = _parse_evaluation_options(request.options)
        contrast = _parse_contrast(request.contrast)
        _require_file(request.protein_matrix, "protein_matrix")
        _require_file(request.sdrf, "sdrf")
        if runtime.ground_truth is not None:
            _require_file(runtime.ground_truth, "ground_truth")
        destination = _output_path(request.output_dir)
        frame = _read_matrix(request.protein_matrix)
        profile, conditions, context = self._profile(
            frame,
            request.sdrf,
            runtime.metadata,
        )
        _validate_contrast_samples(frame, conditions, contrast)
        generated = validate_generated_recommendation(
            request.recommendation,
            context,
        )
        return PreparedEvaluation(
            runtime,
            contrast,
            destination,
            frame,
            profile,
            conditions,
            context,
            generated,
        )

    def _profile(
        self,
        frame: pd.DataFrame,
        sdrf: str,
        metadata: DatasetMetadata,
    ) -> tuple[DataProfile, dict[str, str], BoundContext]:
        """Build the shared profile and policy context for one request."""
        conditions = detect_condition_from_sdrf(sdrf, metadata.factor_column)
        missing = [
            str(sample) for sample in frame.columns[1:] if str(sample) not in conditions
        ]
        if missing:
            raise ValueError(
                "Protein-matrix samples are missing from the SDRF mapping: "
                + ", ".join(missing)
            )
        profile = profile_data(frame, conditions, metadata=metadata.acquisition)
        return profile, conditions, bind_context(profile, self._graph)


def _abstention_payload(
    graph: KnowledgeGraph,
    prepared: PreparedEvaluation,
) -> dict[str, Any]:
    """Build the persisted result for a policy-required abstention."""
    payload = {
        "status": "abstained",
        "abstain_reason": prepared.generated.abstain_reason,
        "ranking_objective": None,
        "results": [],
        "ranking": [],
        "best_config": None,
        "cache": {"hits": 0, "misses": 0, "unique_combos": 0},
    }
    payload.update(
        _audit_fields(
            graph,
            prepared.generated,
            prepared.context,
            input_scale=None,
        )
    )
    return payload


def _run_candidates(
    candidates: list[CandidateConfig],
    protein_df: pd.DataFrame,
    run: EvaluationRun,
) -> dict[str, Any]:
    """Run a bounded candidate set and persist its exact result tables."""
    slugs = _candidate_slugs(candidates)
    run.output_dir.mkdir(parents=True, exist_ok=True)
    cache = PreprocessCache(run.experiment.threads)
    evaluations: list[EvaluationResult] = []
    de_tables: dict[str, pd.DataFrame] = {}
    for config, filename in zip(candidates, slugs):
        measured, de_table = _evaluate_candidate(config, protein_df, run, cache)
        evaluations.append(measured)
        de_tables[config.name] = de_table
        de_table.to_csv(
            run.output_dir / f"{filename}.de.tsv",
            sep="\t",
            index=False,
        )

    sensitivity_table, sensitivity_summary = method_sensitivity(
        de_tables,
        run.expected_direction if run.ground_truth is not None else None,
    )
    sensitivity_table.to_csv(
        run.output_dir / "method_sensitivity.tsv",
        sep="\t",
        index=False,
    )
    rows = [measured.to_dict() for measured in evaluations]
    payload = _ranking_payload(rows, run.ground_truth is not None, cache.stats())
    payload["method_sensitivity"] = sensitivity_summary
    payload["method_sensitivity_artifact"] = "method_sensitivity.tsv"
    return payload


def _evaluate_candidate(
    config: CandidateConfig,
    protein_df: pd.DataFrame,
    run: EvaluationRun,
    cache: PreprocessCache,
) -> tuple[EvaluationResult, pd.DataFrame]:
    """Execute and measure one candidate against the shared round context."""
    de_table = run_experiment(config, protein_df, run.experiment, cache=cache)
    processed = cache.get_or_compute(
        config.normalization,
        config.imputation,
        protein_df,
    )
    measured = evaluate(
        config,
        de_table,
        processed,
        run.experiment.sample_to_condition,
        (
            (run.ground_truth, run.expected_direction)
            if run.ground_truth is not None and run.expected_direction is not None
            else None
        ),
    )
    if run.ground_truth is not None:
        measured = measured.with_score_a(compute_score_ground_truth(measured))
    return measured, de_table


def _candidate_slugs(candidates: list[CandidateConfig]) -> list[str]:
    """Validate candidate identities before creating partial output."""
    names = [candidate.name for candidate in candidates]
    if len(names) != len(set(names)):
        raise ValueError("Candidate names must be unique")
    slugs = [_slug(name) for name in names]
    if len(slugs) != len(set(slugs)):
        raise ValueError("Candidate names must map to distinct output filenames")
    return slugs


def _ranking_payload(
    rows: list[dict[str, Any]],
    has_truth: bool,
    cache: dict[str, int],
) -> dict[str, Any]:
    """Apply the Score A-only ranking contract to measured rows."""
    ranking = (
        sorted(
            (row for row in rows if row["score_a"] is not None),
            key=lambda row: row["score_a"],
            reverse=True,
        )
        if has_truth
        else []
    )
    status = "ranked" if ranking else "ground_truth_unscored"
    if not has_truth:
        status = "exploratory_unranked"
    return {
        "status": status,
        "ranking_objective": "score_a" if has_truth else None,
        "results": rows,
        "ranking": ranking,
        "best_config": ranking[0]["config_name"] if ranking else None,
        "cache": cache,
    }


def _candidate(
    item: dict[str, Any],
    generated: GeneratedRecommendationBlock,
    fdr_threshold: float,
) -> CandidateConfig:
    """Attach block-level provenance to one validated executable config."""
    return CandidateConfig(
        **item,
        fdr_threshold=fdr_threshold,
        evidence_refs=list(generated.evidence_refs),
        confidence=generated.confidence,
        limitations=list(generated.limitations),
        generated_by="host",
    )


def _policy_recommendation(
    profile: DataProfile,
    context: BoundContext,
) -> dict[str, Any]:
    """Return an exact recommendation block from deterministic policy."""
    if not context.scope.generation_allowed:
        errors = [
            item.message for item in context.diagnostics if item.severity == "error"
        ]
        block = {
            "configs": [],
            "evidence_refs": [],
            "confidence": "low",
            "limitations": errors,
            "abstain_reason": "; ".join(errors),
        }
        validate_generated_recommendation(block, context)
        return block

    candidates = rule_propose(profile, context)
    evidence_refs = list(
        dict.fromkeys(
            reference
            for candidate in candidates
            for reference in candidate.evidence_refs
        )
    )
    limitations = list(
        dict.fromkeys(
            limitation
            for candidate in candidates
            for limitation in candidate.limitations
        )
    )
    block = {
        "configs": [
            {field: candidate.to_dict()[field] for field in GENERATED_CONFIG_FIELDS}
            for candidate in candidates
        ],
        "evidence_refs": evidence_refs,
        "confidence": "low",
        "limitations": limitations,
        "abstain_reason": None,
    }
    validate_generated_recommendation(block, context)
    return block


def _parse_dataset_metadata(
    payload: dict[str, Any] | None,
) -> DatasetMetadata:
    """Validate declared acquisition metadata from an MCP object."""
    values = _validated_object(payload, _METADATA_FIELDS, "metadata")
    for name, value in values.items():
        if value is not None and (not isinstance(value, str) or not value.strip()):
            raise ValueError(f"metadata.{name} must be a non-empty string or null")
    return DatasetMetadata(
        acquisition=AcquisitionMetadata(
            values.get("data_type"),
            values.get("quantification"),
            values.get("upstream_engine"),
        ),
        factor_column=values.get("factor_column"),
    )


def _parse_evaluation_options(
    payload: dict[str, Any] | None,
) -> EvaluationOptions:
    """Validate the optional execution block from the MCP boundary."""
    values = _validated_object(payload, _OPTION_FIELDS, "options")
    metadata = _parse_dataset_metadata(
        {name: values[name] for name in _METADATA_FIELDS if name in values}
    )
    ground_truth = values.get("ground_truth")
    if ground_truth is not None and not isinstance(ground_truth, str):
        raise ValueError("options.ground_truth must be a string or null")
    expected_direction = values.get("expected_direction")
    fdr_threshold = values.get("fdr_threshold", 0.05)
    input_scale = values.get("input_scale", "auto")
    threads = values.get("threads", 24)
    if isinstance(fdr_threshold, bool) or not isinstance(fdr_threshold, (int, float)):
        raise ValueError("options.fdr_threshold must be numeric")
    _validate_runtime_options(
        float(fdr_threshold),
        input_scale,
        threads,
        expected_direction,
    )
    if ground_truth is not None and expected_direction is None:
        raise ValueError(
            "options.expected_direction is required when ground_truth is supplied"
        )
    if ground_truth is None and expected_direction is not None:
        raise ValueError(
            "options.expected_direction must be null when ground_truth is null"
        )
    return EvaluationOptions(
        metadata=metadata,
        ground_truth=ground_truth,
        expected_direction=expected_direction,
        fdr_threshold=float(fdr_threshold),
        input_scale=input_scale,
        threads=threads,
    )


def _validated_object(
    payload: dict[str, Any] | None,
    allowed_fields: frozenset[str],
    name: str,
) -> dict[str, Any]:
    """Require an object with no fields outside its public contract."""
    if payload is None:
        return {}
    if not isinstance(payload, dict):
        raise ValueError(f"{name} must be an object or null")
    unknown = set(payload) - allowed_fields
    if unknown:
        raise ValueError(f"Unknown {name} fields: {sorted(unknown)}")
    return payload


def _parse_contrast(contrast: list[str]) -> tuple[str, str]:
    """Require exactly two distinct, non-empty condition labels."""
    if not isinstance(contrast, list) or len(contrast) != 2:
        raise ValueError("contrast must be a two-item list")
    if not all(isinstance(item, str) and item.strip() for item in contrast):
        raise ValueError("contrast entries must be non-empty strings")
    if contrast[0] == contrast[1]:
        raise ValueError("contrast conditions must be distinct")
    return contrast[0], contrast[1]


def _read_matrix(path: str) -> pd.DataFrame:
    """Read a comma- or tab-delimited protein matrix with numeric sample cells."""
    matrix_path = Path(path).expanduser().resolve()
    with matrix_path.open(encoding="utf-8-sig") as handle:
        header = handle.readline()
    separator = "\t" if header.count("\t") > header.count(",") else ","
    columns = next(csv.reader([header], delimiter=separator))
    if any(not column.strip() for column in columns):
        raise ValueError("Protein matrix column names must be non-empty")
    if len(columns) != len(set(columns)):
        raise ValueError("Protein matrix column names must be unique")
    frame = pd.read_csv(matrix_path, sep=separator)
    if frame.shape[1] < 3:
        raise ValueError(
            "Protein matrix requires an identifier and at least two samples"
        )
    identifiers = frame.iloc[:, 0]
    if identifiers.isna().any() or identifiers.astype(str).str.strip().eq("").any():
        raise ValueError("Protein identifiers must be non-empty")
    identifiers = identifiers.astype(str)
    if identifiers.duplicated().any():
        raise ValueError("Protein identifiers must be unique")
    frame[frame.columns[0]] = identifiers
    numeric = frame.iloc[:, 1:].apply(pd.to_numeric, errors="raise")
    values = numeric.to_numpy(dtype=float)
    if np.isinf(values).any():
        raise ValueError("Protein matrix contains an infinite intensity")
    if not np.isfinite(values).any():
        raise ValueError("Protein matrix contains no finite intensities")
    frame[numeric.columns] = numeric
    return frame


def _require_file(path: str, name: str) -> None:
    """Require absolute file paths at the MCP boundary."""
    candidate = Path(path).expanduser()
    if not candidate.is_absolute():
        raise ValueError(f"{name} must be an absolute path")
    if not candidate.is_file():
        raise FileNotFoundError(f"{name} not found: {candidate}")


def _output_path(path: str) -> Path:
    """Resolve one required absolute output directory."""
    destination = Path(path).expanduser()
    if not destination.is_absolute():
        raise ValueError("output_dir must be an absolute path")
    destination = destination.resolve()
    if destination.exists():
        raise FileExistsError(f"output_dir already exists: {destination}")
    return destination


def _validate_contrast_samples(
    frame: pd.DataFrame,
    conditions: dict[str, str],
    contrast: tuple[str, str],
) -> None:
    """Require two mapped matrix replicates for each requested condition."""
    counts = {
        condition: sum(
            conditions[str(sample)] == condition for sample in frame.columns[1:]
        )
        for condition in contrast
    }
    if any(count < 2 for count in counts.values()):
        raise ValueError(
            "Contrast requires at least two matrix samples per condition: "
            + ", ".join(f"{condition}={count}" for condition, count in counts.items())
        )


def _as_linear(
    frame: pd.DataFrame,
    requested: str,
    inferred_log2: bool,
) -> tuple[pd.DataFrame, str]:
    """Return the linear-intensity representation required by the Rust kernel."""
    resolved = "log2" if requested == "auto" and inferred_log2 else requested
    resolved = "linear" if resolved == "auto" else resolved
    if resolved == "linear":
        return frame, resolved
    converted = frame.copy()
    with np.errstate(over="raise", invalid="raise"):
        converted.iloc[:, 1:] = np.exp2(converted.iloc[:, 1:].astype(float))
    return converted, resolved


def _load_ground_truth(path: str | None) -> set[str] | None:
    """Read one protein identifier per line when Score A truth is available."""
    if path is None:
        return None
    values = {
        line.strip()
        for line in Path(path).expanduser().read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    if not values:
        raise ValueError("Ground-truth file contains no protein identifiers")
    return values


def _validate_runtime_options(
    fdr_threshold: float,
    input_scale: Any,
    threads: Any,
    expected_direction: Any,
) -> None:
    """Validate user-controlled execution options before writing output."""
    if not 0.0 < fdr_threshold <= 1.0:
        raise ValueError("options.fdr_threshold must be in (0, 1]")
    if input_scale not in {"auto", "linear", "log2"}:
        raise ValueError("options.input_scale must be auto, linear, or log2")
    if isinstance(threads, bool) or not isinstance(threads, int):
        raise ValueError("options.threads must be an integer between 1 and 256")
    if not 1 <= threads <= 256:
        raise ValueError("options.threads must be an integer between 1 and 256")
    if expected_direction is not None and expected_direction not in {"UP", "DOWN"}:
        raise ValueError("options.expected_direction must be UP, DOWN, or null")


def _audit_fields(
    graph: KnowledgeGraph,
    generated: GeneratedRecommendationBlock,
    context: BoundContext,
    input_scale: str | None,
) -> dict[str, Any]:
    """Return provenance shared by completed and abstained evaluations."""
    return {
        "knowledge_fingerprint": graph.fingerprint,
        "evidence_refs": list(generated.evidence_refs),
        "confidence": generated.confidence,
        "limitations": list(generated.limitations),
        "input_scale": input_scale,
        "diagnostics": [item.to_dict() for item in context.diagnostics],
    }


def _write_evaluation(destination: Path, payload: dict[str, Any]) -> dict[str, Any]:
    """Persist and return one strict JSON-safe evaluation artifact."""
    destination.mkdir(parents=True, exist_ok=True)
    safe_payload = _json_safe(payload)
    (destination / "evaluation.json").write_text(
        json.dumps(safe_payload, indent=2, sort_keys=True, allow_nan=False),
        encoding="utf-8",
    )
    return safe_payload


def _slug(value: str) -> str:
    """Create a portable filename for one candidate result."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._") or "candidate"


def _json_safe(value: Any) -> Any:
    """Replace non-finite numeric diagnostics before JSON/MCP serialization."""
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return None
    if isinstance(value, np.generic):
        return value.item()
    return value
