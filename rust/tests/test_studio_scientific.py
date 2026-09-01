"""Scientific worker and approval contracts for Mokume Studio."""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from mokume.agentic.service import InspectionRequest, RecommendationService
from mokume.studio.jobs import JobConflictError, JobManager
from mokume.studio.models import (
    ProjectRecord,
    RunRecord,
    RunStatus,
    TERMINAL_RUN_STATUSES,
)
from mokume.studio.science import (
    ApprovalRecord,
    ApprovalStatus,
    DatasetInspectionRequest,
    DatasetRecord,
    DatasetStatus,
    ScienceStore,
)
from mokume.studio.scientific import EvaluationPlanRequest, ScientificController
from mokume.studio.state import StateStore


@dataclass
class ScientificHarness:
    """Shared stores, controller, and one mutable project input."""

    state: StateStore
    project: ProjectRecord
    manager: JobManager
    science: ScienceStore
    controller: ScientificController
    matrix: Path


def _project_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    project = tmp_path / "project"
    project.mkdir()
    proteins = [f"P{index}" for index in range(30)]
    baseline = [100.0 + index for index in range(30)]
    matrix = pd.DataFrame({"ProteinName": proteins})
    for sample, multiplier in {
        "A1": 1.0,
        "A2": 1.01,
        "A3": 0.99,
        "B1": 1.4,
        "B2": 1.41,
        "B3": 1.39,
    }.items():
        matrix[sample] = [value * multiplier for value in baseline]
    matrix_path = project / "proteins.tsv"
    matrix.to_csv(matrix_path, sep="\t", index=False)
    sdrf_path = project / "samples.sdrf.tsv"
    pd.DataFrame(
        {
            "source name": ["A1", "A2", "A3", "B1", "B2", "B3"],
            "factor value[condition]": ["A", "A", "A", "B", "B", "B"],
        }
    ).to_csv(sdrf_path, sep="\t", index=False)
    return project, matrix_path, sdrf_path


def _wait_for_terminal(
    store: StateStore, run_id: str, timeout: float = 30
) -> RunRecord:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        record = store.get_run(run_id)
        if record and record.status in TERMINAL_RUN_STATUSES:
            return record
        time.sleep(0.05)
    raise TimeoutError(f"Studio run did not finish: {run_id}")


def _start_after_worker_release(
    controller: ScientificController,
    approval_id: str,
    payload_hash: str,
    project: ProjectRecord,
    timeout: float = 5,
) -> RunRecord:
    """Retry only the brief interval before the prior worker is reaped."""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            return controller.start_approved(approval_id, payload_hash, project)
        except JobConflictError:
            time.sleep(0.05)
    raise TimeoutError("Inspection worker slot was not released")


def _request() -> DatasetInspectionRequest:
    return DatasetInspectionRequest(
        protein_matrix="proteins.tsv",
        sdrf="samples.sdrf.tsv",
        contrast=("A", "B"),
        input_scale="linear",
        data_type="LFQ",
        quantification="directlfq",
    )


def _mapping(value: Any) -> dict[str, Any]:
    assert isinstance(value, dict)
    return value


def _harness(tmp_path: Path) -> ScientificHarness:
    project_root, matrix, _sdrf = _project_inputs(tmp_path)
    state = StateStore(tmp_path / "state")
    project = state.open_project(str(project_root.resolve()))
    manager = JobManager(state)
    science = ScienceStore(state)
    return ScientificHarness(
        state,
        project,
        manager,
        science,
        ScientificController(science, manager),
        matrix,
    )


def _policy_config(dataset: DatasetRecord) -> dict[str, Any]:
    result = _mapping(dataset.result)
    policy = _mapping(result.get("policy_recommendation"))
    return policy["configs"][0]


def _ready_dataset(harness: ScientificHarness) -> DatasetRecord:
    request = _request().model_copy(
        update={
            "protein_matrix": str(harness.matrix),
            "sdrf": str(harness.matrix.parent / "samples.sdrf.tsv"),
        }
    )
    dataset = harness.science.create_dataset(harness.project.id, request)
    result = RecommendationService().inspect_dataset(
        InspectionRequest(
            protein_matrix=request.protein_matrix,
            sdrf=request.sdrf,
            input_scale=request.input_scale,
            contrast=list(request.contrast),
            metadata={"data_type": "LFQ", "quantification": "directlfq"},
        )
    )
    return harness.science.update_dataset(
        dataset.id,
        DatasetStatus.READY,
        result=result,
    )


def _evaluate_policy_config(
    harness: ScientificHarness,
    dataset: DatasetRecord,
    config_name: str,
) -> tuple[ApprovalRecord, RunRecord]:
    approval = harness.controller.prepare_evaluation(
        EvaluationPlanRequest(
            dataset_id=dataset.id,
            config_names=[config_name],
            output_directory="results/evaluation",
        ),
        harness.project,
        provider="openai",
        model="test-model",
    )
    harness.science.decide_approval(
        approval.id,
        payload_hash=approval.payload_hash,
        approved=True,
    )
    submitted = _start_after_worker_release(
        harness.controller,
        approval.id,
        approval.payload_hash,
        harness.project,
    )
    return approval, _wait_for_terminal(harness.state, submitted.id)


def _assert_unranked_evaluation(
    project_root: Path,
    config_name: str,
) -> None:
    evaluation = json.loads(
        (project_root / "results/evaluation/evaluation.json").read_text(
            encoding="utf-8"
        )
    )
    assert evaluation["status"] == "exploratory_unranked"
    assert evaluation["ranking_objective"] is None
    assert evaluation["ranking"] == []
    assert evaluation["best_config"] is None
    assert "winner" not in evaluation
    assert evaluation["results"][0]["config_name"] == config_name
    assert evaluation["results"][0]["score_a"] is None


def test_real_unlabelled_workflow_is_unranked_and_traceable(tmp_path):
    """Run inspect, policy selection, approval, and unlabelled evaluation."""
    harness = _harness(tmp_path)
    try:
        dataset, submitted = harness.controller.inspect(_request(), harness.project)
        inspection_run = _wait_for_terminal(harness.state, submitted.id)
        stored = harness.science.get_dataset(dataset.id)
        assert stored is not None
        config = _policy_config(stored)
        approval, evaluation_run = _evaluate_policy_config(
            harness,
            stored,
            config["name"],
        )
    finally:
        harness.manager.shutdown()

    assert inspection_run.status is RunStatus.SUCCEEDED
    assert evaluation_run.status is RunStatus.SUCCEEDED
    assert stored.status is DatasetStatus.READY
    result = _mapping(stored.result)
    assert _mapping(result.get("profile"))["n_samples"] == 6
    fingerprint = _mapping(result.get("context"))["knowledge_fingerprint"]
    assert fingerprint
    assert len(harness.state.list_artifacts(inspection_run.id)) == 1
    safe = harness.controller.context(harness.project, dataset.id)
    assert safe["disclosure"] == {"metadata": True, "raw_rows": False}
    assert safe["workspace"] == {
        "id": harness.project.id,
        "name": harness.matrix.parent.name,
        "root": harness.project.root,
        "access": "workspace_only",
    }
    safe_dataset = _mapping(safe.get("dataset"))
    assert "missing_per_sample" not in _mapping(safe_dataset.get("profile"))
    _assert_unranked_evaluation(harness.matrix.parent, config["name"])

    provenance = json.loads(
        (Path(evaluation_run.run_directory) / "provenance.json").read_text(
            encoding="utf-8"
        )
    )
    assert provenance["knowledge_fingerprint"] == fingerprint
    assert provenance["plan_source"] == {
        "provider": "openai",
        "model": "test-model",
    }
    linked = harness.science.get_approval(approval.id)
    assert linked is not None
    assert linked.status is ApprovalStatus.CONSUMED
    assert linked.run_id == evaluation_run.id


def test_plan_rejects_unknown_candidate_and_changed_input(tmp_path):
    """A model cannot invent configs or reuse approval after inputs change."""
    harness = _harness(tmp_path)
    dataset = _ready_dataset(harness)

    with pytest.raises(ValueError, match="outside Mokume policy"):
        harness.controller.prepare_evaluation(
            EvaluationPlanRequest(
                dataset_id=dataset.id,
                config_names=["invented"],
                output_directory="results/invented",
            ),
            harness.project,
            provider="openai",
            model="test",
        )

    approval = harness.controller.prepare_evaluation(
        EvaluationPlanRequest(
            dataset_id=dataset.id,
            config_names=[_policy_config(dataset)["name"]],
            output_directory="results/evaluation",
        ),
        harness.project,
        provider="openai",
        model="test",
    )
    harness.science.decide_approval(
        approval.id,
        payload_hash=approval.payload_hash,
        approved=True,
    )
    harness.matrix.write_text("changed after approval\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="changed after the run was approved"):
        harness.controller.start_approved(
            approval.id,
            approval.payload_hash,
            harness.project,
        )
    stored = harness.science.get_approval(approval.id)
    assert stored is not None
    assert stored.status is ApprovalStatus.APPROVED
    harness.manager.shutdown()
