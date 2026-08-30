"""Tests for deterministic Studio science control records."""

from __future__ import annotations

import json
import sqlite3

import pytest
from pydantic import ValidationError

from mokume.studio.science import (
    ApprovalStatus,
    DatasetInspectionRequest,
    DatasetStatus,
    ScienceStore,
)
from mokume.studio.state import StateStore


def _stores(tmp_path):
    state = StateStore(tmp_path / "state")
    return state, ScienceStore(state)


def _request(**overrides) -> DatasetInspectionRequest:
    values = {
        "protein_matrix": "inputs/proteins.tsv",
        "sdrf": "inputs/metadata.sdrf.tsv",
        "contrast": ["treated", "control"],
        "input_scale": "log2",
        "data_type": "LFQ",
    }
    values.update(overrides)
    return DatasetInspectionRequest(**values)


@pytest.mark.parametrize(
    "contrast",
    [
        ["treated"],
        ["treated", "control", "vehicle"],
        ["treated", "treated"],
        ["treated", " treated "],
        ["", "control"],
    ],
)
def test_inspection_request_requires_two_distinct_contrast_values(contrast):
    """Reject missing, extra, blank, and duplicate contrast values."""
    with pytest.raises(ValidationError):
        _request(contrast=contrast)


def test_inspection_request_forbids_unknown_fields():
    """Reject fields outside the deterministic inspection contract."""
    with pytest.raises(ValidationError):
        _request(api_key="must-not-be-stored")


def test_dataset_round_trip_and_latest(tmp_path):
    """Persist canonical requests and results in the shared state database."""
    state, science = _stores(tmp_path)
    project = state.open_project(str(tmp_path / "project"))
    first = science.create_dataset(project.id, _request())

    assert first.status is DatasetStatus.QUEUED
    assert first.request.output_directory == "results/mokume"
    ready = science.update_dataset(
        first.id,
        DatasetStatus.READY,
        result={"z": 2, "a": {"ready": True}},
    )
    assert ready.result == {"a": {"ready": True}, "z": 2}
    assert science.get_dataset(first.id) == ready
    assert science.latest_dataset(project.id) == ready

    with sqlite3.connect(science.database) as connection:
        stored = connection.execute(
            "SELECT result_json FROM datasets WHERE id=?", (first.id,)
        ).fetchone()[0]
    assert stored == '{"a":{"ready":true},"z":2}'


def test_restart_marks_incomplete_dataset_inspections_failed(tmp_path):
    """A restart cannot leave the UI polling a dataset forever."""
    state, science = _stores(tmp_path)
    project = state.open_project(str(tmp_path / "project"))
    queued = science.create_dataset(project.id, _request())
    running = science.create_dataset(project.id, _request())
    science.update_dataset(running.id, DatasetStatus.RUNNING)

    assert science.interrupt_incomplete_datasets() == 2
    for dataset_id in (queued.id, running.id):
        stored = science.get_dataset(dataset_id)
        assert stored is not None
        assert stored.status is DatasetStatus.FAILED
        assert stored.error == "Studio restarted during dataset inspection"


def test_approval_hash_and_one_time_transitions(tmp_path):
    """Bind approval to its payload and permit each transition only once."""
    _, science = _stores(tmp_path)
    approval = science.create_approval(
        "evaluate_recommendation",
        {"project_id": "project-1", "parameters": {"z": 2, "a": 1}},
    )

    with pytest.raises(ValueError, match="hash mismatch"):
        science.decide_approval(
            approval.id,
            payload_hash="0" * 64,
            approved=True,
        )
    approved = science.decide_approval(
        approval.id,
        payload_hash=approval.payload_hash,
        approved=True,
    )
    assert approved.status is ApprovalStatus.APPROVED
    with pytest.raises(ValueError, match="already approved"):
        science.decide_approval(
            approval.id,
            payload_hash=approval.payload_hash,
            approved=True,
        )

    consumed = science.consume_approval(
        approval.id,
        payload_hash=approval.payload_hash,
    )
    assert consumed.status is ApprovalStatus.CONSUMED
    linked = science.link_run(
        approval.id,
        payload_hash=approval.payload_hash,
        run_id="run-1",
    )
    assert linked.run_id == "run-1"
    with pytest.raises(ValueError, match="cannot be linked"):
        science.link_run(
            approval.id,
            payload_hash=approval.payload_hash,
            run_id="run-2",
        )


def test_expired_approval_cannot_be_decided(tmp_path):
    """Persist expiry and refuse a decision after the thirty-minute window."""
    _, science = _stores(tmp_path)
    approval = science.create_approval("evaluate", {"dataset_id": "dataset-1"})
    with sqlite3.connect(science.database) as connection:
        connection.execute(
            "UPDATE approvals SET expires_at=? WHERE id=?",
            ("2000-01-01T00:00:00+00:00", approval.id),
        )

    with pytest.raises(ValueError, match="expired"):
        science.decide_approval(
            approval.id,
            payload_hash=approval.payload_hash,
            approved=True,
        )
    expired = science.get_approval(approval.id)
    assert expired is not None
    assert expired.status is ApprovalStatus.EXPIRED


def test_thread_cannot_cross_project_or_mode(tmp_path):
    """Keep conversation records isolated by project and agent mode."""
    state, science = _stores(tmp_path)
    first = state.open_project(str(tmp_path / "first"))
    second = state.open_project(str(tmp_path / "second"))
    messages = [{"role": "user", "content": "Explain missingness", "z": 2}]
    saved = science.save_thread(
        "thread-1",
        project_id=first.id,
        mode="ask",
        messages=messages,
    )

    assert saved.messages == messages
    assert saved.messages_json == json.dumps(
        messages,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    )
    assert science.get_thread("thread-1", project_id=first.id, mode="ask") == saved
    with pytest.raises(ValueError, match="different project or mode"):
        science.get_thread("thread-1", project_id=second.id, mode="ask")
    with pytest.raises(ValueError, match="different project or mode"):
        science.save_thread(
            "thread-1",
            project_id=first.id,
            mode="plan-run",
            messages=[],
        )
