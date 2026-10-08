"""Preflight review and post-run insight contracts for Mokume Studio."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import httpx
import pytest
from studio_test_support import make_studio_app

from mokume.studio.models import ArtifactRecord, JobSpec, RunStatus, utc_now


PORT = 18767
ORIGIN = f"http://127.0.0.1:{PORT}"
TOKEN = "studio-insights-token"
pytestmark = pytest.mark.anyio


@pytest.fixture(name="studio_client")
async def build_studio_client(tmp_path):
    """Create an isolated Studio app for insight endpoint tests."""
    app = make_studio_app(TOKEN, tmp_path / "state")
    transport = httpx.ASGITransport(app=app)
    async with httpx.AsyncClient(transport=transport, base_url=ORIGIN) as client:
        client.app = app
        yield client


async def _open_project(client: httpx.AsyncClient, project: Path) -> dict[str, str]:
    response = await client.get(f"/?token={TOKEN}", follow_redirects=False)
    assert response.status_code == 303
    session = (await client.get("/api/session")).json()
    headers = {"Origin": ORIGIN, "X-CSRF-Token": session["csrf_token"]}
    opened = await client.post(
        "/api/projects/open",
        json={"path": str(project)},
        headers=headers,
    )
    assert opened.status_code == 200
    return headers


def _review_argv() -> list[str]:
    return [
        "quantify",
        "features2proteins",
        "--parquet",
        "input.parquet",
        "--sdrf",
        "design.sdrf.tsv",
        "--quant-method",
        "sum",
        "--sample-normalization",
        "condition-median",
        "--output",
        "result.csv",
        "--threads",
        "24",
    ]


async def _assert_review_failure_cases(
    client: httpx.AsyncClient, headers: dict[str, str], project: Path
) -> None:
    """Verify existing outputs and missing SDRF conditions are rejected."""
    (project / "result.csv").write_text("already here\n", encoding="utf-8")
    blocked = await client.post(
        "/api/commands/review",
        json={"argv": _review_argv()},
        headers=headers,
    )
    output = next(
        item for item in blocked.json()["outputs"] if item["parameter"] == "output"
    )
    assert blocked.json()["valid"] is False
    assert output["status"] == "error"
    assert output["message"] == "Output already exists"

    (project / "design.sdrf.tsv").write_text(
        "source name\tcomment[data file]\nS1\tplex-1\n",
        encoding="utf-8",
    )
    missing_condition = await client.post(
        "/api/commands/review",
        json={"argv": _review_argv()},
        headers=headers,
    )
    condition_check = next(
        check
        for check in missing_condition.json()["checks"]
        if check["id"] == "sdrf:column:condition"
    )
    assert condition_check["status"] == "error"


async def test_review_reports_mapping_resources_and_existing_output(
    studio_client, tmp_path
):
    """Preflight review stays non-mutating and points to blocking parameters."""
    project = tmp_path / "project"
    project.mkdir()
    (project / "input.parquet").write_bytes(b"PAR1")
    (project / "design.sdrf.tsv").write_text(
        "source name\tfactor value[condition]\t"
        "characteristics[biological replicate]\tcomment[data file]\n"
        "S1\tcontrol\t1\tplex-1\n"
        "S2\ttreated\t2\tplex-1\n",
        encoding="utf-8",
    )
    headers = await _open_project(studio_client, project)

    review = await studio_client.post(
        "/api/commands/review",
        json={"argv": _review_argv()},
        headers=headers,
    )

    assert review.status_code == 200
    payload = review.json()
    assert payload["valid"] is True
    assert payload["resources"]["threads"] == 24
    assert payload["resources"]["available_memory_bytes"] > 0
    assert payload["resources"]["free_disk_bytes"] > 0
    assert payload["planned_steps"] == ["read", "aggregate", "normalize", "export"]
    assert payload["sdrf"]["row_count"] == 2
    assert payload["sdrf"]["columns"]["condition"] == "factor value[condition]"
    assert payload["sdrf"]["rows"][1]["replicate"] == "2"
    assert not (project / "result.csv").exists()
    await _assert_review_failure_cases(studio_client, headers, project)


def _completed_run_files(project: Path) -> tuple[Path, Path, Path]:
    """Create the input, output, and run directory for a completed run."""
    run_directory = project / "results" / "mokume" / "run-1"
    run_directory.mkdir(parents=True)
    input_path = project / "peptides.csv"
    output_path = project / "proteins.csv"
    input_path.write_text("peptide,S1,S2\nAA,1,2\n", encoding="utf-8")
    output_path.write_text(
        "protein_id,S1,S2,S3\nP1,10,12,11\nP2,20,21,19\nP3,5,,6\n",
        encoding="utf-8",
    )
    return run_directory, input_path, output_path


def _completed_run_argv(input_path: Path, output_path: Path) -> list[str]:
    """Build the command recorded by the completed run fixture."""
    return [
        "quantify",
        "peptides2protein",
        "--peptides",
        str(input_path),
        "--quant-method",
        "sum",
        "--output",
        str(output_path),
    ]


def _write_run_audit(
    run_directory: Path,
    parameters: dict,
    input_path: Path,
    artifact: ArtifactRecord,
) -> None:
    """Persist the fixture's parameters, provenance, and logs."""
    (run_directory / "parameters.json").write_text(
        json.dumps(parameters), encoding="utf-8"
    )
    (run_directory / "provenance.json").write_text(
        json.dumps(
            {
                "mokume_version": "test",
                "python_version": "3.test",
                "platform": "test-platform",
                "inputs": [{"path": str(input_path), "sha256": "input-hash"}],
                "artifacts": [artifact.model_dump(mode="json")],
            }
        ),
        encoding="utf-8",
    )
    (run_directory / "stdout.log").write_text("completed\n", encoding="utf-8")
    (run_directory / "stderr.log").write_text("", encoding="utf-8")


async def _seed_completed_run(
    client: httpx.AsyncClient, project: Path
) -> tuple[dict[str, str], str]:
    """Create one succeeded run with a complete local audit record."""
    run_directory, input_path, output_path = _completed_run_files(project)
    headers = await _open_project(client, project)
    store = client.app.state.runtime.store
    project_record = store.active_project()
    argv = _completed_run_argv(input_path, output_path)
    parameters = {"operation": "native", "argv": argv, "threads": 24}
    spec = JobSpec(
        run_id="run-1",
        project_root=str(project),
        run_directory=str(run_directory),
        argv=argv,
        parameters=parameters,
        approved_hash="approved-hash",
        created_at=utc_now(),
    )
    store.create_run(spec, project_record.id, "quantify peptides2protein")
    store.update_run("run-1", RunStatus.RUNNING)
    for stage in ("inputs", "workflow", "artifacts", "provenance"):
        store.add_event(
            "run-1",
            "stage",
            {"stage": stage, "status": "succeeded", "elapsed_seconds": 0.25},
        )
    store.update_run("run-1", RunStatus.SUCCEEDED)
    digest = hashlib.sha256(output_path.read_bytes()).hexdigest()
    artifact = ArtifactRecord(
        id="artifact-1",
        run_id="run-1",
        path=str(output_path),
        media_type="text/csv",
        size=output_path.stat().st_size,
        sha256=digest,
    )
    store.register_artifact(artifact)
    _write_run_audit(run_directory, parameters, input_path, artifact)
    return headers, digest


async def test_run_details_include_qc_stages_logs_and_provenance(
    studio_client, tmp_path
):
    """One successful run exposes its complete workspace-scoped audit record."""
    project = tmp_path / "project"
    headers, digest = await _seed_completed_run(studio_client, project)

    response = await studio_client.get("/api/runs/run-1/details")

    assert response.status_code == 200
    details = response.json()
    assert details["run"]["status"] == "succeeded"
    assert details["template"]["workflow"] == ["quantify", "peptides2protein"]
    assert details["template"]["parameters"]["quant-method"] == "sum"
    assert [stage["status"] for stage in details["stages"]] == ["succeeded"] * 4
    assert details["logs"]["stdout"] == ["completed"]
    assert details["provenance"]["inputs"][0]["sha256"] == "input-hash"
    assert details["artifacts"][0]["sha256"] == digest
    assert details["qc"]["available"] is True
    assert details["qc"]["entity_counts"] == {"proteins": 3}
    assert details["qc"]["sample_count"] == 3
    assert details["qc"]["missing_percent"] > 0

    other = tmp_path / "other"
    other.mkdir()
    switched = await studio_client.post(
        "/api/projects/open",
        json={"path": str(other)},
        headers=headers,
    )
    assert switched.status_code == 200
    assert (await studio_client.get("/api/runs/run-1/details")).status_code == 404
    assert (await studio_client.get("/api/artifacts/artifact-1")).status_code == 404
    assert (await studio_client.get("/api/artifacts")).json() == {"artifacts": []}
