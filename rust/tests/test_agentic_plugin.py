"""Public contract tests for the plugin-owned Mokume agentic workflow."""

from __future__ import annotations

import asyncio
import importlib
import json
from pathlib import Path
import shutil

import pandas as pd
import pytest

from mokume.agentic.evaluator import evaluate, method_sensitivity
from mokume.agentic.knowledge import load_knowledge_graph
from mokume.agentic.service import EvaluationRequest, RecommendationService
from mokume.agentic.state import CandidateConfig

REPOSITORY = Path(__file__).resolve().parents[2]
PLUGIN = REPOSITORY / "plugins" / "mokume"
KNOWLEDGE = PLUGIN / "knowledge" / "knowledge.yaml"
CODEX_MANIFEST = PLUGIN / ".codex-plugin" / "plugin.json"
CLAUDE_MANIFEST = PLUGIN / ".claude-plugin" / "plugin.json"
SERVICE = RecommendationService(str(KNOWLEDGE))


@pytest.fixture(name="lfq_inputs")
def _lfq_inputs(tmp_path: Path) -> dict[str, Path]:
    """Write a small 3-vs-3 LFQ matrix with five known changes."""
    proteins = [f"P{index}" for index in range(30)]
    baseline = [100.0 + index for index in range(30)]
    matrix = pd.DataFrame(
        {
            "ProteinName": proteins,
            "A1": baseline,
            "A2": [value * 1.01 for value in baseline],
            "A3": [value * 0.99 for value in baseline],
            "B1": [
                value * 4.0 if index < 5 else value
                for index, value in enumerate(baseline)
            ],
            "B2": [
                value * 4.04 if index < 5 else value * 1.01
                for index, value in enumerate(baseline)
            ],
            "B3": [
                value * 3.96 if index < 5 else value * 0.99
                for index, value in enumerate(baseline)
            ],
        }
    )
    protein_matrix = tmp_path / "proteins.tsv"
    matrix.to_csv(protein_matrix, sep="\t", index=False)
    sdrf = tmp_path / "samples.sdrf.tsv"
    pd.DataFrame(
        {
            "source name": ["A1", "A2", "A3", "B1", "B2", "B3"],
            "factor value[condition]": ["A", "A", "A", "B", "B", "B"],
        }
    ).to_csv(sdrf, sep="\t", index=False)
    truth = tmp_path / "truth.txt"
    truth.write_text("\n".join(proteins[:5]) + "\n", encoding="utf-8")
    return {"matrix": protein_matrix, "sdrf": sdrf, "truth": truth}


def _inspect(inputs: dict[str, Path]) -> dict:
    return SERVICE.inspect_dataset(
        str(inputs["matrix"]),
        str(inputs["sdrf"]),
        {
            "data_type": "LFQ",
            "quantification": "directlfq",
            "upstream_engine": "quantms",
        },
    )


def test_service_ranks_only_with_ground_truth(
    lfq_inputs: dict[str, Path],
    tmp_path: Path,
) -> None:
    """Score A selects a winner; an unlabelled dataset remains unranked."""
    inspected = _inspect(lfq_inputs)
    assert len(inspected["context"]["knowledge_fingerprint"]) == 64
    diagnostics = next(
        block["content"]
        for block in inspected["context"]["blocks"]
        if block["id"] == "diagnostics"
    )
    assert "OUTSIDE_BENCHMARK_PROFILE" in {
        diagnostic["code"] for diagnostic in diagnostics
    }
    recommendation = inspected["policy_recommendation"]
    assert recommendation["confidence"] == "low"
    recommendation["configs"] = recommendation["configs"][:1]

    ranked = SERVICE.evaluate_recommendation(
        EvaluationRequest(
            str(lfq_inputs["matrix"]),
            str(lfq_inputs["sdrf"]),
            ["A", "B"],
            recommendation,
            str(tmp_path / "ranked"),
            {
                "ground_truth": str(lfq_inputs["truth"]),
                "expected_direction": "UP",
                "data_type": "LFQ",
                "quantification": "directlfq",
                "upstream_engine": "quantms",
                "threads": 24,
            },
        )
    )
    assert ranked["status"] == "ranked"
    assert ranked["ranking_objective"] == "score_a"
    assert ranked["ranking"][0]["score_a"] is not None
    assert ranked["best_config"] == recommendation["configs"][0]["name"]
    assert (tmp_path / "ranked" / "evaluation.json").is_file()

    exploratory = SERVICE.evaluate_recommendation(
        EvaluationRequest(
            str(lfq_inputs["matrix"]),
            str(lfq_inputs["sdrf"]),
            ["A", "B"],
            recommendation,
            str(tmp_path / "exploratory"),
            {
                "data_type": "LFQ",
                "quantification": "directlfq",
                "upstream_engine": "quantms",
                "threads": 24,
            },
        )
    )
    assert exploratory["status"] == "exploratory_unranked"
    assert exploratory["ranking_objective"] is None
    assert exploratory["ranking"] == []
    assert exploratory["results"][0]["score_a"] is None
    assert exploratory["best_config"] is None
    assert exploratory["method_sensitivity"]["comparison_available"] is False
    assert exploratory["method_sensitivity_artifact"] == "method_sensitivity.tsv"
    assert "stability" not in exploratory
    assert (tmp_path / "exploratory" / "method_sensitivity.tsv").is_file()


def test_method_sensitivity_reports_signed_call_agreement() -> None:
    """Exploratory comparison reports sharing without creating a ranking."""
    first = pd.DataFrame(
        {
            "protein": ["P1", "P2", "P3"],
            "significance": ["UP", "DOWN", "NOT_DE"],
        }
    )
    second = pd.DataFrame(
        {
            "protein": ["P1", "P2", "P3"],
            "significance": ["UP", "UP", "DOWN"],
        }
    )

    table, summary = method_sensitivity({"first": first, "second": second})

    assert summary == {
        "comparison_available": True,
        "candidate_count": 2,
        "signed_call_union": 4,
        "shared_signed_calls": 1,
        "method_sensitive_signed_calls": 3,
        "interpretation": (
            "Call sharing describes sensitivity to the tested methods and is not "
            "evidence of biological truth."
        ),
    }
    calls = table.set_index(["protein", "direction"])
    assert calls.loc[("P1", "UP"), "classification"] == "shared"
    assert calls.loc[("P2", "DOWN"), "classification"] == "method_sensitive"
    assert calls.loc[("P2", "UP"), "classification"] == "method_sensitive"


def test_bundled_knowledge_rejects_modified_source_artifacts(tmp_path: Path) -> None:
    """Every bundled benchmark artifact is bound to its declared hash."""
    copied = tmp_path / "knowledge"
    shutil.copytree(PLUGIN / "knowledge", copied)
    report = copied / "sources" / "spike-in-score-a-320" / "benchmark_report.md"
    report.write_text(
        report.read_text(encoding="utf-8") + "\nmodified\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="Knowledge artifact hash mismatch"):
        load_knowledge_graph(copied / "knowledge.yaml")


def test_bundled_knowledge_rejects_nonfinite_profile_ranges(tmp_path: Path) -> None:
    """Applicability ranges cannot contain YAML NaN or infinities."""
    copied = tmp_path / "knowledge"
    shutil.copytree(PLUGIN / "knowledge", copied)
    knowledge = copied / "knowledge.yaml"
    contents = knowledge.read_text(encoding="utf-8").replace(
        "missing_rate: [0.033081, 0.057325]",
        "missing_rate: [.nan, 0.057325]",
        1,
    )
    knowledge.write_text(contents, encoding="utf-8")

    with pytest.raises(ValueError, match="finite.*numeric range"):
        load_knowledge_graph(knowledge)


def test_agentic_catalog_rejects_version_branches(tmp_path: Path) -> None:
    """The runtime catalog is updated in place rather than schema-versioned."""
    copied = tmp_path / "knowledge"
    shutil.copytree(PLUGIN / "knowledge", copied)
    knowledge = copied / "knowledge.yaml"
    knowledge.write_text(
        "version: 1\n" + knowledge.read_text(encoding="utf-8"),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="catalog fields must be exactly"):
        load_knowledge_graph(knowledge)


def test_evaluation_refuses_to_overwrite_existing_output(
    lfq_inputs: dict[str, Path],
    tmp_path: Path,
) -> None:
    """A repeated round must use a new output directory."""
    recommendation = _inspect(lfq_inputs)["policy_recommendation"]
    recommendation["configs"] = recommendation["configs"][:1]
    output = tmp_path / "existing"
    output.mkdir()

    with pytest.raises(FileExistsError, match="output_dir already exists"):
        SERVICE.evaluate_recommendation(
            EvaluationRequest(
                str(lfq_inputs["matrix"]),
                str(lfq_inputs["sdrf"]),
                ["A", "B"],
                recommendation,
                str(output),
                {
                    "data_type": "LFQ",
                    "quantification": "directlfq",
                    "upstream_engine": "quantms",
                    "threads": 24,
                },
            )
        )


def test_inspection_rejects_matrix_samples_absent_from_sdrf(
    lfq_inputs: dict[str, Path],
) -> None:
    """Every matrix sample must bind to authoritative SDRF metadata."""
    matrix = pd.read_csv(lfq_inputs["matrix"], sep="\t").rename(
        columns={"B3": "UNMAPPED"}
    )
    matrix.to_csv(lfq_inputs["matrix"], sep="\t", index=False)

    with pytest.raises(ValueError, match="missing from the SDRF mapping: UNMAPPED"):
        _inspect(lfq_inputs)


def test_generated_recommendation_requires_complete_unique_configs(
    lfq_inputs: dict[str, Path],
    tmp_path: Path,
) -> None:
    """Host-generated candidates need rationale and unique executable settings."""
    recommendation = _inspect(lfq_inputs)["policy_recommendation"]
    first = recommendation["configs"][0]
    first["reasoning"] = ""
    output = tmp_path / "blank-reasoning"

    with pytest.raises(ValueError, match="reasoning must be a non-empty string"):
        SERVICE.evaluate_recommendation(
            EvaluationRequest(
                str(lfq_inputs["matrix"]),
                str(lfq_inputs["sdrf"]),
                ["A", "B"],
                recommendation,
                str(output),
                {
                    "data_type": "LFQ",
                    "quantification": "directlfq",
                    "upstream_engine": "quantms",
                    "threads": 24,
                },
            )
        )
    assert not output.exists()

    recommendation = _inspect(lfq_inputs)["policy_recommendation"]
    duplicate = dict(recommendation["configs"][0])
    duplicate["name"] = "same-settings-under-another-name"
    recommendation["configs"] = [recommendation["configs"][0], duplicate]
    output = tmp_path / "duplicate-settings"

    with pytest.raises(ValueError, match="unique executable settings"):
        SERVICE.evaluate_recommendation(
            EvaluationRequest(
                str(lfq_inputs["matrix"]),
                str(lfq_inputs["sdrf"]),
                ["A", "B"],
                recommendation,
                str(output),
                {
                    "data_type": "LFQ",
                    "quantification": "directlfq",
                    "upstream_engine": "quantms",
                    "threads": 24,
                },
            )
        )
    assert not output.exists()


def test_abstention_contract_is_unambiguous(
    lfq_inputs: dict[str, Path],
    tmp_path: Path,
) -> None:
    """An abstention cannot carry candidates, evidence, or elevated confidence."""
    recommendation = _inspect(lfq_inputs)["policy_recommendation"]
    recommendation.update(
        {
            "configs": [],
            "evidence_refs": [],
            "confidence": "high",
            "abstain_reason": "Host declined to propose a candidate.",
        }
    )

    with pytest.raises(ValueError, match="Abstention requires empty configs"):
        SERVICE.evaluate_recommendation(
            EvaluationRequest(
                str(lfq_inputs["matrix"]),
                str(lfq_inputs["sdrf"]),
                ["A", "B"],
                recommendation,
                str(tmp_path / "invalid-abstention"),
                {
                    "data_type": "LFQ",
                    "quantification": "directlfq",
                    "upstream_engine": "quantms",
                    "threads": 24,
                },
            )
        )


def test_score_a_counts_truth_only_within_the_tested_universe() -> None:
    """Undetected truth identifiers are not false negatives for a tested matrix."""
    de_table = pd.DataFrame(
        {
            "protein": ["P1", "P2"],
            "pvalue": [0.001, 0.8],
            "significance": ["UP", "NOT_DE"],
            "log2FC": [2.0, 0.0],
        }
    )
    protein_matrix = pd.DataFrame(
        {
            "protein": ["P1", "P2"],
            "A1": [100.0, 120.0],
            "A2": [101.0, 121.0],
            "B1": [400.0, 120.0],
            "B2": [404.0, 121.0],
        }
    )

    result = evaluate(
        CandidateConfig(name="tested-universe"),
        de_table,
        protein_matrix,
        {"A1": "A", "A2": "A", "B1": "B", "B2": "B"},
        ({"P1", "NOT_IN_MATRIX"}, "UP"),
    )

    assert result.truth_metrics.tp == 1
    assert result.truth_metrics.fn == 0


def test_ground_truth_requires_expected_direction(
    lfq_inputs: dict[str, Path],
    tmp_path: Path,
) -> None:
    """A labelled benchmark must declare the direction of its truth set."""
    recommendation = _inspect(lfq_inputs)["policy_recommendation"]
    recommendation["configs"] = recommendation["configs"][:1]

    with pytest.raises(ValueError, match="expected_direction is required"):
        SERVICE.evaluate_recommendation(
            EvaluationRequest(
                str(lfq_inputs["matrix"]),
                str(lfq_inputs["sdrf"]),
                ["A", "B"],
                recommendation,
                str(tmp_path / "missing-direction"),
                {
                    "ground_truth": str(lfq_inputs["truth"]),
                    "data_type": "LFQ",
                    "quantification": "directlfq",
                    "upstream_engine": "quantms",
                    "threads": 24,
                },
            )
        )


def test_expected_direction_requires_ground_truth(
    lfq_inputs: dict[str, Path],
    tmp_path: Path,
) -> None:
    """An unlabelled evaluation cannot declare a truth direction."""
    recommendation = _inspect(lfq_inputs)["policy_recommendation"]
    recommendation["configs"] = recommendation["configs"][:1]
    output_dir = tmp_path / "direction-without-truth"

    with pytest.raises(ValueError, match="must be null when ground_truth is null"):
        SERVICE.evaluate_recommendation(
            EvaluationRequest(
                str(lfq_inputs["matrix"]),
                str(lfq_inputs["sdrf"]),
                ["A", "B"],
                recommendation,
                str(output_dir),
                {
                    "ground_truth": None,
                    "expected_direction": "UP",
                    "data_type": "LFQ",
                    "quantification": "directlfq",
                    "upstream_engine": "quantms",
                    "threads": 24,
                },
            )
        )

    assert not output_dir.exists()


def test_inspection_abstains_without_two_conditions(
    lfq_inputs: dict[str, Path],
) -> None:
    """A one-condition SDRF cannot produce a DEA recommendation."""
    pd.DataFrame(
        {
            "source name": ["A1", "A2", "A3", "B1", "B2", "B3"],
            "factor value[condition]": ["A"] * 6,
        }
    ).to_csv(lfq_inputs["sdrf"], sep="\t", index=False)

    inspected = _inspect(lfq_inputs)

    assert inspected["policy_recommendation"]["configs"] == []
    assert inspected["policy_recommendation"]["abstain_reason"]
    diagnostics = next(
        block["content"]
        for block in inspected["context"]["blocks"]
        if block["id"] == "diagnostics"
    )
    assert "INSUFFICIENT_CONDITIONS" in {
        diagnostic["code"] for diagnostic in diagnostics
    }


def test_plugin_manifests_share_one_skill_and_knowledge_tree() -> None:
    """Codex and Claude adapters package the same workflow and evidence."""
    codex = json.loads(CODEX_MANIFEST.read_text(encoding="utf-8"))
    claude = json.loads(CLAUDE_MANIFEST.read_text(encoding="utf-8"))

    assert codex["name"] == claude["name"] == "mokume"
    assert codex["version"] == claude["version"]
    assert codex["skills"] == "./skills/"
    assert (PLUGIN / "skills" / "analyze-proteomics" / "SKILL.md").is_file()
    assert KNOWLEDGE.is_file()

    codex_marketplace = json.loads(
        (REPOSITORY / ".agents" / "plugins" / "marketplace.json").read_text(
            encoding="utf-8"
        )
    )
    claude_marketplace = json.loads(
        (REPOSITORY / ".claude-plugin" / "marketplace.json").read_text(encoding="utf-8")
    )
    assert codex_marketplace["plugins"][0]["name"] == "mokume"
    claude_entry = claude_marketplace["plugins"][0]
    assert claude_entry["name"] == "mokume"
    assert claude_entry["source"] == "./plugins/mokume"


def test_plugin_registers_the_local_mcp_tools(monkeypatch: pytest.MonkeyPatch) -> None:
    """Both host adapters start the same two public local tools."""
    pytest.importorskip("mcp")
    create_server = getattr(
        importlib.import_module("mokume.agentic.mcp_server"),
        "create_server",
    )

    codex_manifest = json.loads(CODEX_MANIFEST.read_text(encoding="utf-8"))
    server_config = codex_manifest["mcpServers"]["mokume"]
    assert server_config["command"] == "mokume"
    assert server_config["args"] == [
        "mcp",
        "serve",
        "--knowledge",
        "./knowledge/knowledge.yaml",
    ]
    assert server_config["cwd"] == "."

    claude = json.loads(CLAUDE_MANIFEST.read_text(encoding="utf-8"))
    claude_server = claude["mcpServers"]["mokume"]
    assert claude_server == {
        "command": "mokume",
        "args": [
            "mcp",
            "serve",
            "--knowledge",
            "${CLAUDE_PLUGIN_ROOT}/knowledge/knowledge.yaml",
        ],
    }

    monkeypatch.chdir(PLUGIN)
    server = create_server(server_config["args"][-1])
    tools = asyncio.run(server.list_tools())
    assert [tool.name for tool in tools] == [
        "inspect_dataset",
        "evaluate_recommendation",
    ]
    evaluation_schema = tools[1].inputSchema
    assert set(evaluation_schema["properties"]) == {
        "protein_matrix",
        "sdrf",
        "contrast",
        "recommendation",
        "options",
    }
    assert set(evaluation_schema["required"]) == set(evaluation_schema["properties"])
