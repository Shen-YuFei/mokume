from __future__ import annotations

import pytest

from mokume.agentic.ranking import build_ranking_payload


def test_ground_truth_ranking_uses_five_metric_mean_rank() -> None:
    """The absolute Score A winner need not be the benchmark mean-rank winner."""
    rows = [
        {
            "config_name": "A",
            "score_a": 0.65,
            "pauc001": 0.90,
            "pauc005": 0.90,
            "pauc": 0.40,
            "nmcc": 0.40,
            "gmean": 0.40,
        },
        {
            "config_name": "B",
            "score_a": 0.64,
            "pauc001": 0.64,
            "pauc005": 0.64,
            "pauc": 0.64,
            "nmcc": 0.64,
            "gmean": 0.64,
        },
        {
            "config_name": "C",
            "score_a": 0.63,
            "pauc001": 0.63,
            "pauc005": 0.63,
            "pauc": 0.63,
            "nmcc": 0.63,
            "gmean": 0.63,
        },
    ]

    payload = build_ranking_payload(rows, True, {})

    assert payload["ranking_objective"] == "benchmark_mean_rank"
    assert payload["best_config"] == "B"
    assert [row["config_name"] for row in payload["ranking"]] == ["B", "A", "C"]
    assert rows[0]["benchmark_mean_rank"] == pytest.approx(2.2)
    assert rows[1]["benchmark_mean_rank"] == pytest.approx(1.4)
    assert rows[2]["benchmark_mean_rank"] == pytest.approx(2.4)
