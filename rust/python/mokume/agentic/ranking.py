"""Candidate ranking for Mokume Plugin evaluation rounds."""

from __future__ import annotations

import math
from typing import Any

import pandas as pd


_BENCHMARK_RANK_METRICS = (
    ("pauc001", "rank_pauc001"),
    ("pauc005", "rank_pauc005"),
    ("pauc", "rank_pauc01"),
    ("nmcc", "rank_nmcc"),
    ("gmean", "rank_gmean"),
)


def build_ranking_payload(
    rows: list[dict[str, Any]],
    has_truth: bool,
    cache: dict[str, int],
) -> dict[str, Any]:
    """Apply the ground-truth mean-rank or exploratory contract."""
    ranking = _rank_ground_truth_rows(rows) if has_truth else []
    status = "ranked" if ranking else "ground_truth_unscored"
    if not has_truth:
        status = "exploratory_unranked"
    return {
        "status": status,
        "ranking_objective": "benchmark_mean_rank" if has_truth else None,
        "results": rows,
        "ranking": ranking,
        "best_config": ranking[0]["config_name"] if ranking else None,
        "cache": cache,
    }


def _rank_ground_truth_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Rank complete candidates per metric and average their five ranks."""
    rank_fields = [rank_field for _, rank_field in _BENCHMARK_RANK_METRICS]
    for row in rows:
        row["benchmark_mean_rank"] = None
        for rank_field in rank_fields:
            row[rank_field] = None

    eligible = [row for row in rows if _has_complete_rank_metrics(row)]
    for metric, rank_field in _BENCHMARK_RANK_METRICS:
        ranks = pd.Series([float(row[metric]) for row in eligible], dtype=float).rank(
            method="average", ascending=False
        )
        for row, rank in zip(eligible, ranks, strict=True):
            row[rank_field] = float(rank)

    for row in eligible:
        row["benchmark_mean_rank"] = sum(
            float(row[rank_field]) for rank_field in rank_fields
        ) / len(rank_fields)
    return sorted(
        eligible,
        key=lambda row: (
            float(row["benchmark_mean_rank"]),
            -float(row["score_a"]),
            str(row["config_name"]),
        ),
    )


def _has_complete_rank_metrics(row: dict[str, Any]) -> bool:
    """Return whether a row has every finite metric required for ranking."""
    values = [row.get(metric) for metric, _ in _BENCHMARK_RANK_METRICS]
    values.append(row.get("score_a"))
    return all(
        isinstance(value, (int, float)) and math.isfinite(float(value))
        for value in values
    )
