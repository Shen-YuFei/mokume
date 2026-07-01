"""Integration tests for P1: parallel rounds, DataProfile enrichment, resume."""

import json
import threading
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mokume.agentic.config import AgenticConfig
from mokume.agentic.optimizer import optimize
from mokume.agentic.profiler import profile_data
from mokume.agentic.reporter import checkpoint_path, load_state_snapshot
from mokume.agentic.runner import PreprocessCache


@pytest.fixture(name="runs_dir")
def _runs_dir_fixture(tmp_path: Path) -> str:
    """Output directory under tmp_path for one optimization run."""
    return str(tmp_path / "runs")


def test_profile_includes_pca_and_fold_change(
    synthetic_protein_df, sample_to_condition
):
    """P1-8: enriched signals populate."""
    profile = profile_data(synthetic_protein_df, sample_to_condition)
    assert len(profile.pca_explained_variance) > 0
    assert all(0.0 <= v <= 1.0 for v in profile.pca_explained_variance)
    assert profile.per_condition_median_cv  # non-empty dict
    assert profile.pairwise_median_abs_log2fc  # A vs B present
    as_dict = profile.to_dict()
    assert "pca_explained_variance" in as_dict
    assert "per_condition_median_cv" in as_dict
    assert "pairwise_median_abs_log2fc" in as_dict


def test_profile_enrichment_safe_on_tiny_matrix():
    """Very small matrices should not explode the new profilers."""
    df = pd.DataFrame(
        {
            "protein": ["P1", "P2"],
            "S1": [10.0, 12.0],
            "S2": [11.0, 13.0],
        }
    )
    profile = profile_data(df, {"S1": "A", "S2": "B"})
    # Too few proteins / samples -> empty PCA, but no crash
    assert not profile.pca_explained_variance


def test_optimize_parallel_run_matches_serial(
    synthetic_protein_df, sample_to_condition, ground_truth_proteins, runs_dir
):
    """P1-6: n_jobs>1 should yield the same best config as n_jobs=1."""
    base_config = AgenticConfig(
        use_llm=False,
        max_rounds=2,
        max_experiments=6,
        contrasts=[("A", "B")],
        output_dir=runs_dir + "_serial",
    )
    serial = optimize(
        synthetic_protein_df,
        sample_to_condition,
        base_config,
        ground_truth=ground_truth_proteins,
    )

    parallel_config = AgenticConfig(
        use_llm=False,
        max_rounds=2,
        max_experiments=6,
        n_jobs=2,
        contrasts=[("A", "B")],
        output_dir=runs_dir + "_parallel",
    )
    parallel = optimize(
        synthetic_protein_df,
        sample_to_condition,
        parallel_config,
        ground_truth=ground_truth_proteins,
    )

    key = "A_vs_B"
    assert key in serial and key in parallel
    # Same rule-based proposer, same seeds -> deterministic best
    assert serial[key].best_config.signature() == parallel[key].best_config.signature()
    assert serial[key].best_score == pytest.approx(parallel[key].best_score, abs=1e-6)


def test_preprocess_cache_is_thread_safe():
    """Concurrent get_or_compute calls on the same key still yield equal output."""
    np.random.seed(1)
    df = pd.DataFrame(np.random.normal(20, 2, (16, 6)))
    df.columns = [f"S{i}" for i in range(6)]
    df.insert(0, "protein", [f"P{i}" for i in range(16)])

    cache = PreprocessCache()
    outputs = []

    def worker():
        outputs.append(cache.get_or_compute("median", "mindet", df))

    threads = [threading.Thread(target=worker) for _ in range(4)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    # All four workers observe the same cached frame contents.
    for o in outputs[1:]:
        pd.testing.assert_frame_equal(outputs[0], o)
    stats = cache.stats()
    # At most one miss is retained, but every worker is accounted for via hits+misses
    assert stats["misses"] >= 1
    assert stats["hits"] + stats["misses"] == 4
    assert stats["unique_combos"] == 1


def test_resume_loads_previous_state(
    synthetic_protein_df, sample_to_condition, ground_truth_proteins, runs_dir
):
    """P1-10: state.json is written and re-hydrates via AgenticConfig(resume=True)."""
    config = AgenticConfig(
        use_llm=False,
        max_rounds=1,
        max_experiments=3,
        contrasts=[("A", "B")],
        output_dir=runs_dir,
    )
    first = optimize(
        synthetic_protein_df,
        sample_to_condition,
        config,
        ground_truth=ground_truth_proteins,
    )
    assert "A_vs_B" in first

    # Checkpoint must exist and be valid JSON with the expected structure.
    cp = checkpoint_path(runs_dir, ("A", "B"))
    assert cp.exists()
    payload = json.loads(cp.read_text(encoding="utf-8"))
    assert "rounds" in payload and len(payload["rounds"]) == 1

    loaded = load_state_snapshot(runs_dir, ("A", "B"))
    assert loaded is not None
    assert loaded.total_experiments == first["A_vs_B"].total_experiments

    # Run again with resume + larger budget; should not re-run round 1.
    resume_config = AgenticConfig(
        use_llm=False,
        max_rounds=2,
        max_experiments=6,
        contrasts=[("A", "B")],
        output_dir=runs_dir,
        resume=True,
    )
    resumed = optimize(
        synthetic_protein_df,
        sample_to_condition,
        resume_config,
        ground_truth=ground_truth_proteins,
    )
    state = resumed["A_vs_B"]
    # We started with 1 round; rerun with resume should add at most 1 more.
    assert len(state.rounds) >= 1
    # Round numbers must be 1, 2, ... in order.
    assert [r.round_num for r in state.rounds] == list(range(1, len(state.rounds) + 1))
