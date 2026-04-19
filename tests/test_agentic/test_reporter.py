"""Tests for mokume.agentic.reporter."""

import json
import tempfile
from pathlib import Path

from mokume.agentic.config import AgenticConfig
from mokume.agentic.profiler import DataProfile
from mokume.agentic.reporter import save_outputs
from mokume.agentic.state import (
    AgenticState,
    CandidateConfig,
    EvaluationResult,
    RoundResult,
)


def test_save_outputs_creates_files():
    """save_outputs writes all expected files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        profile = DataProfile(
            n_proteins=100,
            n_samples=6,
            n_conditions=2,
            samples_per_condition={"A": 3, "B": 3},
            missing_rate=0.1,
            missing_per_sample={},
            median_cv=0.12,
            has_peptide_counts=False,
            data_type="LFQ",
            batch_fields=[],
            intensity_range=(10.0, 25.0),
            is_log_transformed=True,
        )
        best_cfg = CandidateConfig(name="best", de_method="deqms")
        state = AgenticState(
            best_config=best_cfg,
            best_score=0.95,
            total_experiments=5,
            converged=True,
            rounds=[
                RoundResult(
                    round_num=1,
                    configs=[best_cfg],
                    results=[
                        EvaluationResult(
                            config_name="best",
                            config=best_cfg.to_dict(),
                            tp=47,
                            fp=1,
                            auc=0.996,
                            n_de_up=100,
                            n_de_down=20,
                            score=0.95,
                        )
                    ],
                    best_config_name="best",
                )
            ],
        )
        config = AgenticConfig(output_dir=tmpdir, use_llm=False)
        out = save_outputs(profile, state, config)

        assert (Path(out) / "summary.md").exists()
        assert (Path(out) / "best_config.yaml").exists()
        assert (Path(out) / "all_results.tsv").exists()
        assert (Path(out) / "data_profile.json").exists()
        assert (Path(out) / "audit_trail.json").exists()

        # Verify JSON is valid
        profile_json = json.loads((Path(out) / "data_profile.json").read_text())
        assert profile_json["n_proteins"] == 100
