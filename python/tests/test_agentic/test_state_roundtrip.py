"""Tests for AgenticState + children to_dict / from_dict round-trips (--resume)."""

from mokume.agentic.state import (
    AgenticState,
    AuditEntry,
    CandidateConfig,
    EvaluationResult,
    ReflectionResult,
    RoundResult,
)


def _mk_state() -> AgenticState:
    c1 = CandidateConfig(
        name="c1",
        de_method="limma",
        fdr_method="bh",
        normalization="none",
        imputation="none",
        log2fc_threshold=0.5,
        ensemble="none",
        ensemble_k=2,
        reasoning="baseline",
    )
    c2 = CandidateConfig(
        name="c2",
        de_method="ensemble",
        fdr_method="ihw",
        normalization="median",
        imputation="seqknn",
        log2fc_threshold=1.0,
        ensemble="limma,deqms,proda",
        ensemble_k=3,
    )
    r1 = EvaluationResult(
        config_name="c1",
        config=c1.to_dict(),
        tp=4,
        fp=1,
        auc=0.92,
        n_de_up=5,
        n_de_down=0,
        score=0.87,
    )
    r2 = EvaluationResult(
        config_name="c2",
        config=c2.to_dict(),
        score=-1.0,
    )
    reflection = ReflectionResult(
        converged=False,
        next_configs=[c2],
        analysis="explore ensemble",
        adjustments=["switch to ensemble"],
    )
    rnd = RoundResult(
        round_num=1,
        configs=[c1, c2],
        results=[r1, r2],
        best_config_name="c1",
        reflection=reflection,
    )
    audit = [AuditEntry(step="propose", round_num=1, data=["c1", "c2"])]
    return AgenticState(
        rounds=[rnd],
        total_experiments=2,
        converged=False,
        best_config=c1,
        best_score=0.87,
        audit_trail=audit,
    )


def test_candidate_config_roundtrip():
    """CandidateConfig.to_dict -> from_dict is lossless for schema fields."""
    cfg = CandidateConfig(
        name="c",
        de_method="rots",
        ensemble="limma,rots",
        ensemble_k=2,
        log2fc_threshold=1.0,
    )
    back = CandidateConfig.from_dict(cfg.to_dict())
    assert back.signature() == cfg.signature()
    assert back.name == cfg.name


def test_evaluation_result_roundtrip():
    """EvaluationResult survives to_dict -> from_dict."""
    er = EvaluationResult(
        config_name="x",
        config={"name": "x"},
        tp=3,
        fp=0,
        auc=0.99,
        n_de_up=5,
        n_de_down=1,
        score=0.77,
    )
    back = EvaluationResult.from_dict(er.to_dict())
    assert back.config_name == er.config_name
    assert back.tp == 3 and back.fp == 0
    assert back.score == 0.77


def test_reflection_result_roundtrip():
    """ReflectionResult preserves next_configs and adjustments."""
    ref = ReflectionResult(
        converged=True,
        next_configs=[CandidateConfig(name="next")],
        analysis="done",
        adjustments=["none"],
    )
    back = ReflectionResult.from_dict(ref.to_dict())
    assert back.converged is True
    assert back.analysis == "done"
    assert [c.name for c in back.next_configs] == ["next"]


def test_full_state_roundtrip():
    """AgenticState.to_dict -> from_dict preserves the whole optimization state."""
    state = _mk_state()
    restored = AgenticState.from_dict(state.to_dict())
    assert restored.total_experiments == state.total_experiments
    assert restored.best_score == state.best_score
    assert restored.best_config is not None
    assert restored.best_config.signature() == state.best_config.signature()
    assert len(restored.rounds) == 1
    rnd = restored.rounds[0]
    assert rnd.round_num == 1
    assert rnd.best_config_name == "c1"
    assert len(rnd.configs) == 2
    assert len(rnd.results) == 2
    assert rnd.results[0].score == 0.87
    assert rnd.reflection is not None
    assert rnd.reflection.analysis == "explore ensemble"
    assert len(restored.audit_trail) == 1
    assert restored.audit_trail[0].step == "propose"
