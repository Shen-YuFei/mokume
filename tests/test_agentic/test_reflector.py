"""Tests for mokume.agentic.reflector."""

from mokume.agentic.reflector import (
    _collect_all_results,
    _format_results_table,
    _rule_reflect,
)
from mokume.agentic.state import (
    CandidateConfig,
    EvaluationResult,
    RoundResult,
)


def _mk_result(name: str, score: float, tp: int | None = None) -> EvaluationResult:
    return EvaluationResult(
        config_name=name,
        config={"name": name},
        score=score,
        tp=tp,
        fp=0 if tp is not None else None,
        auc=score if tp is not None else None,
        n_de_up=5,
        n_de_down=0,
    )


def _mk_round(num: int, *results: EvaluationResult) -> RoundResult:
    return RoundResult(
        round_num=num,
        configs=[CandidateConfig(name=r.config_name) for r in results],
        results=list(results),
        best_config_name=results[0].config_name if results else "",
    )


def test_format_results_table_empty():
    """Empty rounds produce a graceful message."""
    out = _format_results_table([])
    assert "No experiment results" in out


def test_format_results_table_sorts_and_marks_top3():
    """Top-3 are starred and sorted by score, regardless of input order."""
    rounds = [
        _mk_round(1, _mk_result("low", 0.2, 1), _mk_result("high", 0.9, 5)),
        _mk_round(2, _mk_result("mid", 0.5, 3)),
    ]
    out = _format_results_table(rounds)
    assert "Top-3 globally" in out
    assert "*" in out  # best marker
    lo = out.index("low")
    hi = out.index("high")
    assert hi < lo  # 'high' appears before 'low' in leaderboard


def test_format_results_table_lists_failed():
    """Failed configurations (score == -1) show up in a dedicated block."""
    rounds = [_mk_round(1, _mk_result("ok", 0.8, 3), _mk_result("bad", -1.0))]
    out = _format_results_table(rounds)
    assert "Failed configurations" in out
    assert "bad" in out


def test_format_results_table_per_round_trajectory():
    """Per-round best section shows the trajectory across rounds."""
    rounds = [
        _mk_round(1, _mk_result("v1", 0.5, 2)),
        _mk_round(2, _mk_result("v2", 0.7, 3)),
    ]
    out = _format_results_table(rounds)
    assert "trajectory" in out.lower() or "Per-round best" in out
    assert "R1" in out and "R2" in out


def test_collect_all_results_flattens():
    """_collect_all_results flattens rounds into (round_num, result) pairs."""
    rounds = [
        _mk_round(1, _mk_result("a", 0.1), _mk_result("b", 0.2)),
        _mk_round(2, _mk_result("c", 0.3)),
    ]
    flat = _collect_all_results(rounds)
    assert len(flat) == 3
    assert flat[0][0] == 1 and flat[0][1].config_name == "a"
    assert flat[2][0] == 2 and flat[2][1].config_name == "c"


def test_rule_reflect_converges_when_best_stable_and_flat():
    """Same best name and delta-score <= threshold -> converged."""
    rounds = [
        _mk_round(1, _mk_result("winner", 0.80, 4)),
        _mk_round(2, _mk_result("winner", 0.81, 4)),
    ]
    out = _rule_reflect(rounds, score_gap=0.02)
    assert out.converged is True


def test_rule_reflect_keeps_searching_when_still_improving():
    """Same best name but delta-score > threshold -> keep searching."""
    rounds = [
        _mk_round(1, _mk_result("winner", 0.80, 4)),
        _mk_round(2, _mk_result("winner", 0.90, 4)),
    ]
    out = _rule_reflect(rounds, score_gap=0.02)
    assert out.converged is False


def test_rule_reflect_converges_when_best_changes_but_score_flat():
    """Best name changed but delta-score ~= 0 -> saturation, converge."""
    rounds = [
        _mk_round(1, _mk_result("a", 0.70, 3)),
        _mk_round(2, _mk_result("b", 0.71, 3)),
    ]
    out = _rule_reflect(rounds, score_gap=0.02)
    assert out.converged is True


def test_rule_reflect_keeps_searching_when_best_changes_and_improves():
    """Best changed and delta-score material -> keep searching."""
    rounds = [
        _mk_round(1, _mk_result("first", 0.60, 3)),
        _mk_round(2, _mk_result("second", 0.70, 4)),
    ]
    out = _rule_reflect(rounds, score_gap=0.02)
    assert out.converged is False


def test_rule_reflect_waits_after_single_round():
    """Single round means we cannot judge stability yet."""
    rounds = [_mk_round(1, _mk_result("only", 0.5, 2))]
    out = _rule_reflect(rounds)
    assert out.converged is False


def test_rule_reflect_keeps_searching_when_all_failed():
    """Latest round with no successful results -> keep searching."""
    rounds = [
        _mk_round(1, _mk_result("good", 0.6, 2)),
        _mk_round(2, _mk_result("bad", -1.0)),
    ]
    out = _rule_reflect(rounds)
    assert out.converged is False
