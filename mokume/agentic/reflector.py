"""Result analysis and next-round proposal via LLM or rules."""

import json

from mokume.agentic.config import AgenticConfig
from mokume.agentic.llm_client import (
    REFLECTION_TOOL,
    LLMUnavailableError,
    call_with_tools,
)
from mokume.agentic.profiler import DataProfile
from mokume.agentic.proposer import _build_heuristic_text, _items_to_configs
from mokume.agentic.rules import load_prompts
from mokume.agentic.state import (
    ReflectionResult,
    RoundResult,
)
from mokume.core.logger import get_logger

logger = get_logger("mokume.agentic.reflector")


def _collect_all_results(rounds: list[RoundResult]) -> list[tuple[int, object]]:
    """Flatten (round_num, result) pairs across every round."""
    return [(rnd.round_num, res) for rnd in rounds for res in rnd.results]


def _format_top3_section(successful: list[tuple[int, object]]) -> list[str]:
    """Format the Top-3 globally section with a * marker on the leader."""
    if not successful:
        return []
    lines = ["Top-3 globally (* = current best):"]
    top3 = successful[:3]
    best_score = top3[0][1].score
    for rank, (rnum, res) in enumerate(top3, start=1):
        marker = "*" if rank == 1 else " "
        delta = res.score - best_score
        lines.append(
            f"  {marker} #{rank} R{rnum} {res.config_name} "
            f"score={res.score:.3f} (delta={delta:+.3f}) "
            f"TP={res.tp} FP={res.fp} AUC={res.auc} "
            f"UP={res.n_de_up} DOWN={res.n_de_down}"
        )
    lines.append("")
    return lines


def _format_failed_section(failed: list[tuple[int, object]]) -> list[str]:
    """Format the failed-configurations section (capped at 5 entries)."""
    if not failed:
        return []
    lines = [f"Failed configurations ({len(failed)}, avoid in next round):"]
    for rnum, res in failed[:5]:
        lines.append(f"  X R{rnum} {res.config_name}")
    if len(failed) > 5:
        lines.append(f"  ... and {len(failed) - 5} more")
    lines.append("")
    return lines


def _format_trajectory_section(rounds: list[RoundResult]) -> list[str]:
    """Format the per-round best trajectory section."""
    lines = ["Per-round best (convergence trajectory):"]
    for rnd in rounds:
        if rnd.results and rnd.results[0].score > -1.0:
            top = rnd.results[0]
            lines.append(f"  R{rnd.round_num}: {top.config_name} score={top.score:.3f}")
        else:
            lines.append(f"  R{rnd.round_num}: (all failed)")
    lines.append("")
    return lines


def _format_score_range_section(successful: list[tuple[int, object]]) -> list[str]:
    """Format the score-range / gap summary section."""
    if not successful:
        return []
    score_gap = successful[0][1].score - successful[-1][1].score
    return [
        (
            f"Score range across successful runs: "
            f"{successful[-1][1].score:.3f} - {successful[0][1].score:.3f} "
            f"(gap={score_gap:.3f})"
        ),
        "",
    ]


def _format_leaderboard_section(successful: list[tuple[int, object]]) -> list[str]:
    """Format the full leaderboard section sorted by score descending."""
    lines = [
        "Full leaderboard (score desc):",
        "Rank | Round | Config | Score | TP | FP | AUC | #UP | #DOWN",
        "-" * 70,
    ]
    for rank, (rnum, res) in enumerate(successful, start=1):
        lines.append(
            f"{rank:>4} | R{rnum} | {res.config_name} | "
            f"{res.score:.3f} | {res.tp} | {res.fp} | "
            f"{res.auc} | {res.n_de_up} | {res.n_de_down}"
        )
    return lines


def _format_results_table(rounds: list[RoundResult]) -> str:
    """Format results for the LLM, sorted by score with Top-3 marked.

    Sections produced:
      1. Top-3 configurations globally (* marker)
      2. Failed configurations (score == -1) -- useful for the LLM to avoid
      3. Per-round best (trajectory) -- convergence signal
      4. Score range summary
      5. Full leaderboard
    """
    flat = _collect_all_results(rounds)
    if not flat:
        return "No experiment results yet."

    successful = sorted(
        (p for p in flat if p[1].score > -1.0),
        key=lambda p: p[1].score,
        reverse=True,
    )
    failed = [p for p in flat if p[1].score <= -1.0]

    sections = [
        _format_top3_section(successful),
        _format_failed_section(failed),
        _format_trajectory_section(rounds),
        _format_score_range_section(successful),
        _format_leaderboard_section(successful),
    ]
    return "\n".join(line for section in sections for line in section)


def _parse_reflection_data(data: dict) -> ReflectionResult:
    """Convert parsed tool-call dict to ReflectionResult."""
    next_configs = _items_to_configs(data.get("next_configs", []))
    return ReflectionResult(
        converged=data.get("convergence", False),
        next_configs=next_configs,
        analysis=data.get("analysis", ""),
        adjustments=data.get("adjustments", []),
    )


def _best_successful_score(rnd: RoundResult) -> float | None:
    """Return the top non-failed score of a round, or None if all failed."""
    for res in rnd.results:
        if res.score > -1.0:
            return res.score
    return None


def _rule_reflect(
    rounds: list[RoundResult],
    score_gap: float = 0.02,
) -> ReflectionResult:
    """Rule-based reflection: converge on either name stability or score gap."""
    if len(rounds) < 2:
        return ReflectionResult(
            converged=False,
            analysis="Only 1 round completed, need more data",
        )

    prev, curr = rounds[-2], rounds[-1]
    best_prev, best_curr = prev.best_config_name, curr.best_config_name
    score_prev = _best_successful_score(prev)
    score_curr = _best_successful_score(curr)

    # No successful results at all -> keep searching
    if score_curr is None:
        return ReflectionResult(
            converged=False,
            analysis="Current round has no successful results",
        )

    if best_prev == best_curr:
        delta = (score_curr - score_prev) if score_prev is not None else 0.0
        if delta <= score_gap:
            return ReflectionResult(
                converged=True,
                analysis=(
                    f"Best config '{best_curr}' stable across 2 rounds, "
                    f"score delta={delta:+.3f} <= {score_gap:.3f}"
                ),
            )
        return ReflectionResult(
            converged=False,
            analysis=(
                f"Best '{best_curr}' stable but score still improving "
                f"(delta={delta:+.3f} > {score_gap:.3f})"
            ),
        )

    if score_prev is not None:
        delta = score_curr - score_prev
        if abs(delta) <= score_gap:
            return ReflectionResult(
                converged=True,
                analysis=(
                    f"Best changed {best_prev} -> {best_curr} but score is flat "
                    f"(delta={delta:+.3f}, within {score_gap:.3f})"
                ),
            )

    return ReflectionResult(
        converged=False,
        analysis=f"Best changed: {best_prev} -> {best_curr}",
    )


def reflect(
    profile: DataProfile,
    rounds: list[RoundResult],
    config: AgenticConfig,
) -> ReflectionResult:
    """Analyze results and propose next actions."""
    if config.use_llm:
        try:
            prompts = load_prompts()
            heuristic_text = _build_heuristic_text(profile)

            system_msg = prompts["reflection_system"].format(
                relevant_heuristics=heuristic_text,
            )
            user_msg = prompts["reflection_user"].format(
                data_profile_json=json.dumps(profile.to_dict(), indent=2),
                results_table=_format_results_table(rounds),
            )

            data = call_with_tools(
                system_msg,
                user_msg,
                [REFLECTION_TOOL],
                config,
            )
            result = _parse_reflection_data(data)
            logger.info("LLM reflection: converged=%s", result.converged)
            return result
        except (LLMUnavailableError, ValueError, ConnectionError) as exc:
            logger.warning("LLM reflection failed (%s), using rules", exc)

    return _rule_reflect(rounds, score_gap=config.convergence_score_gap)
