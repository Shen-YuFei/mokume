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


def _format_results_table(rounds: list[RoundResult]) -> str:
    """Format all results as a text table for LLM."""
    lines = ["Round | Config | Score | TP | FP | AUC | #UP | #DOWN"]
    lines.append("-" * 60)
    for rnd in rounds:
        for res in rnd.results:
            lines.append(
                f"{rnd.round_num} | {res.config_name} | "
                f"{res.score:.3f} | {res.tp} | {res.fp} | "
                f"{res.auc} | {res.n_de_up} | {res.n_de_down}"
            )
    return "\n".join(lines)


def _parse_reflection_data(data: dict) -> ReflectionResult:
    """Convert parsed tool-call dict to ReflectionResult."""
    next_configs = _items_to_configs(data.get("next_configs", []))
    return ReflectionResult(
        converged=data.get("convergence", False),
        next_configs=next_configs,
        analysis=data.get("analysis", ""),
        adjustments=data.get("adjustments", []),
    )


def _rule_reflect(rounds: list[RoundResult]) -> ReflectionResult:
    """Rule-based reflection: converge if best unchanged 2 rounds."""
    if len(rounds) < 2:
        return ReflectionResult(
            converged=False,
            analysis="Only 1 round completed, need more data",
        )

    best_prev = rounds[-2].best_config_name
    best_curr = rounds[-1].best_config_name

    if best_prev == best_curr:
        return ReflectionResult(
            converged=True,
            analysis=f"Best config '{best_curr}' stable for 2 rounds",
        )

    return ReflectionResult(
        converged=False,
        analysis=f"Best changed: {best_prev} → {best_curr}",
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

    return _rule_reflect(rounds)
