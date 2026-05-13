"""Configuration proposal via LLM or rule-based fallback."""

import json

from mokume.agentic.config import AgenticConfig
from mokume.agentic.llm_client import (
    PROPOSAL_TOOL,
    call_with_tools,
)
from mokume.agentic.profiler import DataProfile
from mokume.agentic.rules import load_prompts, rule_propose, _load_heuristics
from mokume.agentic.state import CandidateConfig
from mokume.core.logger import get_logger

logger = get_logger("mokume.agentic.proposer")


def _items_to_configs(items: list[dict]) -> list[CandidateConfig]:
    """Convert parsed JSON items to CandidateConfig list."""
    configs = []
    for item in items:
        configs.append(
            CandidateConfig(
                name=item.get("name", "llm_config"),
                de_method=item.get("de_method", "deqms"),
                fdr_method=item.get("fdr_method", "bh"),
                normalization=item.get("normalization", "none"),
                imputation=item.get("imputation", "none"),
                log2fc_threshold=float(item.get("log2fc_threshold", 0.5)),
                ensemble=item.get("ensemble", "none"),
                ensemble_k=int(item.get("ensemble_k", 2)),
                reasoning=item.get("reasoning", ""),
                expected_outcome=item.get("expected_outcome", ""),
            )
        )
    return configs


def _build_heuristic_text(profile: DataProfile) -> str:
    """Format heuristics for the data type as text for the LLM."""
    heuristics = _load_heuristics()
    dt = profile.data_type
    priors = heuristics.get("data_type_priors", {}).get(dt, {})
    setting_specific = heuristics.get("setting_specific", {})
    rules = heuristics.get("condition_rules", [])
    return json.dumps(
        {"priors": priors, "setting_specific": setting_specific, "rules": rules},
        indent=2,
    )


def llm_propose(
    profile: DataProfile,
    config: AgenticConfig,
) -> list[CandidateConfig]:
    """Propose configs using LLM with tool calls."""
    prompts = load_prompts()
    heuristic_text = _build_heuristic_text(profile)

    system_msg = prompts["proposal_system"].format(
        relevant_heuristics=heuristic_text,
    )
    user_msg = prompts["proposal_user"].format(
        data_profile_json=json.dumps(profile.to_dict(), indent=2),
    )

    result = call_with_tools(system_msg, user_msg, [PROPOSAL_TOOL], config)
    items = result.get("configs", [])
    if not items:
        raise ValueError("LLM returned no configs")
    return _items_to_configs(items)


def propose_configs(
    profile: DataProfile,
    config: AgenticConfig,
) -> list[CandidateConfig]:
    """Propose configurations (LLM with rule-based fallback)."""
    if config.use_llm:
        try:
            configs = llm_propose(profile, config)
            logger.info("LLM proposed %d configs", len(configs))
            return configs
        except Exception as exc:  # noqa: BLE001
            logger.warning("LLM unavailable (%s), using rule-based", exc)

    configs = rule_propose(profile)
    return configs
