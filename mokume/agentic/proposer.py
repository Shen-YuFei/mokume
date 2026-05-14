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


def _dedup_configs(
    configs: list[CandidateConfig],
    seen: set[str] | None,
) -> list[CandidateConfig]:
    """Drop configs whose signature is already seen, preserving order."""
    if not seen:
        return configs
    deduped: list[CandidateConfig] = []
    for cfg in configs:
        sig = cfg.signature()
        if sig in seen:
            logger.debug("Skipping duplicate config %s (%s)", cfg.name, sig)
            continue
        deduped.append(cfg)
    return deduped


def _format_seen_block(seen: set[str] | None) -> str:
    """Format previously-tried signatures for the LLM prompt."""
    if not seen:
        return "(none -- this is the first round)"
    return "\n".join(f"- {sig}" for sig in sorted(seen))


def llm_propose(
    profile: DataProfile,
    config: AgenticConfig,
    seen: set[str] | None = None,
) -> list[CandidateConfig]:
    """Propose configs using LLM with tool calls."""
    prompts = load_prompts()
    heuristic_text = _build_heuristic_text(profile)

    system_msg = prompts["proposal_system"].format(
        relevant_heuristics=heuristic_text,
    )
    user_msg = prompts["proposal_user"].format(
        data_profile_json=json.dumps(profile.to_dict(), indent=2),
        previous_signatures=_format_seen_block(seen),
    )

    result = call_with_tools(system_msg, user_msg, [PROPOSAL_TOOL], config)
    items = result.get("configs", [])
    if not items:
        raise ValueError("LLM returned no configs")
    return _items_to_configs(items)


def propose_configs(
    profile: DataProfile,
    config: AgenticConfig,
    seen: set[str] | None = None,
) -> list[CandidateConfig]:
    """Propose configurations (LLM with rule-based fallback).

    ``seen`` is a set of CandidateConfig.signature() strings already tried;
    matching candidates are filtered out and (when LLM is on) also surfaced
    in the prompt so the model proactively avoids them.
    """
    if config.use_llm:
        try:
            configs = llm_propose(profile, config, seen=seen)
            configs = _dedup_configs(configs, seen)
            logger.info("LLM proposed %d configs (after dedup)", len(configs))
            if configs:
                return configs
            logger.warning("LLM returned only duplicates, falling back to rules")
        except Exception as exc:  # noqa: BLE001
            logger.warning("LLM unavailable (%s), using rule-based", exc)

    configs = rule_propose(profile)
    return _dedup_configs(configs, seen)
