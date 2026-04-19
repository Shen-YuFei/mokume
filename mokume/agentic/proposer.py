"""Configuration proposal via LLM or rule-based fallback."""

import json

from mokume.agentic.config import AgenticConfig
from mokume.agentic.profiler import DataProfile
from mokume.agentic.rules import load_prompts, rule_propose, _load_heuristics
from mokume.agentic.state import CandidateConfig
from mokume.core.logger import get_logger

logger = get_logger("mokume.agentic.proposer")


class LLMUnavailableError(Exception):
    """Raised when LLM provider is not available."""


def _get_llm(config: AgenticConfig):
    """Get a ChatOpenAI instance (works with any OpenAI-compatible API)."""
    try:
        from langchain_openai import ChatOpenAI  # pylint: disable=import-outside-toplevel
    except ImportError as exc:
        raise LLMUnavailableError(
            "LLM requires langchain-openai: pip install langchain-openai"
        ) from exc

    kwargs = {
        "model": config.llm_model,
        "temperature": config.llm_temperature,
        "max_tokens": 4096,
    }
    if config.llm_base_url:
        kwargs["base_url"] = config.llm_base_url
    if config.llm_api_key:
        kwargs["api_key"] = config.llm_api_key
    return ChatOpenAI(**kwargs)


def _parse_configs_from_json(raw: str) -> list[CandidateConfig]:
    """Parse LLM JSON response into CandidateConfig list."""
    # Extract JSON array from response
    text = raw.strip()
    start = text.find("[")
    end = text.rfind("]") + 1
    if start < 0 or end <= start:
        raise ValueError("No JSON array found in LLM response")

    items = json.loads(text[start:end])
    configs = []
    for item in items:
        configs.append(CandidateConfig(
            name=item.get("name", "llm_config"),
            de_method=item.get("de_method", "deqms"),
            fdr_method=item.get("fdr_method", "bh"),
            imputation=item.get("imputation", "none"),
            log2fc_threshold=float(item.get("log2fc_threshold", 0.5)),
            reasoning=item.get("reasoning", ""),
            expected_outcome=item.get("expected_outcome", ""),
        ))
    return configs


def _build_proposal_prompt(profile: DataProfile) -> str:
    """Build the proposal prompt from template + data."""
    prompts = load_prompts()
    heuristics = _load_heuristics()

    # Format heuristics for the data type
    dt = profile.data_type
    priors = heuristics.get("data_type_priors", {}).get(dt, {})
    rules = heuristics.get("condition_rules", [])
    heuristic_text = json.dumps(
        {"priors": priors, "rules": rules}, indent=2,
    )

    return prompts["proposal"].format(
        data_profile_json=json.dumps(profile.to_dict(), indent=2),
        relevant_heuristics=heuristic_text,
    )


def llm_propose(
    profile: DataProfile,
    config: AgenticConfig,
) -> list[CandidateConfig]:
    """Propose configs using LLM."""
    llm = _get_llm(config)
    prompt = _build_proposal_prompt(profile)
    response = llm.invoke(prompt)
    text = response.content if hasattr(response, "content") else str(response)
    return _parse_configs_from_json(text)


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
        except (LLMUnavailableError, ValueError, ConnectionError) as exc:
            logger.warning("LLM unavailable (%s), using rule-based", exc)

    configs = rule_propose(profile)
    return configs
