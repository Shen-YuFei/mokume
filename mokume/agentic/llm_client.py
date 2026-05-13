"""Unified LLM client for agentic analysis via OpenAI-compatible API.

Supports DeepSeek V4 Pro thinking mode, Tool Calls (strict), and
automatic fallback (strict → JSON mode → raw text parsing).
"""

import json
from typing import Any

from mokume.agentic.config import AgenticConfig
from mokume.core.logger import get_logger

logger = get_logger("mokume.agentic.llm_client")


class LLMUnavailableError(Exception):
    """Raised when LLM provider is not available."""


# ---------------------------------------------------------------------------
# Tool schemas (strict mode)
# ---------------------------------------------------------------------------

_CONFIG_ITEM_SCHEMA: dict[str, Any] = {
    "type": "object",
    "properties": {
        "name": {"type": "string", "description": "Short descriptive name"},
        "de_method": {
            "type": "string",
            "enum": [
                "limma",
                "rots",
                "deqms",
                "proda",
                "msstats",
                "ensemble",
            ],
            "description": (
                "Differential expression method. Use 'ensemble' to combine "
                "multiple methods via top-k consensus (set ensemble + ensemble_k)."
            ),
        },
        "fdr_method": {
            "type": "string",
            "enum": ["bh", "ihw"],
            "description": "FDR correction method",
        },
        "normalization": {
            "type": "string",
            "enum": [
                "none",
                "median",
                "quantile",
                "mean",
                "rlr",
                "vsn",
                "loess",
                "mbqn",
            ],
            "description": "Normalization method",
        },
        "imputation": {
            "type": "string",
            "enum": [
                "none",
                "minprob",
                "mindet",
                "knn",
                "missforest",
                "seqknn",
                "qrilc",
                "mle",
                "mice",
                "nbavg",
                "gms",
                "bpca",
                "impseq",
                "impseqrob",
            ],
            "description": (
                "Missing value imputation method. "
                "bpca/impseq/impseqrob require an R runtime via rpy2."
            ),
        },
        "ensemble": {
            "type": "string",
            "enum": [
                "none",
                "limma,deqms,proda",
                "limma,rots,deqms",
                "limma,rots,deqms,proda",
            ],
            "description": (
                "Ensemble DE methods (comma-separated) for top-k consensus. "
                "Use 'none' unless de_method == 'ensemble'."
            ),
        },
        "ensemble_k": {
            "type": "integer",
            "minimum": 1,
            "maximum": 5,
            "description": (
                "Minimum number of ensemble methods that must agree "
                "(top-k consensus). Only consumed when ensemble != 'none'."
            ),
        },
        "log2fc_threshold": {
            "type": "number",
            "description": "log2 fold-change threshold",
            "minimum": 0,
            "maximum": 10,
        },
        "reasoning": {
            "type": "string",
            "description": "Why this config is worth testing",
        },
        "expected_outcome": {
            "type": "string",
            "description": "Expected result from this config",
        },
    },
    "required": [
        "name",
        "de_method",
        "fdr_method",
        "normalization",
        "imputation",
        "ensemble",
        "ensemble_k",
        "log2fc_threshold",
        "reasoning",
        "expected_outcome",
    ],
    "additionalProperties": False,
}

PROPOSAL_TOOL: dict[str, Any] = {
    "type": "function",
    "function": {
        "name": "propose_configs",
        "strict": True,
        "description": "Propose 3-5 analysis configurations to test",
        "parameters": {
            "type": "object",
            "properties": {
                "configs": {
                    "type": "array",
                    "items": _CONFIG_ITEM_SCHEMA,
                    "description": "List of candidate configs",
                },
            },
            "required": ["configs"],
            "additionalProperties": False,
        },
    },
}

REFLECTION_TOOL: dict[str, Any] = {
    "type": "function",
    "function": {
        "name": "reflect_on_results",
        "strict": True,
        "description": "Analyze results and decide next steps",
        "parameters": {
            "type": "object",
            "properties": {
                "convergence": {
                    "type": "boolean",
                    "description": "True if optimal config is clear",
                },
                "analysis": {
                    "type": "string",
                    "description": "Analysis of results and patterns",
                },
                "next_configs": {
                    "type": "array",
                    "items": _CONFIG_ITEM_SCHEMA,
                    "description": "Refined configs for next round (empty if converged)",
                },
                "adjustments": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "What changed vs previous round",
                },
            },
            "required": [
                "convergence",
                "analysis",
                "next_configs",
                "adjustments",
            ],
            "additionalProperties": False,
        },
    },
}


# ---------------------------------------------------------------------------
# Client creation
# ---------------------------------------------------------------------------


def create_client(config: AgenticConfig):
    """Create an OpenAI client pointing at the configured base URL.

    Parameters
    ----------
    config : AgenticConfig
        Must have ``llm_base_url`` and ``llm_api_key`` set.
    """
    try:
        from openai import OpenAI  # pylint: disable=import-outside-toplevel
    except ImportError as exc:
        raise LLMUnavailableError("LLM requires openai: pip install openai") from exc

    base = config.llm_base_url.rstrip("/")
    return OpenAI(api_key=config.llm_api_key or "", base_url=base)


# ---------------------------------------------------------------------------
# Call helpers
# ---------------------------------------------------------------------------


def _build_common_kwargs(config: AgenticConfig) -> dict[str, Any]:
    """Build kwargs shared across all chat completion calls."""
    kwargs: dict[str, Any] = {
        "model": config.llm_model,
        "max_tokens": config.llm_max_tokens,
    }
    if config.llm_thinking and config.llm_provider == "deepseek":
        kwargs["extra_body"] = {"thinking": {"type": "enabled"}}
        kwargs["reasoning_effort"] = config.llm_reasoning_effort
    else:
        kwargs["temperature"] = config.llm_temperature
    return kwargs


def call_with_tools(
    system_msg: str,
    user_msg: str,
    tools: list[dict[str, Any]],
    config: AgenticConfig,
) -> dict[str, Any]:
    """Call LLM with tool definitions and return parsed tool-call arguments.

    Falls back: strict tool call → JSON mode → raw text extraction.
    """
    client = create_client(config)

    # --- Attempt 1: tool call ---
    try:
        kwargs = _build_common_kwargs(config)
        kwargs["messages"] = [
            {"role": "system", "content": system_msg},
            {"role": "user", "content": user_msg},
        ]
        kwargs["tools"] = tools
        response = client.chat.completions.create(**kwargs)
        msg = response.choices[0].message
        if msg.tool_calls:
            raw = msg.tool_calls[0].function.arguments
            return json.loads(raw)
        # Model replied with text instead of tool call – fall through
        logger.warning("Tool call returned text, falling back")
    except Exception as exc:  # noqa: BLE001
        logger.warning("Tool call failed (%s), falling back", exc)

    # --- Attempt 2: JSON mode ---
    try:
        kwargs = _build_common_kwargs(config)
        kwargs["messages"] = [
            {"role": "system", "content": system_msg},
            {"role": "user", "content": user_msg + "\n\nRespond with JSON only."},
        ]
        kwargs["response_format"] = {"type": "json_object"}
        response = client.chat.completions.create(**kwargs)
        text = response.choices[0].message.content or ""
        return json.loads(text)
    except Exception as exc:  # noqa: BLE001
        logger.warning("JSON mode failed (%s), falling back to raw parse", exc)

    # --- Attempt 3: raw text extraction ---
    kwargs = _build_common_kwargs(config)
    kwargs["messages"] = [
        {"role": "system", "content": system_msg},
        {"role": "user", "content": user_msg},
    ]
    response = client.chat.completions.create(**kwargs)
    text = response.choices[0].message.content or ""
    return _extract_json(text)


def call_text(
    system_msg: str,
    user_msg: str,
    config: AgenticConfig,
) -> str:
    """Call LLM and return plain text response (for reports)."""
    client = create_client(config)
    kwargs = _build_common_kwargs(config)
    kwargs["messages"] = [
        {"role": "system", "content": system_msg},
        {"role": "user", "content": user_msg},
    ]
    response = client.chat.completions.create(**kwargs)
    return response.choices[0].message.content or ""


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


def _extract_json(text: str) -> dict[str, Any]:
    """Best-effort JSON extraction from raw LLM text."""
    text = text.strip()
    # Try object first, then array wrapped in object
    for start_char, end_char in [("{", "}"), ("[", "]")]:
        start = text.find(start_char)
        end = text.rfind(end_char) + 1
        if start >= 0 and end > start:
            try:
                parsed = json.loads(text[start:end])
                if isinstance(parsed, list):
                    return {"configs": parsed}
                return parsed
            except json.JSONDecodeError:
                continue
    raise ValueError(f"No valid JSON found in LLM response: {text[:200]}")
