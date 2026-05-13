"""Tests for mokume.agentic.config.AgenticConfig defaults."""

from mokume.agentic.config import AgenticConfig


def test_thinking_enabled_by_default():
    """DeepSeek V4 Pro thinking mode is on by default."""
    cfg = AgenticConfig()
    assert cfg.llm_thinking is True
    assert cfg.llm_reasoning_effort == "high"


def test_thinking_can_be_disabled():
    """Users can opt out of thinking mode."""
    cfg = AgenticConfig(llm_thinking=False)
    assert cfg.llm_thinking is False


def test_reasoning_effort_max():
    """reasoning_effort accepts the 'max' tier."""
    cfg = AgenticConfig(llm_reasoning_effort="max")
    assert cfg.llm_reasoning_effort == "max"


def test_deepseek_preset_applied():
    """deepseek provider auto-populates base URL and model."""
    cfg = AgenticConfig()
    assert cfg.llm_provider == "deepseek"
    assert cfg.llm_base_url == "https://api.deepseek.com"
    assert cfg.llm_model == "deepseek-v4-pro"


def test_max_tokens_defaults_to_v4_pro_cap():
    """Default max_tokens matches DeepSeek V4 Pro's 384K output cap."""
    cfg = AgenticConfig()
    assert cfg.llm_max_tokens == 384000


def test_max_tokens_can_be_overridden():
    """Custom-provider users with smaller caps can override programmatically."""
    cfg = AgenticConfig(llm_max_tokens=8192)
    assert cfg.llm_max_tokens == 8192
