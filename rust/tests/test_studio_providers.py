"""Session-local AI provider configuration for Mokume Studio."""

from __future__ import annotations

import uuid

import pytest
from pydantic import ValidationError
from pydantic_ai.models.anthropic import AnthropicModel
from pydantic_ai.models.openai import OpenAIChatModel, OpenAIResponsesModel

from mokume.studio import providers
from mokume.studio.providers import ProviderConfig, ProviderRegistry


CREDENTIAL = uuid.uuid4().hex
PROXY_VARIABLES = (
    "HTTP_PROXY",
    "HTTPS_PROXY",
    "ALL_PROXY",
    "http_proxy",
    "https_proxy",
    "all_proxy",
)


def test_provider_config_validates_contract_without_disclosing_secrets():
    """Config validation is strict and never prints a supplied credential."""
    config = ProviderConfig(provider="openai", model="  gpt-5  ", api_key=CREDENTIAL)
    assert config.model == "gpt-5"
    assert CREDENTIAL not in repr(config)

    with pytest.raises(ValidationError) as error:
        ProviderConfig(
            provider="openai",
            model="gpt-5",
            api_key=CREDENTIAL,
            base_url="http://localhost:8000/v1",
        )
    assert CREDENTIAL not in str(error.value)
    assert "only valid for openai-compatible" in str(error.value)

    with pytest.raises(ValidationError) as error:
        ProviderConfig(
            provider="openai-compatible",
            model="local-model",
            api_key=CREDENTIAL,
        )
    assert CREDENTIAL not in str(error.value)
    assert "base_url is required" in str(error.value)

    with pytest.raises(ValidationError):
        ProviderConfig(provider="anthropic", model="   ")
    with pytest.raises(ValidationError):
        ProviderConfig(provider="openai", model="gpt-5", unexpected=True)


def test_registry_isolates_sessions_and_clears_credentials():
    """One browser session cannot see or reuse another session's settings."""
    registry = ProviderRegistry()
    first = registry.save(
        "session-a",
        ProviderConfig(provider="openai", model="gpt-5", api_key=CREDENTIAL),
    )
    registry.save(
        "session-b",
        ProviderConfig(provider="anthropic", model="claude-sonnet-4-5"),
    )

    assert first.api_key_configured is True
    assert first.model == "gpt-5"
    assert "api_key" not in first.model_dump()
    assert CREDENTIAL not in first.model_dump_json()
    assert CREDENTIAL not in repr(first)
    assert CREDENTIAL not in repr(registry)
    assert registry.summary("session-b").provider == "anthropic"
    assert registry.clear("session-a") is True
    assert registry.summary("session-a") is None
    assert registry.summary("session-b") is not None
    assert registry.clear("session-a") is False


def test_model_for_uses_current_pydantic_ai_model_types(monkeypatch):
    """Each provider maps to the intended PydanticAI 2.36 model class."""
    for variable in PROXY_VARIABLES:
        monkeypatch.delenv(variable, raising=False)
    registry = ProviderRegistry()
    registry.save(
        "openai",
        ProviderConfig(provider="openai", model="gpt-5", api_key="fake-openai"),
    )
    registry.save(
        "compatible",
        ProviderConfig(
            provider="openai-compatible",
            model="local-model",
            base_url="http://127.0.0.1:11434/v1",
        ),
    )
    registry.save(
        "anthropic",
        ProviderConfig(
            provider="anthropic",
            model="claude-sonnet-4-5",
            api_key="fake-anthropic",
        ),
    )

    assert isinstance(registry.model_for("openai"), OpenAIResponsesModel)
    assert isinstance(registry.model_for("compatible"), OpenAIChatModel)
    assert isinstance(registry.model_for("anthropic"), AnthropicModel)


def test_model_for_passes_explicit_keys_and_local_placeholder(monkeypatch):
    """Explicit keys win while keyless compatible endpoints get a placeholder."""
    provider_calls: list[dict[str, str | None]] = []

    def fake_provider(**kwargs):
        provider_calls.append(kwargs)
        return object()

    monkeypatch.setattr(providers, "OpenAIProvider", fake_provider)
    monkeypatch.setattr(
        providers, "OpenAIResponsesModel", lambda model, provider: (model, provider)
    )
    monkeypatch.setattr(
        providers, "OpenAIChatModel", lambda model, provider: (model, provider)
    )
    registry = ProviderRegistry()
    registry.save(
        "official",
        ProviderConfig(provider="openai", model="gpt-5", api_key=CREDENTIAL),
    )
    registry.save(
        "local",
        ProviderConfig(
            provider="openai-compatible",
            model="local-model",
            base_url="http://127.0.0.1:11434/v1/",
        ),
    )

    registry.model_for("official")
    registry.model_for("local")
    assert provider_calls == [
        {"api_key": CREDENTIAL},
        {
            "base_url": "http://127.0.0.1:11434/v1",
            "api_key": providers.OPENAI_COMPATIBLE_PLACEHOLDER_KEY,
        },
    ]


def test_official_provider_delegates_missing_key_to_environment(monkeypatch):
    """A missing official key remains None so the provider can read its env var."""
    calls: list[dict[str, str | None]] = []

    def fake_anthropic_provider(**kwargs):
        calls.append(kwargs)
        return object()

    monkeypatch.setattr(providers, "AnthropicProvider", fake_anthropic_provider)
    monkeypatch.setattr(
        providers, "AnthropicModel", lambda model, provider: (model, provider)
    )
    registry = ProviderRegistry()
    registry.save(
        "session", ProviderConfig(provider="anthropic", model="claude-sonnet-4-5")
    )

    registry.model_for("session")
    assert calls == [{"api_key": None}]


def test_missing_socks_support_has_an_actionable_error(monkeypatch):
    """Proxy dependency failures are surfaced instead of bypassing the proxy."""
    monkeypatch.setattr(
        providers,
        "OpenAIProvider",
        lambda **_kwargs: (_ for _ in ()).throw(
            ImportError("Using SOCKS proxy, but the 'socksio' package is not installed")
        ),
    )
    registry = ProviderRegistry()
    registry.save("session", ProviderConfig(provider="openai", model="gpt-5"))

    with pytest.raises(RuntimeError, match="SOCKS proxy is configured"):
        registry.model_for("session")


def test_model_for_requires_a_configured_nonblank_session():
    """Unconfigured or malformed session identifiers fail explicitly."""
    registry = ProviderRegistry()
    with pytest.raises(LookupError, match="No AI provider"):
        registry.model_for("missing")
    with pytest.raises(ValueError, match="session_id"):
        registry.summary("   ")
