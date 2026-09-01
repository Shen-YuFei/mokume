"""Session and opt-in persistent provider configuration for Mokume Studio."""

from __future__ import annotations

import json
import stat
import uuid
from pathlib import Path

import pytest
from pydantic import ValidationError
from pydantic_ai.models.anthropic import AnthropicModel
from pydantic_ai.models.google import GoogleModel
from pydantic_ai.models.openai import OpenAIChatModel, OpenAIResponsesModel
from pydantic_ai.models.test import TestModel

from mokume.studio import providers
from mokume.studio.providers import (
    ProviderConfig,
    ProviderRegistry,
    default_provider_config_root,
)


CREDENTIAL = uuid.uuid4().hex
PROXY_VARIABLES = (
    "HTTP_PROXY",
    "HTTPS_PROXY",
    "ALL_PROXY",
    "http_proxy",
    "https_proxy",
    "all_proxy",
)


def _persistent_registry(tmp_path):
    config_root = tmp_path / "mokume"
    config_root.mkdir()
    (config_root / ".git").mkdir()
    (config_root / ".gitignore").write_text("*.tmp", encoding="utf-8")
    registry = ProviderRegistry(config_root)
    summary = registry.save(
        "session",
        ProviderConfig(
            provider="anthropic",
            model="model-a",
            api_key=CREDENTIAL,
            base_url="https://api.example.com/anthropic",
            context_tokens=225_000,
            max_output_tokens=12_800,
            thinking_level="high",
            persist=True,
        ),
    )
    return config_root, registry, summary


def test_default_provider_config_root_uses_checkout_or_user_config_directory():
    """Persistent credentials never follow the folder opened for analysis."""
    root = default_provider_config_root()
    module_path = Path(providers.__file__).resolve()
    if (root / ".git").exists():
        assert (
            root / "rust/python/mokume/studio/providers.py"
        ).resolve() == module_path
    else:
        assert root == Path(providers.user_config_dir("mokume", appauthor=False))


def test_pip_install_uses_the_user_config_directory(tmp_path, monkeypatch):
    """Wheel installs never write persistent credentials into site-packages."""
    installed_module = tmp_path / "site-packages/mokume/studio/providers.py"
    user_config = tmp_path / "user-config/mokume"
    monkeypatch.setattr(providers, "__file__", str(installed_module))
    monkeypatch.setattr(
        providers,
        "user_config_dir",
        lambda *_args, **_kwargs: str(user_config),
    )

    assert default_provider_config_root() == user_config


def test_provider_config_validates_contract_without_disclosing_secrets():
    """Config validation is strict and never prints a supplied credential."""
    config = ProviderConfig(
        provider="openai-responses",
        model="  gpt-5  ",
        api_key=CREDENTIAL,
        base_url="https://api.example.com/v1/",
        context_tokens=225_000,
        max_output_tokens=12_800,
        thinking_level="high",
    )
    assert config.model == "gpt-5"
    assert config.base_url == "https://api.example.com/v1"
    assert config.context_tokens == 225_000
    assert config.max_output_tokens == 12_800
    assert config.thinking_level == "high"
    assert CREDENTIAL not in repr(config)

    anthropic = ProviderConfig(
        provider="anthropic",
        model="k3-256k",
        base_url="https://api.kimi.com/coding/",
    )
    assert anthropic.base_url == "https://api.kimi.com/coding"

    with pytest.raises(ValidationError) as error:
        ProviderConfig(
            provider="openai-chat",
            model="local-model",
            api_key=CREDENTIAL,
            base_url="file:///tmp/model",
        )
    assert CREDENTIAL not in str(error.value)
    assert "absolute HTTP(S) URL" in str(error.value)

    with pytest.raises(ValidationError):
        ProviderConfig(provider="anthropic", model="   ")
    with pytest.raises(ValidationError):
        ProviderConfig(provider="openai-responses", model="gpt-5", unexpected=True)
    with pytest.raises(ValidationError, match="smaller than context_tokens"):
        ProviderConfig(
            provider="gemini",
            model="gemini-2.5-flash",
            context_tokens=8_192,
            max_output_tokens=8_192,
        )


def test_registry_isolates_sessions_and_clears_credentials():
    """One browser session cannot see or reuse another session's settings."""
    registry = ProviderRegistry()
    first = registry.save(
        "session-a",
        ProviderConfig(provider="openai-responses", model="gpt-5", api_key=CREDENTIAL),
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
    details = registry.details("session-a")
    assert details.api_key == CREDENTIAL
    assert CREDENTIAL not in repr(details)
    assert CREDENTIAL not in repr(registry)
    assert registry.summary("session-b").provider == "anthropic"
    assert registry.clear("session-a") is True
    assert registry.summary("session-a") is None
    assert registry.summary("session-b") is not None
    assert registry.clear("session-a") is False


def test_registry_keeps_a_key_only_when_the_endpoint_is_unchanged():
    """Blank edits retain a key without leaking it to another API endpoint."""
    registry = ProviderRegistry()
    registry.save(
        "session",
        ProviderConfig(
            provider="anthropic",
            model="model-a",
            api_key=CREDENTIAL,
            base_url="https://api.example.com",
        ),
    )

    retained = registry.save(
        "session",
        ProviderConfig(
            provider="anthropic",
            model="model-b",
            base_url="https://api.example.com",
        ),
    )
    changed = registry.save(
        "session",
        ProviderConfig(
            provider="anthropic",
            model="model-b",
            base_url="https://other.example.com",
        ),
    )

    assert retained.api_key_configured is True
    assert changed.api_key_configured is False


def test_registry_persists_config_only_with_explicit_consent(tmp_path):
    """Write opted-in credentials with restricted access and explicit ownership."""
    config_root, _registry, summary = _persistent_registry(tmp_path)
    config_path = config_root / "mokume-studio-providers.json"
    payload = json.loads(config_path.read_text(encoding="utf-8"))
    assert summary.persistent is True
    assert summary.api_key_configured is True
    assert CREDENTIAL not in summary.model_dump_json()
    assert payload["$schemaVersion"] == 1
    assert payload["application"] == "mokume"
    assert payload["providers"][0]["auth"]["value"] == CREDENTIAL
    assert payload["providers"][0]["models"][0]["apiModel"] == "model-a"
    assert stat.S_IMODE(config_path.stat().st_mode) == 0o600
    assert (config_root / ".gitignore").read_text(encoding="utf-8") == (
        "*.tmp\n/mokume-studio-providers.json\n"
    )


def test_registry_restores_persistent_credentials(tmp_path, monkeypatch):
    """Restore an opted-in credential without exposing it through summaries."""
    config_root, _registry, _summary = _persistent_registry(tmp_path)
    calls = []
    monkeypatch.setattr(
        providers,
        "AnthropicProvider",
        lambda **kwargs: calls.append(kwargs) or object(),
    )
    monkeypatch.setattr(
        providers, "AnthropicModel", lambda model, provider: (model, provider)
    )
    restarted = ProviderRegistry(config_root)
    restored = restarted.summary("new-session")
    restarted.model_for("new-session")
    assert restored.persistent is True
    assert restored.api_key_configured is True
    assert restarted.details("new-session").api_key == CREDENTIAL
    assert calls == [
        {
            "api_key": CREDENTIAL,
            "base_url": "https://api.example.com/anthropic",
        }
    ]


def test_disabling_persistence_keeps_only_the_current_session_key(tmp_path):
    """Remove the file while retaining its key only in the active registry."""
    config_root, _registry, _summary = _persistent_registry(tmp_path)
    config_path = config_root / "mokume-studio-providers.json"
    restarted = ProviderRegistry(config_root)
    memory_only = restarted.save(
        "new-session",
        ProviderConfig(
            provider="anthropic",
            model="model-a",
            base_url="https://api.example.com/anthropic",
            persist=False,
        ),
    )
    assert memory_only.persistent is False
    assert memory_only.api_key_configured is True
    assert restarted.details("new-session").api_key == CREDENTIAL
    assert not config_path.exists()
    assert ProviderRegistry(config_root).summary("after-restart") is None


def test_persistent_provider_requires_a_config_root():
    """Persistence cannot silently select an unspecified credential location."""
    registry = ProviderRegistry()
    with pytest.raises(ValueError, match="storage is unavailable"):
        registry.save(
            "session",
            ProviderConfig(
                provider="anthropic",
                model="model-a",
                api_key=CREDENTIAL,
                persist=True,
            ),
        )


def test_persistence_refuses_to_overwrite_an_unrelated_providers_file(tmp_path):
    """A config-root file with the same name is never replaced or deleted."""
    config_root = tmp_path / "mokume"
    config_root.mkdir()
    config_path = config_root / "mokume-studio-providers.json"
    existing = '{"providers": [{"ownedBy": "another-tool"}]}\n'
    config_path.write_text(existing, encoding="utf-8")
    registry = ProviderRegistry(config_root)

    with pytest.raises(ValueError, match="not a valid Mokume configuration"):
        registry.save(
            "session",
            ProviderConfig(
                provider="anthropic",
                model="model-a",
                api_key=CREDENTIAL,
                persist=True,
            ),
        )

    assert config_path.read_text(encoding="utf-8") == existing
    assert not (config_root / ".gitignore").exists()

    memory_only = registry.save(
        "session",
        ProviderConfig(
            provider="anthropic",
            model="model-a",
            api_key=CREDENTIAL,
            persist=False,
        ),
    )
    assert memory_only.persistent is False
    assert config_path.read_text(encoding="utf-8") == existing


def test_model_for_uses_current_pydantic_ai_model_types(monkeypatch):
    """Each provider maps to the intended PydanticAI 2.36 model class."""
    for variable in PROXY_VARIABLES:
        monkeypatch.delenv(variable, raising=False)
    registry = ProviderRegistry()
    registry.save(
        "responses",
        ProviderConfig(
            provider="openai-responses", model="gpt-5", api_key="fake-openai"
        ),
    )
    registry.save(
        "chat",
        ProviderConfig(
            provider="openai-chat",
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
    registry.save(
        "gemini",
        ProviderConfig(
            provider="gemini",
            model="gemini-2.5-flash",
            api_key="fake-google",
        ),
    )

    assert isinstance(registry.model_for("responses"), OpenAIResponsesModel)
    assert isinstance(registry.model_for("chat"), OpenAIChatModel)
    assert isinstance(registry.model_for("anthropic"), AnthropicModel)
    assert isinstance(registry.model_for("gemini"), GoogleModel)


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
        ProviderConfig(provider="openai-responses", model="gpt-5", api_key=CREDENTIAL),
    )
    registry.save(
        "local",
        ProviderConfig(
            provider="openai-chat",
            model="local-model",
            base_url="http://127.0.0.1:11434/v1/",
        ),
    )
    registry.save(
        "responses-compatible",
        ProviderConfig(
            provider="openai-responses",
            model="responses-model",
            base_url="https://responses.example.com/v1/",
        ),
    )

    registry.model_for("official")
    registry.model_for("local")
    registry.model_for("responses-compatible")
    assert provider_calls == [
        {"api_key": CREDENTIAL},
        {
            "api_key": providers.OPENAI_COMPATIBLE_PLACEHOLDER_KEY,
            "base_url": "http://127.0.0.1:11434/v1",
        },
        {
            "api_key": providers.OPENAI_COMPATIBLE_PLACEHOLDER_KEY,
            "base_url": "https://responses.example.com/v1",
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


def test_anthropic_provider_accepts_a_custom_base_url(monkeypatch):
    """Anthropic-compatible services receive the configured endpoint and key."""
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
        "session",
        ProviderConfig(
            provider="anthropic",
            model="k3-256k",
            api_key=CREDENTIAL,
            base_url="https://api.kimi.com/coding/",
        ),
    )

    registry.model_for("session")
    assert calls == [{"api_key": CREDENTIAL, "base_url": "https://api.kimi.com/coding"}]


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
    registry.save("session", ProviderConfig(provider="openai-responses", model="gpt-5"))

    with pytest.raises(RuntimeError, match="SOCKS proxy is configured"):
        registry.model_for("session")


def test_model_for_requires_a_configured_nonblank_session():
    """Unconfigured or malformed session identifiers fail explicitly."""
    registry = ProviderRegistry()
    with pytest.raises(LookupError, match="No AI provider"):
        registry.model_for("missing")
    with pytest.raises(ValueError, match="session_id"):
        registry.summary("   ")


def test_request_settings_and_context_limit_are_applied_per_session(monkeypatch):
    """Advanced settings reach model requests without becoming global state."""
    monkeypatch.setattr(
        ProviderRegistry,
        "_build_model",
        staticmethod(lambda _config: TestModel()),
    )
    registry = ProviderRegistry()
    summary = registry.save(
        "responses",
        ProviderConfig(
            provider="openai-responses",
            model="gpt-5",
            api_key="fake-openai",
            context_tokens=225_000,
            max_output_tokens=12_800,
            thinking_level="high",
        ),
    )

    assert summary.context_tokens == 225_000
    assert summary.max_output_tokens == 12_800
    assert summary.thinking_level == "high"
    execution = registry.execution_for("responses")
    assert execution.model_settings == {
        "max_tokens": 12_800,
        "thinking": "high",
    }
    limits = execution.usage_limits
    assert limits.per_request_input_tokens_limit == 225_000
    assert limits.count_tokens_before_request is True

    registry.save(
        "chat",
        ProviderConfig(
            provider="openai-chat",
            model="local-model",
            base_url="http://127.0.0.1:11434/v1",
            context_tokens=32_768,
            thinking_level="off",
        ),
    )
    execution = registry.execution_for("chat")
    assert execution.model_settings == {"thinking": False}
    assert execution.usage_limits.count_tokens_before_request is False
