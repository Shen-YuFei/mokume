"""In-memory AI provider configuration for Mokume Studio sessions."""

from __future__ import annotations

from threading import RLock
from typing import Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    SecretStr,
    ValidationInfo,
    field_validator,
)
from pydantic_ai.models import Model
from pydantic_ai.models.anthropic import AnthropicModel
from pydantic_ai.models.openai import OpenAIChatModel, OpenAIResponsesModel
from pydantic_ai.providers.anthropic import AnthropicProvider
from pydantic_ai.providers.openai import OpenAIProvider


ProviderName = Literal["openai", "openai-compatible", "anthropic"]
OPENAI_COMPATIBLE_PLACEHOLDER_KEY = "mokume-local-no-api-key"


class ProviderConfig(BaseModel):
    """Validated provider settings retained only in server memory."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    provider: ProviderName
    model: str = Field(min_length=1)
    api_key: SecretStr | None = Field(default=None, repr=False)
    base_url: str | None = Field(default=None, validate_default=True)

    @field_validator("model")
    @classmethod
    def normalize_model(cls, value: str) -> str:
        """Reject model names containing only whitespace."""
        normalized = value.strip()
        if not normalized:
            raise ValueError("model must not be blank")
        return normalized

    @field_validator("api_key", mode="before")
    @classmethod
    def reject_blank_api_key(cls, value: object) -> object:
        """Treat an explicitly supplied blank credential as invalid."""
        if value is None:
            return None
        secret = (
            value.get_secret_value() if isinstance(value, SecretStr) else str(value)
        )
        if not secret.strip():
            raise ValueError("api_key must not be blank")
        return value

    @field_validator("base_url")
    @classmethod
    def validate_base_url(cls, value: str | None, info: ValidationInfo) -> str | None:
        """Allow custom endpoints only for OpenAI-compatible providers."""
        compatible = info.data.get("provider") == "openai-compatible"
        if value is None:
            if compatible:
                raise ValueError("base_url is required for openai-compatible providers")
            return None
        normalized = value.strip()
        if not normalized:
            raise ValueError("base_url must not be blank")
        if not compatible:
            raise ValueError("base_url is only valid for openai-compatible providers")
        return normalized.rstrip("/")


class ProviderSummary(BaseModel):
    """Credential-free provider settings safe to return to the browser."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    provider: ProviderName
    model: str
    base_url: str | None = None
    api_key_configured: bool


class ProviderRegistry:
    """Thread-safe, process-local provider settings keyed by Studio session."""

    def __init__(self) -> None:
        self._configs: dict[str, ProviderConfig] = {}
        self._lock = RLock()

    def save(self, session_id: str, config: ProviderConfig) -> ProviderSummary:
        """Replace one session's provider configuration and return its summary."""
        self._validate_session_id(session_id)
        with self._lock:
            self._configs[session_id] = config
        return self._summary(config)

    def summary(self, session_id: str) -> ProviderSummary | None:
        """Return credential-free configuration for one session, if present."""
        self._validate_session_id(session_id)
        with self._lock:
            config = self._configs.get(session_id)
        return self._summary(config) if config is not None else None

    def clear(self, session_id: str) -> bool:
        """Forget one session's credentials, reporting whether any existed."""
        self._validate_session_id(session_id)
        with self._lock:
            return self._configs.pop(session_id, None) is not None

    def model_for(self, session_id: str) -> Model:
        """Construct the configured PydanticAI model without persisting secrets."""
        self._validate_session_id(session_id)
        with self._lock:
            config = self._configs.get(session_id)
        if config is None:
            raise LookupError("No AI provider is configured for this Studio session")
        try:
            return self._build_model(config)
        except Exception as exc:
            if "socksio" in str(exc).lower():
                raise RuntimeError(
                    "A SOCKS proxy is configured, but SOCKS support is unavailable; "
                    "install httpx[socks] or remove the SOCKS proxy setting"
                ) from exc
            raise

    @staticmethod
    def _build_model(config: ProviderConfig) -> Model:
        api_key = config.api_key.get_secret_value() if config.api_key else None
        if config.provider == "openai":
            provider = OpenAIProvider(api_key=api_key)
            return OpenAIResponsesModel(config.model, provider=provider)
        if config.provider == "openai-compatible":
            provider = OpenAIProvider(
                base_url=config.base_url,
                api_key=api_key or OPENAI_COMPATIBLE_PLACEHOLDER_KEY,
            )
            return OpenAIChatModel(config.model, provider=provider)
        provider = AnthropicProvider(api_key=api_key)
        return AnthropicModel(config.model, provider=provider)

    @staticmethod
    def _summary(config: ProviderConfig) -> ProviderSummary:
        return ProviderSummary(
            provider=config.provider,
            model=config.model,
            base_url=config.base_url,
            api_key_configured=config.api_key is not None,
        )

    @staticmethod
    def _validate_session_id(session_id: str) -> None:
        if not session_id.strip():
            raise ValueError("session_id must not be blank")
