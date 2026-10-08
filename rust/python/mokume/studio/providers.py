"""Session and opt-in persistent AI provider configuration for Mokume Studio."""

from __future__ import annotations

import json
import logging
import os
import tempfile
from pathlib import Path
from time import perf_counter
from threading import RLock
from typing import Literal, NamedTuple, cast
from urllib.parse import urlsplit

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    SecretStr,
    field_validator,
    model_validator,
)
from pydantic_ai import Agent
from pydantic_ai.models import Model
from pydantic_ai.models.anthropic import AnthropicModel
from pydantic_ai.models.google import GoogleModel
from pydantic_ai.models.openai import OpenAIChatModel, OpenAIResponsesModel
from pydantic_ai.providers.anthropic import AnthropicProvider
from pydantic_ai.providers.google import GoogleProvider
from pydantic_ai.providers.openai import OpenAIProvider
from pydantic_ai.settings import ModelSettings, ThinkingLevel
from pydantic_ai.usage import UsageLimits
from platformdirs import user_config_dir


ProviderName = Literal["anthropic", "openai-chat", "openai-responses", "gemini"]
ThinkingLevelName = str
UNIFIED_THINKING_LEVELS = frozenset({"minimal", "low", "medium", "high", "xhigh"})
OPENAI_COMPATIBLE_PLACEHOLDER_KEY = "mokume-local-no-api-key"
PROVIDER_CONFIG_FILENAME = Path("mokume-studio-providers.json")
PROVIDER_CONFIG_IGNORE_PATTERN = "/mokume-studio-providers.json"
PROVIDER_CONFIG_SCHEMA_VERSION = 1
LOGGER = logging.getLogger(__name__)


def default_provider_config_root() -> Path:
    """Use the Mokume source checkout root or the user's config directory."""
    module_path = Path(__file__).resolve()
    package_root = module_path.parents[1]
    for candidate in package_root.parents:
        source_module = candidate / "rust/python/mokume/studio/providers.py"
        if (
            (candidate / ".git").exists()
            and source_module.is_file()
            and source_module.resolve() == module_path
        ):
            return candidate
    return Path(user_config_dir("mokume", appauthor=False))


class ProviderConfig(BaseModel):
    """Validated provider settings with explicit persistence consent."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    provider: ProviderName
    model: str = Field(min_length=1)
    api_key: SecretStr | None = Field(default=None, repr=False)
    base_url: str | None = Field(default=None, validate_default=True)
    context_tokens: int | None = Field(default=None, ge=1024, le=10_000_000)
    max_output_tokens: int | None = Field(default=None, ge=1, le=1_000_000)
    thinking_level: ThinkingLevelName | None = Field(
        default=None,
        min_length=1,
        max_length=64,
        pattern=r"^[a-z][a-z0-9_-]*$",
    )
    persist: bool = False

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
    def validate_base_url(cls, value: str | None) -> str | None:
        """Validate an optional SDK-compatible HTTP endpoint."""
        if value is None:
            return None
        normalized = value.strip()
        if not normalized:
            raise ValueError("base_url must not be blank")
        parsed = urlsplit(normalized)
        if parsed.scheme not in {"http", "https"} or not parsed.netloc:
            raise ValueError("base_url must be an absolute HTTP(S) URL")
        return normalized.rstrip("/")

    @field_validator("thinking_level", mode="before")
    @classmethod
    def normalize_thinking_level(cls, value: object) -> object:
        """Normalize built-in and provider-specific effort names."""
        return value.strip().lower() if isinstance(value, str) else value

    @model_validator(mode="after")
    def validate_token_limits(self) -> ProviderConfig:
        """Reserve at least one input token when both limits are configured."""
        if (
            self.context_tokens is not None
            and self.max_output_tokens is not None
            and self.max_output_tokens >= self.context_tokens
        ):
            raise ValueError("max_output_tokens must be smaller than context_tokens")
        return self


class ProviderSummary(BaseModel):
    """Credential-free provider settings safe for general Studio state."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    provider: ProviderName
    model: str
    base_url: str | None = None
    api_key_configured: bool
    context_tokens: int | None = None
    max_output_tokens: int | None = None
    thinking_level: ThinkingLevelName | None = None
    persistent: bool = False


class ProviderDetails(ProviderSummary):
    """Provider settings, including the key, for the authenticated config form."""

    api_key: str | None = Field(default=None, repr=False)


class ProviderTestResult(BaseModel):
    """Safe browser result for a live provider capability probe."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    connected: Literal[True] = True
    tool_calling: Literal[True] = True
    latency_ms: int = Field(ge=0)


class ProviderExecution(NamedTuple):
    """One internally consistent snapshot used for an agent run."""

    model: Model
    model_settings: ModelSettings | None
    usage_limits: UsageLimits | None


class ProviderRegistry:
    """Thread-safe session settings with one stable persistence location."""

    def __init__(self, config_root: str | Path | None = None) -> None:
        self._config_path = (
            Path(config_root).resolve() / PROVIDER_CONFIG_FILENAME
            if config_root is not None
            else None
        )
        self._configs: dict[str, ProviderConfig] = {}
        self._lock = RLock()

    def save(self, session_id: str, config: ProviderConfig) -> ProviderSummary:
        """Replace one session's provider configuration and return its summary."""
        self._validate_session_id(session_id)
        config_path = self._config_path
        with self._lock:
            previous = self._load_persistent(config_path) or self._configs.get(
                session_id
            )
            config = self._retain_key(config, previous)
            if config.persist:
                if config_path is None:
                    raise ValueError("Persistent provider storage is unavailable")
                self._write_persistent(config_path, config)
            elif config_path is not None:
                self._remove_persistent(config_path)
            self._configs[session_id] = config.model_copy(update={"persist": False})
        return self._summary(config)

    def summary(self, session_id: str) -> ProviderSummary | None:
        """Return credential-free configuration for one session, if present."""
        self._validate_session_id(session_id)
        config_path = self._config_path
        with self._lock:
            config = self._load_persistent(config_path) or self._configs.get(session_id)
        return self._summary(config) if config is not None else None

    def details(self, session_id: str) -> ProviderDetails | None:
        """Return configuration and its key for the authenticated config form."""
        self._validate_session_id(session_id)
        config_path = self._config_path
        with self._lock:
            config = self._load_persistent(config_path) or self._configs.get(session_id)
        if config is None:
            return None
        return ProviderDetails(
            **self._summary(config).model_dump(),
            api_key=(
                config.api_key.get_secret_value()
                if config.api_key is not None
                else None
            ),
        )

    def clear(self, session_id: str) -> bool:
        """Forget current credentials, reporting whether any existed."""
        self._validate_session_id(session_id)
        config_path = self._config_path
        with self._lock:
            removed = self._configs.pop(session_id, None) is not None
            if config_path is not None and self._remove_persistent(config_path):
                removed = True
            return removed

    def model_for(self, session_id: str) -> Model:
        """Construct the configured PydanticAI model without persisting secrets."""
        return self.execution_for(session_id).model

    def execution_for(self, session_id: str) -> ProviderExecution:
        """Build a model and request options from one session config snapshot."""
        config = self._config_for(session_id)
        try:
            model = self._build_model(config)
        except Exception as exc:
            if "socksio" in str(exc).lower():
                raise RuntimeError(
                    "A SOCKS proxy is configured, but SOCKS support is unavailable; "
                    "install httpx[socks] or remove the SOCKS proxy setting"
                ) from exc
            raise
        return ProviderExecution(
            model=model,
            model_settings=self._model_settings(config),
            usage_limits=self._usage_limits(config),
        )

    @staticmethod
    def _usage_limits(
        config: ProviderConfig,
        *,
        request_limit: int = 50,
    ) -> UsageLimits | None:
        if config.context_tokens is None and request_limit == 50:
            return None
        return UsageLimits(
            request_limit=request_limit,
            per_request_input_tokens_limit=config.context_tokens,
            count_tokens_before_request=(
                config.context_tokens is not None and config.provider != "openai-chat"
            ),
        )

    async def test(
        self,
        session_id: str,
        config: ProviderConfig,
    ) -> ProviderTestResult:
        """Probe authentication, model access, and tool calling without saving config."""
        self._validate_session_id(session_id)
        config_path = self._config_path
        with self._lock:
            previous = self._load_persistent(config_path) or self._configs.get(
                session_id
            )
            config = self._retain_key(config, previous)
        started = perf_counter()
        called = False

        async def confirm_connection() -> str:
            """Confirm that the configured model can invoke a Mokume tool."""
            nonlocal called
            called = True
            return "connected"

        try:
            model = self._build_model(config)
            probe = Agent(
                model=model,
                tools=[confirm_connection],
                instructions=(
                    "This is a connection test. Call confirm_connection exactly once, "
                    "then answer with the single word connected."
                ),
                retries=0,
            )
            settings = dict(self._model_settings(config) or {})
            settings.setdefault("timeout", 30)
            await probe.run(
                "Verify this connection now.",
                model_settings=cast(ModelSettings, settings),
                usage_limits=self._usage_limits(config, request_limit=2),
            )
        except Exception as exc:
            raise RuntimeError(self._safe_error(exc, config)) from exc
        if not called:
            raise RuntimeError(
                "Connection succeeded, but the model did not call the required test tool"
            )
        return ProviderTestResult(
            latency_ms=max(0, round((perf_counter() - started) * 1000))
        )

    @staticmethod
    def _build_model(config: ProviderConfig) -> Model:
        api_key = config.api_key.get_secret_value() if config.api_key else None
        openai_key = api_key or (
            OPENAI_COMPATIBLE_PLACEHOLDER_KEY if config.base_url else None
        )
        openai_kwargs = {"api_key": openai_key}
        if config.base_url:
            openai_kwargs["base_url"] = config.base_url
        if config.provider == "openai-responses":
            provider = OpenAIProvider(**openai_kwargs)
            return OpenAIResponsesModel(config.model, provider=provider)
        if config.provider == "openai-chat":
            provider = OpenAIProvider(**openai_kwargs)
            return OpenAIChatModel(config.model, provider=provider)
        if config.provider == "anthropic":
            anthropic_kwargs = {"api_key": api_key}
            if config.base_url:
                anthropic_kwargs["base_url"] = config.base_url
            provider = AnthropicProvider(**anthropic_kwargs)
            return AnthropicModel(config.model, provider=provider)
        google_kwargs = {"api_key": api_key}
        if config.base_url:
            google_kwargs["base_url"] = config.base_url
        provider = GoogleProvider(**google_kwargs)
        return GoogleModel(config.model, provider=provider)

    @staticmethod
    def _model_settings(config: ProviderConfig) -> ModelSettings | None:
        settings = ModelSettings()
        if config.max_output_tokens is not None:
            settings["max_tokens"] = config.max_output_tokens
        if config.thinking_level is not None:
            level = config.thinking_level
            if level == "off":
                settings["thinking"] = False
            elif level in UNIFIED_THINKING_LEVELS:
                settings["thinking"] = cast(ThinkingLevel, level)
            elif level == "max" and config.provider == "gemini":
                settings["thinking"] = "xhigh"
            elif config.provider.startswith("openai"):
                settings["openai_reasoning_effort"] = level
            elif config.provider == "anthropic":
                settings["anthropic_effort"] = level
            else:
                settings["google_thinking_config"] = {
                    "include_thoughts": True,
                    "thinking_level": level.upper(),
                }
        return settings or None

    @staticmethod
    def _summary(config: ProviderConfig) -> ProviderSummary:
        return ProviderSummary(
            provider=config.provider,
            model=config.model,
            base_url=config.base_url,
            api_key_configured=config.api_key is not None,
            context_tokens=config.context_tokens,
            max_output_tokens=config.max_output_tokens,
            thinking_level=config.thinking_level,
            persistent=config.persist,
        )

    def _config_for(self, session_id: str) -> ProviderConfig:
        self._validate_session_id(session_id)
        config_path = self._config_path
        with self._lock:
            config = self._load_persistent(config_path) or self._configs.get(session_id)
        if config is None:
            raise LookupError("No AI provider is configured for this Studio session")
        return config

    @classmethod
    def _load_persistent(cls, path: Path | None) -> ProviderConfig | None:
        if path is None or not path.is_file():
            return None
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            if payload.get("application") != "mokume":
                return None
            if payload.get("$schemaVersion") != PROVIDER_CONFIG_SCHEMA_VERSION:
                raise ValueError("unsupported schema version")
            provider = payload["providers"][0]
            model = provider["models"][0]
            auth = provider.get("auth")
            api_key = auth.get("value") if auth else None
            return ProviderConfig(
                provider=provider["type"],
                model=model["apiModel"],
                api_key=api_key,
                base_url=provider.get("baseUrl"),
                context_tokens=model.get("contextTokenLimit"),
                max_output_tokens=model.get("maxOutputTokens"),
                thinking_level=model.get("thinkingLevel"),
                persist=True,
            )
        except (
            AttributeError,
            IndexError,
            KeyError,
            OSError,
            TypeError,
            ValueError,
        ):
            LOGGER.warning("Ignoring invalid provider configuration at %s", path)
            return None

    @classmethod
    def _write_persistent(cls, path: Path, config: ProviderConfig) -> None:
        if path.exists() and cls._load_persistent(path) is None:
            raise ValueError(
                "Refusing to overwrite mokume-studio-providers.json because it is "
                "not a valid Mokume configuration"
            )
        cls._ensure_ignored(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        api_key = config.api_key.get_secret_value() if config.api_key else None
        payload = {
            "$schemaVersion": PROVIDER_CONFIG_SCHEMA_VERSION,
            "application": "mokume",
            "providers": [
                {
                    "type": config.provider,
                    "baseUrl": config.base_url,
                    "auth": (
                        {"kind": "apiKey", "value": api_key}
                        if api_key is not None
                        else None
                    ),
                    "models": [
                        {
                            "apiModel": config.model,
                            "thinkingLevel": config.thinking_level,
                            "contextTokenLimit": config.context_tokens,
                            "maxOutputTokens": config.max_output_tokens,
                        }
                    ],
                }
            ],
        }
        file_descriptor, temporary_name = tempfile.mkstemp(
            prefix=f".{path.name}.", dir=path.parent, text=True
        )
        try:
            with os.fdopen(file_descriptor, "w", encoding="utf-8") as handle:
                json.dump(payload, handle, ensure_ascii=False, indent=2)
                handle.write("\n")
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary_name, path)
            path.chmod(0o600)
        finally:
            Path(temporary_name).unlink(missing_ok=True)

    @classmethod
    def _remove_persistent(cls, path: Path) -> bool:
        if not path.exists():
            return False
        if cls._load_persistent(path) is None:
            return False
        path.unlink()
        return True

    @staticmethod
    def _ensure_ignored(path: Path) -> None:
        project_root = path.parent
        if not (project_root / ".git").exists():
            return
        ignore_file = project_root / ".gitignore"
        existing = (
            ignore_file.read_text(encoding="utf-8") if ignore_file.exists() else ""
        )
        if PROVIDER_CONFIG_IGNORE_PATTERN in existing.splitlines():
            return
        separator = "" if not existing or existing.endswith("\n") else "\n"
        with ignore_file.open("a", encoding="utf-8") as handle:
            handle.write(f"{separator}{PROVIDER_CONFIG_IGNORE_PATTERN}\n")

    @staticmethod
    def _retain_key(
        config: ProviderConfig,
        previous: ProviderConfig | None,
    ) -> ProviderConfig:
        if (
            config.api_key is None
            and previous is not None
            and previous.api_key is not None
            and previous.provider == config.provider
            and previous.base_url == config.base_url
        ):
            return config.model_copy(update={"api_key": previous.api_key})
        return config

    @staticmethod
    def _safe_error(exc: Exception, config: ProviderConfig) -> str:
        message = str(exc).strip() or exc.__class__.__name__
        if config.api_key is not None:
            message = message.replace(config.api_key.get_secret_value(), "[redacted]")
        return f"Connection test failed: {message}"

    @staticmethod
    def _validate_session_id(session_id: str) -> None:
        if not session_id.strip():
            raise ValueError("session_id must not be blank")
