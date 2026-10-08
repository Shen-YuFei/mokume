"""Loopback session, Origin, and CSRF protection for Mokume Studio."""

from __future__ import annotations

import secrets
from dataclasses import dataclass
from ipaddress import ip_address
from urllib.parse import urlsplit


SESSION_COOKIE = "mokume_studio_session"
CSRF_HEADER = "x-csrf-token"


@dataclass(frozen=True)
class Session:
    """Server-side session material never persisted to the browser as JSON."""

    id: str
    csrf_token: str


class SessionManager:
    """Exchange one startup token for in-memory browser sessions."""

    def __init__(self, startup_token: str | None = None) -> None:
        self.startup_token = startup_token or secrets.token_urlsafe(32)
        self._startup_token_available = True
        self._sessions: dict[str, Session] = {}

    def exchange_startup_token(self, supplied: str) -> Session | None:
        """Consume the startup token and return a new browser session."""
        valid = self._startup_token_available and secrets.compare_digest(
            supplied, self.startup_token
        )
        if not valid:
            return None
        self._startup_token_available = False
        session = Session(
            id=secrets.token_urlsafe(32), csrf_token=secrets.token_urlsafe(32)
        )
        self._sessions[session.id] = session
        return session

    def get(self, session_id: str | None) -> Session | None:
        """Return a session only for an exact opaque-cookie match."""
        if not session_id:
            return None
        return self._sessions.get(session_id)


def origin_is_allowed(origin: str | None, expected_origin: str) -> bool:
    """Accept only the exact loopback origin for state-changing requests."""
    return bool(origin) and secrets.compare_digest(origin, expected_origin)


def loopback_origin_from_host(host: str | None) -> str | None:
    """Return a canonical origin for one explicit loopback Host authority."""
    parsed = _parse_host_authority(host)
    if parsed is None:
        return None
    hostname, port = parsed
    authority = _loopback_authority(hostname, port)
    return f"http://{authority}" if authority else None


def _parse_host_authority(host: str | None) -> tuple[str, int] | None:
    if not host or host != host.strip() or host.endswith(":"):
        return None
    try:
        parsed = urlsplit(f"//{host}")
        port = parsed.port
    except ValueError:
        return None
    has_extra_components = any(
        (
            parsed.username is not None,
            parsed.password is not None,
            bool(parsed.path),
            bool(parsed.query),
            bool(parsed.fragment),
        )
    )
    if parsed.hostname is None or port is None or port < 1 or has_extra_components:
        return None
    return parsed.hostname.casefold(), port


def _loopback_authority(hostname: str, port: int) -> str | None:
    if hostname == "localhost":
        return f"localhost:{port}"
    if "%" in hostname:
        return None
    try:
        address = ip_address(hostname)
    except ValueError:
        return None
    if not address.is_loopback:
        return None
    literal = f"[{address.compressed}]" if address.version == 6 else address.compressed
    return f"{literal}:{port}"


def csrf_is_allowed(supplied: str | None, session: Session) -> bool:
    """Compare the request CSRF token without timing-dependent early exits."""
    return bool(supplied) and secrets.compare_digest(supplied, session.csrf_token)
