"""Loopback session, Origin, and CSRF protection for Mokume Studio."""

from __future__ import annotations

import secrets
from dataclasses import dataclass


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


def csrf_is_allowed(supplied: str | None, session: Session) -> bool:
    """Compare the request CSRF token without timing-dependent early exits."""
    return bool(supplied) and secrets.compare_digest(supplied, session.csrf_token)
