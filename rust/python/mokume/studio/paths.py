"""Local directory browsing and project-root path guards."""

from __future__ import annotations

import os
from pathlib import Path

from mokume.studio.models import FolderEntry


class PathAccessError(ValueError):
    """Raised when a requested path is unreadable or escapes the project."""


def readable_directory(path: str | Path) -> Path:
    """Resolve a path and require a readable directory."""
    try:
        resolved = Path(path).expanduser().resolve(strict=True)
    except (OSError, RuntimeError) as exc:
        raise PathAccessError(f"directory is unavailable: {path}") from exc
    if not resolved.is_dir() or not os.access(resolved, os.R_OK | os.X_OK):
        raise PathAccessError(f"directory is not readable: {resolved}")
    return resolved


def list_directories(path: str | Path | None = None) -> tuple[Path, list[FolderEntry]]:
    """List readable child directories, starting at the user's home directory."""
    directory = readable_directory(path or Path.home())
    entries: list[FolderEntry] = []
    for child in directory_children(directory):
        try:
            resolved = child.resolve(strict=True)
        except (OSError, RuntimeError):
            continue
        if resolved.is_dir() and os.access(resolved, os.R_OK | os.X_OK):
            entries.append(FolderEntry(name=child.name, path=str(resolved)))
    return directory, entries


def directory_children(directory: Path) -> list[Path]:
    """Return name-sorted children or a stable path-access error."""
    try:
        return sorted(directory.iterdir(), key=lambda item: item.name.casefold())
    except OSError as exc:
        raise PathAccessError(f"cannot list directory: {directory}") from exc


class ProjectPaths:
    """Resolve every file operation beneath one immutable project root."""

    def __init__(self, root: str | Path) -> None:
        self.root = readable_directory(root)

    def resolve_existing(self, requested: str | Path) -> Path:
        """Resolve an existing project-relative path without symlink escape."""
        candidate = self._candidate(requested)
        try:
            resolved = candidate.resolve(strict=True)
        except (OSError, RuntimeError) as exc:
            raise PathAccessError(f"path is unavailable: {requested}") from exc
        self._require_beneath_root(resolved)
        return resolved

    def resolve_output(
        self, requested: str | Path, *, allow_existing: bool = False
    ) -> Path:
        """Resolve a prospective output path and reject project-root escapes."""
        candidate = self._candidate(requested)
        existing_ancestor = next(
            (parent for parent in (candidate, *candidate.parents) if parent.exists()),
            None,
        )
        if existing_ancestor is None:  # pragma: no cover - filesystem roots exist
            raise PathAccessError(f"output parent is unavailable: {requested}")
        try:
            resolved_ancestor = existing_ancestor.resolve(strict=True)
            suffix = candidate.relative_to(existing_ancestor)
        except (OSError, RuntimeError, ValueError) as exc:
            raise PathAccessError(f"output path is unavailable: {requested}") from exc
        self._require_beneath_root(resolved_ancestor)
        resolved = resolved_ancestor.joinpath(suffix)
        if resolved.exists() and not allow_existing:
            raise PathAccessError(f"output already exists: {resolved}")
        return resolved

    def relative(self, path: str | Path) -> str:
        """Return a stable project-relative representation of an existing path."""
        return str(self.resolve_existing(path).relative_to(self.root))

    def _candidate(self, requested: str | Path) -> Path:
        path = Path(requested).expanduser()
        return path if path.is_absolute() else self.root / path

    def _require_beneath_root(self, resolved: Path) -> None:
        try:
            resolved.relative_to(self.root)
        except ValueError as exc:
            raise PathAccessError(f"path escapes project root: {resolved}") from exc
