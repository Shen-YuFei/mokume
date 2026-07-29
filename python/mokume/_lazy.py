"""Shared support for module-level lazy exports."""

import importlib
from typing import Any, Callable


def import_attribute(exports: dict[str, tuple[str, str]], name: str) -> Any:
    """Import one attribute from a lazy export table."""
    target = exports[name]
    target_module, attribute = target
    return getattr(importlib.import_module(target_module), attribute)


def module_api(
    exports: dict[str, tuple[str, str]],
    namespace: dict[str, Any],
    module_name: str,
) -> tuple[Callable[[str], Any], Callable[[], list[str]]]:
    """Build ``__getattr__`` and ``__dir__`` for lazy module exports."""

    def get_attribute(name: str) -> Any:
        if name not in exports:
            raise AttributeError(f"module {module_name!r} has no attribute {name!r}")
        if name in namespace:
            return namespace[name]
        value = import_attribute(exports, name)
        namespace[name] = value
        return value

    def list_attributes() -> list[str]:
        return sorted(set(namespace) | set(exports))

    return get_attribute, list_attributes
