"""
Plugin registry for the mokume package.

This module provides a central registry for all mokume extension types.
Built-in methods register via decorators at import time. Third-party
packages register via Python entry points, discovered on first access.

Extension groups:
    - quantification: Protein quantification algorithms
    - normalization.feature: Feature-level normalization methods
    - normalization.sample: Sample/peptide-level normalization methods
    - harmonization: Batch effect correction methods
    - imputation: Missing value imputation methods
    - filter: Quality control filters

Example — registering a built-in method::

    from mokume.core.registry import PluginRegistry

    @PluginRegistry.register("quantification", "directlfq")
    class DirectLFQQuantification(QuantificationMethod):
        ...

Example — third-party registration via pyproject.toml::

    [project.entry-points."mokume.quantification"]
    spectral_counting = "my_package:SpectralCountingMethod"
"""

import importlib.metadata
import logging
import re
from typing import Any, Dict, List, Optional, Set, Tuple, Type

from mokume import _plugins

logger = logging.getLogger(__name__)

# Sentinel for TopN pattern matching
_TOPN_PATTERN = re.compile(r"^top(\d+)$")

# Valid input_level values for quantification methods.
# Must match keys in FLOW_DISPATCH (runner.py) and QpxDataset._VALID_LEVELS.
VALID_INPUT_LEVELS: Set[str] = {"peptides", "psms", "peptides_raw", "features"}


class PluginRegistry:
    """Central registry for all mokume extension types.

    Manages registration and discovery of plugins across five extension
    groups. Supports both decorator-based registration (built-in) and
    entry-point discovery (third-party packages).
    """

    _stores: Dict[str, Dict[str, Any]] = {
        "quantification": {},
        "normalization.feature": {},
        "normalization.sample": {},
        "harmonization": {},
        "imputation": {},
        "filter": {},
    }

    _discovered: Set[Tuple[str, Optional[str]]] = set()
    _builtins: Dict[Tuple[str, str], Any] = {}
    _entry_points_cache: Dict[str, List[Any]] = {}

    @classmethod
    def register(cls, group: str, name: str):
        """Decorator to register a plugin class.

        Parameters
        ----------
        group : str
            Extension group (e.g., "quantification", "normalization.feature").
        name : str
            Name to register the plugin under (e.g., "maxlfq").

        Returns
        -------
        Callable
            Decorator that registers the class and returns it unchanged.

        Raises
        ------
        ValueError
            If the group is not recognized.

        Examples
        --------
        >>> @PluginRegistry.register("quantification", "my_method")
        ... class MyMethod(QuantificationMethod):
        ...     ...
        """
        if group not in cls._stores:
            raise ValueError(
                f"Unknown plugin group: '{group}'. "
                f"Available groups: {list(cls._stores.keys())}"
            )

        def decorator(klass: Type) -> Type:
            key = (group, name.lower())
            store = cls._stores[group]
            existing = store.get(key[1])
            if key in _plugins.plugins_for_module(klass.__module__):
                previous_builtin = cls._builtins.get(key)
                cls._builtins[key] = klass
                if existing is not None and existing is not previous_builtin:
                    return klass
            store[key[1]] = klass
            return klass

        return decorator

    @classmethod
    def register_instance_factory(cls, group: str, name: str, factory):
        """Register a callable factory for a plugin.

        Useful for registering aliases like top3/top5 that create
        instances with specific parameters.

        Parameters
        ----------
        group : str
            Extension group.
        name : str
            Name to register under.
        factory : callable
            A callable that accepts **kwargs and returns a plugin instance.
        """
        if group not in cls._stores:
            raise ValueError(
                f"Unknown plugin group: '{group}'. "
                f"Available groups: {list(cls._stores.keys())}"
            )
        cls._stores[group][name.lower()] = factory

    @classmethod
    def get(cls, group: str, name: str, **kwargs: Any) -> Any:
        """Get a plugin instance by group and name.

        Handles special patterns like topN (top3, top5, top10) by
        parsing the numeric suffix.

        Parameters
        ----------
        group : str
            Extension group.
        name : str
            Plugin name.
        **kwargs
            Arguments passed to the plugin constructor.

        Returns
        -------
        Any
            An instance of the requested plugin.

        Raises
        ------
        ValueError
            If the plugin is not found.
        """
        name_lower = name.lower()
        cls._ensure_registered(group, name_lower)

        # Check direct match first
        entry = cls._stores.get(group, {}).get(name_lower)
        instance = None
        if entry is not None:
            if isinstance(entry, type):
                instance = entry(**kwargs)
            else:
                # It's a factory callable
                instance = entry(**kwargs)
        else:
            # Handle topN pattern: top3, top5, top10, etc.
            if group == "quantification":
                match = _TOPN_PATTERN.match(name_lower)
                if match:
                    cls._ensure_registered(group, "topn")
                    topn_cls = cls._stores.get(group, {}).get("topn")
                    if topn_cls is not None:
                        n = int(match.group(1))
                        instance = topn_cls(n=n, **kwargs)

        if instance is None:
            available = cls.available(group)
            raise ValueError(
                f"Unknown {group} method: '{name}'. Available: {available}"
            )

        # Validate input_level for quantification methods
        if group == "quantification" and hasattr(instance, "input_level"):
            level = instance.input_level
            if level not in VALID_INPUT_LEVELS:
                raise ValueError(
                    f"Quantification method '{name}' declares "
                    f"input_level='{level}', which is not valid. "
                    f"Must be one of: {sorted(VALID_INPUT_LEVELS)}"
                )

        return instance

    @classmethod
    def get_class(cls, group: str, name: str) -> Optional[Type]:
        """Get the registered class (not an instance) for a plugin.

        Parameters
        ----------
        group : str
            Extension group.
        name : str
            Plugin name.

        Returns
        -------
        Type or None
            The registered class, or None if not found.
        """
        cls._ensure_registered(group, name.lower())
        return cls._stores.get(group, {}).get(name.lower())

    @classmethod
    def available(cls, group: str) -> List[str]:
        """List registered plugin names for a group.

        Parameters
        ----------
        group : str
            Extension group.

        Returns
        -------
        list[str]
            Sorted list of available plugin names.
        """
        for module_name in _plugins.load_all_plugins(group):
            cls._restore_module(module_name)
        cls._ensure_discovered(group)
        return sorted(cls._stores.get(group, {}).keys())

    @classmethod
    def is_registered(cls, group: str, name: str) -> bool:
        """Check if a plugin is registered.

        Parameters
        ----------
        group : str
            Extension group.
        name : str
            Plugin name.

        Returns
        -------
        bool
        """
        cls._ensure_registered(group, name.lower())
        return name.lower() in cls._stores.get(group, {})

    @classmethod
    def _ensure_registered(cls, group: str, name: str) -> None:
        """Load one matching built-in and third-party plugin."""
        module_name = _plugins.load_plugin(group, name)
        if module_name is not None:
            cls._restore_module(module_name)
        cls._ensure_discovered(group, name)

    @classmethod
    def _restore_module(cls, module_name: str) -> None:
        """Remember or restore built-ins from one imported module."""
        for group, name in _plugins.plugins_for_module(module_name):
            entry = cls._stores[group].get(name)
            if entry is not None and getattr(entry, "__module__", None) == module_name:
                cls._builtins[(group, name)] = entry
            elif entry is None and (group, name) in cls._builtins:
                cls._stores[group][name] = cls._builtins[(group, name)]

    @classmethod
    def _entry_points(cls, group: str) -> List[Any]:
        """Return the installed entry points for one plugin group."""
        if group in cls._entry_points_cache:
            return cls._entry_points_cache[group]

        ep_group = f"mokume.{group}"
        try:
            entry_points = importlib.metadata.entry_points(group=ep_group)
        except TypeError:
            all_entry_points = importlib.metadata.entry_points()
            if isinstance(all_entry_points, dict):
                entry_points = all_entry_points.get(ep_group, [])
            else:
                entry_points = [
                    entry_point
                    for entry_point in all_entry_points
                    if entry_point.group == ep_group
                ]

        cls._entry_points_cache[group] = list(entry_points)
        return cls._entry_points_cache[group]

    @classmethod
    def _ensure_discovered(
        cls,
        group: str,
        name: Optional[str] = None,
    ) -> None:
        """Discover matching entry-point plugins once on first access."""
        key = (group, name)
        if (group, None) in cls._discovered or key in cls._discovered:
            return
        cls._discovered.add(key)

        ep_group = f"mokume.{group}"
        try:
            group_eps = cls._entry_points(group)
            if name is not None:
                group_eps = [ep for ep in group_eps if ep.name.lower() == name]
            for ep in group_eps:
                try:
                    klass = ep.load()
                    cls._stores[group][ep.name.lower()] = klass
                    logger.debug(
                        "Discovered plugin: %s.%s -> %s",
                        group,
                        ep.name,
                        klass,
                    )
                except Exception as exc:
                    logger.warning(
                        "Failed to load plugin '%s' from group '%s': %s",
                        ep.name,
                        ep_group,
                        exc,
                    )
        except Exception as exc:
            logger.debug(
                "Entry point discovery failed for group '%s': %s",
                ep_group,
                exc,
            )

    @classmethod
    def reset(cls):
        """Reset the registry. Mainly useful for testing."""
        for group in cls._stores:
            cls._stores[group].clear()
        cls._discovered.clear()
        cls._entry_points_cache.clear()
