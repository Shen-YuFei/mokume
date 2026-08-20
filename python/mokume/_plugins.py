"""Lazy loading for Mokume's built-in plugins.

Built-in classes keep their existing ``@PluginRegistry.register`` decorators.
The registry calls :func:`load_plugin` before resolving one name, so importing
or running a lightweight method does not load unrelated implementations.
Listing a group calls :func:`load_all_plugins` for that group and preserves the
existing availability contract.

Each owning module is imported in isolation. A missing optional dependency
skips that module as before, while any other import-time exception still
signals a plugin bug.
"""

import importlib
import logging

logger = logging.getLogger(__name__)

# Importing an owning module runs its existing @PluginRegistry.register
# decorators. Several related names can share one module.
_PLUGIN_MODULES = {
    ("quantification", "sum"): "mokume.quantification.all_peptides",
    ("quantification", "directlfq"): "mokume.quantification.directlfq",
    ("quantification", "maxlfq"): "mokume.quantification.maxlfq",
    ("quantification", "median"): "mokume.quantification.median",
    ("quantification", "ratio"): "mokume.quantification.ratio",
    ("quantification", "spectral_count"): "mokume.quantification.spectral_count",
    ("quantification", "tmt_abundance"): "mokume.quantification.tmt_abundance",
    ("quantification", "tmt_reporter"): "mokume.quantification.tmt_reporter",
    ("quantification", "topn"): "mokume.quantification.topn",
    ("quantification", "pibaq"): "mokume.quantification._dev_registrations",
    ("imputation", "mean"): "mokume.imputation.simple",
    ("imputation", "median"): "mokume.imputation.simple",
    ("imputation", "most_frequent"): "mokume.imputation.simple",
    ("imputation", "constant"): "mokume.imputation.simple",
    ("imputation", "knn"): "mokume.imputation.knn",
    ("imputation", "bpca"): "mokume.imputation._dev_registrations",
    ("imputation", "qrilc"): "mokume.imputation._dev_registrations",
    ("imputation", "seqknn"): "mokume.imputation._dev_registrations",
    ("imputation", "impseq"): "mokume.imputation._dev_registrations",
    ("imputation", "impseqrob"): "mokume.imputation._dev_registrations",
    ("imputation", "gms"): "mokume.imputation._dev_registrations",
    ("imputation", "missforest"): "mokume.imputation._dev_registrations",
    ("normalization.feature", "none"): "mokume.normalization.feature_normalizers",
    ("normalization.feature", "mean"): "mokume.normalization.feature_normalizers",
    ("normalization.feature", "median"): "mokume.normalization.feature_normalizers",
    ("normalization.feature", "max"): "mokume.normalization.feature_normalizers",
    ("normalization.feature", "global"): "mokume.normalization.feature_normalizers",
    ("normalization.feature", "max_min"): "mokume.normalization.feature_normalizers",
    ("normalization.feature", "iqr"): "mokume.normalization.feature_normalizers",
    ("normalization.feature", "loess"): "mokume.normalization._dev_registrations",
    ("normalization.feature", "quantile"): "mokume.normalization._dev_registrations",
    ("normalization.sample", "none"): "mokume.normalization.sample_normalizers",
    (
        "normalization.sample",
        "globalmedian",
    ): "mokume.normalization.sample_normalizers",
    (
        "normalization.sample",
        "conditionmedian",
    ): "mokume.normalization.sample_normalizers",
    (
        "normalization.sample",
        "hierarchical",
    ): "mokume.normalization.sample_normalizers",
    ("normalization.sample", "tmm"): "mokume.normalization.sample_normalizers",
    ("normalization.sample", "irs"): "mokume.normalization.sample_normalizers",
    ("normalization.sample", "rlr"): "mokume.normalization._dev_registrations",
    ("harmonization", "combat"): "mokume.harmonization.combat",
}


def _load_module(module_name: str) -> str | None:
    """Import one plugin module, allowing its decorators to register classes."""
    try:
        importlib.import_module(module_name)
    except ModuleNotFoundError as exc:
        missing_name = exc.name or ""
        if missing_name == module_name or module_name.startswith(f"{missing_name}."):
            logger.error(
                "Plugin module '%s' could not be found: %s",
                module_name,
                exc,
            )
        else:
            logger.debug(
                "Skipping plugin module '%s' (optional dependency missing): %s",
                module_name,
                exc,
            )
        return None
    except ImportError as exc:
        logger.debug(
            "Skipping plugin module '%s' (optional dependency missing): %s",
            module_name,
            exc,
        )
        return None
    return module_name


def load_plugin(group: str, name: str) -> str | None:
    """Load the built-in module that owns one plugin name."""
    module_name = _PLUGIN_MODULES.get((group, name.lower()))
    return _load_module(module_name) if module_name is not None else None


def load_all_plugins(group: str | None = None) -> list[str]:
    """Load every built-in plugin in one group, or in every group."""
    modules = dict.fromkeys(
        module
        for (plugin_group, _), module in _PLUGIN_MODULES.items()
        if group is None or plugin_group == group
    )
    return [
        loaded
        for module_name in modules
        if (loaded := _load_module(module_name)) is not None
    ]


def plugins_for_module(module_name: str) -> tuple[tuple[str, str], ...]:
    """Return the plugin keys owned by one implementation module."""
    return tuple(
        plugin for plugin, owner in _PLUGIN_MODULES.items() if owner == module_name
    )
