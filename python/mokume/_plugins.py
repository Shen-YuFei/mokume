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
    ("quantification", "ibaq"): "mokume.quantification._dev_registrations",
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


def _load_module(module_name: str) -> bool:
    """Import one plugin module, allowing its decorators to register classes."""
    try:
        importlib.import_module(module_name)
    except ImportError as exc:
        logger.debug(
            "Skipping plugin module '%s' (optional dependency missing): %s",
            module_name,
            exc,
        )
        return False

    from mokume.core.registry import PluginRegistry

    for (group, name), owner in _PLUGIN_MODULES.items():
        if owner != module_name:
            continue
        PluginRegistry._restore_builtin(group, name, module_name)
    return True


def load_plugin(group: str, name: str) -> bool:
    """Load the built-in module that owns one plugin name."""
    module_name = _PLUGIN_MODULES.get((group, name.lower()))
    return _load_module(module_name) if module_name is not None else False


def load_all_plugins(group: str | None = None) -> None:
    """Load every built-in plugin in one group, or in every group."""
    modules = dict.fromkeys(
        module
        for (plugin_group, _), module in _PLUGIN_MODULES.items()
        if group is None or plugin_group == group
    )
    for module_name in modules:
        _load_module(module_name)
