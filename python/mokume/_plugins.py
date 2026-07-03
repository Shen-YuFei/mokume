"""Plugin auto-loading for the mokume package.

Importing this module registers every built-in plugin with
:class:`~mokume.core.registry.PluginRegistry`. Registration happens as an
import side effect: each submodule listed below carries
``@PluginRegistry.register`` decorators (or defines adapter classes that do),
so importing the module runs those decorators exactly once.

The top-level :mod:`mokume` package imports this module on first import
(inside a ``try``/``except`` so an optional-dependency ``ImportError`` in any
one plugin family never aborts ``import mokume``). Callers who only need a
subset can also import the relevant submodule directly.

Three kinds of registration sources are pulled in:

* Per-file ``@register`` class decorators that live on editable class modules
  (e.g. :mod:`mokume.imputation.simple`,
  :mod:`mokume.normalization.feature_normalizers`,
  :mod:`mokume.quantification.median`). These are *not* imported by their
  package ``__init__`` and would otherwise stay unregistered.
* The three ``_dev_registrations`` adapter modules that expose dev-authored
  algorithms (iBAQ, the advanced imputers, RLR/LOESS/quantile normalizers)
  through thin registry adapters.
* Modules that are already imported by their package ``__init__`` (e.g. the
  quantification methods, ComBat harmonization) are intentionally *not*
  re-listed here; importing them again is a no-op but adds noise.
"""

import importlib
import logging

logger = logging.getLogger(__name__)

# Modules whose import side effect registers one or more plugins. Each entry is
# imported in isolation so that a missing optional dependency in one family
# (e.g. sklearn for the simple imputers) does not prevent the others from
# registering.
_PLUGIN_MODULES = (
    # Quantification -----------------------------------------------------
    # median.py is not imported by quantification/__init__.py.
    "mokume.quantification.median",
    # iBAQ adapter delegating to the locked ibaq.py core.
    "mokume.quantification._dev_registrations",
    # Imputation ---------------------------------------------------------
    # P2 base imputers (sklearn-backed) and the KNN imputer.
    "mokume.imputation.simple",
    "mokume.imputation.knn",
    # Dev-authored advanced imputers (bpca/qrilc/seqknn/...).
    "mokume.imputation._dev_registrations",
    # Normalization ------------------------------------------------------
    # P2 feature and sample normalizers.
    "mokume.normalization.feature_normalizers",
    "mokume.normalization.sample_normalizers",
    # Dev-authored RLR / LOESS / quantile normalizers.
    "mokume.normalization._dev_registrations",
    # Harmonization ------------------------------------------------------
    "mokume.harmonization.combat",
)


def load_all_plugins() -> None:
    """Import every plugin module, registering its plugins as a side effect.

    Each module is imported independently; if one raises ``ImportError``
    (typically a missing optional dependency) the failure is logged at debug
    level and the remaining modules are still imported. Any other exception is
    re-raised, since it signals a real bug rather than an absent extra.
    """
    for module_name in _PLUGIN_MODULES:
        try:
            importlib.import_module(module_name)
        except ImportError as exc:
            logger.debug(
                "Skipping plugin module '%s' (optional dependency missing): %s",
                module_name,
                exc,
            )


# Load on import so that ``import mokume._plugins`` is enough to populate the
# registry, mirroring the "import for side effects" idiom used by the
# individual registration modules.
load_all_plugins()
