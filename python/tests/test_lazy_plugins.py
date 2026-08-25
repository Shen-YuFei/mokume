"""Import-isolation tests for built-in plugins."""

import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from mokume import _lazy
from mokume import _plugins


PYTHON_ROOT = Path(__file__).parents[1]
PROBE_PYTHON = os.environ.get("MOKUME_TEST_PYTHON", sys.executable)


def _run_import_probe(source: str) -> dict:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(PYTHON_ROOT)
    result = subprocess.run(
        [PROBE_PYTHON, "-c", source],
        check=True,
        capture_output=True,
        env=env,
        text=True,
    )
    return json.loads(result.stdout.splitlines()[-1])


def test_import_mokume_does_not_load_plugin_implementations():
    result = _run_import_probe(
        """
import json
import sys

import mokume

plugin_modules = (
    "mokume.harmonization.combat",
    "mokume.imputation._dev_registrations",
    "mokume.normalization._dev_registrations",
    "mokume.quantification._dev_registrations",
    "mokume.quantification.all_peptides",
)
print(json.dumps({
    "helpers_listed": all(
        name in dir(mokume)
        for name in ("is_batch_correction_available", "is_directlfq_available")
    ),
    "loaded_plugins": [name for name in plugin_modules if name in sys.modules],
    "scipy": "scipy" in sys.modules,
    "sklearn": "sklearn" in sys.modules,
}))
"""
    )

    assert result == {
        "helpers_listed": True,
        "loaded_plugins": [],
        "scipy": False,
        "sklearn": False,
    }


def test_top_level_helpers_keep_their_original_function_identity():
    result = _run_import_probe(
        """
import json
import sys

import mokume

loaded_before = {
    "batch": "mokume.postprocessing.batch_correction" in sys.modules,
    "directlfq": "mokume.quantification.directlfq" in sys.modules,
}
batch_helper = mokume.is_batch_correction_available
directlfq_helper = mokume.is_directlfq_available
from mokume.postprocessing.batch_correction import (
    is_batch_correction_available as original_batch_helper,
)
from mokume.quantification.directlfq import (
    is_directlfq_available as original_directlfq_helper,
)
print(json.dumps({
    "batch_identity": batch_helper is original_batch_helper,
    "directlfq_identity": directlfq_helper is original_directlfq_helper,
    "loaded_after": {
        "batch": "mokume.postprocessing.batch_correction" in sys.modules,
        "directlfq": "mokume.quantification.directlfq" in sys.modules,
    },
    "loaded_before": loaded_before,
}))
"""
    )

    assert result == {
        "batch_identity": True,
        "directlfq_identity": True,
        "loaded_after": {"batch": True, "directlfq": True},
        "loaded_before": {"batch": False, "directlfq": False},
    }


def test_sum_registry_lookup_loads_only_its_implementation():
    result = _run_import_probe(
        """
import json
import sys

from mokume.core.registry import PluginRegistry

method = PluginRegistry.get("quantification", "sum")
unrelated_prefixes = (
    "mokume.harmonization.",
    "mokume.imputation.",
    "mokume.normalization.",
)
unrelated_modules = (
    "mokume.quantification._dev_registrations",
    "mokume.quantification.directlfq",
    "mokume.quantification.maxlfq",
    "mokume.quantification.median",
)
print(json.dumps({
    "method_module": type(method).__module__,
    "scipy": "scipy" in sys.modules,
    "sklearn": "sklearn" in sys.modules,
    "unrelated": sorted(
        name for name in sys.modules
        if name in unrelated_modules or name.startswith(unrelated_prefixes)
    ),
}))
"""
    )

    assert result == {
        "method_module": "mokume.quantification.all_peptides",
        "scipy": False,
        "sklearn": False,
        "unrelated": [],
    }


def test_quantification_exports_resolve_without_eager_imports():
    result = _run_import_probe(
        """
import json
import sys

import mokume.quantification as quantification

loaded_before = "mokume.quantification.all_peptides" in sys.modules
method_class = quantification.AllPeptidesQuantification
print(json.dumps({
    "class_module": method_class.__module__,
    "listed": "AllPeptidesQuantification" in dir(quantification),
    "loaded_before": loaded_before,
}))
"""
    )

    assert result == {
        "class_module": "mokume.quantification.all_peptides",
        "listed": True,
        "loaded_before": False,
    }


def test_quantification_availability_and_reset_are_preserved():
    result = _run_import_probe(
        """
import json

from mokume.core.registry import PluginRegistry

names = PluginRegistry.available("quantification")
PluginRegistry.reset()
method = PluginRegistry.get("quantification", "sum")
print(json.dumps({
    "method_module": type(method).__module__,
    "names": names,
}))
"""
    )

    assert result == {
        "method_module": "mokume.quantification.all_peptides",
        "names": [
            "directlfq",
            "maxlfq",
            "median",
            "peptide_count",
            "pibaq",
            "ratio",
            "spectral_count",
            "sum",
            "tmt_abundance",
            "tmt_reporter",
            "topn",
        ],
    }


def test_registry_class_registration_and_topn_paths_are_preserved():
    result = _run_import_probe(
        """
import json

from mokume.core.registry import PluginRegistry

@PluginRegistry.register("quantification", "sum")
class CustomSum:
    pass

first = PluginRegistry.get("quantification", "sum")
names = PluginRegistry.available("quantification")
second = PluginRegistry.get("quantification", "sum")
PluginRegistry.reset()

method_class = PluginRegistry.get_class("quantification", "sum")
registered = PluginRegistry.is_registered("quantification", "median")
top5 = PluginRegistry.get("quantification", "top5")
print(json.dumps({
    "class_module": method_class.__module__,
    "custom_first_module": type(first).__module__,
    "custom_listed": "sum" in names,
    "custom_second_module": type(second).__module__,
    "registered": registered,
    "topn_module": type(top5).__module__,
    "topn_name": top5.name,
}))
"""
    )

    assert result == {
        "class_module": "mokume.quantification.all_peptides",
        "custom_first_module": "__main__",
        "custom_listed": True,
        "custom_second_module": "__main__",
        "registered": True,
        "topn_module": "mokume.quantification.topn",
        "topn_name": "Top5",
    }


def test_registry_loads_only_the_matching_third_party_entry_point():
    result = _run_import_probe(
        """
import importlib.metadata
import json

from mokume.core.registry import PluginRegistry

loads = []
scans = 0

class ThirdPartySum:
    pass

class EntryPoint:
    def __init__(self, name):
        self.name = name

    def load(self):
        loads.append(self.name)
        return ThirdPartySum

def entry_points(*, group):
    global scans
    scans += 1
    return [EntryPoint("sum"), EntryPoint("unrelated")]

importlib.metadata.entry_points = entry_points
PluginRegistry.reset()
method = PluginRegistry.get("quantification", "sum")
other_method = PluginRegistry.get("quantification", "unrelated")
PluginRegistry.reset()
reset_method = PluginRegistry.get("quantification", "sum")
print(json.dumps({
    "loads": loads,
    "method_module": type(method).__module__,
    "other_method_module": type(other_method).__module__,
    "reset_method_module": type(reset_method).__module__,
    "scans": scans,
}))
"""
    )

    assert result == {
        "loads": ["sum", "unrelated", "sum"],
        "method_module": "__main__",
        "other_method_module": "__main__",
        "reset_method_module": "__main__",
        "scans": 2,
    }


def test_plugin_loader_distinguishes_missing_owner_from_dependency(monkeypatch):
    module_name = "mokume.missing_plugin"
    debug_calls = []
    error_calls = []
    monkeypatch.setattr(
        _plugins.logger, "debug", lambda *args: debug_calls.append(args)
    )
    monkeypatch.setattr(
        _plugins.logger, "error", lambda *args: error_calls.append(args)
    )

    def raise_missing_owner(_):
        raise ModuleNotFoundError("missing owner", name=module_name)

    monkeypatch.setattr(_plugins.importlib, "import_module", raise_missing_owner)
    assert _plugins._load_module(module_name) is None
    assert len(error_calls) == 1
    assert not debug_calls

    error_calls.clear()

    def raise_missing_dependency(_):
        raise ModuleNotFoundError("missing dependency", name="optional_backend")

    monkeypatch.setattr(_plugins.importlib, "import_module", raise_missing_dependency)
    assert _plugins._load_module(module_name) is None
    assert len(debug_calls) == 1
    assert not error_calls


def test_lazy_export_preserves_errors_from_target_modules(monkeypatch):
    def raise_target_error(_):
        raise KeyError("inside target")

    monkeypatch.setattr(_lazy.importlib, "import_module", raise_target_error)
    get_attribute, _ = _lazy.module_api(
        {"known": ("mokume.target", "value")},
        {},
        "mokume.probe",
    )

    with pytest.raises(KeyError, match="inside target"):
        get_attribute("known")
    with pytest.raises(
        AttributeError,
        match="module 'mokume.probe' has no attribute 'unknown'",
    ):
        get_attribute("unknown")
