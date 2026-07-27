"""Import-isolation tests for built-in plugins."""

import json
import os
from pathlib import Path
import subprocess
import sys


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

import mokume

print(json.dumps({
    "batch_module": mokume.is_batch_correction_available.__module__,
    "directlfq_module": mokume.is_directlfq_available.__module__,
}))
"""
    )

    assert result == {
        "batch_module": "mokume.postprocessing.batch_correction",
        "directlfq_module": "mokume.quantification.directlfq",
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
            "ibaq",
            "maxlfq",
            "median",
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

method_class = PluginRegistry.get_class("quantification", "sum")
registered = PluginRegistry.is_registered("quantification", "median")
top5 = PluginRegistry.get("quantification", "top5")
print(json.dumps({
    "class_module": method_class.__module__,
    "registered": registered,
    "topn_module": type(top5).__module__,
    "topn_name": top5.name,
}))
"""
    )

    assert result == {
        "class_module": "mokume.quantification.all_peptides",
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

class ThirdPartySum:
    pass

class EntryPoint:
    def __init__(self, name):
        self.name = name

    def load(self):
        loads.append(self.name)
        return ThirdPartySum

def entry_points(*, group):
    return [EntryPoint("sum"), EntryPoint("unrelated")]

importlib.metadata.entry_points = entry_points
PluginRegistry.reset()
method = PluginRegistry.get("quantification", "sum")
print(json.dumps({
    "loads": loads,
    "method_module": type(method).__module__,
}))
"""
    )

    assert result == {
        "loads": ["sum"],
        "method_module": "__main__",
    }
