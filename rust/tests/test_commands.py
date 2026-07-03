"""The relocated periphery command modules import and expose ``main(argv)``."""

import importlib

import pytest

COMMAND_MODULES = [
    "visualize",
    "tissuemap",
    "de_plots",
    "interactive_report",
    "peptides2protein_qc",
    "peptides2protein_ibaq",
]


@pytest.mark.parametrize("name", COMMAND_MODULES)
def test_command_module_imports_with_main(name):
    module = importlib.import_module(f"mokume.commands.{name}")
    assert callable(getattr(module, "main", None))
