from __future__ import annotations

import pytest

from mokume._command_flags import flags_for


def test_rejects_unknown_wrapper_keyword() -> None:
    with pytest.raises(TypeError, match="unexpected keyword"):
        flags_for("features2proteins", {"not_a_real_option": 1})


def test_csv_and_repeatable_options_use_their_real_shapes() -> None:
    assert flags_for(
        "features2proteins",
        {
            "de_contrasts": ["A vs B", "C vs D"],
            "irs_reference_sample": ["Pool A", "Pool B"],
        },
    ) == [
        "--de-contrasts",
        "A vs B,C vs D",
        "--irs-reference-sample",
        "Pool A",
        "--irs-reference-sample",
        "Pool B",
    ]


def test_reverse_boolean_can_disable_default_contaminant_removal() -> None:
    assert flags_for("features2proteins", {"remove_contaminants": False}) == [
        "--keep-contaminants"
    ]
    assert flags_for("features2proteins", {"remove_contaminants": True}) == []


def test_tissuemap_plural_destination_maps_to_repeatable_singular_flag() -> None:
    assert flags_for("tissuemap", {"tmt_datasets": ["PXD1", "PXD2"]}) == [
        "--tmt-dataset",
        "PXD1",
        "--tmt-dataset",
        "PXD2",
    ]


def test_nonboolean_value_option_rejects_bool() -> None:
    with pytest.raises(TypeError, match="expects a value"):
        flags_for("features2proteins", {"min_unique": False})


def test_python_name_maps_to_real_families_flag() -> None:
    assert flags_for("peptides2protein", {"families_yaml": "families.yaml"}) == [
        "--families",
        "families.yaml",
    ]
