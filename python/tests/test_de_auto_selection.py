"""Tests for the data-aware ``de.method = "auto"`` rule in the pipeline."""

import pytest

from mokume.pipeline.stages import _auto_select_de_method


def _sample_map(samples_per_condition: dict[str, int]) -> dict[str, str]:
    """Build sample metadata with the requested group sizes."""
    sample_to_condition = {}
    for condition, n in samples_per_condition.items():
        for i in range(n):
            name = f"{condition}_{i}"
            sample_to_condition[name] = condition
    return sample_to_condition


@pytest.mark.parametrize("quant", ["directlfq", "maxlfq", "pibaq"])
def test_tiny_groups_get_the_permutation_test(quant):
    """Below 3 replicates the moderated-t methods are underpowered."""
    assert _auto_select_de_method(_sample_map({"A": 2, "B": 3}), quant) == "rots"


def test_directlfq_keeps_deqms_when_groups_are_large_enough():
    assert _auto_select_de_method(_sample_map({"A": 3, "B": 3}), "directlfq") == "deqms"


def test_other_quantifications_keep_limrots():
    assert _auto_select_de_method(_sample_map({"A": 4, "B": 4}), "maxlfq") == "limrots"


def test_the_smallest_group_decides_not_the_average():
    """A 2-vs-20 design is still a 2-replicate problem."""
    assert _auto_select_de_method(_sample_map({"A": 2, "B": 20}), "maxlfq") == "rots"
