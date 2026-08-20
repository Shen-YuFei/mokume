"""Tests for :class:`mokume.pipeline.config.QuantificationConfig` validation.

``ratio_fraction_merge`` is read by the Python and Rust implementations, which
disagree on case, so the config layer folds it. See the ``__post_init__``
comment for why that has to happen here rather than in either reader.
"""

from __future__ import annotations

import pytest

from mokume.pipeline.config import QuantificationConfig


class TestRatioFractionMerge:
    def test_default_is_mean(self) -> None:
        assert QuantificationConfig().ratio_fraction_merge == "mean"

    @pytest.mark.parametrize(
        ("given", "expected"),
        [
            ("mean", "mean"),
            ("max", "max"),
            ("Max", "max"),
            ("MEAN", "mean"),
        ],
    )
    def test_case_is_folded(self, given: str, expected: str) -> None:
        config = QuantificationConfig(ratio_fraction_merge=given)
        assert config.ratio_fraction_merge == expected

    @pytest.mark.parametrize("bad", ["median", "sum", "", "maximum"])
    def test_rejects_unknown_methods(self, bad: str) -> None:
        with pytest.raises(ValueError, match="must be 'mean' or 'max'"):
            QuantificationConfig(ratio_fraction_merge=bad)

    def test_error_names_the_value_the_user_wrote(self) -> None:
        with pytest.raises(ValueError, match="'Median'"):
            QuantificationConfig(ratio_fraction_merge="Median")
