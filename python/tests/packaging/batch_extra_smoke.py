"""Exercise Python batch correction after installing an advertised extra."""

from importlib.metadata import requires

import numpy as np
import pandas as pd
from inmoose.pycombat import pycombat_norm
from mokume.harmonization.correction import (
    apply_batch_correction as apply_harmonization_correction,
)
from mokume.postprocessing.batch_correction import (
    apply_batch_correction,
    is_batch_correction_available,
)
from packaging.requirements import Requirement


def _requirements_for_extra(extra: str) -> set[str]:
    """Return normalized dependency names selected by an extra."""
    selected = set()
    for requirement_text in requires("mokume") or []:
        requirement = Requirement(requirement_text)
        if requirement.marker and requirement.marker.evaluate({"extra": extra}):
            selected.add(requirement.name.lower())
    return selected


def _check_extra_metadata() -> None:
    """Verify that each ComBat consumer retains its own backend."""
    batch_requirements = _requirements_for_extra("batch-correction")
    assert "inmoose" in batch_requirements
    assert "combat" not in batch_requirements

    alias_requirements = _requirements_for_extra("inmoose")
    assert "inmoose" in alias_requirements
    assert "combat" not in alias_requirements

    tissuemap_requirements = _requirements_for_extra("tissuemap")
    assert "combat" in tissuemap_requirements
    assert "inmoose" not in tissuemap_requirements

    all_requirements = _requirements_for_extra("all")
    assert {"combat", "inmoose"} <= all_requirements


def main() -> None:
    """Run a small finite matrix through the backend imported by Mokume."""
    _check_extra_metadata()
    assert callable(pycombat_norm)
    assert is_batch_correction_available()

    data = pd.DataFrame(
        [
            [10.0, 11.0, 9.0, 14.0, 15.0, 13.0],
            [20.0, 19.0, 21.0, 25.0, 24.0, 26.0],
            [5.0, 6.0, 5.5, 8.0, 9.0, 8.5],
            [30.0, 31.0, 29.0, 34.0, 35.0, 33.0],
        ],
        columns=[f"sample-{index}" for index in range(6)],
    )
    batches = [0, 0, 0, 1, 1, 1]
    before_gap = abs(data.iloc[:, :3].mean().mean() - data.iloc[:, 3:].mean().mean())

    for apply_correction in (
        apply_harmonization_correction,
        apply_batch_correction,
    ):
        corrected = apply_correction(data, batches)

        assert corrected.shape == data.shape
        assert np.isfinite(corrected.to_numpy()).all()

        after_gap = abs(
            corrected.iloc[:, :3].mean().mean() - corrected.iloc[:, 3:].mean().mean()
        )
        assert after_gap < before_gap


if __name__ == "__main__":
    main()
