"""The pure-Python DEqMS reproduces the R oracle that pins the Rust kernel.

The fixture is the Rust crate's (DEqMS 1.28.0 under R 4.5.3): count-moderated
two-group tests with positive varying counts, including a small n and repeated
counts.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mokume.analysis._deqms_loess import loess_predict
from mokume.analysis.deqms import _spectra_count_ebayes, run_deqms
from mokume.analysis.limma import _contrasts_fit, _ebayes, _lm_fit_with_na

FIXTURE = (
    Path(__file__).resolve().parents[2]
    / "rust/crates/mokume-stats/tests/data/deqms_official.json"
)
CASES = json.loads(FIXTURE.read_text())["cases"]


@pytest.mark.parametrize("case", CASES)
def test_loess_and_prior_df_match_deqms(case):
    matrix = np.array(case["matrix"])
    design = np.zeros((matrix.shape[1], 2))
    design[: case["n_a"], 0] = 1.0
    design[case["n_a"] :, 1] = 1.0
    proteins = [str(index) for index in range(len(matrix))]
    fit = _ebayes(
        _contrasts_fit(
            _lm_fit_with_na(matrix, design, proteins, ["A", "B"]),
            np.array([[1.0], [-1.0]]),
        )
    )
    counts = np.array(case["counts"])
    log_var = np.log(fit.sigma**2)

    fitted = loess_predict(log_var, np.log2(counts), fit.df_residual > 0)

    np.testing.assert_allclose(fitted, case["loess_fitted"], rtol=0, atol=1e-10)
    assert _spectra_count_ebayes(fit, counts)[4] == case["prior_df"]


@pytest.mark.parametrize("case", CASES)
def test_deqms_matches_r(case):
    matrix = pd.DataFrame(case["matrix"]).rename(columns=str)
    samples = list(matrix.columns)

    result = run_deqms(
        matrix,
        samples[: case["n_a"]],
        samples[case["n_a"] :],
        ("A", "B"),
        peptide_counts=pd.Series(case["counts"]),
    ).set_index("ProteinName")

    for column, expected in (
        ("log2FC", "log2fc"),
        ("sca_t", "statistic"),
        ("pvalue", "pvalue"),
        ("adj_pvalue", "bh"),
    ):
        np.testing.assert_allclose(
            result.loc[matrix.index, column], case[expected], rtol=1e-8, atol=1e-10
        )
