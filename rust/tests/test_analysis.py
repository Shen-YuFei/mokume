"""Migrated analysis features: imputation fallback and the QC /
workflow-comparison HTML reports. Skipped unless the ``analysis`` extra
(scipy + scikit-learn) is installed.
"""

import os

import pytest

import mokume

pytest.importorskip("scipy")
pytest.importorskip("sklearn")
pd = pytest.importorskip("pandas")


def _matrix():
    """A 12-protein x 6-sample wide matrix (proteins as the index)."""
    values = [[float((i + 1) * (j + 2)) for j in range(6)] for i in range(12)]
    return pd.DataFrame(
        values,
        columns=[f"S{k + 1}" for k in range(6)],
        index=[f"P{i + 1}" for i in range(12)],
    )


def _conditions():
    return {f"S{k + 1}": ("A" if k < 3 else "B") for k in range(6)}


def test_impute_fills_missing(tmp_path):
    frame = _matrix()
    frame.iloc[0, 0] = float("nan")
    result = mokume.impute(frame, method="knn", output=str(tmp_path / "imputed.csv"))
    assert not result.isna().any().any()


def test_qc_report_writes_html(tmp_path):
    from mokume.reports.qc_report import generate_qc_report

    protein_df = _matrix().reset_index().rename(columns={"index": "ProteinName"})
    out = tmp_path / "qc.html"
    path = generate_qc_report(
        protein_df, _conditions(), output_html=str(out), title="t"
    )
    assert os.path.getsize(path) > 0


def test_workflow_comparison_writes_html(tmp_path):
    protein_df = _matrix().reset_index().rename(columns={"index": "ProteinName"})
    conditions = _conditions()
    out = tmp_path / "cmp.html"
    path = mokume.workflow_comparison(
        [
            {
                "name": "wfA",
                "protein_df": protein_df,
                "sample_to_condition": conditions,
            },
            {
                "name": "wfB",
                "protein_df": protein_df,
                "sample_to_condition": conditions,
            },
        ],
        output=str(out),
    )
    assert os.path.getsize(path) > 0
