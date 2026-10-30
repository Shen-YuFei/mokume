"""Matrix DEqMS accepts per-cell quantification precision weights."""

import pytest

import mokume


def test_matrix_deqms_precision_weights():
    """Precision 2 floors to unit weights; other methods and shapes are rejected."""
    matrix = [
        [2.0, 4.0, 8.0, 32.0, 64.0, 128.0],
        [8.0, 16.0, 32.0, 8.0, 16.0, 32.0],
        [32.0, 30.0, 34.0, 4.0, 8.0, 16.0],
        [5.0, 6.0, 7.0, 5.5, 6.5, None],
    ]
    proteins = ["P1", "P2", "P3", "P4"]
    counts = [1, 2, 3, 4]
    plain = mokume.differential_expression(
        proteins, matrix, 3, 3, "deqms", peptide_counts=counts
    )
    unit = mokume.differential_expression(
        proteins,
        matrix,
        3,
        3,
        "deqms",
        peptide_counts=counts,
        precision=[[2.0] * 6 for _ in matrix],
    )
    assert unit == plain

    uneven = [[1.0, 8.0, 2.0, 4.0, 0.0, 3.0] for _ in matrix]
    weighted = mokume.differential_expression(
        proteins, matrix, 3, 3, "deqms", peptide_counts=counts, precision=uneven
    )
    assert [row["pvalue"] for row in weighted] != [row["pvalue"] for row in plain]

    with pytest.raises(RuntimeError, match="precision weights require DEqMS"):
        mokume.differential_expression(
            proteins, matrix, 3, 3, "limma", precision=uneven
        )
    with pytest.raises(RuntimeError, match="shape of the matrix"):
        mokume.differential_expression(
            proteins, matrix, 3, 3, "deqms", peptide_counts=counts, precision=uneven[:3]
        )
