"""Periphery rendering tests: the folded ``mokume.plotting`` draws real files.

Skipped unless the ``plotting`` extra (matplotlib + pandas) is installed.
"""

import glob
import os

import pytest

import mokume

pytest.importorskip("matplotlib")
pytest.importorskip("pandas")

DATA = os.path.join(os.path.dirname(__file__), "data")


def test_peptides2protein_qc_writes_pdf(tmp_path):
    pdf = tmp_path / "qc.pdf"
    mokume.peptides2protein_qc(
        protein_table=os.path.join(DATA, "proteins_qc.tsv"),
        qc_report=str(pdf),
        plot_column="PiBAQ",
    )
    assert pdf.exists()
    assert pdf.stat().st_size > 0


def test_de_plots_writes_volcano_png(tmp_path):
    plot_dir = tmp_path / "plots"
    mokume.de_plots(
        [
            "--protein-matrix",
            os.path.join(DATA, "proteins_matrix.csv"),
            "--plot-dir",
            str(plot_dir),
            "--volcano",
            "--contrast",
            "c1",
            "CondA",
            "CondB",
            os.path.join(DATA, "de.csv"),
        ]
    )
    assert glob.glob(str(plot_dir / "*.png"))
