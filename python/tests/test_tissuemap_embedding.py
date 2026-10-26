"""Regression tests for optional TissueMap embedding dependencies."""

from __future__ import annotations

import logging
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import numpy as np
import pandas as pd
import pytest

from mokume.tissuemap.embedding import _impute_for_pca
from mokume.plotting.style import CONDITION_COLORS
from mokume.tissuemap.plotting.atlas import (
    _draw_dendrogram_panel,
    _group_styles,
    _label_groups,
    _name_cluster,
    plot_slide_atlas_dendrogram,
)
from mokume.tissuemap.plotting import markers as marker_plots
from mokume.tissuemap.plotting.markers import _protein_label
from mokume.tissuemap.plotting.markers import _significant_markers
from mokume.tissuemap.tissue_specificity import _compute_ts_vectorized_mad


PYTHON_ROOT = Path(__file__).parents[1]


def test_tsne_module_import_does_not_load_optional_umap(tmp_path):
    """A broken optional UMAP install must not block the independent t-SNE path."""
    (tmp_path / "umap.py").write_text(
        "raise RuntimeError('UMAP should be loaded only when requested')\n",
        encoding="utf-8",
    )
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join((str(tmp_path), str(PYTHON_ROOT)))
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; import mokume.tissuemap.embedding; "
            "assert 'umap' not in sys.modules",
        ],
        capture_output=True,
        check=False,
        env=env,
        text=True,
    )

    assert result.returncode == 0, result.stderr


def test_embedding_mindet_uses_each_samples_protein_distribution():
    """Embedding adapts its sample-by-protein matrix to the imputer contract."""
    matrix = np.array(
        [
            [1.0, 10.0, np.nan],
            [2.0, np.nan, 100.0],
            [np.nan, 30.0, 1000.0],
        ]
    )

    result = _impute_for_pca(matrix, "mindet", n_neighbors=2)

    assert result[0, 2] == pytest.approx(np.quantile([1.0, 10.0], 0.01))
    assert result[1, 1] == pytest.approx(np.quantile([2.0, 100.0], 0.01))
    assert result[2, 0] == pytest.approx(np.quantile([30.0, 1000.0], 0.01))


def test_pure_mad_pi_uses_observed_values_as_denominator():
    """Missing values are not counted as population-model outliers."""
    matrix = np.array([[0.0], [1.0], [2.0], [np.nan], [np.nan]])

    _scores, params = _compute_ts_vectorized_mad(
        matrix,
        np.array(["A", "A", "A", "B", "B"]),
        ["A", "B"],
        sigma_floor=0.01,
    )

    assert params[0].pi == pytest.approx(1.0)


def test_marker_display_requires_effect_and_adjusted_significance():
    """Top-ranked non-significant proteins are not labelled as markers."""
    result = pd.DataFrame(
        {
            "names": ["kept", "weak_fdr", "weak_effect"],
            "pvals_adj": [0.01, 0.1, 0.01],
            "logfoldchanges": [0.6, 1.0, 0.2],
        }
    )

    assert _significant_markers(result)["names"].tolist() == ["kept"]


def test_atlas_dendrogram_skips_single_tissue(tmp_path, caplog):
    """A single-tissue dataset has no meaningful dendrogram to plot."""
    adata = SimpleNamespace(
        obsm={"X_tsne": np.array([[0.0, 0.0], [1.0, 1.0]])},
        obs=pd.DataFrame({"tissue": ["liver", "liver"]}),
    )

    with caplog.at_level(logging.WARNING):
        plot_slide_atlas_dendrogram(adata, tmp_path)

    assert "requires at least two tissues" in caplog.text
    assert not any(tmp_path.iterdir())


@pytest.mark.parametrize(
    ("members", "used", "expected"),
    [
        (["brain, cortex", "brain, cerebellum", "lung"], set(), "Brain"),
        (["brain, cortex", "lung"], set(), "Group 3"),
        (["brain, cortex", "brain, cerebellum"], {"Brain"}, "Group 3"),
    ],
)
def test_atlas_names_require_a_unique_strict_organ_majority(members, used, expected):
    """A minority family or a previously used organ name cannot name a group."""
    assert _name_cluster(members, 3, used) == expected


def test_atlas_group_codes_use_circles_then_unique_color_shape_pairs():
    """Atlas and dendrogram group styles must remain unambiguous beyond six."""
    small = _group_styles({f"Group {i}": [] for i in range(6)})
    large = _group_styles({f"Group {i}": [] for i in range(10)})

    assert list(small.values()) == [(color, "o") for color in CONDITION_COLORS]
    assert len(set(large.values())) == 10
    assert len({marker for _, marker in large.values()}) > 1


@pytest.mark.parametrize(
    ("identifier", "label"),
    [
        ("sp|P12345|ABC_HUMAN", "ABC"),
        ("tr|A0A123|ABC_SUB_HUMAN", "ABC_SUB"),
        ("db|ACC|ENTRY", "db|ACC|ENTRY"),
        ("P12345", "P12345"),
        ("ABC_HUMAN", "ABC_HUMAN"),
    ],
)
def test_marker_labels_shorten_only_uniprot_entries(identifier, label):
    """Plain identifiers retain their original labels."""
    assert _protein_label(identifier) == label


def test_marker_tsne_shares_a_protein_relative_scale_without_changing_values(
    tmp_path, monkeypatch
):
    """Differently scaled proteins share one colorbar; missing samples stay missing."""
    matrix = np.array([[1, 100], [2, 200], [3, 300], [4, 400], [np.nan, np.nan]])
    adata = ad.AnnData(
        matrix.copy(),
        obs=pd.DataFrame({"tissue": ["a"] * 5}, index=list("abcde")),
        var=pd.DataFrame(index=["p1", "p2"]),
    )
    adata.obsm["X_tsne"] = np.arange(10).reshape(5, 2)
    adata.uns["tissue_markers"] = {}
    monkeypatch.setattr(
        marker_plots, "_select_showcase_markers", lambda *_: [("p1", "a"), ("p2", "a")]
    )
    figures = []
    monkeypatch.setattr(
        Figure, "savefig", lambda fig, *_args, **_kwargs: figures.append(fig)
    )

    marker_plots.plot_marker_tsne(adata, ["a"], tmp_path, save_pdf=False)

    fig = figures[0]
    visible = [ax for ax in fig.axes if ax.get_visible()]
    assert len(visible) == 3
    for ax in visible[:2]:
        detected = ax.collections[1]
        np.testing.assert_allclose(
            detected.get_array(), [0, 0.97 / 2.94, 1.97 / 2.94, 1]
        )
        assert detected.get_clim() == (0, 1)
    assert "within each protein" in visible[-1].get_ylabel()
    np.testing.assert_array_equal(adata.X, matrix)


def test_atlas_group_labels_stay_inside_axes_without_overlapping():
    """Labels near an embedding boundary must not cross the title or each other."""
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set(xlim=(0, 1), ylim=(0, 1))
    _label_groups(ax, {f"Group {i}": np.array([0.95, 0.95]) for i in range(10)})
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bounds = ax.get_window_extent()
    boxes = [text.get_window_extent(renderer) for text in ax.texts]

    assert all(bounds.contains(*box.min) and bounds.contains(*box.max) for box in boxes)
    assert not any(a.overlaps(b) for i, a in enumerate(boxes) for b in boxes[i + 1 :])
    plt.close(fig)


def test_dendrogram_group_markers_do_not_cover_leaf_labels():
    """Marker spacing in points must not depend on the dendrogram width."""
    fig, ax = plt.subplots(figsize=(10, 6))
    _draw_dendrogram_panel(
        ax,
        np.array([[0, 1, 1.0, 2]]),
        ["brain", "skin"],
        {"brain": "Group 1", "skin": "Group 1"},
        {"Group 1": (CONDITION_COLORS[0], "o")},
    )
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for label, marker in zip(ax.get_yticklabels(), ax.collections[1:]):
        center = marker.get_offset_transform().transform(marker.get_offsets())[0]
        radius = (np.sqrt(marker.get_sizes()[0]) + marker.get_linewidths()[0]) / 2
        assert label.get_window_extent(renderer).x1 < center[0] - radius * fig.dpi / 72
    plt.close(fig)
