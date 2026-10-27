"""Local plot styling and matrix preparation shared by diagnostic plots."""

from contextlib import contextmanager
from textwrap import fill

import numpy as np
import pandas as pd

from mokume.plotting.style import (
    CONDITION_MARKERS,
    condition_style,
    MATPLOTLIB_FONTS,
    MATPLOTLIB_STYLE,
    SAMPLE_POINT_COLOR,
    TEXT_COLOR,
)


@contextmanager
def static_style():
    """Scope static-plot defaults to creation and export, restoring caller settings.

    Prefer an installed Latin sans-serif, followed by the caller's installed
    font families for glyphs such as Chinese sample names. Import matplotlib
    only when plotting, since matrix helpers also serve HTML-only reports.
    """
    from matplotlib import font_manager, rc_context, rcParams

    installed = {
        font.name.casefold(): font.name for font in font_manager.fontManager.ttflist
    }
    primary = next(
        (name for name in MATPLOTLIB_FONTS if name.casefold() in installed), None
    )
    configured = []
    for family in rcParams["font.family"]:
        configured.extend(rcParams.get(f"font.{family}", [family]))
    families = [primary] if primary else []
    families.extend(
        installed[name.casefold()]
        for name in configured
        if name.casefold() in installed
    )
    with rc_context({**MATPLOTLIB_STYLE, "font.family": list(dict.fromkeys(families))}):
        yield


def style_axes(ax):
    """Style axes locally; static_style scopes defaults without global changes."""
    ax.set_axisbelow(True)
    ax.tick_params(colors=TEXT_COLOR, direction="out")
    ax.xaxis.label.set_color(TEXT_COLOR)
    ax.yaxis.label.set_color(TEXT_COLOR)
    ax.title.set_color(TEXT_COLOR)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_color(TEXT_COLOR)
        ax.spines[side].set_linewidth(0.8)


def condition_markers(conditions):
    """Use circles through three conditions, then the shared shape codes."""
    return {
        condition: (
            CONDITION_MARKERS[condition_style(index)[1]] if len(conditions) > 3 else "o"
        )
        for index, condition in enumerate(conditions)
    }


def nearest_anchor(points):
    """Choose the observed point nearest the coordinate-wise group median."""
    points = np.asarray(points)
    return points[np.argmin(np.linalg.norm(points - np.median(points, axis=0), axis=1))]


def _overlap_count(box, boxes):
    """Count overlapping display-space rectangles."""
    rectangles = np.asarray(boxes).reshape(-1, 4)
    return np.count_nonzero(
        (box[0] < rectangles[:, 2])
        & (box[2] > rectangles[:, 0])
        & (box[1] < rectangles[:, 3])
        & (box[3] > rectangles[:, 1])
    )


def _label_candidates(size, center, bounds, radii):
    """Try eight directions at increasing distances within the label bounds."""
    half = np.asarray(size) / 2
    for radius in radii:
        for angle in np.arange(8) * np.pi / 4:
            direction = np.array([np.cos(angle), np.sin(angle)])
            point = center + direction * (radius + half)
            point = np.clip(point, bounds[:2] + half, bounds[2:] - half)
            yield point, np.r_[point - half, point + half]


def _label_position(size, center, placed, layout):
    """Prefer nearby clear positions; far or compromised labels get a leader."""
    bounds, points, scale = layout["bounds"], layout["points"], layout["scale"]
    candidates = list(
        _label_candidates(size, center, bounds, np.arange(4, 25, 4) * scale)
    )
    available = [
        (point, box) for point, box in candidates if not _overlap_count(box, placed)
    ]
    clear = [
        (point, box) for point, box in available if not _overlap_count(box, points)
    ]
    if clear:
        point, box = clear[0]
        # As in ggrepel, a label pushed past the nearest rings keeps a leader.
        gap = np.linalg.norm(center - np.clip(center, box[:2], box[2:]))
        return point, box, gap > 8 * scale
    if available:
        point, box = min(available, key=lambda item: _overlap_count(item[1], points))
        return point, box, True
    # Prefer a longer leader to overlap, keeping a nearby fallback for dense plots.
    radii = np.arange(28 * scale, np.linalg.norm(bounds[2:] - bounds[:2]), 4 * scale)
    point, box = next(
        (
            (point, box)
            for point, box in _label_candidates(size, center, bounds, radii)
            if not _overlap_count(box, placed)
        ),
        min(candidates, key=lambda item: np.linalg.norm(item[0] - center)),
    )
    return point, box, True


def _place_one_label(ax, label, anchor, placed, layout):
    """Measure and position one label in the already settled figure layout."""
    text = ax.annotate(
        label,
        anchor,
        xytext=(0, 0),
        textcoords="offset points",
        ha="center",
        va="center",
        color=TEXT_COLOR,
        fontsize=layout["fontsize"],
        arrowprops={"arrowstyle": "-", "color": SAMPLE_POINT_COLOR, "lw": 0.5},
        zorder=5,
    )
    extent = text.get_window_extent(ax.figure.canvas.get_renderer())
    center = ax.transData.transform(anchor)
    point, box, leader = _label_position(
        (extent.width + 2 * layout["scale"], extent.height + 2 * layout["scale"]),
        center,
        placed,
        layout,
    )
    text.xyann = (point - center) / layout["scale"]
    text.arrow_patch.set_visible(leader)
    placed.append(box)


def place_labels(ax, anchors, points=None, fontsize=7):
    """Keep neutral labels inside their axes, avoiding overlaps when space allows."""
    ax.figure.canvas.draw()
    scale = ax.figure.dpi / 72
    if points is None:
        points = [
            xy for collection in ax.collections for xy in collection.get_offsets()
        ]
    centers = ax.transData.transform(np.asarray(points).reshape(-1, 2))
    bounds = ax.get_window_extent().extents + np.array([2, 2, -2, -2]) * scale
    layout = {
        "bounds": bounds,
        "scale": scale,
        "fontsize": fontsize,
        "points": [np.r_[center - 4 * scale, center + 4 * scale] for center in centers],
    }
    placed = [
        text.get_window_extent().extents
        for text in [
            ax.xaxis.label,
            ax.yaxis.label,
            *ax.get_xticklabels(),
            *ax.get_yticklabels(),
        ]
    ]
    for label in sorted(anchors, key=lambda name: -len(str(name))):
        _place_one_label(ax, label, anchors[label], placed, layout)


def label_categories(ax, frame, coordinates, group_column):
    """Label each category at its median-nearest sample."""
    anchors = {
        fill(str(group), 28): nearest_anchor(part[list(coordinates)].to_numpy())
        for group, part in frame.groupby(group_column, sort=True)
    }
    place_labels(ax, anchors)


def expand_category_axes(ax):
    """Give crowded embeddings a page-width canvas and 6.4 by 4.5 inch axes."""
    ax.figure.set_layout_engine(None)
    ax.figure.set_size_inches(7.2, 5.5)
    ax.set_position([0.65 / 7.2, 0.55 / 5.5, 6.4 / 7.2, 4.5 / 5.5])


def _legend_columns(fig, handles, labels, title):
    """Use the largest column count whose rendered width fits 7.2 inches."""
    renderer = fig.canvas.get_renderer()
    labels = [fill(str(label), 38) for label in labels]
    options = {
        "title": title,
        "loc": "lower center",
        "bbox_to_anchor": (0.5, 0.02),
        "frameon": False,
        "labelcolor": TEXT_COLOR,
        "columnspacing": 1,
        "handletextpad": 0.4,
    }
    for columns in range(len(labels), 1, -1):
        legend = fig.legend(handles, labels, ncol=columns, **options)
        if legend.get_window_extent(renderer).width / fig.dpi <= 6.9:
            return legend
        legend.remove()
    return fig.legend(handles, labels, ncol=1, **options)


def bottom_legend(ax, handles=None, title=None):
    """Grow the canvas below a legend while retaining all axes' physical sizes."""
    fig = ax.figure
    fig.canvas.draw()
    if handles is None:
        handles, labels = ax.get_legend_handles_labels()
    else:
        labels = [handle.get_label() for handle in handles]
    if ax.get_legend() is not None:
        ax.get_legend().remove()
    legend = _legend_columns(fig, handles, labels, title)
    bounds = legend.get_window_extent(fig.canvas.get_renderer())
    width, height = fig.get_size_inches()
    positions = [
        axis.get_position().bounds * np.array([width, height, width, height])
        for axis in fig.axes
    ]
    extra = bounds.height / fig.dpi + 0.35
    new_width = max(width, bounds.width / fig.dpi + 0.3)
    fig.set_layout_engine(None)
    fig.set_size_inches(new_width, height + extra)
    for axis, position in zip(fig.axes, positions):
        position[:2] += ((new_width - width) / 2, extra)
        axis.set_position(
            position / [new_width, height + extra, new_width, height + extra]
        )


def panel_label(ax, letter):
    """Keep a lowercase bold panel letter separate from the regular title."""
    ax.text(
        -0.025,
        1.025,
        letter.lower(),
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontweight="bold",
        color=TEXT_COLOR,
    )


def read_frame(source, sep=","):
    """Accept an in-memory table or a table path."""
    return source if isinstance(source, pd.DataFrame) else pd.read_csv(source, sep=sep)


def log2_matrix(frame, is_log2=False):
    """Index a wide protein table and mask nonpositive linear intensities."""
    values = frame.set_index(frame.columns[0])
    return values if is_log2 else np.log2(values.where(values > 0))


def ordered_samples(samples, conditions):
    """Order samples by condition and then by sample name."""
    return sorted(samples, key=lambda sample: (conditions.get(sample, ""), sample))


def finite_list(values):
    """Return compact, standards-compliant numeric JSON values."""
    return [round(float(value), 2) if np.isfinite(value) else None for value in values]
