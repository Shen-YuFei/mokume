"""Local plot styling and matrix preparation shared by diagnostic plots."""

import numpy as np
import pandas as pd

from mokume.plotting.style import AXIS_COLOR, GRID_COLOR, TEXT_COLOR


def style_axes(ax):
    """Style a single axes without modifying matplotlib defaults."""
    ax.set_axisbelow(True)
    ax.grid(color=GRID_COLOR, linestyle="-", linewidth=0.6)
    ax.tick_params(colors=TEXT_COLOR)
    ax.xaxis.label.set_color(TEXT_COLOR)
    ax.yaxis.label.set_color(TEXT_COLOR)
    ax.title.set_color(TEXT_COLOR)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_color(AXIS_COLOR)
        ax.spines[side].set_linewidth(0.6)


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
