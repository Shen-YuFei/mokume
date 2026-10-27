"""Local plot styling and matrix preparation shared by diagnostic plots."""

from contextlib import contextmanager

import numpy as np
import pandas as pd

from mokume.plotting.style import MATPLOTLIB_FONTS, MATPLOTLIB_STYLE, TEXT_COLOR


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
