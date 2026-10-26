"""Render data-driven Plotly cards using the shared report presentation."""

import html
import json

from mokume.plotting.profile import condition_groups
from mokume.plotting.style import PEPTIDE_COLOR, PLOTLY_TEMPLATE, TEXT_COLOR


def compact_json(value):
    """Serialize chart data compactly and escape HTML script terminators."""
    return json.dumps(value, separators=(",", ":")).replace("</", "<\\/")


def condition_guides(samples, conditions):
    """Return Plotly separators and top labels for condition-grouped samples."""
    boundaries, labels = condition_groups(samples, conditions)
    line = {"color": PEPTIDE_COLOR, "width": 1}
    return {
        "shapes": [
            {
                "type": "line",
                "x0": x,
                "x1": x,
                "y0": 0,
                "y1": 1,
                "yref": "paper",
                "line": line,
            }
            for x in boundaries
        ],
        "annotations": [
            {
                "x": x,
                "y": 1,
                "yref": "paper",
                "yanchor": "bottom",
                "text": text,
                "showarrow": False,
                "font": {"color": TEXT_COLOR},
            }
            for x, text in labels
        ],
    }


def plot_cards(plots):
    """Return card markup and a short renderer for standalone reports."""
    cards = "".join(
        f'<div class="card{" wide" if plot.get("wide") else ""}">'
        f"<h3>{html.escape(plot['title'])}</h3>"
        f'<div id="{plot["id"]}" style="height:{plot.get("height", 420)}px"></div>'
        f'<p class="info-text">{html.escape(plot.get("caption", ""))}</p></div>'
        for plot in plots
    )
    specs = [
        {
            "id": plot["id"],
            "data": plot["data"],
            "layout": {
                "template": PLOTLY_TEMPLATE,
                "margin": {"t": 35, "b": 100, "l": 70, "r": 30},
                **plot["layout"],
            },
        }
        for plot in plots
    ]
    payload = compact_json(specs)
    script = f"""
(function() {{
    const specs = {payload};
    requestAnimationFrame(() => requestAnimationFrame(() => {{
        specs.forEach(spec => Plotly.newPlot(spec.id, spec.data, spec.layout,
            {{responsive: true, displaylogo: false, displayModeBar: 'hover'}}));
    }}));
}})();
"""
    return f'<div class="plot-grid">{cards}</div>', script
