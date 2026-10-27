"""Render data-driven Plotly cards using the shared report presentation."""

import html
import json

from mokume.plotting.profile import condition_groups
from mokume.plotting.style import PEPTIDE_COLOR, PLOTLY_TEMPLATE, TEXT_COLOR


def compact_json(value):
    """Serialize chart data compactly and escape HTML script terminators."""
    return json.dumps(value, separators=(",", ":")).replace("</", "<\\/")


# Each marker trace gets a text twin in its legend group, drawn after every
# marker so no point covers a label. Outliers keep their labels; one that would
# cover a placed label or another sample's point is hidden, and hover still
# names the sample. Visibility is redone after every redraw (zoom, resize,
# legend toggle).
EMBEDDING_LABEL_JS = """
function labelTraces(traces) {
    return traces.map(trace => {
        trace.legendgroup = trace.legendgroup || trace.name;
        // Transparent copies of the markers keep the usual offset from the point.
        return {type: 'scatter', mode: 'markers+text', x: trace.x, y: trace.y, text: trace.text,
            marker: {size: trace.marker.size, symbol: trace.marker.symbol, opacity: 0},
            legendgroup: trace.legendgroup, textposition: 'top right', textfont: {size: 8},
            cliponaxis: false, hoverinfo: 'skip', showlegend: false};
    });
}
function labelEmbedding(gd) {
    const overlaps = (a, b) => a.left < b.right && b.left < a.right
        && a.top < b.bottom && b.top < a.bottom;
    const thin = () => {
        const traces = [...gd.querySelectorAll('.scatterlayer .trace')];
        const markers = traces.filter(trace => !trace.querySelector('.textpoint'));
        const dots = markers.flatMap(trace => [...trace.querySelectorAll('.point')])
            .map(node => ({node, box: node.getBoundingClientRect()}));
        // Text twins follow their marker traces in order, point for point.
        const labels = traces.filter(trace => trace.querySelector('.textpoint'))
            .flatMap((trace, j) => {
                const own = markers[j].querySelectorAll('.point');
                return [...trace.querySelectorAll('.textpoint')].map((node, i) => {
                    const box = node.getBoundingClientRect();
                    return {node, box, own: own[i],
                        x: box.x + box.width / 2, y: box.y + box.height / 2};
                });
            });
        const cx = labels.reduce((sum, label) => sum + label.x, 0) / labels.length;
        const cy = labels.reduce((sum, label) => sum + label.y, 0) / labels.length;
        labels.sort((a, b) => Math.hypot(b.x - cx, b.y - cy) - Math.hypot(a.x - cx, a.y - cy));
        const placed = [];
        labels.forEach(({node, box, own}) => {
            const hidden = placed.some(other => overlaps(box, other))
                || dots.some(dot => dot.node !== own && overlaps(box, dot.box));
            node.style.visibility = hidden ? 'hidden' : '';
            if (!hidden) placed.push(box);
        });
    };
    gd.on('plotly_afterplot', thin);
    thin();
}
"""


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
