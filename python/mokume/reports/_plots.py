"""Render data-driven Plotly cards using the shared report presentation."""

import html
import json

from mokume.plotting.profile import condition_groups
from mokume.plotting.style import (
    PEPTIDE_COLOR,
    PLOTLY_TEMPLATE,
    SAMPLE_POINT_COLOR,
    TEXT_COLOR,
)


def compact_json(value):
    """Serialize chart data compactly and escape HTML script terminators."""
    return json.dumps(value, separators=(",", ":")).replace("</", "<\\/")


LABEL_LAYOUT_JS = """
const labelOverlap = (a, b) => a.left < b.right && b.left < a.right
    && a.top < b.bottom && b.top < a.bottom;
function labelCandidates(center, size, bounds, radii) {
    const candidates = [];
    radii.forEach(radius => {
        for (let direction = 0; direction < 8; direction++) {
            const angle = direction * Math.PI / 4;
            const x = Math.max(bounds.left + size.width / 2,
                Math.min(bounds.right - size.width / 2,
                    center.x + Math.cos(angle) * (radius + size.width / 2)));
            const y = Math.max(bounds.top + size.height / 2,
                Math.min(bounds.bottom - size.height / 2,
                    center.y + Math.sin(angle) * (radius + size.height / 2)));
            candidates.push({x, y, left: x - size.width / 2, right: x + size.width / 2,
                top: y - size.height / 2, bottom: y + size.height / 2});
        }
    });
    return candidates;
}
function nearbyLabel(center, size, placed, dots, bounds) {
    const candidates = labelCandidates(center, size, bounds, [4, 8, 12, 16, 20, 24].map(r => r * 4 / 3));
    const available = candidates.filter(box => !placed.some(other => labelOverlap(box, other)));
    const clear = available.find(box => !dots.some(dot => labelOverlap(box, dot)));
    if (clear) {
        const gap = Math.hypot(center.x - Math.max(clear.left, Math.min(clear.right, center.x)),
            center.y - Math.max(clear.top, Math.min(clear.bottom, center.y)));
        return {...clear, leader: gap > 8 * 4 / 3};
    }
    if (available.length) {
        const hits = box => dots.filter(dot => labelOverlap(box, dot)).length;
        available.sort((a, b) => hits(a) - hits(b));
        return {...available[0], leader: true};
    }
    const radii = [];
    for (let radius = 38; radius < Math.hypot(bounds.right - bounds.left,
            bounds.bottom - bounds.top); radius += 6) radii.push(radius);
    return {...labelCandidates(center, size, bounds, radii)
        .find(box => !placed.some(other => labelOverlap(box, other))), leader: true};
}
function plotBounds(gd) {
    const rect = gd.getBoundingClientRect(), size = gd._fullLayout._size;
    return {left: rect.left + size.l + 3, right: rect.left + size.l + size.w - 3,
        top: rect.top + size.t + 3, bottom: rect.top + size.t + size.h - 3};
}
function labelLeader(gd, center, box) {
    if (!box.leader) return;
    const svg = gd.querySelector('svg.main-svg'), rect = svg.getBoundingClientRect();
    const x = Math.max(box.left, Math.min(box.right, center.x));
    const y = Math.max(box.top, Math.min(box.bottom, center.y));
    const edges = [{x: box.left, y}, {x: box.right, y},
        {x, y: box.top}, {x, y: box.bottom}];
    const distance = point => Math.hypot(point.x - center.x, point.y - center.y);
    const end = edges.reduce((a, b) => distance(a) <= distance(b) ? a : b);
    const line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
    Object.entries({x1: center.x - rect.left, y1: center.y - rect.top,
        x2: end.x - rect.left, y2: end.y - rect.top, stroke: '__LEADER_COLOR__',
        'stroke-width': 2 / 3, class: 'mokume-label-leader'}).forEach(([key, value]) => line.setAttribute(key, value));
    svg.appendChild(line);
}
function categoryLayout(target, traces, height) {
    const count = new Set(traces.map(trace => trace.meta.condition)).size;
    if (count <= 12) return {};
    const extra = Math.ceil(count / 3) * 20;
    document.getElementById(target).style.height = (height + extra) + 'px';
    return {height: height + extra,
        legend: {orientation: 'h', y: -0.2, x: 0.5, xanchor: 'center', font: {size: 10}},
        margin: {t: 55, b: 95 + extra, l: 60, r: 45}};
}
""".replace("__LEADER_COLOR__", SAMPLE_POINT_COLOR)


# Each marker trace gets a text twin in its legend group, drawn after every
# marker so no point covers a label. Outliers keep their labels; one that would
# cover a placed label or another sample's point is hidden, and hover still
# names the sample. Visibility is redone after every redraw (zoom, resize,
# legend toggle).
EMBEDDING_LABEL_JS = (
    LABEL_LAYOUT_JS
    + """
function categoryTraces(traces) {
    const groups = [...new Set(traces.map(trace => trace.meta.condition))];
    const median = values => {
        const sorted = [...values].sort((a, b) => a - b), mid = Math.floor(sorted.length / 2);
        return (sorted[mid] + sorted[Math.floor((sorted.length - 1) / 2)]) / 2;
    };
    return groups.map(condition => {
        const own = traces.filter(trace => trace.meta.condition === condition);
        const points = own.flatMap(trace => trace.x.map((x, i) => [x, trace.y[i]]));
        const cx = median(points.map(point => point[0])), cy = median(points.map(point => point[1]));
        points.sort((a, b) => Math.hypot(a[0] - cx, a[1] - cy) - Math.hypot(b[0] - cx, b[1] - cy));
        return {type: 'scatter', mode: 'text', x: [points[0][0]], y: [points[0][1]],
            text: [condition], textposition: 'middle center', textfont: {size: 9},
            legendgroup: condition, showlegend: false, hoverinfo: 'skip', cliponaxis: false,
            meta: {categoryLabel: true}};
    });
}
function placeCategoryLabels(gd) {
    gd.querySelectorAll('.mokume-label-leader').forEach(node => node.remove());
    const bounds = plotBounds(gd), placed = [];
    const dots = [...gd.querySelectorAll('.scatterlayer .point')].map(node => node.getBoundingClientRect());
    gd._fullData.filter(trace => trace.meta?.categoryLabel && trace.visible === true).forEach(trace => {
        const node = gd.querySelector('.trace' + trace.uid + ' .textpoint');
        if (!node) return;
        node.setAttribute('transform', 'translate(0,0)');
        const rect = node.getBoundingClientRect();
        const center = {x: (rect.left + rect.right) / 2, y: (rect.top + rect.bottom) / 2};
        const box = nearbyLabel(center, {width: rect.width + 4, height: rect.height + 4}, placed, dots, bounds);
        node.setAttribute('transform', `translate(${box.x - center.x},${box.y - center.y})`);
        placed.push(box);
        labelLeader(gd, center, box);
    });
}
function labelTraces(traces) {
    if (new Set(traces.map(trace => trace.meta.condition)).size > 12) return categoryTraces(traces);
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
        if (gd.data.some(trace => trace.meta?.categoryLabel)) return placeCategoryLabels(gd);
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
)


VOLCANO_LABEL_JS = (
    LABEL_LAYOUT_JS
    + """
function volcanoLabelChanges(gd) {
    const rect = gd.getBoundingClientRect(), layout = gd._fullLayout;
    const nodes = [...gd.querySelectorAll('.annotation')];
    const centers = layout.annotations.map(annotation => ({
        x: rect.left + layout.xaxis._offset + layout.xaxis.l2p(annotation.x),
        y: rect.top + layout.yaxis._offset + layout.yaxis.l2p(annotation.y)}));
    const dots = centers.map(point => ({left: point.x - 8, right: point.x + 8,
        top: point.y - 8, bottom: point.y + 8}));
    const placed = [], changes = {};
    layout.annotations.forEach((annotation, i) => {
        const text = nodes.find(node => Number(node.getAttribute('data-index')) === i)
            ?.querySelector('.annotation-text');
        if (!text) return;
        const size = text.getBoundingClientRect(), center = centers[i];
        const box = nearbyLabel(center, {width: size.width + 6, height: size.height + 6},
            placed, dots, plotBounds(gd));
        placed.push(box);
        const ax = box.x - center.x, ay = box.y - center.y;
        if (Math.abs(annotation.ax - ax) > 0.5 || Math.abs(annotation.ay - ay) > 0.5) {
            changes[`annotations[${i}].ax`] = ax;
            changes[`annotations[${i}].ay`] = ay;
        }
    });
    return changes;
}
function labelVolcano(gd) {
    let updating = false;
    const update = () => {
        if (updating) return;
        const changes = volcanoLabelChanges(gd);
        if (!Object.keys(changes).length) return;
        updating = true;
        Plotly.relayout(gd, changes).finally(() => {updating = false;});
    };
    gd.on('plotly_afterplot', update);
    update();
}
"""
)


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
