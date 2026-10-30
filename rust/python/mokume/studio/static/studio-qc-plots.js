// Small SVG diagnostics drawn from the bounded post-run QC summary.
import { setTranslatedText, translate, translatedElement } from "./studio-i18n.js";

import { categoryColors, categoryOrder } from "./studio-qc-colors.js";

const SVG_NS = "http://www.w3.org/2000/svg";
const AREA = { left: 70, top: 26, width: 240, height: 174 };
const SEQUENTIAL_STEPS = 6;

function svgNode(name, attributes = {}, text = null) {
  const node = document.createElementNS(SVG_NS, name);
  Object.entries(attributes).forEach(([key, value]) => node.setAttribute(key, value));
  if (text !== null) node.textContent = text;
  return node;
}

function label(svg, text, x, y, attributes = {}) {
  const node = svgNode("text", { x, y, ...attributes }, text);
  svg.append(node);
  return node;
}

function sampleLabel(svg, sample, x, y, attributes = {}) {
  const node = label(svg, sample.length > 14 ? `${sample.slice(0, 13)}…` : sample, x, y, attributes);
  node.append(svgNode("title", {}, sample));
}

function plotFrame(title, xLabel, yLabel, { xLabelY = 308, height = 320 } = {}) {
  const figure = document.createElement("figure");
  figure.className = "qc-plot";
  figure.append(translatedElement("figcaption", "", title));
  const svg = svgNode("svg", { viewBox: `0 0 360 ${height}`, role: "img" });
  const heading = svgNode("title");
  setTranslatedText(heading, title);
  svg.append(heading, svgNode("path", {
    d: "M70 26V200H310", class: "qc-plot-axis",
  }));
  setTranslatedText(label(svg, "", 190, xLabelY, { "text-anchor": "middle" }), xLabel);
  setTranslatedText(label(svg, "", 14, 113, {
    "text-anchor": "middle", transform: "rotate(-90 14 113)",
  }), yLabel);
  figure.append(svg);
  return { figure, svg };
}

function gridLine(svg, y) {
  svg.append(svgNode("path", { d: `M${AREA.left} ${y}H${AREA.left + AREA.width}`,
    class: "qc-plot-grid" }));
}

function unavailable(plot) {
  setTranslatedText(label(plot.svg, "", 190, 110, { "text-anchor": "middle" }), "Plot unavailable");
  return plot.figure;
}

export function missingnessPlot(qc) {
  const plot = plotFrame("Missing values by sample", "Sample", "Missing values (%)");
  if (!qc.missing_by_sample?.length) return unavailable(plot);
  const step = AREA.width / qc.samples.length;
  const width = Math.min(step * 0.7, 18);
  [0, 50, 100].forEach((value) => {
    if (value) gridLine(plot.svg, 200 - value * AREA.height / 100);
    label(plot.svg, value, 64, 204 - value * AREA.height / 100, { "text-anchor": "end" });
  });
  qc.samples.forEach((sample, index) => {
    const value = qc.missing_by_sample[index];
    const height = value * AREA.height / 100;
    const x = AREA.left + step * (index + 0.5);
    const radius = Math.min(3, width / 2, height);
    // Rounded at the data end, square at the baseline.
    const bar = svgNode("path", {
      d: `M${x - width / 2} 200V${200 - height + radius}`
        + `Q${x - width / 2} ${200 - height} ${x - width / 2 + radius} ${200 - height}`
        + `H${x + width / 2 - radius}`
        + `Q${x + width / 2} ${200 - height} ${x + width / 2} ${200 - height + radius}V200Z`,
      fill: "var(--qc-condition-1)",
    });
    bar.append(svgNode("title", {}, `${sample}: ${value}%`));
    plot.svg.append(bar);
    sampleLabel(plot.svg, sample, x, 211, {
      "text-anchor": "end", transform: `rotate(-55 ${x} 211)`,
    });
  });
  return plot.figure;
}

function coordinateRange(values) {
  const minimum = Math.min(...values);
  const maximum = Math.max(...values);
  const padding = (maximum - minimum || 1) * 0.1;
  return [minimum - padding, maximum + padding];
}

function scatterTicks(svg, xRange, yRange) {
  [0, 0.5, 1].forEach((fraction) => {
    const x = AREA.left + AREA.width * fraction;
    const y = 200 - AREA.height * fraction;
    if (fraction) gridLine(svg, y);
    label(svg, (xRange[0] + fraction * (xRange[1] - xRange[0])).toFixed(1), x, 216,
      { "text-anchor": "middle" });
    label(svg, (yRange[0] + fraction * (yRange[1] - yRange[0])).toFixed(1), 64, y + 3,
      { "text-anchor": "end" });
  });
}

function project(point, [xRange, yRange]) {
  return [AREA.left + (point[0] - xRange[0]) / (xRange[1] - xRange[0]) * AREA.width,
    200 - (point[1] - yRange[0]) / (yRange[1] - yRange[0]) * AREA.height];
}

// Category marks carry their slot so a theme change can repaint them.
function categoryPaint(node, group, paint) {
  node.dataset.category = group;
  node.dataset.paint = paint;
  return node;
}

function pcaPoint(svg, point, sample, ranges, [color, group]) {
  const [x, y] = project(point, ranges);
  const marker = categoryPaint(svgNode("circle", { cx: x, cy: y, r: 3.3, fill: color,
    "fill-opacity": 0.85, stroke: "none" }), group, "fill");
  marker.append(svgNode("title", {}, `${sample}: ${point[0].toFixed(2)}, ${point[1].toFixed(2)}`));
  svg.append(marker);
  return [x, y];
}

// 95% normal data ellipse, as ggplot2 stat_ellipse(type = "norm") and the
// reports: the covariance scaled by 2 F(0.95; 2, n - 1), in closed form.
function confidenceEllipse(points) {
  const dof = points.length - 1;
  const mean = [0, 1].map((axis) =>
    points.reduce((sum, point) => sum + point[axis], 0) / points.length);
  const [xx, yy, xy] = [[0, 0], [1, 1], [0, 1]].map(([a, b]) => points.reduce((sum, point) =>
    sum + (point[a] - mean[a]) * (point[b] - mean[b]), 0) / dof);
  const radius = Math.sqrt(dof * (0.05 ** (-2 / dof) - 1));
  const centre = (xx + yy) / 2, spread = Math.hypot((xx - yy) / 2, xy);
  const angle = Math.atan2(2 * xy, xx - yy) / 2;
  const major = radius * Math.sqrt(centre + spread);
  const minor = radius * Math.sqrt(Math.max(centre - spread, 0));
  return Array.from({ length: 73 }, (_, step) => {
    const u = major * Math.cos(step * Math.PI / 36), v = minor * Math.sin(step * Math.PI / 36);
    return [mean[0] + u * Math.cos(angle) - v * Math.sin(angle),
      mean[1] + u * Math.sin(angle) + v * Math.cos(angle)];
  });
}

// Groups of three or more samples, for twelve conditions or fewer.
function groupEllipses(points, groups, count) {
  if (count > 12) return [];
  return [...Array(count).keys()].map((group) => [group,
    points.filter((_, index) => groups[index] === group)])
    .filter(([, own]) => own.length >= 3)
    .map(([group, own]) => [group, confidenceEllipse(own)]);
}

function ellipsePath(svg, outline, ranges, [color, group]) {
  const path = outline.map((point, index) =>
    `${index ? "L" : "M"}${project(point, ranges).map((value) => value.toFixed(1)).join(" ")}`);
  svg.append(categoryPaint(svgNode("path", { d: path.join(""), fill: "none", stroke: color,
    "stroke-opacity": 0.8, "stroke-width": 1 }), group, "stroke"));
}

function overlaps(a, b) {
  return a[0] < b[2] && b[0] < a[2] && a[1] < b[3] && b[1] < a[3];
}

// Outlying samples are labelled first; a label that would collide is left
// to the marker's hover title instead of stacking unreadable text.
function pcaLabels(svg, positions, samples) {
  const centre = [0, 1].map((axis) =>
    positions.reduce((sum, position) => sum + position[axis], 0) / positions.length);
  const spread = (index) => Math.hypot(positions[index][0] - centre[0],
    positions[index][1] - centre[1]);
  const boxes = positions.map(([x, y]) => [x - 5, y - 5, x + 5, y + 5]);
  [...positions.keys()].sort((a, b) => spread(b) - spread(a)).forEach((index) => {
    const [x, y] = positions[index];
    const box = [x + 6, y - 13, x + 6 + Math.min(samples[index].length, 14) * 5, y - 3];
    if (boxes.some((other) => overlaps(box, other))) return;
    boxes.push(box);
    sampleLabel(svg, samples[index], x + 6, y - 5);
  });
}

function median(values) {
  const sorted = [...values].sort((a, b) => a - b);
  return (sorted[Math.floor((sorted.length - 1) / 2)] + sorted[Math.floor(sorted.length / 2)]) / 2;
}

// Nearest clear spot in eight directions, then farther out; as in ggrepel,
// a label more than 8 px from its anchor keeps a leader.
// As in the reports, a condition label with no spot free of other labels is
// left to the legend rather than stacked on another label.
function labelBox(anchor, size, placed, dots) {
  const radii = [4, 8, 12, 16, 20, 24];
  for (let radius = 28; radius < Math.hypot(AREA.width, AREA.height); radius += 4) radii.push(radius);
  const candidates = radii.flatMap((radius) => [...Array(8).keys()].map((step) => {
    const angle = step * Math.PI / 4;
    const x = Math.max(AREA.left + size[0] / 2, Math.min(AREA.left + AREA.width - size[0] / 2,
      anchor[0] + Math.cos(angle) * (radius + size[0] / 2)));
    const y = Math.max(AREA.top + size[1] / 2, Math.min(200 - size[1] / 2,
      anchor[1] + Math.sin(angle) * (radius + size[1] / 2)));
    return { radius, box: [x - size[0] / 2, y - size[1] / 2, x + size[0] / 2, y + size[1] / 2] };
  })).filter(({ box }) => !placed.some((other) => overlaps(box, other)));
  const hits = (box) => dots.filter((dot) => overlaps(box, dot)).length;
  const near = candidates.filter(({ radius }) => radius <= 24);
  const clear = near.find(({ box }) => !hits(box));
  if (clear) {
    const gap = Math.hypot(anchor[0] - Math.max(clear.box[0], Math.min(clear.box[2], anchor[0])),
      anchor[1] - Math.max(clear.box[1], Math.min(clear.box[3], anchor[1])));
    return { box: clear.box, leader: gap > 8 };
  }
  if (near.length) {
    return { box: near.reduce((best, item) => (hits(item.box) < hits(best.box) ? item : best)).box,
      leader: true };
  }
  return candidates.length ? { box: candidates[0].box, leader: true } : null;
}

// Above twelve conditions colours repeat, so each condition is named once
// beside its sample nearest the condition median instead of naming samples.
function conditionLabels(svg, positions, groups, conditions) {
  const dots = positions.map(([x, y]) => [x - 4, y - 4, x + 4, y + 4]);
  const placed = [];
  conditions.forEach((condition, group) => {
    const own = positions.filter((_, index) => groups[index] === group);
    const centre = [0, 1].map((axis) => median(own.map((position) => position[axis])));
    const anchor = own.reduce((best, position) => (Math.hypot(position[0] - centre[0],
      position[1] - centre[1]) < Math.hypot(best[0] - centre[0], best[1] - centre[1]) ? position : best));
    const name = condition || translate("Unknown condition");
    const spot = labelBox(anchor, [Math.min(name.length, 14) * 5 + 2, 12], placed, dots);
    if (!spot) return;
    const { box, leader } = spot;
    placed.push(box);
    if (leader) {
      const end = [Math.max(box[0], Math.min(box[2], anchor[0])), Math.max(box[1], Math.min(box[3], anchor[1]))];
      svg.append(svgNode("path", { d: `M${anchor[0]} ${anchor[1]}L${end[0]} ${end[1]}`, class: "qc-plot-leader" }));
    }
    sampleLabel(svg, name, (box[0] + box[2]) / 2, box[3] - 3, { "text-anchor": "middle" });
  });
}

function conditionLegend(figure, conditions, colors) {
  const legend = document.createElement("div");
  legend.className = "qc-plot-legend";
  conditions.forEach((condition, index) => {
    const entry = document.createElement("span");
    const swatch = svgNode("svg", { viewBox: "0 0 10 10", "aria-hidden": "true" });
    swatch.append(categoryPaint(svgNode("circle", { cx: 5, cy: 5, r: 3.3, fill: colors[index],
      "fill-opacity": 0.85, stroke: "none" }), index, "fill"));
    entry.append(swatch, document.createTextNode(condition || translate("Unknown condition")));
    legend.append(entry);
  });
  figure.append(legend);
}

export function pcaPlot(qc) {
  const variance = qc.pca?.variance_percent || [];
  const plot = plotFrame("PCA by sample", "PC1", "PC2", { xLabelY: 232, height: 262 });
  const points = qc.pca?.scores || [];
  if (!points.length) return unavailable(plot);
  const conditions = [...new Set(qc.samples.map((sample) => qc.conditions?.[sample] || ""))]
    .sort(categoryOrder);
  const colors = categoryColors(conditions.length);
  const groups = qc.samples.map((sample) => conditions.indexOf(qc.conditions?.[sample] || ""));
  const ellipses = groupEllipses(points, groups, conditions.length);
  const extent = [...points, ...ellipses.flatMap(([, outline]) => outline)];
  const ranges = [0, 1].map((axis) => coordinateRange(extent.map((point) => point[axis])));
  scatterTicks(plot.svg, ...ranges);
  plot.figure.dataset.categoryCount = conditions.length;
  ellipses.forEach(([group, outline]) => ellipsePath(plot.svg, outline, ranges,
    [colors[group], group]));
  const positions = points.map((point, index) => pcaPoint(plot.svg, point, qc.samples[index],
    ranges, [colors[groups[index]], groups[index]]));
  if (conditions.length > 12) conditionLabels(plot.svg, positions, groups, conditions);
  else pcaLabels(plot.svg, positions, qc.samples);
  if (conditions.length > 1) conditionLegend(plot.figure, conditions, colors);
  label(plot.svg, translate("Explained variance: {pc1}% / {pc2}%", {
    pc1: variance[0] || 0, pc2: variance[1] || 0,
  }), 190, 252, { "text-anchor": "middle" });
  return plot.figure;
}

function heatColor(value, minimum) {
  if (value === null) return "var(--border-soft)";
  const fraction = minimum === 1 ? 1 : (value - minimum) / (1 - minimum);
  const index = Math.max(0, Math.min(SEQUENTIAL_STEPS, Math.round(fraction * SEQUENTIAL_STEPS)));
  return `var(--qc-seq-${index})`;
}

function heatmapLabels(svg, samples) {
  const width = AREA.width / samples.length;
  const height = AREA.height / samples.length;
  samples.forEach((sample, index) => {
    const x = AREA.left + (index + 0.5) * width;
    sampleLabel(svg, sample, x, 211, {
      "text-anchor": "end", transform: `rotate(-55 ${x} 211)`,
    });
    sampleLabel(svg, sample, 64, AREA.top + (index + 0.5) * height + 3,
      { "text-anchor": "end" });
  });
}

function heatmapLegend(svg, minimum) {
  const width = 220 / (SEQUENTIAL_STEPS + 1);
  for (let index = 0; index <= SEQUENTIAL_STEPS; index += 1) {
    svg.append(svgNode("rect", { x: 80 + index * width, y: 272, width, height: 8,
      fill: `var(--qc-seq-${index})` }));
  }
  label(svg, minimum.toFixed(2), 80, 294);
  label(svg, "1.00", 300, 294, { "text-anchor": "end" });
}

export function correlationPlot(qc) {
  const plot = plotFrame("Sample correlation", "Sample", "Sample");
  const matrix = qc.correlation?.matrix;
  if (!matrix?.length) return unavailable(plot);
  const minimum = matrix.reduce((lowest, row) => row.reduce((value, cell) =>
    cell === null ? value : Math.min(value, cell), lowest), 1);
  const width = AREA.width / qc.samples.length;
  const height = AREA.height / qc.samples.length;
  const gap = Math.min(1, width / 10);
  matrix.forEach((row, y) => row.forEach((value, x) => {
    const cell = svgNode("rect", { x: AREA.left + x * width, y: AREA.top + y * height,
      width: width - gap, height: height - gap, fill: heatColor(value, minimum) });
    cell.append(svgNode("title", {}, `${qc.samples[y]} / ${qc.samples[x]}: ${value ?? "—"}`));
    plot.svg.append(cell);
  }));
  heatmapLabels(plot.svg, qc.samples);
  heatmapLegend(plot.svg, minimum);
  return plot.figure;
}

// Repaint category marks with the palette of the theme now in use.
function repaintCategories() {
  document.querySelectorAll(".qc-plot[data-category-count]").forEach((figure) => {
    const colors = categoryColors(Number(figure.dataset.categoryCount));
    figure.querySelectorAll("[data-category]").forEach((node) => {
      node.setAttribute(node.dataset.paint, colors[Number(node.dataset.category)]);
    });
  });
}
document.addEventListener("studio:theme", repaintCategories);

export function qcPlots(qc) {
  const container = document.createElement("div");
  container.className = "qc-plots";
  container.append(missingnessPlot(qc), pcaPlot(qc), correlationPlot(qc));
  return container;
}
