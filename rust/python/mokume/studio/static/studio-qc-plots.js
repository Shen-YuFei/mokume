// Small SVG diagnostics drawn from the bounded post-run QC summary.
import { setTranslatedText, translate, translatedElement } from "./studio-i18n.js";

const SVG_NS = "http://www.w3.org/2000/svg";
const AREA = { left: 70, top: 26, width: 240, height: 174 };
const CONDITION_SLOTS = 6;
const SEQUENTIAL_STEPS = 6;
// Shapes join colour once a scatter shows more than three conditions,
// so every condition keeps a distinct colour-shape pair.
const SHAPES = [
  (x, y) => `M${x - 4} ${y}a4 4 0 1 0 8 0a4 4 0 1 0 -8 0Z`,
  (x, y) => `M${x - 3.5} ${y - 3.5}h7v7h-7Z`,
  (x, y) => `M${x} ${y - 4.6}L${x + 4.2} ${y + 3}H${x - 4.2}Z`,
  (x, y) => `M${x} ${y - 5}L${x + 5} ${y}L${x} ${y + 5}L${x - 5} ${y}Z`,
  (x, y) => `M${x} ${y + 4.6}L${x + 4.2} ${y - 3}H${x - 4.2}Z`,
];

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

function conditionStyle(index, count) {
  return {
    fill: `var(--qc-condition-${index % CONDITION_SLOTS + 1})`,
    shape: SHAPES[count > 3 ? index % SHAPES.length : 0],
  };
}

function pcaPoint(svg, point, sample, ranges, style) {
  const [xRange, yRange] = ranges;
  const x = AREA.left + (point[0] - xRange[0]) / (xRange[1] - xRange[0]) * AREA.width;
  const y = 200 - (point[1] - yRange[0]) / (yRange[1] - yRange[0]) * AREA.height;
  const marker = svgNode("path", { d: style.shape(x, y), fill: style.fill });
  marker.append(svgNode("title", {}, `${sample}: ${point[0].toFixed(2)}, ${point[1].toFixed(2)}`));
  svg.append(marker);
  return [x, y];
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

function conditionLegend(figure, conditions) {
  const legend = document.createElement("div");
  legend.className = "qc-plot-legend";
  conditions.forEach((condition, index) => {
    const entry = document.createElement("span");
    const style = conditionStyle(index, conditions.length);
    const swatch = svgNode("svg", { viewBox: "0 0 10 10", "aria-hidden": "true" });
    swatch.append(svgNode("path", { d: style.shape(5, 5), fill: style.fill }));
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
  const xRange = coordinateRange(points.map((point) => point[0]));
  const yRange = coordinateRange(points.map((point) => point[1]));
  scatterTicks(plot.svg, xRange, yRange);
  const conditions = [...new Set(qc.samples.map((sample) => qc.conditions?.[sample] || ""))];
  const positions = points.map((point, index) => pcaPoint(plot.svg, point, qc.samples[index],
    [xRange, yRange], conditionStyle(
      conditions.indexOf(qc.conditions?.[qc.samples[index]] || ""), conditions.length)));
  pcaLabels(plot.svg, positions, qc.samples);
  if (conditions.length > 1) conditionLegend(plot.figure, conditions);
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

export function qcPlots(qc) {
  const container = document.createElement("div");
  container.className = "qc-plots";
  container.append(missingnessPlot(qc), pcaPlot(qc), correlationPlot(qc));
  return container;
}
