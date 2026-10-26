// Small SVG diagnostics drawn from the bounded post-run QC summary.
import { setTranslatedText, translate, translatedElement } from "./studio-i18n.js";

const SVG_NS = "http://www.w3.org/2000/svg";
const AREA = { left: 70, top: 26, width: 240, height: 174 };

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

function plotFrame(title, xLabel, yLabel) {
  const figure = document.createElement("figure");
  figure.className = "qc-plot";
  figure.append(translatedElement("figcaption", "", title));
  const svg = svgNode("svg", { viewBox: "0 0 360 320", role: "img" });
  const heading = svgNode("title");
  setTranslatedText(heading, title);
  svg.append(heading, svgNode("path", {
    d: "M70 26V200H310", class: "qc-plot-axis",
  }));
  setTranslatedText(label(svg, "", 190, 308, { "text-anchor": "middle" }), xLabel);
  setTranslatedText(label(svg, "", 14, 113, {
    "text-anchor": "middle", transform: "rotate(-90 14 113)",
  }), yLabel);
  figure.append(svg);
  return { figure, svg };
}

function unavailable(plot) {
  setTranslatedText(label(plot.svg, "", 190, 110, { "text-anchor": "middle" }), "Plot unavailable");
  return plot.figure;
}

export function missingnessPlot(qc) {
  const plot = plotFrame("Missing values by sample", "Sample", "Missing values (%)");
  if (!qc.missing_by_sample?.length) return unavailable(plot);
  const step = AREA.width / qc.samples.length;
  [0, 50, 100].forEach((value) => label(plot.svg, value, 64,
    204 - value * AREA.height / 100, { "text-anchor": "end" }));
  qc.samples.forEach((sample, index) => {
    const value = qc.missing_by_sample[index];
    const height = value * AREA.height / 100;
    const x = AREA.left + step * (index + 0.5);
    const bar = svgNode("rect", {
      x: x - step * 0.4, y: 200 - height, width: step * 0.8, height,
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
    label(svg, (xRange[0] + fraction * (xRange[1] - xRange[0])).toFixed(1), x, 216,
      { "text-anchor": "middle" });
    label(svg, (yRange[0] + fraction * (yRange[1] - yRange[0])).toFixed(1), 64, y + 3,
      { "text-anchor": "end" });
  });
}

function pcaPoint(svg, point, sample, ranges, colorIndex) {
  const [xRange, yRange] = ranges;
  const x = AREA.left + (point[0] - xRange[0]) / (xRange[1] - xRange[0]) * AREA.width;
  const y = 200 - (point[1] - yRange[0]) / (yRange[1] - yRange[0]) * AREA.height;
  const marker = svgNode("circle", { cx: x, cy: y, r: 4,
    fill: `var(--qc-condition-${colorIndex + 1})` });
  marker.append(svgNode("title", {}, `${sample}: ${point[0].toFixed(2)}, ${point[1].toFixed(2)}`));
  svg.append(marker);
  sampleLabel(svg, sample, x + 6, y - 5);
}

function conditionLegend(figure, conditions) {
  const legend = document.createElement("div");
  legend.className = "qc-plot-legend";
  conditions.forEach((condition, index) => {
    const entry = document.createElement("span");
    const swatch = document.createElement("i");
    swatch.className = `qc-condition-${index % 8 + 1}`;
    entry.append(swatch, document.createTextNode(condition || translate("Unknown condition")));
    legend.append(entry);
  });
  figure.append(legend);
}

export function pcaPlot(qc) {
  const variance = qc.pca?.variance_percent || [];
  const plot = plotFrame("PCA by sample", "PC1", "PC2");
  const points = qc.pca?.scores || [];
  if (!points.length) return unavailable(plot);
  const xRange = coordinateRange(points.map((point) => point[0]));
  const yRange = coordinateRange(points.map((point) => point[1]));
  scatterTicks(plot.svg, xRange, yRange);
  const conditions = [...new Set(qc.samples.map((sample) => qc.conditions?.[sample] || ""))];
  points.forEach((point, index) => pcaPoint(plot.svg, point, qc.samples[index],
    [xRange, yRange], conditions.indexOf(qc.conditions?.[qc.samples[index]] || "") % 8));
  if (conditions.some(Boolean)) conditionLegend(plot.figure, conditions);
  label(plot.svg, translate("Explained variance: {pc1}% / {pc2}%", {
    pc1: variance[0] || 0, pc2: variance[1] || 0,
  }), 190, 245, { "text-anchor": "middle" });
  return plot.figure;
}

function heatColor(value, minimum) {
  if (value === null) return "var(--border-soft)";
  const fraction = minimum === 1 ? 1 : (value - minimum) / (1 - minimum);
  const index = Math.max(0, Math.min(10, Math.round(fraction * 10)));
  // High correlation is red, matching the reversed RdBu scale of the HTML QC report.
  return `var(--qc-rdbu-${10 - index})`;
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
  for (let index = 0; index <= 10; index += 1) {
    svg.append(svgNode("rect", { x: 80 + index * 20, y: 272, width: 20, height: 8,
      fill: `var(--qc-rdbu-${10 - index})` }));
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
  matrix.forEach((row, y) => row.forEach((value, x) => {
    const cell = svgNode("rect", { x: AREA.left + x * width, y: AREA.top + y * height,
      width, height, fill: heatColor(value, minimum) });
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
