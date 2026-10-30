// Deterministic OKLab maximin extension of the six theme-compatible seeds.
// Each theme extends from its own lightness band and surface, so dark-theme
// colors are not limited to the lightness range that also suits the light theme.
function multiply(vector, matrix) {
  return matrix.map((row) => row.reduce((sum, value, index) => sum + value * vector[index], 0));
}

function linearRgb(lab) {
  const lms = multiply(lab, [
    [1, 0.3963377774, 0.2158037573],
    [1, -0.1055613458, -0.0638541728],
    [1, -0.0894841775, -1.2914855480],
  ]).map((value) => value ** 3);
  return multiply(lms, [
    [4.0767416621, -3.3077115913, 0.2309699292],
    [-1.2684380046, 2.6097574011, -0.3413193965],
    [-0.0041960863, -0.7034186147, 1.7076147010],
  ]);
}

function oklab(color) {
  const rgb = [1, 3, 5].map((offset) => parseInt(color.slice(offset, offset + 2), 16) / 255)
    .map((value) => value <= 0.04045 ? value / 12.92 : ((value + 0.055) / 1.055) ** 2.4);
  const lms = multiply(rgb, [
    [0.4122214708, 0.5363325363, 0.0514459929],
    [0.2119034982, 0.6806995451, 0.1073969566],
    [0.0883024619, 0.2817188376, 0.6299787005],
  ]).map(Math.cbrt);
  return multiply(lms, [
    [0.2104542553, 0.7936177850, -0.0040720468],
    [1.9779984951, -2.4285922050, 0.4505937099],
    [0.0259040371, 0.7827717662, -0.8086757660],
  ]);
}

function luminance(color) {
  return [1, 3, 5].map((offset) => parseInt(color.slice(offset, offset + 2), 16) / 255)
    .map((value) => value <= 0.04045 ? value / 12.92 : ((value + 0.055) / 1.055) ** 2.4)
    .reduce((sum, value, index) => sum + value * [0.2126, 0.7152, 0.0722][index], 0);
}

function contrast(a, b) {
  const [light, dark] = [luminance(a), luminance(b)].sort((x, y) => y - x);
  return (light + 0.05) / (dark + 0.05);
}

// The theme's lightness band, kept to colors with 3:1 contrast on its surface;
// chroma stays within the seeds' range.
const POOLS = new Map();
function candidates(band, surface) {
  const key = `${band} ${surface}`;
  if (POOLS.has(key)) return POOLS.get(key);
  const [low, high] = band.split(/\s+/).map(Number);
  const colors = [];
  for (let light = 0; low + light * 0.01 <= high + 1e-9; light += 1) {
    for (let chroma = 0; chroma < 8; chroma += 1) {
      for (let hue = 0; hue < 72; hue += 1) {
        const angle = hue * Math.PI / 36, c = 0.10 + chroma * 0.01;
        const lab = [low + light * 0.01, c * Math.cos(angle), c * Math.sin(angle)];
        const rgb = linearRgb(lab);
        if (rgb.some((value) => value < 0 || value > 1)) continue;
        const color = "#" + rgb.map((value) => value <= 0.0031308
          ? 12.92 * value : 1.055 * value ** (1 / 2.4) - 0.055)
          .map((value) => Math.round(value * 255).toString(16).padStart(2, "0")).join("");
        if (contrast(color, surface) >= 3) colors.push({ lab, color });
      }
    }
  }
  POOLS.set(key, colors);
  return colors;
}

function distance(a, b) {
  return a.reduce((sum, value, index) => sum + (value - b[index]) ** 2, 0);
}

// Case-insensitive category order, as the reports assign their colors.
export function categoryOrder(a, b) {
  const x = String(a).toLowerCase(), y = String(b).toLowerCase();
  return x < y ? -1 : x > y ? 1 : a < b ? -1 : a > b ? 1 : 0;
}

export function categoryColors(count) {
  const css = getComputedStyle(document.documentElement);
  const colors = Array.from({ length: 6 }, (_, index) =>
    css.getPropertyValue(`--qc-condition-${index + 1}`).trim());
  const seeds = colors.map(oklab);
  const pool = candidates(css.getPropertyValue("--qc-extension-lightness").trim(),
    css.getPropertyValue("--surface").trim());
  const distances = pool.map(({ lab }) => Math.min(...seeds.map((seed) => distance(lab, seed))));
  for (let next = colors.length; next < count; next += 1) {
    const index = distances.reduce((best, value, i) => value > distances[best] ? i : best, 0);
    colors.push(pool[index].color);
    pool.forEach(({ lab }, i) => { distances[i] = Math.min(distances[i], distance(lab, pool[index].lab)); });
    distances[index] = -Infinity;
  }
  return colors.slice(0, count);
}
