// Appearance, panel visibility and resizing, and the system memory indicator.
import {
  APPEARANCE_STORAGE_KEY,
  DEFAULT_PANEL_SIZES,
  api,
  byId,
  formatMemorySize,
  state,
} from "./studio-core.js";
import { setTranslatedAttribute, translate } from "./studio-i18n.js";
import { closeMenus } from "./studio-menus.js";

export const SYSTEM_DARK_QUERY = window.matchMedia("(prefers-color-scheme: dark)");

export const PANEL_SIZE_LIMITS = Object.freeze({
  sidebar: Object.freeze({ min: 180 }),
  assistant: Object.freeze({ min: 280 }),
  bottom: Object.freeze({ min: 120, max: 560 }),
});

export const SIDE_PANEL_MAX_VIEWPORT_RATIO = 0.5;

export const SIDE_PANEL_COLLAPSE_RATIO = 0.5;

export const SYSTEM_MEMORY_REFRESH_MS = 5000;

export const MIN_WORKFLOW_WIDTH = 420;

export const HIGH_MEMORY_PERCENT = 85;

export function applyAppearance(appearance, persist = true) {
  state.appearance = ["light", "dark"].includes(appearance) ? appearance : "system";
  const theme = state.appearance === "system"
    ? (SYSTEM_DARK_QUERY.matches ? "dark" : "light")
    : state.appearance;
  document.documentElement.dataset.theme = theme;
  document.querySelectorAll("[data-appearance]").forEach((button) => {
    button.setAttribute("aria-checked", String(button.dataset.appearance === state.appearance));
  });
  if (persist) localStorage.setItem(APPEARANCE_STORAGE_KEY, state.appearance);
}

export function renderSystemMemory() {
  const indicator = byId("system-memory");
  const memory = state.systemMemory;
  if (!memory) {
    indicator.style.setProperty("--memory-percent", "0");
    indicator.setAttribute("aria-valuenow", "0");
    setTranslatedAttribute(indicator, "aria-label", "System memory");
    indicator.title = translate("System memory");
    indicator.classList.add("unavailable");
    indicator.classList.remove("high");
    byId("system-memory-label").textContent = "–";
    return;
  }
  const percent = Math.min(100, Math.max(0, Number(memory.percent) || 0));
  const label = translate("System memory: {percent}% ({used} / {total})", {
    percent: percent.toFixed(1),
    used: formatMemorySize(memory.used_bytes),
    total: formatMemorySize(memory.total_bytes),
  });
  indicator.style.setProperty("--memory-percent", String(percent));
  indicator.setAttribute("aria-valuenow", percent.toFixed(1));
  indicator.setAttribute("aria-label", label);
  indicator.title = label;
  indicator.classList.remove("unavailable");
  indicator.classList.toggle("high", percent >= HIGH_MEMORY_PERCENT);
  byId("system-memory-label").textContent = `${Math.round(percent)}%`;
}

export async function refreshSystemMemory() {
  try {
    state.systemMemory = (await api("/api/system")).memory;
  } catch (_) {
    state.systemMemory = null;
  }
  renderSystemMemory();
}

export function toggleBottom() {
  byId("bottom-panel").classList.toggle("collapsed");
  resizeWorkspace();
}

export function toggleSidebar() {
  const workspace = document.querySelector(".workspace");
  if (window.matchMedia("(max-width: 680px)").matches) {
    if (!workspace.classList.contains("sidebar-mobile-visible")) {
      workspace.classList.remove("assistant-mobile-visible");
    }
    workspace.classList.toggle("sidebar-mobile-visible");
    applyPanelSizes();
    return;
  }
  workspace.classList.toggle("sidebar-hidden");
  applyPanelSizes();
}

export function toggleAssistant() {
  const workspace = document.querySelector(".workspace");
  if (window.matchMedia("(max-width: 1050px)").matches) {
    if (!workspace.classList.contains("assistant-mobile-visible")) {
      workspace.classList.remove("sidebar-mobile-visible");
    }
    workspace.classList.toggle("assistant-mobile-visible");
    applyPanelSizes();
    return;
  }
  workspace.classList.toggle("assistant-hidden");
  applyPanelSizes();
}

export function setSidePanelCollapsed(panel, collapsed) {
  const workspace = document.querySelector(".workspace");
  if (panel === "sidebar") {
    if (window.matchMedia("(max-width: 680px)").matches) {
      workspace.classList.toggle("sidebar-mobile-visible", !collapsed);
    } else {
      workspace.classList.toggle("sidebar-hidden", collapsed);
    }
  }
  if (panel === "assistant") {
    if (window.matchMedia("(max-width: 1050px)").matches) {
      workspace.classList.toggle("assistant-mobile-visible", !collapsed);
    } else {
      workspace.classList.toggle("assistant-hidden", collapsed);
    }
  }
  applyPanelSizes();
}

export function panelSizeBounds(panel) {
  const limits = PANEL_SIZE_LIMITS[panel];
  if (panel === "bottom") {
    const headerHeight = document.querySelector(".app-header").getBoundingClientRect().height;
    const available = window.innerHeight - headerHeight - 240;
    return { min: limits.min, max: Math.max(limits.min, Math.min(limits.max, available)) };
  }

  const workspace = document.querySelector(".workspace");
  const viewportMax = Math.floor(window.innerWidth * SIDE_PANEL_MAX_VIEWPORT_RATIO);
  if (panel === "assistant" && window.matchMedia("(max-width: 1050px)").matches) {
    return { min: Math.min(limits.min, viewportMax), max: viewportMax };
  }
  if (panel === "sidebar" && window.matchMedia("(max-width: 680px)").matches) {
    return { min: Math.min(limits.min, viewportMax), max: viewportMax };
  }

  let occupied = 0;
  if (panel === "sidebar" && !window.matchMedia("(max-width: 1050px)").matches
      && !workspace.classList.contains("assistant-hidden")) {
    occupied = state.panelSizes.assistant;
  }
  if (panel === "assistant" && !workspace.classList.contains("sidebar-hidden")) {
    occupied = state.panelSizes.sidebar;
  }
  const available = workspace.clientWidth - MIN_WORKFLOW_WIDTH - occupied;
  const max = Math.max(0, Math.min(viewportMax, available));
  return { min: Math.min(limits.min, max), max };
}

export function setPanelSize(panel, requestedSize) {
  const bounds = panelSizeBounds(panel);
  const size = Math.round(Math.min(Math.max(requestedSize, bounds.min), bounds.max));
  state.panelSizes[panel] = size;
  if (panel === "bottom") {
    if (!byId("bottom-panel").classList.contains("collapsed")) {
      document.documentElement.style.setProperty("--bottom-height", `${size}px`);
    }
  } else {
    document.documentElement.style.setProperty(`--${panel}-width`, `${size}px`);
  }
  const handle = document.querySelector(`[data-resize-panel="${panel}"]`);
  handle?.setAttribute("aria-valuemin", String(Math.round(bounds.min)));
  handle?.setAttribute("aria-valuemax", String(Math.round(bounds.max)));
  handle?.setAttribute("aria-valuenow", String(size));
}

export function applyPanelSizes() {
  setPanelSize("assistant", state.panelSizes.assistant);
  setPanelSize("sidebar", state.panelSizes.sidebar);
  setPanelSize("assistant", state.panelSizes.assistant);
  setPanelSize("bottom", state.panelSizes.bottom);
  resizeWorkspace();
}

export function bindPanelResizers() {
  document.querySelectorAll("[data-resize-panel]").forEach((handle) => {
    const panel = handle.dataset.resizePanel;
    handle.addEventListener("pointerdown", (event) => {
      if (event.button !== 0) return;
      event.preventDefault();
      closeMenus();
      if (panel === "bottom") {
        byId("bottom-panel").classList.remove("collapsed");
        resizeWorkspace();
      }
      const startX = event.clientX;
      const startY = event.clientY;
      const startSize = state.panelSizes[panel];
      const cursorClass = panel === "bottom" ? "resize-rows" : "resize-columns";
      const dragSurface = document.querySelector(".workspace");
      const collapseThreshold = panel === "bottom"
        ? null
        : panelSizeBounds(panel).min * SIDE_PANEL_COLLAPSE_RATIO;
      let collapsedDuringDrag = false;
      handle.setPointerCapture(event.pointerId);
      handle.classList.add("active");
      document.body.classList.add("resizing", cursorClass);

      const move = (moveEvent) => {
        let size = startSize;
        if (panel === "sidebar") size += moveEvent.clientX - startX;
        if (panel === "assistant") size -= moveEvent.clientX - startX;
        if (panel === "bottom") size -= moveEvent.clientY - startY;
        const shouldCollapse = collapseThreshold !== null && size <= collapseThreshold;
        if (shouldCollapse !== collapsedDuringDrag) {
          collapsedDuringDrag = shouldCollapse;
          if (shouldCollapse) {
            state.panelSizes[panel] = startSize;
            dragSurface.setPointerCapture(moveEvent.pointerId);
          }
          setSidePanelCollapsed(panel, shouldCollapse);
        }
        if (shouldCollapse) return;
        setPanelSize(panel, size);
      };
      const finish = (finishEvent) => {
        if (handle.hasPointerCapture(finishEvent.pointerId)) {
          handle.releasePointerCapture(finishEvent.pointerId);
        }
        if (dragSurface.hasPointerCapture(finishEvent.pointerId)) {
          dragSurface.releasePointerCapture(finishEvent.pointerId);
        }
        window.removeEventListener("pointermove", move);
        window.removeEventListener("pointerup", finish);
        window.removeEventListener("pointercancel", finish);
        handle.classList.remove("active");
        document.body.classList.remove("resizing", cursorClass);
      };
      window.addEventListener("pointermove", move);
      window.addEventListener("pointerup", finish);
      window.addEventListener("pointercancel", finish);
    });
    handle.addEventListener("dblclick", () => {
      if (panel === "bottom") byId("bottom-panel").classList.remove("collapsed");
      setPanelSize(panel, DEFAULT_PANEL_SIZES[panel]);
      resizeWorkspace();
    });
    handle.addEventListener("keydown", (event) => {
      const bounds = panelSizeBounds(panel);
      let size = state.panelSizes[panel];
      if (event.key === "Home") size = bounds.min;
      else if (event.key === "End") size = bounds.max;
      else if (panel === "sidebar" && event.key === "ArrowLeft") size -= 10;
      else if (panel === "sidebar" && event.key === "ArrowRight") size += 10;
      else if (panel === "assistant" && event.key === "ArrowLeft") size += 10;
      else if (panel === "assistant" && event.key === "ArrowRight") size -= 10;
      else if (panel === "bottom" && event.key === "ArrowUp") size += 10;
      else if (panel === "bottom" && event.key === "ArrowDown") size -= 10;
      else return;
      event.preventDefault();
      if (panel === "bottom") byId("bottom-panel").classList.remove("collapsed");
      setPanelSize(panel, size);
      resizeWorkspace();
    });
  });
  window.addEventListener("resize", applyPanelSizes);
}

export function resizeWorkspace() {
  const collapsed = byId("bottom-panel").classList.contains("collapsed");
  const height = collapsed ? 39 : state.panelSizes.bottom;
  document.documentElement.style.setProperty("--bottom-height", `${height}px`);
}
