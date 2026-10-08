// Shared state, constants, HTTP helpers, notifications, and value formatting.
import { translate } from "./studio-i18n.js";

export const ACTIVE_RUN_STATUSES = new Set(["queued", "starting", "running", "cancelling"]);

export const WORKFLOW_TEMPLATE_SCHEMA_VERSION = 1;

export const LANGUAGE_STORAGE_KEY = "mokume:language";

export const APPEARANCE_STORAGE_KEY = "mokume:appearance";

export const QUEUE_STORAGE_PREFIX = "mokume:queue:";

export const DEFAULT_PANEL_SIZES = Object.freeze({ sidebar: 235, assistant: 350, bottom: 188 });

export function savedLanguage() {
  return localStorage.getItem(LANGUAGE_STORAGE_KEY) === "zh-CN" ? "zh-CN" : "en";
}

export function savedAppearance() {
  const appearance = localStorage.getItem(APPEARANCE_STORAGE_KEY);
  return ["system", "light", "dark"].includes(appearance) ? appearance : "system";
}

export const state = {
  csrf: null,
  version: null,
  project: null,
  folderPath: null,
  folderParent: null,
  workflowTemplateMode: null,
  workflowTemplateDirectory: ".",
  workflowTemplateParent: null,
  workflowTemplateEntries: [],
  workflowTemplateSelected: null,
  commands: [],
  selectedCommand: null,
  runs: [],
  activeRunId: null,
  eventSource: null,
  eventRunId: null,
  logs: [],
  logRenderPending: false,
  artifacts: [],
  bottomTab: "runs",
  provider: null,
  dataset: null,
  projectId: null,
  agentBusy: false,
  agentAbort: null,
  pendingApproval: null,
  systemMemory: null,
  threads: { ask: crypto.randomUUID(), agent: crypto.randomUUID() },
  language: savedLanguage(),
  appearance: savedAppearance(),
  panelSizes: { ...DEFAULT_PANEL_SIZES },
  commandReview: null,
  selectedRunDetails: null,
  runComparisonDetails: null,
  qcDetails: null,
  compareRunIds: new Set(),
  batchQueue: [],
  queueBusy: false,
};

export const byId = (id) => document.getElementById(id);

export const originHeaders = () => ({
  "Content-Type": "application/json",
  "X-CSRF-Token": state.csrf,
});

export const agentHeaders = () => ({
  ...originHeaders(),
  Accept: "text/event-stream",
});

export async function api(path, options = {}) {
  const response = await fetch(path, options);
  if (!response.ok) {
    let message = `${response.status} ${response.statusText}`;
    try {
      const body = await response.json();
      message = body.detail || message;
    } catch (_) {
      // The status text is sufficient for non-JSON failures.
    }
    throw new Error(message);
  }
  return response.status === 204 ? null : response.json();
}

export function toast(message, error = false) {
  const element = byId("toast");
  element.textContent = message;
  element.classList.toggle("error", error);
  element.classList.add("visible");
  window.clearTimeout(toast.timer);
  toast.timer = window.setTimeout(() => element.classList.remove("visible"), 3200);
}

export function formatMemorySize(bytes) {
  const gibibytes = bytes / (1024 ** 3);
  return `${gibibytes.toFixed(gibibytes >= 100 ? 0 : 1)} GiB`;
}

export function titleCase(value) {
  return value.replace(/\b\w/g, (letter) => letter.toUpperCase());
}

export function formatBytes(value) {
  const bytes = Number(value);
  if (!Number.isFinite(bytes) || bytes < 0) return "—";
  const units = ["B", "KiB", "MiB", "GiB", "TiB"];
  let amount = bytes;
  let unit = 0;
  while (amount >= 1024 && unit < units.length - 1) {
    amount /= 1024;
    unit += 1;
  }
  const digits = unit === 0 || amount >= 100 ? 0 : amount >= 10 ? 1 : 2;
  return `${amount.toFixed(digits)} ${units[unit]}`;
}

export function formatDate(value) {
  if (!value) return "—";
  const date = new Date(value);
  if (Number.isNaN(date.getTime())) return value;
  return new Intl.DateTimeFormat(state.language === "zh-CN" ? "zh-CN" : "en", {
    dateStyle: "medium",
    timeStyle: "medium",
  }).format(date);
}

export function formatDuration(seconds) {
  if (!Number.isFinite(Number(seconds))) return "—";
  const total = Math.max(0, Number(seconds));
  if (total < 1) return `${Math.round(total * 1000)} ms`;
  if (total < 60) return `${total.toFixed(total < 10 ? 1 : 0)} s`;
  const minutes = Math.floor(total / 60);
  const remainder = Math.round(total % 60);
  return `${minutes} min ${remainder} s`;
}

export function formatRunDuration(run) {
  if (!run.started_at) return translate("Not started");
  if (!run.finished_at) return translate("In progress");
  const seconds = (new Date(run.finished_at) - new Date(run.started_at)) / 1000;
  return formatDuration(seconds);
}

export function shortRunId(runId) {
  return String(runId || "").slice(0, 8);
}

export function reportError(error) {
  if (error?.assistantDisplayed) return;
  toast(error.message || String(error), true);
}
