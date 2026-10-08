// Run history, the batch queue, live run events, artifacts, and the bottom panel.
import {
  ACTIVE_RUN_STATUSES,
  QUEUE_STORAGE_PREFIX,
  api,
  byId,
  formatDate,
  formatRunDuration,
  originHeaders,
  reportError,
  shortRunId,
  state,
  toast,
} from "./studio-core.js";
import { translate, translatedElement } from "./studio-i18n.js";
import { updateActionStates } from "./studio-menus.js";
import { resizeWorkspace } from "./studio-layout.js";
import {
  loadLatestQc,
  openRunComparison,
  openRunDetails,
  refreshOpenRunDetails,
  renderQcPanel,
  showCompletedRunQc,
} from "./studio-run-details.js";

export const MAX_LOG_LINES = 1000;

export function queueStorageKey(projectId = state.projectId) {
  return projectId ? `${QUEUE_STORAGE_PREFIX}${projectId}` : null;
}

export function loadBatchQueue(projectId) {
  const key = queueStorageKey(projectId);
  if (!key) return [];
  try {
    const saved = JSON.parse(localStorage.getItem(key) || "[]");
    return Array.isArray(saved) ? saved.filter(validQueueItem) : [];
  } catch (_) {
    return [];
  }
}

export function validQueueItem(item) {
  return item && typeof item.id === "string" && Array.isArray(item.argv)
    && ["pending", "queued", "starting", "running", "cancelling", "succeeded", "failed", "cancelled", "interrupted"].includes(item.status);
}

export function saveBatchQueue() {
  const key = queueStorageKey();
  if (key) localStorage.setItem(key, JSON.stringify(state.batchQueue));
  renderTabCounts();
}

export function queueOutputConflict(review) {
  const outputs = new Set(review.outputs.map((item) => item.path));
  return state.batchQueue.some((item) => ["pending", "queued", "starting", "running", "cancelling"].includes(item.status)
    && (item.outputs || []).some((path) => outputs.has(path)));
}

export async function addReviewedToQueue() {
  const review = state.commandReview;
  if (!review?.valid) return;
  if (queueOutputConflict(review)) throw new Error(translate("A queued analysis already writes to this output"));
  state.batchQueue.push({
    id: crypto.randomUUID(),
    argv: review.argv,
    command: review.command,
    workflow: review.display_name,
    outputs: review.outputs.map((item) => item.path),
    status: "pending",
    runId: null,
    error: null,
    createdAt: new Date().toISOString(),
  });
  saveBatchQueue();
  byId("command-review-dialog").close();
  showBottomTab("queue");
  toast(translate("Added to queue"));
  await processBatchQueue();
}

export function syncBatchQueue() {
  const runs = new Map(state.runs.map((run) => [run.id, run]));
  let changed = false;
  state.batchQueue.forEach((item) => {
    if (!item.runId) return;
    const run = runs.get(item.runId);
    if (run && item.status !== run.status) {
      item.status = run.status;
      item.error = run.error;
      changed = true;
    }
  });
  if (changed) saveBatchQueue();
}

export async function processBatchQueue() {
  if (!state.project || state.queueBusy || state.activeRunId) return;
  const item = state.batchQueue.find((candidate) => candidate.status === "pending");
  if (!item) return;
  state.queueBusy = true;
  try {
    const run = await api("/api/runs", {
      method: "POST",
      headers: originHeaders(),
      body: JSON.stringify({ argv: item.argv }),
    });
    item.status = run.status;
    item.runId = run.id;
    state.activeRunId = run.id;
  } catch (error) {
    item.status = "failed";
    item.error = error.message || String(error);
  } finally {
    state.queueBusy = false;
    saveBatchQueue();
    if (state.bottomTab === "queue") renderBottom();
  }
  await refreshRuns();
  if (!state.activeRunId) window.setTimeout(() => void processBatchQueue(), 0);
}

export function moveQueueItem(id, direction) {
  const index = state.batchQueue.findIndex((item) => item.id === id);
  const target = index + direction;
  if (index < 0 || target < 0 || target >= state.batchQueue.length) return;
  if (state.batchQueue[index].status !== "pending" || state.batchQueue[target].status !== "pending") return;
  [state.batchQueue[index], state.batchQueue[target]] = [state.batchQueue[target], state.batchQueue[index]];
  saveBatchQueue();
  renderBottom();
}

export async function cancelQueueItem(id) {
  const item = state.batchQueue.find((candidate) => candidate.id === id);
  if (!item) return;
  if (item.runId && ACTIVE_RUN_STATUSES.has(item.status)) {
    await api(`/api/runs/${encodeURIComponent(item.runId)}/cancel`, {
      method: "POST",
      headers: originHeaders(),
    });
  }
  item.status = "cancelled";
  saveBatchQueue();
  await refreshRuns();
}

export function clearFinishedQueueItems() {
  state.batchQueue = state.batchQueue.filter((item) => ["pending", "running", "starting", "queued"].includes(item.status));
  saveBatchQueue();
  renderBottom();
}

export async function refreshRuns() {
  const payload = await api("/api/runs");
  state.runs = payload.runs;
  const active = state.runs.find((run) => ACTIVE_RUN_STATUSES.has(run.status)) || null;
  state.activeRunId = active?.id || null;
  syncBatchQueue();
  connectRunEvents(active?.id || null);
  updateActionStates();
  renderTabCounts();
  if (["runs", "queue"].includes(state.bottomTab)) renderBottom();
  if (!state.activeRunId) window.setTimeout(() => void processBatchQueue(), 0);
}

export function connectRunEvents(runId) {
  if (state.eventRunId === runId) return;
  state.eventSource?.close();
  state.eventSource = null;
  state.eventRunId = runId;
  if (!runId) return;
  const source = new EventSource(`/api/runs/${encodeURIComponent(runId)}/events`);
  state.eventSource = source;
  source.addEventListener("log", (event) => {
    const payload = JSON.parse(event.data);
    state.logs.push(`[${payload.stream}] ${payload.line}`);
    if (state.logs.length > MAX_LOG_LINES) state.logs.splice(0, state.logs.length - MAX_LOG_LINES);
    if (state.bottomTab === "logs") scheduleLogRender();
  });
  source.addEventListener("artifact", () => refreshArtifacts().catch(reportError));
  source.addEventListener("stage", () => {
    if (state.selectedRunDetails?.run?.id === runId) refreshOpenRunDetails(runId).catch(reportError);
  });
  source.addEventListener("status", (event) => {
    const payload = JSON.parse(event.data);
    if (!ACTIVE_RUN_STATUSES.has(payload.status)) {
      source.close();
      state.eventSource = null;
      state.eventRunId = null;
    }
    refreshRuns().then(() => {
      if (payload.status === "succeeded") return showCompletedRunQc(runId);
      return null;
    }).catch(reportError);
  });
}

export async function cancelRun() {
  if (!state.activeRunId) throw new Error(translate("No run is active"));
  const runId = state.activeRunId;
  await api(`/api/runs/${encodeURIComponent(runId)}/cancel`, {
    method: "POST",
    headers: originHeaders(),
  });
  toast(translate("Cancelled run {id}", { id: runId }));
  await refreshRuns();
}

export async function refreshArtifacts() {
  const payload = await api("/api/artifacts");
  state.artifacts = payload.artifacts;
  renderTabCounts();
  if (state.bottomTab === "artifacts") renderBottom();
}

export function renderTabCounts() {
  const counts = {
    runs: state.runs.length,
    queue: state.batchQueue.filter((item) => item.status === "pending").length,
    artifacts: state.artifacts.length,
  };
  Object.entries(counts).forEach(([tab, count]) => {
    const badge = document.querySelector(`[data-bottom-tab="${tab}"] .tab-count`);
    if (!badge) return;
    badge.textContent = String(count);
    badge.hidden = count === 0;
  });
}

// Coalesce bursts of log events into one render per frame.
export function scheduleLogRender() {
  if (state.logRenderPending) return;
  state.logRenderPending = true;
  window.requestAnimationFrame(() => {
    state.logRenderPending = false;
    if (state.bottomTab === "logs") renderLogs(byId("bottom-content"));
  });
}

export function renderLogs(content) {
  // Follow new output only while the reader is already at the end.
  const following = content.scrollHeight - content.scrollTop - content.clientHeight < 24;
  content.textContent = state.logs.length ? state.logs.join("\n") : translate("No log events yet.");
  if (following) content.scrollTop = content.scrollHeight;
}

export function renderBottom() {
  const content = byId("bottom-content");
  content.removeAttribute("data-i18n");
  content.dataset.tab = state.bottomTab;
  content.classList.remove("qc-panel");
  content.replaceChildren();
  if (state.bottomTab === "runs") {
    renderRunList(content);
    return;
  }
  if (state.bottomTab === "queue") {
    renderQueue(content);
    return;
  }
  if (state.bottomTab === "logs") {
    renderLogs(content);
    return;
  }
  if (state.bottomTab === "qc") {
    renderQcPanel(content, state.qcDetails?.qc);
    return;
  }
  if (!state.artifacts.length) {
    content.textContent = translate("No artifacts yet.");
    return;
  }
  state.artifacts.forEach((artifact) => {
    const link = document.createElement("a");
    link.className = "artifact-link";
    link.href = `/api/artifacts/${encodeURIComponent(artifact.id)}`;
    link.target = "_blank";
    link.rel = "noopener";
    link.textContent = artifact.path;
    content.append(link);
  });
}

export function renderRunList(content) {
  if (!state.runs.length) {
    content.append(translatedElement("p", "empty-state", "No runs yet."));
    return;
  }
  const toolbar = document.createElement("div");
  toolbar.className = "bottom-toolbar";
  const compare = translatedElement("button", "", "Compare selected");
  compare.type = "button";
  compare.disabled = state.compareRunIds.size !== 2;
  compare.addEventListener("click", () => openRunComparison().catch(reportError));
  toolbar.append(compare);
  const list = document.createElement("div");
  list.className = "run-list";
  state.runs.forEach((run) => list.append(runListRow(run)));
  content.append(toolbar, list);
}

export function runListRow(run) {
  const row = document.createElement("div");
  row.className = "run-row";
  const compare = document.createElement("input");
  compare.type = "checkbox";
  compare.checked = state.compareRunIds.has(run.id);
  compare.setAttribute("aria-label", translate("Select for comparison"));
  compare.addEventListener("change", () => toggleRunComparison(run.id, compare));
  const open = document.createElement("button");
  open.type = "button";
  open.className = "run-open";
  open.addEventListener("click", () => openRunDetails(run.id).catch(reportError));
  const statusBadge = translatedElement("span", `status-badge ${run.status}`, run.status);
  const command = document.createElement("strong");
  command.textContent = run.command;
  const meta = document.createElement("span");
  meta.textContent = `${formatDate(run.created_at)} · ${shortRunId(run.id)} · ${formatRunDuration(run)}`;
  open.append(statusBadge, command, meta);
  row.append(compare, open, miniRunProgress(run));
  return row;
}

export function miniRunProgress(run) {
  const bar = document.createElement("span");
  bar.className = `mini-run-progress ${run.status}`;
  for (let index = 0; index < 4; index += 1) bar.append(document.createElement("i"));
  return bar;
}

export function toggleRunComparison(runId, checkbox) {
  if (checkbox.checked && state.compareRunIds.size >= 2) {
    checkbox.checked = false;
    toast(translate("Select exactly two runs to compare"), true);
    return;
  }
  if (checkbox.checked) state.compareRunIds.add(runId);
  else state.compareRunIds.delete(runId);
  renderBottom();
}

export function renderQueue(content) {
  if (!state.batchQueue.length) {
    content.append(translatedElement("p", "empty-state", "No queued analyses."));
    return;
  }
  const complete = state.batchQueue.filter((item) => !["pending", "queued", "starting", "running"].includes(item.status)).length;
  const toolbar = document.createElement("div");
  toolbar.className = "bottom-toolbar queue-summary";
  const label = translatedElement("span", "", "Queue progress");
  const progress = document.createElement("progress");
  progress.max = state.batchQueue.length;
  progress.value = complete;
  const count = document.createElement("span");
  count.textContent = `${complete}/${state.batchQueue.length}`;
  const clear = translatedElement("button", "", "Clear finished");
  clear.type = "button";
  clear.addEventListener("click", clearFinishedQueueItems);
  toolbar.append(label, progress, count, clear);
  const list = document.createElement("div");
  list.className = "queue-list";
  state.batchQueue.forEach((item, index) => list.append(queueListRow(item, index)));
  const note = translatedElement(
    "small",
    "queue-note",
    "Queue runs while Studio is open and resumes after a refresh.",
  );
  content.append(toolbar, note, list);
}

export function queueListRow(item, index) {
  const row = document.createElement("div");
  row.className = "queue-row";
  const statusBadge = translatedElement("span", `status-badge ${item.status}`, item.status);
  const text = document.createElement("span");
  const name = document.createElement("strong");
  name.textContent = item.workflow;
  const detail = document.createElement("small");
  detail.textContent = item.error || item.command;
  text.append(name, detail);
  const actions = document.createElement("span");
  actions.className = "queue-actions";
  actions.append(
    queueActionButton("↑", "Move up", () => moveQueueItem(item.id, -1), index === 0 || item.status !== "pending"),
    queueActionButton("↓", "Move down", () => moveQueueItem(item.id, 1), index === state.batchQueue.length - 1 || item.status !== "pending"),
    queueActionButton("×", "Cancel queued item", () => cancelQueueItem(item.id).catch(reportError), !["pending", "queued", "starting", "running"].includes(item.status)),
  );
  row.append(statusBadge, text, actions);
  return row;
}

export function queueActionButton(text, label, handler, disabled) {
  const button = document.createElement("button");
  button.type = "button";
  button.textContent = text;
  button.title = translate(label);
  button.setAttribute("aria-label", translate(label));
  button.disabled = disabled;
  button.addEventListener("click", handler);
  return button;
}

export function navigateBottomTabs(event) {
  const tabs = [...document.querySelectorAll("[data-bottom-tab]")];
  const index = tabs.indexOf(document.activeElement);
  if (index < 0) return;
  const targets = { ArrowRight: index + 1, ArrowLeft: index - 1, Home: 0, End: tabs.length - 1 };
  if (!(event.key in targets)) return;
  event.preventDefault();
  const next = tabs[(targets[event.key] + tabs.length) % tabs.length];
  next.focus();
  showBottomTab(next.dataset.bottomTab);
}

export function showBottomTab(tab) {
  state.bottomTab = tab;
  byId("bottom-panel").classList.remove("collapsed");
  document.querySelectorAll("[data-bottom-tab]").forEach((button) => {
    const selected = button.dataset.bottomTab === tab;
    button.classList.toggle("active", selected);
    button.setAttribute("aria-selected", String(selected));
    button.tabIndex = selected ? 0 : -1;
    if (selected) byId("bottom-content").setAttribute("aria-labelledby", button.id);
  });
  if (tab === "artifacts") refreshArtifacts().catch(reportError);
  if (tab === "qc") loadLatestQc().catch(reportError);
  renderBottom();
  resizeWorkspace();
}
