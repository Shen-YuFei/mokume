const ACTIVE_RUN_STATUSES = new Set(["queued", "starting", "running", "cancelling"]);
const state = {
  csrf: null,
  version: null,
  project: null,
  folderPath: null,
  folderParent: null,
  commands: [],
  selectedCommand: null,
  runs: [],
  activeRunId: null,
  eventSource: null,
  eventRunId: null,
  logs: [],
  artifacts: [],
  bottomTab: "runs",
  provider: null,
  dataset: null,
  projectId: null,
  agentBusy: false,
  agentAbort: null,
  pendingApproval: null,
  threads: { ask: crypto.randomUUID(), plan: crypto.randomUUID() },
};

const byId = (id) => document.getElementById(id);
const originHeaders = () => ({
  "Content-Type": "application/json",
  "X-CSRF-Token": state.csrf,
});

const agentHeaders = () => ({
  ...originHeaders(),
  Accept: "text/event-stream",
});

async function api(path, options = {}) {
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

function toast(message, error = false) {
  const element = byId("toast");
  element.textContent = message;
  element.classList.toggle("error", error);
  element.classList.add("visible");
  window.clearTimeout(toast.timer);
  toast.timer = window.setTimeout(() => element.classList.remove("visible"), 3200);
}

function menuItems(menu) {
  return [...menu.querySelectorAll('[role="menuitem"]:not([disabled])')];
}

function closeMenus({ restoreFocus = false } = {}) {
  const activeTrigger = document.querySelector(".menu-trigger.active");
  document.querySelectorAll(".menu-popover.open").forEach((menu) => {
    menu.classList.remove("open");
    menuItems(menu).forEach((item) => { item.tabIndex = -1; });
  });
  document.querySelectorAll(".menu-trigger.active").forEach((trigger) => {
    trigger.classList.remove("active");
    trigger.setAttribute("aria-expanded", "false");
  });
  if (restoreFocus) activeTrigger?.focus();
}

function openMenu(trigger, { focusFirst = false } = {}) {
  const menu = byId(trigger.dataset.menu);
  const alreadyOpen = menu.classList.contains("open");
  closeMenus();
  if (alreadyOpen) return;
  menu.classList.add("open");
  const bounds = trigger.getBoundingClientRect();
  const maximum = Math.max(6, window.innerWidth - menu.offsetWidth - 6);
  menu.style.left = `${Math.min(Math.max(6, bounds.left), maximum)}px`;
  trigger.classList.add("active");
  trigger.setAttribute("aria-expanded", "true");
  if (focusFirst) focusMenuItem(menu, 0);
}

function focusMenuItem(menu, index) {
  const items = menuItems(menu);
  items.forEach((item, itemIndex) => { item.tabIndex = itemIndex === index ? 0 : -1; });
  items[index]?.focus();
}

function showInfo(title, content) {
  byId("info-title").textContent = title;
  byId("info-content").textContent = content;
  byId("info-dialog").showModal();
}

function setActionDisabled(action, disabled) {
  document.querySelectorAll(`[data-action="${action}"]`).forEach((element) => {
    if (element instanceof HTMLButtonElement) element.disabled = disabled;
  });
}

function actionEnabled(action) {
  const controls = [...document.querySelectorAll(`[data-action="${action}"]`)]
    .filter((element) => element instanceof HTMLButtonElement);
  return !controls.length || controls.some((control) => !control.disabled);
}

function updateActionStates() {
  const canConfigure = Boolean(state.project && state.selectedCommand);
  setActionDisabled("open-folder", Boolean(state.activeRunId));
  setActionDisabled("close-project", !state.project || Boolean(state.activeRunId));
  setActionDisabled("inspect-dataset", !state.project || Boolean(state.activeRunId));
  setActionDisabled("validate-command", !canConfigure);
  setActionDisabled("run-command", !canConfigure || Boolean(state.activeRunId));
  setActionDisabled("cancel-run", !state.activeRunId);
  setActionDisabled("exit-studio", Boolean(state.activeRunId));
  byId("assistant-send").disabled = !state.project || !state.provider || state.agentBusy;
}

function threadsForProject(projectId) {
  if (!projectId) return { ask: crypto.randomUUID(), plan: crypto.randomUUID() };
  return Object.fromEntries(["ask", "plan"].map((mode) => {
    const key = `mokume:thread:${projectId}:${mode}`;
    let value = sessionStorage.getItem(key);
    if (!value) {
      value = crypto.randomUUID();
      sessionStorage.setItem(key, value);
    }
    return [mode, value];
  }));
}

function resetAssistantConversation() {
  const messages = [...byId("assistant-messages").querySelectorAll(".assistant-message")];
  messages.slice(1).forEach((message) => message.remove());
}

function appendAssistantMessage(role, text) {
  const messages = byId("assistant-messages");
  const article = document.createElement("article");
  article.className = `assistant-message ${role}`;
  if (role === "assistant") {
    const orb = document.createElement("div");
    orb.className = "assistant-orb";
    orb.setAttribute("aria-hidden", "true");
    article.append(orb);
  }
  const body = document.createElement("div");
  body.className = "message-body";
  const heading = document.createElement("strong");
  heading.textContent = role === "assistant" ? "Mokume Assistant" : "You";
  const paragraph = document.createElement("p");
  paragraph.textContent = text;
  body.append(heading, paragraph);
  article.append(body);
  messages.insertBefore(article, messages.querySelector(".assistant-notice"));
  messages.scrollTop = messages.scrollHeight;
  return paragraph;
}

function renderProviderState() {
  const notice = document.querySelector(".assistant-notice");
  if (state.provider) {
    notice.classList.add("configured");
    notice.querySelector("strong").textContent = `${state.provider.provider} · ${state.provider.model}`;
    notice.querySelector("p").textContent = "Provider credentials are held in memory for this Studio session only.";
  } else {
    notice.classList.remove("configured");
    notice.querySelector("strong").textContent = "AI is not configured";
    notice.querySelector("p").textContent = "Native workflows remain available. Add a provider for Ask and Plan; credentials stay in this server session.";
  }
  updateActionStates();
}

function toggleProviderBaseUrl() {
  const compatible = byId("provider-kind").value === "openai-compatible";
  const input = byId("provider-base-url");
  input.closest(".dialog-field").classList.toggle("hidden", !compatible);
  input.required = compatible;
  if (!compatible) input.value = "";
}

async function openProviderDialog() {
  closeMenus();
  const summary = await api("/api/ai/config");
  state.provider = summary;
  byId("provider-kind").value = summary?.provider || "openai";
  byId("provider-model").value = summary?.model || "";
  byId("provider-base-url").value = summary?.base_url || "";
  byId("provider-api-key").value = "";
  toggleProviderBaseUrl();
  byId("provider-dialog").showModal();
}

async function saveProvider(event) {
  event.preventDefault();
  const provider = byId("provider-kind").value;
  const apiKey = byId("provider-api-key").value.trim();
  const baseUrl = byId("provider-base-url").value.trim();
  state.provider = await api("/api/ai/config", {
    method: "POST",
    headers: originHeaders(),
    body: JSON.stringify({
      provider,
      model: byId("provider-model").value.trim(),
      api_key: apiKey || null,
      base_url: provider === "openai-compatible" ? baseUrl : null,
    }),
  });
  byId("provider-dialog").close();
  renderProviderState();
  toast(`Configured ${state.provider.provider} for this session`);
}

async function clearProvider() {
  await api("/api/ai/config", { method: "DELETE", headers: originHeaders() });
  state.provider = null;
  byId("provider-api-key").value = "";
  byId("provider-dialog").close();
  renderProviderState();
  toast("Provider credentials cleared");
}

async function refreshDataset() {
  state.dataset = await api("/api/datasets/latest");
  return state.dataset;
}

function openDatasetDialog() {
  closeMenus();
  if (!state.project) throw new Error("Open a folder before inspecting a dataset");
  byId("dataset-dialog").showModal();
}

function datasetRequest() {
  return {
    protein_matrix: byId("dataset-matrix").value.trim(),
    sdrf: byId("dataset-sdrf").value.trim(),
    contrast: [
      byId("dataset-condition-a").value.trim(),
      byId("dataset-condition-b").value.trim(),
    ],
    input_scale: byId("dataset-input-scale").value,
    peptide_counts: byId("dataset-peptide-counts").value.trim() || null,
    data_type: byId("dataset-data-type").value || null,
    quantification: byId("dataset-quantification").value.trim() || null,
    upstream_engine: byId("dataset-engine").value.trim() || null,
    factor_column: byId("dataset-factor-column").value.trim() || null,
  };
}

async function inspectDataset(event) {
  event.preventDefault();
  const payload = await api("/api/datasets/inspect", {
    method: "POST",
    headers: originHeaders(),
    body: JSON.stringify(datasetRequest()),
  });
  state.dataset = payload.dataset;
  byId("dataset-dialog").close();
  appendAssistantMessage("assistant", `Dataset inspection queued as run ${payload.run.id}.`);
  await refreshRuns();
  showBottomTab("runs");
  pollDataset(payload.dataset.id, payload.run.id).catch(reportError);
}

async function pollDataset(datasetId, runId) {
  while (state.project && state.dataset?.id === datasetId) {
    const record = await api(`/api/datasets/${encodeURIComponent(datasetId)}`);
    state.dataset = record;
    if (record.status === "ready") {
      const configs = record.result?.policy_recommendation?.configs || [];
      const names = configs.map((config) => config.name).join(", ");
      appendAssistantMessage("assistant", names
        ? `Inspection complete. Policy candidates: ${names}.`
        : "Inspection complete. The dataset profile is ready for Ask or Plan.");
      toast("Dataset inspection complete");
      return;
    }
    if (record.status === "failed") throw new Error(record.error || "Dataset inspection failed");
    const run = await api(`/api/runs/${encodeURIComponent(runId)}`);
    if (!ACTIVE_RUN_STATUSES.has(run.status)) {
      throw new Error(`Dataset inspection ended with status ${run.status}`);
    }
    await new Promise((resolve) => window.setTimeout(resolve, 750));
  }
}

function setAgentBusy(busy) {
  state.agentBusy = busy;
  byId("assistant-input").disabled = busy;
  byId("assistant-mode").disabled = busy;
  byId("assistant-send").textContent = busy ? "Working…" : "Send";
  updateActionStates();
}

function newAgentBody(
  mode,
  { message = null, resume = null, datasetId = state.dataset?.id || null } = {},
) {
  return {
    threadId: state.threads[mode],
    runId: crypto.randomUUID(),
    state: null,
    messages: message ? [{ id: crypto.randomUUID(), role: "user", content: message }] : [],
    tools: [],
    context: [],
    forwardedProps: { mode, datasetId },
    resume,
  };
}

function streamMessage(stream, event) {
  const id = event.messageId;
  if (!stream.messages.has(id)) stream.messages.set(id, appendAssistantMessage("assistant", ""));
  return stream.messages.get(id);
}

function parseToolArguments(tool) {
  try {
    return JSON.parse(tool.args || "{}");
  } catch (_) {
    return {};
  }
}

async function presentApproval(interrupt, stream, mode) {
  const toolId = interrupt.toolCallId;
  const tool = stream.tools.get(toolId);
  const args = parseToolArguments(tool || {});
  if (tool?.name !== "run_approved_evaluation" || !args.approval_id || !args.payload_hash) {
    throw new Error("The analysis approval did not contain a valid server plan");
  }
  const record = await api(`/api/approvals/${encodeURIComponent(args.approval_id)}`);
  if (record.payload_hash !== args.payload_hash) throw new Error("Approval hash mismatch");
  state.pendingApproval = {
    record,
    interruptId: interrupt.id,
    mode,
    datasetId: state.dataset?.id || null,
    decision: null,
  };
  persistPendingApproval();
  renderApproval(record);
  byId("approval-dialog").showModal();
}

function pendingApprovalKey() {
  return state.projectId ? `mokume:approval:${state.projectId}` : null;
}

function persistPendingApproval() {
  const key = pendingApprovalKey();
  const pending = state.pendingApproval;
  if (!key || !pending) return;
  sessionStorage.setItem(key, JSON.stringify({
    approvalId: pending.record.id,
    interruptId: pending.interruptId,
    mode: pending.mode,
    datasetId: pending.datasetId,
    decision: pending.decision,
  }));
}

function clearPendingApproval() {
  const key = pendingApprovalKey();
  if (key) sessionStorage.removeItem(key);
  state.pendingApproval = null;
}

async function restorePendingApproval() {
  const key = pendingApprovalKey();
  if (!key) return;
  const saved = JSON.parse(sessionStorage.getItem(key) || "null");
  if (!saved?.approvalId || !saved?.interruptId || !saved?.mode) return;
  const record = await api(`/api/approvals/${encodeURIComponent(saved.approvalId)}`);
  if (["consumed", "expired"].includes(record.status)) {
    sessionStorage.removeItem(key);
    return;
  }
  const decision = record.status === "approved" ? true
    : record.status === "rejected" ? false : null;
  state.pendingApproval = {
    record,
    interruptId: saved.interruptId,
    mode: saved.mode,
    datasetId: saved.datasetId || null,
    decision,
  };
  renderApproval(record);
  byId("approval-dialog").showModal();
}

function renderApproval(record) {
  const content = byId("approval-content");
  content.replaceChildren();
  const section = document.createElement("section");
  section.className = "approval-card";
  const title = document.createElement("div");
  title.className = "approval-card-title";
  const label = document.createElement("span");
  label.textContent = "Canonical parameters";
  const status = document.createElement("span");
  status.className = "approval-status";
  status.textContent = record.status;
  title.append(label, status);
  const list = document.createElement("dl");
  list.className = "approval-parameters";
  Object.entries(record.payload.card).forEach(([name, value]) => {
    const term = document.createElement("dt");
    term.textContent = name.replaceAll("_", " ");
    const detail = document.createElement("dd");
    detail.textContent = typeof value === "string" ? value : JSON.stringify(value);
    list.append(term, detail);
  });
  section.append(title, list);
  content.append(section);
  const decided = state.pendingApproval?.decision;
  byId("approval-approve").disabled = decided === false;
  byId("approval-reject").disabled = decided === true;
  byId("approval-approve").textContent = decided === true ? "Resume Approved Run" : "Approve and Run";
  byId("approval-reject").textContent = decided === false ? "Resume Rejection" : "Reject";
}

async function decideApproval(approved) {
  const pending = state.pendingApproval;
  if (!pending) throw new Error("No analysis plan is awaiting approval");
  if (pending.decision !== null && pending.decision !== approved) {
    throw new Error("This analysis plan already has the opposite decision");
  }
  if (pending.decision === null) {
    pending.record = await api(`/api/approvals/${encodeURIComponent(pending.record.id)}`, {
      method: "POST",
      headers: originHeaders(),
      body: JSON.stringify({ approved, payload_hash: pending.record.payload_hash }),
    });
    pending.decision = approved;
    persistPendingApproval();
  }
  byId("approval-dialog").close();
  const body = newAgentBody(pending.mode, {
    datasetId: pending.datasetId,
    resume: [{
      interruptId: pending.interruptId,
      status: "resolved",
      payload: { approved, reason: approved ? null : "Rejected in Mokume Studio" },
    }],
  });
  try {
    await runAgentRequest(body, pending.mode);
    clearPendingApproval();
    if (approved) {
      await refreshRuns();
      showBottomTab("runs");
    }
  } catch (error) {
    state.pendingApproval = pending;
    renderApproval(pending.record);
    byId("approval-dialog").showModal();
    throw error;
  }
}

async function handleAgentEvent(event, stream, mode) {
  if (event.type === "TEXT_MESSAGE_START") streamMessage(stream, event);
  if (event.type === "TEXT_MESSAGE_CONTENT") {
    const paragraph = streamMessage(stream, event);
    paragraph.textContent += event.delta;
    paragraph.closest(".assistant-messages")?.scrollTo(0, paragraph.scrollHeight);
  }
  if (event.type === "TOOL_CALL_START") {
    stream.tools.set(event.toolCallId, { name: event.toolCallName, args: "" });
  }
  if (event.type === "TOOL_CALL_ARGS") {
    const tool = stream.tools.get(event.toolCallId);
    if (tool) tool.args += event.delta;
  }
  if (event.type === "RUN_ERROR") throw new Error(event.message || "Assistant run failed");
  if (event.type === "RUN_FINISHED" && event.outcome?.type === "interrupt") {
    const [interrupt] = event.outcome.interrupts || [];
    if (!interrupt) throw new Error("Assistant paused without an approval request");
    await presentApproval(interrupt, stream, mode);
  }
}

async function consumeAgentStream(response, mode) {
  const stream = { buffer: "", messages: new Map(), tools: new Map() };
  const reader = response.body.getReader();
  const decoder = new TextDecoder();
  while (true) {
    const { value, done } = await reader.read();
    stream.buffer += decoder.decode(value || new Uint8Array(), { stream: !done });
    const blocks = stream.buffer.split("\n\n");
    stream.buffer = blocks.pop() || "";
    for (const block of blocks) {
      const data = block.split("\n").filter((line) => line.startsWith("data:"))
        .map((line) => line.slice(5).trimStart()).join("\n");
      if (data) await handleAgentEvent(JSON.parse(data), stream, mode);
    }
    if (done) break;
  }
}

async function runAgentRequest(body, mode) {
  if (!state.provider) throw new Error("Configure an AI provider first");
  if (!state.project) throw new Error("Open a folder first");
  state.agentAbort?.abort();
  state.agentAbort = new AbortController();
  setAgentBusy(true);
  try {
    const response = await fetch("/api/agent/run", {
      method: "POST",
      headers: agentHeaders(),
      body: JSON.stringify(body),
      signal: state.agentAbort.signal,
    });
    if (!response.ok) {
      let detail = `${response.status} ${response.statusText}`;
      try { detail = (await response.json()).detail || detail; } catch (_) { /* status is enough */ }
      throw new Error(detail);
    }
    await consumeAgentStream(response, mode);
  } finally {
    state.agentAbort = null;
    setAgentBusy(false);
  }
}

async function sendAssistantMessage() {
  const input = byId("assistant-input");
  const message = input.value.trim();
  if (!message) return;
  const mode = byId("assistant-mode").value;
  if (mode === "plan" && state.dataset?.status !== "ready") {
    throw new Error("Inspect a dataset before requesting an analysis plan");
  }
  appendAssistantMessage("user", message);
  input.value = "";
  await runAgentRequest(newAgentBody(mode, { message }), mode);
}

async function refreshProject() {
  state.project = await api("/api/project");
  const hasProject = Boolean(state.project);
  const projectId = state.project?.id || null;
  if (projectId !== state.projectId) {
    if (byId("approval-dialog").open) byId("approval-dialog").close();
    state.projectId = projectId;
    state.dataset = null;
    state.pendingApproval = null;
    state.threads = threadsForProject(state.projectId);
    resetAssistantConversation();
  }
  byId("project-chip").textContent = hasProject ? state.project.root : "No project";
  byId("welcome").classList.toggle("hidden", hasProject);
  byId("workflow").classList.toggle("hidden", !hasProject);
  if (hasProject) {
    await Promise.all([refreshFiles(), refreshCommands(), refreshDataset()]);
  } else {
    state.selectedCommand = null;
    state.commands = [];
    byId("project-files").textContent = "Open a folder to browse input files.";
    byId("command-catalog").replaceChildren();
    byId("command-form").textContent = "Select a workflow to configure its parameters.";
  }
  updateActionStates();
}

async function loadFolders(path = null) {
  const query = path ? `?path=${encodeURIComponent(path)}` : "";
  const payload = await api(`/api/folders${query}`);
  state.folderPath = payload.path;
  state.folderParent = payload.parent;
  byId("folder-current").value = payload.path;
  byId("folder-parent").disabled = !payload.parent;
  const list = byId("folder-list");
  list.replaceChildren();
  if (!payload.directories.length) {
    const empty = document.createElement("div");
    empty.className = "empty-state";
    empty.textContent = "No readable subfolders.";
    list.append(empty);
  }
  payload.directories.forEach((directory) => {
    const row = document.createElement("button");
    row.className = "folder-row";
    row.textContent = directory.name;
    row.dataset.path = directory.path;
    row.addEventListener("dblclick", () => loadFolders(directory.path).catch(reportError));
    row.addEventListener("click", () => {
      list.querySelectorAll(".selected").forEach((item) => item.classList.remove("selected"));
      row.classList.add("selected");
    });
    list.append(row);
  });
}

async function openFolderDialog() {
  closeMenus();
  await loadFolders(state.project?.root || null);
  byId("folder-dialog").showModal();
}

async function selectFolder() {
  const selected = byId("folder-list").querySelector(".folder-row.selected");
  const path = selected?.dataset.path || state.folderPath;
  await api("/api/projects/open", {
    method: "POST",
    headers: originHeaders(),
    body: JSON.stringify({ path }),
  });
  byId("folder-dialog").close();
  await refreshProject();
  toast(`Opened ${path}`);
}

async function refreshFiles() {
  const payload = await api("/api/files");
  const tree = byId("project-files");
  tree.replaceChildren();
  payload.entries.forEach((entry) => {
    const row = document.createElement("div");
    row.className = "file-entry";
    const icon = document.createElement("span");
    icon.className = "kind";
    icon.textContent = entry.kind === "directory" ? "▸" : "·";
    const name = document.createElement("span");
    name.textContent = entry.name;
    row.append(icon, name);
    tree.append(row);
  });
  if (!payload.entries.length) tree.textContent = "This folder is empty.";
}

async function refreshCommands() {
  const payload = await api("/api/commands");
  state.commands = payload.commands;
  state.selectedCommand = null;
  const catalog = byId("command-catalog");
  catalog.replaceChildren();
  state.commands.forEach((command, index) => {
    const button = document.createElement("button");
    button.className = "command-card";
    button.textContent = command.path.join(" ");
    button.title = command.help || "";
    button.addEventListener("click", () => selectCommand(index));
    catalog.append(button);
  });
  updateActionStates();
}

function selectCommand(index) {
  state.selectedCommand = state.commands[index];
  document.querySelectorAll(".command-card").forEach((card, cardIndex) => {
    card.classList.toggle("active", cardIndex === index);
  });
  renderCommandForm(state.selectedCommand);
  updateActionStates();
}

function renderCommandForm(command) {
  const form = byId("command-form");
  form.classList.remove("empty-state");
  form.replaceChildren();
  command.flags.filter((flag) => !flag.global).forEach((flag) => {
    const field = document.createElement("div");
    field.className = "form-field";
    const heading = document.createElement("span");
    heading.className = "field-heading";
    heading.textContent = `--${flag.long || flag.id}${flag.required ? " *" : ""}`;
    const control = buildArgumentControl(flag);
    const help = document.createElement("small");
    help.textContent = flag.help || valueHint(flag);
    field.append(heading, control, help);
    form.append(field);
  });
}

function buildArgumentControl(flag) {
  const arity = flag.value_arity || { min: 0, max: 0 };
  const control = document.createElement("div");
  control.className = "argument-control";
  control.dataset.flag = flag.long || flag.id;
  control.dataset.boolean = String(arity.max === 0);
  control.dataset.required = String(Boolean(flag.required));
  control.dataset.repeat = String(Boolean(flag.repeat));
  control.dataset.valueCount = String(arity.max || arity.min || 0);
  if (arity.max === 0) {
    const checkbox = document.createElement("input");
    checkbox.type = "checkbox";
    checkbox.setAttribute("aria-label", `--${control.dataset.flag}`);
    control.append(checkbox);
    return control;
  }
  addValueRow(control, flag);
  if (flag.repeat) {
    const add = document.createElement("button");
    add.type = "button";
    add.className = "add-value";
    add.textContent = "Add another";
    add.addEventListener("click", () => addValueRow(control, flag, true));
    control.append(add);
  }
  return control;
}

function addValueRow(control, flag, removable = false) {
  const count = Number(control.dataset.valueCount);
  const row = document.createElement("div");
  row.className = "value-row";
  row.style.setProperty("--value-count", String(count));
  for (let index = 0; index < count; index += 1) row.append(buildValueInput(flag, index));
  if (removable) {
    const remove = document.createElement("button");
    remove.type = "button";
    remove.className = "remove-value";
    remove.textContent = "Remove";
    remove.addEventListener("click", () => row.remove());
    row.append(remove);
  }
  control.insertBefore(row, control.querySelector(".add-value"));
}

function buildValueInput(flag, index) {
  const name = flag.value_names?.[index] || flag.value_names?.at(-1) || "VALUE";
  let input;
  if (flag.possible_values?.length && index === 0) {
    input = document.createElement("select");
    const blank = document.createElement("option");
    blank.value = "";
    blank.textContent = flag.default?.length ? `Default (${flag.default.join(", ")})` : "Select…";
    input.append(blank);
    flag.possible_values.forEach((value) => {
      const option = document.createElement("option");
      option.value = value;
      option.textContent = value;
      input.append(option);
    });
  } else {
    input = document.createElement("input");
    input.type = numericValue(flag) ? "number" : "text";
    input.placeholder = flag.default?.[index] ? `Default: ${flag.default[index]}` : name;
  }
  input.dataset.value = "";
  input.setAttribute("aria-label", `--${flag.long || flag.id} ${name}`);
  return input;
}

function valueHint(flag) {
  const names = flag.value_names || [];
  const hint = names.length ? names.map((name) => `<${name}>`).join(" ") : "Optional switch";
  return flag.repeat ? `${hint}; may be repeated` : hint;
}

function numericValue(flag) {
  return (flag.value_names || []).some((name) => ["N", "VALUE", "FRACTION", "CORRELATION"].includes(name));
}

function commandArgv() {
  if (!state.selectedCommand) throw new Error("Select a workflow first");
  const argv = [...state.selectedCommand.path];
  byId("command-form").querySelectorAll(".argument-control").forEach((control) => {
    const option = `--${control.dataset.flag}`;
    if (control.dataset.boolean === "true") {
      if (control.querySelector("input").checked) argv.push(option);
      return;
    }
    let occurrences = 0;
    control.querySelectorAll(".value-row").forEach((row) => {
      const values = [...row.querySelectorAll("[data-value]")].map((input) => input.value.trim());
      if (values.every((value) => !value)) return;
      if (values.some((value) => !value)) throw new Error(`${option} requires ${values.length} values`);
      argv.push(option, ...values);
      occurrences += 1;
    });
    if (control.dataset.required === "true" && occurrences === 0) throw new Error(`${option} is required`);
  });
  return argv;
}

async function validateCommand() {
  const payload = await api("/api/commands/validate", {
    method: "POST",
    headers: originHeaders(),
    body: JSON.stringify({ argv: commandArgv() }),
  });
  showBottomTab("logs");
  state.logs.push(`Validated: ${payload.argv.join(" ")}`);
  renderBottom();
  toast("Parameters are valid");
  return payload.argv;
}

async function runCommand() {
  const argv = await validateCommand();
  const payload = await api("/api/runs", {
    method: "POST",
    headers: originHeaders(),
    body: JSON.stringify({ argv }),
  });
  state.activeRunId = payload.id;
  toast(`Run ${payload.id} queued`);
  await refreshRuns();
  showBottomTab("runs");
}

function updateRunIndicator(active) {
  const dot = document.querySelector(".status-dot");
  dot.classList.remove("idle", "running", "failed");
  if (active) {
    dot.classList.add("running");
    byId("run-state").textContent = active.status;
    return;
  }
  const latest = state.runs[0];
  dot.classList.add(latest?.status === "failed" ? "failed" : "idle");
  byId("run-state").textContent = latest ? latest.status : "Idle";
}

async function refreshRuns() {
  const payload = await api("/api/runs");
  state.runs = payload.runs;
  const active = state.runs.find((run) => ACTIVE_RUN_STATUSES.has(run.status)) || null;
  state.activeRunId = active?.id || null;
  updateRunIndicator(active);
  connectRunEvents(active?.id || null);
  updateActionStates();
  if (state.bottomTab === "runs") renderBottom();
}

function connectRunEvents(runId) {
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
    state.logs = state.logs.slice(-1000);
    if (state.bottomTab === "logs") renderBottom();
  });
  source.addEventListener("artifact", () => refreshArtifacts().catch(reportError));
  source.addEventListener("status", (event) => {
    const payload = JSON.parse(event.data);
    if (!ACTIVE_RUN_STATUSES.has(payload.status)) {
      source.close();
      state.eventSource = null;
      state.eventRunId = null;
    }
    refreshRuns().catch(reportError);
  });
}

async function cancelRun() {
  if (!state.activeRunId) throw new Error("No run is active");
  const runId = state.activeRunId;
  await api(`/api/runs/${encodeURIComponent(runId)}/cancel`, {
    method: "POST",
    headers: originHeaders(),
  });
  toast(`Cancelled run ${runId}`);
  await refreshRuns();
}

async function refreshArtifacts() {
  const payload = await api("/api/artifacts");
  state.artifacts = payload.artifacts;
  if (state.bottomTab === "artifacts") renderBottom();
}

function renderBottom() {
  const content = byId("bottom-content");
  content.replaceChildren();
  if (state.bottomTab === "runs") {
    content.textContent = state.runs.length
      ? state.runs.map((run) => `${run.status.padEnd(12)} ${run.command}  ${run.id}`).join("\n")
      : "No runs yet.";
    return;
  }
  if (state.bottomTab === "logs") {
    content.textContent = state.logs.length ? state.logs.join("\n") : "No log events yet.";
    content.scrollTop = content.scrollHeight;
    return;
  }
  if (!state.artifacts.length) {
    content.textContent = "No artifacts yet.";
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

function showBottomTab(tab) {
  state.bottomTab = tab;
  byId("bottom-panel").classList.remove("collapsed");
  document.querySelectorAll("[data-bottom-tab]").forEach((button) => {
    button.classList.toggle("active", button.dataset.bottomTab === tab);
  });
  if (tab === "artifacts") refreshArtifacts().catch(reportError);
  renderBottom();
  resizeWorkspace();
}

function toggleBottom() {
  byId("bottom-panel").classList.toggle("collapsed");
  resizeWorkspace();
}

function toggleSidebar() {
  const workspace = document.querySelector(".workspace");
  if (window.matchMedia("(max-width: 680px)").matches) {
    if (!workspace.classList.contains("sidebar-mobile-visible")) {
      workspace.classList.remove("assistant-mobile-visible");
    }
    workspace.classList.toggle("sidebar-mobile-visible");
    return;
  }
  workspace.classList.toggle("sidebar-hidden");
}

function toggleAssistant() {
  const workspace = document.querySelector(".workspace");
  if (window.matchMedia("(max-width: 1050px)").matches) {
    if (!workspace.classList.contains("assistant-mobile-visible")) {
      workspace.classList.remove("sidebar-mobile-visible");
    }
    workspace.classList.toggle("assistant-mobile-visible");
    return;
  }
  workspace.classList.toggle("assistant-hidden");
}

function resizeWorkspace() {
  const collapsed = byId("bottom-panel").classList.contains("collapsed");
  document.documentElement.style.setProperty("--bottom-height", collapsed ? "39px" : "188px");
  document.querySelector(".workspace").style.height = `calc(100vh - var(--header-height) - ${collapsed ? 39 : 188}px)`;
}

async function performAction(action) {
  closeMenus({ restoreFocus: true });
  if (!actionEnabled(action)) return;
  const handlers = {
    "open-folder": openFolderDialog,
    "close-project": async () => {
      await api("/api/projects/close", { method: "POST", headers: originHeaders() });
      await refreshProject();
    },
    "exit-studio": async () => {
      await api("/api/studio/exit", { method: "POST", headers: originHeaders() });
      document.body.textContent = "Mokume Studio is stopping. You can close this tab.";
    },
    "validate-command": validateCommand,
    "run-command": runCommand,
    "cancel-run": cancelRun,
    "show-runs": async () => {
      await refreshRuns();
      showBottomTab("runs");
    },
    "toggle-sidebar": async () => toggleSidebar(),
    "toggle-assistant": async () => toggleAssistant(),
    "toggle-bottom": async () => toggleBottom(),
    "show-artifacts": async () => showBottomTab("artifacts"),
    fullscreen: async () => document.fullscreenElement ? document.exitFullscreen() : document.documentElement.requestFullscreen(),
    shortcuts: async () => showInfo("Keyboard Shortcuts", "Ctrl+O  Open Folder\nAlt+1  Toggle Sidebar\nAlt+2  Toggle Bottom Panel\nAlt+3  Show Artifacts\nAlt+4  Toggle Assistant\nCtrl+Shift+Enter  Validate\nCtrl+Enter  Run\nCtrl+.  Cancel Run\nCtrl+/  Keyboard Shortcuts\nF11  Full Screen"),
    "system-status": async () => showInfo("System Status", JSON.stringify(await api("/api/system"), null, 2)),
    about: async () => showInfo("About Mokume", `Mokume Studio ${state.version}\nA local, Rust-backed proteomics analysis workspace.`),
    "assistant-settings": openProviderDialog,
    "inspect-dataset": async () => openDatasetDialog(),
  };
  if (handlers[action]) await handlers[action]();
}

function reportError(error) {
  toast(error.message || String(error), true);
}

function bindMenuEvents() {
  const triggers = [...document.querySelectorAll(".menu-trigger")];
  triggers.forEach((trigger, index) => {
    trigger.addEventListener("click", (event) => {
      event.stopPropagation();
      openMenu(trigger, { focusFirst: true });
    });
    trigger.addEventListener("keydown", (event) => {
      if (["ArrowDown", "Enter", " "].includes(event.key)) {
        event.preventDefault();
        openMenu(trigger, { focusFirst: true });
      } else if (["ArrowLeft", "ArrowRight"].includes(event.key)) {
        event.preventDefault();
        const delta = event.key === "ArrowRight" ? 1 : -1;
        triggers[(index + delta + triggers.length) % triggers.length].focus();
      }
    });
  });
  document.querySelectorAll(".menu-popover").forEach((menu) => {
    menuItems(menu).forEach((item) => { item.tabIndex = -1; });
    menu.addEventListener("keydown", (event) => navigateMenu(event, menu, triggers));
  });
  document.querySelectorAll(".menu-popover a").forEach((link) => {
    link.addEventListener("click", () => closeMenus({ restoreFocus: true }));
  });
  document.addEventListener("click", () => closeMenus());
  document.addEventListener("focusin", (event) => {
    const activeTrigger = document.querySelector(".menu-trigger.active");
    const insideMenu = event.target.closest?.(".menu-popover");
    if (!insideMenu && event.target !== activeTrigger) closeMenus();
  });
  window.addEventListener("resize", () => {
    closeMenus();
    document.querySelector(".workspace").classList.remove("sidebar-mobile-visible", "assistant-mobile-visible");
  });
}

function navigateMenu(event, menu, triggers) {
  const items = menuItems(menu);
  const current = items.indexOf(document.activeElement);
  if (["ArrowDown", "ArrowUp", "Home", "End"].includes(event.key)) {
    event.preventDefault();
    let next = event.key === "Home" ? 0 : items.length - 1;
    if (event.key === "ArrowDown") next = (current + 1) % items.length;
    if (event.key === "ArrowUp") next = (current - 1 + items.length) % items.length;
    focusMenuItem(menu, next);
    return;
  }
  if (event.key === "Escape") {
    event.preventDefault();
    closeMenus({ restoreFocus: true });
    return;
  }
  if (["ArrowLeft", "ArrowRight"].includes(event.key)) {
    event.preventDefault();
    const trigger = document.querySelector(`[data-menu="${menu.id}"]`);
    const index = triggers.indexOf(trigger);
    const delta = event.key === "ArrowRight" ? 1 : -1;
    openMenu(triggers[(index + delta + triggers.length) % triggers.length], { focusFirst: true });
  }
}

function bindEvents() {
  bindMenuEvents();
  document.querySelectorAll("[data-action]").forEach((element) => element.addEventListener("click", (event) => {
    event.stopPropagation();
    performAction(element.dataset.action).catch(reportError);
  }));
  document.querySelectorAll("[data-close-dialog]").forEach((button) => button.addEventListener("click", () => button.closest("dialog").close()));
  document.querySelectorAll("[data-bottom-tab]").forEach((button) => button.addEventListener("click", () => showBottomTab(button.dataset.bottomTab)));
  byId("folder-parent").addEventListener("click", () => loadFolders(state.folderParent).catch(reportError));
  byId("select-folder").addEventListener("click", () => selectFolder().catch(reportError));
  byId("provider-kind").addEventListener("change", toggleProviderBaseUrl);
  byId("provider-form").addEventListener("submit", (event) => saveProvider(event).catch(reportError));
  byId("clear-provider").addEventListener("click", () => clearProvider().catch(reportError));
  byId("dataset-form").addEventListener("submit", (event) => inspectDataset(event).catch(reportError));
  byId("assistant-send").addEventListener("click", () => sendAssistantMessage().catch(reportError));
  byId("assistant-input").addEventListener("keydown", (event) => {
    if (event.key === "Enter" && !event.shiftKey) {
      event.preventDefault();
      sendAssistantMessage().catch(reportError);
    }
  });
  byId("approval-approve").addEventListener("click", () => decideApproval(true).catch(reportError));
  byId("approval-reject").addEventListener("click", () => decideApproval(false).catch(reportError));
  document.addEventListener("keydown", keyboardShortcut);
}

function keyboardShortcut(event) {
  if (event.defaultPrevented) return;
  if (event.target.closest?.("input, textarea, select, dialog")) return;
  const control = event.ctrlKey || event.metaKey;
  let action = null;
  if (control && event.key.toLowerCase() === "o") action = "open-folder";
  if (control && event.key === "Enter") action = event.shiftKey ? "validate-command" : "run-command";
  if (control && event.key === ".") action = "cancel-run";
  if (control && event.key === "/") action = "shortcuts";
  if (event.altKey && event.key === "1") action = "toggle-sidebar";
  if (event.altKey && event.key === "2") action = "toggle-bottom";
  if (event.altKey && event.key === "3") action = "show-artifacts";
  if (event.altKey && event.key === "4") action = "toggle-assistant";
  if (event.key === "F11") action = "fullscreen";
  if (!action) return;
  event.preventDefault();
  performAction(action).catch(reportError);
}

async function boot() {
  const session = await api("/api/session");
  state.csrf = session.csrf_token;
  state.version = session.version;
  state.provider = session.ai_provider;
  bindEvents();
  await Promise.all([refreshProject(), refreshRuns(), refreshArtifacts()]);
  renderProviderState();
  await restorePendingApproval();
  renderBottom();
}

boot().catch(reportError);
