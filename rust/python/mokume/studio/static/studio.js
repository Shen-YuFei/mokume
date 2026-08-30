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
};

const byId = (id) => document.getElementById(id);
const originHeaders = () => ({
  "Content-Type": "application/json",
  "X-CSRF-Token": state.csrf,
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
  document.querySelectorAll(".menu-popover.open").forEach((menu) => menu.classList.remove("open"));
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
  if (focusFirst) menuItems(menu)[0]?.focus();
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

function updateActionStates() {
  const canConfigure = Boolean(state.project && state.selectedCommand);
  setActionDisabled("close-project", !state.project);
  setActionDisabled("validate-command", !canConfigure);
  setActionDisabled("run-command", !canConfigure || Boolean(state.activeRunId));
  setActionDisabled("cancel-run", !state.activeRunId);
}

async function refreshProject() {
  state.project = await api("/api/project");
  const hasProject = Boolean(state.project);
  byId("project-chip").textContent = hasProject ? state.project.root : "No project";
  byId("welcome").classList.toggle("hidden", hasProject);
  byId("workflow").classList.toggle("hidden", !hasProject);
  if (hasProject) {
    await Promise.all([refreshFiles(), refreshCommands()]);
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

function resizeWorkspace() {
  const collapsed = byId("bottom-panel").classList.contains("collapsed");
  document.documentElement.style.setProperty("--bottom-height", collapsed ? "39px" : "188px");
  document.querySelector(".workspace").style.height = `calc(100vh - var(--header-height) - ${collapsed ? 39 : 188}px)`;
}

async function performAction(action) {
  closeMenus();
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
    "toggle-sidebar": async () => document.querySelector(".workspace").classList.toggle("sidebar-hidden"),
    "toggle-bottom": async () => toggleBottom(),
    "show-artifacts": async () => showBottomTab("artifacts"),
    fullscreen: async () => document.fullscreenElement ? document.exitFullscreen() : document.documentElement.requestFullscreen(),
    shortcuts: async () => showInfo("Keyboard Shortcuts", "Ctrl+O  Open Folder\nAlt+1  Toggle Sidebar\nAlt+2  Toggle Bottom Panel\nAlt+3  Show Artifacts\nCtrl+Shift+Enter  Validate\nCtrl+Enter  Run\nCtrl+.  Cancel Run\nCtrl+/  Keyboard Shortcuts\nF11  Full Screen"),
    "system-status": async () => showInfo("System Status", JSON.stringify(await api("/api/system"), null, 2)),
    about: async () => showInfo("About Mokume", `Mokume Studio ${state.version}\nA local, Rust-backed proteomics analysis workspace.`),
    "assistant-settings": async () => showInfo("Assistant", "Provider configuration will be available in the AI phase. Mokume compute remains fully usable without AI."),
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
      openMenu(trigger);
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
    menu.addEventListener("keydown", (event) => navigateMenu(event, menu, triggers));
  });
  document.addEventListener("click", () => closeMenus());
  window.addEventListener("resize", () => closeMenus());
}

function navigateMenu(event, menu, triggers) {
  const items = menuItems(menu);
  const current = items.indexOf(document.activeElement);
  if (["ArrowDown", "ArrowUp", "Home", "End"].includes(event.key)) {
    event.preventDefault();
    let next = event.key === "Home" ? 0 : items.length - 1;
    if (event.key === "ArrowDown") next = (current + 1) % items.length;
    if (event.key === "ArrowUp") next = (current - 1 + items.length) % items.length;
    items[next]?.focus();
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
  document.addEventListener("keydown", keyboardShortcut);
}

function keyboardShortcut(event) {
  if (event.defaultPrevented) return;
  const control = event.ctrlKey || event.metaKey;
  let action = null;
  if (control && event.key.toLowerCase() === "o") action = "open-folder";
  if (control && event.key === "Enter") action = event.shiftKey ? "validate-command" : "run-command";
  if (control && event.key === ".") action = "cancel-run";
  if (control && event.key === "/") action = "shortcuts";
  if (event.altKey && event.key === "1") action = "toggle-sidebar";
  if (event.altKey && event.key === "2") action = "toggle-bottom";
  if (event.altKey && event.key === "3") action = "show-artifacts";
  if (event.key === "F11") action = "fullscreen";
  if (!action) return;
  event.preventDefault();
  performAction(action).catch(reportError);
}

async function boot() {
  const session = await api("/api/session");
  state.csrf = session.csrf_token;
  state.version = session.version;
  bindEvents();
  await Promise.all([refreshProject(), refreshRuns(), refreshArtifacts()]);
  renderBottom();
}

boot().catch(reportError);
