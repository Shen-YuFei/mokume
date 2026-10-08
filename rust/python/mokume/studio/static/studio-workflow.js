// Workflow catalog, command arguments, workflow templates, and Agent parameter patches.
import {
  WORKFLOW_TEMPLATE_SCHEMA_VERSION,
  api,
  byId,
  originHeaders,
  reportError,
  state,
  toast,
} from "./studio-core.js";
import { setTranslatedAttribute, setTranslatedText, translate } from "./studio-i18n.js";
import { updateActionStates } from "./studio-menus.js";
import { createFileIcon, filePresentation, refreshFiles } from "./studio-files.js";
import {
  clearArgumentControl,
  fillTemplateControl,
  highlightAgentFields,
  renderCommandForm,
} from "./studio-form.js";
import { configureConditionalParameters } from "./studio-workflow-rules.js";

export function agentWorkflowFormState() {
  if (!state.selectedCommand) return null;
  const parameters = {};
  byId("command-form").querySelectorAll(".argument-control").forEach((control) => {
    try {
      const value = argumentTemplateValue(control);
      if (value !== null) parameters[control.dataset.flag] = value;
    } catch (_) {
      // Keep an unfinished field local until the user completes it.
    }
  });
  return { workflow: [...state.selectedCommand.path], parameters };
}

export async function refreshCommands() {
  const payload = await api("/api/commands");
  state.commands = payload.commands;
  state.selectedCommand = null;
  const catalog = byId("command-catalog");
  catalog.classList.remove("compact");
  catalog.replaceChildren();
  const categories = new Map();
  state.commands.forEach((command, index) => {
    const category = command.category || "Workflow";
    let actions = categories.get(category);
    if (!actions) {
      const group = document.createElement("section");
      group.className = "command-family";
      const label = document.createElement("span");
      label.className = "command-family-label";
      setTranslatedText(label, category);
      actions = document.createElement("div");
      actions.className = "command-family-actions";
      group.append(label, actions);
      catalog.append(group);
      categories.set(category, actions);
    }
    const button = document.createElement("button");
    button.className = "command-card";
    button.type = "button";
    button.setAttribute("aria-pressed", "false");
    button.textContent = command.display_name;
    const { summary, requirement } = commandSummary(command.help);
    if (summary) {
      const description = document.createElement("span");
      description.className = "command-card-description";
      setTranslatedText(description, summary);
      button.append(description);
      setTranslatedAttribute(button, "title", summary);
    }
    if (requirement) {
      const extra = document.createElement("span");
      extra.className = "command-card-requirement";
      extra.textContent = requirement;
      button.append(extra);
    }
    button.addEventListener("click", () => selectCommand(index));
    actions.append(button);
  });
  updateActionStates();
}

// Split "Build a report [requires: mokume[analysis]]" into its summary and optional extra.
export function commandSummary(help) {
  const match = /^(.*?)\s*\[requires: (.+)\]$/.exec(help || "");
  return match ? { summary: match[1], requirement: match[2] } : { summary: help || "", requirement: null };
}

export function selectCommand(index) {
  state.selectedCommand = state.commands[index];
  // Once a workflow is chosen, the catalog shrinks to a row of names.
  byId("command-catalog").classList.add("compact");
  document.querySelectorAll(".command-card").forEach((card, cardIndex) => {
    card.classList.toggle("active", cardIndex === index);
    card.setAttribute("aria-pressed", String(cardIndex === index));
  });
  renderCommandForm(state.selectedCommand);
  updateActionStates();
}

export function commandArgv() {
  if (!state.selectedCommand) throw new Error(translate("Select a workflow first"));
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
      if (values.some((value) => !value)) {
        throw new Error(translate("{option} requires {count} values", { option, count: values.length }));
      }
      argv.push(option, ...values);
      occurrences += 1;
    });
    if (control.dataset.required === "true" && occurrences === 0) {
      throw new Error(translate("{option} is required", { option }));
    }
  });
  return argv;
}

export function workflowTemplateParameters() {
  const parameters = {};
  byId("command-form").querySelectorAll(".argument-control").forEach((control) => {
    const value = argumentTemplateValue(control);
    if (value !== null) parameters[control.dataset.flag] = value;
  });
  return parameters;
}

export function argumentTemplateValue(control) {
  if (control.dataset.boolean === "true") {
    return control.querySelector('input[type="checkbox"]')?.checked ? true : null;
  }
  const rows = [...control.querySelectorAll(".value-row")].map((row) =>
    [...row.querySelectorAll("[data-value]")].map((input) => input.value.trim()),
  ).filter((values) => values.some(Boolean));
  const incomplete = rows.find((values) => values.some((value) => !value));
  if (incomplete) {
    throw new Error(translate("{option} requires {count} values", {
      option: `--${control.dataset.flag}`,
      count: incomplete.length,
    }));
  }
  if (!rows.length) return null;
  const repeated = control.dataset.repeat === "true";
  const singleValue = Number(control.dataset.valueCount) === 1;
  if (singleValue) return repeated ? rows.map(([value]) => value) : rows[0][0];
  return repeated ? rows : rows[0];
}

export function currentWorkflowTemplate() {
  if (!state.selectedCommand) throw new Error(translate("Select a workflow first"));
  return {
    $schemaVersion: WORKFLOW_TEMPLATE_SCHEMA_VERSION,
    workflow: [...state.selectedCommand.path],
    parameters: workflowTemplateParameters(),
  };
}

export function workflowTemplatePath(directory, name) {
  return directory === "." ? name : `${directory}/${name}`;
}

export function workflowTemplateDirectoryLabel(directory) {
  if (directory === ".") return state.project.root;
  return `${state.project.root.replace(/[\\/]$/, "")}/${directory}`;
}

export function updateWorkflowTemplateConfirm() {
  const confirm = byId("workflow-template-confirm");
  confirm.disabled = state.workflowTemplateMode === "import"
    ? !state.workflowTemplateSelected
    : !byId("workflow-template-name").value.trim();
}

export function renderWorkflowTemplateEntries() {
  const list = byId("workflow-template-list");
  list.replaceChildren();
  const entries = state.workflowTemplateEntries.filter(
    (entry) => entry.kind === "directory" || entry.name.toLowerCase().endsWith(".json"),
  );
  if (!entries.length) {
    const empty = document.createElement("div");
    empty.className = "empty-state";
    setTranslatedText(empty, "No JSON templates in this folder.");
    list.append(empty);
    return;
  }
  entries.forEach((entry) => {
    const row = document.createElement("button");
    row.type = "button";
    row.className = "workflow-template-row";
    row.dataset.path = entry.path;
    const icon = createFileIcon(entry.kind === "directory"
      ? { icon: "folder", tone: "folder" }
      : filePresentation(entry.name));
    const name = document.createElement("span");
    name.textContent = entry.name;
    row.append(icon, name);
    if (entry.kind === "directory") {
      row.classList.add("directory");
      row.addEventListener("click", () => loadWorkflowTemplateDirectory(entry.path).catch(reportError));
    } else {
      row.setAttribute("aria-pressed", String(entry.path === state.workflowTemplateSelected));
      row.addEventListener("click", () => selectWorkflowTemplateFile(entry));
    }
    list.append(row);
  });
}

export function selectWorkflowTemplateFile(entry) {
  state.workflowTemplateSelected = entry.path;
  byId("workflow-template-name").value = entry.name;
  byId("workflow-template-list").querySelectorAll(".workflow-template-row:not(.directory)")
    .forEach((row) => row.setAttribute("aria-pressed", String(row.dataset.path === entry.path)));
  updateWorkflowTemplateConfirm();
}

export async function loadWorkflowTemplateDirectory(path = ".") {
  const payload = await api(`/api/files?path=${encodeURIComponent(path)}`);
  state.workflowTemplateDirectory = payload.path;
  state.workflowTemplateParent = payload.parent;
  state.workflowTemplateEntries = payload.entries;
  state.workflowTemplateSelected = null;
  byId("workflow-template-parent").disabled = payload.parent === null;
  byId("workflow-template-current").value = workflowTemplateDirectoryLabel(payload.path);
  if (state.workflowTemplateMode === "import") byId("workflow-template-name").value = "";
  renderWorkflowTemplateEntries();
  updateWorkflowTemplateConfirm();
}

export async function openWorkflowTemplateDialog(mode) {
  if (!state.project) throw new Error(translate("Open a folder first"));
  if (mode === "export" && !state.selectedCommand) {
    throw new Error(translate("Select a workflow first"));
  }
  state.workflowTemplateMode = mode;
  const importing = mode === "import";
  const name = byId("workflow-template-name");
  name.readOnly = importing;
  name.value = importing ? "" : `mokume-${state.selectedCommand.path.join("-")}-template.json`;
  setTranslatedAttribute(name, "placeholder", importing
    ? "Select a JSON template"
    : "Enter a JSON file name");
  setTranslatedText(byId("workflow-template-title"), importing
    ? "Import workflow template"
    : "Export workflow template");
  setTranslatedText(byId("workflow-template-confirm"), importing
    ? "Import from workspace"
    : "Save to workspace");
  await loadWorkflowTemplateDirectory(".");
  byId("workflow-template-dialog").showModal();
  if (!importing) name.select();
}

export async function submitWorkflowTemplate(event) {
  event.preventDefault();
  if (state.workflowTemplateMode === "import") {
    const path = state.workflowTemplateSelected;
    if (!path) return;
    const payload = await api(`/api/workflow-template?path=${encodeURIComponent(path)}`);
    applyWorkflowTemplate(payload.template);
    byId("workflow-template-dialog").close();
    toast(translate("Workflow template imported"));
    return;
  }
  const name = byId("workflow-template-name").value.trim();
  if (!name.toLowerCase().endsWith(".json") || name.includes("/") || name.includes("\\")) {
    throw new Error(translate("Template file name must end with .json"));
  }
  const path = workflowTemplatePath(state.workflowTemplateDirectory, name);
  const exists = state.workflowTemplateEntries.some(
    (entry) => entry.kind === "file" && entry.path === path,
  );
  if (exists && !window.confirm(translate("Replace the existing workflow template?"))) return;
  await api("/api/workflow-template", {
    method: "PUT",
    headers: originHeaders(),
    body: JSON.stringify({
      path,
      overwrite: exists,
      template: currentWorkflowTemplate(),
    }),
  });
  byId("workflow-template-dialog").close();
  await refreshFiles();
  toast(translate("Workflow template exported"));
}

export function applyWorkflowTemplate(template) {
  if (!template || typeof template !== "object" || Array.isArray(template)) {
    throw new Error(translate("Invalid workflow template"));
  }
  if (template.$schemaVersion !== WORKFLOW_TEMPLATE_SCHEMA_VERSION) {
    throw new Error(translate("Unsupported workflow template version"));
  }
  if (!Array.isArray(template.workflow) || !template.workflow.length
      || !template.workflow.every((part) => typeof part === "string")) {
    throw new Error(translate("Invalid workflow template"));
  }
  if (!template.parameters || typeof template.parameters !== "object"
      || Array.isArray(template.parameters)) {
    throw new Error(translate("Invalid workflow template"));
  }
  const workflow = template.workflow.join(" ");
  const commandIndex = state.commands.findIndex((command) => command.path.join(" ") === workflow);
  if (commandIndex < 0) {
    throw new Error(translate("Template workflow is unavailable: {workflow}", { workflow }));
  }
  const command = state.commands[commandIndex];
  const flags = new Map(command.flags
    .filter((flag) => !flag.global && !flag.studio_hidden)
    .map((flag) => [flag.long || flag.id, flag]));
  const values = new Map();
  Object.entries(template.parameters).forEach(([parameter, value]) => {
    const flag = flags.get(parameter);
    if (!flag) {
      throw new Error(translate("Template parameter is unavailable: --{parameter}", { parameter }));
    }
    values.set(parameter, normalizeTemplateParameter(flag, value));
  });
  selectCommand(commandIndex);
  const form = byId("command-form");
  form.querySelectorAll(".argument-control").forEach(clearArgumentControl);
  values.forEach((value, parameter) => {
    const control = [...form.querySelectorAll(".argument-control")]
      .find((candidate) => candidate.dataset.flag === parameter);
    fillTemplateControl(control, flags.get(parameter), value);
  });
  configureConditionalParameters(command, form);
}

export function applyWorkflowParameterPatch(patch) {
  if (!patch || patch.type !== "workflow_parameter_patch"
      || !Array.isArray(patch.workflow)
      || !patch.parameters || typeof patch.parameters !== "object"
      || Array.isArray(patch.parameters)) {
    throw new Error(translate("Invalid workflow template"));
  }
  const workflow = patch.workflow.join(" ");
  if (!state.selectedCommand
      || state.selectedCommand.path.join(" ") !== workflow) {
    throw new Error(translate("The Agent response no longer matches the current workflow"));
  }
  const command = state.selectedCommand;
  const flags = new Map(command.flags
    .filter((flag) => !flag.global && !flag.studio_hidden)
    .map((flag) => [flag.long || flag.id, flag]));
  const flagsById = new Map(command.flags.map((flag) => [flag.id, flag]));
  const values = new Map();
  Object.entries(patch.parameters).forEach(([parameter, value]) => {
    const flag = flags.get(parameter);
    if (!flag) {
      throw new Error(translate("Template parameter is unavailable: --{parameter}", { parameter }));
    }
    values.set(parameter, value === null ? null : normalizeTemplateParameter(flag, value));
  });
  const form = byId("command-form");
  const controlFor = (parameter) => [...form.querySelectorAll(".argument-control")]
    .find((candidate) => candidate.dataset.flag === parameter);
  values.forEach((_value, parameter) => clearArgumentControl(controlFor(parameter)));
  values.forEach((value, parameter) => {
    if (value === null || value === false) return;
    const flag = flags.get(parameter);
    (flag.conflicts || []).forEach((conflict) => {
      const conflictFlag = flagsById.get(conflict) || flags.get(conflict);
      if (!conflictFlag) return;
      const conflictName = conflictFlag.long || conflictFlag.id;
      if (!values.has(conflictName)) clearArgumentControl(controlFor(conflictName));
    });
  });
  values.forEach((value, parameter) => {
    if (value !== null) fillTemplateControl(controlFor(parameter), flags.get(parameter), value);
  });
  configureConditionalParameters(command, form);
  const changedFields = [...values.keys()].map((parameter) =>
    [...form.querySelectorAll(".form-field")]
      .find((field) => field.dataset.flag === parameter),
  ).filter(Boolean);
  highlightAgentFields(changedFields);
  toast(translate("Agent updated {count} parameters", { count: changedFields.length }));
}

export function normalizeTemplateParameter(flag, value) {
  const parameter = flag.long || flag.id;
  const count = Number(flag.value_arity?.max || flag.value_arity?.min || 0);
  if (count === 0) {
    if (typeof value !== "boolean") throw invalidTemplateValue(parameter);
    return value;
  }
  let rows;
  if (count === 1) {
    const values = flag.repeat && Array.isArray(value) ? value : [value];
    rows = values.map((item) => [templateScalar(item, parameter)]);
  } else {
    const values = flag.repeat ? value : [value];
    if (!Array.isArray(values)) throw invalidTemplateValue(parameter);
    rows = values.map((row) => {
      if (!Array.isArray(row) || row.length !== count) throw invalidTemplateValue(parameter);
      return row.map((item) => templateScalar(item, parameter));
    });
  }
  const choices = flag.possible_values || [];
  if (choices.length && rows.some(([item]) => !choices.includes(item))) {
    throw invalidTemplateValue(parameter);
  }
  return rows;
}

export function templateScalar(value, parameter) {
  if (!(["string", "number"].includes(typeof value))) throw invalidTemplateValue(parameter);
  const text = String(value).trim();
  if (!text) throw invalidTemplateValue(parameter);
  return text;
}

export function invalidTemplateValue(parameter) {
  return new Error(translate("Invalid template value for --{parameter}", { parameter }));
}
