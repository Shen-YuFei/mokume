// Workflow parameter form: fields, value rows, conditional visibility, and field focus.
import { byId, state } from "./studio-core.js";
import { setTranslatedAttribute, setTranslatedText } from "./studio-i18n.js";
import {
  closeParameterSelectors,
  navigateParameterSelector,
  openParameterSelector,
} from "./studio-menus.js";
import { configureConditionalParameters } from "./studio-workflow-rules.js";
import { commandSummary, invalidTemplateValue } from "./studio-workflow.js";

export function renderCommandForm(command) {
  const form = byId("command-form");
  form.classList.remove("empty-state");
  form.removeAttribute("data-i18n");
  form.replaceChildren();
  const flags = command.flags.filter((flag) => !flag.global && !flag.studio_hidden);
  const flagsById = new Map(flags.map((flag) => [flag.id, flag]));
  form.append(commandIntro(command));
  parameterGroups(command, flags).forEach((group) => {
    form.append(buildParameterGroup(group, flagsById));
  });
  configureConditionalParameters(command, form);
}

export function commandIntro(command) {
  const intro = document.createElement("header");
  intro.className = "command-form-intro";
  const title = document.createElement("h3");
  title.textContent = command.display_name;
  const path = document.createElement("code");
  path.className = "command-path";
  path.textContent = `mokume ${command.path.join(" ")}`;
  intro.append(title, path);
  const { summary, requirement } = commandSummary(command.help);
  if (requirement) {
    const extra = document.createElement("code");
    extra.className = "command-card-requirement";
    extra.textContent = requirement;
    intro.append(extra);
  }
  if (summary) {
    const description = document.createElement("p");
    setTranslatedText(description, summary);
    intro.append(description);
  }
  return intro;
}

export function parameterGroups(command, flags) {
  const groups = command.presentation?.groups;
  if (groups?.length) return groups;
  return [{
    id: "parameters",
    title: "Parameters",
    flags: flags.map((flag) => flag.id),
    common: flags.map((flag) => flag.id),
  }];
}

export function buildParameterGroup(group, flagsById) {
  const section = document.createElement("section");
  section.className = "parameter-group";
  section.dataset.group = group.id;
  const heading = document.createElement("h3");
  heading.className = "parameter-group-title";
  setTranslatedText(heading, group.title);
  section.append(heading);

  const commonIds = new Set(group.common || []);
  const flags = group.flags.map((id) => flagsById.get(id)).filter(Boolean);
  const common = flags.filter((flag) => flag.required || commonIds.has(flag.id));
  const advanced = flags.filter((flag) => !flag.required && !commonIds.has(flag.id));
  if (common.length) {
    const primary = buildParameterFields(common);
    primary.classList.add("parameter-primary");
    section.append(primary);
  }
  if (advanced.length) section.append(buildAdvancedParameters(advanced));
  return section;
}

export function buildParameterFields(flags) {
  const fields = document.createElement("div");
  fields.className = "parameter-fields";
  flags.forEach((flag) => fields.append(buildArgumentField(flag)));
  return fields;
}

export function buildAdvancedParameters(flags) {
  const details = document.createElement("details");
  details.className = "parameter-advanced";
  const summary = document.createElement("summary");
  setTranslatedText(summary, "Advanced");
  details.append(summary, buildParameterFields(flags));
  return details;
}

export function buildArgumentField(flag) {
  const field = document.createElement("div");
  field.className = "form-field";
  field.dataset.flag = flag.long || flag.id;
  const heading = document.createElement("span");
  heading.className = "field-heading";
  renderFieldHeading(heading, flag.long || flag.id, Boolean(flag.required));
  const control = buildArgumentControl(flag);
  const help = document.createElement("small");
  if (flag.help) setTranslatedText(help, flag.help);
  else setValueHint(help, flag);
  field.append(heading, control, help);
  return field;
}

export function renderFieldHeading(heading, flag, requirement) {
  const parameterSelector = heading.querySelector(".parameter-selector");
  heading.replaceChildren(parameterSelector || document.createTextNode(`--${flag}`));
  if (!requirement) return;
  const label = document.createElement("span");
  label.className = "required-label";
  setTranslatedText(label, requirement === true ? "Required" : requirement);
  heading.append(label);
}

export function buildArgumentControl(flag) {
  const arity = flag.value_arity || { min: 0, max: 0 };
  const control = document.createElement("div");
  control.className = "argument-control";
  control.dataset.flag = flag.long || flag.id;
  control.dataset.boolean = String(arity.max === 0);
  control.dataset.required = String(Boolean(flag.required));
  control.dataset.schemaRequired = String(Boolean(flag.required));
  control.dataset.repeat = String(Boolean(flag.repeat));
  control.dataset.valueCount = String(arity.max || arity.min || 0);
  control.dataset.defaultValue = flag.default?.[0] || "";
  if (arity.max === 0) {
    const checkboxControl = document.createElement("span");
    checkboxControl.className = "boolean-checkbox";
    const checkbox = document.createElement("input");
    checkbox.type = "checkbox";
    checkbox.setAttribute("aria-label", `--${control.dataset.flag}`);
    const mark = document.createElement("span");
    mark.className = "boolean-checkbox-mark";
    mark.setAttribute("aria-hidden", "true");
    checkboxControl.append(checkbox, mark);
    control.append(checkboxControl);
    return control;
  }
  addValueRow(control, flag);
  if (flag.repeat) {
    const add = document.createElement("button");
    add.type = "button";
    add.className = "add-value";
    setTranslatedText(add, "Add another");
    add.addEventListener("click", () => addValueRow(control, flag, true));
    control.append(add);
  }
  return control;
}

export function addValueRow(control, flag, removable = false) {
  const count = Number(control.dataset.valueCount);
  const row = document.createElement("div");
  row.className = "value-row";
  row.style.setProperty("--value-count", String(count));
  for (let index = 0; index < count; index += 1) row.append(buildValueInput(flag, index));
  if (removable) {
    const remove = document.createElement("button");
    remove.type = "button";
    remove.className = "remove-value";
    setTranslatedText(remove, "Remove");
    remove.addEventListener("click", () => row.remove());
    row.append(remove);
  }
  control.insertBefore(row, control.querySelector(".add-value"));
}

export function buildValueInput(flag, index) {
  const name = flag.value_names?.[index] || flag.value_names?.at(-1) || "VALUE";
  let input;
  if (flag.possible_values?.length && index === 0) {
    input = document.createElement("select");
    const blank = document.createElement("option");
    blank.value = "";
    if (flag.default?.length) {
      setTranslatedText(blank, "Default ({value})", { value: flag.default.join(", ") });
    } else {
      setTranslatedText(blank, "Select…");
    }
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
    if (flag.default?.[index]) {
      setTranslatedAttribute(input, "placeholder", "Default: {value}", { value: flag.default[index] });
    } else {
      input.placeholder = name;
    }
  }
  input.dataset.value = "";
  input.setAttribute("aria-label", `--${flag.long || flag.id} ${name}`);
  input.setAttribute("aria-required", String(Boolean(flag.required)));
  return input;
}

export function configureAlternativeParameterSelector(form, flags, stateKey, label) {
  const importedType = flags.find((flag) => argumentControlValue(form, flag));
  if (!form.querySelector(`[data-parameter-selector="${stateKey}"]`)) {
    flags.forEach((flag) => {
      const field = [...form.querySelectorAll(".form-field")]
        .find((candidate) => candidate.dataset.flag === flag);
      if (!field) return;
      const selector = document.createElement("span");
      selector.className = "parameter-selector";
      selector.dataset.parameterSelector = stateKey;
      const trigger = document.createElement("button");
      trigger.type = "button";
      trigger.className = "parameter-selector-trigger";
      trigger.setAttribute("aria-haspopup", "listbox");
      trigger.setAttribute("aria-expanded", "false");
      setTranslatedAttribute(trigger, "aria-label", label);
      const valueLabel = document.createElement("span");
      valueLabel.className = "parameter-selector-value";
      const chevron = document.createElement("span");
      chevron.className = "parameter-selector-chevron";
      chevron.setAttribute("aria-hidden", "true");
      trigger.append(valueLabel, chevron);
      const menu = document.createElement("span");
      menu.className = "parameter-selector-menu";
      menu.id = `parameter-selector-${stateKey}-${flag}`;
      menu.hidden = true;
      menu.setAttribute("role", "listbox");
      setTranslatedAttribute(menu, "aria-label", label);
      trigger.setAttribute("aria-controls", menu.id);
      flags.forEach((value) => {
        const option = document.createElement("button");
        option.type = "button";
        option.className = "parameter-selector-option";
        option.dataset.value = value;
        option.setAttribute("role", "option");
        option.setAttribute("aria-selected", "false");
        option.tabIndex = -1;
        option.textContent = `--${value}`;
        option.addEventListener("click", (event) => {
          event.stopPropagation();
          setAlternativeParameter(form, flags, stateKey, value);
          form.dispatchEvent(new Event("change", { bubbles: true }));
          closeParameterSelectors({ restoreFocus: true });
        });
        menu.append(option);
      });
      trigger.addEventListener("click", (event) => {
        event.stopPropagation();
        if (selector.classList.contains("open")) closeParameterSelectors();
        else openParameterSelector(selector);
      });
      trigger.addEventListener("keydown", (event) => {
        if (!["ArrowDown", "ArrowUp", "Enter", " "].includes(event.key)) return;
        event.preventDefault();
        openParameterSelector(selector, true);
      });
      menu.addEventListener("keydown", (event) => navigateParameterSelector(event, selector));
      selector.append(trigger, menu);
      field.querySelector(".field-heading").replaceChildren(selector);
    });
  }
  setAlternativeParameter(
    form,
    flags,
    stateKey,
    importedType || form.dataset[stateKey] || flags[0],
  );
}

export function setAlternativeParameter(form, flags, stateKey, selected) {
  form.dataset[stateKey] = selected;
  flags.forEach((flag) => {
    const field = [...form.querySelectorAll(".form-field")]
      .find((candidate) => candidate.dataset.flag === flag);
    if (!field) return;
    const active = flag === selected;
    const control = field.querySelector(".argument-control");
    if (!active) clearArgumentControl(control);
    control.querySelectorAll("input, select, button").forEach((element) => {
      element.disabled = !active;
    });
    field.hidden = !active;
    const selector = field.querySelector(`[data-parameter-selector="${stateKey}"]`);
    selector.querySelector(".parameter-selector-value").textContent = `--${selected}`;
    selector.querySelectorAll(".parameter-selector-option").forEach((option) => {
      option.setAttribute("aria-selected", String(option.dataset.value === selected));
    });
  });
}

export function argumentControlValue(form, flag, useDefault = false) {
  const control = [...form.querySelectorAll(".argument-control")]
    .find((candidate) => candidate.dataset.flag === flag);
  if (!control) return "";
  if (control.dataset.boolean === "true") {
    return control.querySelector('input[type="checkbox"]')?.checked ? "true" : "";
  }
  const value = [...control.querySelectorAll("[data-value]")]
    .map((input) => input.value.trim())
    .find(Boolean);
  return value || (useDefault ? control.dataset.defaultValue : "") || "";
}

export function setConditionalFields(form, flags, visible, promote = false) {
  flags.forEach((flag) => {
    const field = [...form.querySelectorAll(".form-field")]
      .find((candidate) => candidate.dataset.flag === flag);
    if (!field) return;
    if (visible && promote) promoteParameterField(field);
    const control = field.querySelector(".argument-control");
    if (!visible) clearArgumentControl(control);
    control.querySelectorAll("input, select, button").forEach((element) => {
      element.disabled = !visible;
    });
    field.hidden = !visible;
  });
}

export function setArgumentFieldsEnabled(form, flags, enabled) {
  flags.forEach((flag) => {
    const field = [...form.querySelectorAll(".form-field")]
      .find((candidate) => candidate.dataset.flag === flag);
    if (!field) return;
    const control = field.querySelector(".argument-control");
    if (!enabled) clearArgumentControl(control);
    control.querySelectorAll("input, select, button").forEach((element) => {
      element.disabled = !enabled;
    });
  });
}

export function setArgumentFieldRequired(form, flag, required) {
  const field = [...form.querySelectorAll(".form-field")]
    .find((candidate) => candidate.dataset.flag === flag);
  if (!field) return;
  const control = field.querySelector(".argument-control");
  const effective = required || control.dataset.schemaRequired === "true";
  control.dataset.required = String(effective);
  control.querySelectorAll("[data-value]").forEach((input) => {
    input.setAttribute("aria-required", String(effective));
  });
  renderFieldHeading(field.querySelector(".field-heading"), control.dataset.flag, effective);
}

export function promoteParameterField(field) {
  const section = field.closest(".parameter-group");
  const primary = [...section.children]
    .find((child) => child.classList.contains("parameter-primary"));
  if (primary && field.parentElement !== primary) primary.append(field);
}

export function clearArgumentControl(control) {
  control.querySelectorAll('input[type="checkbox"]').forEach((input) => {
    input.checked = false;
  });
  control.querySelectorAll("[data-value]").forEach((input) => {
    input.value = "";
  });
  [...control.querySelectorAll(".value-row")].slice(1).forEach((row) => row.remove());
}

export function updateAdvancedVisibility(form) {
  form.querySelectorAll(".parameter-advanced").forEach((details) => {
    const visible = [...details.querySelectorAll(".form-field")].some((field) => !field.hidden);
    details.hidden = !visible;
    if (!visible) details.open = false;
  });
  form.querySelectorAll(".parameter-group").forEach((section) => {
    const fields = [...section.querySelectorAll(".form-field")];
    const visible = fields.some((field) => !field.hidden);
    const required = fields.some(
      (field) => field.querySelector(".argument-control")?.dataset.required === "true",
    );
    section.hidden = !visible && !required;
  });
}

export function setValueHint(element, flag) {
  const names = flag.value_names || [];
  if (!names.length) {
    setTranslatedText(element, "Optional switch");
    return;
  }
  const hint = names.map((name) => `<${name}>`).join(" ");
  if (flag.repeat) setTranslatedText(element, "{hint}; may be repeated", { hint });
  else element.textContent = hint;
}

export const NUMERIC_VALUE_FLAGS = new Set([
  "impute-tune-sigma",
  "filter-min-intensity",
  "filter-cv-threshold",
  "cpc",
]);

export function numericValue(flag) {
  const names = flag.value_names || [];
  return names.some((name) => ["N", "FRACTION", "CORRELATION"].includes(name))
    || (names.includes("VALUE") && NUMERIC_VALUE_FLAGS.has(flag.long || flag.id));
}

export function highlightAgentFields(fields) {
  window.clearTimeout(highlightAgentFields.timer);
  document.querySelectorAll(".form-field.agent-filled")
    .forEach((field) => field.classList.remove("agent-filled"));
  fields.forEach((field) => {
    field.closest(".parameter-advanced")?.setAttribute("open", "");
    field.classList.add("agent-filled");
  });
  fields.find((field) => !field.hidden)?.scrollIntoView({ block: "nearest" });
  highlightAgentFields.timer = window.setTimeout(() => {
    fields.forEach((field) => field.classList.remove("agent-filled"));
  }, 3200);
}

export function fillTemplateControl(control, flag, value) {
  if (control.dataset.boolean === "true") {
    control.querySelector('input[type="checkbox"]').checked = value;
    return;
  }
  value.forEach((values, index) => {
    if (index > 0) addValueRow(control, flag, true);
    const row = control.querySelectorAll(".value-row")[index];
    [...row.querySelectorAll("[data-value]")].forEach((input, valueIndex) => {
      input.value = values[valueIndex];
      if (!input.value) throw invalidTemplateValue(control.dataset.flag);
    });
  });
}

export function focusWorkflowParameter(parameter) {
  byId("command-review-dialog").close();
  const field = [...byId("command-form").querySelectorAll(".form-field")]
    .find((candidate) => candidate.dataset.flag === parameter);
  if (!field) return;
  field.closest(".parameter-advanced")?.setAttribute("open", "");
  field.classList.add("review-target");
  field.scrollIntoView({ behavior: "smooth", block: "center" });
  field.querySelector("input, select, button")?.focus();
  window.setTimeout(() => field.classList.remove("review-target"), 2400);
}
