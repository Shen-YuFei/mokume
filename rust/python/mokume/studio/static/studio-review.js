// Preflight review dialog and submission of a reviewed analysis.
import { api, byId, formatBytes, originHeaders, state, titleCase, toast } from "./studio-core.js";
import { setTranslatedText, translate, translatedElement } from "./studio-i18n.js";
import { focusWorkflowParameter } from "./studio-form.js";
import { commandArgv } from "./studio-workflow.js";
import { refreshRuns, showBottomTab } from "./studio-runs.js";

export async function requestCommandReview() {
  return api("/api/commands/review", {
    method: "POST",
    headers: originHeaders(),
    body: JSON.stringify({ argv: commandArgv() }),
  });
}

export async function openCommandReview(intent = "validate") {
  const review = await requestCommandReview();
  state.commandReview = { ...review, intent };
  renderCommandReview();
  byId("command-review-dialog").showModal();
}

export async function validateCommand() {
  await openCommandReview("validate");
}

export async function runCommand() {
  await openCommandReview("run");
}

export async function queueCommand() {
  await openCommandReview("queue");
}

export function workbenchSection(title) {
  const section = document.createElement("section");
  section.className = "workbench-section";
  const heading = document.createElement("h3");
  setTranslatedText(heading, title);
  section.append(heading);
  return section;
}

export function renderCommandReview() {
  const review = state.commandReview;
  const content = byId("command-review-content");
  content.replaceChildren();
  if (!review) return;
  content.append(
    reviewHeadline(review),
    reviewChecksSection(review),
    reviewPathsSection(review),
    reviewResourcesSection(review),
    workflowStepsSection(review.planned_steps),
  );
  if (review.sdrf) content.append(sdrfSection(review.sdrf));
  const canSubmit = review.valid && !state.activeRunId;
  byId("review-run").disabled = !canSubmit;
  byId("review-queue").disabled = !review.valid;
}

export function reviewHeadline(review) {
  const headline = document.createElement("div");
  headline.className = `review-headline ${review.valid ? "ok" : "error"}`;
  const title = translatedElement("strong", "", review.valid
    ? "Ready to run"
    : "Resolve the highlighted checks before running");
  const command = document.createElement("code");
  command.textContent = review.command;
  headline.append(title, command);
  return headline;
}

export function reviewChecksSection(review) {
  const section = workbenchSection("Status");
  const list = document.createElement("div");
  list.className = "review-checks";
  review.checks.forEach((check) => {
    const row = document.createElement(check.parameter ? "button" : "div");
    row.className = `review-check ${check.status}`;
    if (check.parameter) {
      row.type = "button";
      row.addEventListener("click", () => focusWorkflowParameter(check.parameter));
    }
    const marker = document.createElement("span");
    marker.className = "review-check-marker";
    const text = document.createElement("span");
    const label = document.createElement("strong");
    label.textContent = translate(check.label);
    const message = document.createElement("small");
    message.textContent = translate(check.message);
    text.append(label, message);
    row.append(marker, text);
    list.append(row);
  });
  section.append(list);
  return section;
}

export function reviewPathsSection(review) {
  const section = workbenchSection("Inputs & outputs");
  const list = document.createElement("div");
  list.className = "review-paths";
  [...review.inputs.map((item) => ["Input", item]), ...review.outputs.map((item) => ["Output", item])]
    .forEach(([kind, item]) => {
      const row = document.createElement("button");
      row.type = "button";
      row.className = `review-path ${item.status}`;
      row.addEventListener("click", () => focusWorkflowParameter(item.parameter));
      const type = translatedElement("span", "review-path-kind", kind);
      const body = document.createElement("span");
      const option = document.createElement("strong");
      option.textContent = `--${item.parameter}`;
      const path = document.createElement("code");
      path.textContent = item.path;
      body.append(option, path);
      const size = document.createElement("span");
      size.textContent = item.size_bytes === null ? "—" : formatBytes(item.size_bytes);
      row.append(type, body, size);
      list.append(row);
    });
  if (!list.children.length) list.append(translatedElement("p", "empty-state", "No artifacts yet."));
  section.append(list);
  return section;
}

export function reviewResourcesSection(review) {
  const section = workbenchSection("Resource estimate");
  const grid = document.createElement("dl");
  grid.className = "metric-grid compact";
  appendMetric(grid, "Threads", String(review.resources.threads));
  appendMetric(grid, "Configured memory", review.resources.configured_memory || "—");
  appendMetric(grid, "Suggested peak memory", formatBytes(review.resources.suggested_peak_memory_bytes));
  appendMetric(grid, "Available memory", formatBytes(review.resources.available_memory_bytes));
  appendMetric(grid, "Suggested free disk", formatBytes(review.resources.suggested_free_disk_bytes));
  appendMetric(grid, "Available disk", formatBytes(review.resources.free_disk_bytes));
  section.append(grid, translatedElement("small", "estimate-note", "Estimate only"));
  return section;
}

export function workflowStepsSection(steps, title = "Planned workflow") {
  const section = workbenchSection(title);
  const list = document.createElement("ol");
  list.className = "planned-steps";
  (steps || []).forEach((step) => list.append(translatedElement("li", "", step)));
  section.append(list);
  return section;
}

export function sdrfSection(sdrf) {
  const section = workbenchSection("SDRF sample mapping");
  if (!sdrf.rows?.length) {
    const message = document.createElement("p");
    message.className = "empty-state";
    message.textContent = sdrf.issues?.join("; ") || translate("No SDRF was selected for this workflow.");
    section.append(message);
    return section;
  }
  const mapping = document.createElement("div");
  mapping.className = "sdrf-column-map";
  Object.entries(sdrf.columns || {}).forEach(([name, column]) => {
    const chip = document.createElement("span");
    chip.textContent = `${translate(titleCase(name.replaceAll("_", " ")))}: ${column || "—"}`;
    mapping.append(chip);
  });
  const table = document.createElement("table");
  table.className = "workbench-table";
  const headers = ["Sample", "Condition", "Batch", "Replicate", "Plex", "Label", "Reference"];
  const head = document.createElement("thead");
  const headRow = document.createElement("tr");
  headers.forEach((label) => headRow.append(translatedElement("th", "", label)));
  head.append(headRow);
  const body = document.createElement("tbody");
  sdrf.rows.forEach((row) => {
    const line = document.createElement("tr");
    [row.sample, row.condition, row.batch, row.replicate, row.plex, row.label,
      translate(row.reference ? "Yes" : "No")].forEach((value) => {
      const cell = document.createElement("td");
      cell.textContent = value || "—";
      line.append(cell);
    });
    body.append(line);
  });
  table.append(head, body);
  const note = translatedElement("small", "estimate-note", "Showing {shown} of {total} rows", {
    shown: sdrf.rows.length,
    total: sdrf.row_count,
  });
  section.append(mapping, table, note);
  if (sdrf.issues?.length) {
    const issues = document.createElement("ul");
    issues.className = "review-issues";
    sdrf.issues.forEach((issue) => {
      const item = document.createElement("li");
      item.textContent = issue;
      issues.append(item);
    });
    section.append(issues);
  }
  return section;
}

export async function submitReviewedRun() {
  const review = state.commandReview;
  if (!review?.valid) return;
  const payload = await api("/api/runs", {
    method: "POST",
    headers: originHeaders(),
    body: JSON.stringify({ argv: review.argv }),
  });
  byId("command-review-dialog").close();
  state.activeRunId = payload.id;
  toast(translate("Run {id} queued", { id: payload.id }));
  await refreshRuns();
  showBottomTab("runs");
}

export function appendMetric(list, label, value) {
  const term = translatedElement("dt", "", label);
  const description = document.createElement("dd");
  description.textContent = value ?? "—";
  list.append(term, description);
}
