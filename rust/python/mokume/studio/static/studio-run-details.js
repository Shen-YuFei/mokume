// Run details, provenance, run comparison, and the QC summary.
import {
  api,
  byId,
  formatBytes,
  formatDate,
  formatDuration,
  formatRunDuration,
  shortRunId,
  state,
  titleCase,
  toast,
} from "./studio-core.js";
import { translate, translatedElement } from "./studio-i18n.js";
import { applyWorkflowTemplate } from "./studio-workflow.js";
import {
  appendMetric,
  openCommandReview,
  workbenchSection,
  workflowStepsSection,
} from "./studio-review.js";
import { renderBottom, showBottomTab } from "./studio-runs.js";

export function appendQcMetric(list, label, value) {
  const card = document.createElement("div");
  const term = translatedElement("dt", "", label);
  const description = document.createElement("dd");
  description.textContent = value ?? "—";
  card.append(term, description);
  list.append(card);
}

export async function openRunDetails(runId) {
  const details = await api(`/api/runs/${encodeURIComponent(runId)}/details`);
  const index = state.runs.findIndex((run) => run.id === runId);
  const previousRecord = index >= 0 ? state.runs[index + 1] : null;
  let previous = null;
  if (previousRecord) {
    try {
      previous = await api(`/api/runs/${encodeURIComponent(previousRecord.id)}/details`);
    } catch (_) {
      previous = null;
    }
  }
  state.selectedRunDetails = { ...details, previous };
  renderRunDetails();
  const dialog = byId("run-details-dialog");
  if (!dialog.open) dialog.showModal();
}

export async function refreshOpenRunDetails(runId) {
  if (state.selectedRunDetails?.run?.id !== runId) return;
  const details = await api(`/api/runs/${encodeURIComponent(runId)}/details`);
  state.selectedRunDetails = { ...details, previous: state.selectedRunDetails.previous || null };
  renderRunDetails();
}

export function renderRunDetails() {
  const details = state.selectedRunDetails;
  const content = byId("run-details-content");
  content.replaceChildren();
  if (!details) return;
  byId("run-details-title").removeAttribute("data-i18n");
  byId("run-details-title").textContent = `${translate("Run details")} · ${shortRunId(details.run.id)}`;
  content.append(runOverview(details), executionStageSection(details), workflowStepsSection(
    details.planned_steps,
    "Scientific plan",
  ));
  const qc = workbenchSection("QC summary");
  const qcBody = document.createElement("div");
  renderQcPanel(qcBody, details.qc);
  qc.append(qcBody);
  content.append(qc, runArtifactSection(details), runProvenanceSection(details));
  if (details.previous) content.append(previousRunSection(details, details.previous));
  content.append(runLogsSection(details.logs));
  const selected = state.compareRunIds.has(details.run.id);
  byId("run-compare-choice").checked = selected;
  byId("run-load-parameters").disabled = !details.template;
  byId("run-rerun").disabled = !details.template || Boolean(state.activeRunId);
}

export function runOverview(details) {
  const wrapper = document.createElement("section");
  wrapper.className = "run-overview";
  const heading = document.createElement("div");
  heading.className = "run-overview-heading";
  const badge = translatedElement("span", `status-badge ${details.run.status}`, details.run.status);
  const workflow = document.createElement("strong");
  workflow.textContent = details.template?.workflow?.join(" ") || details.run.argv.join(" ");
  heading.append(badge, workflow);
  const command = document.createElement("code");
  command.textContent = details.command;
  const metrics = document.createElement("dl");
  metrics.className = "metric-grid compact";
  appendMetric(metrics, "Created", formatDate(details.run.created_at));
  appendMetric(metrics, "Started", formatDate(details.run.started_at));
  appendMetric(metrics, "Finished", formatDate(details.run.finished_at));
  appendMetric(metrics, "Duration", details.duration_seconds === null
    ? formatRunDuration(details.run)
    : formatDuration(details.duration_seconds));
  wrapper.append(heading, command, metrics);
  if (details.run.error) {
    const error = document.createElement("p");
    error.className = "run-error";
    error.textContent = details.run.error;
    wrapper.append(error);
  }
  return wrapper;
}

export function executionStageSection(details) {
  const section = workbenchSection("Execution progress");
  const list = document.createElement("ol");
  list.className = "stage-stepper";
  (details.stages || []).forEach((stage) => {
    const item = document.createElement("li");
    item.className = `stage-step ${stage.status}`;
    const marker = document.createElement("span");
    marker.className = "stage-marker";
    const text = document.createElement("span");
    const label = translatedElement("strong", "", stage.label);
    const meta = document.createElement("small");
    meta.textContent = `${translate(stage.status)} · ${formatDuration(stage.elapsed_seconds)}`;
    text.append(label, meta);
    if (stage.error) {
      const error = document.createElement("small");
      error.className = "stage-error";
      error.textContent = stage.error;
      text.append(error);
    }
    item.append(marker, text);
    list.append(item);
  });
  section.append(list);
  return section;
}

export function runArtifactSection(details) {
  const section = workbenchSection("Run artifacts");
  const list = document.createElement("div");
  list.className = "detail-file-list";
  if (!details.artifacts?.length) {
    list.append(translatedElement("p", "empty-state", "No artifacts were recorded."));
  }
  (details.artifacts || []).forEach((artifact) => {
    const link = document.createElement("a");
    link.href = `/api/artifacts/${encodeURIComponent(artifact.id)}`;
    link.target = "_blank";
    link.rel = "noopener";
    const path = document.createElement("code");
    path.textContent = artifact.path;
    const meta = document.createElement("small");
    meta.textContent = `${formatBytes(artifact.size)} · sha256 ${artifact.sha256}`;
    link.append(path, meta);
    list.append(link);
  });
  section.append(list);
  return section;
}

export function runProvenanceSection(details) {
  const section = workbenchSection("Parameters & provenance");
  const provenance = details.provenance || {};
  const software = document.createElement("dl");
  software.className = "metric-grid compact";
  appendMetric(software, "Software", provenance.mokume_version ? `Mokume ${provenance.mokume_version}` : "—");
  appendMetric(software, "Python", provenance.python_version || "—");
  appendMetric(software, "Platform", provenance.platform || "—");
  appendMetric(software, "Canonical command", details.command);
  section.append(software);
  section.append(hashList("Input hashes", provenance.inputs || [], "input"));
  section.append(hashList("Artifact hashes", provenance.artifacts || details.artifacts || [], "artifact"));
  section.append(jsonDisclosure("Standard parameters", details.parameters || {}));
  return section;
}

export function hashList(title, records, kind) {
  const block = document.createElement("div");
  block.className = "hash-block";
  block.append(translatedElement("h4", "", title));
  const rows = [];
  records.forEach((record) => {
    if (record.sha256) rows.push([record.path || record.resolved_path, record.sha256]);
    (record.entries || []).forEach((entry) => {
      if (entry.sha256) rows.push([entry.path || entry.resolved_path, entry.sha256]);
    });
  });
  if (!rows.length) {
    block.append(translatedElement("p", "empty-state", kind === "input"
      ? "No input hashes were recorded."
      : "No artifacts were recorded."));
    return block;
  }
  rows.forEach(([path, hash]) => {
    const row = document.createElement("div");
    row.className = "hash-row";
    const name = document.createElement("code");
    name.textContent = path || "—";
    const digest = document.createElement("code");
    digest.textContent = hash;
    row.append(name, digest);
    block.append(row);
  });
  return block;
}

export function jsonDisclosure(title, payload) {
  const details = document.createElement("details");
  details.className = "json-disclosure";
  const summary = translatedElement("summary", "", title);
  const content = document.createElement("pre");
  content.textContent = JSON.stringify(payload, null, 2);
  details.append(summary, content);
  return details;
}

export function runLogsSection(logs) {
  const section = workbenchSection("Run logs");
  const entries = ["stdout", "stderr"].filter((name) => logs?.[name]?.length);
  if (!entries.length) {
    section.append(translatedElement("p", "empty-state", "No logs were written."));
    return section;
  }
  entries.forEach((name) => {
    const details = document.createElement("details");
    details.className = "log-disclosure";
    if (name === "stderr") details.open = true;
    const summary = translatedElement("summary", "", name);
    const content = document.createElement("pre");
    content.textContent = logs[name].join("\n");
    details.append(summary, content);
    section.append(details);
  });
  return section;
}

export function comparisonValues(details) {
  const values = new Map();
  values.set("$workflow", details.template?.workflow || details.run.argv.slice(0, 2));
  Object.entries(details.template?.parameters || {}).forEach(([name, value]) => {
    values.set(`--${name}`, value);
  });
  values.set("$mokume", details.provenance?.mokume_version || null);
  values.set("$python", details.provenance?.python_version || null);
  values.set("$inputs", (details.provenance?.inputs || []).map((item) => ({
    path: item.path,
    sha256: item.sha256,
    entries: item.entries?.map((entry) => ({ path: entry.path, sha256: entry.sha256 })),
  })));
  values.set("$outputs", (details.artifacts || []).map((item) => ({
    path: item.path,
    sha256: item.sha256,
  })));
  return values;
}

export function comparisonRows(first, second) {
  const left = comparisonValues(first);
  const right = comparisonValues(second);
  return [...new Set([...left.keys(), ...right.keys()])].sort().flatMap((name) => {
    const leftValue = left.has(name) ? left.get(name) : null;
    const rightValue = right.has(name) ? right.get(name) : null;
    return JSON.stringify(leftValue) === JSON.stringify(rightValue)
      ? []
      : [{ name, left: leftValue, right: rightValue }];
  });
}

export function comparisonLabel(name) {
  const labels = {
    $workflow: "Workflow",
    $mokume: "Software",
    $python: "Python",
    $inputs: "Input files",
    $outputs: "Output files",
  };
  return translate(labels[name] || name);
}

export function comparisonValue(value) {
  if (value === null || value === undefined) return "—";
  return typeof value === "string" ? value : JSON.stringify(value, null, 2);
}

export function comparisonTable(rows, leftLabel, rightLabel) {
  if (!rows.length) return translatedElement("p", "empty-state", "No parameter differences.");
  const table = document.createElement("table");
  table.className = "workbench-table comparison-table";
  const head = document.createElement("thead");
  const header = document.createElement("tr");
  ["Parameter", leftLabel, rightLabel].forEach((label, index) => {
    const cell = document.createElement("th");
    cell.textContent = index === 0 ? translate(label) : label;
    header.append(cell);
  });
  head.append(header);
  const body = document.createElement("tbody");
  rows.forEach((row) => {
    const line = document.createElement("tr");
    const name = document.createElement("th");
    name.scope = "row";
    name.textContent = comparisonLabel(row.name);
    const left = document.createElement("td");
    const right = document.createElement("td");
    left.textContent = comparisonValue(row.left);
    right.textContent = comparisonValue(row.right);
    line.append(name, left, right);
    body.append(line);
  });
  table.append(head, body);
  return table;
}

export function previousRunSection(details, previous) {
  const section = workbenchSection("Previous run changes");
  const rows = comparisonRows(previous, details);
  if (!rows.length) {
    section.append(translatedElement("p", "empty-state", "No changes from the previous run."));
    return section;
  }
  section.append(comparisonTable(rows, shortRunId(previous.run.id), shortRunId(details.run.id)));
  return section;
}

export function loadSelectedRunParameters() {
  const details = state.selectedRunDetails;
  if (!details?.template) return;
  applyWorkflowTemplate(details.template);
  byId("run-details-dialog").close();
  toast(translate("Parameters loaded"));
}

export async function copySelectedRunCommand() {
  const command = state.selectedRunDetails?.command;
  if (!command) return;
  await navigator.clipboard.writeText(command);
  toast(translate("Command copied"));
}

export async function prepareSelectedRerun() {
  const details = state.selectedRunDetails;
  if (!details?.template) return;
  applyWorkflowTemplate(details.template);
  byId("run-details-dialog").close();
  toast(translate("Parameters loaded. Review output paths before running again."));
  await openCommandReview("run");
}

export function updateSelectedRunComparison() {
  const runId = state.selectedRunDetails?.run?.id;
  if (!runId) return;
  if (byId("run-compare-choice").checked) {
    if (state.compareRunIds.size >= 2 && !state.compareRunIds.has(runId)) {
      byId("run-compare-choice").checked = false;
      toast(translate("Select exactly two runs to compare"), true);
      return;
    }
    state.compareRunIds.add(runId);
  } else {
    state.compareRunIds.delete(runId);
  }
  if (state.bottomTab === "runs") renderBottom();
}

export async function openRunComparison() {
  if (state.compareRunIds.size !== 2) {
    throw new Error(translate("Select exactly two runs to compare"));
  }
  state.runComparisonDetails = await Promise.all([...state.compareRunIds].map((runId) =>
    api(`/api/runs/${encodeURIComponent(runId)}/details`),
  ));
  renderRunComparison();
  byId("run-compare-dialog").showModal();
}

export function renderRunComparison() {
  const content = byId("run-compare-content");
  content.replaceChildren();
  if (!state.runComparisonDetails?.length) return;
  const [first, second] = state.runComparisonDetails;
  content.append(translatedElement("p", "estimate-note", "Comparison includes canonical workflow parameters."));
  content.append(comparisonTable(
    comparisonRows(first, second),
    `${translate("Run A")} · ${shortRunId(first.run.id)}`,
    `${translate("Run B")} · ${shortRunId(second.run.id)}`,
  ));
}

export async function loadLatestQc() {
  const succeeded = state.runs.find((run) => run.status === "succeeded");
  if (!succeeded) {
    state.qcDetails = null;
    if (state.bottomTab === "qc") renderBottom();
    return;
  }
  state.qcDetails = await api(`/api/runs/${encodeURIComponent(succeeded.id)}/details`);
  if (state.bottomTab === "qc") renderBottom();
}

export async function showCompletedRunQc(runId) {
  state.qcDetails = await api(`/api/runs/${encodeURIComponent(runId)}/details`);
  showBottomTab("qc");
}

export function renderQcPanel(container, qc) {
  container.replaceChildren();
  container.classList.add("qc-panel");
  if (!qc?.available) {
    const reason = qc?.reason || "No completed run with QC is selected.";
    container.append(translatedElement("p", "empty-state", reason));
    return;
  }
  const metrics = document.createElement("dl");
  metrics.className = "qc-metrics";
  appendQcMetric(metrics, "Samples", String(qc.sample_count));
  Object.entries(qc.entity_counts || { [qc.entity_label]: qc.entity_count }).forEach(([label, count]) => {
    appendQcMetric(metrics, titleCase(label), String(count));
  });
  appendQcMetric(metrics, "Missing values", `${qc.missing_percent}%`);
  appendQcMetric(metrics, "Median CV", qc.median_cv_percent === null ? "—" : `${qc.median_cv_percent}%`);
  const quartiles = qc.log2_intensity_quartiles || [];
  appendQcMetric(metrics, "Log2 intensity", quartiles.length === 3
    ? `Q1 ${quartiles[0]} · median ${quartiles[1]} · Q3 ${quartiles[2]}`
    : "—");
  container.append(metrics);

  const diagnostics = document.createElement("div");
  diagnostics.className = "qc-diagnostics";
  diagnostics.append(qcDiagnostic(
    "Sample correlation",
    qc.correlation?.median === null ? "—" : String(qc.correlation?.median),
    [
      `${translate("Lowest-correlation sample")}: ${qc.correlation?.lowest_sample || "—"}`,
      outlierText(qc.correlation?.outliers),
    ],
  ));
  const variance = qc.pca?.variance_percent || [];
  diagnostics.append(qcDiagnostic(
    "PCA variance",
    variance.length ? variance.map((value) => `${value}%`).join(" / ") : "—",
    [`${translate("PC1 / PC2")}`, outlierText(qc.pca?.outliers)],
  ));
  diagnostics.append(normalizationCard(qc.normalization));
  container.append(diagnostics);

  const footer = document.createElement("div");
  footer.className = "qc-footer";
  const artifact = document.createElement("code");
  artifact.textContent = `${translate("Primary QC artifact")}: ${qc.artifact_path}`;
  footer.append(artifact);
  (qc.reports || []).forEach((report) => {
    const link = translatedElement("a", "", "Open detailed report");
    link.href = `/api/artifacts/${encodeURIComponent(report.id)}`;
    link.target = "_blank";
    link.rel = "noopener";
    footer.append(link);
  });
  container.append(footer);
}

export function qcDiagnostic(title, value, notes) {
  const card = document.createElement("article");
  card.className = "qc-diagnostic";
  card.append(translatedElement("h4", "", title));
  const metric = document.createElement("strong");
  metric.textContent = value;
  card.append(metric);
  notes.filter(Boolean).forEach((note) => {
    const line = document.createElement("small");
    line.textContent = note;
    card.append(line);
  });
  return card;
}

export function outlierText(outliers) {
  return outliers?.length
    ? `${translate("Potential outliers")}: ${outliers.join(", ")}`
    : translate("No outliers detected");
}

export function normalizationCard(normalization) {
  const methods = Object.entries(normalization?.methods || {});
  const card = qcDiagnostic(
    "Normalization",
    methods.length ? methods.map(([name, value]) => `${name}: ${value}`).join(" · ") : "—",
    [],
  );
  const note = document.createElement("small");
  note.className = "normalization-note";
  note.textContent = translate(normalization?.message || "Normalization comparison unavailable");
  card.append(note);
  return card;
}
