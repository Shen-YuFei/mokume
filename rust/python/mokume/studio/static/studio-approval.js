// Approval of Agent-proposed analyses, including approvals restored after a refresh.
import { api, byId, originHeaders, state } from "./studio-core.js";
import { setTranslatedText, translate } from "./studio-i18n.js";
import { updateActionStates } from "./studio-menus.js";
import { refreshRuns, showBottomTab } from "./studio-runs.js";
import { newAgentBody, parseToolArguments, runAgentRequest } from "./studio-assistant.js";

export async function presentApproval(interrupt, stream, mode) {
  const toolId = interrupt.toolCallId;
  const tool = stream.tools.get(toolId);
  const args = parseToolArguments(tool || {});
  if (tool?.name !== "run_approved_evaluation" || !args.approval_id || !args.payload_hash) {
    throw new Error(translate("The analysis approval did not contain a valid server plan"));
  }
  const record = await api(`/api/approvals/${encodeURIComponent(args.approval_id)}`);
  if (record.payload_hash !== args.payload_hash) throw new Error(translate("Approval hash mismatch"));
  state.pendingApproval = {
    record,
    interruptId: interrupt.id,
    mode,
    datasetId: state.dataset?.id || null,
    decision: null,
  };
  updateActionStates();
  persistPendingApproval();
  renderApproval(record);
  byId("approval-dialog").showModal();
}

export function pendingApprovalKey() {
  return state.projectId ? `mokume:approval:${state.projectId}` : null;
}

export function persistPendingApproval() {
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

export function clearPendingApproval() {
  const key = pendingApprovalKey();
  if (key) sessionStorage.removeItem(key);
  state.pendingApproval = null;
  updateActionStates();
}

export async function restorePendingApproval() {
  const key = pendingApprovalKey();
  if (!key) return;
  const saved = JSON.parse(sessionStorage.getItem(key) || "null");
  if (!saved?.approvalId || !saved?.interruptId || !saved?.mode) return;
  if (saved.mode !== "agent") {
    sessionStorage.removeItem(key);
    return;
  }
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
  updateActionStates();
  renderApproval(record);
  byId("approval-dialog").showModal();
}

export function renderApproval(record) {
  const content = byId("approval-content");
  content.replaceChildren();
  const section = document.createElement("section");
  section.className = "approval-card";
  const title = document.createElement("div");
  title.className = "approval-card-title";
  const label = document.createElement("span");
  setTranslatedText(label, "Canonical parameters");
  const status = document.createElement("span");
  status.className = "approval-status";
  status.textContent = translate(record.status);
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
  setTranslatedText(byId("approval-approve"), decided === true ? "Resume Approved Run" : "Approve and Run");
  setTranslatedText(byId("approval-reject"), decided === false ? "Resume Rejection" : "Reject");
}

export async function decideApproval(approved) {
  const pending = state.pendingApproval;
  if (!pending) throw new Error(translate("No analysis plan is awaiting approval"));
  if (pending.decision !== null && pending.decision !== approved) {
    throw new Error(translate("This analysis plan already has the opposite decision"));
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
