// Assistant conversations: streaming replies, reasoning, Markdown, and stored threads.
import {
  agentHeaders,
  api,
  byId,
  originHeaders,
  reportError,
  state,
  toast,
} from "./studio-core.js";
import { setTranslatedAttribute, setTranslatedText, translate } from "./studio-i18n.js";
import {
  closeAssistantModeMenu,
  selectAssistantMode,
  updateActionStates,
} from "./studio-menus.js";
import { agentWorkflowFormState, applyWorkflowParameterPatch } from "./studio-workflow.js";
import { presentApproval } from "./studio-approval.js";

export const CONVERSATION_TITLE_MAX_LENGTH = 120;

export const CONVERSATION_DELETE_CONFIRM_MS = 4000;

export function threadsForProject(projectId) {
  if (!projectId) return { ask: crypto.randomUUID(), agent: crypto.randomUUID() };
  return Object.fromEntries(["ask", "agent"].map((mode) => {
    const key = `mokume:thread:${projectId}:${mode}`;
    let value = sessionStorage.getItem(key);
    if (!value) {
      value = crypto.randomUUID();
      sessionStorage.setItem(key, value);
    }
    return [mode, value];
  }));
}

export function resetAssistantConversation() {
  byId("assistant-messages").querySelectorAll(".assistant-message").forEach((message) => message.remove());
}

export function startNewAssistantChat() {
  if (state.agentBusy || state.pendingApproval) return;
  state.threads = { ask: crypto.randomUUID(), agent: crypto.randomUUID() };
  if (state.projectId) {
    Object.entries(state.threads).forEach(([mode, threadId]) => {
      sessionStorage.setItem(`mokume:thread:${state.projectId}:${mode}`, threadId);
    });
  }
  resetAssistantConversation();
  const input = byId("assistant-input");
  input.value = "";
  resizeAssistantInput();
  input.focus();
}

export function conversationTimestamp(value) {
  const locale = state.language === "zh-CN" ? "zh-CN" : "en";
  return new Date(value).toLocaleString(locale, { dateStyle: "medium", timeStyle: "short" });
}

export function renderConversationHistory(threads, workspace) {
  const list = byId("conversation-list");
  list.replaceChildren();
  if (!threads.length) {
    const empty = document.createElement("div");
    empty.className = "empty-state";
    setTranslatedText(empty, "No conversations in this workspace.");
    list.append(empty);
    return;
  }
  threads.forEach((thread) => {
    const row = document.createElement("div");
    row.className = "conversation-row";
    const openButton = document.createElement("button");
    openButton.className = "conversation-open";
    openButton.type = "button";
    const title = document.createElement("span");
    title.className = "conversation-title";
    title.textContent = thread.title;
    const meta = document.createElement("span");
    meta.className = "conversation-meta";
    meta.textContent = `${translate(thread.mode === "agent" ? "Agent" : "Ask")} · ${conversationTimestamp(thread.updated_at)}`;
    const owner = document.createElement("span");
    owner.className = "conversation-workspace";
    setTranslatedText(owner, "Workspace: {path}", { path: workspace.root });
    owner.title = workspace.root;
    openButton.append(title, meta, owner);
    openButton.addEventListener("click", () => openStoredConversation(thread).catch(reportError));

    const renameButton = document.createElement("button");
    renameButton.className = "conversation-rename";
    renameButton.type = "button";
    setTranslatedAttribute(renameButton, "aria-label", "Rename conversation: {title}", { title: thread.title });
    setTranslatedAttribute(renameButton, "title", "Rename conversation: {title}", { title: thread.title });
    const renameIcon = document.createElementNS("http://www.w3.org/2000/svg", "svg");
    renameIcon.setAttribute("viewBox", "0 0 24 24");
    renameIcon.setAttribute("aria-hidden", "true");
    const renamePath = document.createElementNS("http://www.w3.org/2000/svg", "path");
    renamePath.setAttribute("d", "M12 20h9M16.5 3.5a2.1 2.1 0 0 1 3 3L8 18l-4 1 1-4Z");
    renameIcon.append(renamePath);
    renameButton.append(renameIcon);
    renameButton.addEventListener("click", () => beginConversationRename(row, thread));

    const deleteButton = document.createElement("button");
    deleteButton.className = "conversation-delete";
    deleteButton.type = "button";
    setTranslatedAttribute(deleteButton, "aria-label", "Delete conversation: {title}", { title: thread.title });
    setTranslatedAttribute(deleteButton, "title", "Delete conversation: {title}", { title: thread.title });
    const icon = document.createElementNS("http://www.w3.org/2000/svg", "svg");
    icon.setAttribute("viewBox", "0 0 24 24");
    icon.setAttribute("aria-hidden", "true");
    const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
    path.setAttribute("d", "M4 7h16M9 7V4h6v3m-8 0 1 13h8l1-13M10 11v5M14 11v5");
    icon.append(path);
    deleteButton.append(icon);
    deleteButton.addEventListener("click", () => {
      if (deleteButton.classList.contains("confirming")) deleteStoredConversation(thread).catch(reportError);
      else armConversationDelete(deleteButton);
    });
    row.append(openButton, renameButton, deleteButton);
    list.append(row);
  });
}

export async function refreshConversationHistory() {
  const payload = await api("/api/agent/threads");
  if (payload.project_id !== state.projectId) {
    throw new Error("Conversation history belongs to a different workspace");
  }
  renderConversationHistory(payload.threads, payload.workspace);
}

export async function openConversationHistory() {
  if (!state.projectId) throw new Error(translate("Open a folder first"));
  await refreshConversationHistory();
  byId("conversation-dialog").showModal();
}

export async function openStoredConversation(summary) {
  const thread = await api(`/api/agent/threads/${encodeURIComponent(summary.id)}?mode=${encodeURIComponent(summary.mode)}`);
  if (thread.project_id !== state.projectId) {
    throw new Error("Conversation belongs to a different workspace");
  }
  state.threads[thread.mode] = thread.id;
  sessionStorage.setItem(`mokume:thread:${state.projectId}:${thread.mode}`, thread.id);
  byId("conversation-dialog").close();
  selectAssistantMode(thread.mode);
  resetAssistantConversation();
  thread.conversation.forEach((message) => {
    if (message.role === "reasoning") {
      appendAssistantReasoning(message.text);
    } else if (["user", "assistant"].includes(message.role)) {
      appendAssistantMessage(message.role, message.text);
    }
  });
}

// Edit the title in place; Enter or leaving the field saves, Escape cancels.
export function beginConversationRename(row, summary) {
  const openButton = row.querySelector(".conversation-open");
  const input = document.createElement("input");
  input.className = "conversation-title-input";
  input.type = "text";
  input.value = summary.title;
  input.maxLength = CONVERSATION_TITLE_MAX_LENGTH;
  setTranslatedAttribute(input, "aria-label", "Rename conversation");
  openButton.hidden = true;
  row.prepend(input);
  input.focus();
  input.select();
  let finished = false;
  const finish = (save) => {
    if (finished) return;
    finished = true;
    const title = input.value.trim();
    const restore = () => {
      input.remove();
      openButton.hidden = false;
    };
    if (!save || !title || title === summary.title) {
      restore();
      return;
    }
    renameStoredConversation(summary, title).catch((error) => {
      restore();
      reportError(error);
    });
  };
  input.addEventListener("keydown", (event) => {
    if (event.key === "Enter") {
      event.preventDefault();
      finish(true);
    } else if (event.key === "Escape") {
      event.preventDefault();
      event.stopPropagation();
      finish(false);
      row.querySelector(".conversation-rename")?.focus();
    }
  });
  input.addEventListener("blur", () => finish(true));
}

// The first click arms the button; a second click within a few seconds deletes.
export function armConversationDelete(button) {
  const labelKey = button.getAttribute("data-i18n-aria-label");
  const labelValues = JSON.parse(button.getAttribute("data-i18n-aria-label-values") || "{}");
  const text = document.createElement("span");
  setTranslatedText(text, "Delete?");
  button.append(text);
  button.classList.add("confirming");
  setTranslatedAttribute(button, "aria-label", "Confirm delete");
  const reset = () => {
    window.clearTimeout(timer);
    if (!button.isConnected || !button.classList.contains("confirming")) return;
    button.classList.remove("confirming");
    text.remove();
    setTranslatedAttribute(button, "aria-label", labelKey, labelValues);
  };
  const timer = window.setTimeout(reset, CONVERSATION_DELETE_CONFIRM_MS);
  button.addEventListener("blur", reset, { once: true });
}

export async function renameStoredConversation(summary, title) {
  await api(`/api/agent/threads/${encodeURIComponent(summary.id)}?mode=${encodeURIComponent(summary.mode)}`, {
    method: "PATCH",
    headers: originHeaders(),
    body: JSON.stringify({ title }),
  });
  await refreshConversationHistory();
  toast(translate("Conversation renamed"));
}

export async function deleteStoredConversation(summary) {
  await api(`/api/agent/threads/${encodeURIComponent(summary.id)}?mode=${encodeURIComponent(summary.mode)}`, {
    method: "DELETE",
    headers: originHeaders(),
  });
  if (state.threads[summary.mode] === summary.id) {
    const threadId = crypto.randomUUID();
    state.threads[summary.mode] = threadId;
    sessionStorage.setItem(`mokume:thread:${state.projectId}:${summary.mode}`, threadId);
    if (byId("assistant-mode").value === summary.mode) resetAssistantConversation();
  }
  await refreshConversationHistory();
  toast(translate("Conversation deleted"));
}

export function appendAssistantMessage(role, text) {
  const messages = byId("assistant-messages");
  const article = document.createElement("article");
  article.className = `assistant-message ${role}`;
  const body = document.createElement("div");
  body.className = "message-body";
  const paragraph = document.createElement("p");
  paragraph.textContent = text;
  body.append(paragraph);
  article.append(body);
  messages.append(article);
  if (role === "assistant" && text) void renderAssistantMarkdown(body, text);
  messages.scrollTop = messages.scrollHeight;
  return paragraph;
}

export function appendAssistantThinking() {
  const messages = byId("assistant-messages");
  const article = document.createElement("article");
  article.className = "assistant-message assistant assistant-thinking";
  article.setAttribute("role", "status");
  setTranslatedAttribute(article, "aria-label", "Thinking…");
  const body = document.createElement("div");
  body.className = "message-body";
  const label = document.createElement("span");
  label.className = "assistant-thinking-label";
  setTranslatedText(label, "Thinking");
  const dots = document.createElement("span");
  dots.className = "assistant-thinking-dots";
  dots.setAttribute("aria-hidden", "true");
  for (let index = 0; index < 3; index += 1) {
    const dot = document.createElement("span");
    dot.className = "assistant-thinking-dot";
    dots.append(dot);
  }
  body.append(label, dots);
  article.append(body);
  messages.append(article);
  messages.scrollTop = messages.scrollHeight;
  return article;
}

export function clearAssistantThinking(stream) {
  stream.thinking?.remove();
  stream.thinking = null;
}

export function appendAssistantReasoning(text = "", active = false) {
  const messages = byId("assistant-messages");
  const article = document.createElement("article");
  article.className = "assistant-message assistant assistant-reasoning";
  const details = document.createElement("details");
  details.className = "assistant-reasoning-details";
  details.classList.toggle("active", active);
  const summary = document.createElement("summary");
  const label = document.createElement("span");
  setTranslatedText(label, active ? "Thinking…" : "Thought process");
  const content = document.createElement("div");
  content.className = "assistant-reasoning-content";
  content.textContent = text;
  summary.append(label);
  details.append(summary, content);
  article.append(details);
  messages.append(article);
  messages.scrollTop = messages.scrollHeight;
  return { details, label, content };
}

export function adjacentAssistantReasoning() {
  const article = byId("assistant-messages").lastElementChild;
  if (!article?.classList.contains("assistant-reasoning")) return null;
  const details = article.querySelector(".assistant-reasoning-details");
  const label = details?.querySelector("summary span");
  const content = details?.querySelector(".assistant-reasoning-content");
  return details && label && content ? { details, label, content } : null;
}

export function reasoningMessage(stream, event) {
  if (!stream.reasoning) {
    stream.reasoning = adjacentAssistantReasoning()
      || appendAssistantReasoning("", true);
  }
  if (event.messageId && !stream.reasoningMessageIds.has(event.messageId)) {
    if (stream.reasoning.content.textContent
        && !stream.reasoning.content.textContent.endsWith("\n")) {
      stream.reasoning.content.textContent += "\n";
    }
    stream.reasoningMessageIds.add(event.messageId);
  }
  stream.reasoning.details.classList.add("active");
  setTranslatedText(stream.reasoning.label, "Thinking…");
  return stream.reasoning;
}

export function finishReasoning(stream) {
  const reasoning = stream.reasoning;
  if (!reasoning) return;
  reasoning.details.classList.remove("active");
  setTranslatedText(reasoning.label, "Thought process");
}

export function appendAssistantError(error) {
  const messages = byId("assistant-messages");
  const article = document.createElement("article");
  article.className = "assistant-message assistant assistant-error";
  article.setAttribute("role", "alert");
  const body = document.createElement("div");
  body.className = "message-body";
  const title = document.createElement("strong");
  title.className = "assistant-error-title";
  setTranslatedText(title, "Request failed");
  const detail = document.createElement("p");
  detail.className = "assistant-error-detail";
  detail.textContent = error?.message || String(error);
  body.append(title, detail);
  article.append(body);
  messages.append(article);
  messages.scrollTop = messages.scrollHeight;
}

export async function renderAssistantMarkdown(body, text) {
  try {
    const rendered = await api("/api/agent/markdown", {
      method: "POST",
      headers: originHeaders(),
      body: JSON.stringify({ text }),
    });
    if (!body.isConnected) return;
    body.innerHTML = rendered.html;
    body.classList.add("markdown-body");
    body.querySelectorAll("a").forEach((link) => {
      link.target = "_blank";
      link.rel = "noopener noreferrer";
    });
  } catch (_) {
    // Plain text remains readable if local Markdown rendering is unavailable.
  }
}

export async function refreshDataset() {
  state.dataset = await api("/api/datasets/latest");
  return state.dataset;
}

export function setAgentBusy(busy) {
  state.agentBusy = busy;
  byId("assistant-input").disabled = busy;
  byId("assistant-mode").disabled = busy;
  if (busy) closeAssistantModeMenu();
  const send = byId("assistant-send");
  send.classList.toggle("busy", busy);
  setTranslatedAttribute(send, "aria-label", busy ? "Working…" : "Send");
  setTranslatedAttribute(send, "title", busy ? "Working…" : "Send");
  updateActionStates();
}

export function resizeAssistantInput() {
  const input = byId("assistant-input");
  input.style.height = "auto";
  input.style.height = `${Math.min(input.scrollHeight, 180)}px`;
}

export function newAgentBody(
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
    forwardedProps: {
      mode,
      datasetId,
      projectId: state.projectId,
      formState: agentWorkflowFormState(),
    },
    resume,
  };
}

export function streamMessage(stream, event) {
  const id = event.messageId;
  if (!stream.messages.has(id)) stream.messages.set(id, appendAssistantMessage("assistant", ""));
  return stream.messages.get(id);
}

export function parseToolArguments(tool) {
  try {
    return JSON.parse(tool.args || "{}");
  } catch (_) {
    return {};
  }
}

export function handleReasoningEvent(event, stream) {
  if (event.type === "REASONING_MESSAGE_START") {
    clearAssistantThinking(stream);
    reasoningMessage(stream, event);
  }
  if (event.type === "REASONING_MESSAGE_CONTENT") {
    clearAssistantThinking(stream);
    const reasoning = reasoningMessage(stream, event);
    reasoning.content.textContent += event.delta;
    const messages = reasoning.details.closest(".assistant-messages");
    messages?.scrollTo(0, messages.scrollHeight);
  }
  if (["REASONING_MESSAGE_END", "REASONING_END"].includes(event.type)) {
    finishReasoning(stream);
  }
}

export async function handleTextEvent(event, stream) {
  if (event.type === "TEXT_MESSAGE_START") {
    clearAssistantThinking(stream);
    streamMessage(stream, event);
  }
  if (event.type === "TEXT_MESSAGE_CONTENT") {
    clearAssistantThinking(stream);
    const paragraph = streamMessage(stream, event);
    paragraph.textContent += event.delta;
    const messages = paragraph.closest(".assistant-messages");
    messages?.scrollTo(0, messages.scrollHeight);
  }
  if (event.type === "TEXT_MESSAGE_END") {
    const paragraph = streamMessage(stream, event);
    await renderAssistantMarkdown(paragraph.parentElement, paragraph.textContent);
  }
}

export function handleToolEvent(event, stream) {
  if (event.type === "TOOL_CALL_START") {
    stream.tools.set(event.toolCallId, { name: event.toolCallName, args: "" });
  }
  if (event.type === "TOOL_CALL_ARGS") {
    const tool = stream.tools.get(event.toolCallId);
    if (tool) tool.args += event.delta;
  }
  if (event.type === "TOOL_CALL_RESULT") {
    const tool = stream.tools.get(event.toolCallId);
    if (tool?.name === "fill_workflow_parameters"
        && !stream.appliedTools.has(event.toolCallId)) {
      let result = event.content;
      if (typeof result === "string") {
        try {
          result = JSON.parse(result);
        } catch (_) {
          result = null;
        }
      }
      if (result?.type === "workflow_parameter_patch") {
        applyWorkflowParameterPatch(result);
        stream.appliedTools.add(event.toolCallId);
      }
    }
  }
}

export async function handleAgentEvent(event, stream, mode) {
  handleReasoningEvent(event, stream);
  await handleTextEvent(event, stream);
  handleToolEvent(event, stream);
  if (event.type === "RUN_ERROR") throw new Error(event.message || translate("Assistant run failed"));
  if (event.type === "RUN_FINISHED" && event.outcome?.type === "interrupt") {
    const [interrupt] = event.outcome.interrupts || [];
    if (!interrupt) throw new Error(translate("Assistant paused without an approval request"));
    await presentApproval(interrupt, stream, mode);
  }
}

export async function consumeAgentStream(response, mode, thinking) {
  const stream = {
    buffer: "",
    messages: new Map(),
    reasoning: null,
    reasoningMessageIds: new Set(),
    tools: new Map(),
    appliedTools: new Set(),
    thinking,
  };
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

export async function runAgentRequest(body, mode) {
  if (!state.provider) throw new Error(translate("Configure an AI provider first"));
  if (!state.project) throw new Error(translate("Open a folder first"));
  state.agentAbort?.abort();
  state.agentAbort = new AbortController();
  setAgentBusy(true);
  const thinking = appendAssistantThinking();
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
    await consumeAgentStream(response, mode, thinking);
  } catch (error) {
    if (error?.name !== "AbortError") appendAssistantError(error);
    if (error && typeof error === "object") error.assistantDisplayed = true;
    throw error;
  } finally {
    thinking.remove();
    state.agentAbort = null;
    setAgentBusy(false);
  }
}

export async function sendAssistantMessage() {
  const input = byId("assistant-input");
  const message = input.value.trim();
  if (!message) return;
  const mode = byId("assistant-mode").value;
  appendAssistantMessage("user", message);
  input.value = "";
  resizeAssistantInput();
  await runAgentRequest(newAgentBody(mode, { message }), mode);
}
