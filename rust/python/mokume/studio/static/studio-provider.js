// Model provider settings dialog.
import { api, byId, originHeaders, state, toast } from "./studio-core.js";
import { setTranslatedAttribute, setTranslatedText, translate } from "./studio-i18n.js";
import { closeMenus, updateActionStates } from "./studio-menus.js";

export const BUILTIN_THINKING_LEVELS = new Set(["", "low", "medium", "high", "xhigh", "max"]);

export function renderProviderState() {
  const button = byId("assistant-model");
  const name = byId("assistant-model-name");
  if (state.provider) {
    button.classList.add("configured");
    name.removeAttribute("data-i18n");
    name.removeAttribute("data-i18n-values");
    name.textContent = state.provider.model;
  } else {
    button.classList.remove("configured");
    setTranslatedText(name, "Configure model");
  }
  byId("assistant-panel").classList.toggle("provider-configured", Boolean(state.provider));
  updateActionStates();
}

export function renderProviderKeyNote() {
  let key = "The key is held only in server memory. It is not written to the project or Studio database.";
  if (byId("provider-persist").checked) {
    key = "The API key will be stored in Mokume's mokume-studio-providers.json.";
  } else if (state.provider?.api_key_configured) {
    key = "The configured API key is loaded into this form. Editing it replaces the current key.";
  }
  setTranslatedText(byId("provider-key-note"), key);
}

export function optionalNumber(id) {
  const input = byId(id);
  return input.value === "" ? null : input.valueAsNumber;
}

export function providerThinkingLevel() {
  const selected = byId("provider-thinking-level").value;
  if (selected !== "custom") return selected || null;
  return byId("provider-thinking-custom").value.trim().toLowerCase() || null;
}

export function updateProviderThinkingControl() {
  const input = byId("provider-thinking-custom");
  const custom = byId("provider-thinking-level").value === "custom";
  input.hidden = !custom;
  input.disabled = !custom;
  input.required = custom;
}

export function restoreProviderThinkingLevel(level) {
  const normalized = level || "";
  const builtin = BUILTIN_THINKING_LEVELS.has(normalized);
  byId("provider-thinking-level").value = builtin ? normalized : "custom";
  byId("provider-thinking-custom").value = builtin ? "" : normalized;
  updateProviderThinkingControl();
}

export function providerPayload() {
  const apiKey = byId("provider-api-key").value.trim();
  const baseUrl = byId("provider-base-url").value.trim();
  return {
    provider: byId("provider-kind").value,
    model: byId("provider-model").value.trim(),
    api_key: apiKey || null,
    base_url: baseUrl || null,
    context_tokens: optionalNumber("provider-context-tokens"),
    max_output_tokens: optionalNumber("provider-max-output-tokens"),
    thinking_level: providerThinkingLevel(),
    persist: byId("provider-persist").checked,
  };
}

export function clearProviderTestStatus() {
  const status = byId("provider-test-status");
  status.className = "provider-test-status";
  status.removeAttribute("data-i18n");
  status.removeAttribute("data-i18n-values");
  status.textContent = "";
}

export function setProviderTestError(error) {
  const status = byId("provider-test-status");
  status.className = "provider-test-status error";
  status.removeAttribute("data-i18n");
  status.removeAttribute("data-i18n-values");
  status.textContent = error.message || String(error);
}

export function resetProviderTestFeedback() {
  byId("test-provider").dataset.state = "idle";
  clearProviderTestStatus();
}

export function setProviderKeyVisibility(visible) {
  const input = byId("provider-api-key");
  const button = byId("toggle-provider-key");
  const label = visible ? "Hide API key" : "Show API key";
  input.classList.toggle("secret-masked", !visible);
  button.setAttribute("aria-pressed", String(visible));
  setTranslatedAttribute(button, "aria-label", label);
  setTranslatedAttribute(button, "title", label);
}

export async function openProviderDialog() {
  closeMenus();
  const summary = await api("/api/ai/config");
  state.provider = summary;
  byId("provider-kind").value = summary?.provider || "openai-responses";
  byId("provider-model").value = summary?.model || "";
  byId("provider-base-url").value = summary?.base_url || "";
  byId("provider-api-key").value = summary?.api_key || "";
  byId("provider-context-tokens").value = summary?.context_tokens ?? "";
  byId("provider-max-output-tokens").value = summary?.max_output_tokens ?? "";
  restoreProviderThinkingLevel(summary?.thinking_level);
  byId("provider-persist").checked = Boolean(summary?.persistent);
  byId("provider-advanced").open = false;
  setProviderKeyVisibility(false);
  resetProviderTestFeedback();
  renderProviderKeyNote();
  byId("provider-dialog").showModal();
}

export async function saveProvider(event) {
  event.preventDefault();
  state.provider = await api("/api/ai/config", {
    method: "POST",
    headers: originHeaders(),
    body: JSON.stringify(providerPayload()),
  });
  byId("provider-dialog").close();
  renderProviderState();
  toast(translate("Provider configuration saved"));
}

export async function testProviderConnection() {
  const form = byId("provider-form");
  if (!form.reportValidity()) return;
  const button = byId("test-provider");
  const label = byId("provider-test-label");
  button.disabled = true;
  button.dataset.state = "testing";
  setTranslatedText(label, "Testing connection…");
  clearProviderTestStatus();
  try {
    await api("/api/ai/config/test", {
      method: "POST",
      headers: originHeaders(),
      body: JSON.stringify(providerPayload()),
    });
    button.dataset.state = "success";
  } catch (error) {
    button.dataset.state = "error";
    setProviderTestError(error);
  } finally {
    button.disabled = false;
    setTranslatedText(label, "Test service");
  }
}
