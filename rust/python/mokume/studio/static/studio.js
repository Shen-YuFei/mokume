// Mokume Studio entry point: startup, project switching, actions, events, and shortcuts.
import {
  LANGUAGE_STORAGE_KEY,
  api,
  byId,
  originHeaders,
  reportError,
  state,
} from "./studio-core.js";
import { setTranslatedText, translate } from "./studio-i18n.js";
import {
  actionEnabled,
  bindMenuEvents,
  closeAssistantModeMenu,
  closeMenus,
  navigateAssistantModeMenu,
  openAssistantModeMenu,
  selectAssistantMode,
  showInfo,
  updateActionStates,
} from "./studio-menus.js";
import {
  SYSTEM_DARK_QUERY,
  SYSTEM_MEMORY_REFRESH_MS,
  applyAppearance,
  applyPanelSizes,
  bindPanelResizers,
  refreshSystemMemory,
  renderSystemMemory,
  toggleAssistant,
  toggleBottom,
  toggleSidebar,
} from "./studio-layout.js";
import {
  bindFileDrops,
  collapseFileTree,
  loadFolders,
  openFolderDialog,
  refreshFiles,
  selectFolder,
} from "./studio-files.js";
import {
  loadWorkflowTemplateDirectory,
  openWorkflowTemplateDialog,
  refreshCommands,
  submitWorkflowTemplate,
  updateWorkflowTemplateConfirm,
} from "./studio-workflow.js";
import {
  queueCommand,
  renderCommandReview,
  runCommand,
  submitReviewedRun,
  validateCommand,
} from "./studio-review.js";
import {
  addReviewedToQueue,
  cancelRun,
  loadBatchQueue,
  navigateBottomTabs,
  processBatchQueue,
  refreshArtifacts,
  refreshRuns,
  renderBottom,
  showBottomTab,
} from "./studio-runs.js";
import {
  copySelectedRunCommand,
  loadSelectedRunParameters,
  prepareSelectedRerun,
  renderRunComparison,
  renderRunDetails,
  updateSelectedRunComparison,
} from "./studio-run-details.js";
import {
  openProviderDialog,
  renderProviderKeyNote,
  renderProviderState,
  resetProviderTestFeedback,
  saveProvider,
  setProviderKeyVisibility,
  testProviderConnection,
  updateProviderThinkingControl,
} from "./studio-provider.js";
import { decideApproval, renderApproval, restorePendingApproval } from "./studio-approval.js";
import {
  openConversationHistory,
  refreshDataset,
  resetAssistantConversation,
  resizeAssistantInput,
  sendAssistantMessage,
  setAgentBusy,
  startNewAssistantChat,
  threadsForProject,
} from "./studio-assistant.js";

export function translateInterface(language, persist = true) {
  state.language = language === "zh-CN" ? "zh-CN" : "en";
  document.documentElement.lang = state.language;
  document.querySelectorAll("[data-i18n]").forEach((element) => {
    const values = JSON.parse(element.dataset.i18nValues || "{}");
    element.textContent = translate(element.dataset.i18n, values);
  });
  for (const attribute of ["placeholder", "title", "aria-label"]) {
    document.querySelectorAll(`[data-i18n-${attribute}]`).forEach((element) => {
      const values = JSON.parse(element.getAttribute(`data-i18n-${attribute}-values`) || "{}");
      element.setAttribute(attribute, translate(element.getAttribute(`data-i18n-${attribute}`), values));
    });
  }
  document.querySelectorAll("[data-language]").forEach((button) => {
    button.setAttribute("aria-checked", String(button.dataset.language === state.language));
  });
  if (persist) localStorage.setItem(LANGUAGE_STORAGE_KEY, state.language);
  renderProviderState();
  setAgentBusy(state.agentBusy);
  if (state.pendingApproval) renderApproval(state.pendingApproval.record);
  renderSystemMemory();
  if (byId("command-review-dialog")?.open) renderCommandReview();
  if (byId("run-details-dialog")?.open) renderRunDetails();
  if (byId("run-compare-dialog")?.open) renderRunComparison();
  renderBottom();
}

export async function refreshProject() {
  state.project = await api("/api/project");
  const hasProject = Boolean(state.project);
  const projectId = state.project?.id || null;
  if (projectId !== state.projectId) {
    state.agentAbort?.abort();
    if (byId("approval-dialog").open) byId("approval-dialog").close();
    if (byId("conversation-dialog").open) byId("conversation-dialog").close();
    if (byId("workflow-template-dialog").open) byId("workflow-template-dialog").close();
    if (byId("command-review-dialog").open) byId("command-review-dialog").close();
    if (byId("run-details-dialog").open) byId("run-details-dialog").close();
    if (byId("run-compare-dialog").open) byId("run-compare-dialog").close();
    state.projectId = projectId;
    state.dataset = null;
    state.pendingApproval = null;
    state.commandReview = null;
    state.selectedRunDetails = null;
    state.runComparisonDetails = null;
    state.qcDetails = null;
    state.compareRunIds.clear();
    state.batchQueue = loadBatchQueue(projectId);
    state.threads = threadsForProject(state.projectId);
    resetAssistantConversation();
  }
  const chip = byId("project-chip");
  if (hasProject) {
    chip.removeAttribute("data-i18n");
    chip.textContent = state.project.root.split(/[\\/]/).filter(Boolean).at(-1) || state.project.root;
    chip.title = state.project.root;
  } else {
    setTranslatedText(chip, "No project");
    chip.removeAttribute("title");
  }
  byId("welcome").classList.toggle("hidden", hasProject);
  byId("workflow").classList.toggle("hidden", !hasProject);
  if (hasProject) {
    await Promise.all([refreshFiles(), refreshCommands(), refreshDataset()]);
  } else {
    state.selectedCommand = null;
    state.commands = [];
    setTranslatedText(byId("project-files"), "Open a folder to browse input files.");
    byId("command-catalog").replaceChildren();
    setTranslatedText(byId("command-form"), "Select a workflow to configure its parameters.");
  }
  updateActionStates();
}

export function shortcutHelp() {
  return [
    `Ctrl+O  ${translate("Open Folder")}`,
    `Alt+1  ${translate("Toggle Sidebar")}`,
    `Alt+2  ${translate("Toggle Bottom Panel")}`,
    `Alt+3  ${translate("Artifacts")}`,
    `Alt+4  ${translate("Toggle Assistant")}`,
    `Ctrl+Shift+Enter  ${translate("Validate")}`,
    `Ctrl+Enter  ${translate("Run")}`,
    `Ctrl+.  ${translate("Cancel Run")}`,
    `Ctrl+/  ${translate("Keyboard Shortcuts")}`,
    `F11  ${translate("Full Screen")}`,
  ].join("\n");
}

export async function performAction(action) {
  closeMenus({ restoreFocus: true });
  if (!actionEnabled(action)) return;
  const handlers = {
    "open-folder": openFolderDialog,
    "refresh-files": refreshFiles,
    "collapse-folders": async () => collapseFileTree(),
    "assistant-history": openConversationHistory,
    "close-project": async () => {
      await api("/api/projects/close", { method: "POST", headers: originHeaders() });
      await refreshProject();
      state.logs = [];
      await refreshRuns();
    },
    "exit-studio": async () => {
      await api("/api/studio/exit", { method: "POST", headers: originHeaders() });
      document.body.textContent = translate("Mokume Studio is stopping. You can close this tab.");
    },
    "import-workflow-template": async () => openWorkflowTemplateDialog("import"),
    "export-workflow-template": async () => openWorkflowTemplateDialog("export"),
    "validate-command": validateCommand,
    "queue-command": queueCommand,
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
    shortcuts: async () => showInfo(translate("Keyboard Shortcuts"), shortcutHelp()),
    "system-status": async () => showInfo(translate("System Status"), JSON.stringify(await api("/api/system"), null, 2)),
    about: async () => showInfo(translate("About Mokume"), `Mokume Studio ${state.version}\n${translate("A local, Rust-backed proteomics analysis workspace.")}`),
    "new-assistant-chat": async () => startNewAssistantChat(),
    "assistant-settings": openProviderDialog,
  };
  if (handlers[action]) await handlers[action]();
}

export function bindEvents() {
  bindMenuEvents();
  bindPanelResizers();
  SYSTEM_DARK_QUERY.addEventListener("change", () => {
    if (state.appearance === "system") applyAppearance("system", false);
  });
  document.querySelectorAll("[data-appearance]").forEach((button) => {
    button.addEventListener("click", (event) => {
      event.stopPropagation();
      applyAppearance(button.dataset.appearance);
      closeMenus({ restoreFocus: true });
    });
  });
  document.querySelectorAll("[data-language]").forEach((button) => {
    button.addEventListener("click", (event) => {
      event.stopPropagation();
      translateInterface(button.dataset.language);
      closeMenus({ restoreFocus: true });
    });
  });
  byId("assistant-mode").addEventListener("click", (event) => {
    event.stopPropagation();
    if (byId("assistant-mode-menu").classList.contains("open")) closeAssistantModeMenu();
    else openAssistantModeMenu();
  });
  byId("assistant-mode").addEventListener("keydown", (event) => {
    if (!["ArrowDown", "ArrowUp"].includes(event.key)) return;
    event.preventDefault();
    openAssistantModeMenu();
  });
  byId("assistant-mode-menu").addEventListener("keydown", navigateAssistantModeMenu);
  document.querySelectorAll("[data-assistant-mode]").forEach((option) => {
    option.addEventListener("click", (event) => {
      event.stopPropagation();
      selectAssistantMode(option.dataset.assistantMode);
    });
  });
  document.querySelectorAll("[data-action]").forEach((element) => element.addEventListener("click", (event) => {
    event.stopPropagation();
    performAction(element.dataset.action).catch(reportError);
  }));
  document.querySelectorAll("[data-close-dialog]").forEach((button) => button.addEventListener("click", () => button.closest("dialog").close()));
  document.querySelectorAll("[data-bottom-tab]").forEach((button) => button.addEventListener("click", () => showBottomTab(button.dataset.bottomTab)));
  document.querySelector(".bottom-tab-list").addEventListener("keydown", navigateBottomTabs);
  document.querySelectorAll("[data-assistant-prompt]").forEach((button) => button.addEventListener("click", () => {
    const input = byId("assistant-input");
    input.value = translate(button.dataset.assistantPrompt);
    resizeAssistantInput();
    input.focus();
  }));
  bindFileDrops();
  byId("folder-parent").addEventListener("click", () => loadFolders(state.folderParent).catch(reportError));
  byId("select-folder").addEventListener("click", () => selectFolder().catch(reportError));
  byId("workflow-template-parent").addEventListener("click", () => {
    loadWorkflowTemplateDirectory(state.workflowTemplateParent).catch(reportError);
  });
  byId("workflow-template-name").addEventListener("input", () => {
    state.workflowTemplateSelected = null;
    byId("workflow-template-list").querySelectorAll(".workflow-template-row:not(.directory)")
      .forEach((row) => row.setAttribute("aria-pressed", "false"));
    updateWorkflowTemplateConfirm();
  });
  byId("workflow-template-form").addEventListener("submit", (event) => {
    submitWorkflowTemplate(event).catch(reportError);
  });
  byId("provider-form").addEventListener("submit", (event) => saveProvider(event).catch(reportError));
  byId("provider-form").addEventListener("input", resetProviderTestFeedback);
  byId("provider-form").addEventListener("change", resetProviderTestFeedback);
  byId("provider-persist").addEventListener("change", renderProviderKeyNote);
  byId("provider-thinking-level").addEventListener("change", updateProviderThinkingControl);
  byId("provider-thinking-custom").addEventListener("input", (event) => {
    event.target.value = event.target.value.toLowerCase();
  });
  byId("test-provider").addEventListener("click", () => testProviderConnection());
  byId("toggle-provider-key").addEventListener("click", () => {
    setProviderKeyVisibility(byId("provider-api-key").classList.contains("secret-masked"));
  });
  byId("assistant-send").addEventListener("click", () => sendAssistantMessage().catch(reportError));
  byId("assistant-input").addEventListener("input", resizeAssistantInput);
  byId("assistant-input").addEventListener("keydown", (event) => {
    if (event.key === "Enter" && !event.shiftKey) {
      event.preventDefault();
      sendAssistantMessage().catch(reportError);
    }
  });
  byId("approval-approve").addEventListener("click", () => decideApproval(true).catch(reportError));
  byId("approval-reject").addEventListener("click", () => decideApproval(false).catch(reportError));
  byId("review-run").addEventListener("click", () => submitReviewedRun().catch(reportError));
  byId("review-queue").addEventListener("click", () => addReviewedToQueue().catch(reportError));
  byId("run-load-parameters").addEventListener("click", loadSelectedRunParameters);
  byId("run-copy-command").addEventListener("click", () => copySelectedRunCommand().catch(reportError));
  byId("run-rerun").addEventListener("click", () => prepareSelectedRerun().catch(reportError));
  byId("run-compare-choice").addEventListener("change", updateSelectedRunComparison);
  document.addEventListener("keydown", keyboardShortcut);
}

export function keyboardShortcut(event) {
  if (event.defaultPrevented) return;
  const control = event.ctrlKey || event.metaKey;
  // Ctrl+Enter still validates or runs while a parameter field has focus.
  const formSubmit = control && event.key === "Enter" && Boolean(event.target.closest?.("#command-form"));
  if (!formSubmit && event.target.closest?.("input, textarea, select, dialog")) return;
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

export async function boot() {
  applyAppearance(state.appearance, false);
  translateInterface(state.language, false);
  applyPanelSizes();
  const session = await api("/api/session");
  state.csrf = session.csrf_token;
  state.version = session.version;
  state.provider = session.ai_provider;
  bindEvents();
  await Promise.all([refreshProject(), refreshRuns(), refreshArtifacts(), refreshSystemMemory()]);
  await processBatchQueue();
  renderProviderState();
  await restorePendingApproval();
  renderBottom();
  window.setInterval(() => void refreshSystemMemory(), SYSTEM_MEMORY_REFRESH_MS);
}

boot().catch(reportError);
