// Menu bar, submenus, assistant mode menu, parameter selectors, and action availability.
import { byId, state } from "./studio-core.js";
import { setTranslatedText } from "./studio-i18n.js";

export function menuItems(menu) {
  return [...menu.querySelectorAll('[role="menuitem"]:not([disabled]), [role="menuitemradio"]:not([disabled])')];
}

export function closeParameterSelectors({ restoreFocus = false } = {}) {
  const openSelector = document.querySelector(".parameter-selector.open");
  const stateKey = openSelector?.dataset.parameterSelector;
  document.querySelectorAll(".parameter-selector.open").forEach((selector) => {
    selector.classList.remove("open");
    selector.querySelector(".parameter-selector-trigger")?.setAttribute("aria-expanded", "false");
    const menu = selector.querySelector(".parameter-selector-menu");
    if (menu) menu.hidden = true;
    selector.querySelectorAll(".parameter-selector-option").forEach((option) => {
      option.tabIndex = -1;
    });
  });
  if (!restoreFocus || !stateKey) return;
  const visibleSelector = [...document.querySelectorAll(".parameter-selector")]
    .find((selector) => selector.dataset.parameterSelector === stateKey
      && !selector.closest(".form-field")?.hidden);
  visibleSelector?.querySelector(".parameter-selector-trigger")?.focus();
}

export function openParameterSelector(selector, focusOption = false) {
  const trigger = selector.querySelector(".parameter-selector-trigger");
  const menu = selector.querySelector(".parameter-selector-menu");
  if (trigger.disabled) return;
  closeMenus();
  selector.classList.add("open");
  trigger.setAttribute("aria-expanded", "true");
  menu.hidden = false;
  const bounds = trigger.getBoundingClientRect();
  menu.style.minWidth = `${Math.max(136, bounds.width)}px`;
  const maximumLeft = Math.max(6, window.innerWidth - menu.offsetWidth - 6);
  menu.style.left = `${Math.min(Math.max(6, bounds.left), maximumLeft)}px`;
  const below = bounds.bottom + 4;
  menu.style.top = `${below + menu.offsetHeight <= window.innerHeight - 6
    ? below
    : Math.max(6, bounds.top - menu.offsetHeight - 4)}px`;
  if (!focusOption) return;
  const options = [...menu.querySelectorAll(".parameter-selector-option")];
  const selected = options.find((option) => option.getAttribute("aria-selected") === "true") || options[0];
  options.forEach((option) => { option.tabIndex = option === selected ? 0 : -1; });
  selected?.focus();
}

export function navigateParameterSelector(event, selector) {
  const options = [...selector.querySelectorAll(".parameter-selector-option")];
  const current = options.indexOf(document.activeElement);
  if (["ArrowDown", "ArrowUp", "Home", "End"].includes(event.key)) {
    event.preventDefault();
    let next = event.key === "Home" ? 0 : options.length - 1;
    if (event.key === "ArrowDown") next = (current + 1) % options.length;
    if (event.key === "ArrowUp") next = (current - 1 + options.length) % options.length;
    options.forEach((option, index) => { option.tabIndex = index === next ? 0 : -1; });
    options[next]?.focus();
  } else if (event.key === "Escape") {
    event.preventDefault();
    closeParameterSelectors({ restoreFocus: true });
  } else if (event.key === "Tab") {
    closeParameterSelectors();
  }
}

export let submenuCloseTimer = null;

export function cancelSubmenuClose() {
  window.clearTimeout(submenuCloseTimer);
  submenuCloseTimer = null;
}

export function closeSubmenus({ restoreFocus = false } = {}) {
  cancelSubmenuClose();
  const activeTrigger = document.querySelector('.submenu-trigger[aria-expanded="true"]');
  document.querySelectorAll(".submenu-popover.open").forEach((menu) => {
    menu.classList.remove("open");
    menuItems(menu).forEach((item) => { item.tabIndex = -1; });
  });
  document.querySelectorAll('.submenu-trigger[aria-expanded="true"]').forEach((trigger) => {
    trigger.setAttribute("aria-expanded", "false");
  });
  if (restoreFocus) activeTrigger?.focus();
}

export function openSubmenu(trigger, { focusFirst = false } = {}) {
  const menu = byId(trigger.dataset.submenu);
  if (!menu) return;
  closeSubmenus();
  menu.classList.add("open");
  trigger.setAttribute("aria-expanded", "true");
  const bounds = trigger.getBoundingClientRect();
  const parentBounds = trigger.closest(".menu-popover")?.getBoundingClientRect() || bounds;
  const gap = 0;
  const opensRight = parentBounds.right + gap + menu.offsetWidth <= window.innerWidth - 6;
  const left = opensRight
    ? parentBounds.right + gap
    : Math.max(6, parentBounds.left - menu.offsetWidth - gap);
  const maximumTop = Math.max(6, window.innerHeight - menu.offsetHeight - 6);
  menu.style.left = `${left}px`;
  menu.style.top = `${Math.min(Math.max(6, bounds.top), maximumTop)}px`;
  if (focusFirst) focusMenuItem(menu, 0);
}

export function scheduleSubmenuClose() {
  cancelSubmenuClose();
  submenuCloseTimer = window.setTimeout(() => closeSubmenus(), 140);
}

export function closeMenus({ restoreFocus = false } = {}) {
  const activeTrigger = document.querySelector(".menu-trigger.active");
  closeParameterSelectors();
  closeAssistantModeMenu();
  closeSubmenus();
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

export function openMenu(trigger, { focusFirst = false } = {}) {
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

export function focusMenuItem(menu, index) {
  const items = menuItems(menu);
  items.forEach((item, itemIndex) => { item.tabIndex = itemIndex === index ? 0 : -1; });
  items[index]?.focus();
}

export function closeAssistantModeMenu({ restoreFocus = false } = {}) {
  const trigger = byId("assistant-mode");
  const menu = byId("assistant-mode-menu");
  menu.classList.remove("open");
  trigger.setAttribute("aria-expanded", "false");
  menu.querySelectorAll('[role="option"]').forEach((option) => { option.tabIndex = -1; });
  if (restoreFocus) trigger.focus();
}

export function openAssistantModeMenu() {
  closeMenus();
  const trigger = byId("assistant-mode");
  const menu = byId("assistant-mode-menu");
  menu.classList.add("open");
  trigger.setAttribute("aria-expanded", "true");
  const options = [...menu.querySelectorAll('[role="option"]')];
  const selected = options.find((option) => option.dataset.assistantMode === trigger.value) || options[0];
  options.forEach((option) => { option.tabIndex = option === selected ? 0 : -1; });
  selected?.focus();
}

export function selectAssistantMode(mode) {
  const trigger = byId("assistant-mode");
  trigger.value = mode;
  setTranslatedText(byId("assistant-mode-value"), mode === "agent" ? "Agent" : "Ask");
  byId("assistant-mode-icon-use").setAttribute("href", `#assistant-icon-${mode}`);
  document.querySelectorAll("[data-assistant-mode]").forEach((option) => {
    option.setAttribute("aria-selected", String(option.dataset.assistantMode === mode));
  });
  closeAssistantModeMenu({ restoreFocus: true });
}

export function navigateAssistantModeMenu(event) {
  const options = [...byId("assistant-mode-menu").querySelectorAll('[role="option"]')];
  const current = options.indexOf(document.activeElement);
  if (["ArrowDown", "ArrowUp", "Home", "End"].includes(event.key)) {
    event.preventDefault();
    let next = event.key === "Home" ? 0 : options.length - 1;
    if (event.key === "ArrowDown") next = (current + 1) % options.length;
    if (event.key === "ArrowUp") next = (current - 1 + options.length) % options.length;
    options.forEach((option, index) => { option.tabIndex = index === next ? 0 : -1; });
    options[next]?.focus();
  }
  if (event.key === "Escape") {
    event.preventDefault();
    closeAssistantModeMenu({ restoreFocus: true });
  }
}

export function showInfo(title, content) {
  byId("info-title").textContent = title;
  byId("info-content").textContent = content;
  byId("info-dialog").showModal();
}

export function setActionDisabled(action, disabled) {
  document.querySelectorAll(`[data-action="${action}"]`).forEach((element) => {
    if (element instanceof HTMLButtonElement) element.disabled = disabled;
  });
}

export function actionEnabled(action) {
  const controls = [...document.querySelectorAll(`[data-action="${action}"]`)]
    .filter((element) => element instanceof HTMLButtonElement);
  return !controls.length || controls.some((control) => !control.disabled);
}

export function updateActionStates() {
  const canConfigure = Boolean(state.project && state.selectedCommand);
  const workspaceLocked = Boolean(state.activeRunId || state.agentBusy || state.pendingApproval);
  setActionDisabled("open-folder", workspaceLocked);
  setActionDisabled("refresh-files", !state.project);
  setActionDisabled("collapse-folders", !state.project);
  setActionDisabled("close-project", !state.project || workspaceLocked);
  setActionDisabled("import-workflow-template", !state.project);
  setActionDisabled("export-workflow-template", !canConfigure);
  setActionDisabled("validate-command", !canConfigure);
  setActionDisabled("queue-command", !canConfigure);
  setActionDisabled("run-command", !canConfigure || Boolean(state.activeRunId));
  setActionDisabled("cancel-run", !state.activeRunId);
  setActionDisabled("exit-studio", Boolean(state.activeRunId));
  setActionDisabled("assistant-history", !state.project || state.agentBusy || Boolean(state.pendingApproval));
  setActionDisabled("new-assistant-chat", state.agentBusy || Boolean(state.pendingApproval));
  byId("assistant-send").disabled = !state.project || !state.provider || state.agentBusy;
  byId("run-indicator").hidden = !state.activeRunId;
}

export function bindMenuEvents() {
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
  document.querySelectorAll(".submenu-trigger").forEach((trigger) => {
    trigger.addEventListener("mouseenter", () => openSubmenu(trigger));
    trigger.addEventListener("mouseleave", scheduleSubmenuClose);
    trigger.addEventListener("click", (event) => {
      event.stopPropagation();
      openSubmenu(trigger, { focusFirst: true });
    });
  });
  document.querySelectorAll(".submenu-popover").forEach((menu) => {
    menuItems(menu).forEach((item) => { item.tabIndex = -1; });
    menu.addEventListener("mouseenter", cancelSubmenuClose);
    menu.addEventListener("mouseleave", scheduleSubmenuClose);
    menu.addEventListener("keydown", (event) => navigateSubmenu(event, menu));
  });
  document.querySelectorAll(".menu-popover [role='menuitem']:not(.submenu-trigger)").forEach((item) => {
    item.addEventListener("mouseenter", closeSubmenus);
  });
  document.querySelectorAll(".menu-popover a").forEach((link) => {
    link.addEventListener("click", () => closeMenus({ restoreFocus: true }));
  });
  document.addEventListener("click", () => closeMenus());
  document.addEventListener("focusin", (event) => {
    const activeTrigger = document.querySelector(".menu-trigger.active");
    const insideMenu = event.target.closest?.(".menu-popover, .submenu-popover, .assistant-mode-control, .parameter-selector");
    if (!insideMenu && event.target !== activeTrigger) closeMenus();
  });
  document.querySelector(".workflow-panel").addEventListener("scroll", () => {
    closeParameterSelectors();
  }, { passive: true });
  window.addEventListener("resize", () => {
    closeMenus();
    document.querySelector(".workspace").classList.remove("sidebar-mobile-visible", "assistant-mobile-visible");
  });
}

export function navigateMenu(event, menu, triggers) {
  const items = menuItems(menu);
  const current = items.indexOf(document.activeElement);
  const currentItem = items[current];
  if (currentItem?.dataset.submenu && ["ArrowRight", "Enter", " "].includes(event.key)) {
    event.preventDefault();
    openSubmenu(currentItem, { focusFirst: true });
    return;
  }
  if (["ArrowDown", "ArrowUp", "Home", "End"].includes(event.key)) {
    event.preventDefault();
    let next = event.key === "Home" ? 0 : items.length - 1;
    if (event.key === "ArrowDown") next = (current + 1) % items.length;
    if (event.key === "ArrowUp") next = (current - 1 + items.length) % items.length;
    closeSubmenus();
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
    closeSubmenus();
    const trigger = document.querySelector(`[data-menu="${menu.id}"]`);
    const index = triggers.indexOf(trigger);
    const delta = event.key === "ArrowRight" ? 1 : -1;
    openMenu(triggers[(index + delta + triggers.length) % triggers.length], { focusFirst: true });
  }
}

export function navigateSubmenu(event, menu) {
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
  if (["ArrowLeft", "Escape"].includes(event.key)) {
    event.preventDefault();
    closeSubmenus({ restoreFocus: true });
  }
}
