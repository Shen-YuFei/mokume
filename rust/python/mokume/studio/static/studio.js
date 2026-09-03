const ACTIVE_RUN_STATUSES = new Set(["queued", "starting", "running", "cancelling"]);
const WORKFLOW_TEMPLATE_SCHEMA_VERSION = 1;
const LANGUAGE_STORAGE_KEY = "mokume:language";
const APPEARANCE_STORAGE_KEY = "mokume:appearance";
const SYSTEM_DARK_QUERY = window.matchMedia("(prefers-color-scheme: dark)");
const DEFAULT_PANEL_SIZES = Object.freeze({ sidebar: 235, assistant: 350, bottom: 188 });
const PANEL_SIZE_LIMITS = Object.freeze({
  sidebar: Object.freeze({ min: 180 }),
  assistant: Object.freeze({ min: 280 }),
  bottom: Object.freeze({ min: 120, max: 560 }),
});
const SIDE_PANEL_MAX_VIEWPORT_RATIO = 0.5;
const SIDE_PANEL_COLLAPSE_RATIO = 0.5;
const SYSTEM_MEMORY_REFRESH_MS = 5000;
const BUILTIN_THINKING_LEVELS = new Set(["", "low", "medium", "high", "xhigh", "max"]);
const PROTEOMICS_FILE_RULES = Object.freeze([
  Object.freeze({ pattern: /(^|[._-])sdrf([._-]|$)/, icon: "sdrf", tone: "sdrf" }),
  Object.freeze({ pattern: /(^|[._-])msstats([._-]|$)/, icon: "msstats", tone: "stats" }),
]);
const FILE_ICON_RULES = Object.freeze([
  Object.freeze({ suffixes: [".parquet", ".arrow", ".feather"], icon: "vscode-parquet", tone: "parquet", filled: true }),
  Object.freeze({ suffixes: [".mzml", ".mzxml", ".mgf", ".raw", ".wiff", ".mzml.gz", ".mzxml.gz", ".mgf.gz"], icon: "mass-spectrum", tone: "spectrum" }),
  Object.freeze({ suffixes: [".featurexml", ".consensusxml", ".osw"], icon: "feature-map", tone: "features" }),
  Object.freeze({ suffixes: [".mzid", ".mzidentml", ".idxml", ".pepxml"], icon: "identification", tone: "identification" }),
  Object.freeze({ suffixes: [".fasta", ".fa", ".faa", ".fna", ".fastq", ".fq", ".fasta.gz", ".fastq.gz"], icon: "dna", tone: "sequence" }),
  Object.freeze({ suffixes: [".nf"], icon: "workflow", tone: "workflow" }),
  Object.freeze({ suffixes: [".matrix"], icon: "table-properties", tone: "matrix" }),
  Object.freeze({ suffixes: [".csv", ".tsv", ".tab", ".xls", ".xlsx", ".mztab"], icon: "file-spreadsheet", tone: "table" }),
  Object.freeze({ suffixes: [".json", ".jsonl", ".geojson"], icon: "braces", tone: "structured" }),
  Object.freeze({ suffixes: [".yaml", ".yml", ".toml", ".ini", ".cfg", ".conf", ".config", ".env", ".quantms", ".diann"], icon: "vscode-config", tone: "config", filled: true }),
  Object.freeze({ suffixes: [".sh", ".bash", ".zsh", ".fish", ".ps1", ".bat", ".cmd"], icon: "vscode-shell", tone: "script", filled: true }),
  Object.freeze({ suffixes: [".log"], icon: "vscode-log", tone: "log", filled: true }),
  Object.freeze({ suffixes: [".py", ".r", ".rs", ".js", ".jsx", ".ts", ".tsx", ".sql"], icon: "file-code", tone: "code" }),
  Object.freeze({ suffixes: [".png", ".jpg", ".jpeg", ".gif", ".webp", ".svg", ".tif", ".tiff"], icon: "image", tone: "image" }),
  Object.freeze({ suffixes: [".zip", ".tar", ".gz", ".bz2", ".xz", ".7z"], icon: "file-archive", tone: "archive" }),
  Object.freeze({ suffixes: [".txt", ".md", ".rst", ".pdf"], icon: "file-text", tone: "document" }),
]);
const MIN_WORKFLOW_WIDTH = 420;
const CHINESE_TRANSLATIONS = Object.freeze({
  File: "文件",
  Analysis: "分析",
  View: "视图",
  Help: "帮助",
  "Application menu": "应用菜单",
  "No project": "未打开项目",
  Idle: "空闲",
  "Open Folder": "打开文件夹",
  "Open folder": "打开文件夹",
  "Refresh files": "刷新文件",
  "Collapse all folders": "收起所有文件夹",
  "Conversation history": "历史对话",
  "Close conversation history": "关闭历史对话",
  "No conversations in this workspace.": "当前工作区暂无历史对话。",
  "Workspace: {path}": "工作区：{path}",
  "Rename conversation: {title}": "重命名对话：{title}",
  "Rename conversation": "重命名对话",
  "Conversation renamed": "对话已重命名",
  "Delete conversation: {title}": "删除对话：{title}",
  "Delete this conversation? This cannot be undone.": "确定删除这条对话吗？此操作无法撤销。",
  "Conversation deleted": "对话已删除",
  "New chat": "新建对话",
  "Collapse assistant": "收起助手",
  "Expand assistant": "展开助手",
  "Resize data panel": "调整数据面板大小",
  "Resize assistant panel": "调整助手面板大小",
  "Resize bottom panel": "调整底部面板大小",
  "Collapse sidebar": "收起侧边栏",
  "Expand sidebar": "展开侧边栏",
  "Collapse bottom panel": "收起底部面板",
  "Close Project": "关闭项目",
  "Exit Mokume Studio": "退出 Mokume Studio",
  "Validate Parameters": "验证参数",
  "Run Analysis": "运行分析",
  "Cancel Run": "取消运行",
  "Run History": "运行历史",
  "Toggle Sidebar": "切换侧边栏",
  "Toggle Assistant": "切换助手",
  "Toggle Bottom Panel": "切换底部面板",
  Artifacts: "产物",
  "Full Screen": "全屏",
  Appearance: "外观",
  System: "跟随系统",
  Light: "浅色",
  Dark: "深色",
  Language: "语言",
  Documentation: "文档",
  "Keyboard Shortcuts": "快捷键",
  "System Status": "系统状态",
  "System memory": "系统内存",
  "System memory: {percent}% ({used} / {total})": "系统内存：{percent}%（{used} / {total}）",
  "About Mokume": "关于 Mokume",
  DATA: "数据",
  "Open a folder to browse input files.": "打开文件夹以浏览输入文件。",
  "Start a Mokume analysis": "开始 Mokume 分析",
  "Open a local folder, configure a reproducible workflow, and review its results.": "打开本地文件夹、配置可复现的工作流并查看结果。",
  "Files stay on this computer. Inputs are opened read-only.": "文件保留在此计算机上，输入文件以只读方式打开。",
  WORKFLOW: "工作流",
  "Configure analysis": "配置分析",
  "Import template": "导入模板",
  "Export template": "导出模板",
  "Import workflow template": "导入工作流模板",
  "Export workflow template": "导出工作流模板",
  "Close workflow template": "关闭工作流模板",
  "Workspace folder": "工作区文件夹",
  "Template file": "模板文件",
  "Select a JSON template": "选择一个 JSON 模板",
  "Enter a JSON file name": "输入 JSON 文件名",
  "No JSON templates in this folder.": "此文件夹中没有 JSON 模板。",
  "Import from workspace": "从工作区导入",
  "Save to workspace": "保存到工作区",
  "Template file name must end with .json": "模板文件名必须以 .json 结尾",
  "Replace the existing workflow template?": "替换现有工作流模板吗？",
  "Workflow template imported": "工作流模板已导入",
  "Workflow template exported": "工作流模板已导出",
  "Invalid workflow template": "工作流模板无效",
  "Unsupported workflow template version": "不支持该工作流模板版本",
  "Template workflow is unavailable: {workflow}": "模板中的工作流不可用：{workflow}",
  "Template parameter is unavailable: --{parameter}": "模板中的参数不可用：--{parameter}",
  "Invalid template value for --{parameter}": "参数 --{parameter} 的模板值无效",
  "The Agent response no longer matches the current workflow": "Agent 返回的参数已不属于当前工作流",
  "Agent updated {count} parameters": "Agent 已填写 {count} 个参数",
  Validate: "验证",
  Run: "运行",
  "Select a workflow to configure its parameters.": "选择工作流以配置参数。",
  "Input & output": "输入与输出",
  Quantification: "定量",
  "Normalization & correction": "归一化与校正",
  "Imputation & differential expression": "缺失值填补与差异表达",
  Runtime: "运行资源",
  Filtering: "过滤",
  "Normalization & aggregation": "归一化与聚合",
  "Internal reference scaling": "内部参照缩放（IRS）",
  "Absolute abundance": "绝对丰度",
  "QC & runtime": "质控与运行资源",
  "Column mapping & export": "列映射与导出",
  "Batch correction": "批次校正",
  "Tissue atlas": "组织图谱",
  Visualization: "可视化",
  Reports: "报告",
  Embedding: "降维",
  "Plots & contrasts": "图表与对比",
  Contrasts: "对比",
  "Thresholds & highlights": "阈值与高亮",
  Parameters: "参数",
  Required: "必填",
  "Input file type": "输入文件类型",
  "IRS reference type": "IRS 参照选择方式",
  "Feature-level QPX input for protein quantification; also supplies protein-group assignments for spectral-count.": "用于蛋白定量的 feature-level QPX 输入；同时为 spectral-count 提供 protein-group 归属。",
  "Feature-level MSstats input for protein quantification; requires SDRF metadata.": "用于蛋白定量的 feature-level MSstats 输入；需要 SDRF 元数据。",
  ASSISTANT: "助手",
  Optional: "可选",
  Mode: "模式",
  Ask: "问答",
  Agent: "Agent",
  Provider: "模型服务",
  "Configure model provider": "配置模型服务",
  "Assistant conversation": "助手对话",
  "Configure model": "配置模型",
  "Message Mokume Assistant": "向 Mokume 助手发送消息",
  "Ask about this workspace or ask Agent to help with an analysis": "询问当前工作区，或让 Agent 协助分析",
  Send: "发送",
  "Working…": "处理中…",
  Thinking: "正在思考",
  "Thinking…": "正在思考…",
  "Thought process": "思考中",
  "Request failed": "请求失败",
  Runs: "运行记录",
  Logs: "日志",
  "No runs yet.": "暂无运行记录。",
  PROJECT: "项目",
  "Parent folder": "上级文件夹",
  "Current folder": "当前文件夹",
  "No readable subfolders.": "没有可读取的子文件夹。",
  Cancel: "取消",
  "Open this folder": "打开此文件夹",
  "This folder is empty.": "此文件夹为空。",
  Done: "完成",
  "Configure Provider": "配置模型服务",
  "Close provider settings": "关闭模型服务设置",
  "API format": "API 格式",
  "Model ID": "模型 ID",
  "e.g. gpt-5, claude-sonnet-4-5, k3-256k": "例如 gpt-5、claude-sonnet-4-5、k3-256k",
  "The selected model must support tool calling.": "所选模型必须支持 tool calling。",
  "Base URL": "基础 URL",
  "Leave empty to use the SDK default endpoint.": "留空时使用 SDK 默认端点。",
  "API key": "API 密钥",
  "Show API key": "显示 API 密钥",
  "Hide API key": "隐藏 API 密钥",
  "Optional when supplied by the server environment": "服务器环境已提供时可留空",
  "The key is held only in server memory. It is not written to the project or Studio database.": "密钥仅保留在服务器内存中，不会写入项目或 Studio 数据库。",
  "The configured API key is loaded into this form. Editing it replaces the current key.": "已配置的 API 密钥会载入此表单；修改后将替换当前密钥。",
  "The API key will be stored in Mokume's mokume-studio-providers.json.": "API 密钥将保存到 Mokume 的 mokume-studio-providers.json。",
  "Persist Studio configuration": "持久化保存 Studio 配置",
  Advanced: "高级设置",
  "Context token limit": "上下文 token 上限",
  "Max output tokens": "最大输出 token",
  "Thinking level": "思考等级",
  "Provider default": "模型服务默认值",
  "Maximum input context per model request.": "单次模型请求允许的最大输入上下文。",
  "Maximum tokens generated in one model response.": "单次模型响应允许生成的最大 token 数。",
  "Supported models use the closest available level; other models may ignore this setting.": "支持思考配置的模型会使用最接近的可用等级；其他模型可能忽略此设置。",
  "Test service": "点击测试服务",
  "Testing connection…": "正在测试连接…",
  Save: "保存",
  "FINAL REVIEW": "最终审核",
  "Approve Analysis": "审批分析",
  "Close analysis approval": "关闭分析审批",
  "Review the canonical parameters below. Computation cannot start until this exact plan is approved.": "请审核以下标准参数。只有批准此确切方案后才能开始计算。",
  "Canonical parameters": "标准参数",
  Pending: "待处理",
  "No analysis plan is awaiting approval.": "当前没有等待审批的分析方案。",
  Reject: "拒绝",
  "Approve and Run": "批准并运行",
  "Resume Approved Run": "继续已批准的运行",
  "Resume Rejection": "继续拒绝",
  "Provider configuration saved": "模型服务配置已保存",
  "The analysis approval did not contain a valid server plan": "分析审批中没有有效的服务器方案",
  "Approval hash mismatch": "审批哈希不匹配",
  "No analysis plan is awaiting approval": "当前没有等待审批的分析方案",
  "This analysis plan already has the opposite decision": "此分析方案已作出相反决定",
  "Assistant run failed": "助手运行失败",
  "Assistant paused without an approval request": "助手已暂停，但没有提供审批请求",
  "Configure an AI provider first": "请先配置 AI 模型服务",
  "Open a folder first": "请先打开文件夹",
  "Opened {path}": "已打开 {path}",
  "Add another": "继续添加",
  Remove: "移除",
  "Default ({value})": "默认值（{value}）",
  "Select…": "选择…",
  "Default: {value}": "默认值：{value}",
  "Optional switch": "可选开关",
  "{hint}; may be repeated": "{hint}；可重复输入",
  "{option} requires {count} values": "{option} 需要 {count} 个值",
  "{option} is required": "{option} 为必填项",
  "Select a workflow first": "请先选择工作流",
  "Validated: {command}": "已验证：{command}",
  "Parameters are valid": "参数有效",
  "Run {id} queued": "运行 {id} 已加入队列",
  "No run is active": "当前没有正在运行的任务",
  "Cancelled run {id}": "已取消运行 {id}",
  "No log events yet.": "暂无日志事件。",
  "No artifacts yet.": "暂无产物。",
  "Mokume Studio is stopping. You can close this tab.": "Mokume Studio 正在停止，可以关闭此标签页。",
  "A local, Rust-backed proteomics analysis workspace.": "基于 Rust 的本地蛋白质组分析工作区。",
  queued: "等待中",
  starting: "启动中",
  running: "运行中",
  cancelling: "正在取消",
  succeeded: "成功",
  failed: "失败",
  cancelled: "已取消",
  interrupted: "已中断",
  approved: "已批准",
  rejected: "已拒绝",
  consumed: "已使用",
  expired: "已过期",
});

function savedLanguage() {
  return localStorage.getItem(LANGUAGE_STORAGE_KEY) === "zh-CN" ? "zh-CN" : "en";
}

function savedAppearance() {
  const appearance = localStorage.getItem(APPEARANCE_STORAGE_KEY);
  return ["system", "light", "dark"].includes(appearance) ? appearance : "system";
}

const state = {
  csrf: null,
  version: null,
  project: null,
  folderPath: null,
  folderParent: null,
  workflowTemplateMode: null,
  workflowTemplateDirectory: ".",
  workflowTemplateParent: null,
  workflowTemplateEntries: [],
  workflowTemplateSelected: null,
  commands: [],
  selectedCommand: null,
  runs: [],
  activeRunId: null,
  eventSource: null,
  eventRunId: null,
  logs: [],
  artifacts: [],
  bottomTab: "runs",
  provider: null,
  dataset: null,
  projectId: null,
  agentBusy: false,
  agentAbort: null,
  pendingApproval: null,
  systemMemory: null,
  threads: { ask: crypto.randomUUID(), agent: crypto.randomUUID() },
  language: savedLanguage(),
  appearance: savedAppearance(),
  panelSizes: { ...DEFAULT_PANEL_SIZES },
};

const byId = (id) => document.getElementById(id);
function translate(text, values = {}) {
  let result = state.language === "zh-CN" ? CHINESE_TRANSLATIONS[text] || text : text;
  Object.entries(values).forEach(([name, value]) => {
    result = result.replaceAll(`{${name}}`, String(value));
  });
  return result;
}

function setTranslatedText(element, text, values = {}) {
  element.dataset.i18n = text;
  element.dataset.i18nValues = JSON.stringify(values);
  element.textContent = translate(text, values);
}

function setTranslatedAttribute(element, attribute, text, values = {}) {
  element.setAttribute(`data-i18n-${attribute}`, text);
  element.setAttribute(`data-i18n-${attribute}-values`, JSON.stringify(values));
  element.setAttribute(attribute, translate(text, values));
}

function applyAppearance(appearance, persist = true) {
  state.appearance = ["light", "dark"].includes(appearance) ? appearance : "system";
  const theme = state.appearance === "system"
    ? (SYSTEM_DARK_QUERY.matches ? "dark" : "light")
    : state.appearance;
  document.documentElement.dataset.theme = theme;
  document.querySelectorAll("[data-appearance]").forEach((button) => {
    button.setAttribute("aria-checked", String(button.dataset.appearance === state.appearance));
  });
  if (persist) localStorage.setItem(APPEARANCE_STORAGE_KEY, state.appearance);
}

function translateInterface(language, persist = true) {
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
  renderBottom();
}

const originHeaders = () => ({
  "Content-Type": "application/json",
  "X-CSRF-Token": state.csrf,
});

const agentHeaders = () => ({
  ...originHeaders(),
  Accept: "text/event-stream",
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
  return [...menu.querySelectorAll('[role="menuitem"]:not([disabled]), [role="menuitemradio"]:not([disabled])')];
}

function closeParameterSelectors({ restoreFocus = false } = {}) {
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

function openParameterSelector(selector, focusOption = false) {
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

function navigateParameterSelector(event, selector) {
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

let submenuCloseTimer = null;

function cancelSubmenuClose() {
  window.clearTimeout(submenuCloseTimer);
  submenuCloseTimer = null;
}

function closeSubmenus({ restoreFocus = false } = {}) {
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

function openSubmenu(trigger, { focusFirst = false } = {}) {
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

function scheduleSubmenuClose() {
  cancelSubmenuClose();
  submenuCloseTimer = window.setTimeout(() => closeSubmenus(), 140);
}

function closeMenus({ restoreFocus = false } = {}) {
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
  if (focusFirst) focusMenuItem(menu, 0);
}

function focusMenuItem(menu, index) {
  const items = menuItems(menu);
  items.forEach((item, itemIndex) => { item.tabIndex = itemIndex === index ? 0 : -1; });
  items[index]?.focus();
}

function closeAssistantModeMenu({ restoreFocus = false } = {}) {
  const trigger = byId("assistant-mode");
  const menu = byId("assistant-mode-menu");
  menu.classList.remove("open");
  trigger.setAttribute("aria-expanded", "false");
  menu.querySelectorAll('[role="option"]').forEach((option) => { option.tabIndex = -1; });
  if (restoreFocus) trigger.focus();
}

function openAssistantModeMenu() {
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

function selectAssistantMode(mode) {
  const trigger = byId("assistant-mode");
  trigger.value = mode;
  setTranslatedText(byId("assistant-mode-value"), mode === "agent" ? "Agent" : "Ask");
  byId("assistant-mode-icon-use").setAttribute("href", `#assistant-icon-${mode}`);
  document.querySelectorAll("[data-assistant-mode]").forEach((option) => {
    option.setAttribute("aria-selected", String(option.dataset.assistantMode === mode));
  });
  closeAssistantModeMenu({ restoreFocus: true });
}

function navigateAssistantModeMenu(event) {
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

function actionEnabled(action) {
  const controls = [...document.querySelectorAll(`[data-action="${action}"]`)]
    .filter((element) => element instanceof HTMLButtonElement);
  return !controls.length || controls.some((control) => !control.disabled);
}

function updateActionStates() {
  const canConfigure = Boolean(state.project && state.selectedCommand);
  const workspaceLocked = Boolean(state.activeRunId || state.agentBusy || state.pendingApproval);
  setActionDisabled("open-folder", workspaceLocked);
  setActionDisabled("refresh-files", !state.project);
  setActionDisabled("collapse-folders", !state.project);
  setActionDisabled("close-project", !state.project || workspaceLocked);
  setActionDisabled("import-workflow-template", !state.project);
  setActionDisabled("export-workflow-template", !canConfigure);
  setActionDisabled("validate-command", !canConfigure);
  setActionDisabled("run-command", !canConfigure || Boolean(state.activeRunId));
  setActionDisabled("cancel-run", !state.activeRunId);
  setActionDisabled("exit-studio", Boolean(state.activeRunId));
  setActionDisabled("assistant-history", !state.project || state.agentBusy || Boolean(state.pendingApproval));
  setActionDisabled("new-assistant-chat", state.agentBusy || Boolean(state.pendingApproval));
  byId("assistant-send").disabled = !state.project || !state.provider || state.agentBusy;
}

function threadsForProject(projectId) {
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

function resetAssistantConversation() {
  byId("assistant-messages").querySelectorAll(".assistant-message").forEach((message) => message.remove());
}

function startNewAssistantChat() {
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

function conversationTimestamp(value) {
  const locale = state.language === "zh-CN" ? "zh-CN" : "en";
  return new Date(value).toLocaleString(locale, { dateStyle: "medium", timeStyle: "short" });
}

function renderConversationHistory(threads, workspace) {
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
    renameButton.addEventListener("click", () => renameStoredConversation(thread).catch(reportError));

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
    deleteButton.addEventListener("click", () => deleteStoredConversation(thread).catch(reportError));
    row.append(openButton, renameButton, deleteButton);
    list.append(row);
  });
}

async function refreshConversationHistory() {
  const payload = await api("/api/agent/threads");
  if (payload.project_id !== state.projectId) {
    throw new Error("Conversation history belongs to a different workspace");
  }
  renderConversationHistory(payload.threads, payload.workspace);
}

async function openConversationHistory() {
  if (!state.projectId) throw new Error(translate("Open a folder first"));
  await refreshConversationHistory();
  byId("conversation-dialog").showModal();
}

async function openStoredConversation(summary) {
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

async function renameStoredConversation(summary) {
  const title = window.prompt(translate("Rename conversation"), summary.title);
  if (title === null || title.trim() === summary.title) return;
  await api(`/api/agent/threads/${encodeURIComponent(summary.id)}?mode=${encodeURIComponent(summary.mode)}`, {
    method: "PATCH",
    headers: originHeaders(),
    body: JSON.stringify({ title }),
  });
  await refreshConversationHistory();
  toast(translate("Conversation renamed"));
}

async function deleteStoredConversation(summary) {
  if (!window.confirm(translate("Delete this conversation? This cannot be undone."))) return;
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

function appendAssistantMessage(role, text) {
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

function appendAssistantThinking() {
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

function clearAssistantThinking(stream) {
  stream.thinking?.remove();
  stream.thinking = null;
}

function appendAssistantReasoning(text = "", active = false) {
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

function reasoningMessage(stream, event) {
  if (!stream.reasoning.has(event.messageId)) {
    stream.reasoning.set(event.messageId, appendAssistantReasoning("", true));
  }
  return stream.reasoning.get(event.messageId);
}

function finishReasoning(stream, event) {
  const reasoning = stream.reasoning.get(event.messageId);
  if (!reasoning) return;
  reasoning.details.classList.remove("active");
  setTranslatedText(reasoning.label, "Thought process");
}

function appendAssistantError(error) {
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

function formatMemorySize(bytes) {
  const gibibytes = bytes / (1024 ** 3);
  return `${gibibytes.toFixed(gibibytes >= 100 ? 0 : 1)} GiB`;
}

function renderSystemMemory() {
  const indicator = byId("system-memory");
  const memory = state.systemMemory;
  if (!memory) {
    indicator.style.setProperty("--memory-percent", "0");
    indicator.setAttribute("aria-valuenow", "0");
    setTranslatedAttribute(indicator, "aria-label", "System memory");
    indicator.title = translate("System memory");
    indicator.classList.add("unavailable");
    return;
  }
  const percent = Math.min(100, Math.max(0, Number(memory.percent) || 0));
  const label = translate("System memory: {percent}% ({used} / {total})", {
    percent: percent.toFixed(1),
    used: formatMemorySize(memory.used_bytes),
    total: formatMemorySize(memory.total_bytes),
  });
  indicator.style.setProperty("--memory-percent", String(percent));
  indicator.setAttribute("aria-valuenow", percent.toFixed(1));
  indicator.setAttribute("aria-label", label);
  indicator.title = label;
  indicator.classList.remove("unavailable");
}

async function refreshSystemMemory() {
  try {
    state.systemMemory = (await api("/api/system")).memory;
  } catch (_) {
    state.systemMemory = null;
  }
  renderSystemMemory();
}

async function renderAssistantMarkdown(body, text) {
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

function renderProviderState() {
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
  updateActionStates();
}

function renderProviderKeyNote() {
  let key = "The key is held only in server memory. It is not written to the project or Studio database.";
  if (byId("provider-persist").checked) {
    key = "The API key will be stored in Mokume's mokume-studio-providers.json.";
  } else if (state.provider?.api_key_configured) {
    key = "The configured API key is loaded into this form. Editing it replaces the current key.";
  }
  setTranslatedText(byId("provider-key-note"), key);
}

function optionalNumber(id) {
  const input = byId(id);
  return input.value === "" ? null : input.valueAsNumber;
}

function providerThinkingLevel() {
  const selected = byId("provider-thinking-level").value;
  if (selected !== "custom") return selected || null;
  return byId("provider-thinking-custom").value.trim().toLowerCase() || null;
}

function updateProviderThinkingControl() {
  const input = byId("provider-thinking-custom");
  const custom = byId("provider-thinking-level").value === "custom";
  input.hidden = !custom;
  input.disabled = !custom;
  input.required = custom;
}

function restoreProviderThinkingLevel(level) {
  const normalized = level || "";
  const builtin = BUILTIN_THINKING_LEVELS.has(normalized);
  byId("provider-thinking-level").value = builtin ? normalized : "custom";
  byId("provider-thinking-custom").value = builtin ? "" : normalized;
  updateProviderThinkingControl();
}

function providerPayload() {
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

function clearProviderTestStatus() {
  const status = byId("provider-test-status");
  status.className = "provider-test-status";
  status.removeAttribute("data-i18n");
  status.removeAttribute("data-i18n-values");
  status.textContent = "";
}

function setProviderTestError(error) {
  const status = byId("provider-test-status");
  status.className = "provider-test-status error";
  status.removeAttribute("data-i18n");
  status.removeAttribute("data-i18n-values");
  status.textContent = error.message || String(error);
}

function resetProviderTestFeedback() {
  byId("test-provider").dataset.state = "idle";
  clearProviderTestStatus();
}

function setProviderKeyVisibility(visible) {
  const input = byId("provider-api-key");
  const button = byId("toggle-provider-key");
  const label = visible ? "Hide API key" : "Show API key";
  input.classList.toggle("secret-masked", !visible);
  button.setAttribute("aria-pressed", String(visible));
  setTranslatedAttribute(button, "aria-label", label);
  setTranslatedAttribute(button, "title", label);
}

async function openProviderDialog() {
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

async function saveProvider(event) {
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

async function testProviderConnection() {
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

async function refreshDataset() {
  state.dataset = await api("/api/datasets/latest");
  return state.dataset;
}

function setAgentBusy(busy) {
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

function resizeAssistantInput() {
  const input = byId("assistant-input");
  input.style.height = "auto";
  input.style.height = `${Math.min(input.scrollHeight, 180)}px`;
}

function newAgentBody(
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

function agentWorkflowFormState() {
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

function streamMessage(stream, event) {
  const id = event.messageId;
  if (!stream.messages.has(id)) stream.messages.set(id, appendAssistantMessage("assistant", ""));
  return stream.messages.get(id);
}

function parseToolArguments(tool) {
  try {
    return JSON.parse(tool.args || "{}");
  } catch (_) {
    return {};
  }
}

async function presentApproval(interrupt, stream, mode) {
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

function pendingApprovalKey() {
  return state.projectId ? `mokume:approval:${state.projectId}` : null;
}

function persistPendingApproval() {
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

function clearPendingApproval() {
  const key = pendingApprovalKey();
  if (key) sessionStorage.removeItem(key);
  state.pendingApproval = null;
  updateActionStates();
}

async function restorePendingApproval() {
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

function renderApproval(record) {
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

async function decideApproval(approved) {
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

function handleReasoningEvent(event, stream) {
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
    finishReasoning(stream, event);
  }
}

async function handleTextEvent(event, stream) {
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

function handleToolEvent(event, stream) {
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

async function handleAgentEvent(event, stream, mode) {
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

async function consumeAgentStream(response, mode, thinking) {
  const stream = {
    buffer: "",
    messages: new Map(),
    reasoning: new Map(),
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

async function runAgentRequest(body, mode) {
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

async function sendAssistantMessage() {
  const input = byId("assistant-input");
  const message = input.value.trim();
  if (!message) return;
  const mode = byId("assistant-mode").value;
  appendAssistantMessage("user", message);
  input.value = "";
  resizeAssistantInput();
  await runAgentRequest(newAgentBody(mode, { message }), mode);
}

async function refreshProject() {
  state.project = await api("/api/project");
  const hasProject = Boolean(state.project);
  const projectId = state.project?.id || null;
  if (projectId !== state.projectId) {
    state.agentAbort?.abort();
    if (byId("approval-dialog").open) byId("approval-dialog").close();
    if (byId("conversation-dialog").open) byId("conversation-dialog").close();
    if (byId("workflow-template-dialog").open) byId("workflow-template-dialog").close();
    state.projectId = projectId;
    state.dataset = null;
    state.pendingApproval = null;
    state.threads = threadsForProject(state.projectId);
    resetAssistantConversation();
  }
  if (hasProject) {
    byId("project-chip").removeAttribute("data-i18n");
    byId("project-chip").textContent = state.project.root;
  } else {
    setTranslatedText(byId("project-chip"), "No project");
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
    setTranslatedText(empty, "No readable subfolders.");
    list.append(empty);
  }
  payload.directories.forEach((directory) => {
    const row = document.createElement("button");
    row.className = "folder-row";
    row.textContent = directory.name;
    row.addEventListener("click", () => loadFolders(directory.path).catch(reportError));
    list.append(row);
  });
}

async function openFolderDialog() {
  closeMenus();
  const projectRoot = state.project?.root || null;
  try {
    await loadFolders(projectRoot);
  } catch (error) {
    if (!projectRoot) throw error;
    await loadFolders();
  }
  byId("folder-dialog").showModal();
}

async function selectFolder() {
  const path = state.folderPath;
  await api("/api/projects/open", {
    method: "POST",
    headers: originHeaders(),
    body: JSON.stringify({ path }),
  });
  byId("folder-dialog").close();
  await refreshProject();
  state.logs = [];
  await refreshRuns();
  toast(translate("Opened {path}", { path }));
}

async function refreshFiles() {
  const tree = byId("project-files");
  const expandedPaths = expandedDirectoryPaths(tree);
  const payload = await api("/api/files");
  renderFileEntries(payload.entries, tree);
  await restoreExpandedDirectories(tree, expandedPaths);
}

function expandedDirectoryPaths(container) {
  return [...container.querySelectorAll('.file-entry.directory[aria-expanded="true"]')]
    .map((row) => row.dataset.path)
    .filter(Boolean);
}

async function restoreExpandedDirectories(container, paths) {
  for (const path of paths) {
    const row = [...container.querySelectorAll(".file-entry.directory")]
      .find((candidate) => candidate.dataset.path === path);
    if (!row) continue;
    const disclosure = row.querySelector(".file-disclosure");
    const icon = row.querySelector(".file-kind-icon");
    const children = row.nextElementSibling;
    if (disclosure && icon && children?.classList.contains("file-children")) {
      await toggleDirectory(row, disclosure, icon, children, path);
    }
  }
}

function renderFileEntries(entries, container) {
  container.classList.remove("empty-state");
  container.removeAttribute("data-i18n");
  container.replaceChildren();
  entries.forEach((entry) => container.append(createFileNode(entry)));
  if (!entries.length) {
    const empty = document.createElement("div");
    empty.className = "file-empty";
    setTranslatedText(empty, "This folder is empty.");
    container.append(empty);
  }
}

function createFileNode(entry) {
  const node = document.createElement("div");
  node.className = "file-node";
  const isDirectory = entry.kind === "directory";
  const row = document.createElement(isDirectory ? "button" : "div");
  row.className = `file-entry${isDirectory ? " directory" : ""}`;
  if (isDirectory) row.type = "button";
  const disclosure = document.createElement("span");
  disclosure.className = "file-disclosure";
  disclosure.textContent = isDirectory ? "▸" : "";
  const icon = createFileIcon(isDirectory ? { icon: "folder", tone: "folder" } : filePresentation(entry.name));
  const name = document.createElement("span");
  name.className = "file-name";
  name.textContent = entry.name;
  row.append(disclosure, icon, name);
  node.append(row);

  if (isDirectory) {
    row.dataset.path = entry.path;
    const children = document.createElement("div");
    children.className = "file-children";
    children.hidden = true;
    children.setAttribute("role", "group");
    row.setAttribute("aria-expanded", "false");
    row.addEventListener("click", () => {
      toggleDirectory(row, disclosure, icon, children, entry.path).catch(reportError);
    });
    node.append(children);
  }
  return node;
}

function filePresentation(name) {
  const normalized = name.toLowerCase();
  const proteomics = PROTEOMICS_FILE_RULES.find((rule) => rule.pattern.test(normalized));
  if (proteomics) return proteomics;
  const special = {
    dockerfile: { icon: "vscode-config", tone: "config", filled: true },
    makefile: { icon: "file-code", tone: "code" },
    snakefile: { icon: "file-code", tone: "code" },
  }[normalized];
  if (special) return special;
  return FILE_ICON_RULES.find((rule) => rule.suffixes.some((suffix) => normalized.endsWith(suffix)))
    || { icon: "file", tone: "default" };
}

function createFileIcon(presentation) {
  const icon = document.createElementNS("http://www.w3.org/2000/svg", "svg");
  icon.classList.add("file-kind-icon", `file-kind-${presentation.tone}`);
  if (presentation.filled) icon.classList.add("file-icon-filled");
  icon.setAttribute("viewBox", "0 0 24 24");
  icon.setAttribute("aria-hidden", "true");
  const use = document.createElementNS("http://www.w3.org/2000/svg", "use");
  icon.append(use);
  setFileIcon(icon, presentation.icon);
  return icon;
}

function setFileIcon(icon, name) {
  icon.firstElementChild.setAttribute("href", `#file-icon-${name}`);
}

function collapseDirectory(row, disclosure, icon, children) {
  row.setAttribute("aria-expanded", "false");
  disclosure.textContent = "▸";
  setFileIcon(icon, "folder");
  children.hidden = true;
}

function collapseFileTree() {
  byId("project-files").querySelectorAll('.file-entry.directory[aria-expanded="true"]').forEach((row) => {
    const disclosure = row.querySelector(".file-disclosure");
    const icon = row.querySelector(".file-kind-icon");
    const children = row.nextElementSibling;
    if (disclosure && icon && children?.classList.contains("file-children")) {
      collapseDirectory(row, disclosure, icon, children);
    }
  });
}

async function toggleDirectory(row, disclosure, icon, children, path) {
  const expanded = row.getAttribute("aria-expanded") === "true";
  if (expanded) {
    collapseDirectory(row, disclosure, icon, children);
    return;
  }
  if (row.getAttribute("aria-busy") === "true") return;
  if (children.dataset.loaded !== "true") {
    row.setAttribute("aria-busy", "true");
    disclosure.textContent = "…";
    try {
      const payload = await api(`/api/files?path=${encodeURIComponent(path)}`);
      renderFileEntries(payload.entries, children);
      children.dataset.loaded = "true";
    } catch (error) {
      disclosure.textContent = "▸";
      throw error;
    } finally {
      row.removeAttribute("aria-busy");
    }
  }
  row.setAttribute("aria-expanded", "true");
  disclosure.textContent = "▾";
  setFileIcon(icon, "folder-open");
  children.hidden = false;
}

async function refreshCommands() {
  const payload = await api("/api/commands");
  state.commands = payload.commands;
  state.selectedCommand = null;
  const catalog = byId("command-catalog");
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
    button.textContent = command.display_name;
    button.title = command.help || "";
    button.addEventListener("click", () => selectCommand(index));
    actions.append(button);
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
  form.removeAttribute("data-i18n");
  form.replaceChildren();
  const flags = command.flags.filter((flag) => !flag.global && !flag.studio_hidden);
  const flagsById = new Map(flags.map((flag) => [flag.id, flag]));
  parameterGroups(command, flags).forEach((group) => {
    form.append(buildParameterGroup(group, flagsById));
  });
  configureConditionalParameters(command, form);
}

function parameterGroups(command, flags) {
  const groups = command.presentation?.groups;
  if (groups?.length) return groups;
  return [{
    id: "parameters",
    title: "Parameters",
    flags: flags.map((flag) => flag.id),
    common: flags.map((flag) => flag.id),
  }];
}

function buildParameterGroup(group, flagsById) {
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

function buildParameterFields(flags) {
  const fields = document.createElement("div");
  fields.className = "parameter-fields";
  flags.forEach((flag) => fields.append(buildArgumentField(flag)));
  return fields;
}

function buildAdvancedParameters(flags) {
  const details = document.createElement("details");
  details.className = "parameter-advanced";
  const summary = document.createElement("summary");
  setTranslatedText(summary, "Advanced");
  details.append(summary, buildParameterFields(flags));
  return details;
}

function buildArgumentField(flag) {
  const field = document.createElement("div");
  field.className = "form-field";
  field.dataset.flag = flag.long || flag.id;
  const heading = document.createElement("span");
  heading.className = "field-heading";
  renderFieldHeading(heading, flag.long || flag.id, Boolean(flag.required));
  const control = buildArgumentControl(flag);
  const help = document.createElement("small");
  if (flag.help) help.textContent = flag.help;
  else setValueHint(help, flag);
  field.append(heading, control, help);
  return field;
}

function renderFieldHeading(heading, flag, requirement) {
  const parameterSelector = heading.querySelector(".parameter-selector");
  heading.replaceChildren(parameterSelector || document.createTextNode(`--${flag}`));
  if (!requirement) return;
  const label = document.createElement("span");
  label.className = "required-label";
  setTranslatedText(label, requirement === true ? "Required" : requirement);
  heading.append(label);
}

function buildArgumentControl(flag) {
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
    setTranslatedText(remove, "Remove");
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

function configureConditionalParameters(command, form) {
  const path = command.path.join(" ");
  form.oninput = null;
  form.onchange = null;
  form.onclick = null;
  let update;
  if (path === "quantify features2proteins") {
    configureFeatures2ProteinsInputSelector(form);
    update = () => updateFeatures2ProteinsParameters(form);
  } else if (path === "quantify features2peptides") {
    configureFeatures2PeptidesIrsSelector(form);
    update = () => updateFeatures2PeptidesParameters(form);
  } else if (path === "quantify peptides2protein") {
    update = () => updatePeptides2ProteinParameters(form);
  } else if (path === "plot de") {
    update = () => updateDifferentialExpressionPlotParameters(form);
  } else {
    return;
  }
  form.oninput = update;
  form.onchange = update;
  form.onclick = (event) => {
    if (event.target.closest(".add-value, .remove-value")) update();
  };
  update();
}

function updateDifferentialExpressionPlotParameters(form) {
  const volcano = argumentControlValue(form, "volcano") === "true";
  const heatmap = argumentControlValue(form, "heatmap") === "true";
  setConditionalFields(form, ["sdrf"], heatmap, true);
  setArgumentFieldRequired(form, "sdrf", heatmap);
  setConditionalFields(form, ["highlight-protein"], volcano);
  updateAdvancedVisibility(form);
}

function configureFeatures2ProteinsInputSelector(form) {
  configureAlternativeParameterSelector(
    form,
    ["parquet", "msstats"],
    "featureInputType",
    "Input file type",
  );
  const helpByFlag = {
    parquet: "Feature-level QPX input for protein quantification; also supplies protein-group assignments for spectral-count.",
    msstats: "Feature-level MSstats input for protein quantification; requires SDRF metadata.",
  };
  Object.entries(helpByFlag).forEach(([flag, help]) => {
    const field = [...form.querySelectorAll(".form-field")]
      .find((candidate) => candidate.dataset.flag === flag);
    if (field) setTranslatedText(field.querySelector("small"), help);
  });
}

function configureFeatures2PeptidesIrsSelector(form) {
  configureAlternativeParameterSelector(
    form,
    ["irs-channel", "irs-autodetect-regex"],
    "peptideIrsType",
    "IRS reference type",
  );
}

function configureAlternativeParameterSelector(form, flags, stateKey, label) {
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

function setAlternativeParameter(form, flags, stateKey, selected) {
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

function updateFeatures2ProteinsParameters(form) {
  const quantMethod = argumentControlValue(form, "quant-method", true).toLowerCase();
  updateFeatures2ProteinsQuantification(form, quantMethod);
  updateFeatures2ProteinsCorrections(form, quantMethod);
  updateFeatures2ProteinsDownstream(form, quantMethod);
  updateFeatures2ProteinsRequirements(form, quantMethod);
  updateAdvancedVisibility(form);
}

function updateFeatures2ProteinsQuantification(form, quantMethod) {
  const pibaq = quantMethod === "pibaq";
  const ratio = quantMethod === "ratio";
  const spectralCount = quantMethod === "spectral-count";
  const qpxOnly = ratio || spectralCount;
  if (qpxOnly) {
    setAlternativeParameter(
      form,
      ["parquet", "msstats"],
      "featureInputType",
      "parquet",
    );
  }
  form.querySelectorAll('[data-parameter-selector="featureInputType"]').forEach((selector) => {
    selector.querySelector(".parameter-selector-trigger").disabled = qpxOnly;
    selector.classList.toggle("disabled", qpxOnly);
    if (qpxOnly && selector.classList.contains("open")) closeParameterSelectors();
  });
  setConditionalFields(form, ["fasta", "pibaq-enzyme"], pibaq, true);
  setConditionalFields(form, [
    "pibaq-max-aa",
    "pibaq-min-shared",
    "pibaq-families",
    "pibaq-min-anchors",
  ], pibaq);
  setArgumentFieldRequired(form, "fasta", pibaq);
  setConditionalFields(form, ["min-unique"], !pibaq);
  setConditionalFields(form, ["psm"], spectralCount);
  setConditionalFields(form, ["ratio-fraction-merge"], ratio);
  setConditionalFields(form, ["directlfq-min-nonan"], quantMethod === "directlfq");
  setConditionalFields(form, ["export-ions"], quantMethod === "directlfq");
  setConditionalFields(
    form,
    ["export-peptides"],
    !["directlfq", "ratio", "spectral-count"].includes(quantMethod),
  );
  setConditionalFields(
    form,
    ["directlfq-num-samples-quadratic"],
    ["directlfq", "maxlfq"].includes(quantMethod),
  );
  setConditionalFields(
    form,
    ["normalization-proteins"],
    !["directlfq", "ratio", "peptide-count", "spectral-count"].includes(quantMethod),
  );
  setConditionalFields(
    form,
    ["run-normalization", "sample-normalization"],
    !["directlfq", "ratio", "peptide-count", "spectral-count"].includes(quantMethod),
  );
}

function updateFeatures2ProteinsCorrections(form, quantMethod) {
  const ratio = quantMethod === "ratio";
  const batchCorrection = argumentControlValue(form, "batch-correction") === "true";
  const batchMethod = argumentControlValue(form, "batch-method", true).toLowerCase() || "sample-prefix";
  const hasBatchCovariate = Boolean(argumentControlValue(form, "batch-covariate"));
  setConditionalFields(form, ["batch-method", "batch-covariate"], batchCorrection, true);
  setConditionalFields(form, ["batch-column"], batchCorrection && batchMethod === "column", true);
  setArgumentFieldRequired(form, "batch-column", batchCorrection && batchMethod === "column");
  setConditionalFields(form, ["batch-nonparametric", "batch-mean-only", "batch-ref"], batchCorrection);

  const irsAvailable = !["ratio", "peptide-count", "spectral-count"].includes(quantMethod);
  const irs = irsAvailable && argumentControlValue(form, "irs") === "true";
  setConditionalFields(form, ["irs"], irsAvailable);
  setConditionalFields(form, ["irs-reference-sample", "irs-reference-regex"], ratio || irs, true);
  setConditionalFields(form, ["irs-sdrf-column", "irs-sdrf-value"], irs, true);
  setConditionalFields(form, ["irs-stat", "irs-remove-reference"], irs);
}

function updateFeatures2ProteinsDownstream(form, quantMethod) {
  const imputeMethod = argumentControlValue(form, "impute-method").toLowerCase();
  setConditionalFields(form, ["impute-quantile"], ["mindet", "minprob"].includes(imputeMethod), true);
  setConditionalFields(form, ["impute-shift", "impute-scale"], imputeMethod === "minprob", true);
  setConditionalFields(form, ["impute-n-neighbors"], ["knn", "seqknn"].includes(imputeMethod), true);

  const differentialExpression = Boolean(
    argumentControlValue(form, "de-contrast") || argumentControlValue(form, "de-contrast-file"),
  );
  setConditionalFields(form, ["de-method", "de-log2fc", "de-fdr", "de-output"], differentialExpression, true);
  setConditionalFields(form, ["de-effect-size-gate"], differentialExpression);
  const deMethod = argumentControlValue(form, "de-method", true).toLowerCase() || "auto";
  const effectiveDeMethod = deMethod === "auto"
    ? (quantMethod === "directlfq" ? "deqms" : "limrots")
    : deMethod;
  setConditionalFields(
    form,
    ["de-fdr-method"],
    differentialExpression && !["limrots", "rots"].includes(effectiveDeMethod),
  );
  setConditionalFields(
    form,
    ["de-ensemble-method", "de-ensemble-min-k"],
    differentialExpression && deMethod === "ensemble",
  );
}

function updateFeatures2ProteinsRequirements(form, quantMethod) {
  const spectralCount = quantMethod === "spectral-count";
  const qpxOnly = quantMethod === "ratio" || spectralCount;
  const inputType = form.dataset.featureInputType || "parquet";
  const hasMsstats = inputType === "msstats";
  const sampleNormalization = argumentControlValue(form, "sample-normalization", true).toLowerCase();
  const batchCorrection = argumentControlValue(form, "batch-correction") === "true";
  const batchMethod = argumentControlValue(form, "batch-method", true).toLowerCase() || "sample-prefix";
  const hasBatchCovariate = Boolean(argumentControlValue(form, "batch-covariate"));
  const differentialExpression = Boolean(
    argumentControlValue(form, "de-contrast") || argumentControlValue(form, "de-contrast-file"),
  );
  const requiresSdrf = hasMsstats
    || qpxOnly
    || argumentControlValue(form, "irs") === "true"
    || Boolean(argumentControlValue(form, "coverage-threshold"))
    || Boolean(argumentControlValue(form, "min-sample-correlation"))
    || differentialExpression
    || sampleNormalization === "condition-median"
    || (batchCorrection && (batchMethod === "column" || hasBatchCovariate));

  setArgumentFieldRequired(form, "parquet", inputType === "parquet");
  setArgumentFieldRequired(form, "msstats", hasMsstats);
  setArgumentFieldRequired(form, "psm", spectralCount);
  setArgumentFieldRequired(form, "sdrf", requiresSdrf);
  setArgumentFieldRequired(form, "de-output", differentialExpression);
}

function updateFeatures2PeptidesParameters(form) {
  setArgumentFieldRequired(form, "parquet", true);
  setArgumentFieldRequired(form, "output", true);
  const skipNormalization = argumentControlValue(form, "skip-normalization") === "true";
  setArgumentFieldsEnabled(form, ["run-normalization", "sample-normalization"], !skipNormalization);
  const sampleNormalization = argumentControlValue(form, "sample-normalization", true).toLowerCase();
  const irsAutodetect = argumentControlValue(form, "irs-autodetect-regex");
  setArgumentFieldRequired(
    form,
    "sdrf",
    Boolean(irsAutodetect) || sampleNormalization === "condition-median",
  );
  const irs = Boolean(
    argumentControlValue(form, "irs-channel") || argumentControlValue(form, "irs-autodetect-regex"),
  );
  setConditionalFields(form, ["irs-stat", "irs-scope"], irs);
  updateAdvancedVisibility(form);
}

function updatePeptides2ProteinParameters(form) {
  const quantMethod = argumentControlValue(form, "quant-method", true).toLowerCase();
  const pibaq = quantMethod === "pibaq";
  setConditionalFields(form, ["fasta", "enzyme"], pibaq, true);
  setArgumentFieldRequired(form, "fasta", pibaq);
  setConditionalFields(form, [
    "min-aa",
    "max-aa",
    "families",
    "min-shared",
    "min-anchors",
    "high-anchor-threshold",
  ], pibaq);
  setConditionalFields(form, ["tpa"], pibaq);
  const tpa = pibaq && argumentControlValue(form, "tpa") === "true";
  setConditionalFields(form, ["ruler"], tpa);
  const ruler = tpa && argumentControlValue(form, "ruler") === "true";
  setConditionalFields(form, ["ploidy", "organism", "cpc"], ruler);
  setConditionalFields(form, ["qc-report"], pibaq);
  setConditionalFields(form, ["threads"], ["directlfq", "maxlfq"].includes(quantMethod));
  setConditionalFields(form, ["directlfq-min-nonan"], quantMethod === "directlfq");
  updateAdvancedVisibility(form);
}

function argumentControlValue(form, flag, useDefault = false) {
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

function setConditionalFields(form, flags, visible, promote = false) {
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

function setArgumentFieldsEnabled(form, flags, enabled) {
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

function setArgumentFieldRequired(form, flag, required) {
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

function promoteParameterField(field) {
  const section = field.closest(".parameter-group");
  const primary = [...section.children]
    .find((child) => child.classList.contains("parameter-primary"));
  if (primary && field.parentElement !== primary) primary.append(field);
}

function clearArgumentControl(control) {
  control.querySelectorAll('input[type="checkbox"]').forEach((input) => {
    input.checked = false;
  });
  control.querySelectorAll("[data-value]").forEach((input) => {
    input.value = "";
  });
  [...control.querySelectorAll(".value-row")].slice(1).forEach((row) => row.remove());
}

function updateAdvancedVisibility(form) {
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

function setValueHint(element, flag) {
  const names = flag.value_names || [];
  if (!names.length) {
    setTranslatedText(element, "Optional switch");
    return;
  }
  const hint = names.map((name) => `<${name}>`).join(" ");
  if (flag.repeat) setTranslatedText(element, "{hint}; may be repeated", { hint });
  else element.textContent = hint;
}

function numericValue(flag) {
  return (flag.value_names || []).some((name) => ["N", "VALUE", "FRACTION", "CORRELATION"].includes(name));
}

function commandArgv() {
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

function workflowTemplateParameters() {
  const parameters = {};
  byId("command-form").querySelectorAll(".argument-control").forEach((control) => {
    const value = argumentTemplateValue(control);
    if (value !== null) parameters[control.dataset.flag] = value;
  });
  return parameters;
}

function argumentTemplateValue(control) {
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

function currentWorkflowTemplate() {
  if (!state.selectedCommand) throw new Error(translate("Select a workflow first"));
  return {
    $schemaVersion: WORKFLOW_TEMPLATE_SCHEMA_VERSION,
    workflow: [...state.selectedCommand.path],
    parameters: workflowTemplateParameters(),
  };
}

function workflowTemplatePath(directory, name) {
  return directory === "." ? name : `${directory}/${name}`;
}

function workflowTemplateDirectoryLabel(directory) {
  if (directory === ".") return state.project.root;
  return `${state.project.root.replace(/[\\/]$/, "")}/${directory}`;
}

function updateWorkflowTemplateConfirm() {
  const confirm = byId("workflow-template-confirm");
  confirm.disabled = state.workflowTemplateMode === "import"
    ? !state.workflowTemplateSelected
    : !byId("workflow-template-name").value.trim();
}

function renderWorkflowTemplateEntries() {
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

function selectWorkflowTemplateFile(entry) {
  state.workflowTemplateSelected = entry.path;
  byId("workflow-template-name").value = entry.name;
  byId("workflow-template-list").querySelectorAll(".workflow-template-row:not(.directory)")
    .forEach((row) => row.setAttribute("aria-pressed", String(row.dataset.path === entry.path)));
  updateWorkflowTemplateConfirm();
}

async function loadWorkflowTemplateDirectory(path = ".") {
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

async function openWorkflowTemplateDialog(mode) {
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

async function submitWorkflowTemplate(event) {
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

function applyWorkflowTemplate(template) {
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

function applyWorkflowParameterPatch(patch) {
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

function highlightAgentFields(fields) {
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

function normalizeTemplateParameter(flag, value) {
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

function templateScalar(value, parameter) {
  if (!(["string", "number"].includes(typeof value))) throw invalidTemplateValue(parameter);
  const text = String(value).trim();
  if (!text) throw invalidTemplateValue(parameter);
  return text;
}

function invalidTemplateValue(parameter) {
  return new Error(translate("Invalid template value for --{parameter}", { parameter }));
}

function fillTemplateControl(control, flag, value) {
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

async function validateCommand() {
  const payload = await api("/api/commands/validate", {
    method: "POST",
    headers: originHeaders(),
    body: JSON.stringify({ argv: commandArgv() }),
  });
  showBottomTab("logs");
  state.logs.push(translate("Validated: {command}", { command: payload.argv.join(" ") }));
  renderBottom();
  toast(translate("Parameters are valid"));
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
  toast(translate("Run {id} queued", { id: payload.id }));
  await refreshRuns();
  showBottomTab("runs");
}

async function refreshRuns() {
  const payload = await api("/api/runs");
  state.runs = payload.runs;
  const active = state.runs.find((run) => ACTIVE_RUN_STATUSES.has(run.status)) || null;
  state.activeRunId = active?.id || null;
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
  if (!state.activeRunId) throw new Error(translate("No run is active"));
  const runId = state.activeRunId;
  await api(`/api/runs/${encodeURIComponent(runId)}/cancel`, {
    method: "POST",
    headers: originHeaders(),
  });
  toast(translate("Cancelled run {id}", { id: runId }));
  await refreshRuns();
}

async function refreshArtifacts() {
  const payload = await api("/api/artifacts");
  state.artifacts = payload.artifacts;
  if (state.bottomTab === "artifacts") renderBottom();
}

function renderBottom() {
  const content = byId("bottom-content");
  content.removeAttribute("data-i18n");
  content.replaceChildren();
  if (state.bottomTab === "runs") {
    content.textContent = state.runs.length
      ? state.runs.map((run) => `${translate(run.status).padEnd(12)} ${run.command}  ${run.id}`).join("\n")
      : translate("No runs yet.");
    return;
  }
  if (state.bottomTab === "logs") {
    content.textContent = state.logs.length ? state.logs.join("\n") : translate("No log events yet.");
    content.scrollTop = content.scrollHeight;
    return;
  }
  if (!state.artifacts.length) {
    content.textContent = translate("No artifacts yet.");
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

function toggleSidebar() {
  const workspace = document.querySelector(".workspace");
  if (window.matchMedia("(max-width: 680px)").matches) {
    if (!workspace.classList.contains("sidebar-mobile-visible")) {
      workspace.classList.remove("assistant-mobile-visible");
    }
    workspace.classList.toggle("sidebar-mobile-visible");
    applyPanelSizes();
    return;
  }
  workspace.classList.toggle("sidebar-hidden");
  applyPanelSizes();
}

function toggleAssistant() {
  const workspace = document.querySelector(".workspace");
  if (window.matchMedia("(max-width: 1050px)").matches) {
    if (!workspace.classList.contains("assistant-mobile-visible")) {
      workspace.classList.remove("sidebar-mobile-visible");
    }
    workspace.classList.toggle("assistant-mobile-visible");
    applyPanelSizes();
    return;
  }
  workspace.classList.toggle("assistant-hidden");
  applyPanelSizes();
}

function setSidePanelCollapsed(panel, collapsed) {
  const workspace = document.querySelector(".workspace");
  if (panel === "sidebar") {
    if (window.matchMedia("(max-width: 680px)").matches) {
      workspace.classList.toggle("sidebar-mobile-visible", !collapsed);
    } else {
      workspace.classList.toggle("sidebar-hidden", collapsed);
    }
  }
  if (panel === "assistant") {
    if (window.matchMedia("(max-width: 1050px)").matches) {
      workspace.classList.toggle("assistant-mobile-visible", !collapsed);
    } else {
      workspace.classList.toggle("assistant-hidden", collapsed);
    }
  }
  applyPanelSizes();
}

function panelSizeBounds(panel) {
  const limits = PANEL_SIZE_LIMITS[panel];
  if (panel === "bottom") {
    const headerHeight = document.querySelector(".app-header").getBoundingClientRect().height;
    const available = window.innerHeight - headerHeight - 240;
    return { min: limits.min, max: Math.max(limits.min, Math.min(limits.max, available)) };
  }

  const workspace = document.querySelector(".workspace");
  const viewportMax = Math.floor(window.innerWidth * SIDE_PANEL_MAX_VIEWPORT_RATIO);
  if (panel === "assistant" && window.matchMedia("(max-width: 1050px)").matches) {
    return { min: Math.min(limits.min, viewportMax), max: viewportMax };
  }
  if (panel === "sidebar" && window.matchMedia("(max-width: 680px)").matches) {
    return { min: Math.min(limits.min, viewportMax), max: viewportMax };
  }

  let occupied = 0;
  if (panel === "sidebar" && !window.matchMedia("(max-width: 1050px)").matches
      && !workspace.classList.contains("assistant-hidden")) {
    occupied = state.panelSizes.assistant;
  }
  if (panel === "assistant" && !workspace.classList.contains("sidebar-hidden")) {
    occupied = state.panelSizes.sidebar;
  }
  const available = workspace.clientWidth - MIN_WORKFLOW_WIDTH - occupied;
  const max = Math.max(0, Math.min(viewportMax, available));
  return { min: Math.min(limits.min, max), max };
}

function setPanelSize(panel, requestedSize) {
  const bounds = panelSizeBounds(panel);
  const size = Math.round(Math.min(Math.max(requestedSize, bounds.min), bounds.max));
  state.panelSizes[panel] = size;
  if (panel === "bottom") {
    if (!byId("bottom-panel").classList.contains("collapsed")) {
      document.documentElement.style.setProperty("--bottom-height", `${size}px`);
    }
  } else {
    document.documentElement.style.setProperty(`--${panel}-width`, `${size}px`);
  }
  const handle = document.querySelector(`[data-resize-panel="${panel}"]`);
  handle?.setAttribute("aria-valuemin", String(Math.round(bounds.min)));
  handle?.setAttribute("aria-valuemax", String(Math.round(bounds.max)));
  handle?.setAttribute("aria-valuenow", String(size));
}

function applyPanelSizes() {
  setPanelSize("assistant", state.panelSizes.assistant);
  setPanelSize("sidebar", state.panelSizes.sidebar);
  setPanelSize("assistant", state.panelSizes.assistant);
  setPanelSize("bottom", state.panelSizes.bottom);
  resizeWorkspace();
}

function bindPanelResizers() {
  document.querySelectorAll("[data-resize-panel]").forEach((handle) => {
    const panel = handle.dataset.resizePanel;
    handle.addEventListener("pointerdown", (event) => {
      if (event.button !== 0) return;
      event.preventDefault();
      closeMenus();
      if (panel === "bottom") {
        byId("bottom-panel").classList.remove("collapsed");
        resizeWorkspace();
      }
      const startX = event.clientX;
      const startY = event.clientY;
      const startSize = state.panelSizes[panel];
      const cursorClass = panel === "bottom" ? "resize-rows" : "resize-columns";
      const dragSurface = document.querySelector(".workspace");
      const collapseThreshold = panel === "bottom"
        ? null
        : panelSizeBounds(panel).min * SIDE_PANEL_COLLAPSE_RATIO;
      let collapsedDuringDrag = false;
      handle.setPointerCapture(event.pointerId);
      handle.classList.add("active");
      document.body.classList.add("resizing", cursorClass);

      const move = (moveEvent) => {
        let size = startSize;
        if (panel === "sidebar") size += moveEvent.clientX - startX;
        if (panel === "assistant") size -= moveEvent.clientX - startX;
        if (panel === "bottom") size -= moveEvent.clientY - startY;
        const shouldCollapse = collapseThreshold !== null && size <= collapseThreshold;
        if (shouldCollapse !== collapsedDuringDrag) {
          collapsedDuringDrag = shouldCollapse;
          if (shouldCollapse) {
            state.panelSizes[panel] = startSize;
            dragSurface.setPointerCapture(moveEvent.pointerId);
          }
          setSidePanelCollapsed(panel, shouldCollapse);
        }
        if (shouldCollapse) return;
        setPanelSize(panel, size);
      };
      const finish = (finishEvent) => {
        if (handle.hasPointerCapture(finishEvent.pointerId)) {
          handle.releasePointerCapture(finishEvent.pointerId);
        }
        if (dragSurface.hasPointerCapture(finishEvent.pointerId)) {
          dragSurface.releasePointerCapture(finishEvent.pointerId);
        }
        window.removeEventListener("pointermove", move);
        window.removeEventListener("pointerup", finish);
        window.removeEventListener("pointercancel", finish);
        handle.classList.remove("active");
        document.body.classList.remove("resizing", cursorClass);
      };
      window.addEventListener("pointermove", move);
      window.addEventListener("pointerup", finish);
      window.addEventListener("pointercancel", finish);
    });
    handle.addEventListener("dblclick", () => {
      if (panel === "bottom") byId("bottom-panel").classList.remove("collapsed");
      setPanelSize(panel, DEFAULT_PANEL_SIZES[panel]);
      resizeWorkspace();
    });
    handle.addEventListener("keydown", (event) => {
      const bounds = panelSizeBounds(panel);
      let size = state.panelSizes[panel];
      if (event.key === "Home") size = bounds.min;
      else if (event.key === "End") size = bounds.max;
      else if (panel === "sidebar" && event.key === "ArrowLeft") size -= 10;
      else if (panel === "sidebar" && event.key === "ArrowRight") size += 10;
      else if (panel === "assistant" && event.key === "ArrowLeft") size += 10;
      else if (panel === "assistant" && event.key === "ArrowRight") size -= 10;
      else if (panel === "bottom" && event.key === "ArrowUp") size += 10;
      else if (panel === "bottom" && event.key === "ArrowDown") size -= 10;
      else return;
      event.preventDefault();
      if (panel === "bottom") byId("bottom-panel").classList.remove("collapsed");
      setPanelSize(panel, size);
      resizeWorkspace();
    });
  });
  window.addEventListener("resize", applyPanelSizes);
}

function resizeWorkspace() {
  const collapsed = byId("bottom-panel").classList.contains("collapsed");
  const height = collapsed ? 39 : state.panelSizes.bottom;
  document.documentElement.style.setProperty("--bottom-height", `${height}px`);
}

function shortcutHelp() {
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

async function performAction(action) {
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

function reportError(error) {
  if (error?.assistantDisplayed) return;
  toast(error.message || String(error), true);
}

function bindMenuEvents() {
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

function navigateMenu(event, menu, triggers) {
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

function navigateSubmenu(event, menu) {
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

function bindEvents() {
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
  document.addEventListener("keydown", keyboardShortcut);
}

function keyboardShortcut(event) {
  if (event.defaultPrevented) return;
  if (event.target.closest?.("input, textarea, select, dialog")) return;
  const control = event.ctrlKey || event.metaKey;
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

async function boot() {
  applyAppearance(state.appearance, false);
  translateInterface(state.language, false);
  applyPanelSizes();
  const session = await api("/api/session");
  state.csrf = session.csrf_token;
  state.version = session.version;
  state.provider = session.ai_provider;
  bindEvents();
  await Promise.all([refreshProject(), refreshRuns(), refreshArtifacts(), refreshSystemMemory()]);
  renderProviderState();
  await restorePendingApproval();
  renderBottom();
  window.setInterval(() => void refreshSystemMemory(), SYSTEM_MEMORY_REFRESH_MS);
}

boot().catch(reportError);
