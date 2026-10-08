// Folder picker, project file tree, and file-type icons.
import { api, byId, formatBytes, originHeaders, reportError, state, toast } from "./studio-core.js";
import { setTranslatedAttribute, setTranslatedText, translate } from "./studio-i18n.js";
import { closeMenus } from "./studio-menus.js";
import { refreshRuns } from "./studio-runs.js";
import { refreshProject } from "./studio.js";

export const FILE_PATH_TYPE = "application/x-mokume-path";

export const PROTEOMICS_FILE_RULES = Object.freeze([
  Object.freeze({ pattern: /(^|[._-])sdrf([._-]|$)/, icon: "sdrf", tone: "sdrf" }),
  Object.freeze({ pattern: /(^|[._-])msstats([._-]|$)/, icon: "msstats", tone: "stats" }),
]);

export const FILE_ICON_RULES = Object.freeze([
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

export async function loadFolders(path = null) {
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

export async function openFolderDialog() {
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

export async function selectFolder() {
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

export async function refreshFiles() {
  const tree = byId("project-files");
  const expandedPaths = expandedDirectoryPaths(tree);
  const payload = await api("/api/files");
  renderFileEntries(payload.entries, tree);
  await restoreExpandedDirectories(tree, expandedPaths);
}

export function expandedDirectoryPaths(container) {
  return [...container.querySelectorAll('.file-entry.directory[aria-expanded="true"]')]
    .map((row) => row.dataset.path)
    .filter(Boolean);
}

export async function restoreExpandedDirectories(container, paths) {
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

export function renderFileEntries(entries, container) {
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

export function createFileNode(entry) {
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
  } else {
    row.draggable = true;
    setTranslatedAttribute(row, "title", "{path} ({size}). Drag onto a parameter to use this path.", {
      path: entry.path,
      size: formatBytes(entry.size),
    });
    row.addEventListener("dragstart", (event) => {
      event.dataTransfer.setData(FILE_PATH_TYPE, entry.path);
      event.dataTransfer.setData("text/plain", entry.path);
      event.dataTransfer.effectAllowed = "copy";
    });
  }
  return node;
}

// A file dropped on a text parameter replaces its value with the project-relative path.
export function bindFileDrops() {
  const form = byId("command-form");
  const dropTarget = (event) => {
    const input = event.target.closest?.('input[type="text"]');
    if (!input || input.disabled || !event.dataTransfer.types.includes(FILE_PATH_TYPE)) return null;
    return input;
  };
  form.addEventListener("dragover", (event) => {
    const input = dropTarget(event);
    if (!input) return;
    event.preventDefault();
    event.dataTransfer.dropEffect = "copy";
    input.classList.add("drop-target");
  });
  form.addEventListener("dragleave", (event) => event.target.classList?.remove("drop-target"));
  form.addEventListener("drop", (event) => {
    const input = dropTarget(event);
    if (!input) return;
    event.preventDefault();
    input.classList.remove("drop-target");
    input.value = event.dataTransfer.getData(FILE_PATH_TYPE);
    input.dispatchEvent(new Event("input", { bubbles: true }));
    input.dispatchEvent(new Event("change", { bubbles: true }));
    input.focus();
  });
}

export function filePresentation(name) {
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

export function createFileIcon(presentation) {
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

export function setFileIcon(icon, name) {
  icon.firstElementChild.setAttribute("href", `#file-icon-${name}`);
}

export function collapseDirectory(row, disclosure, icon, children) {
  row.setAttribute("aria-expanded", "false");
  disclosure.textContent = "▸";
  setFileIcon(icon, "folder");
  children.hidden = true;
}

export function collapseFileTree() {
  byId("project-files").querySelectorAll('.file-entry.directory[aria-expanded="true"]').forEach((row) => {
    const disclosure = row.querySelector(".file-disclosure");
    const icon = row.querySelector(".file-kind-icon");
    const children = row.nextElementSibling;
    if (disclosure && icon && children?.classList.contains("file-children")) {
      collapseDirectory(row, disclosure, icon, children);
    }
  });
}

export async function toggleDirectory(row, disclosure, icon, children, path) {
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
