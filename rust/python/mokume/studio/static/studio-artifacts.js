// Shared artifact links, image thumbnails, and sandboxed HTML previews.
import { formatBytes, state } from "./studio-core.js";
import { setTranslatedAttribute, translatedElement } from "./studio-i18n.js";

export function relativeArtifactPath(path) {
  const root = state.project?.root?.replace(/\/$/, "");
  return root && path.startsWith(`${root}/`) ? path.slice(root.length + 1) : path;
}

export function artifactItem(artifact) {
  const item = document.createElement("article");
  item.className = "artifact-item";
  const link = document.createElement("a");
  link.className = "artifact-link";
  link.href = `/api/artifacts/${encodeURIComponent(artifact.id)}`;
  link.target = "_blank";
  link.rel = "noopener";
  link.title = artifact.path;
  const relativePath = relativeArtifactPath(artifact.path);
  if (/\.(png|jpe?g|gif|webp|svg)$/i.test(artifact.path)) {
    const image = document.createElement("img");
    image.className = "artifact-thumbnail";
    image.src = link.href;
    image.alt = relativePath;
    image.loading = "lazy";
    link.append(image);
  }
  const description = document.createElement("span");
  description.className = "artifact-description";
  const path = document.createElement("code");
  path.textContent = relativePath;
  const meta = document.createElement("small");
  meta.textContent = `${formatBytes(artifact.size)} · sha256 ${artifact.sha256}`;
  description.append(path, meta);
  link.append(description);
  item.append(link);
  if (artifact.media_type === "text/html") item.append(htmlPreview(artifact));
  return item;
}

export function htmlPreview(artifact) {
  const preview = document.createElement("details");
  preview.className = "artifact-preview";
  preview.append(translatedElement("summary", "", "Preview report"));
  preview.addEventListener("toggle", () => {
    if (!preview.open || preview.querySelector("iframe")) return;
    const frame = document.createElement("iframe");
    frame.setAttribute("sandbox", "allow-scripts");
    setTranslatedAttribute(frame, "title", "HTML report preview");
    frame.src = `/api/artifacts/${encodeURIComponent(artifact.id)}/preview`;
    preview.append(frame);
  });
  return preview;
}
