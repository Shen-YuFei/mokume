// Workflow-specific rules that show, hide, and require parameters as values change.
import { setTranslatedText } from "./studio-i18n.js";
import { closeParameterSelectors } from "./studio-menus.js";
import {
  argumentControlValue,
  configureAlternativeParameterSelector,
  setAlternativeParameter,
  setArgumentFieldRequired,
  setArgumentFieldsEnabled,
  setConditionalFields,
  updateAdvancedVisibility,
} from "./studio-form.js";

export function configureConditionalParameters(command, form) {
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

export function updateDifferentialExpressionPlotParameters(form) {
  const volcano = argumentControlValue(form, "volcano") === "true";
  const heatmap = argumentControlValue(form, "heatmap") === "true";
  setConditionalFields(form, ["sdrf"], heatmap, true);
  setArgumentFieldRequired(form, "sdrf", heatmap);
  setConditionalFields(form, ["highlight-protein"], volcano);
  updateAdvancedVisibility(form);
}

export function configureFeatures2ProteinsInputSelector(form) {
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

export function configureFeatures2PeptidesIrsSelector(form) {
  configureAlternativeParameterSelector(
    form,
    ["irs-channel", "irs-autodetect-regex"],
    "peptideIrsType",
    "IRS reference type",
  );
}

export function updateFeatures2ProteinsParameters(form) {
  const quantMethod = argumentControlValue(form, "quant-method", true).toLowerCase();
  updateFeatures2ProteinsQuantification(form, quantMethod);
  updateFeatures2ProteinsCorrections(form, quantMethod);
  updateFeatures2ProteinsDownstream(form, quantMethod);
  updateFeatures2ProteinsRequirements(form, quantMethod);
  updateAdvancedVisibility(form);
}

export function updateFeatures2ProteinsQuantification(form, quantMethod) {
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

export function updateFeatures2ProteinsCorrections(form, quantMethod) {
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

export function updateFeatures2ProteinsDownstream(form, quantMethod) {
  const imputeMethod = argumentControlValue(form, "impute-method").toLowerCase();
  setConditionalFields(form, ["impute-quantile"], ["mindet", "minprob"].includes(imputeMethod), true);
  setConditionalFields(form, ["impute-seed", "impute-tune-sigma"], ["minprob", "qrilc"].includes(imputeMethod), true);
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

export function updateFeatures2ProteinsRequirements(form, quantMethod) {
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

export function updateFeatures2PeptidesParameters(form) {
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

export function updatePeptides2ProteinParameters(form) {
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
