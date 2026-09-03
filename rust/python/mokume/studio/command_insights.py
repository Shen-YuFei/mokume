"""Parse Studio commands for review, templates, and scientific step labels."""

from __future__ import annotations

from typing import Any


def planned_workflow_steps(argv: list[str]) -> list[str]:
    """Describe configured scientific steps without claiming runtime callbacks."""
    path = tuple(argv[:2]) if argv[:1] == ["quantify"] else tuple(argv[:1])
    options = set(token.split("=", 1)[0] for token in argv if token.startswith("--"))
    if path == ("quantify", "features2proteins"):
        steps = _features_to_proteins_steps(argv, options)
    elif path == ("quantify", "features2peptides"):
        middle = [] if "--skip-normalization" in options else ["normalize"]
        steps = ["read", "filter", *middle, "aggregate", "export"]
    elif path == ("quantify", "peptides2protein"):
        middle = ["normalize"] if "--normalize" in options else []
        steps = ["read", "aggregate", *middle, "export"]
    elif path == ("correct-batches",):
        steps = ["read", "correct", "export"]
    elif path == ("tissuemap",):
        steps = ["read", "correct", "impute", "embed", "export"]
    elif path in {("plot", "pca"), ("plot", "tsne")}:
        steps = ["read", "embed", "export"]
    elif path in {("plot", "de"), ("interactive-report",)}:
        steps = ["read", "differential", "visualize", "export"]
    else:
        steps = ["read", "analyze", "export"]
    return steps


def parse_occurrences(argv: list[str], spec: dict[str, Any]) -> list[dict[str, Any]]:
    """Return supplied command arguments paired with their catalog flags."""
    options = {}
    for flag in spec.get("flags", ()):
        if flag.get("long"):
            options[f"--{flag['long']}"] = flag
        if flag.get("short"):
            options[f"-{flag['short']}"] = flag
    occurrences = []
    index = len(spec["path"])
    while index < len(argv):
        option, inline = _split_option(argv[index])
        flag = options.get(option)
        if flag is None:
            index += 1
            continue
        count = value_count(flag)
        values = ([inline] if inline is not None else []) + argv[
            index + 1 : index + 1 + count - (inline is not None)
        ]
        occurrences.append({"flag": flag, "option": option, "values": values})
        index += 1 + count - (inline is not None)
    return occurrences


def value_count(flag: dict[str, Any]) -> int:
    """Return the number of values consumed by a catalog flag."""
    arity = flag.get("value_arity") or {}
    maximum = arity.get("max")
    minimum = int(arity.get("min", 0))
    if maximum == 0 or not flag.get("value_names"):
        return 0
    return max(1, minimum) if maximum is None or maximum != minimum else int(maximum)


def template_from_argv(
    argv: list[str], spec: dict[str, Any] | None
) -> dict[str, Any] | None:
    """Convert canonical argv into a reusable workflow template."""
    if spec is None:
        return None
    parameters: dict[str, Any] = {}
    grouped: dict[str, list[list[str]]] = {}
    booleans: set[str] = set()
    flags: dict[str, dict[str, Any]] = {}
    for occurrence in parse_occurrences(argv, spec):
        name = str(occurrence["flag"].get("long") or occurrence["flag"].get("id"))
        flags[name] = occurrence["flag"]
        if value_count(occurrence["flag"]) == 0:
            booleans.add(name)
        else:
            grouped.setdefault(name, []).append(occurrence["values"])
    parameters.update({name: True for name in booleans})
    for name, rows in grouped.items():
        repeated = bool(flags[name].get("repeat"))
        single = value_count(flags[name]) == 1
        if single:
            parameters[name] = [row[0] for row in rows] if repeated else rows[-1][0]
        else:
            parameters[name] = rows if repeated else rows[-1]
    return {
        "$schemaVersion": 1,
        "workflow": list(spec["path"]),
        "parameters": parameters,
    }


def _features_to_proteins_steps(argv: list[str], options: set[str]) -> list[str]:
    steps = ["read", "aggregate"]
    quant_method = _argv_value(argv, "quant-method") or "maxlfq"
    run_normalization = _argv_value(argv, "run-normalization")
    sample_normalization = _argv_value(argv, "sample-normalization")
    explicit_normalization = any(
        method and method.casefold() != "none"
        for method in (run_normalization, sample_normalization)
    )
    defaults_to_normalization = quant_method not in {
        "ratio",
        "peptide-count",
        "spectral-count",
    }
    disabled = run_normalization == "none" and sample_normalization == "none"
    if explicit_normalization or (defaults_to_normalization and not disabled):
        steps.append("normalize")
    if "--batch-correction" in options or "--irs" in options:
        steps.append("correct")
    if "--impute-method" in options:
        steps.append("impute")
    if "--de-contrast" in options or "--de-contrast-file" in options:
        steps.append("differential")
    return [*steps, "export"]


def _argv_value(argv: list[str], name: str) -> str | None:
    option = f"--{name}"
    for index in range(len(argv) - 2, -1, -1):
        if argv[index] == option:
            return argv[index + 1]
        if argv[index].startswith(f"{option}="):
            return argv[index].split("=", 1)[1]
    return None


def _split_option(token: str) -> tuple[str, str | None]:
    if token.startswith("--") and "=" in token:
        option, value = token.split("=", 1)
        return option, value
    return token, None
