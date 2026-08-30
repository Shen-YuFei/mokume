"""Machine-readable native command catalog and safe argv canonicalization."""

from __future__ import annotations

import importlib
from collections.abc import Iterable
from pathlib import Path
from typing import Any

from mokume.studio.paths import ProjectPaths


MVP_COMMAND_PATHS = {
    ("quantify", "features2proteins"),
    ("quantify", "features2peptides"),
    ("quantify", "peptides2protein"),
    ("correct-batches",),
}

OUTPUT_ARGUMENTS = {
    "de-output",
    "export-ions",
    "export-peptides",
    "generate-filter-config",
    "log-file",
    "output",
    "qc-report",
}


class CommandValidationError(ValueError):
    """Raised when a Studio command does not match the installed CLI contract."""


def command_schema() -> list[dict[str, Any]]:
    """Return the installed native schema, restricted to Studio MVP commands."""
    package = importlib.import_module("mokume")
    schema = getattr(package, "command_schema")()
    commands = schema.get("commands", schema) if isinstance(schema, dict) else schema
    return [
        item for item in commands if tuple(item.get("path", ())) in MVP_COMMAND_PATHS
    ]


def validate_and_canonicalize(
    argv: Iterable[str], project_root: str | Path
) -> list[str]:
    """Validate native argv and turn every path argument into a guarded absolute path."""
    original = list(argv)
    spec = _command_spec(original)
    canonical = _canonicalize_paths(original, spec, ProjectPaths(project_root))
    package = importlib.import_module("mokume")
    getattr(package, "validate_args")(canonical)
    return canonical


def command_paths(argv: Iterable[str]) -> tuple[list[Path], list[Path]]:
    """Return explicit input and output paths from already-canonical argv."""
    tokens = list(argv)
    spec = _command_spec(tokens)
    options = _option_index(spec.get("flags", ()))
    inputs: list[Path] = []
    outputs: list[Path] = []
    index = len(spec["path"])
    while index < len(tokens):
        option, inline_value = _split_option(tokens[index])
        argument = options.get(option)
        if argument is None:
            index += 1
            continue
        count = _value_count(argument)
        values = ([inline_value] if inline_value is not None else []) + tokens[
            index + 1 : index + 1 + count - (inline_value is not None)
        ]
        target = outputs if str(argument.get("long")) in OUTPUT_ARGUMENTS else inputs
        target.extend(_path_values(argument, values))
        index += 1 + count - (inline_value is not None)
    return inputs, outputs


def _path_values(argument: dict[str, Any], values: list[str]) -> list[Path]:
    names = argument.get("value_names") or []
    return [
        Path(value)
        for position, value in enumerate(values)
        if names and names[min(position, len(names) - 1)] in {"FILE", "DIR"}
    ]


def _command_spec(argv: list[str]) -> dict[str, Any]:
    for spec in command_schema():
        path = list(spec["path"])
        if argv[: len(path)] == path:
            return spec
    raise CommandValidationError("command is not available in Mokume Studio")


def _canonicalize_paths(
    argv: list[str], spec: dict[str, Any], paths: ProjectPaths
) -> list[str]:
    command_length = len(spec["path"])
    options = _option_index(spec.get("flags", ()))
    canonical = argv[:command_length]
    index = command_length
    while index < len(argv):
        token = argv[index]
        option, inline_value = _split_option(token)
        argument = options.get(option)
        if argument is None:
            canonical.append(token)
            index += 1
            continue
        canonical.append(option)
        value_count = _value_count(argument)
        values = ([inline_value] if inline_value is not None else []) + argv[
            index + 1 : index + 1 + value_count - (inline_value is not None)
        ]
        if len(values) != value_count:
            raise CommandValidationError(f"missing value for {option}")
        canonical.extend(_canonical_values(argument, option, values, paths))
        index += 1 + value_count - (inline_value is not None)
    return canonical


def _option_index(arguments: Iterable[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    options: dict[str, dict[str, Any]] = {}
    for argument in arguments:
        if argument.get("long"):
            options[f"--{argument['long']}"] = argument
        if argument.get("short"):
            options[f"-{argument['short']}"] = argument
    return options


def _split_option(token: str) -> tuple[str, str | None]:
    if token.startswith("--") and "=" in token:
        return tuple(token.split("=", 1))  # type: ignore[return-value]
    return token, None


def _value_count(argument: dict[str, Any]) -> int:
    arity = argument.get("value_arity") or {}
    maximum = arity.get("max")
    minimum = int(arity.get("min", 0))
    if maximum == 0 or not argument.get("value_names"):
        return 0
    if maximum is None or maximum != minimum:
        return max(1, int(minimum))
    return int(maximum)


def _canonical_values(
    argument: dict[str, Any],
    option: str,
    values: list[str],
    paths: ProjectPaths,
) -> list[str]:
    names = argument.get("value_names") or []
    canonical: list[str] = []
    for position, value in enumerate(values):
        value_name = names[min(position, len(names) - 1)] if names else ""
        if value_name not in {"FILE", "DIR"}:
            canonical.append(value)
            continue
        long_name = str(argument.get("long") or option.lstrip("-"))
        resolved = (
            paths.resolve_output(value)
            if long_name in OUTPUT_ARGUMENTS
            else paths.resolve_existing(value)
        )
        canonical.append(str(resolved))
    return canonical
