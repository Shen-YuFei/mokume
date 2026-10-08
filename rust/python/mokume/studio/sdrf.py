"""Bounded SDRF previews for Studio preflight review."""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any


_REFERENCE_PATTERN = re.compile(
    r"(?:^|[^a-z])(bridge|pool(?:ed)?|reference)(?:[^a-z]|$)", re.I
)
_EMPTY_VALUES = {"", "na", "n/a", "not applicable", "not available", "unknown"}


def read_sdrf(path: Path) -> dict[str, Any]:
    """Read a small mapped preview while counting every SDRF row."""
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        columns = _sdrf_columns(reader.fieldnames or [])
        rows = []
        issues: list[str] = []
        seen = set()
        row_count = 0
        reference_count = 0
        for row_count, source in enumerate(reader, start=1):
            mapped = _map_sdrf_row(source, columns)
            reference_count += int(mapped["reference"])
            key = (mapped["sample"], mapped["plex"], mapped["label"])
            if not mapped["sample"]:
                issues.append(f"row {row_count}: sample is missing")
            if key in seen and any(key):
                issues.append(f"row {row_count}: duplicate sample mapping")
            seen.add(key)
            if len(rows) < 24:
                rows.append(mapped)
    return {
        "path": str(path),
        "row_count": row_count,
        "headers": reader.fieldnames or [],
        "columns": columns,
        "rows": rows,
        "reference_count": reference_count,
        "issues": issues[:50],
        "truncated": row_count > len(rows),
    }


def _sdrf_columns(headers: list[str]) -> dict[str, str | None]:
    lowered = {header.casefold(): header for header in headers}

    def exact(*names: str) -> str | None:
        return next((lowered[name] for name in names if name in lowered), None)

    factor = next(
        (header for header in headers if header.casefold().startswith("factor value[")),
        None,
    )
    batch = next((header for header in headers if "batch" in header.casefold()), None)
    return {
        "sample": exact("source name", "assay name"),
        "condition": factor
        or exact("characteristics[disease]", "characteristics[organism part]"),
        "batch": batch,
        "biological_replicate": exact("characteristics[biological replicate]"),
        "technical_replicate": exact("comment[technical replicate]"),
        "plex": exact("comment[data file]", "assay name"),
        "label": exact("comment[label]"),
    }


def _clean_value(value: Any) -> str:
    text = "" if value is None else str(value).strip()
    return "" if text.casefold() in _EMPTY_VALUES else text


def _map_sdrf_row(
    source: dict[str, Any], columns: dict[str, str | None]
) -> dict[str, Any]:
    def value(name: str) -> str:
        column = columns.get(name)
        return _clean_value(source.get(column)) if column else ""

    values = [value(name) for name in columns]
    replicate = "/".join(
        filter(None, (value("biological_replicate"), value("technical_replicate")))
    )
    return {
        "sample": value("sample"),
        "condition": value("condition"),
        "batch": value("batch"),
        "replicate": replicate,
        "plex": value("plex"),
        "label": value("label"),
        "reference": any(_REFERENCE_PATTERN.search(item) for item in values),
    }
