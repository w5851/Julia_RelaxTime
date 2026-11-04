"""Store shared provenance once while retaining every figure-specific record."""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any

from scripts.plotting.plot_manifest import sha256_file


BUNDLE_SCHEMA = "plot_manifest_bundle_v1"


def build_bundle(index: dict[str, Any], records: list[dict[str, Any]]) -> dict[str, Any]:
    charts = index.get("charts", [])
    if not records or len(charts) != len(records):
        raise ValueError("a bundle requires one manifest record for each chart")
    common_keys = set.intersection(*(set(record) for record in records))
    shared = {key: copy.deepcopy(records[0][key]) for key in records[0]
              if key in common_keys and all(record[key] == records[0][key] for record in records)}
    figures = []
    for number, (chart, record) in enumerate(zip(charts, records)):
        figure_id = record.get("asset_id") or f"figure-{number + 1:03d}"
        summary = {key: copy.deepcopy(value) for key, value in chart.items()
                   if key not in {"manifest", "manifest_sha256", "manifest_record"}}
        figures.append({"figure_id": figure_id, "chart": summary,
                        "record": {key: copy.deepcopy(value) for key, value in record.items() if key not in shared}})
    if len({figure["figure_id"] for figure in figures}) != len(figures):
        raise ValueError("bundle figure IDs must be unique")
    return {"schema_version": BUNDLE_SCHEMA,
            "package": {key: copy.deepcopy(value) for key, value in index.items() if key != "charts"},
            "shared": shared, "figures": figures}


def expand_bundle(bundle: dict[str, Any]) -> list[tuple[dict[str, Any], dict[str, Any]]]:
    if bundle.get("schema_version") != BUNDLE_SCHEMA:
        raise ValueError("unsupported plot bundle schema")
    if not isinstance(bundle.get("package"), dict):
        raise ValueError("bundle requires package metadata")
    shared, figures = bundle.get("shared"), bundle.get("figures")
    if not isinstance(shared, dict) or not isinstance(figures, list) or not figures:
        raise ValueError("bundle requires shared provenance and a non-empty figures list")
    seen = set()
    expanded = []
    for entry in figures:
        identifier = entry.get("figure_id")
        if not isinstance(identifier, str) or not identifier or identifier in seen:
            raise ValueError("bundle figure IDs must be non-empty and unique")
        seen.add(identifier)
        if not isinstance(entry.get("record"), dict) or not isinstance(entry.get("chart"), dict):
            raise ValueError("each bundle figure requires chart and record objects")
        if set(shared) & set(entry["record"]):
            raise ValueError("a figure must not override shared provenance")
        expanded.append((copy.deepcopy(entry["chart"]), {**copy.deepcopy(shared), **copy.deepcopy(entry["record"])}))
    return expanded


def load_chart_records(path: Path, *, root: Path) -> tuple[dict, list[tuple[dict, dict]]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema_version") == BUNDLE_SCHEMA:
        return payload["package"], expand_bundle(payload)
    records = []
    for chart in payload["charts"]:
        source = root / chart["manifest"]
        if sha256_file(source) != chart["manifest_sha256"]:
            raise ValueError(f"chart manifest changed: {source}")
        records.append((chart, json.loads(source.read_text(encoding="utf-8"))))
    return payload, records
