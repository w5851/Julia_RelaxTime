#!/usr/bin/env python3
"""Freeze immutable validation windows from the three v2 cascade jobs."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

XIS = {-0.5, 0.0, 0.5}


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _find_proposals(input_dir: Path) -> list[dict[str, Any]]:
    proposals = []
    for path in sorted(input_dir.rglob("window_proposal.json")):
        payload = _json(path)
        if payload.get("schema_version") != "cep_narrow_pilot_v2_window_proposal":
            raise ValueError(f"unexpected proposal schema: {path}")
        payload["source_path"] = str(path)
        proposals.append(payload)
    return proposals


def freeze(input_dir: Path, output_path: Path) -> dict[str, Any]:
    proposals = _find_proposals(input_dir)
    if len(proposals) != 3:
        raise ValueError(f"expected exactly 3 cascade proposals, found {len(proposals)}")
    seen = {float(row["xi"]) for row in proposals}
    if seen != XIS:
        raise ValueError(f"cascade proposals must cover {sorted(XIS)}, found {sorted(seen)}")
    shas = {str(row.get("calculation_sha", "")) for row in proposals}
    if len(shas) != 1 or not next(iter(shas)):
        raise ValueError(f"cascade proposals must share one non-empty calculation SHA: {sorted(shas)}")

    windows = []
    for row in sorted(proposals, key=lambda value: float(value["xi"])):
        anchors = sorted({float(value) for value in row.get("required_T_anchors", [])})
        windows.append(
            {
                "xi": float(row["xi"]),
                "T_min_MeV": float(row["T_min_MeV"]),
                "T_max_MeV": float(row["T_max_MeV"]),
                "required_T_anchors": anchors,
                "discovery_result_status": row.get("cep", {}).get("result_status", "ambiguous"),
                "discovery_T_last_first_order_MeV": row.get("cep", {}).get("T_last_first_order_MeV"),
                "discovery_T_first_monotone_MeV": row.get("cep", {}).get("T_first_monotone_MeV"),
            }
        )

    payload: dict[str, Any] = {
        "schema_version": "cep_narrow_pilot_v2_validation_windows",
        "generated_at_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "calculation_sha": next(iter(shas)),
        "label_policy": "bounds_and_common_anchors_only; discovery labels are diagnostics",
        "windows": windows,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    payload["sha256"] = hashlib.sha256(output_path.read_bytes()).hexdigest()
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = freeze(args.input_dir, args.output)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
