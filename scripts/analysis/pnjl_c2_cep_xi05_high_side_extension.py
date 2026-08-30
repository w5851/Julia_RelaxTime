"""Validate and plot the xi=0.5 high-side CEP extension artifact.

This collector is solver-free. It accepts only the three fixed equal-step
temperatures and leaves the physical bracket decision to author review.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
from typing import Any


_MANUAL_COLLECTOR_PATH = Path(__file__).with_name("pnjl_c2_cep_manual_bisection.py")
_MANUAL_COLLECTOR_SPEC = importlib.util.spec_from_file_location(
    "pnjl_c2_cep_manual_bisection", _MANUAL_COLLECTOR_PATH
)
if _MANUAL_COLLECTOR_SPEC is None or _MANUAL_COLLECTOR_SPEC.loader is None:
    raise ImportError(f"cannot load sibling collector: {_MANUAL_COLLECTOR_PATH}")
_MANUAL_COLLECTOR = importlib.util.module_from_spec(_MANUAL_COLLECTOR_SPEC)
_MANUAL_COLLECTOR_SPEC.loader.exec_module(_MANUAL_COLLECTOR)
plot_slice = _MANUAL_COLLECTOR.plot_slice
slug = _MANUAL_COLLECTOR.slug


CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
JOB_SCHEMA = "pnjl_c2_cep_xi05_high_side_extension_job_v1"
AUDIT_SCHEMA = "pnjl_c2_cep_xi05_high_side_extension_audit_v1"
XI = 0.5
ANCHOR_T = 107.0625
STEP = 0.0625
TEMPERATURES = (107.125, 107.1875, 107.25)
RHO_STEP = 0.003125
RHO_COUNT = 1281


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def finite(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _declared_file(manifest: dict[str, Any], root: Path, name: str) -> Path:
    path = root / name
    if not path.is_file():
        raise ValueError(f"missing {name}: {root}")
    if manifest.get("files", {}).get(name) != sha256(path):
        raise ValueError(f"hash mismatch for {name}: {root}")
    return path


def _validate_trace(trace: list[dict[str, str]]) -> None:
    if len(trace) != len(TEMPERATURES):
        raise ValueError("high-side trace must contain exactly three rows")
    trace.sort(key=lambda row: int(row["sequence"]))
    if [int(row["sequence"]) for row in trace] != [1, 2, 3]:
        raise ValueError("high-side trace sequence must be 1,2,3")
    temperatures = tuple(float(row["T_MeV"]) for row in trace)
    if temperatures != TEMPERATURES:
        raise ValueError(f"unexpected high-side temperatures: {temperatures}")
    for row in trace:
        if float(row["xi"]) != XI:
            raise ValueError("high-side trace xi mismatch")
        if row.get("direction") != "high":
            raise ValueError("high-side trace direction mismatch")
        if not row.get("role", "").startswith("high_extension_"):
            raise ValueError("high-side trace role is not versioned")
        if float(row["anchor_T_MeV"]) != ANCHOR_T:
            raise ValueError("high-side anchor mismatch")
        if float(row["step_MeV"]) != STEP:
            raise ValueError("high-side step mismatch")
        if row.get("oracle_labels_used_for_routing", "").lower() == "true":
            raise ValueError("oracle labels leaked into high-side routing")
    if any(abs(temperatures[i] - temperatures[i - 1] - STEP) > 1e-12
           for i in range(1, len(temperatures))):
        raise ValueError("high-side temperatures are not equally spaced")


def validate_artifact(root: Path, expected_source_run: str, expected_head: str) -> dict[str, Any]:
    manifest = read_json(root / "manifest.json")
    if manifest.get("schema_version") != JOB_SCHEMA:
        raise ValueError(f"unexpected job schema: {root}")
    if manifest.get("calculation_sha", "").lower() != CALCULATION_SHA:
        raise ValueError(f"calculation SHA mismatch: {root}")
    if expected_head and manifest.get("workflow_head_sha") != expected_head:
        raise ValueError(f"workflow head mismatch: {root}")
    if str(manifest.get("source_run_id")) != str(expected_source_run):
        raise ValueError(f"source run mismatch: {root}")
    if manifest.get("solver_called") is not True:
        raise ValueError("numerical artifact must record solver_called=true")
    if manifest.get("manual_decision_required") is not True:
        raise ValueError("manual decision contract missing")
    if manifest.get("extension_direction") != "high":
        raise ValueError("high-side extension direction missing")
    if float(manifest.get("extension_anchor_T_MeV")) != ANCHOR_T:
        raise ValueError("high-side extension anchor missing")
    if float(manifest.get("extension_step_MeV")) != STEP:
        raise ValueError("high-side extension step missing")
    fine_path = _declared_file(manifest, root, "fine_pool.csv")
    trace_path = _declared_file(manifest, root, "high_side_extension_trace.csv")
    cost_path = _declared_file(manifest, root, "method_costs.csv")
    fine = read_csv(fine_path)
    trace = read_csv(trace_path)
    costs = read_csv(cost_path)
    _validate_trace(trace)
    if len(fine) != len(TEMPERATURES) * RHO_COUNT:
        raise ValueError(f"unexpected fine-pool rows: {len(fine)}")
    keys = set()
    for row in fine:
        key = (round(float(row["xi"]), 8), round(float(row["T_MeV"]), 8),
               round(float(row["rho"]) / RHO_STEP))
        if key in keys:
            raise ValueError(f"duplicate fine-pool key: {key}")
        keys.add(key)
        if key[0] != XI or key[1] not in TEMPERATURES:
            raise ValueError(f"unexpected fine-pool coordinate: {key}")
        if not finite(row.get("muq_MeV")) or row.get("converged", "").lower() != "true" or row.get("finite", "").lower() != "true":
            raise ValueError("non-finite/non-converged fine-pool row")
        if row.get("calculation_sha", "").lower() != CALCULATION_SHA:
            raise ValueError("fine-pool calculation SHA mismatch")
    for row in costs:
        if int(row.get("failed_points", 0)) != 0:
            raise ValueError("failed solver points in cost table")
        if int(row["point_requests"]) != int(row["unique_solves"]) + int(row["cache_hits"]):
            raise ValueError("high-side cache cost does not reconcile")
    return {"manifest": manifest, "fine": fine, "trace": trace, "costs": costs}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--source-run-id", required=True)
    parser.add_argument("--expected-workflow-head", default="")
    args = parser.parse_args()

    artifact = validate_artifact(args.input_dir, args.source_run_id, args.expected_workflow_head)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figures = args.output_dir / "figures"
    audit_rows = []
    for row in artifact["trace"]:
        audit_rows.append({**row, **plot_slice(args.input_dir, row, figures)})
    fields = sorted({key for row in audit_rows for key in row})
    with (args.output_dir / "high_side_extension_audit.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows({key: row.get(key, "") for key in fields} for row in audit_rows)
    cost_fields = sorted({key for row in artifact["costs"] for key in row})
    with (args.output_dir / "method_costs.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=cost_fields)
        writer.writeheader()
        writer.writerows(artifact["costs"])
    manifest = {
        "schema_version": AUDIT_SCHEMA,
        "source_run_id": str(args.source_run_id),
        "calculation_sha": CALCULATION_SHA,
        "workflow_head_sha": args.expected_workflow_head,
        "xi": XI,
        "extension_direction": "high",
        "extension_anchor_T_MeV": ANCHOR_T,
        "extension_step_MeV": STEP,
        "temperatures": list(TEMPERATURES),
        "manual_decision_required": True,
        "oracle_labels_used_for_routing": False,
        "solver_called": False,
        "input_artifact_manifest_sha256": sha256(args.input_dir / "manifest.json"),
        "files": {},
        "plots": len(list(figures.glob("*.png"))),
    }
    for path in (args.output_dir / "high_side_extension_audit.csv", args.output_dir / "method_costs.csv"):
        manifest["files"][path.name] = sha256(path)
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    (args.output_dir / "claim_ledger.md").write_text(
        "# C2 xi=0.5 high-side extension audit v1\n\n"
        "- The three temperatures are fixed equal-step observations above the author-reviewed anchor.\n"
        "- The numerical runner did not use oracle labels or author decisions to select a slice.\n"
        "- This aggregate is solver-free and remains diagnostic-only until author review.\n"
        "- If all three slices remain first-order, the next action must extend the same 0.0625 MeV step.\n",
        encoding="utf-8",
    )
    (args.output_dir / "README.md").write_text(
        "# C2 xi=0.5 high-side extension audit v1\n\n"
        "Review `high_side_extension_audit.csv` and the local rho-mu figures. "
        "The three slices extend the open high-temperature side of the xi=0.5 bracket.\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
