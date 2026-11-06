import csv
import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "pnjl_c2_cep_xi05_high_side_extension.py"


def load_module():
    spec = importlib.util.spec_from_file_location("c2_xi05_high_side", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_minimal_artifact(module, root: Path, *, routing_leak: bool = False):
    root.mkdir()
    rows = []
    for temperature in module.TEMPERATURES:
        for index in range(module.RHO_COUNT):
            rows.append({
                "xi": "0.5", "T_MeV": str(temperature), "rho": str(index * module.RHO_STEP),
                "muq_MeV": "300.0", "residual_norm": "0.0", "iterations": "1",
                "converged": "true", "finite": "true", "calculation_sha": module.CALCULATION_SHA,
            })
    with (root / "fine_pool.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    trace_rows = [{
        "xi": "0.5", "sequence": str(index), "role": f"high_extension_{index}",
        "direction": "high", "T_MeV": str(temperature), "anchor_T_MeV": "107.0625",
        "step_MeV": "0.0625", "oracle_labels_used_for_routing": str(routing_leak).lower(),
    } for index, temperature in enumerate(module.TEMPERATURES, start=1)]
    with (root / "high_side_extension_trace.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(trace_rows[0]))
        writer.writeheader()
        writer.writerows(trace_rows)
    costs = [{"point_requests": "1", "unique_solves": "1", "cache_hits": "0", "failed_points": "0"}]
    with (root / "method_costs.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(costs[0]))
        writer.writeheader()
        writer.writerows(costs)
    manifest = {
        "schema_version": module.JOB_SCHEMA, "calculation_sha": module.CALCULATION_SHA,
        "workflow_head_sha": "a" * 40, "source_run_id": "123", "solver_called": True,
        "manual_decision_required": True, "extension_direction": "high",
        "extension_anchor_T_MeV": module.ANCHOR_T, "extension_step_MeV": module.STEP,
        "files": {name: module.sha256(root / name) for name in (
            "fine_pool.csv", "high_side_extension_trace.csv", "method_costs.csv")},
    }
    (root / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")


def test_validator_rejects_oracle_routing_leak(tmp_path):
    module = load_module()
    root = tmp_path / "artifact"
    _write_minimal_artifact(module, root, routing_leak=True)
    with pytest.raises(ValueError, match="oracle labels leaked"):
        module.validate_artifact(root, "123", "a" * 40)
