import csv
import importlib.util
import json
from argparse import Namespace
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "collect_pnjl_issue130_crossover_mu_endpoint_pilot.py"
TARGET_LIST = ROOT / "docs" / "analysis" / "pnjl" / "issue130_endpoint_refinement_preflight_v1" / "crossover_mu" / "tables" / "target_list.csv"
CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"


def module():
    spec = importlib.util.spec_from_file_location("issue130_crossover_mu_endpoint_pilot", SCRIPT)
    assert spec is not None and spec.loader is not None
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


def write_target_artifact(
    root: Path,
    target_id: str,
    target_hash: str,
    *,
    found: bool = True,
    schema_version: str = "pnjl_issue130_crossover_mu_endpoint_pilot_v1",
    target_selection: str | None = None,
    workflow_sha: str = "f" * 40,
):
    target_dir = root / target_id
    target_dir.mkdir(parents=True)
    summary = {
        "schema_version": schema_version,
        "target_schema": "pnjl_issue130_endpoint_refinement_preflight_v1",
        "target_id": target_id,
        "calculation_sha": CALCULATION_SHA,
        "workflow_head_sha": workflow_sha,
        "target_list_sha256": target_hash,
        "target": {
            "xi": 0.0,
            "target_mu_MeV": 290.0,
            "mu_CEP_proxy_MeV": 300.0,
            "physical_side": "crossover_mu_lt_CEP_proxy",
        },
        "result": {
            "found": found,
            "status": "crossover_candidate" if found else "not_found",
            "T_crossover_MeV": 120.0 if found else None,
            "rho_crossover": 1.0 if found else None,
            "derivative_peak": 2.0 if found else None,
        },
        "finite_and_converged": found,
        "runner_seconds": 1.0,
        "solver": {"solver_called": True, "estimated_detector_calls": 24},
    }
    if target_selection is not None:
        summary["target_selection"] = target_selection
    (target_dir / "target_summary.json").write_text(json.dumps(summary), encoding="utf-8")
    (target_dir / "provenance.json").write_text(json.dumps({
        "calculation_sha": CALCULATION_SHA,
        "workflow_head_sha": workflow_sha,
        "reference_write": False,
    }), encoding="utf-8")
    with (target_dir / "temperature_response.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["T_MeV", "derivative", "rho", "finite", "solver_ok", "error", "scan_derivative"])
        writer.writeheader()
        writer.writerow({"T_MeV": "120.0", "derivative": "2.0", "rho": "1.0", "finite": "true", "solver_ok": "true", "error": "", "scan_derivative": "2.0"})


def test_aggregate_requires_exact_eight_targets(tmp_path):
    collector = module()
    expected, target_hash = collector.target_index(TARGET_LIST)
    input_dir = tmp_path / "inputs"
    for target_id in expected:
        write_target_artifact(input_dir, target_id, target_hash)
    output_dir = tmp_path / "aggregate"
    status = collector.aggregate(Namespace(
        input_dir=str(input_dir),
        output_dir=str(output_dir),
        target_list=str(TARGET_LIST),
        calculation_sha=CALCULATION_SHA,
        postprocess_sha="e" * 40,
        run_mode="numerical",
        source_run_id="",
    ))
    assert status == 0
    assert json.loads((output_dir / "verdict.json").read_text(encoding="utf-8"))["verdict"] == "pilot_candidate"
    assert len(list(csv.DictReader((output_dir / "pilot_summary.csv").open(encoding="utf-8")))) == 8


def test_aggregate_stops_on_missing_target(tmp_path):
    collector = module()
    expected, target_hash = collector.target_index(TARGET_LIST)
    input_dir = tmp_path / "inputs"
    for target_id in list(expected)[:-1]:
        write_target_artifact(input_dir, target_id, target_hash)
    output_dir = tmp_path / "aggregate"
    status = collector.aggregate(Namespace(
        input_dir=str(input_dir),
        output_dir=str(output_dir),
        target_list=str(TARGET_LIST),
        calculation_sha=CALCULATION_SHA,
        postprocess_sha="e" * 40,
        run_mode="numerical",
        source_run_id="",
    ))
    assert status == 2
    assert json.loads((output_dir / "verdict.json").read_text(encoding="utf-8"))["verdict"] == "pilot_artifact_invalid"


def test_full_expansion_selection_and_replay_provenance(tmp_path):
    collector = module()
    expected, target_hash = collector.target_index(TARGET_LIST, "full")
    assert len(expected) == 186
    input_dir = tmp_path / "inputs"
    for target_id in expected:
        write_target_artifact(
            input_dir,
            target_id,
            target_hash,
            schema_version="pnjl_issue130_crossover_mu_endpoint_expansion_v1",
            target_selection="full",
        )
    output_dir = tmp_path / "aggregate"
    status = collector.aggregate(Namespace(
        input_dir=str(input_dir),
        output_dir=str(output_dir),
        target_list=str(TARGET_LIST),
        calculation_sha=CALCULATION_SHA,
        postprocess_sha="e" * 40,
        source_workflow_sha="f" * 40,
        schema_version="pnjl_issue130_crossover_mu_endpoint_expansion_v1",
        selection="full",
        run_mode="aggregate_replay",
        source_run_id="12345",
    ))
    assert status == 0
    manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["target_selection"] == "full"
    assert manifest["solver_called"] is False
    assert json.loads((output_dir / "verdict.json").read_text(encoding="utf-8"))["verdict"] == "expansion_candidate"
    assert len(list(csv.DictReader((output_dir / "expansion_summary.csv").open(encoding="utf-8")))) == 186
