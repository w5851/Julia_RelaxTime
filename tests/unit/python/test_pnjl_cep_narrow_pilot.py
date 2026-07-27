import importlib.util
import json
from pathlib import Path

import pytest
import yaml


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "collect_pnjl_cep_narrow_pilot.py"
WORKFLOW = ROOT / ".github" / "workflows" / "pnjl-cep-narrow-pilot.yml"


def load_collector():
    spec = importlib.util.spec_from_file_location("cep_pilot_collector", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_workflow_has_immutable_matrix_and_aggregate_contract():
    text = WORKFLOW.read_text(encoding="utf-8")
    workflow = yaml.load(text, Loader=yaml.BaseLoader)
    inputs = workflow["on"]["workflow_dispatch"]["inputs"]
    assert set(inputs) == {"tag", "calculation_ref"}
    assert "^\u005b0-9a-fA-F\u005d{40}$" in text
    matrix = workflow["jobs"]["pilot"]["strategy"]["matrix"]
    assert len(matrix["xi"]) == 3
    assert len(matrix["method"]) == 3
    assert workflow["jobs"]["aggregate"]["needs"] == ["pilot"]
    assert "collect_pnjl_cep_narrow_pilot.py" in text
    assert "timeout-minutes: 330" in text


def test_collector_rejects_incomplete_matrix(tmp_path):
    module = load_collector()
    result = module.collect(tmp_path, tmp_path / "out")
    assert result["gate"]["status"] == "diagnostic_only"
    assert any("missing matrix jobs" in error for error in result["validation_errors"])


def test_collector_accepts_schema_smoke(tmp_path):
    module = load_collector()
    methods = sorted(module.METHODS)
    for xi in sorted(module.XIS):
        for method in methods:
            job_dir = tmp_path / f"job-{xi}-{method}"
            job_dir.mkdir()
            summary = {
                "xi": xi,
                "method": method,
                "calculation_sha": "a" * 40,
                "finite_and_converged": True,
                "targeted_additions": 0,
            }
            (job_dir / "job_summary.json").write_text(json.dumps(summary), encoding="utf-8")
            (job_dir / "slice_metrics.csv").write_text("xi,method,T_MeV\n0,method,130\n", encoding="utf-8")
            (job_dir / "method_costs.csv").write_text("xi,method,runner_seconds\n0,method,1\n", encoding="utf-8")
            (job_dir / "cep_accuracy.csv").write_text(
                "xi,method,canonical_T_CEP_MeV,canonical_muq_CEP_MeV,estimated_T_CEP_MeV,estimated_muq_CEP_MeV,delta_T_MeV,delta_muq_MeV,bracket_low_MeV,bracket_high_MeV,bracket_width_MeV,found,oracle_refine,oracle_refine_stable\n"
                f"{xi},{method},130,292,130,292,0,0,129.9,130.1,0.2,true,{{}},true\n",
                encoding="utf-8",
            )
    result = module.collect(tmp_path, tmp_path / "out")
    assert (tmp_path / "out" / "manifest.json").is_file()
    assert (tmp_path / "out" / "actions_costs.csv").is_file()
    assert result["gate"]["status"] == "pass"
