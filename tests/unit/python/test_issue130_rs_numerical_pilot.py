from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import yaml


ROOT = Path(__file__).parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "relaxtime" / "collect_issue130_rs_numerical_pilot.py"
WORKFLOW = ROOT / ".github" / "workflows" / "relaxtime-issue130-rs-numerical-pilot-v1.yml"
CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
MODE = "mode_a_fixed_muB_phase_scaled"
FIELDS = [
    "muB_MeV",
    "xi",
    "mode",
    "alpha_T",
    "T_MeV",
    "converged",
    "quality_flag",
    "quality_reason",
    "phase_reference_kind",
    "phase_structure",
    "eta_over_s",
    "sigma_over_T",
    "zeta_over_s",
    "tau_u",
    "tau_d",
    "tau_s",
    "tau_ubar",
    "tau_dbar",
    "tau_sbar",
]
IDENTITY = {"muB_MeV": "450.0", "xi": "0.0", "mode": MODE, "alpha_T": "1.0"}


def _write_csv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _row(*, temperature: str = "180.0", sigma: str = "0.4", xi: str | None = None) -> dict[str, str]:
    row = {
        **IDENTITY,
        "T_MeV": temperature,
        "converged": "true",
        "quality_flag": "false",
        "quality_reason": "ok",
        "phase_reference_kind": "crossover",
        "phase_structure": "crossover_possible",
        "eta_over_s": "0.2",
        "sigma_over_T": sigma,
        "zeta_over_s": "0.03",
        "tau_u": "1.0",
        "tau_d": "1.1",
        "tau_s": "1.2",
        "tau_ubar": "1.3",
        "tau_dbar": "1.4",
        "tau_sbar": "1.5",
    }
    if xi is not None:
        row["xi"] = xi
    return row


def _prepare_pair(tmp_path: Path, *, candidate_rows: list[dict[str, str]] | None = None, legacy_rows: list[dict[str, str]] | None = None) -> Path:
    root = tmp_path / "pilot"
    for reference, rows in (("candidate_runtime", candidate_rows), ("legacy", legacy_rows)):
        rows = rows if rows is not None else [_row(temperature="180.0" if reference == "candidate_runtime" else "179.9")]
        result_dir = root / "results" / reference
        _write_csv(result_dir / "phase_guided_transport_scan.csv", FIELDS, rows)
        plan_fields = ["T_MeV", "muB_MeV", "xi", "mode", "alpha_T"]
        _write_csv(result_dir / "sampling_plan.csv", plan_fields, [{field: row[field] for field in plan_fields} for row in rows])
        _write_csv(result_dir / "failed_points.csv", ["T_MeV", "muB_MeV", "xi"], [])
        (result_dir / "effective_config.json").write_text(
            json.dumps({"options": {"phase_reference_mode": "runtime" if reference == "candidate_runtime" else "legacy", "phase_reference_source": "candidate" if reference == "candidate_runtime" else "legacy"}}),
            encoding="utf-8",
        )
        (result_dir / "run_manifest.json").write_text(json.dumps({"run_id": f"{reference}-run"}), encoding="utf-8")
        (result_dir / "pilot_metadata.json").write_text(
            json.dumps({"reference_mode": reference, "calculation_sha": CALCULATION_SHA, "solver_called": True}),
            encoding="utf-8",
        )
    (root / "status").mkdir(parents=True, exist_ok=True)
    (root / "status" / "candidate_runtime.txt").write_text("0\n", encoding="utf-8")
    (root / "status" / "legacy.txt").write_text("0\n", encoding="utf-8")
    return root


def _run(root: Path, output: Path, *, expected_points: int = 1, fail: bool = False) -> subprocess.CompletedProcess[str]:
    cmd = [
        sys.executable,
        str(SCRIPT),
        "--pilot-root",
        str(root),
        "--output",
        str(output),
        "--repo-sha",
        "workflow-sha",
        "--calculation-sha",
        CALCULATION_SHA,
        "--expected-points",
        str(expected_points),
    ]
    if fail:
        cmd.append("--fail-on-hard-failure")
    return subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, check=False)


def test_paired_success_pairs_by_request_not_anchor_temperature(tmp_path: Path) -> None:
    root = _prepare_pair(tmp_path)
    output = tmp_path / "aggregate"
    result = _run(root, output, fail=True)
    assert result.returncode == 0, result.stderr
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["verdict"] == "pilot_pair_complete_diagnostic_only"
    assert manifest["common_row_count"] == 1
    assert manifest["reference_modes"]["candidate_runtime"]["result_path"] == "results/candidate_runtime/phase_guided_transport_scan.csv"
    comparison = (output / "transport_comparison.csv").read_text(encoding="utf-8")
    assert "T_MeV_candidate" in comparison
    assert "T_MeV_legacy" in comparison


def test_common_tau_ratio_warning_is_diagnostic_not_solver_failure(tmp_path: Path) -> None:
    candidate = _row()
    legacy = _row(temperature="179.9")
    for row in (candidate, legacy):
        row.update(
            {
                "quality_flag": "true",
                "quality_reason": "tau_u_ubar_ratio_high",
                "tau_u": "7.0",
                "tau_ubar": "1.0",
            }
        )
    root = _prepare_pair(tmp_path, candidate_rows=[candidate], legacy_rows=[legacy])
    result = _run(root, tmp_path / "aggregate", fail=True)
    assert result.returncode == 0, result.stderr
    manifest = json.loads((tmp_path / "aggregate" / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["verdict"] == "pilot_pair_complete_with_common_quality_warnings_diagnostic_only"
    assert manifest["common_quality_warning_count"] == 1


def test_quality_warning_only_on_one_reference_is_hard_failure(tmp_path: Path) -> None:
    candidate = _row()
    candidate.update({"quality_flag": "true", "quality_reason": "tau_u_ubar_ratio_high"})
    root = _prepare_pair(tmp_path, candidate_rows=[candidate])
    result = _run(root, tmp_path / "aggregate", fail=True)
    assert result.returncode == 2
    verdict = json.loads((tmp_path / "aggregate" / "verdict.json").read_text(encoding="utf-8"))
    assert "candidate_legacy_quality_warning_key_mismatch" in verdict["hard_failures"]


def test_unknown_quality_reason_remains_hard_failure(tmp_path: Path) -> None:
    candidate = _row()
    legacy = _row(temperature="179.9")
    for row in (candidate, legacy):
        row.update({"quality_flag": "true", "quality_reason": "unexpected_quality_reason"})
    root = _prepare_pair(tmp_path, candidate_rows=[candidate], legacy_rows=[legacy])
    result = _run(root, tmp_path / "aggregate", fail=True)
    assert result.returncode == 2
    verdict = json.loads((tmp_path / "aggregate" / "verdict.json").read_text(encoding="utf-8"))
    assert "candidate_runtime:invalid_quality_hard_1" in verdict["hard_failures"]
    assert "legacy:invalid_quality_hard_1" in verdict["hard_failures"]


def test_missing_failed_point_is_hard_failure(tmp_path: Path) -> None:
    root = _prepare_pair(tmp_path)
    _write_csv(root / "results" / "legacy" / "failed_points.csv", ["T_MeV", "muB_MeV", "xi"], [{"T_MeV": "1", "muB_MeV": "2", "xi": "0"}])
    result = _run(root, tmp_path / "aggregate", fail=True)
    assert result.returncode == 2
    verdict = json.loads((tmp_path / "aggregate" / "verdict.json").read_text(encoding="utf-8"))
    assert "legacy:failed_points_1" in verdict["hard_failures"]


def test_nonfinite_transport_is_hard_failure(tmp_path: Path) -> None:
    root = _prepare_pair(tmp_path, candidate_rows=[_row(sigma="nan")])
    result = _run(root, tmp_path / "aggregate", fail=True)
    assert result.returncode == 2
    verdict = json.loads((tmp_path / "aggregate" / "verdict.json").read_text(encoding="utf-8"))
    assert "candidate_runtime:nonfinite_1" in verdict["hard_failures"]


def test_duplicate_request_key_is_hard_failure(tmp_path: Path) -> None:
    rows = [_row(temperature="180.0"), _row(temperature="181.0")]
    root = _prepare_pair(tmp_path, candidate_rows=rows, legacy_rows=rows)
    result = _run(root, tmp_path / "aggregate", expected_points=2, fail=True)
    assert result.returncode == 2
    verdict = json.loads((tmp_path / "aggregate" / "verdict.json").read_text(encoding="utf-8"))
    assert "candidate_runtime:duplicate_keys_1" in verdict["hard_failures"]


def test_no_common_request_keys_is_hard_failure(tmp_path: Path) -> None:
    root = _prepare_pair(tmp_path, candidate_rows=[_row()], legacy_rows=[_row(xi="0.5")])
    result = _run(root, tmp_path / "aggregate", fail=True)
    assert result.returncode == 2
    verdict = json.loads((tmp_path / "aggregate" / "verdict.json").read_text(encoding="utf-8"))
    assert "candidate_legacy_plan_key_mismatch" in verdict["hard_failures"]
    assert "no_common_transport_rows" in verdict["hard_failures"]


def test_workflow_freezes_paired_scope_and_solver_provenance() -> None:
    text = WORKFLOW.read_text(encoding="utf-8")
    payload = yaml.load(text, Loader=yaml.BaseLoader)
    assert payload["jobs"]["numerical_pilot"]["timeout-minutes"] == "720"
    assert payload["env"]["CALCULATION_SHA"] == CALCULATION_SHA
    assert payload["env"]["EXPECTED_POINTS"] == "19"
    assert "tag contains unsupported characters" in text
    assert "run_reference candidate_runtime runtime" in text
    assert "run_reference legacy legacy" in text
    assert "--phase-reference-mode \"$phase_mode\"" in text
    assert "legacy_phase_reference_v1/RETIREMENT_MANIFEST.json" in text
    assert "--fail-on-hard-failure" in text
    assert '"solver_called": true' in text
    assert '"production_write": false' in text
