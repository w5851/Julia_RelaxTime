import csv
import hashlib
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
COLLECTOR = ROOT / "scripts" / "analysis" / "collect_pnjl_cep_hybrid_production_shadow.py"
WORKFLOW = ROOT / ".github" / "workflows" / "pnjl-cep-hybrid-production-shadow.yml"
ENDPOINT_LOCAL_WORKFLOW = ROOT / ".github" / "workflows" / "pnjl-maxwell-endpoint-local-production-shadow-v4.yml"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _status(xi: float, t: float):
    first_order = {
        -0.5: {5.0, 20.0, 60.0, 100.0, 130.0, 147.0947265625},
        0.0: {5.0, 20.0, 60.0, 100.0, 120.0, 130.9619140625},
        0.5: {5.0, 20.0, 60.0, 90.0, 100.0, 106.9599609375},
    }
    monotone = {
        -0.5: {147.2197265625, 160.0},
        0.0: {131.0869140625, 145.0},
        0.5: {107.0849609375, 120.0},
    }
    if t in first_order[xi]:
        return "confirmed_first_order"
    if t in monotone[xi]:
        return "confirmed_monotone"
    return "ambiguous_near_critical"


def write_job(root: Path, xi: float, method: str, *, cost: int = 100, scope: str = "full"):
    collector = load_module(COLLECTOR, f"hybrid_collector_constants_{xi}_{method}")
    job = root / f"job-{method}-{xi}"
    job.mkdir()
    anchors = (
        collector.FOCUSED_ANCHORS[xi]
        if scope == "focused"
        else collector.TARGETED_ANCHORS[xi]
        if scope == "targeted"
        else collector.ANCHORS[xi]
    )
    sha = "a" * 40
    with (job / "curve_points.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["xi", "method", "T_MeV", "rho_level", "rho", "muq_MeV", "converged", "finite"])
        writer.writeheader()
        for index, t in enumerate(anchors):
            writer.writerow({"xi": xi, "method": method, "T_MeV": t, "rho_level": 0, "rho": 1.0, "muq_MeV": 300.0 + index, "converged": "true", "finite": "true"})
    with (job / "slice_metrics.csv").open("w", newline="", encoding="utf-8") as handle:
        fields = ["xi", "method", "T_MeV", "stage_a_status", "stage_b_status", "stage_c_status", "stage_used", "upgrade_reason", "result_status", "geometry_converged", "position_error_MeV", "density_error", "maxwell_area_gate", "area_residual", "mu_transition_MeV", "rho_hadron", "rho_quark", "mu_spinodal_hadron_MeV", "mu_spinodal_quark_MeV", "rho_spinodal_hadron", "rho_spinodal_quark", "support_low", "support_high", "support_point_count", "point_ranking_version", "stage_c_stop_reason", "stage_c_actual_cap", "stage_c_selected_points_json", "stage_c_component_geometry_json", "stage_c_refinement_trace_json", "targeted_additions", "solver_failure_count"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for t in anchors:
            status = _status(xi, t)
            first = status == "confirmed_first_order"
            writer.writerow({
                "xi": xi, "method": method, "T_MeV": t,
                "stage_a_status": status, "stage_b_status": "not_run", "stage_c_status": "not_run",
                "stage_used": "stage_a", "upgrade_reason": "synthetic", "result_status": status,
                "geometry_converged": "true", "position_error_MeV": 0.01, "density_error": 0.001,
                "maxwell_area_gate": 0.00001, "area_residual": 0.00001,
                "mu_transition_MeV": 300.0 if first else "", "rho_hadron": 1.0 if first else "",
                "rho_quark": 2.0 if first else "", "mu_spinodal_hadron_MeV": 301.0 if first else "",
                "mu_spinodal_quark_MeV": 299.0 if first else "", "rho_spinodal_hadron": 0.9 if first else "",
                "rho_spinodal_quark": 2.1 if first else "", "support_low": "", "support_high": "",
                "support_point_count": 0,
                "point_ranking_version": "stage_b_features_v1" if method == "production_hybrid" else "",
                "stage_c_stop_reason": "geometry_certificate_closed" if method == "production_hybrid" else "",
                "stage_c_actual_cap": 9 if method == "production_hybrid" else 0,
                "stage_c_selected_points_json": '[{"feature":"stage_b_features_v1","rank":1,"batch":1,"rho":0.1}]' if method == "production_hybrid" else "[]",
                "stage_c_component_geometry_json": '[{"component":"rho_hadron","status":"pass","pass":true}]' if method == "production_hybrid" else "[]",
                "stage_c_refinement_trace_json": "[]",
                "targeted_additions": 9 if method == "production_hybrid" else 0,
                "solver_failure_count": 0,
            })
    (job / "method_costs.csv").write_text(
        "xi,method,unique_solves,equilibrium_requests,fixedrho_requests,residual_calls,jacobian_calls,newton_iterations,runner_seconds,fallback_count,retry_count\n"
        f"{xi},{method},{cost},{cost},{cost},100,50,120,10,0,0\n", encoding="utf-8")
    (job / "cep_accuracy.csv").write_text(
        "xi,method,anchor_T_MeV,result_status,T_last_first_order_MeV,muq_last_first_order_MeV,T_first_monotone_MeV,ambiguity_width_T_MeV,temperature_resolution_target_MeV\n"
        + "".join(f"{xi},{method},{t},{_status(xi, t)},20,300,120,100,0.125\n" for t in anchors), encoding="utf-8")
    hashes = {name: hashlib.sha256((job / name).read_bytes()).hexdigest() for name in ("curve_points.csv", "slice_metrics.csv", "method_costs.csv", "cep_accuracy.csv")}
    (job / "job_summary.json").write_text(json.dumps({
        "schema_version": "cep_cascade_production_shadow_v2", "xi": xi, "method": method,
        "calculation_sha": sha, "workflow_head_sha": sha, "anchors": list(anchors),
        "scope": scope,
        "finite_and_converged_final": True, "curve_file_sha256": hashes,
        "provenance": {"calculation_sha": sha, "reference_write": False, "anchor_run_count": len(anchors)},
    }), encoding="utf-8")


def _matrix(tmp_path: Path):
    for xi in (-0.5, 0.0, 0.5):
        write_job(tmp_path, xi, "production_hybrid", cost=90)
        write_job(tmp_path, xi, "memoized_dense", cost=100)
        write_job(tmp_path, xi, "independent_oracle", cost=110)


def _targeted_matrix(tmp_path: Path):
    for xi in (-0.5, 0.0, 0.5):
        write_job(tmp_path, xi, "production_hybrid", cost=90, scope="targeted")
        write_job(tmp_path, xi, "memoized_dense", cost=100, scope="targeted")
        write_job(tmp_path, xi, "independent_oracle", cost=110, scope="targeted")


def _focused_matrix(tmp_path: Path):
    for xi in (-0.5, 0.0, 0.5):
        write_job(tmp_path, xi, "production_hybrid", cost=90, scope="focused")
        write_job(tmp_path, xi, "memoized_dense", cost=100, scope="focused")
        write_job(tmp_path, xi, "independent_oracle", cost=110, scope="focused")


def test_hybrid_collector_accepts_complete_matrix(tmp_path):
    module = load_module(COLLECTOR, "hybrid_collector_complete")
    _matrix(tmp_path)
    gate = module.collect(tmp_path, tmp_path / "aggregate", None, "a" * 40)
    assert gate["verdict"] == "full_hybrid_candidate"
    assert (tmp_path / "aggregate" / "geometry_accuracy.csv").is_file()
    assert (tmp_path / "aggregate" / "manifest.json").is_file()
    rows = list(csv.DictReader((tmp_path / "aggregate" / "geometry_accuracy.csv").open(newline="", encoding="utf-8")))
    production_row = next(row for row in rows if row["method"] == "production_hybrid")
    assert production_row["point_ranking_version"] == "stage_b_features_v1"
    assert production_row["stage_c_stop_reason"] == "geometry_certificate_closed"
    assert production_row["stage_c_actual_cap"] == "9"


def test_hybrid_collector_accepts_targeted_matrix_and_uses_targeted_verdict(tmp_path):
    module = load_module(COLLECTOR, "hybrid_collector_targeted")
    _targeted_matrix(tmp_path)
    gate = module.collect(tmp_path, tmp_path / "aggregate", None, "a" * 40, scope="targeted")
    assert gate["verdict"] == "targeted_hybrid_candidate"
    assert gate["scope"] == "targeted"
    assert gate["expected_anchor_count"] == 18


def test_hybrid_collector_rejects_confirmation_on_oracle_ambiguous(tmp_path):
    module = load_module(COLLECTOR, "hybrid_collector_unsupported")
    _matrix(tmp_path)
    oracle = tmp_path / "job-independent_oracle--0.5" / "slice_metrics.csv"
    text = oracle.read_text(encoding="utf-8").replace(
        "stage_a,synthetic,confirmed_first_order,true",
        "stage_a,synthetic,ambiguous_near_critical,true",
        1,
    )
    oracle.write_text(text, encoding="utf-8")
    summary_path = oracle.parent / "job_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["curve_file_sha256"]["slice_metrics.csv"] = hashlib.sha256(oracle.read_bytes()).hexdigest()
    summary_path.write_text(json.dumps(summary), encoding="utf-8")
    gate = module.collect(tmp_path, tmp_path / "aggregate", None, "a" * 40)
    assert gate["verdict"] == "oracle_inconclusive"
    assert gate["classification_errors"]


def test_replay_repairs_legacy_mu_transition_field_and_records_correction(tmp_path):
    module = load_module(COLLECTOR, "hybrid_collector_repair")
    _matrix(tmp_path)
    slice_path = tmp_path / "job-production_hybrid--0.5" / "slice_metrics.csv"
    rows = list(csv.DictReader(slice_path.open(newline="", encoding="utf-8")))
    rows[0]["mu_transition_MeV"] = ""
    with slice_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    summary_path = slice_path.parent / "job_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["curve_file_sha256"]["slice_metrics.csv"] = hashlib.sha256(slice_path.read_bytes()).hexdigest()
    summary_path.write_text(json.dumps(summary), encoding="utf-8")

    output = tmp_path / "aggregate"
    gate = module.collect(tmp_path, output, None, "a" * 40)
    assert gate["verdict"] == "full_hybrid_candidate"
    repaired = list(csv.DictReader((output / "slice_metrics.csv").open(newline="", encoding="utf-8")))
    row = next(row for row in repaired if row["method"] == "production_hybrid" and row["xi"] == "-0.5" and row["T_MeV"] == "5.0")
    assert float(row["mu_transition_MeV"]) == 300.0
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["aggregation_corrections"]


def test_required_oracle_anchor_is_explicitly_gated(tmp_path):
    module = load_module(COLLECTOR, "hybrid_collector_required_anchor")
    _matrix(tmp_path)
    oracle_path = tmp_path / "job-independent_oracle--0.5" / "slice_metrics.csv"
    rows = list(csv.DictReader(oracle_path.open(newline="", encoding="utf-8")))
    next(row for row in rows if row["T_MeV"] == "60.0")["result_status"] = "ambiguous_near_critical"
    with oracle_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    summary_path = oracle_path.parent / "job_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["curve_file_sha256"]["slice_metrics.csv"] = hashlib.sha256(oracle_path.read_bytes()).hexdigest()
    summary_path.write_text(json.dumps(summary), encoding="utf-8")
    gate = module.collect(tmp_path, tmp_path / "aggregate", None, "a" * 40)
    assert gate["verdict"] == "oracle_inconclusive"
    assert any("required first-order anchor" in error for error in gate["oracle_errors"])


def test_legacy_replay_keeps_stage_c_support_mismatch_as_warning(tmp_path):
    module = load_module(COLLECTOR, "hybrid_collector_legacy_replay")
    _matrix(tmp_path)
    job = tmp_path / "job-production_hybrid--0.5"

    curve_path = job / "curve_points.csv"
    curve_rows = list(csv.DictReader(curve_path.open(newline="", encoding="utf-8")))
    for row in curve_rows:
        row["sampling_role"] = "stage_c_support"
    with curve_path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = list(curve_rows[0])
        "sampling_role" not in fieldnames and fieldnames.append("sampling_role")
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(curve_rows)

    slice_path = job / "slice_metrics.csv"
    slice_rows = list(csv.DictReader(slice_path.open(newline="", encoding="utf-8")))
    slice_rows[0]["stage_c_status"] = "invalid"
    slice_rows[0]["support_low"] = "0.0"
    slice_rows[0]["support_high"] = "0.5"
    with slice_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(slice_rows[0]))
        writer.writeheader()
        writer.writerows(slice_rows)

    summary_path = job / "job_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["curve_file_sha256"]["curve_points.csv"] = hashlib.sha256(curve_path.read_bytes()).hexdigest()
    summary["curve_file_sha256"]["slice_metrics.csv"] = hashlib.sha256(slice_path.read_bytes()).hexdigest()
    summary_path.write_text(json.dumps(summary), encoding="utf-8")

    strict = module.collect(tmp_path, tmp_path / "strict", None, "a" * 40)
    assert strict["verdict"] == "workflow_failure"
    replay = module.collect(tmp_path, tmp_path / "replay", None, "a" * 40, legacy_replay=True)
    assert replay["verdict"] == "full_hybrid_candidate"
    assert replay["compatibility_warnings"]
    manifest = json.loads((tmp_path / "replay" / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["compatibility_warnings"]


def test_stage_c_guard_matches_six_decimal_curve_temperature(tmp_path):
    module = load_module(COLLECTOR, "hybrid_collector_curve_temperature_rounding")
    _matrix(tmp_path)
    job = tmp_path / "job-production_hybrid--0.5"
    curve_path = job / "curve_points.csv"
    curve_rows = list(csv.DictReader(curve_path.open(newline="", encoding="utf-8")))
    for row in curve_rows:
        if row["T_MeV"] == "147.0947265625":
            row["T_MeV"] = "147.094727"
            row["sampling_role"] = "stage_c_guard"
    with curve_path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = list(curve_rows[0])
        "sampling_role" not in fieldnames and fieldnames.append("sampling_role")
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(curve_rows)
    slice_path = job / "slice_metrics.csv"
    slice_rows = list(csv.DictReader(slice_path.open(newline="", encoding="utf-8")))
    target = next(row for row in slice_rows if row["T_MeV"] == "147.0947265625")
    target["stage_c_status"] = "invalid"
    target["support_low"], target["support_high"] = "0.0", "1.5"
    with slice_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(slice_rows[0]))
        writer.writeheader()
        writer.writerows(slice_rows)
    summary_path = job / "job_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["curve_file_sha256"]["curve_points.csv"] = hashlib.sha256(curve_path.read_bytes()).hexdigest()
    summary["curve_file_sha256"]["slice_metrics.csv"] = hashlib.sha256(slice_path.read_bytes()).hexdigest()
    summary_path.write_text(json.dumps(summary), encoding="utf-8")
    gate = module.collect(tmp_path, tmp_path / "aggregate", None, "a" * 40)
    assert not gate["workflow_contract_errors"]


def test_hybrid_workflow_contract_has_immutable_matrix_and_replay():
    module = load_module(COLLECTOR, "hybrid_collector_workflow")
    text = WORKFLOW.read_text(encoding="utf-8")
    assert "calculation_ref" in text
    assert "production_hybrid" in text and "independent_oracle" in text
    assert "timeout-minutes: 180" in text
    assert "aggregate_replay" in text and "source_run_id" in text
    assert "deep_oracle_run_id" in text
    assert "scope" in text and "targeted" in text and "full" in text
    assert "numeric_jobs" in text and "numeric_failures" in text
    assert "matplotlib==3.9.2" in text
    assert "--schema-version cep_maxwell_endpoint_production_shadow_v3" in text
    assert "--endpoint-mode" in text
    assert "--endpoint-policy bounded_zero_density_v1" in text
    assert "--candidate-policy unique_three_crossing_topology_v1" in text
    assert "--deep-input-dir deep-artifacts" in text
    assert "--deep-run-id \"$DEEP_RUN_ID\"" in text
    assert "--pattern '*aggregate'" in text
    assert "--run-mode aggregate_replay" in text
    assert len(module.XIS) * len(module.METHODS) == 9


def test_replay_accepts_numeric_success_when_source_aggregate_failed():
    text = WORKFLOW.read_text(encoding="utf-8")
    provenance_block = text.split("- name: Verify source run provenance", 1)[1].split("- name: Download source numerical artifacts", 1)[0]
    assert 'test "$(echo "$payload" | jq -r .status)" = completed' in provenance_block
    assert 'numeric_failures="$(echo "$payload"' in provenance_block
    assert 'test "$(echo "$payload" | jq -r .conclusion)" = success' not in provenance_block


def test_replay_source_accepts_legacy_hybrid_numeric_job_prefix():
    module = load_module(COLLECTOR, "hybrid_collector_source_replay_prefix")
    payload = {
        "status": "completed",
        "conclusion": "failure",
        "jobs": [
            {"name": f"hybrid shadow xi={xi} method={method}", "conclusion": "success"}
            for xi in ("-0.5", "0.0", "0.5")
            for method in ("production_hybrid", "memoized_dense", "independent_oracle")
        ],
    }
    assert module._source_run_completed_success(payload)


def test_endpoint_v3_schema_and_policy_gate_are_supported(tmp_path):
    module = load_module(COLLECTOR, "hybrid_collector_endpoint_v3")
    _matrix(tmp_path)
    for summary_path in tmp_path.glob("job-*/job_summary.json"):
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        summary["schema_version"] = "cep_maxwell_endpoint_production_shadow_v3"
        summary.setdefault("parameters", {})["rho_hybrid_candidate_policy"] = "unique_three_crossing_topology_v1"
        summary.setdefault("parameters", {})["rho_hybrid_endpoint_policy"] = "bounded_zero_density_v1"
        summary_path.write_text(json.dumps(summary), encoding="utf-8")
    gate = module.collect(
        tmp_path, tmp_path / "aggregate", None, "a" * 40,
        schema_version="cep_maxwell_endpoint_production_shadow_v3",
        endpoint_mode=True,
    )
    assert gate["verdict"] == "full_hybrid_candidate"
    assert not gate["endpoint_errors"]


def test_endpoint_local_v4_focused_schema_and_policy_gate_are_supported(tmp_path):
    module = load_module(COLLECTOR, "hybrid_collector_endpoint_local_v4")
    _focused_matrix(tmp_path)
    for summary_path in tmp_path.glob("job-*/job_summary.json"):
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        summary["schema_version"] = "cep_maxwell_endpoint_local_production_shadow_v4"
        summary.setdefault("parameters", {})["rho_hybrid_candidate_policy"] = "unique_three_crossing_sign_change_v2"
        summary.setdefault("parameters", {})["rho_hybrid_endpoint_policy"] = "three_crossing_endpoint_local_v2"
        summary_path.write_text(json.dumps(summary), encoding="utf-8")
    gate = module.collect(
        tmp_path,
        tmp_path / "aggregate",
        None,
        "a" * 40,
        scope="focused",
        schema_version="cep_maxwell_endpoint_local_production_shadow_v4",
        endpoint_mode=True,
        endpoint_policy="three_crossing_endpoint_local_v2",
        candidate_policy="unique_three_crossing_sign_change_v2",
    )
    assert gate["verdict"] == "focused_hybrid_candidate"
    assert gate["scope"] == "focused"
    assert gate["expected_anchor_count"] == 9
    assert not gate["endpoint_errors"]


def test_endpoint_local_support_envelope_covers_initial_bracket(tmp_path):
    module = load_module(COLLECTOR, "hybrid_collector_endpoint_support_envelope")
    _focused_matrix(tmp_path)
    job = tmp_path / "job-production_hybrid--0.5"
    slice_path = job / "slice_metrics.csv"
    slice_rows = list(csv.DictReader(slice_path.open(newline="", encoding="utf-8")))
    row = next(row for row in slice_rows if row["T_MeV"] == "5.0")
    row.update({
        "stage_c_status": "confirmed_first_order",
        "certificate_type": "endpoint_local_geometry_first_order",
        "result_status": "confirmed_first_order",
        "support_low": "0.0", "support_high": "0.00625",
        "endpoint_anchor_rho": "0.003125",
        "endpoint_left_bracket_low": "0.0", "endpoint_left_bracket_high": "0.00625",
        "endpoint_right_bracket_low": "2.8", "endpoint_right_bracket_high": "2.80625",
        "rho_hadron": "0.001", "rho_quark": "2.8",
        "maxwell_candidate_count": "1", "maxwell_crossing_count": "3",
    })
    with slice_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(slice_rows[0]))
        writer.writeheader()
        writer.writerows(slice_rows)
    curve_path = job / "curve_points.csv"
    curve_rows = list(csv.DictReader(curve_path.open(newline="", encoding="utf-8")))
    guard = next(row for row in curve_rows if row["T_MeV"] == "5.0")
    guard["rho"] = "0.003125"
    guard["sampling_role"] = "stage_c_guard"
    with curve_path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = list(curve_rows[0])
        if "sampling_role" not in fieldnames:
            fieldnames.append("sampling_role")
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(curve_rows)
    summary_path = job / "job_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    for name, path in (("slice_metrics.csv", slice_path), ("curve_points.csv", curve_path)):
        summary["curve_file_sha256"][name] = hashlib.sha256(path.read_bytes()).hexdigest()
    summary_path.write_text(json.dumps(summary), encoding="utf-8")
    for summary_path in tmp_path.glob("job-*/job_summary.json"):
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        summary["schema_version"] = "cep_maxwell_endpoint_local_production_shadow_v4"
        summary.setdefault("parameters", {})["rho_hybrid_endpoint_policy"] = "three_crossing_endpoint_local_v2"
        summary_path.write_text(json.dumps(summary), encoding="utf-8")
    gate = module.collect(
        tmp_path, tmp_path / "aggregate", None, "a" * 40,
        scope="focused", schema_version="cep_maxwell_endpoint_local_production_shadow_v4",
        endpoint_mode=True, endpoint_policy="three_crossing_endpoint_local_v2",
    )
    assert not any("support envelope" in error for error in gate["endpoint_errors"])


def test_endpoint_gate_finite_helper_accepts_serialized_values():
    module = load_module(COLLECTOR, "hybrid_collector_finite_helper")
    assert module._finite("0.125")
    assert not module._finite("NaN")
    assert not module._finite(None)


def test_endpoint_local_v4_workflow_contract_is_versioned_and_scoped():
    text = ENDPOINT_LOCAL_WORKFLOW.read_text(encoding="utf-8")
    assert "cep_maxwell_endpoint_local_production_shadow_v4" in text
    assert "three_crossing_endpoint_local_v2" in text
    assert "options: [focused, targeted, full]" in text
    assert "aggregate_replay" in text
    assert "rerun_failed_only" in text
    assert "timeout-minutes: 180" in text
    assert "GH_TOKEN: ${{ github.token }}" in text
    assert "--pattern '*required_three*'" in text
    assert "source_workflow_head_sha" in text
    assert "job summaries remain the authoritative calculation-SHA check" in text


def test_actions_cost_snapshot_is_provisional_until_authenticated_replay(tmp_path, monkeypatch):
    module = load_module(COLLECTOR, "hybrid_collector_actions_modes")
    payload = {
        "headSha": "a" * 40,
        "url": "https://example.invalid/run/1",
        "status": "completed",
        "conclusion": "success",
        "jobs": [{
            "name": "job",
            "databaseId": 1,
            "status": "completed",
            "conclusion": "success",
            "startedAt": "2026-08-05T00:00:00Z",
            "completedAt": "2026-08-05T00:01:00Z",
        }],
    }
    monkeypatch.setenv("GH_TOKEN", "token")
    monkeypatch.setattr(module.subprocess, "check_output", lambda *args, **kwargs: json.dumps(payload))
    provisional = module._actions("1", tmp_path / "provisional", run_mode="numerical")
    replay = module._actions("1", tmp_path / "replay", run_mode="aggregate_replay", source_run_id="1")
    assert provisional["snapshot_phase"] == "provisional"
    assert provisional["cost_snapshot_is_final"] is False
    assert replay["snapshot_phase"] == "final"
    assert replay["cost_snapshot_is_final"] is True


def test_replay_accepts_source_with_successful_numeric_jobs_and_failed_aggregate(tmp_path, monkeypatch):
    module = load_module(COLLECTOR, "hybrid_collector_numeric_source_replay")
    xis = ("-0.5", "0.0", "0.5")
    methods = ("production_hybrid", "memoized_dense", "independent_oracle")
    payload = {
        "headSha": "a" * 40,
        "url": "https://example.invalid/run/2",
        "status": "completed",
        "conclusion": "failure",
        "jobs": [
            {
                "name": f"endpoint-local v4 xi={xi} method={method}",
                "databaseId": index,
                "status": "completed",
                "conclusion": "success",
                "startedAt": "2026-08-05T00:00:00Z",
                "completedAt": "2026-08-05T00:01:00Z",
            }
            for index, (xi, method) in enumerate((
                (xi, method) for xi in xis for method in methods
            ), start=1)
        ],
    }
    monkeypatch.setenv("GH_TOKEN", "token")
    monkeypatch.setattr(module.subprocess, "check_output", lambda *args, **kwargs: json.dumps(payload))
    replay = module._actions("2", tmp_path / "replay", run_mode="aggregate_replay", source_run_id="2")
    assert replay["source_run_overall_success"] is False
    assert replay["source_run_numeric_jobs"] == 9
    assert replay["source_run_completed_success"] is True
    assert replay["cost_snapshot_is_final"] is True
