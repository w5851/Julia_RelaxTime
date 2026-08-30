#!/usr/bin/env python3
"""Build the solver-free PNJL Maxwell endpoint-limit evidence contract.

The source Actions artifact remains immutable.  This postprocessor validates
its hashes and numerical trace, then records the author-approved interpretation
that a shrinking left coexistence interval ending at rho=0 is sufficient
first-order evidence.  It never calls the equilibrium or Maxwell solver.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any, Iterable


SCHEMA_VERSION = "pnjl_maxwell_endpoint_limit_contract_v1"
SOURCE_SCHEMA = "pnjl_maxwell_endpoint_refinement_v1"
SOURCE_RUN_ID = "30980094983"
SOURCE_CALCULATION_SHA = "3217bed3635574f00c04cbee75e843b4c49451db"
SOURCE_WORKFLOW_HEAD_SHA = "df1fcdece7bc0888dd57c465c0015828743596c5"
SOURCE_VERDICT = "candidate_only_endpoint_inconclusive"
SOURCE_JOB_MANIFEST_SHA256 = (
    "ee2c74821d11d4b79cde932ccffc247e4602b7246697ebd6fe8ab4f937ffced5"
)
SOURCE_AGGREGATE_MANIFEST_SHA256 = (
    "e91e3cd3fb783213f7d6e5ec179e8c25f90763c20a292aa88326c2dbdc661af1"
)
EXPECTED_XI = -0.5
EXPECTED_T_MEV = 5.0
EXPECTED_BASE_STEP = 0.003125
EXPECTED_TARGETED_CAP = 12
AREA_SOLVER_TOL = 5.0e-6
OUTER_AREA_GATE = 5.0e-5
POSITION_TOL_MEV = 0.025
DENSITY_TOL = 0.0025
AUTHOR_DECISION_DATE = "2026-08-05"
AUTHOR_DECISION = (
    "The bounded endpoint interval 0 <= rho_h < 7.63e-7 is sufficient "
    "first-order evidence; a strictly positive lower bound is not required."
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def as_float(value: Any) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"expected finite numeric value, got {value!r}")
    return result


def as_int(value: Any) -> int:
    return int(value)


def as_bool(value: Any) -> bool:
    return str(value).strip().lower() == "true"


def git_head(project_root: Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(project_root), "rev-parse", "HEAD"],
            text=True,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def verify_file_hashes(directory: Path, expected: dict[str, str]) -> list[str]:
    errors: list[str] = []
    for name, digest in expected.items():
        path = directory / name
        if not path.is_file():
            errors.append(f"missing hashed source file: {name}")
            continue
        actual = sha256(path)
        if actual != str(digest):
            errors.append(f"source hash mismatch: {name}")
    return errors


def _strictly_decreasing(values: Iterable[float]) -> bool:
    sequence = list(values)
    return all(right < left for left, right in zip(sequence, sequence[1:]))


def evaluate_endpoint_limit_contract(
    metrics: list[dict[str, str]],
    costs: list[dict[str, str]],
    curves: list[dict[str, str]],
    source_manifest: dict[str, Any],
) -> tuple[dict[str, Any], list[str]]:
    """Evaluate the versioned endpoint-limit contract without solver calls."""

    errors: list[str] = []
    ordered = sorted(metrics, key=lambda row: as_int(row["level"]))
    expected_levels = list(range(EXPECTED_TARGETED_CAP + 1))
    levels = [as_int(row["level"]) for row in ordered]
    if levels != expected_levels:
        errors.append(f"refinement levels mismatch: {levels}")

    if not ordered:
        return {}, ["empty refinement trace"]

    unique_candidate_gate = all(
        row["status"] == "first_order"
        and row["reason"] == "unique_three_crossing_candidate"
        and as_int(row["candidate_count"]) == 1
        for row in ordered
    )
    if not unique_candidate_gate:
        errors.append("unique three-crossing candidate gate failed")

    endpoint_gate = all(
        as_bool(row["endpoint_dependent"])
        and not as_bool(row["positive_rho_bracket"])
        and math.isclose(as_float(row["bracket_low"]), 0.0, abs_tol=1.0e-15)
        and as_float(row["bracket_high"]) > 0.0
        and math.isclose(
            as_float(row["bracket_width"]),
            as_float(row["bracket_high"]) - as_float(row["bracket_low"]),
            rel_tol=0.0,
            abs_tol=1.0e-15,
        )
        for row in ordered
    )
    if not endpoint_gate:
        errors.append("endpoint-dependent zero-lower-bound gate failed")

    bracket_high = [as_float(row["bracket_high"]) for row in ordered]
    bracket_width = [as_float(row["bracket_width"]) for row in ordered]
    halving_gate = _strictly_decreasing(bracket_high) and all(
        math.isclose(right, left / 2.0, rel_tol=0.0, abs_tol=1.0e-15)
        for left, right in zip(bracket_high, bracket_high[1:])
    )
    if not halving_gate:
        errors.append("endpoint bracket does not shrink deterministically by bisection")

    expected_upper = EXPECTED_BASE_STEP / (2**EXPECTED_TARGETED_CAP)
    final_upper = bracket_high[-1]
    final_bound_gate = math.isclose(
        final_upper, expected_upper, rel_tol=0.0, abs_tol=1.0e-15
    )
    if not final_bound_gate:
        errors.append(
            f"final endpoint upper bound {final_upper} != expected {expected_upper}"
        )

    area_residuals = [abs(as_float(row["candidate_area"])) for row in ordered]
    area_solver_gate = all(value <= AREA_SOLVER_TOL for value in area_residuals)
    if not area_solver_gate:
        errors.append("candidate area exceeds the internal solver tolerance")

    geometry_rows = ordered[1:]
    geometry_gate = bool(geometry_rows) and all(
        as_bool(row["geometry_converged"])
        and as_float(row["position_error_MeV"]) <= POSITION_TOL_MEV
        and as_float(row["density_error"]) <= DENSITY_TOL
        and abs(as_float(row["candidate_area"])) <= OUTER_AREA_GATE
        for row in geometry_rows
    )
    if not geometry_gate:
        errors.append("cross-level Maxwell geometry gate failed")

    tail = ordered[-4:]
    tail_position = [as_float(row["position_error_MeV"]) for row in tail]
    tail_density = [as_float(row["density_error"]) for row in tail]
    convergence_trend_gate = (
        _strictly_decreasing(tail_position)
        and _strictly_decreasing(tail_density)
        and tail_position[-1] <= POSITION_TOL_MEV
        and tail_density[-1] <= DENSITY_TOL
    )
    if not convergence_trend_gate:
        errors.append("final Maxwell coordinate convergence trend failed")

    rho_hadron = [as_float(row["rho_hadron"]) for row in ordered]
    crossing_bound_gate = all(
        0.0 <= crossing <= upper
        for crossing, upper in zip(rho_hadron, bracket_high)
    )
    if not crossing_bound_gate:
        errors.append("interpolated left crossing escapes its endpoint interval")

    trace_accounting_gate = all(
        as_int(row["targeted_additions"]) == level
        and as_int(row["unique_solves"]) == 1281 + level
        and as_int(row["solver_failures"]) == 0
        for level, row in enumerate(ordered)
    )
    if not trace_accounting_gate:
        errors.append("per-level targeted-point or solver-work accounting gate failed")

    curve_keys: set[tuple[float, float, float]] = set()
    curve_gate = True
    for row in curves:
        try:
            key = (
                as_float(row["xi"]),
                as_float(row["T_MeV"]),
                as_float(row["rho"]),
            )
            as_float(row["muq_MeV"])
            as_float(row["pressure_fm4"])
            as_float(row["residual_norm"])
        except (KeyError, ValueError):
            curve_gate = False
            continue
        if key in curve_keys:
            curve_gate = False
        curve_keys.add(key)
        curve_gate = curve_gate and as_bool(row["converged"]) and as_bool(row["finite"])
        curve_gate = curve_gate and math.isclose(
            key[0], EXPECTED_XI, abs_tol=1.0e-12
        )
        curve_gate = curve_gate and math.isclose(
            key[1], EXPECTED_T_MEV, abs_tol=1.0e-12
        )
        curve_gate = curve_gate and row["calculation_sha"] == SOURCE_CALCULATION_SHA
        curve_gate = curve_gate and row["workflow_head_sha"] == SOURCE_WORKFLOW_HEAD_SHA
    if len(curves) != 1281 + EXPECTED_TARGETED_CAP or not curve_gate:
        errors.append("curve finiteness, provenance, count, or unique-key gate failed")

    if len(costs) != 1:
        errors.append("expected exactly one method-cost row")
        cost = {}
    else:
        cost = costs[0]
    zero_fields = (
        "nonconverged_attempts",
        "fallback_count",
        "governed_rescue_count",
        "retry_count",
        "exception_count",
    )
    cost_gate = bool(cost)
    if cost:
        cost_gate = cost_gate and as_int(cost["targeted_additions"]) == EXPECTED_TARGETED_CAP
        cost_gate = cost_gate and as_int(cost["unique_solves"]) == len(curves)
        cost_gate = cost_gate and as_int(cost["point_requests"]) == len(curves)
        cost_gate = cost_gate and as_int(cost["fixedrho_requests"]) == len(curves)
        cost_gate = cost_gate and as_int(cost["cache_hits"]) == 0
        cost_gate = cost_gate and all(as_int(cost[field]) == 0 for field in zero_fields)
    if not cost_gate:
        errors.append("solver-work conservation or failure gate failed")

    manifest_gate = (
        source_manifest.get("schema_version") == SOURCE_SCHEMA
        and source_manifest.get("verdict") == SOURCE_VERDICT
        and source_manifest.get("calculation_sha") == SOURCE_CALCULATION_SHA
        and source_manifest.get("workflow_head_sha") == SOURCE_WORKFLOW_HEAD_SHA
        and source_manifest.get("reference_write") is False
        and math.isclose(
            as_float(source_manifest.get("base_rho_step")),
            EXPECTED_BASE_STEP,
            abs_tol=1.0e-15,
        )
    )
    if not manifest_gate:
        errors.append("source manifest contract failed")

    passed = not errors
    certificate = {
        "contract_version": "endpoint_limited_first_order_v1",
        "contract_pass": passed,
        "physical_classification": (
            "confirmed_first_order" if passed else "ambiguous_near_critical"
        ),
        "certificate_kind": (
            "endpoint_limited_first_order" if passed else "none"
        ),
        "source_verdict": SOURCE_VERDICT,
        "xi": EXPECTED_XI,
        "T_MeV": EXPECTED_T_MEV,
        "candidate_mu_MeV": as_float(ordered[-1]["candidate_mu"]),
        "candidate_area_residual": as_float(ordered[-1]["candidate_area"]),
        "rho_hadron_representative": 0.0,
        "rho_hadron_lower_bound": 0.0,
        "rho_hadron_upper_bound": final_upper,
        "rho_hadron_interpolated_diagnostic": as_float(ordered[-1]["rho_hadron"]),
        "rho_quark": as_float(ordered[-1]["rho_quark"]),
        "targeted_cap": EXPECTED_TARGETED_CAP,
        "base_rho_step": EXPECTED_BASE_STEP,
        "final_position_error_MeV": as_float(ordered[-1]["position_error_MeV"]),
        "final_density_error": as_float(ordered[-1]["density_error"]),
        "max_abs_area_residual": max(area_residuals),
        "unique_candidate_gate": unique_candidate_gate,
        "endpoint_gate": endpoint_gate,
        "bracket_halving_gate": halving_gate,
        "final_bound_gate": final_bound_gate,
        "area_solver_gate": area_solver_gate,
        "geometry_gate": geometry_gate,
        "convergence_trend_gate": convergence_trend_gate,
        "crossing_bound_gate": crossing_bound_gate,
        "trace_accounting_gate": trace_accounting_gate,
        "curve_gate": curve_gate,
        "cost_gate": cost_gate,
        "strict_positive_rho_claimed": False,
        "author_decision_status": "accepted",
        "author_decision_date": AUTHOR_DECISION_DATE,
    }
    return certificate, errors


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        with path.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("\n")
        return
    fieldnames = fields or list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(value)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_text(path, json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def finite_or_empty(value: str) -> str:
    """Keep finite values and normalize unavailable source sentinels to empty CSV cells."""

    try:
        return value if math.isfinite(float(value)) else ""
    except (TypeError, ValueError):
        return ""


def plot_trace(metrics: list[dict[str, str]], output: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ordered = sorted(metrics, key=lambda row: as_int(row["level"]))
    levels = [as_int(row["level"]) for row in ordered]
    upper = [as_float(row["bracket_high"]) for row in ordered]
    mu = [as_float(row["candidate_mu"]) for row in ordered]
    rho_q = [as_float(row["rho_quark"]) for row in ordered]
    area = [abs(as_float(row["candidate_area"])) for row in ordered]
    final_mu, final_rho_q = mu[-1], rho_q[-1]

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.2), constrained_layout=True)
    ax = axes[0, 0]
    ax.semilogy(levels, upper, marker="o", color="#264653")
    ax.set_xlabel("targeted endpoint points")
    ax.set_ylabel("rho_h upper bound")
    ax.set_title("Endpoint bracket contraction")
    ax.grid(alpha=0.25)

    ax = axes[0, 1]
    mu_delta = [max(abs(value - final_mu), 1.0e-16) for value in mu]
    ax.semilogy(levels, mu_delta, marker="o", color="#e76f51")
    ax.set_xlabel("targeted endpoint points")
    ax.set_ylabel("abs(mu_M - final mu_M) [MeV]")
    ax.set_title("Maxwell chemical-potential convergence")
    ax.grid(alpha=0.25)

    ax = axes[1, 0]
    rho_delta = [max(abs(value - final_rho_q), 1.0e-16) for value in rho_q]
    ax.semilogy(levels, rho_delta, marker="o", color="#2a9d8f")
    ax.set_xlabel("targeted endpoint points")
    ax.set_ylabel("abs(rho_q - final rho_q)")
    ax.set_title("Right coexistence-density convergence")
    ax.grid(alpha=0.25)

    ax = axes[1, 1]
    ax.semilogy(levels, area, marker="o", color="#f4a261", label="area residual")
    ax.axhline(AREA_SOLVER_TOL, color="#9b2226", linestyle="--", label="solver tolerance")
    ax.axhline(OUTER_AREA_GATE, color="#6c757d", linestyle=":", label="outer gate")
    ax.set_xlabel("targeted endpoint points")
    ax.set_ylabel("absolute area residual")
    ax.set_title("Maxwell area contract")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def build_evidence(
    input_dir: Path,
    aggregate_dir: Path,
    output_dir: Path,
    project_root: Path,
    source_run_id: str,
) -> dict[str, Any]:
    errors: list[str] = []
    if source_run_id != SOURCE_RUN_ID:
        errors.append(f"source run must be {SOURCE_RUN_ID}, got {source_run_id}")

    source_manifest = load_json(input_dir / "manifest.json")
    aggregate_manifest = load_json(aggregate_dir / "manifest.json")
    selected_source = load_json(input_dir / "selected_policy.json")
    if sha256(input_dir / "manifest.json") != SOURCE_JOB_MANIFEST_SHA256:
        errors.append("source job manifest SHA256 mismatch")
    if sha256(aggregate_dir / "manifest.json") != SOURCE_AGGREGATE_MANIFEST_SHA256:
        errors.append("source aggregate manifest SHA256 mismatch")
    errors.extend(verify_file_hashes(input_dir, source_manifest.get("files", {})))
    errors.extend(verify_file_hashes(input_dir, aggregate_manifest.get("files", {})))

    aggregate_gate = (
        aggregate_manifest.get("schema_version") == SOURCE_SCHEMA
        and aggregate_manifest.get("verdict") == SOURCE_VERDICT
        and aggregate_manifest.get("workflow_contract_errors") == []
        and aggregate_manifest.get("calculation_sha") == SOURCE_CALCULATION_SHA
        and aggregate_manifest.get("workflow_head_sha") == SOURCE_WORKFLOW_HEAD_SHA
        and aggregate_manifest.get("reference_write") is False
        and selected_source.get("verdict") == SOURCE_VERDICT
    )
    if not aggregate_gate:
        errors.append("aggregate/source-policy provenance gate failed")

    curves = read_rows(input_dir / "curve_points.csv")
    metrics = read_rows(input_dir / "slice_metrics.csv")
    frontier = read_rows(input_dir / "policy_frontier.csv")
    costs = read_rows(input_dir / "method_costs.csv")
    if len(curves) != as_int(aggregate_manifest.get("curve_rows", -1)):
        errors.append("aggregate curve-row count mismatch")
    if len(metrics) != as_int(aggregate_manifest.get("metric_rows", -1)):
        errors.append("aggregate metric-row count mismatch")
    if len(frontier) != as_int(aggregate_manifest.get("frontier_rows", -1)):
        errors.append("aggregate frontier-row count mismatch")

    certificate, contract_errors = evaluate_endpoint_limit_contract(
        metrics, costs, curves, source_manifest
    )
    errors.extend(contract_errors)
    verdict = "endpoint_limited_first_order_candidate" if not errors else "contract_inconclusive"

    output_dir.mkdir(parents=True, exist_ok=True)
    trace_fields = [
        "level",
        "targeted_additions",
        "status",
        "reason",
        "candidate_count",
        "candidate_mu",
        "candidate_area",
        "rho_hadron",
        "rho_quark",
        "endpoint_dependent",
        "bracket_low",
        "bracket_high",
        "bracket_width",
        "positive_rho_bracket",
        "position_error_MeV",
        "density_error",
        "geometry_converged",
        "unique_solves",
        "solver_failures",
    ]
    trace = [{field: row[field] for field in trace_fields} for row in metrics]
    for row in trace:
        row["position_error_MeV"] = finite_or_empty(row["position_error_MeV"])
        row["density_error"] = finite_or_empty(row["density_error"])
    write_csv(output_dir / "tables" / "refinement_trace.csv", trace, trace_fields)
    write_csv(output_dir / "tables" / "endpoint_limit_certificate.csv", [certificate])
    write_csv(output_dir / "tables" / "cost_summary.csv", costs)

    inventory = []
    for name in sorted(aggregate_manifest.get("files", {})):
        inventory.append(
            {
                "source_file": name,
                "sha256": aggregate_manifest["files"][name],
                "stored_in_repository": "false",
                "role": "immutable Actions input",
            }
        )
    write_csv(output_dir / "tables" / "source_inventory.csv", inventory)

    claims = [
        {
            "claim_id": "source_integrity",
            "claim": "The endpoint artifact is hash-complete, finite, converged, and key-unique.",
            "status": "pass" if not errors else "inconclusive",
            "evidence": "provenance.json; tables/cost_summary.csv",
            "fields": "source_*_manifest_sha256; unique_solves; failure counters",
            "boundary": "run 30980094983 only",
        },
        {
            "claim_id": "unique_candidate",
            "claim": "All refinement levels retain one stable three-crossing Maxwell candidate.",
            "status": "pass" if certificate.get("unique_candidate_gate") else "inconclusive",
            "evidence": "tables/refinement_trace.csv",
            "fields": "level; status; reason; candidate_count",
            "boundary": "single xi=-0.5, T=5 MeV anchor",
        },
        {
            "claim_id": "endpoint_limit",
            "claim": "The hadronic coexistence density is bounded by 0 <= rho_h < 7.63e-7.",
            "status": "pass" if certificate.get("contract_pass") else "inconclusive",
            "evidence": "tables/endpoint_limit_certificate.csv",
            "fields": "rho_hadron_lower_bound; rho_hadron_upper_bound; bracket_halving_gate",
            "boundary": "endpoint-limit certificate; strict positive rho_h is not claimed",
        },
        {
            "claim_id": "physical_classification",
            "claim": "The author-approved endpoint-limit contract supports confirmed first order.",
            "status": "author_confirmed" if certificate.get("contract_pass") else "inconclusive",
            "evidence": "selected_policy.json; provenance.json",
            "fields": "final_physical_classification; author_decision",
            "boundary": "policy decision dated 2026-08-05",
        },
        {
            "claim_id": "production_integration",
            "claim": "The endpoint-limit certificate is implemented in production.",
            "status": "not_claimed",
            "evidence": "README.md",
            "fields": "evidence boundary",
            "boundary": "requires a separate public-core/hybrid production PR and shadow validation",
        },
        {
            "claim_id": "cep_or_reference",
            "claim": "This single anchor resolves the CEP or promotes a phase reference.",
            "status": "not_claimed",
            "evidence": "README.md",
            "fields": "evidence boundary",
            "boundary": "single-temperature diagnostic only",
        },
    ]
    write_csv(output_dir / "tables" / "claim_ledger.csv", claims)

    selected_policy = {
        "schema_version": SCHEMA_VERSION,
        "verdict": verdict,
        "certificate_policy": "endpoint_limited_first_order_v1",
        "final_physical_classification": certificate.get("physical_classification"),
        "rho_hadron_scalar_semantics": "endpoint_limit_representative_zero",
        "rho_hadron_representative": certificate.get("rho_hadron_representative"),
        "rho_hadron_interval": [
            certificate.get("rho_hadron_lower_bound"),
            certificate.get("rho_hadron_upper_bound"),
        ],
        "base_rho_step": EXPECTED_BASE_STEP,
        "targeted_cap": EXPECTED_TARGETED_CAP,
        "candidate_policy": "unique_three_crossing_topology_v1",
        "endpoint_policy": "bounded_zero_density_v1",
        "require_deterministic_bracket_halving": True,
        "require_cross_level_geometry": True,
        "require_solver_failures_zero": True,
        "strict_positive_rho_claimed": False,
        "author_decision_status": "accepted",
        "author_decision_date": AUTHOR_DECISION_DATE,
    }
    write_json(output_dir / "selected_policy.json", selected_policy)

    provenance = {
        "schema_version": SCHEMA_VERSION,
        "source_run_id": SOURCE_RUN_ID,
        "source_run_url": f"https://github.com/w5851/Julia_RelaxTime/actions/runs/{SOURCE_RUN_ID}",
        "source_calculation_sha": SOURCE_CALCULATION_SHA,
        "source_workflow_head_sha": SOURCE_WORKFLOW_HEAD_SHA,
        "source_verdict": SOURCE_VERDICT,
        "source_job_manifest_sha256": sha256(input_dir / "manifest.json"),
        "source_aggregate_manifest_sha256": sha256(aggregate_dir / "manifest.json"),
        "producer_head_sha": git_head(project_root),
        "producer_script": "scripts/analysis/build_pnjl_maxwell_endpoint_limit_contract.py",
        "producer_script_sha256": sha256(Path(__file__).resolve()),
        "solver_called": False,
        "reference_write": False,
        "raw_curves_external": True,
        "author_decision": {
            "status": "accepted",
            "date": AUTHOR_DECISION_DATE,
            "statement": AUTHOR_DECISION,
            "scope": "endpoint_limited_first_order_v1",
        },
    }
    write_json(output_dir / "provenance.json", provenance)

    figure_path = output_dir / "figures" / "endpoint_limit_convergence.png"
    plot_trace(trace, figure_path)
    plot_manifest = {
        "schema_version": SCHEMA_VERSION,
        "figures": [
            {
                "path": "figures/endpoint_limit_convergence.png",
                "sha256": sha256(figure_path),
                "inputs": ["tables/refinement_trace.csv"],
                "claim_boundary": "single-anchor numerical convergence only",
            }
        ],
        "command": (
            "python scripts/analysis/build_pnjl_maxwell_endpoint_limit_contract.py "
            "--input-dir <job-artifact> --aggregate-dir <aggregate-artifact> "
            "--output-dir docs/analysis/pnjl/cep_maxwell/maxwell_contracts/pnjl_maxwell_endpoint_limit_contract_v1"
        ),
    }
    write_json(output_dir / "figures" / "plot_manifest.json", plot_manifest)

    readme = f"""# PNJL Maxwell endpoint-limit contract v1

本分析包固定读取 GitHub Actions run `{SOURCE_RUN_ID}` 的单点 endpoint refinement
artifact，并在不调用 equilibrium/Maxwell solver 的前提下形成端点极限证书。

## 结论

- source artifact 原 verdict 保持 `{SOURCE_VERDICT}`，没有改写历史产物。
- 13 个 refinement levels 均保持唯一三交点 Maxwell candidate。
- 左侧 bracket 从 `[0, {EXPECTED_BASE_STEP}]` 确定性收缩至
  `[0, {certificate.get('rho_hadron_upper_bound')}]`。
- 最终 `mu_M={certificate.get('candidate_mu_MeV')} MeV`，
  `rho_q={certificate.get('rho_quark')}`，面积残差
  `{certificate.get('candidate_area_residual')}`；geometry 与 solver-work gate 均通过。
- 作者已确认该有界端点区间足以作为一阶相变证据。因此 derived verdict 为
  `{verdict}`，三态分类为 `{certificate.get('physical_classification')}`。

## 合同语义

`endpoint_limited_first_order` 不声称存在严格正的 `rho_h` 下界。兼容 scalar 字段在后续
production 实现中应以 `rho_hadron=0.0` 表示端点极限，同时把真实区间
`[0, rho_hadron_upper_bound]` 写入 diagnostics/CSV。该证书不增加第四种 CEP 物理状态；
最终仍映射到 `confirmed_first_order`。

## 证据边界

本包只覆盖 `(xi,T)=(-0.5,5 MeV)`，不构成 CEP resolved、24-anchor shadow 通过、
phase-reference 晋升或 transport 输入。公共 Maxwell 核心和 hybrid production 尚未修改；
完整原始 `curve_points.csv` 仅保留在 Actions/local artifact，通过 provenance/hash 追溯。
"""
    write_text(output_dir / "README.md", readme)

    audit = f"""# 审计记录

- source run：`{SOURCE_RUN_ID}`
- calculation SHA：`{SOURCE_CALCULATION_SHA}`
- workflow head SHA：`{SOURCE_WORKFLOW_HEAD_SHA}`
- curve rows：`{len(curves)}`；metric rows：`{len(metrics)}`；frontier rows：`{len(frontier)}`
- targeted additions：`{EXPECTED_TARGETED_CAP}`
- source job manifest SHA256：`{provenance['source_job_manifest_sha256']}`
- source aggregate manifest SHA256：`{provenance['source_aggregate_manifest_sha256']}`
- solver called：`false`；reference write：`false`
- contract errors：`{errors}`

作者裁决只接受 `rho_h` 的有界零密度端点极限，不把该区间伪装成严格正密度根，也不
替代后续 production integration 与 targeted/full shadow gate。
"""
    write_text(output_dir / "AUDIT.md", audit)

    files: dict[str, str] = {}
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and path.name != "manifest.json":
            files[path.relative_to(output_dir).as_posix()] = sha256(path)
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "verdict": verdict,
        "physical_classification": certificate.get("physical_classification"),
        "certificate_kind": certificate.get("certificate_kind"),
        "source_run_id": SOURCE_RUN_ID,
        "source_calculation_sha": SOURCE_CALCULATION_SHA,
        "source_workflow_head_sha": SOURCE_WORKFLOW_HEAD_SHA,
        "source_verdict": SOURCE_VERDICT,
        "solver_called": False,
        "reference_write": False,
        "raw_curves_external": True,
        "contract_errors": errors,
        "author_decision_status": "accepted",
        "files": files,
    }
    write_json(output_dir / "manifest.json", manifest)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--aggregate-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--source-run-id", default=SOURCE_RUN_ID)
    args = parser.parse_args()
    manifest = build_evidence(
        args.input_dir.resolve(),
        args.aggregate_dir.resolve(),
        args.output_dir.resolve(),
        args.project_root.resolve(),
        args.source_run_id,
    )
    print(json.dumps(manifest, ensure_ascii=False))
    return 0 if manifest["verdict"] == "endpoint_limited_first_order_candidate" else 2


if __name__ == "__main__":
    raise SystemExit(main())
