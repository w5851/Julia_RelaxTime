"""Validate and plot the fixed low/mid/high CEP manual-bisection artifacts.

This collector is solver-free. It deliberately records manual decisions as
``manual_pending`` and never uses oracle labels to select a temperature.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
SCHEMA_VERSION = "pnjl_c2_cep_manual_bisection_job_v1"
AUDIT_SCHEMA_VERSION = "pnjl_c2_cep_manual_bisection_audit_v1"
XIS = (0.125, 0.39375, 0.5)
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


def float_or_none(value: Any) -> float | None:
    return float(value) if finite(value) else None


def slug(value: float) -> str:
    return str(value).replace("-", "m").replace(".", "p")


def _declared_file(manifest: dict[str, Any], root: Path, name: str) -> Path:
    path = root / name
    if not path.is_file():
        raise ValueError(f"missing {name}: {root}")
    declared = manifest.get("files", {}).get(name)
    if declared != sha256(path):
        raise ValueError(f"hash mismatch for {name}: {root}")
    return path


def _phase_summary(slice_dir: Path, method: str) -> dict[str, Any]:
    path = slice_dir / method / "phase_summary.json"
    return read_json(path) if path.is_file() else {}


def _sweep_record(summary: dict[str, Any]) -> dict[str, Any]:
    records = summary.get("stats", {}).get("sweep_records", [])
    return records[0] if records else {}


def _candidate_fields(summary: dict[str, Any]) -> dict[str, Any]:
    record = _sweep_record(summary)
    return {
        "crossing_count": record.get("maxwell_crossing_count"),
        "candidate_count": record.get("maxwell_candidate_count"),
        "mu_transition_MeV": record.get("mu_transition_MeV"),
        "rho_hadron": record.get("rho_hadron"),
        "rho_quark": record.get("rho_quark"),
        "area_residual": record.get("area_residual"),
        "area_gate": record.get("maxwell_area_gate"),
        "position_error_MeV": record.get("position_error_MeV"),
        "density_error": record.get("density_error"),
        "geometry_converged": record.get("geometry_converged"),
        "phase_status": summary.get("conclusion", {}).get("phase_structure"),
        "scan_failure": summary.get("stats", {}).get("scan_failure"),
        "scan_success": summary.get("stats", {}).get("scan_success"),
        "hybrid_stage_b_status": record.get("hybrid_stage_b_status"),
        "hybrid_stage_c_status": record.get("hybrid_stage_c_status"),
        "hybrid_detailed_reason": record.get("hybrid_upgrade_reason"),
    }


def validate_artifact(root: Path, expected_source_run: str, expected_head: str) -> dict[str, Any]:
    manifest = read_json(root / "manifest.json")
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"unexpected schema: {root}")
    if manifest.get("calculation_sha", "").lower() != CALCULATION_SHA:
        raise ValueError(f"calculation SHA mismatch: {root}")
    if expected_head and manifest.get("workflow_head_sha") != expected_head:
        raise ValueError(f"workflow head mismatch: {root}")
    if str(manifest.get("source_run_id")) != str(expected_source_run):
        raise ValueError(f"source run mismatch: {root}")
    if manifest.get("solver_called") is not True:
        raise ValueError(f"numerical artifact must record solver_called=true: {root}")
    if manifest.get("manual_decision_required") is not True:
        raise ValueError(f"manual decision contract missing: {root}")
    if manifest.get("oracle_labels_used_for_routing") is not False:
        raise ValueError(f"oracle routing leakage: {root}")
    fine_path = _declared_file(manifest, root, "fine_pool.csv")
    trace_path = _declared_file(manifest, root, "manual_bisection_trace.csv")
    cost_path = _declared_file(manifest, root, "method_costs.csv")
    fine = read_csv(fine_path)
    trace = read_csv(trace_path)
    costs = read_csv(cost_path)
    xi = float(manifest["xi"])
    if xi not in XIS:
        raise ValueError(f"unexpected xi: {xi}")
    if len(trace) != 3 or [row["role"] for row in trace] != ["low", "midpoint", "high"]:
        raise ValueError(f"manual trace must contain low/midpoint/high: {root}")
    if len(fine) != 3 * RHO_COUNT:
        raise ValueError(f"unexpected fine pool rows {len(fine)}: {root}")
    keys = set()
    for row in fine:
        key = (round(float(row["xi"]), 8), round(float(row["T_MeV"]), 8),
               round(float(row["rho"]) / RHO_STEP))
        if key in keys:
            raise ValueError(f"duplicate fine-pool key {key}: {root}")
        keys.add(key)
        if not finite(row.get("muq_MeV")) or row.get("converged", "").lower() != "true" or row.get("finite", "").lower() != "true":
            raise ValueError(f"non-finite/non-converged fine-pool row: {root}")
        if row.get("calculation_sha", "").lower() != CALCULATION_SHA:
            raise ValueError(f"fine-pool SHA mismatch: {root}")
    if any(int(row.get("failed_points", "0")) != 0 for row in costs):
        raise ValueError(f"failed solver points in cost table: {root}")
    for row in costs:
        if int(row["point_requests"]) != int(row["unique_solves"]) + int(row["cache_hits"]):
            raise ValueError(f"cache cost does not reconcile: {root}")
    return {"manifest": manifest, "fine": fine, "trace": trace, "costs": costs}


def _hybrid_rows(slice_dir: Path) -> list[dict[str, str]]:
    path = slice_dir / "hybrid" / "trho_scan.csv"
    if not path.is_file():
        return []
    return read_csv(path)


def _mu(row: dict[str, str]) -> float | None:
    return float_or_none(row.get("mu_avg_MeV", row.get("muq_MeV")))


def plot_slice(root: Path, row: dict[str, str], output_dir: Path) -> dict[str, Any]:
    xi = float(row["xi"])
    temperature = float(row["T_MeV"])
    slice_dir = root / "slices" / f"T_{str(temperature).replace('.', 'p')}"
    oracle = _phase_summary(slice_dir, "oracle")
    candidate = _candidate_fields(oracle)
    oracle_rows = [r for r in read_csv(root / "fine_pool.csv")
                   if abs(float(r["T_MeV"]) - temperature) < 1e-8]
    oracle_rows.sort(key=lambda item: float(item["rho"]))
    x = [float(item["rho"]) for item in oracle_rows]
    y = [_mu(item) for item in oracle_rows]
    if any(value is None for value in y):
        raise ValueError(f"non-finite curve for {root} T={temperature}")
    y_float = [float(value) for value in y if value is not None]
    hybrid_rows = _hybrid_rows(slice_dir)
    hx = [float(item["rho"]) for item in hybrid_rows if _mu(item) is not None]
    hy = [float(_mu(item)) for item in hybrid_rows if _mu(item) is not None]
    rho_h = float_or_none(candidate.get("rho_hadron"))
    rho_q = float_or_none(candidate.get("rho_quark"))
    mu_m = float_or_none(candidate.get("mu_transition_MeV"))
    if rho_h is not None and rho_q is not None:
        xmin = max(0.0, min(rho_h, rho_q) - 0.22)
        xmax = min(4.0, max(rho_h, rho_q) + 0.22)
    else:
        xmin, xmax = 0.0, 4.0
    local_values = [value for rho, value in zip(x, y_float) if xmin <= rho <= xmax]
    center = mu_m if mu_m is not None else 0.5 * (min(y_float) + max(y_float))
    span = max(local_values) - min(local_values) if local_values else max(y_float) - min(y_float)
    ypad = max(0.12 * span, 0.05, 0.0005)
    ymin, ymax = center - ypad, center + ypad

    fig, (ax, inset) = plt.subplots(1, 2, figsize=(12.0, 4.7),
                                    gridspec_kw={"width_ratios": [1.35, 1.0]})
    ax.plot(x, y_float, color="#1f4e79", linewidth=1.0, label="full-fine oracle")
    if hx:
        ax.plot(hx, hy, "o", color="#c44e52", markersize=2.5, alpha=0.7,
                label="hybrid Stage-A/B points")
    if mu_m is not None:
        ax.axhline(mu_m, color="#333333", linestyle="--", linewidth=0.8,
                   label="oracle Maxwell μ")
    if rho_h is not None and rho_q is not None and mu_m is not None:
        ax.scatter([rho_h, rho_q], [mu_m, mu_m], color="#111111", s=20, zorder=4)
    ax.set_xlabel(r"$\rho$")
    ax.set_ylabel(r"$\mu$ (MeV)")
    ax.set_title(f"xi={xi:g}, T={temperature:g} MeV")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8, loc="best")

    inset.plot(x, y_float, color="#1f4e79", linewidth=1.0)
    if hx:
        inset.plot(hx, hy, "o", color="#c44e52", markersize=3.0, alpha=0.7)
    if mu_m is not None:
        inset.axhline(mu_m, color="#333333", linestyle="--", linewidth=0.8)
    if rho_h is not None and rho_q is not None and mu_m is not None:
        inset.scatter([rho_h, rho_q], [mu_m, mu_m], color="#111111", s=24, zorder=4)
    inset.set_xlim(xmin, xmax)
    inset.set_ylim(ymin, ymax)
    inset.set_xlabel(r"$\rho$")
    inset.set_ylabel(r"$\mu$ (MeV)")
    inset.set_title("local Maxwell region")
    inset.grid(alpha=0.2)
    fig.suptitle(f"manual CEP bisection audit: {row['role']}", y=1.02)
    fig.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    out = output_dir / f"rho_mu_xi_{slug(xi)}_T_{slug(temperature)}_{row['role']}.png"
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return {
        "xi": xi, "T_MeV": temperature, "role": row["role"],
        "crossing_count": candidate["crossing_count"],
        "candidate_count": candidate["candidate_count"],
        "mu_transition_MeV": candidate["mu_transition_MeV"],
        "rho_hadron": candidate["rho_hadron"], "rho_quark": candidate["rho_quark"],
        "hybrid_stage_b_status": candidate["hybrid_stage_b_status"],
        "hybrid_stage_c_status": candidate["hybrid_stage_c_status"],
        "hybrid_detailed_reason": candidate["hybrid_detailed_reason"],
        "plot": str(out.name),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--source-run-id", required=True)
    parser.add_argument("--expected-workflow-head", default="")
    args = parser.parse_args()

    artifacts = []
    for xi in XIS:
        root = args.input_dir / f"xi_{xi}"
        if not root.is_dir():
            raise SystemExit(f"missing artifact directory: {root}")
        artifacts.append(validate_artifact(root, args.source_run_id, args.expected_workflow_head))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    figures = args.output_dir / "figures"
    audit_rows: list[dict[str, Any]] = []
    candidate_rows: list[dict[str, Any]] = []
    total_cost_rows: list[dict[str, str]] = []
    for artifact in artifacts:
        manifest = artifact["manifest"]
        root = args.input_dir / f"xi_{manifest['xi']}"
        for row in artifact["trace"]:
            plotted = plot_slice(root, row, figures)
            audit_rows.append({**row, **plotted})
        total_cost_rows.extend(artifact["costs"])
        candidate_rows.extend([row for row in audit_rows if row["xi"] == float(manifest["xi"])])

    fields = sorted({key for row in audit_rows for key in row})
    with (args.output_dir / "manual_bisection_audit.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows({key: row.get(key, "") for key in fields} for row in audit_rows)
    cost_fields = sorted({key for row in total_cost_rows for key in row})
    with (args.output_dir / "method_costs.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=cost_fields)
        writer.writeheader()
        writer.writerows({key: row.get(key, "") for key in cost_fields} for row in total_cost_rows)

    manifest = {
        "schema_version": AUDIT_SCHEMA_VERSION,
        "source_run_id": str(args.source_run_id),
        "calculation_sha": CALCULATION_SHA,
        "workflow_head_sha": args.expected_workflow_head,
        "manual_decision_required": True,
        "oracle_labels_used_for_routing": False,
        "solver_called": False,
        "input_artifacts": [
            {"xi": item["manifest"]["xi"], "plan_sha256": item["manifest"].get("plan_sha256"),
             "manifest_sha256": sha256(args.input_dir / f"xi_{item['manifest']['xi']}" / "manifest.json")}
            for item in artifacts
        ],
        "files": {},
        "plots": len(list(figures.glob("*.png"))),
    }
    for path in [args.output_dir / "manual_bisection_audit.csv", args.output_dir / "method_costs.csv"]:
        manifest["files"][path.name] = sha256(path)
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    (args.output_dir / "claim_ledger.md").write_text(
        "# C2 CEP manual bisection audit v1\n\n"
        "- This package is solver-free after the numerical source run.\n"
        "- The nine slices are fixed low/midpoint/high observations of the C2 brackets.\n"
        "- `manual_pending` is intentional; no oracle label selected a temperature.\n"
        "- A unique oracle Maxwell candidate is diagnostic evidence, not a hybrid production certificate.\n"
        "- CEP acceptance still requires author bracket decisions and width <= 0.1 MeV.\n",
        encoding="utf-8",
    )
    (args.output_dir / "README.md").write_text(
        "# C2 CEP manual bisection audit v1\n\n"
        "This package contains full-fine oracle curves and hybrid overlays for the low, midpoint, "
        "and high temperatures of the three C2 CEP brackets. Review the local panels and update "
        "the bracket trace only after author classification.\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
