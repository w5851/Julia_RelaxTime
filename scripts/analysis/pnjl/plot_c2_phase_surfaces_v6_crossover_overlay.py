#!/usr/bin/env python3
"""Build the v6 diagnostic phase view from the v5 postprocessed package.

The v5 package is the immutable baseline.  This script only appends the
solver-produced crossover endpoint expansion and never rebuilds a boundary,
Maxwell candidate, CEP bracket, or spinodal from raw C2 rows.  Plotting uses
native ordered segments; gaps remain visible and no interpolation or
triangulation is performed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import shutil
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


SCHEMA = "pnjl_c2_phase_surfaces_diagnostic_v6_crossover_overlay"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v5-root", type=Path, required=True)
    parser.add_argument("--expansion-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--spinodals", type=Path, default=None)
    parser.add_argument("--max-mu-gap", type=float, default=55.0)
    parser.add_argument("--max-temperature-gap", type=float, default=8.0)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"CSV has no rows: {path}")
    return rows


def write_csv(path: Path, fields: Iterable[str], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields))
        writer.writeheader()
        writer.writerows(rows)


def copy_input(path: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(path, destination)


def finite(row: dict[str, str], field: str) -> float:
    raw = row.get(field, "")
    try:
        value = float(raw)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid numeric {field}: {raw!r}") from exc
    if not math.isfinite(value):
        raise ValueError(f"non-finite {field}: {value}")
    return value


def is_true(value: object) -> bool:
    return str(value).strip().lower() == "true"


def coordinate_key(xi: float, mu: float) -> tuple[float, float]:
    return round(xi, 10), round(mu, 10)


def xi_key(value: float) -> float:
    return round(value, 10)


def source_record(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "sha256": sha256(path), "bytes": path.stat().st_size}


def infer_replay_run_id(path: Path) -> str | None:
    match = re.search(r"_(\d+)_replay$", path.name)
    return match.group(1) if match else None


def require_fields(rows: list[dict[str, str]], fields: set[str], label: str) -> None:
    missing = fields.difference(rows[0])
    if missing:
        raise ValueError(f"{label} schema missing: {sorted(missing)}")


def native_segments(
    axis,
    rows: list[dict[str, object]],
    *,
    coordinate: str,
    max_gap: float,
    color: str,
    alpha: float,
    linewidth: float,
    linestyle: str = "-",
    marker: str | None = None,
    marker_size: float = 8.0,
    marker_face: str | None = None,
) -> None:
    """Draw ordered native support and leave larger gaps blank."""
    grouped: defaultdict[float, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[xi_key(float(row["xi"]))].append(row)
    for xi, group in sorted(grouped.items()):
        ordered = sorted(group, key=lambda row: float(row[coordinate]))
        if marker is not None:
            axis.scatter(
                [float(row["mu_MeV"]) for row in ordered],
                [xi] * len(ordered),
                [float(row["T_MeV"]) for row in ordered],
                marker=marker,
                s=marker_size,
                facecolors=marker_face if marker_face is not None else color,
                edgecolors=color,
                linewidths=0.55,
                alpha=alpha,
                depthshade=False,
            )
        segment = [ordered[0]]
        for row in ordered[1:]:
            previous = segment[-1]
            if float(row[coordinate]) - float(previous[coordinate]) <= max_gap:
                segment.append(row)
                continue
            if len(segment) >= 2:
                axis.plot(
                    [float(item["mu_MeV"]) for item in segment],
                    [xi] * len(segment),
                    [float(item["T_MeV"]) for item in segment],
                    color=color,
                    alpha=alpha,
                    linewidth=linewidth,
                    linestyle=linestyle,
                )
            segment = [row]
        if len(segment) >= 2:
            axis.plot(
                [float(item["mu_MeV"]) for item in segment],
                [xi] * len(segment),
                [float(item["T_MeV"]) for item in segment],
                color=color,
                alpha=alpha,
                linewidth=linewidth,
                linestyle=linestyle,
            )


def main() -> int:
    args = parse_args()
    v5_root = args.v5_root.resolve()
    expansion_root = args.expansion_root.resolve()
    output_root = args.output_root.resolve()
    figures = output_root / "figures"
    tables = output_root / "tables"
    figures.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)

    v5_manifest_path = v5_root / "manifest.json"
    v5_manifest = json.loads(v5_manifest_path.read_text(encoding="utf-8"))
    if v5_manifest.get("schema_version") != "pnjl_c2_phase_surfaces_diagnostic_v5":
        raise ValueError("v5-root does not contain the expected v5 manifest")
    calculation_sha = str(v5_manifest.get("calculation_sha", ""))
    if len(calculation_sha) != 40:
        raise ValueError("v5 calculation SHA is not a 40-character commit")

    expansion_manifest_path = expansion_root / "manifest.json"
    expansion_summary_path = expansion_root / "expansion_summary.csv"
    expansion_verdict_path = expansion_root / "verdict.json"
    expansion_manifest = json.loads(expansion_manifest_path.read_text(encoding="utf-8"))
    expansion_verdict = json.loads(expansion_verdict_path.read_text(encoding="utf-8"))
    replay_run_id = infer_replay_run_id(expansion_root)
    if expansion_manifest.get("schema_version") != "pnjl_issue130_crossover_mu_endpoint_expansion_v1":
        raise ValueError("unexpected endpoint expansion schema")
    if expansion_manifest.get("calculation_sha") != calculation_sha:
        raise ValueError("endpoint expansion calculation SHA differs from v5 baseline")
    if expansion_manifest.get("run_mode") != "aggregate_replay":
        raise ValueError("v6 requires the solver-free aggregate replay manifest")
    if expansion_manifest.get("solver_called") is not False or expansion_manifest.get("reference_write") is not False:
        raise ValueError("endpoint replay provenance is not solver-free/reference-free")
    if expansion_manifest.get("oracle_labels_consumed") is not False:
        raise ValueError("endpoint replay consumed oracle labels")
    if expansion_manifest.get("verdict") != "expansion_candidate" or expansion_verdict.get("verdict") != "expansion_candidate":
        raise ValueError("endpoint expansion is not an accepted diagnostic candidate")

    v5_tables = {
        "maxwell": v5_root / "tables" / "maxwell_surface_point_status.csv",
        "maxwell_status": v5_root / "tables" / "maxwell_grid_status_counts.csv",
        "grid_unresolved": v5_root / "tables" / "grid_unresolved_diagnostics.csv",
        "crossover_filter": v5_root / "tables" / "crossover_physical_filter.csv",
        "crossover_gap": v5_root / "tables" / "crossover_cep_sampling_gap.csv",
        "claim_ledger": v5_root / "tables" / "claim_ledger.csv",
    }
    for path in v5_tables.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    maxwell_rows = read_csv(v5_tables["maxwell"])
    crossover_filter_rows = read_csv(v5_tables["crossover_filter"])
    crossover_gap_rows = read_csv(v5_tables["crossover_gap"])
    require_fields(maxwell_rows, {"xi", "T_MeV", "mu_MeV", "grid_unresolved"}, "v5 Maxwell")
    require_fields(crossover_filter_rows, {"xi", "mu_MeV", "T_MeV", "rho", "physical_status"}, "v5 crossover")
    require_fields(crossover_gap_rows, {"xi", "muq_CEP_proxy_MeV", "T_CEP_bracket_low_MeV", "T_CEP_bracket_high_MeV"}, "v5 CEP")

    spinodals_path = args.spinodals.resolve() if args.spinodals else None
    if spinodals_path is None:
        manifest_spinodal = v5_manifest.get("inputs", {}).get("spinodals", {}).get("path")
        if manifest_spinodal:
            spinodals_path = Path(manifest_spinodal)
    if spinodals_path is None or not spinodals_path.is_file():
        raise FileNotFoundError("v5 spinodal input is not available; pass --spinodals explicitly")
    spinodal_rows = read_csv(spinodals_path)
    require_fields(
        spinodal_rows,
        {"xi", "T_MeV", "mu_spinodal_hadron_MeV", "mu_spinodal_quark_MeV"},
        "v5 spinodals",
    )

    expansion_rows = read_csv(expansion_summary_path)
    require_fields(
        expansion_rows,
        {
            "target_id",
            "xi",
            "target_mu_MeV",
            "mu_CEP_proxy_MeV",
            "T_crossover_MeV",
            "rho_crossover",
            "status",
            "finite_and_converged",
        },
        "endpoint expansion",
    )
    if len(expansion_rows) != int(expansion_manifest.get("materialized_target_count", len(expansion_rows))):
        raise ValueError("endpoint summary/materialized count mismatch")

    cep_by_xi: dict[float, dict[str, float]] = {}
    for row in crossover_gap_rows:
        xi = finite(row, "xi")
        k = xi_key(xi)
        if k in cep_by_xi:
            raise ValueError(f"duplicate v5 CEP xi: {xi}")
        cep_by_xi[k] = {
            "mu": finite(row, "muq_CEP_proxy_MeV"),
            "low": finite(row, "T_CEP_bracket_low_MeV"),
            "high": finite(row, "T_CEP_bracket_high_MeV"),
        }

    baseline_rows: list[dict[str, object]] = []
    baseline_keys: set[tuple[float, float]] = set()
    for row in crossover_filter_rows:
        if row.get("physical_status") != "retained_crossover":
            continue
        xi = finite(row, "xi")
        mu = finite(row, "mu_MeV")
        key = coordinate_key(xi, mu)
        if key in baseline_keys:
            raise ValueError(f"duplicate v5 crossover key: {key}")
        baseline_keys.add(key)
        if key[0] not in cep_by_xi or mu > cep_by_xi[key[0]]["mu"] + 1e-8:
            raise ValueError(f"v5 retained crossover violates CEP physical filter: {key}")
        baseline_rows.append(
            {
                "xi": xi,
                "mu_MeV": mu,
                "T_MeV": finite(row, "T_MeV"),
                "rho": finite(row, "rho"),
                "source_layer": "v5_baseline",
                "source_id": "v5_crossover_physical_filter.csv",
                "physical_status": "retained_crossover",
                "plot_status": "retained",
                "mu_CEP_proxy_MeV": cep_by_xi[key[0]]["mu"],
            }
        )

    endpoint_rows: list[dict[str, object]] = []
    excluded_endpoint_rows: list[dict[str, object]] = []
    endpoint_keys: set[tuple[float, float]] = set()
    for row in expansion_rows:
        xi = finite(row, "xi")
        mu = finite(row, "target_mu_MeV")
        key = coordinate_key(xi, mu)
        target_id = row["target_id"]
        proxy = finite(row, "mu_CEP_proxy_MeV")
        if key[0] not in cep_by_xi:
            raise ValueError(f"endpoint xi missing from v5 CEP table: {xi}")
        if abs(proxy - cep_by_xi[key[0]]["mu"]) > 1e-8:
            raise ValueError(f"endpoint/v5 CEP proxy mismatch at xi={xi}")
        status = "eligible"
        reason = "finite_converged_candidate_below_CEP_proxy"
        if not is_true(row["finite_and_converged"]) or row["status"] != "crossover_candidate":
            status = "excluded"
            reason = "endpoint_not_finite_converged_candidate"
        elif mu > proxy + 1e-8:
            status = "excluded"
            reason = "target_mu_above_CEP_proxy"
        elif key in baseline_keys:
            status = "excluded"
            reason = "duplicate_v5_baseline_coordinate"
        elif key in endpoint_keys:
            raise ValueError(f"duplicate endpoint expansion key: {key}")
        endpoint_record = {
            "xi": xi,
            "mu_MeV": mu,
            "T_MeV": finite(row, "T_crossover_MeV"),
            "rho": finite(row, "rho_crossover"),
            "source_layer": "endpoint_expansion",
            "source_id": target_id,
            "physical_status": "retained_crossover" if status == "eligible" else "excluded",
            "plot_status": "retained" if status == "eligible" else "excluded",
            "mu_CEP_proxy_MeV": proxy,
            "endpoint_status": row["status"],
            "finite_and_converged": row["finite_and_converged"],
            "exclusion_reason": "" if status == "eligible" else reason,
        }
        if status == "eligible":
            endpoint_keys.add(key)
            endpoint_rows.append(endpoint_record)
        else:
            excluded_endpoint_rows.append(endpoint_record)

    overlay_rows = baseline_rows + endpoint_rows
    overlay_rows.sort(key=lambda row: (float(row["xi"]), float(row["mu_MeV"]), str(row["source_layer"])))
    overlay_fields = (
        "xi",
        "mu_MeV",
        "T_MeV",
        "rho",
        "source_layer",
        "source_id",
        "physical_status",
        "plot_status",
        "mu_CEP_proxy_MeV",
        "endpoint_status",
        "finite_and_converged",
        "exclusion_reason",
    )
    write_csv(tables / "crossover_overlay.csv", overlay_fields, overlay_rows)
    write_csv(tables / "crossover_endpoint_expansion_summary.csv", expansion_rows[0].keys(), expansion_rows)
    write_csv(tables / "crossover_endpoint_excluded.csv", overlay_fields, excluded_endpoint_rows)

    overlay_gap_rows: list[dict[str, object]] = []
    for xi in sorted({xi_key(float(row["xi"])) for row in overlay_rows}):
        base = sorted((row for row in baseline_rows if xi_key(float(row["xi"])) == xi), key=lambda row: float(row["mu_MeV"]))
        endpoint = sorted((row for row in endpoint_rows if xi_key(float(row["xi"])) == xi), key=lambda row: float(row["mu_MeV"]))
        first_excluded = next((row for row in crossover_filter_rows if xi_key(finite(row, "xi")) == xi and row.get("physical_status") == "excluded_response_peak_above_CEP"), None)
        base_last = base[-1] if base else None
        ep_first = endpoint[0] if endpoint else None
        ep_last = endpoint[-1] if endpoint else None
        overlay_gap_rows.append(
            {
                "xi": xi,
                "v5_baseline_count": len(base),
                "endpoint_count": len(endpoint),
                "v5_last_mu_MeV": float(base_last["mu_MeV"]) if base_last else "",
                "endpoint_first_mu_MeV": float(ep_first["mu_MeV"]) if ep_first else "",
                "endpoint_last_mu_MeV": float(ep_last["mu_MeV"]) if ep_last else "",
                "gap_v5_last_to_endpoint_first_MeV": (float(ep_first["mu_MeV"]) - float(base_last["mu_MeV"])) if base_last and ep_first else "",
                "v5_first_excluded_mu_MeV": finite(first_excluded, "mu_MeV") if first_excluded else "",
                "cep_proxy_mu_MeV": cep_by_xi[xi]["mu"],
                "endpoint_gap_to_CEP_proxy_MeV": (cep_by_xi[xi]["mu"] - float(ep_last["mu_MeV"])) if ep_last else "",
                "endpoint_proxy_relation": "below_CEP_proxy" if ep_last else "no_endpoint_candidate",
                "interpretation": "overlay extends native crossover support; no implicit connector",
            }
        )
    write_csv(tables / "crossover_overlay_gap.csv", overlay_gap_rows[0].keys(), overlay_gap_rows)
    write_csv(
        tables / "crossover_source_counts.csv",
        ("source_layer", "row_count", "interpretation"),
        [
            {"source_layer": "v5_baseline", "row_count": len(baseline_rows), "interpretation": "copied from v5 retained crossover support"},
            {"source_layer": "endpoint_expansion", "row_count": len(endpoint_rows), "interpretation": "new finite/converged numerical endpoint candidates"},
            {"source_layer": "endpoint_excluded", "row_count": len(excluded_endpoint_rows), "interpretation": "retained in audit only; not plotted"},
        ],
    )

    # Preserve v5 tables byte-for-byte in the v6 package for local inspection.
    copied_inputs = {
        "v5_source_manifest.json": v5_manifest_path,
        "v5_maxwell_surface_point_status.csv": v5_tables["maxwell"],
        "v5_maxwell_grid_status_counts.csv": v5_tables["maxwell_status"],
        "v5_grid_unresolved_diagnostics.csv": v5_tables["grid_unresolved"],
        "v5_crossover_physical_filter.csv": v5_tables["crossover_filter"],
        "v5_crossover_cep_sampling_gap.csv": v5_tables["crossover_gap"],
        "v5_claim_ledger.csv": v5_tables["claim_ledger"],
        "v5_spinodals.csv": spinodals_path,
        "endpoint_replay_manifest.json": expansion_manifest_path,
    }
    for name, source in copied_inputs.items():
        copy_input(source, tables / name)

    # Normalize Maxwell and spinodal rows only for plotting.  The copied v5
    # files above remain the authoritative retained evidence.
    maxwell_plot: list[dict[str, object]] = []
    for row in maxwell_rows:
        maxwell_plot.append(
            {
                "xi": finite(row, "xi"),
                "T_MeV": finite(row, "T_MeV"),
                "mu_MeV": finite(row, "mu_MeV"),
                "rho": finite(row, "rho_hadron") if row.get("rho_hadron", "") else 0.0,
            }
        )
    spinodal_plot: list[dict[str, object]] = []
    for row in spinodal_rows:
        xi = finite(row, "xi")
        T = finite(row, "T_MeV")
        for branch, field in (("hadron", "mu_spinodal_hadron_MeV"), ("quark", "mu_spinodal_quark_MeV")):
            spinodal_plot.append({"xi": xi, "T_MeV": T, "mu_MeV": finite(row, field), "branch": branch})

    figure = plt.figure(figsize=(13.2, 9.0))
    axis = figure.add_subplot(111, projection="3d")
    native_segments(
        axis,
        maxwell_plot,
        coordinate="T_MeV",
        max_gap=args.max_temperature_gap,
        color="#1f5a91",
        alpha=0.72,
        linewidth=0.48,
    )
    for branch, color in (("hadron", "#526b82"), ("quark", "#7c8c9b")):
        native_segments(
            axis,
            [row for row in spinodal_plot if row["branch"] == branch],
            coordinate="T_MeV",
            max_gap=args.max_temperature_gap,
            color=color,
            alpha=0.42,
            linewidth=0.32,
            linestyle=":",
        )
    native_segments(
        axis,
        baseline_rows,
        coordinate="mu_MeV",
        max_gap=args.max_mu_gap,
        color="#d9792b",
        alpha=0.38,
        linewidth=0.48,
    )
    native_segments(
        axis,
        endpoint_rows,
        coordinate="mu_MeV",
        max_gap=args.max_mu_gap,
        color="#d9792b",
        alpha=0.95,
        linewidth=1.05,
        marker="o",
        marker_size=13,
        marker_face="none",
    )

    for row in crossover_gap_rows:
        xi = finite(row, "xi")
        mu = finite(row, "muq_CEP_proxy_MeV")
        low = finite(row, "T_CEP_bracket_low_MeV")
        high = finite(row, "T_CEP_bracket_high_MeV")
        lo, hi = min(low, high), max(low, high)
        axis.plot([mu, mu], [xi, xi], [lo, hi], color="#333333", linestyle="--", linewidth=0.55, alpha=0.48)
        axis.scatter([mu], [xi], [0.5 * (lo + hi)], color="#333333", marker="D", s=8, facecolors="none", depthshade=False)

    legend = [
        Line2D([0], [0], color="#1f5a91", lw=2.0, alpha=0.72, label="Maxwell v5 native support"),
        Line2D([0], [0], color="#526b82", lw=1.2, linestyle=":", alpha=0.55, label="spinodal v5 input support"),
        Line2D([0], [0], color="#d9792b", lw=2.0, alpha=0.38, label="crossover v5 baseline"),
        Line2D([0], [0], color="#d9792b", lw=2.0, marker="o", markerfacecolor="none", label="crossover endpoint expansion"),
        Line2D([0], [0], color="#333333", lw=1.0, linestyle="--", marker="D", markerfacecolor="none", label="CEP bracket / midpoint"),
    ]
    axis.legend(handles=legend, loc="upper left", bbox_to_anchor=(0.01, 0.98), framealpha=0.92, fontsize=8)
    axis.set_xlabel(r"$\mu_q$ [MeV]", labelpad=10)
    axis.set_ylabel(r"$\xi$", labelpad=10)
    axis.set_zlabel(r"$T$ [MeV]", labelpad=10)
    axis.set_title(r"C2 v6 phase curves: v5 baseline + crossover endpoint overlay", pad=12)
    axis.view_init(elev=25, azim=-62)
    try:
        axis.set_box_aspect((1.55, 1.0, 1.25))
    except AttributeError:
        pass
    axis.text2D(
        0.01,
        0.01,
        f"diagnostic-only | Maxwell(v5)={len(maxwell_plot)} | spinodal(v5)={len(spinodal_rows)} | crossover(v5)={len(baseline_rows)} + endpoint={len(endpoint_rows)} | CEP brackets={len(crossover_gap_rows)} | no interpolation/triangulation",
        transform=figure.transFigure,
        fontsize=7,
        color="#444444",
    )
    figure.tight_layout()
    figure_path = figures / "c2_phase_surfaces_mu_xi_T_v6_crossover_overlay.png"
    figure.savefig(figure_path, dpi=220, bbox_inches="tight")
    plt.close(figure)

    output_inputs = {
        "v5_manifest": v5_manifest_path,
        "v5_maxwell": v5_tables["maxwell"],
        "v5_maxwell_status": v5_tables["maxwell_status"],
        "v5_grid_unresolved": v5_tables["grid_unresolved"],
        "v5_crossover_filter": v5_tables["crossover_filter"],
        "v5_crossover_gap": v5_tables["crossover_gap"],
        "v5_claim_ledger": v5_tables["claim_ledger"],
        "v5_spinodals": spinodals_path,
        "endpoint_manifest": expansion_manifest_path,
        "endpoint_summary": expansion_summary_path,
        "endpoint_verdict": expansion_verdict_path,
    }
    script_path = Path(__file__).resolve()
    plot_manifest = {
        "schema_version": SCHEMA,
        "figure_mode": "diagnostic_no_interpolation_no_triangulation",
        "orientation": "x=mu_q, y=xi, z=T",
        "calculation_sha": calculation_sha,
        "solver_called": False,
        "reference_write": False,
        "generator": {"script": str(script_path), "script_sha256": sha256(script_path)},
        "source_layers": {
            "v5_baseline": "all Maxwell, spinodal, CEP and retained crossover processing is copied from v5 evidence",
            "endpoint_expansion": "186 aggregate-replay-validated candidates appended only to crossover support",
        },
        "endpoint_expansion_provenance": {
            "source_run_id": expansion_manifest.get("source_run_id"),
            "replay_run_id": replay_run_id,
            "postprocess_sha": expansion_manifest.get("postprocess_sha"),
        },
        "surface_policy": {
            "maxwell": "v5 maxwell_surface_point_status copied unchanged; no C2 raw recomputation",
            "spinodal": "v5 spinodal input copied unchanged and drawn as native dotted support",
            "crossover_baseline": "v5 retained crossover rows copied unchanged",
            "crossover_overlay": "finite/converged endpoint candidates with target_mu <= the same v5 CEP proxy",
            "cep": "v5 bracket projections and midpoint markers retained; midpoint is not a certified CEP",
            "interpolation": "disabled; native ordered segments only",
            "triangulation": "disabled",
            "gap_policy": "segments are not connected when their native coordinate gap exceeds the explicit limit",
        },
        "limits": {"max_mu_gap_MeV": args.max_mu_gap, "max_temperature_gap_MeV": args.max_temperature_gap},
        "inputs": {name: source_record(path) for name, path in output_inputs.items()},
        "counts": {
            "maxwell_rows": len(maxwell_plot),
            "spinodal_rows": len(spinodal_rows),
            "v5_crossover_rows": len(baseline_rows),
            "endpoint_rows": len(endpoint_rows),
            "endpoint_excluded_rows": len(excluded_endpoint_rows),
            "cep_brackets": len(crossover_gap_rows),
        },
        "output": {"path": str(figure_path), "sha256": sha256(figure_path)},
    }
    (figures / "plot_manifest.json").write_text(json.dumps(plot_manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    claim_rows = [
        {"claim_id": "v5_baseline_preserved", "claim": "Maxwell, spinodal, CEP and retained crossover evidence are inherited from the v5 postprocessed package", "status": "supported", "evidence": "tables/v5_source_manifest.json; tables/v5_maxwell_surface_point_status.csv; tables/v5_spinodals.csv; tables/v5_crossover_cep_sampling_gap.csv", "boundary": "no raw C2 recomputation"},
        {"claim_id": "endpoint_overlay_materialized", "claim": "The 186 endpoint expansion candidates are appended as a separate crossover source layer", "status": "supported", "evidence": "tables/crossover_endpoint_expansion_summary.csv; tables/crossover_overlay.csv", "boundary": "diagnostic-only aggregate replay provenance"},
        {"claim_id": "physical_crossover_filter_preserved", "claim": "No response peak with mu_q above the per-xi CEP proxy is plotted as crossover", "status": "supported", "evidence": "tables/v5_crossover_physical_filter.csv; tables/crossover_endpoint_excluded.csv", "boundary": "CEP proxy remains bracket-derived"},
        {"claim_id": "native_gaps_preserved", "claim": "v5 and endpoint support are connected only across native gaps within the explicit plotting limit", "status": "supported", "evidence": "tables/crossover_overlay_gap.csv; figures/plot_manifest.json", "boundary": "no interpolation or triangulation"},
        {"claim_id": "phase_reference_promotion", "claim": "The v6 figure authorizes phase-reference promotion", "status": "blocked", "evidence": "decision.json", "boundary": "this is a diagnostic overlay, not a production/reference artifact"},
    ]
    write_csv(tables / "claim_ledger.csv", claim_rows[0].keys(), claim_rows)

    decision = {
        "schema_version": SCHEMA,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "calculation_sha": calculation_sha,
        "solver_called": False,
        "reference_write": False,
        "verdict": "diagnostic_v6_overlay_ready_for_author_review" if not excluded_endpoint_rows else "diagnostic_v6_overlay_with_excluded_endpoint_rows",
        "promotion_effect": "none",
        "baseline": {"source": "v5_postprocessed_package", "maxwell_rows": len(maxwell_plot), "spinodal_rows": len(spinodal_rows), "crossover_rows": len(baseline_rows), "cep_brackets": len(crossover_gap_rows)},
        "endpoint_expansion": {"source_run_id": expansion_manifest.get("source_run_id"), "replay_run_id": replay_run_id, "replay_artifact_root": str(expansion_root), "postprocess_sha": expansion_manifest.get("postprocess_sha"), "rows": len(endpoint_rows), "excluded_rows": len(excluded_endpoint_rows), "target_xi_count": len({xi_key(finite(row, 'xi')) for row in expansion_rows})},
        "policy": "v5 baseline plus endpoint crossover overlay; no Maxwell extension in this artifact",
    }
    (output_root / "decision.json").write_text(json.dumps(decision, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    (output_root / "README.md").write_text(
        "# C2 phase surfaces v6: v5 baseline + crossover endpoint overlay\n\n"
        "本包明确以 v5 后处理结果为基线，而不是直接从 C2 原始 CSV 重建。v5 的 Maxwell boundary、spinodal、CEP bracket、物理 crossover 筛选、unresolved 诊断和无三角化规则均保留；本包只追加 endpoint expansion 的 crossover 结果。\n\n"
        f"固定 calculation SHA 为 `{calculation_sha}`。endpoint expansion 来自 numerical run `{expansion_manifest.get('source_run_id')}` 的 aggregate replay，replay run provenance 位于 `tables/endpoint_replay_manifest.json`；本地物化过程 `solver_called=false`、`reference_write=false`、`oracle_labels_consumed=false`。\n\n"
        f"v5 基线包含 Maxwell `{len(maxwell_plot)}` 行、spinodal `{len(spinodal_rows)}` 行、保留 crossover `{len(baseline_rows)}` 行和 `{len(crossover_gap_rows)}` 个 CEP bracket。新增 overlay 为 `{len(endpoint_rows)}` 行，覆盖 `{len({xi_key(finite(row, 'xi')) for row in expansion_rows})}` 个非均匀 xi 切片；它不是完整均匀 xi 网格，也不代表 Maxwell 区已补全。\n\n"
        "图中不使用插值或三角化；超过显式 gap 上限的 native support 保持空白。endpoint 点以空心圆标识，便于审核其相对 v5 端点的延伸，但该标识不表示新的物理证书。`mu_q > mu_CEP` 的响应峰仍不绘制为 crossover。\n\n"
        "本包 verdict 为 diagnostic-only，不能触发 phase-reference 晋升、正式 production 或 transport。Maxwell 区数值扩展仍是独立任务。\n",
        encoding="utf-8",
    )
    (output_root / "AUDIT.md").write_text(
        "## Evidence boundary\n\n"
        "- v5 tables are copied byte-for-byte into `tables/v5_*` and remain the baseline.\n"
        "- Endpoint expansion summary is copied from the solver-free replay artifact; the full response curves remain in the external artifact referenced by the replay manifest.\n"
        "- The plotting script does not call the PNJL solver and does not write a reference.\n"
        "- The figure is a visual diagnostic. It does not convert unresolved v5 geometry into a certificate and does not infer missing xi slices.\n",
        encoding="utf-8",
    )

    output_files: dict[str, str] = {}
    for path in sorted(output_root.rglob("*")):
        if path.is_file() and path.name != "manifest.json":
            output_files[path.relative_to(output_root).as_posix()] = sha256(path)
    manifest = {
        "schema_version": SCHEMA,
        "generated_at_utc": decision["generated_at_utc"],
        "calculation_sha": calculation_sha,
        "solver_called": False,
        "reference_write": False,
        "generator": {"script": str(script_path), "script_sha256": sha256(script_path)},
        "verdict": decision["verdict"],
        "source_manifests": {"v5": source_record(v5_manifest_path), "endpoint_replay": source_record(expansion_manifest_path)},
        "counts": decision["baseline"] | {"endpoint_rows": len(endpoint_rows), "endpoint_excluded_rows": len(excluded_endpoint_rows)},
        "output_files": output_files,
    }
    (output_root / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"verdict": decision["verdict"], "output_root": str(output_root), "counts": manifest["counts"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
