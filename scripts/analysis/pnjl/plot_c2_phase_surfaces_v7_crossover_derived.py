#!/usr/bin/env python3
"""Build the v7 solver-free crossover-derived phase-surface package.

The v6 package is immutable input.  This script materializes two explicitly
non-certified derived layers:

* piecewise-linear points inside each native ``xi`` slice;
* midpoint ``xi`` slices, but only over the common native crossover support of
  the two neighbouring slices.

The selected CEP boundary estimate is allowed as a physical crossover
endpoint.  It is recorded as a boundary-constrained interpolation and is not
represented as a new strict CEP solve.  Maxwell and spinodal rows are copied
unchanged from v6; no Maxwell values are interpolated here.

No equilibrium solver is called and no reference artifact is written.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


SCHEMA = "pnjl_c2_phase_surfaces_diagnostic_v7_crossover_derived"
DEFAULT_SUBDIVISIONS = 4
DEFAULT_MAX_XI_GAP = 0.025


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v6-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--mu-subdivisions", type=int, default=DEFAULT_SUBDIVISIONS)
    parser.add_argument("--max-xi-gap", type=float, default=DEFAULT_MAX_XI_GAP)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_record(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "sha256": sha256(path), "bytes": path.stat().st_size}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"CSV has no rows: {path}")
    return rows


def write_csv(path: Path, fields: Iterable[str], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def finite(row: dict[str, object], field: str) -> float:
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


def xi_key(value: float) -> float:
    return round(value, 10)


def coordinate_key(xi: float, mu: float) -> tuple[float, float]:
    return xi_key(xi), round(mu, 10)


def source_key(xi: float, mu: float) -> str:
    return f"xi={xi:.10g};mu_MeV={mu:.12g}"


def lerp(x0: float, y0: float, x1: float, y1: float, x: float) -> float:
    if math.isclose(x0, x1, rel_tol=0.0, abs_tol=1e-12):
        return y0
    weight = (x - x0) / (x1 - x0)
    return y0 + weight * (y1 - y0)


def bracket_value(rows: list[dict[str, object]], mu: float) -> tuple[float, float, float, float]:
    """Interpolate T/rho and return the two source rows' mu values."""
    ordered = sorted(rows, key=lambda row: float(row["mu_MeV"]))
    if mu < float(ordered[0]["mu_MeV"]) - 1e-8 or mu > float(ordered[-1]["mu_MeV"]) + 1e-8:
        raise ValueError(f"mu={mu} is outside native support")
    for left, right in zip(ordered, ordered[1:]):
        mu_left = float(left["mu_MeV"])
        mu_right = float(right["mu_MeV"])
        if mu <= mu_right + 1e-10:
            return (
                lerp(mu_left, float(left["T_MeV"]), mu_right, float(right["T_MeV"]), mu),
                lerp(mu_left, float(left["rho"]), mu_right, float(right["rho"]), mu),
                mu_left,
                mu_right,
            )
    last = ordered[-1]
    return float(last["T_MeV"]), float(last["rho"]), float(last["mu_MeV"]), float(last["mu_MeV"])


COMMON_FIELDS = (
    "surface",
    "xi",
    "mu_MeV",
    "T_MeV",
    "rho",
    "layer",
    "status",
    "interpolation_method",
    "source_key_left",
    "source_key_right",
    "source_keys",
    "distance_to_source_MeV",
    "mu_CEP_proxy_MeV",
    "T_CEP_bracket_low_MeV",
    "T_CEP_bracket_high_MeV",
    "boundary_relation",
    "physical_region",
    "error_estimate",
    "source_layer",
    "calculation_sha",
    "v6_manifest_hash",
)


def make_row(
    *,
    xi: float,
    mu: float,
    T: float,
    rho: float,
    layer: str,
    status: str,
    method: str,
    source_left: str,
    source_right: str,
    source_keys: str,
    distance: float,
    cep_mu: float,
    cep_low: float,
    cep_high: float,
    boundary_relation: str,
    physical_region: str,
    source_layer: str,
    calculation_sha: str,
    v6_manifest_hash: str,
) -> dict[str, object]:
    values = {
        "surface": "crossover",
        "xi": xi,
        "mu_MeV": mu,
        "T_MeV": T,
        "rho": rho,
        "layer": layer,
        "status": status,
        "interpolation_method": method,
        "source_key_left": source_left,
        "source_key_right": source_right,
        "source_keys": source_keys,
        "distance_to_source_MeV": distance,
        "mu_CEP_proxy_MeV": cep_mu,
        "T_CEP_bracket_low_MeV": cep_low,
        "T_CEP_bracket_high_MeV": cep_high,
        "boundary_relation": boundary_relation,
        "physical_region": physical_region,
        "error_estimate": "not_estimated_linear_interpolation",
        "source_layer": source_layer,
        "calculation_sha": calculation_sha,
        "v6_manifest_hash": v6_manifest_hash,
    }
    for field in ("xi", "mu_MeV", "T_MeV", "rho", "distance_to_source_MeV", "mu_CEP_proxy_MeV", "T_CEP_bracket_low_MeV", "T_CEP_bracket_high_MeV"):
        if not math.isfinite(float(values[field])):
            raise ValueError(f"non-finite derived field {field}")
    return values


def plot_segments(axis, rows: list[dict[str, object]], *, color: str, alpha: float, linewidth: float, linestyle: str) -> None:
    grouped: defaultdict[float, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[xi_key(float(row["xi"]))].append(row)
    for xi, group in sorted(grouped.items()):
        ordered = sorted(group, key=lambda row: float(row["mu_MeV"]))
        if len(ordered) < 2:
            continue
        axis.plot(
            [float(row["mu_MeV"]) for row in ordered],
            [xi] * len(ordered),
            [float(row["T_MeV"]) for row in ordered],
            color=color,
            alpha=alpha,
            linewidth=linewidth,
            linestyle=linestyle,
        )


def make_figure(
    path: Path,
    *,
    maxwell: list[dict[str, object]],
    spinodal: list[dict[str, object]],
    certified: list[dict[str, object]],
    same_xi: list[dict[str, object]],
    xi_derived: list[dict[str, object]],
    boundary: list[dict[str, object]],
    layer_audit: bool,
) -> None:
    figure = plt.figure(figsize=(13.6, 9.0))
    axis = figure.add_subplot(111, projection="3d")
    plot_segments(axis, maxwell, color="#1f5a91", alpha=0.78, linewidth=0.5, linestyle="-")
    plot_segments(axis, spinodal, color="#6f7f8c", alpha=0.30, linewidth=0.28, linestyle=":")
    plot_segments(axis, certified, color="#d9792b", alpha=0.90, linewidth=0.65, linestyle="-")
    plot_segments(axis, same_xi, color="#d9792b", alpha=0.36 if layer_audit else 0.52, linewidth=0.45, linestyle="--")
    plot_segments(axis, xi_derived, color="#a95422" if layer_audit else "#d9792b", alpha=0.60 if layer_audit else 0.42, linewidth=0.42, linestyle="-." if layer_audit else "--")
    if boundary:
        axis.scatter(
            [float(row["mu_MeV"]) for row in boundary],
            [float(row["xi"]) for row in boundary],
            [float(row["T_MeV"]) for row in boundary],
            marker="D",
            s=11,
            facecolors="none",
            edgecolors="#222222",
            linewidths=0.45,
            alpha=0.75,
            depthshade=False,
        )
    legend = [
        Line2D([0], [0], color="#1f5a91", lw=2, label="Maxwell v6 native support"),
        Line2D([0], [0], color="#6f7f8c", lw=1, linestyle=":", alpha=0.65, label="spinodal v6 native support"),
        Line2D([0], [0], color="#d9792b", lw=2, label="crossover native support"),
        Line2D([0], [0], color="#d9792b", lw=1.5, linestyle="--", alpha=0.48, label="same-xi derived"),
        Line2D([0], [0], color="#a95422" if layer_audit else "#d9792b", lw=1.5, linestyle="-." if layer_audit else "--", alpha=0.58, label="xi-derived common-support"),
        Line2D([0], [0], color="#222222", marker="D", markerfacecolor="none", lw=0, label="CEP boundary estimate"),
    ]
    axis.legend(handles=legend, loc="upper left", bbox_to_anchor=(0.01, 0.98), framealpha=0.92, fontsize=8)
    axis.set_xlabel(r"$\mu_q$ [MeV]", labelpad=10)
    axis.set_ylabel(r"$\xi$", labelpad=10)
    axis.set_zlabel(r"$T$ [MeV]", labelpad=10)
    mode = "layer audit" if layer_audit else "complete derived display"
    axis.set_title(f"C2 v7 crossover-derived phase surfaces ({mode})", pad=12)
    axis.view_init(elev=25, azim=-62)
    try:
        axis.set_box_aspect((1.55, 1.0, 1.25))
    except AttributeError:
        pass
    axis.text2D(
        0.01,
        0.01,
        f"diagnostic-only | native Maxwell={len(maxwell)} | native crossover={len(certified)} | same-xi derived={len(same_xi)} | xi-derived={len(xi_derived)} | boundary={len(boundary)} | no triangulation",
        transform=figure.transFigure,
        fontsize=7,
        color="#444444",
    )
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def main() -> int:
    args = parse_args()
    if args.mu_subdivisions < 2:
        raise ValueError("--mu-subdivisions must be at least 2")
    if args.max_xi_gap <= 0:
        raise ValueError("--max-xi-gap must be positive")
    v6_root = args.v6_root.resolve()
    output_root = args.output_root.resolve()
    v6_manifest_path = v6_root / "manifest.json"
    if not v6_manifest_path.is_file():
        raise FileNotFoundError(v6_manifest_path)
    v6_manifest = json.loads(v6_manifest_path.read_text(encoding="utf-8"))
    if v6_manifest.get("schema_version") != "pnjl_c2_phase_surfaces_diagnostic_v6_crossover_overlay":
        raise ValueError("v6-root does not contain the expected v6 manifest")
    calculation_sha = str(v6_manifest.get("calculation_sha", ""))
    if len(calculation_sha) != 40 or any(char not in "0123456789abcdef" for char in calculation_sha.lower()):
        raise ValueError("v6 calculation SHA is not a 40-character hexadecimal commit")
    v6_manifest_hash = sha256(v6_manifest_path)

    v6_tables = v6_root / "tables"
    crossover_path = v6_tables / "crossover_overlay.csv"
    gap_path = v6_tables / "crossover_overlay_gap.csv"
    maxwell_path = v6_tables / "v5_maxwell_surface_point_status.csv"
    spinodal_path = v6_tables / "v5_spinodals.csv"
    for path in (crossover_path, gap_path, maxwell_path, spinodal_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    crossover_rows = read_csv(crossover_path)
    gap_rows = read_csv(gap_path)
    maxwell_rows = read_csv(maxwell_path)
    spinodal_rows = read_csv(spinodal_path)

    cep_by_xi: dict[float, dict[str, float]] = {}
    for row in gap_rows:
        xi = finite(row, "xi")
        key = xi_key(xi)
        if key in cep_by_xi:
            raise ValueError(f"duplicate CEP xi in v6 gap table: {xi}")
        cep_by_xi[key] = {
            "mu": finite(row, "cep_proxy_mu_MeV"),
            "low": math.nan,
            "high": math.nan,
        }
    # v6 overlay_gap carries the endpoint geometry but the bracket values live
    # in the copied v5 gap table.  Join them explicitly instead of guessing a
    # bracket from the endpoint points.
    v5_gap_path = v6_tables / "v5_crossover_cep_sampling_gap.csv"
    v5_gap_rows = read_csv(v5_gap_path)
    for row in v5_gap_rows:
        key = xi_key(finite(row, "xi"))
        if key not in cep_by_xi:
            raise ValueError(f"v5 CEP xi missing from v6 overlay gap: {key}")
        cep_by_xi[key]["low"] = finite(row, "T_CEP_bracket_low_MeV")
        cep_by_xi[key]["high"] = finite(row, "T_CEP_bracket_high_MeV")

    native_by_xi: defaultdict[float, list[dict[str, object]]] = defaultdict(list)
    native_keys: set[tuple[float, float]] = set()
    for row in crossover_rows:
        if row.get("physical_status") != "retained_crossover" or row.get("plot_status") != "retained":
            continue
        xi = finite(row, "xi")
        mu = finite(row, "mu_MeV")
        key = coordinate_key(xi, mu)
        if key in native_keys:
            raise ValueError(f"duplicate native crossover key: {key}")
        native_keys.add(key)
        if xi_key(xi) not in cep_by_xi:
            raise ValueError(f"native crossover xi has no CEP boundary: {xi}")
        cep = cep_by_xi[xi_key(xi)]
        if mu > cep["mu"] + 1e-8:
            raise ValueError(f"native crossover crosses CEP boundary: {key}")
        native_by_xi[xi_key(xi)].append({
            "xi": xi,
            "mu_MeV": mu,
            "T_MeV": finite(row, "T_MeV"),
            "rho": finite(row, "rho"),
            "source_layer": row.get("source_layer", "v6_native"),
        })
    if not native_by_xi:
        raise ValueError("v6 has no retained crossover support")

    native_rows: list[dict[str, object]] = []
    for xi, rows in sorted(native_by_xi.items()):
        cep = cep_by_xi[xi]
        for row in sorted(rows, key=lambda item: float(item["mu_MeV"])):
            mu = float(row["mu_MeV"])
            native_rows.append(make_row(
                xi=float(row["xi"]), mu=mu, T=float(row["T_MeV"]), rho=float(row["rho"]),
                layer="certified_layer", status="native_support", method="none_native_v6",
                source_left=source_key(float(row["xi"]), mu), source_right=source_key(float(row["xi"]), mu),
                source_keys=source_key(float(row["xi"]), mu), distance=0.0,
                cep_mu=cep["mu"], cep_low=cep["low"], cep_high=cep["high"],
                boundary_relation="below_CEP_boundary", physical_region="crossover_native_support",
                source_layer=str(row["source_layer"]), calculation_sha=calculation_sha, v6_manifest_hash=v6_manifest_hash,
            ))

    same_xi_rows: list[dict[str, object]] = []
    boundary_rows: list[dict[str, object]] = []
    for xi, rows in sorted(native_by_xi.items()):
        ordered = sorted(rows, key=lambda item: float(item["mu_MeV"]))
        cep = cep_by_xi[xi]
        for left, right in zip(ordered, ordered[1:]):
            mu_left = float(left["mu_MeV"])
            mu_right = float(right["mu_MeV"])
            for index in range(1, args.mu_subdivisions):
                mu = mu_left + (mu_right - mu_left) * index / args.mu_subdivisions
                T = lerp(mu_left, float(left["T_MeV"]), mu_right, float(right["T_MeV"]), mu)
                rho = lerp(mu_left, float(left["rho"]), mu_right, float(right["rho"]), mu)
                same_xi_rows.append(make_row(
                    xi=float(left["xi"]), mu=mu, T=T, rho=rho,
                    layer="interpolated_noncertified", status="derived_internal_mu",
                    method="piecewise_linear_same_xi", source_left=source_key(float(left["xi"]), mu_left),
                    source_right=source_key(float(right["xi"]), mu_right),
                    source_keys=f"{source_key(float(left['xi']), mu_left)}|{source_key(float(right['xi']), mu_right)}",
                    distance=min(mu - mu_left, mu_right - mu), cep_mu=cep["mu"], cep_low=cep["low"], cep_high=cep["high"],
                    boundary_relation="inside_CEP_boundary", physical_region="crossover_native_interval",
                    source_layer="v6_native_support", calculation_sha=calculation_sha, v6_manifest_hash=v6_manifest_hash,
                ))
        last = ordered[-1]
        previous = ordered[-2] if len(ordered) >= 2 else None
        if previous is not None and cep["mu"] > float(last["mu_MeV"]) + 1e-8:
            boundary_mu = cep["mu"]
            boundary_T = lerp(float(previous["mu_MeV"]), float(previous["T_MeV"]), float(last["mu_MeV"]), float(last["T_MeV"]), boundary_mu)
            boundary_rho = lerp(float(previous["mu_MeV"]), float(previous["rho"]), float(last["mu_MeV"]), float(last["rho"]), boundary_mu)
            boundary_rows.append(make_row(
                xi=float(last["xi"]), mu=boundary_mu, T=boundary_T, rho=boundary_rho,
                layer="boundary_constrained_endpoint_interpolated_noncertified",
                status="derived_cep_boundary", method="boundary_constrained_linear_extension",
                source_left=source_key(float(previous["xi"]), float(previous["mu_MeV"])),
                source_right=source_key(float(last["xi"]), float(last["mu_MeV"])),
                source_keys=f"{source_key(float(previous['xi']), float(previous['mu_MeV']))}|{source_key(float(last['xi']), float(last['mu_MeV']))}",
                distance=boundary_mu - float(last["mu_MeV"]), cep_mu=cep["mu"], cep_low=cep["low"], cep_high=cep["high"],
                boundary_relation="terminates_at_CEP_boundary_estimate", physical_region="crossover_CEP_boundary",
                source_layer="v6_native_support", calculation_sha=calculation_sha, v6_manifest_hash=v6_manifest_hash,
            ))

    xi_derived_rows: list[dict[str, object]] = []
    xi_mask_rows: list[dict[str, object]] = []
    xi_values = sorted(native_by_xi)
    for index, (xi_left, xi_right) in enumerate(zip(xi_values, xi_values[1:])):
        xi_gap = xi_right - xi_left
        left_rows = sorted(native_by_xi[xi_left], key=lambda item: float(item["mu_MeV"]))
        right_rows = sorted(native_by_xi[xi_right], key=lambda item: float(item["mu_MeV"]))
        left_cep = cep_by_xi[xi_left]
        right_cep = cep_by_xi[xi_right]
        common_upper = min(
            float(left_rows[-1]["mu_MeV"]), float(right_rows[-1]["mu_MeV"]),
            left_cep["mu"], right_cep["mu"],
        )
        common_lower = max(float(left_rows[0]["mu_MeV"]), float(right_rows[0]["mu_MeV"]))
        target_xi = 0.5 * (xi_left + xi_right)
        interval_points = {common_lower, common_upper}
        for rows in (left_rows, right_rows):
            interval_points.update(float(row["mu_MeV"]) for row in rows if common_lower <= float(row["mu_MeV"]) <= common_upper)
        ordered_points = sorted(interval_points)
        target_points: set[float] = set(ordered_points)
        for left_mu, right_mu in zip(ordered_points, ordered_points[1:]):
            for sub in range(1, args.mu_subdivisions):
                target_points.add(left_mu + (right_mu - left_mu) * sub / args.mu_subdivisions)
        if xi_gap <= args.max_xi_gap + 1e-10 and common_upper >= common_lower - 1e-8:
            for mu in sorted(target_points):
                left_T, left_rho, left_bracket_left, left_bracket_right = bracket_value(left_rows, mu)
                right_T, right_rho, right_bracket_left, right_bracket_right = bracket_value(right_rows, mu)
                T = 0.5 * (left_T + right_T)
                rho = 0.5 * (left_rho + right_rho)
                left_source = f"{source_key(xi_left, left_bracket_left)}|{source_key(xi_left, left_bracket_right)}"
                right_source = f"{source_key(xi_right, right_bracket_left)}|{source_key(xi_right, right_bracket_right)}"
                xi_derived_rows.append(make_row(
                    xi=target_xi, mu=mu, T=T, rho=rho,
                    layer="interpolated_noncertified", status="derived_internal_xi",
                    method="piecewise_linear_common_xi_support", source_left=left_source, source_right=right_source,
                    source_keys=f"{left_source}|{right_source}", distance=0.5 * xi_gap,
                    cep_mu=0.5 * (left_cep["mu"] + right_cep["mu"]),
                    cep_low=0.5 * (left_cep["low"] + right_cep["low"]), cep_high=0.5 * (left_cep["high"] + right_cep["high"]),
                    boundary_relation="inside_common_CEP_boundary", physical_region="crossover_common_native_support",
                    source_layer="v6_adjacent_xi_slices", calculation_sha=calculation_sha, v6_manifest_hash=v6_manifest_hash,
                ))
        xi_mask_rows.append({
            "xi": target_xi,
            "xi_role": "derived_midpoint",
            "xi_left": xi_left,
            "xi_right": xi_right,
            "xi_gap": xi_gap,
            "native_mu_low_left_MeV": float(left_rows[0]["mu_MeV"]),
            "native_mu_low_right_MeV": float(right_rows[0]["mu_MeV"]),
            "common_mu_high_MeV": common_upper,
            "cep_mu_left_MeV": left_cep["mu"],
            "cep_mu_right_MeV": right_cep["mu"],
            "materialized_rows": sum(1 for row in xi_derived_rows if xi_key(float(row["xi"])) == xi_key(target_xi)),
            "unresolved_tail_to_nearest_CEP_MeV": max(0.0, min(left_cep["mu"], right_cep["mu"]) - common_upper),
            "coverage_status": "common_support_only" if xi_gap <= args.max_xi_gap + 1e-10 else "unresolved_xi_gap_exceeds_policy",
        })

    # Native Maxwell/spinodal data are retained as a separate, unchanged
    # surface layer.  They are not interpolated by v7.
    maxwell_v7: list[dict[str, object]] = []
    maxwell_fields = ("xi", "T_MeV", "mu_MeV", "rho_hadron", "rho_quark", "area_residual", "grid_status", "grid_unresolved", "layer", "status", "calculation_sha", "v6_manifest_hash")
    for row in maxwell_rows:
        maxwell_v7.append({
            **{field: row.get(field, "") for field in maxwell_fields[:8]},
            "layer": "certified_layer",
            "status": "native_support_unresolved_geometry" if is_true(row.get("grid_unresolved", "")) else "native_support",
            "calculation_sha": calculation_sha,
            "v6_manifest_hash": v6_manifest_hash,
        })
    spinodal_v7: list[dict[str, object]] = []
    spinodal_plot: list[dict[str, object]] = []
    spinodal_fields = ("xi", "T_MeV", "mu_spinodal_hadron_MeV", "mu_spinodal_quark_MeV", "layer", "status", "calculation_sha", "v6_manifest_hash")
    for row in spinodal_rows:
        for field in spinodal_fields[:4]:
            finite(row, field)
        spinodal_v7.append({**{field: row[field] for field in spinodal_fields[:4]}, "layer": "certified_layer", "status": "native_support", "calculation_sha": calculation_sha, "v6_manifest_hash": v6_manifest_hash})
        for branch, field in (("hadron", "mu_spinodal_hadron_MeV"), ("quark", "mu_spinodal_quark_MeV")):
            spinodal_plot.append({
                "xi": finite(row, "xi"),
                "T_MeV": finite(row, "T_MeV"),
                "mu_MeV": finite(row, field),
                "branch": branch,
            })

    output_root.mkdir(parents=True, exist_ok=True)
    tables = output_root / "tables"
    figures = output_root / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    for name, source in (("v6_manifest.json", v6_manifest_path), ("v6_crossover_overlay.csv", crossover_path), ("v6_crossover_overlay_gap.csv", gap_path), ("v6_maxwell_surface_point_status.csv", maxwell_path), ("v6_spinodals.csv", spinodal_path), ("v6_crossover_cep_sampling_gap.csv", v5_gap_path)):
        shutil.copyfile(source, tables / name)

    crossover_surface = native_rows + same_xi_rows + boundary_rows + xi_derived_rows
    crossover_surface.sort(key=lambda row: (float(row["xi"]), float(row["mu_MeV"]), str(row["layer"]), str(row["status"])))
    crossover_fields = COMMON_FIELDS
    write_csv(tables / "crossover_certified_v6.csv", crossover_fields, native_rows)
    write_csv(tables / "crossover_derived_same_xi.csv", crossover_fields, same_xi_rows)
    write_csv(tables / "crossover_derived_xi.csv", crossover_fields, xi_derived_rows)
    write_csv(tables / "crossover_derived_boundary.csv", crossover_fields, boundary_rows)
    write_csv(tables / "crossover_surface_v7.csv", crossover_fields, crossover_surface)
    write_csv(tables / "maxwell_surface_v7.csv", maxwell_fields, maxwell_v7)
    write_csv(tables / "spinodal_surface_v7.csv", spinodal_fields, spinodal_v7)

    cep_boundary_rows = []
    for xi in sorted(cep_by_xi):
        cep = cep_by_xi[xi]
        cep_boundary_rows.append({
            "xi": xi,
            "mu_CEP_proxy_MeV": cep["mu"],
            "T_CEP_bracket_low_MeV": cep["low"],
            "T_CEP_bracket_high_MeV": cep["high"],
            "T_CEP_estimated_midpoint_MeV": 0.5 * (cep["low"] + cep["high"]),
            "boundary_mode": "estimated_midpoint",
            "strict_status": "bracket_preserved_not_strict_single_value",
            "use_in_v7": "crossover_boundary_only",
        })
    write_csv(tables / "cep_boundary_estimates.csv", cep_boundary_rows[0].keys(), cep_boundary_rows)

    for xi in sorted(native_by_xi):
        rows = sorted(native_by_xi[xi], key=lambda row: float(row["mu_MeV"]))
        cep = cep_by_xi[xi]
        same_count = sum(1 for row in same_xi_rows if xi_key(float(row["xi"])) == xi)
        boundary_count = sum(1 for row in boundary_rows if xi_key(float(row["xi"])) == xi)
        xi_mask_rows.append({
            "xi": xi,
            "xi_role": "native_slice",
            "xi_left": xi,
            "xi_right": xi,
            "xi_gap": 0.0,
            "native_mu_low_left_MeV": float(rows[0]["mu_MeV"]),
            "native_mu_low_right_MeV": float(rows[0]["mu_MeV"]),
            "common_mu_high_MeV": float(rows[-1]["mu_MeV"]),
            "cep_mu_left_MeV": cep["mu"],
            "cep_mu_right_MeV": cep["mu"],
            "materialized_rows": len(rows) + same_count + boundary_count,
            "unresolved_tail_to_nearest_CEP_MeV": 0.0,
            "coverage_status": "native_support_plus_boundary_constrained_completion",
        })
    xi_mask_rows.sort(key=lambda row: float(row["xi"]))
    write_csv(tables / "crossover_coverage_mask.csv", xi_mask_rows[0].keys(), xi_mask_rows)

    summary = [{
        "layer": "certified_layer",
        "surface": "crossover",
        "rows": len(native_rows),
        "xi_count": len(native_by_xi),
        "interpretation": "v6 native finite/converged support; not phase-reference promotion",
    }, {
        "layer": "interpolated_noncertified",
        "surface": "crossover",
        "rows": len(same_xi_rows) + len(xi_derived_rows),
        "xi_count": len(native_by_xi) + len(xi_mask_rows) - len(native_by_xi),
        "interpretation": "piecewise-linear derived values with source provenance",
    }, {
        "layer": "boundary_constrained_endpoint_interpolated_noncertified",
        "surface": "crossover",
        "rows": len(boundary_rows),
        "xi_count": len({xi_key(float(row["xi"])) for row in boundary_rows}),
        "interpretation": "terminates at the selected CEP boundary estimate; not a strict CEP solve",
    }, {
        "layer": "certified_layer",
        "surface": "maxwell",
        "rows": len(maxwell_v7),
        "xi_count": len({xi_key(finite(row, "xi")) for row in maxwell_v7}),
        "interpretation": "v6 Maxwell rows copied unchanged; unresolved geometry status retained",
    }]
    write_csv(tables / "v7_layer_summary.csv", summary[0].keys(), summary)

    make_figure(figures / "c2_phase_surfaces_mu_xi_T_v7_complete_display.png", maxwell=maxwell_v7, spinodal=spinodal_plot, certified=native_rows, same_xi=same_xi_rows + boundary_rows, xi_derived=xi_derived_rows, boundary=boundary_rows, layer_audit=False)
    make_figure(figures / "c2_phase_surfaces_mu_xi_T_v7_layer_audit.png", maxwell=maxwell_v7, spinodal=spinodal_plot, certified=native_rows, same_xi=same_xi_rows + boundary_rows, xi_derived=xi_derived_rows, boundary=boundary_rows, layer_audit=True)

    figure_files = {}
    for path in sorted(figures.glob("*.png")):
        figure_files[path.name] = {"sha256": sha256(path), "bytes": path.stat().st_size, "mode": "complete_display" if "complete_display" in path.name else "layer_audit"}
    plot_manifest = {
        "schema_version": SCHEMA,
        "calculation_sha": calculation_sha,
        "v6_manifest_hash": v6_manifest_hash,
        "solver_called": False,
        "reference_write": False,
        "triangulation": False,
        "orientation": "x=mu_q, y=xi, z=T",
        "generator": {"script": str(Path(__file__).resolve()), "script_sha256": sha256(Path(__file__).resolve())},
        "figures": figure_files,
        "visual_policy": {
            "complete_display": "same orange hue for native and derived crossover, with derived support translucent/dashed; Maxwell remains one uniform blue layer",
            "layer_audit": "native, same-xi derived and xi-derived support use distinct line styles/alpha",
            "maxwell": "copied native v6 rows only; no interpolation or geometry closure",
        },
    }
    (figures / "plot_manifest.json").write_text(json.dumps(plot_manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    claim_rows = [
        {"claim_id": "v6_native_preserved", "claim": "v6 native crossover, Maxwell and spinodal rows are preserved as separate native support layers", "status": "supported", "evidence": "tables/crossover_certified_v6.csv; tables/maxwell_surface_v7.csv; tables/spinodal_surface_v7.csv", "boundary": "native support is diagnostic and not phase-reference promotion"},
        {"claim_id": "same_xi_interpolation", "claim": "Internal crossover gaps are filled only by piecewise-linear interpolation within an existing xi slice", "status": "supported", "evidence": "tables/crossover_derived_same_xi.csv; tables/crossover_surface_v7.csv", "boundary": "error estimate is not independently calibrated"},
        {"claim_id": "xi_common_support_interpolation", "claim": "Midpoint xi rows are generated only over the common native crossover support of adjacent xi slices", "status": "supported", "evidence": "tables/crossover_derived_xi.csv; tables/crossover_coverage_mask.csv", "boundary": "tails not shared by both slices remain unresolved"},
        {"claim_id": "cep_boundary_semantics", "claim": "The selected CEP boundary estimate can terminate the derived crossover layer without reassigning any mu_q above it to crossover", "status": "supported", "evidence": "tables/cep_boundary_estimates.csv; tables/crossover_derived_boundary.csv", "boundary": "bracket remains the strict evidence; midpoint is an estimated boundary"},
        {"claim_id": "maxwell_not_filled", "claim": "v7 does not interpolate or repair the Maxwell coexistence surface", "status": "supported", "evidence": "tables/maxwell_surface_v7.csv; README.md", "boundary": "Maxwell near-CEP completion remains a separate task"},
        {"claim_id": "reference_promotion", "claim": "v7 authorizes phase-reference promotion", "status": "blocked", "evidence": "decision.json", "boundary": "derived values are non-certified and no reference write occurred"},
    ]
    write_csv(tables / "claim_ledger.csv", claim_rows[0].keys(), claim_rows)

    decision = {
        "schema_version": SCHEMA,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "calculation_sha": calculation_sha,
        "v6_manifest_hash": v6_manifest_hash,
        "solver_called": False,
        "reference_write": False,
        "verdict": "v7_crossover_derived_candidate",
        "promotion_effect": "none",
        "policy": {
            "base": "v6_postprocessed_package",
            "mu_interpolation": "piecewise_linear_same_xi",
            "xi_interpolation": "piecewise_linear_common_support_midpoint_xi",
            "cep_boundary": "estimated_midpoint boundary estimate; bracket retained",
            "maxwell": "native v6 only; no derived values",
            "outside_support": "unresolved; no ordinary extrapolation",
        },
        "parameters": {"mu_subdivisions": args.mu_subdivisions, "max_xi_gap": args.max_xi_gap},
        "counts": {
            "native_crossover_rows": len(native_rows),
            "same_xi_derived_rows": len(same_xi_rows),
            "xi_derived_rows": len(xi_derived_rows),
            "boundary_rows": len(boundary_rows),
            "maxwell_rows": len(maxwell_v7),
            "spinodal_rows": len(spinodal_v7),
            "cep_slices": len(cep_by_xi),
        },
    }
    (output_root / "decision.json").write_text(json.dumps(decision, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (output_root / "README.md").write_text(
        "# C2 phase surfaces v7：crossover 派生补全层\n\n"
        "本包以 v6 为唯一输入，生成第七版后处理结果。它不重跑 C2、不调用 equilibrium solver，也不写入 `data/reference`。v6 的 native crossover、Maxwell 和 spinodal 行均保留；本包新增的行明确标记为 `interpolated_noncertified` 或 `boundary_constrained_endpoint_interpolated_noncertified`。\n\n"
        f"固定 calculation SHA 为 `{calculation_sha}`，v6 manifest SHA-256 为 `{v6_manifest_hash}`。同一 xi 切片内的原生支持间隔用 `{args.mu_subdivisions}` 等分的分段线性插值填充；相邻 xi 只在共同 native crossover 支持区间内生成 midpoint xi。超过共同支持的尾部保留为 unresolved，不做普通外推。\n\n"
        "CEP 采用 `estimated_midpoint` 作为派生层的边界估计：它只用于把 crossover 终止在物理 CEP 边界，不等同于 strict CEP 单值求解；完整 bracket 保存在 `tables/cep_boundary_estimates.csv`。任何 `mu_q > mu_CEP` 的响应峰都不会被重新标成 crossover。\n\n"
        "Maxwell 面在本包只复制 v6 native rows，包含原有 `grid_unresolved`/geometry 状态；没有用插值伪造 `maxwell_area` 或 geometry certificate。Maxwell 近 CEP 补点是独立后续任务。\n\n"
        "两张图分别用于完整派生显示和层级审计。完整显示图用同一橙色系表示 crossover，但派生层用半透明/虚线；审计图额外区分同 xi 与 xi 方向派生。两张图均不使用三角化。\n",
        encoding="utf-8",
    )
    (output_root / "AUDIT.md").write_text(
        "## Evidence boundary\n\n"
        "- Input is the immutable v6 postprocessed package; input files are copied under `tables/v6_*`.\n"
        "- Native rows are not overwritten. Derived rows carry source keys, source layer, calculation SHA and v6 manifest hash.\n"
        "- The CEP bracket is retained as an interval. The estimated midpoint is a boundary estimate for the derived crossover display, not a strict CEP result.\n"
        "- Maxwell rows are copied only; no Maxwell interpolation, candidate selection or geometry closure is performed.\n"
        "- No ordinary extrapolation crosses a missing native support interval, a neighbouring-xi common-support boundary, or the physical `mu_q > mu_CEP` region.\n"
        "- v7 is a diagnostic/derived candidate and does not promote phase-reference or unlock transport.\n",
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
        "v6_manifest_hash": v6_manifest_hash,
        "solver_called": False,
        "reference_write": False,
        "verdict": decision["verdict"],
        "generator": {"script": str(Path(__file__).resolve()), "script_sha256": sha256(Path(__file__).resolve())},
        "input_records": {"v6_manifest": source_record(v6_manifest_path), "v6_crossover": source_record(crossover_path), "v6_gap": source_record(gap_path), "v6_maxwell": source_record(maxwell_path), "v6_spinodals": source_record(spinodal_path), "v5_cep_gap": source_record(v5_gap_path)},
        "parameters": decision["parameters"],
        "counts": decision["counts"],
        "output_files": output_files,
    }
    (output_root / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"verdict": decision["verdict"], "output_root": str(output_root), "counts": decision["counts"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
