#!/usr/bin/env python3
"""Build the review-only publication_clean_v4 display derivative.

v4 inherits publication_clean_v3 and adds three mode-B composite-curve
interpolations plus one left-branch endpoint adjustment for mode-A tau_sbar.
The endpoint adjustment is deliberately display-only: it does not fill the
first-order gap, alter raw values, regularize a propagator, or call a solver.
"""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import importlib.util
import json
import math
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[3]
V3_SCRIPT = ROOT / "scripts" / "analysis" / "relaxtime" / "build_phase_guided_publication_clean_v3.py"
V3_SPEC = importlib.util.spec_from_file_location("phase_guided_publication_clean_v3_v4", V3_SCRIPT)
if V3_SPEC is None or V3_SPEC.loader is None:  # pragma: no cover
    raise RuntimeError(f"unable to load v3 builder: {V3_SCRIPT}")
V3 = importlib.util.module_from_spec(V3_SPEC)
V3_SPEC.loader.exec_module(V3)


PARENT_DIR = (
    ROOT
    / "docs"
    / "analysis"
    / "relaxtime"
    / "phase_guided_transport"
    / "phase_guided_transport_publication_clean_v3"
)
PARENT_TABLE_DIR = PARENT_DIR / "tables"
PARENT_POINTS = PARENT_TABLE_DIR / "publication_clean_points.csv"
PARENT_MANIFEST = PARENT_DIR / "manifest.json"
PARENT_PLOT_MANIFEST = PARENT_DIR / "figures" / "plot_manifest.json"
RECIPE = (
    ROOT
    / "docs"
    / "analysis"
    / "relaxtime"
    / "phase_guided_transport"
    / "phase_guided_transport_v4_residual_smoothing"
    / "tables"
    / "publication_clean_v4_display_adjustments.csv"
)
MECHANISM_DIR = (
    ROOT
    / "docs"
    / "analysis"
    / "relaxtime"
    / "phase_guided_transport"
    / "phase_guided_transport_publication_clean_v4_mechanism_review"
)
MECHANISM_TABLE_DIR = MECHANISM_DIR / "tables"
OUT_DIR = (
    ROOT
    / "docs"
    / "analysis"
    / "relaxtime"
    / "phase_guided_transport"
    / "phase_guided_transport_publication_clean_v4"
)
TABLE_DIR = OUT_DIR / "tables"
FIGURE_DIR = OUT_DIR / "figures"
AUDIT_DIR = OUT_DIR / "audit"
OBSERVABLES = list(V3.OBSERVABLES)
CANONICAL_XI = V3.V2.canonical_xi
POINT_FIELDS = list(V3.POINT_FIELDS)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relpath(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def read_csv(path: Path) -> list[dict[str, str]]:
    return V3.read_csv(path)


def write_csv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    V3.write_csv(path, rows, fields)


def write_json(path: Path, payload: Any) -> None:
    V3.write_json(path, payload)


def finite(row: dict[str, Any], field: str) -> float:
    return V3.finite(row, field)


def fmt(value: float) -> str:
    return format(value, ".17g")


def point_key(row: dict[str, Any]) -> tuple[str, str, str, str, str]:
    return (
        str(row["mode_key"]),
        str(row["plot_panel"]),
        str(row["plot_series"]),
        CANONICAL_XI(row["xi"]),
        str(row["observable"]),
    )


def curve_key(row: dict[str, Any]) -> tuple[str, str, str, str]:
    return (
        str(row["mode_key"]),
        str(row["plot_panel"]),
        str(row["plot_series"]),
        str(row["observable"]),
    )


def local_log_residual(rows: list[dict[str, Any]], index: int) -> float:
    return V3.local_log_residual(rows, index)


def load_parent() -> tuple[list[dict[str, str]], dict[str, Any]]:
    if not PARENT_POINTS.is_file() or not PARENT_MANIFEST.is_file():
        raise FileNotFoundError("publication_clean_v3 must exist before v4")
    manifest = json.loads(PARENT_MANIFEST.read_text(encoding="utf-8"))
    if manifest.get("schema") != "phase_guided_transport_publication_clean_manifest_v3":
        raise ValueError("unexpected publication_clean_v3 parent manifest schema")
    expected = int(manifest["derived_counts"]["publication_clean_point_rows"])
    rows = read_csv(PARENT_POINTS)
    if len(rows) != expected:
        raise ValueError(f"parent point rows {len(rows)} != manifest count {expected}")
    keys: set[tuple[str, str, str, str, str]] = set()
    for row in rows:
        key = point_key(row)
        if key in keys:
            raise ValueError(f"duplicate parent point key: {key}")
        keys.add(key)
        for field in ("raw_value", "v2_clean_value", "clean_value"):
            finite(row, field)
    return rows, manifest


def load_recipe() -> list[dict[str, Any]]:
    rows = read_csv(RECIPE)
    required = {
        "adjustment_id",
        "mode_key",
        "plot_panel",
        "plot_series",
        "observable",
        "target_xi",
        "anchor_left_xi",
        "anchor_right_xi",
        "phase_policy",
        "method",
        "fit_space",
        "mechanism_status",
        "reason",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"v4 recipe is missing required columns: {sorted(required)}")
    out: list[dict[str, Any]] = []
    ids: set[tuple[str, str, str]] = set()
    for row in rows:
        key = (row["adjustment_id"], row["observable"], row["mode_key"])
        if key in ids:
            raise ValueError(f"duplicate v4 recipe entry: {key}")
        ids.add(key)
        if row["observable"] not in OBSERVABLES:
            raise ValueError(f"unknown observable in v4 recipe: {row['observable']}")
        item = dict(row)
        for field in ("target_xi", "anchor_left_xi", "anchor_right_xi"):
            item[field] = float(row[field])
        if item["method"] == "linear_interpolation":
            if not item["anchor_left_xi"] < item["target_xi"] < item["anchor_right_xi"]:
                raise ValueError(f"linear target is not interior to anchors: {row}")
        elif item["method"] == "log_linear_extrapolation":
            if not item["anchor_left_xi"] < item["anchor_right_xi"] < item["target_xi"]:
                raise ValueError(f"endpoint target is not to the right of anchors: {row}")
        else:
            raise ValueError(f"unsupported v4 method: {item['method']}")
        out.append(item)
    return out


def load_raw_context(parent_manifest: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    loaded, inventory, gaps = V3.load_raw_context(parent_manifest)
    return loaded, inventory, gaps


def gaps_for_curve(gaps: list[dict[str, Any]], mode_key: str, panel: str, series: str) -> list[dict[str, Any]]:
    return V3.gaps_for_curve(gaps, mode_key, panel, series)


def interval_crosses_gap(left: float, right: float, gaps: list[dict[str, Any]]) -> bool:
    return V3.interval_crosses_gap(left, right, gaps)


def point_inside_gap(xi: float, gaps: list[dict[str, Any]]) -> bool:
    return V3.point_inside_gap(xi, gaps)


def load_mechanism_review() -> dict[str, dict[str, Any]]:
    summary_path = MECHANISM_TABLE_DIR / "mechanism_window_summary.csv"
    upstream_path = MECHANISM_TABLE_DIR / "upstream_branch_smoothness_summary.csv"
    denominator_path = MECHANISM_TABLE_DIR / "denominator_chain_summary.csv"
    if not all(path.is_file() for path in (summary_path, upstream_path, denominator_path)):
        raise FileNotFoundError("v4 sbar mechanism review tables are required before building v4")
    summary_rows = read_csv(summary_path)
    upstream_rows = read_csv(upstream_path)
    denominator_rows = read_csv(denominator_path)
    endpoint_id = "mode_a_muB900p0_alpha1p0_xim0p003_sbar_review"
    control_id = "mode_a_muB900p0_alpha1p0_xim0p01_sbar_review"
    by_id = {str(row["window_id"]): row for row in summary_rows}
    upstream_by_id = {str(row["window_id"]): row for row in upstream_rows}
    endpoint = by_id.get(endpoint_id)
    control = by_id.get(control_id)
    endpoint_upstream = upstream_by_id.get(endpoint_id)
    control_upstream = upstream_by_id.get(control_id)
    if endpoint is None or control is None or endpoint_upstream is None or control_upstream is None:
        raise ValueError("v4 mechanism review is missing endpoint or control rows")
    if endpoint["mechanism_verdict"] != "small_denominator_supported":
        raise ValueError(f"unexpected endpoint mechanism verdict: {endpoint['mechanism_verdict']}")
    if endpoint["dominant_denominator_branch"] != "simple_1m4KPi":
        raise ValueError("endpoint is not supported by the expected simple denominator branch")
    if endpoint["denominator_sigma_alignment"].lower() != "true":
        raise ValueError("endpoint denominator and sigma peaks are not aligned")
    if float(endpoint["max_rate_reproduction_rel_error"]) > 0.05:
        raise ValueError("endpoint rate reproduction error exceeds the mechanism gate")
    if endpoint_upstream["upstream_branch_flag"].lower() != "true":
        raise ValueError("endpoint upstream sensitivity flag is missing")
    if control["mechanism_verdict"] != "small_denominator_supported":
        raise ValueError("left-branch control did not pass the denominator review")
    if control_upstream["upstream_branch_flag"].lower() != "false":
        raise ValueError("left-branch control unexpectedly has an upstream branch flag")
    target_denominators = [
        row
        for row in denominator_rows
        if row["window_id"] == endpoint_id
        and row["xi"] == "-0.003"
        and row["exchange_channel"] == "s"
    ]
    if not target_denominators:
        raise ValueError("endpoint s-channel denominator rows are missing")
    return {
        endpoint_id: {
            "mechanism_verdict": endpoint["mechanism_verdict"],
            "dominant_channels": endpoint["dominant_channels"],
            "dominant_denominator_branch": endpoint["dominant_denominator_branch"],
            "max_rate_reproduction_rel_error": endpoint["max_rate_reproduction_rel_error"],
            "denominator_sigma_alignment": endpoint["denominator_sigma_alignment"],
            "upstream_branch_flag": endpoint_upstream["upstream_branch_flag"],
            "upstream_max_rel_step": endpoint_upstream["max_rel_step"],
            "upstream_max_rel_curvature": endpoint_upstream["max_rel_curvature"],
            "control_window_id": control_id,
            "evidence_files": ";".join(
                (
                    relpath(summary_path),
                    relpath(upstream_path),
                    relpath(denominator_path),
                )
            ),
        }
    }


def endpoint_log_slope(left_xi: float, left_value: float, right_xi: float, right_value: float) -> float:
    if min(left_value, right_value) <= 0.0 or right_xi == left_xi:
        return math.nan
    return (math.log(right_value) - math.log(left_value)) / (right_xi - left_xi)


def phase_gate(
    adjustment: dict[str, Any],
    rows: list[dict[str, str]],
    left_index: int,
    target_index: int,
    right_index: int,
    curve_gaps: list[dict[str, Any]],
) -> str:
    left_xi = finite(rows[left_index], "xi")
    target_xi = finite(rows[target_index], "xi")
    right_xi = finite(rows[right_index], "xi")
    if point_inside_gap(target_xi, curve_gaps):
        raise ValueError(f"target lies inside a rendered phase gap: {adjustment['adjustment_id']}")
    if adjustment["phase_policy"] == "crossover_only":
        if interval_crosses_gap(left_xi, right_xi, curve_gaps):
            raise ValueError(f"crossover adjustment crosses a phase gap: {adjustment['adjustment_id']}")
        window_rows = rows[left_index : right_index + 1]
        if any(
            row.get("phase_reference_kind") != "crossover"
            or str(row.get("phase_structure", "")).startswith("first_order")
            for row in window_rows
        ):
            raise ValueError(f"crossover gate failed: {adjustment['adjustment_id']}")
        return "passed_crossover_only"
    if adjustment["phase_policy"] == "left_branch_endpoint_only":
        if not left_xi < right_xi < target_xi:
            raise ValueError(f"endpoint anchors are not ordered: {adjustment['adjustment_id']}")
        matching_gap = [
            gap
            for gap in curve_gaps
            if abs(float(gap["gap_xi_low"]) - target_xi) <= 1.0e-9
        ]
        if len(matching_gap) != 1:
            raise ValueError(f"endpoint is not the audited left gap endpoint: {adjustment['adjustment_id']}")
        if interval_crosses_gap(left_xi, target_xi, curve_gaps):
            raise ValueError(f"endpoint adjustment crosses a phase gap: {adjustment['adjustment_id']}")
        window_rows = rows[left_index : target_index + 1]
        if any(
            row.get("phase_reference_kind") != "first_order"
            or not str(row.get("phase_structure", "")).startswith("first_order")
            for row in window_rows
        ):
            raise ValueError(f"left-branch endpoint gate failed: {adjustment['adjustment_id']}")
        return "passed_left_branch_endpoint_only"
    raise ValueError(f"unsupported v4 phase policy: {adjustment['phase_policy']}")


def derive_display_value(
    method: str,
    fit_space: str,
    target_xi: float,
    left_xi: float,
    right_xi: float,
    left_value: float,
    right_value: float,
) -> float:
    if min(left_value, right_value) <= 0.0:
        raise ValueError("display smoothing requires positive anchor values")
    if method == "linear_interpolation":
        fraction = (target_xi - left_xi) / (right_xi - left_xi)
        if fit_space == "linear":
            return left_value + fraction * (right_value - left_value)
        if fit_space == "log":
            return math.exp(math.log(left_value) + fraction * (math.log(right_value) - math.log(left_value)))
    elif method == "log_linear_extrapolation":
        if fit_space != "log":
            raise ValueError("endpoint extrapolation must use log fit space")
        fraction = (target_xi - right_xi) / (right_xi - left_xi)
        return math.exp(math.log(right_value) + fraction * (math.log(right_value) - math.log(left_value)))
    raise ValueError(f"unsupported v4 method/fit space: {method}/{fit_space}")


def apply_adjustments(
    points: list[dict[str, str]],
    adjustments: list[dict[str, Any]],
    gaps: list[dict[str, Any]],
    mechanism_review: dict[str, dict[str, Any]] | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, Any]]]:
    grouped: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in points:
        grouped[curve_key(row)].append(row)
    for rows in grouped.values():
        rows.sort(key=lambda row: finite(row, "xi"))

    mechanism_review = mechanism_review or {}
    modified: set[tuple[str, str, str, str, str]] = set()
    audit: list[dict[str, Any]] = []
    for adjustment in adjustments:
        key = (
            adjustment["mode_key"],
            adjustment["plot_panel"],
            adjustment["plot_series"],
            adjustment["observable"],
        )
        rows = grouped.get(key)
        if rows is None:
            raise ValueError(f"recipe curve is missing: {key}")
        by_xi = {CANONICAL_XI(row["xi"]): index for index, row in enumerate(rows)}
        target_key = CANONICAL_XI(adjustment["target_xi"])
        left_key = CANONICAL_XI(adjustment["anchor_left_xi"])
        right_key = CANONICAL_XI(adjustment["anchor_right_xi"])
        if target_key not in by_xi or left_key not in by_xi or right_key not in by_xi:
            raise ValueError(f"recipe points are missing for {key}: {left_key}/{target_key}/{right_key}")
        target_index = by_xi[target_key]
        left_index = by_xi[left_key]
        right_index = by_xi[right_key]
        if adjustment["method"] == "linear_interpolation" and not left_index < target_index < right_index:
            raise ValueError(f"linear target is not an interior point: {key}")
        if adjustment["method"] == "log_linear_extrapolation" and not left_index < right_index < target_index:
            raise ValueError(f"endpoint target is not to the right of anchors: {key}")

        curve_gaps = gaps_for_curve(gaps, *key[:3])
        gate = phase_gate(adjustment, rows, left_index, target_index, right_index, curve_gaps)
        point_key_value = point_key(rows[target_index])
        if point_key_value in modified:
            raise ValueError(f"overlapping v4 adjustment: {point_key_value}")
        modified.add(point_key_value)

        target_xi = finite(rows[target_index], "xi")
        left_xi = finite(rows[left_index], "xi")
        right_xi = finite(rows[right_index], "xi")
        old_value = finite(rows[target_index], "clean_value")
        left_value = finite(rows[left_index], "clean_value")
        right_value = finite(rows[right_index], "clean_value")
        derived_value = derive_display_value(
            str(adjustment["method"]),
            str(adjustment["fit_space"]),
            target_xi,
            left_xi,
            right_xi,
            left_value,
            right_value,
        )
        if not math.isfinite(derived_value) or derived_value <= 0.0:
            raise ValueError(f"invalid v4 display value for {point_key_value}: {derived_value}")

        before_residual = ""
        after_residual = ""
        endpoint_slope_before = ""
        endpoint_slope_after = ""
        if adjustment["method"] == "linear_interpolation":
            before = local_log_residual(rows, target_index)
            rows[target_index]["clean_value"] = fmt(derived_value)
            after = local_log_residual(rows, target_index)
            before_residual = fmt(before)
            after_residual = fmt(after)
            adjustment_type = "interior_crossover_display_interpolation"
        else:
            endpoint_slope_before = fmt(endpoint_log_slope(right_xi, right_value, target_xi, old_value))
            endpoint_slope_after = fmt(endpoint_log_slope(left_xi, left_value, right_xi, right_value))
            rows[target_index]["clean_value"] = fmt(derived_value)
            adjustment_type = "left_branch_endpoint_display_extrapolation"

        rows[target_index]["display_status"] = "v4_residual_smoothed"
        rows[target_index]["value_source"] = "v3 display value with explicit v4 local display adjustment"
        rows[target_index]["canonical_data_modified"] = "False"
        rows[target_index]["smoothing_window"] = adjustment["adjustment_id"]
        rows[target_index]["smoothing_anchor_points"] = (
            f"left={CANONICAL_XI(adjustment['anchor_left_xi'])};"
            f"right={CANONICAL_XI(adjustment['anchor_right_xi'])}"
        )
        rows[target_index]["smoothing_fit_method"] = str(adjustment["method"])
        rows[target_index]["protected_by_phase_gate"] = "True"

        evidence = mechanism_review.get(str(adjustment["adjustment_id"]), {})
        if not evidence and adjustment["mechanism_status"].startswith("small_denominator"):
            evidence = mechanism_review.get("mode_a_muB900p0_alpha1p0_xim0p003_sbar_review", {})
        audit.append(
            {
                "adjustment_id": adjustment["adjustment_id"],
                "mode_key": adjustment["mode_key"],
                "plot_panel": adjustment["plot_panel"],
                "plot_series": adjustment["plot_series"],
                "observable": adjustment["observable"],
                "target_xi": CANONICAL_XI(adjustment["target_xi"]),
                "anchor_left_xi": CANONICAL_XI(adjustment["anchor_left_xi"]),
                "anchor_right_xi": CANONICAL_XI(adjustment["anchor_right_xi"]),
                "raw_value": rows[target_index]["raw_value"],
                "v2_clean_value": rows[target_index]["v2_clean_value"],
                "parent_v3_display_value": fmt(old_value),
                "derived_v4_display_value": fmt(derived_value),
                "relative_change_from_parent": fmt((derived_value - old_value) / old_value if old_value else math.nan),
                "local_log_residual_before": before_residual,
                "local_log_residual_after": after_residual,
                "endpoint_log_slope_before": endpoint_slope_before,
                "endpoint_log_slope_after": endpoint_slope_after,
                "phase_policy": adjustment["phase_policy"],
                "phase_gate": gate,
                "method": adjustment["method"],
                "fit_space": adjustment["fit_space"],
                "adjustment_type": adjustment_type,
                "mechanism_status": adjustment["mechanism_status"],
                "mechanism_evidence": evidence.get("evidence_files", "not_assessed"),
                "canonical_data_modified": "False",
                "reason": adjustment["reason"],
            }
        )
    return points, audit


def build_curve_index(points: list[dict[str, str]], audit: list[dict[str, Any]]) -> list[dict[str, Any]]:
    changed = defaultdict(int)
    for row in audit:
        changed[(row["mode_key"], row["plot_panel"], row["plot_series"], row["observable"])] += 1
    grouped: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in points:
        grouped[curve_key(row)].append(row)
    out: list[dict[str, Any]] = []
    for key, rows in sorted(grouped.items()):
        out.append(
            {
                "mode_key": key[0],
                "plot_panel": key[1],
                "plot_series": key[2],
                "observable": key[3],
                "point_count": len(rows),
                "v2_replacement_count": sum(row.get("display_status") == "replacement" for row in rows),
                "residual_smoothing_count": sum(row.get("display_status") == "residual_smoothed" for row in rows),
                "v3_adjustment_count": sum(row.get("display_status") == "v3_residual_smoothed" for row in rows),
                "v4_adjustment_count": changed[key],
                "xi_min": min(float(row["xi"]) for row in rows),
                "xi_max": max(float(row["xi"]) for row in rows),
                "canonical_data_modified": False,
            }
        )
    return out


def render_audit_figures(parent: list[dict[str, str]], final: list[dict[str, str]], audit: list[dict[str, Any]]) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    parent_by_curve: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    final_by_curve: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in parent:
        parent_by_curve[curve_key(row)].append(row)
    for row in final:
        final_by_curve[curve_key(row)].append(row)

    paths: list[Path] = []
    groups = [
        (
            "mode_b_composite_xi036_before_after",
            [
                row
                for row in audit
                if row["mode_key"] == "mode_b"
                and row["plot_panel"] == "T200.0"
                and row["plot_series"] == "muB0.0"
            ],
            (0.30, 0.42),
            "mode-B composite curves near xi=0.36",
        ),
        (
            "mode_a_tau_sbar_endpoint_before_after",
            [
                row
                for row in audit
                if row["mode_key"] == "mode_a"
                and row["plot_panel"] == "muB900.0"
                and row["plot_series"] == "alpha1.0"
            ],
            (-0.04, 0.001),
            "mode-A tau_sbar left-branch endpoint",
        ),
    ]
    for stem, group, xlim, title in groups:
        if not group:
            continue
        fig, axes = plt.subplots(1, len(group), figsize=(5.0 * len(group), 4.2), squeeze=False)
        for axis, item in zip(axes.flat, group):
            key = (
                item["mode_key"],
                item["plot_panel"],
                item["plot_series"],
                item["observable"],
            )
            before = sorted(parent_by_curve[key], key=lambda row: float(row["xi"]))
            after = sorted(final_by_curve[key], key=lambda row: float(row["xi"]))
            before = [row for row in before if xlim[0] <= float(row["xi"]) <= xlim[1]]
            after = [row for row in after if xlim[0] <= float(row["xi"]) <= xlim[1]]
            axis.plot(
                [float(row["xi"]) for row in before],
                [float(row["clean_value"]) for row in before],
                marker="o",
                ms=3,
                lw=1.2,
                label="parent v3",
                color="#4477AA",
            )
            axis.plot(
                [float(row["xi"]) for row in after],
                [float(row["clean_value"]) for row in after],
                marker="o",
                ms=3,
                lw=1.2,
                linestyle="--",
                label="v4 display",
                color="#CC3311",
            )
            axis.scatter(
                [float(item["target_xi"])],
                [float(item["derived_v4_display_value"])],
                color="#228833",
                zorder=5,
                label="adjusted point",
            )
            axis.axvline(float(item["target_xi"]), color="0.4", ls=":", lw=0.8)
            axis.set_title(item["observable"])
            axis.set_xlabel("xi")
            axis.set_ylabel("display value")
            axis.set_xlim(*xlim)
            axis.grid(alpha=0.25)
        axes.flat[0].legend(loc="best", fontsize=8)
        fig.suptitle(title)
        fig.tight_layout()
        path = AUDIT_DIR / f"{stem}.png"
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=450, bbox_inches="tight", pad_inches=0.08)
        plt.close(fig)
        paths.append(path)
    return paths


def mechanism_summary_rows(mechanism_review: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for window_id, item in sorted(mechanism_review.items()):
        rows.append(
            {
                "window_id": window_id,
                "mechanism_verdict": item["mechanism_verdict"],
                "dominant_channels": item["dominant_channels"],
                "dominant_denominator_branch": item["dominant_denominator_branch"],
                "max_rate_reproduction_rel_error": item["max_rate_reproduction_rel_error"],
                "denominator_sigma_alignment": item["denominator_sigma_alignment"],
                "upstream_branch_flag": item["upstream_branch_flag"],
                "upstream_max_rel_step": item["upstream_max_rel_step"],
                "upstream_max_rel_curvature": item["upstream_max_rel_curvature"],
                "control_window_id": item["control_window_id"],
                "evidence_files": item["evidence_files"],
                "high_rate_convergence_gate": "not_run",
            }
        )
    return rows


def claim_ledger(audit: list[dict[str, Any]]) -> list[dict[str, str]]:
    return [
        {
            "claim_id": "PC-V4-001",
            "status": "supported_with_scope_limit",
            "claim_zh": "publication_clean_v4 以 publication_clean_v3 为父级，保留 raw_value、v2_clean_value 和 v3 显示值；所有变换仍是 solver-free display-only 派生。",
            "evidence": "manifest.json; tables/publication_clean_points.csv; tables/v4_display_adjustment_map.csv",
            "scope_limit": "不修改 raw/production 数据，不重新求解，不代表正式产物晋升。",
        },
        {
            "claim_id": "PC-V4-002",
            "status": "supported_with_scope_limit",
            "claim_zh": "mode-B T=200, muB=0 的 zeta、sigma/T、sigma 在 xi=0.36 使用同一 crossover 区间的显式线性显示插值；v3 已处理的四条 tau 曲线保持继承。",
            "evidence": "tables/v4_display_adjustment_map.csv; tables/publication_clean_points.csv",
            "scope_limit": "这是图形残余肩部的显示处理，不是对底层 channel 或传播子做正则化。",
        },
        {
            "claim_id": "PC-V4-003",
            "status": "supported_with_scope_limit",
            "claim_zh": "mode-A muB=900, alpha_T=1 的 tau_sbar 左支端点 xi=-0.003 使用只依赖左支锚点的 log-linear display extrapolation；端点 gap 保留为硬断线，不填补另一相分支。",
            "evidence": "tables/v4_display_adjustment_map.csv; tables/boundary_gap_map.csv; audit/mode_a_tau_sbar_endpoint_before_after.png",
            "scope_limit": "端点显示值不再是原始 branch solution 的定量值，不能用于精确导数或相变跃迁幅度。",
        },
        {
            "claim_id": "PC-V4-004",
            "status": "supported_with_scope_limit",
            "claim_zh": "同一点机制审计支持 tau_sbar 主导 usbar_to_usbar channel 的 simple 1-4KPi 分母敏感性，并与 sigma 峰同 band 对齐；端点同时存在上游 m_u 曲率，因此不能把整段结构归因于单一分母机制。",
            "evidence": "tables/mechanism_review_summary.csv; " + relpath(MECHANISM_TABLE_DIR / "denominator_chain_summary.csv"),
            "scope_limit": "没有额外 high-rate convergence gate；该结论是机制补证，不是物理正则化或 production-grade 收敛证明。",
        },
        {
            "claim_id": "PC-V4-005",
            "status": "author_check",
            "claim_zh": "v4 图层仅供作者审查，manuscript_eligible=false；采用前应确认是否接受对一阶左端点显示值的单点平滑。",
            "evidence": "README.md; manifest.json; audit/",
            "scope_limit": "作者确认前继续以 v3/raw 作为定量数据源。",
        },
    ]


def render_readme(
    parent_manifest: dict[str, Any],
    recipe: list[dict[str, Any]],
    audit: list[dict[str, Any]],
    mechanism_review: dict[str, dict[str, Any]],
    figure_count: int,
) -> str:
    rows = "\n".join(
        f"| {item['observable']} | {item['target_xi']} | {item['parent_v3_display_value']} | {item['derived_v4_display_value']} | {float(item['relative_change_from_parent']):.6g} | {item['adjustment_type']} |"
        for item in audit
    )
    endpoint = mechanism_review["mode_a_muB900p0_alpha1p0_xim0p003_sbar_review"]
    return f"""# Issue #130 RS publication_clean_v4 review candidate

## Purpose and boundary

This package is an independent display-only derivative of
publication_clean_v3. It keeps the v3 four-tau adjustment and adds three
mode-B composite-curve adjustments at T=200 MeV, muB=0 MeV, xi=0.36.
It also adds one one-sided left-branch endpoint adjustment to mode-A
tau_sbar at muB=900 MeV, alpha_T=1.0, xi=-0.003.

No raw CSV, production registry, solver output, or v3 package is modified.
The endpoint adjustment does not fill the rendered first-order gap
[-0.003, +0.003], and it is not a propagator finite-width prescription,
denominator clipping rule, or numerical re-solve.

## Display transformations

| scope | rule |
| --- | --- |
| mode-B composite curves | Linear interpolation between xi=0.35 and xi=0.37 using the v3 display anchors. |
| mode-A tau_sbar endpoint | Log-linear extrapolation from the left-branch anchors xi=-0.02 and xi=-0.01 to xi=-0.003. |
| phase gate | The first rule is crossover-only. The endpoint rule accepts only the audited left endpoint, uses anchors strictly on the left branch, and keeps the gap hard. |

## Point audit

| observable | target xi | parent v3 value | v4 display value | relative change | adjustment type |
| --- | ---: | ---: | ---: | ---: | --- |
{rows}

The endpoint raw value remains in tables/publication_clean_points.csv as
raw_value, while clean_value is the v4 display value. Do not use the v4
endpoint value for exact branch derivatives, jump amplitudes, or quantitative
phase-transition fitting.

## Mechanism evidence for tau_sbar

The dedicated same-branch diagnostic is recorded in
tables/mechanism_review_summary.csv. At xi=-0.003 the dominant
usbar_to_usbar channel covers about 80.8 percent of the sbar rate; its simple
1-4KPi denominator peak and the full-process sigma peak occupy the same
near-threshold band, and the production rate is reproduced with relative
error {endpoint['max_rate_reproduction_rel_error']}. The control point
xi=-0.01 has the same denominator verdict but no upstream branch flag.

The endpoint itself also has an upstream branch-sensitivity flag
(max relative step {endpoint['upstream_max_rel_step']}, max local curvature
{endpoint['upstream_max_rel_curvature']}). The defensible interpretation is
therefore combined denominator sensitivity plus left-branch endpoint response,
not a claim that the visible feature is purely numerical noise. The local
high-rate convergence gate was not run in this review package.

## Provenance

- parent artifact: {relpath(PARENT_DIR)}
- parent manifest SHA256: {sha256_file(PARENT_MANIFEST)}
- parent plot manifest SHA256: {sha256_file(PARENT_PLOT_MANIFEST)}
- adjustment recipe: {relpath(RECIPE)}
- adjustment recipe SHA256: {sha256_file(RECIPE)}
- mechanism review directory: {relpath(MECHANISM_DIR)}
- solver called for this derivative: false
- canonical/raw data modified: false
- production write: false
- manuscript eligible: false
- publication figures: {figure_count}

## Reproduction

    python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v4.py
    python scripts/analysis/relaxtime/sync_phase_guided_publication_clean_figure_layer_v4.py
    python -m pytest tests/unit/python/test_phase_guided_publication_clean_v4.py
"""


def main() -> None:
    if (OUT_DIR / "manifest.json").exists():
        raise FileExistsError(f"refusing to overwrite completed v4 analysis package: {OUT_DIR}")
    parent_points, parent_manifest = load_parent()
    recipe = load_recipe()
    _, inventory, gaps = load_raw_context(parent_manifest)
    mechanism_review = load_mechanism_review()
    final_points = [dict(row) for row in parent_points]
    final_points, audit = apply_adjustments(final_points, recipe, gaps, mechanism_review)
    if len(audit) != len(recipe):
        raise ValueError(f"applied audit rows {len(audit)} != recipe rows {len(recipe)}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    for source_table in PARENT_TABLE_DIR.glob("*.csv"):
        shutil.copy2(source_table, TABLE_DIR / source_table.name)

    V3.V2.FIGURE_DIR = FIGURE_DIR
    figure_paths, figure_specs = V3.V2.render_figures(final_points, gaps, OBSERVABLES)
    audit_paths = render_audit_figures(parent_points, final_points, audit)

    write_csv(TABLE_DIR / "publication_clean_points.csv", final_points, POINT_FIELDS)
    curve_fields = [
        "mode_key",
        "plot_panel",
        "plot_series",
        "observable",
        "point_count",
        "v2_replacement_count",
        "residual_smoothing_count",
        "v3_adjustment_count",
        "v4_adjustment_count",
        "xi_min",
        "xi_max",
        "canonical_data_modified",
    ]
    write_csv(TABLE_DIR / "curve_index.csv", build_curve_index(final_points, audit), curve_fields)
    adjustment_fields = [
        "adjustment_id",
        "mode_key",
        "plot_panel",
        "plot_series",
        "observable",
        "target_xi",
        "anchor_left_xi",
        "anchor_right_xi",
        "raw_value",
        "v2_clean_value",
        "parent_v3_display_value",
        "derived_v4_display_value",
        "relative_change_from_parent",
        "local_log_residual_before",
        "local_log_residual_after",
        "endpoint_log_slope_before",
        "endpoint_log_slope_after",
        "phase_policy",
        "phase_gate",
        "method",
        "fit_space",
        "adjustment_type",
        "mechanism_status",
        "mechanism_evidence",
        "canonical_data_modified",
        "reason",
    ]
    write_csv(TABLE_DIR / "v4_display_adjustment_map.csv", audit, adjustment_fields)
    write_csv(
        TABLE_DIR / "mechanism_review_summary.csv",
        mechanism_summary_rows(mechanism_review),
        [
            "window_id",
            "mechanism_verdict",
            "dominant_channels",
            "dominant_denominator_branch",
            "max_rate_reproduction_rel_error",
            "denominator_sigma_alignment",
            "upstream_branch_flag",
            "upstream_max_rel_step",
            "upstream_max_rel_curvature",
            "control_window_id",
            "evidence_files",
            "high_rate_convergence_gate",
        ],
    )
    write_csv(
        TABLE_DIR / "claim_ledger.csv",
        claim_ledger(audit),
        ["claim_id", "status", "claim_zh", "evidence", "scope_limit"],
    )

    generated_at = dt.datetime.now(dt.timezone.utc).isoformat()
    axis_counts: dict[str, int] = defaultdict(int)
    figure_assets: list[dict[str, Any]] = []
    for spec in figure_specs:
        axis_counts[str(spec["axis_scale"])] += 1
        path = spec["path"]
        figure_assets.append(
            {
                "path": relpath(path),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
                "mode_key": spec["mode_key"],
                "plot_panel": spec["plot_panel"],
                "observable": spec["observable"],
                "axis_scale": spec["axis_scale"],
                "axis_scale_reason": spec["axis_scale_reason"],
                "data_min": spec["data_min"],
                "data_max": spec["data_max"],
                "dynamic_range_ratio": spec["dynamic_range_ratio"],
                "first_order_gap_present": spec["first_order_gap_present"],
            }
        )

    plot_manifest = {
        "schema": "phase_guided_transport_publication_clean_plot_manifest_v4",
        "marker_contract": "phase_endpoint_gap_v2",
        "case": parent_manifest["case"],
        "generated_at": generated_at,
        "base_git_commit": V3.V2.git_head(),
        "generator": relpath(Path(__file__).resolve()),
        "generator_sha256": sha256_file(Path(__file__).resolve()),
        "source_parent_artifact": relpath(PARENT_DIR),
        "source_parent_manifest_sha256": sha256_file(PARENT_MANIFEST),
        "source_parent_plot_manifest_sha256": sha256_file(PARENT_PLOT_MANIFEST),
        "adjustment_recipe": relpath(RECIPE),
        "adjustment_recipe_sha256": sha256_file(RECIPE),
        "mechanism_review_artifact": relpath(MECHANISM_DIR),
        "observables": OBSERVABLES,
        "boundary_gap_count": len(gaps),
        "axis_scale_policy": "inherited from parent v2 residual-smoothed renderer",
        "axis_scale_counts": dict(sorted(axis_counts.items())),
        "manuscript_eligible": False,
        "canonical_data_modified": False,
        "production_write": False,
        "solver_called": False,
        "rendering_semantics": "v3 display values plus three mode-B composite interpolations and one one-sided mode-A tau_sbar endpoint extrapolation; raw and v2 values retained; first-order gaps remain split",
        "changed_point_count": len(audit),
        "figures": figure_assets,
    }
    write_json(FIGURE_DIR / "plot_manifest.json", plot_manifest)

    readme_path = OUT_DIR / "README.md"
    readme_path.write_text(
        render_readme(parent_manifest, recipe, audit, mechanism_review, len(figure_paths)),
        encoding="utf-8",
    )

    output_paths = [
        readme_path,
        *sorted(TABLE_DIR.glob("*.csv")),
        FIGURE_DIR / "plot_manifest.json",
        *figure_paths,
        *audit_paths,
    ]
    calculation_sha = parent_manifest.get("calculation_sha")
    if calculation_sha is None:
        calculation_sha = parent_manifest["source_inputs"][0].get("calculation_sha")
    workflow_head_sha = parent_manifest.get("workflow_head_sha")
    if workflow_head_sha is None:
        workflow_head_sha = parent_manifest["source_inputs"][0].get("workflow_head_sha")
    package_manifest = {
        "schema": "phase_guided_transport_publication_clean_manifest_v4",
        "marker_contract": "phase_endpoint_gap_v2",
        "case": parent_manifest["case"],
        "generated_at": generated_at,
        "base_git_commit": V3.V2.git_head(),
        "generator": relpath(Path(__file__).resolve()),
        "generator_sha256": sha256_file(Path(__file__).resolve()),
        "source_parent_artifact": relpath(PARENT_DIR),
        "source_parent_manifest_sha256": sha256_file(PARENT_MANIFEST),
        "source_parent_plot_manifest_sha256": sha256_file(PARENT_PLOT_MANIFEST),
        "adjustment_recipe": relpath(RECIPE),
        "adjustment_recipe_sha256": sha256_file(RECIPE),
        "mechanism_review_artifact": relpath(MECHANISM_DIR),
        "mechanism_review_summary_sha256": sha256_file(TABLE_DIR / "mechanism_review_summary.csv"),
        "status": "derived_author_review_required",
        "manuscript_eligible": False,
        "canonical_data_modified": False,
        "production_write": False,
        "solver_called": False,
        "source_solver_called": parent_manifest.get("source_solver_called", True),
        "calculation_sha": calculation_sha,
        "workflow_head_sha": workflow_head_sha,
        "source_inputs": inventory,
        "derived_counts": {
            "parent_point_rows": len(parent_points),
            "publication_clean_point_rows": len(final_points),
            "curve_rows": len(build_curve_index(final_points, audit)),
            "adjustment_recipe_rows": len(recipe),
            "adjusted_point_rows": len(audit),
            "publication_figure_count": len(figure_paths),
            "audit_figure_count": len(audit_paths),
        },
        "adjustment_summary": {
            "targets": [
                "mode_b/T200.0/muB0.0/xi=0.36/zeta",
                "mode_b/T200.0/muB0.0/xi=0.36/sigma_over_T",
                "mode_b/T200.0/muB0.0/xi=0.36/sigma",
                "mode_a/muB900.0/alpha1.0/xi=-0.003/tau_sbar",
            ],
            "methods": {
                "mode_b_composite": "linear_interpolation_between_v3_display_anchors",
                "mode_a_tau_sbar_endpoint": "left_branch_log_linear_extrapolation",
            },
            "phase_gates": {
                "mode_b_composite": "crossover_only_passed",
                "mode_a_tau_sbar_endpoint": "left_branch_endpoint_only_passed",
            },
        },
        "known_boundaries": [
            "publication_clean_v3 remains unchanged",
            "publication_clean_v2 and publication_clean_v2_residual_smoothed remain unchanged",
            "v4 is a review-only display derivative, not a solver or production rerun",
            "raw_value and v2_clean_value remain available for quantitative use",
            "v3 display values and v3 adjustment provenance remain inherited",
            "the mode-A first-order gap and both endpoint rows remain rendered as a hard split",
            "the tau_sbar endpoint display value is not a replacement equilibrium branch solution",
            "denominator-chain evidence is mechanism context, not a finite-width regularization",
            "the local high-rate convergence gate was not run",
        ],
        "outputs": [
            {"path": relpath(path), "sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in output_paths
        ],
    }
    write_json(OUT_DIR / "manifest.json", package_manifest)
    print(
        json.dumps(
            {
                "output": relpath(OUT_DIR),
                "manifest": relpath(OUT_DIR / "manifest.json"),
                "publication_figures": len(figure_paths),
                "audit_figures": len(audit_paths),
                "adjusted_points": len(audit),
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
