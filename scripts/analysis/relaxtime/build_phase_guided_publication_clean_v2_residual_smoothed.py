#!/usr/bin/env python3
"""Build a second, display-only residual smoothing derivative of v2.

The source is the existing publication_clean_v2 point table.  No raw result,
solver output, production registry, or existing v2 artifact is overwritten.
Only explicit windows from the residual-smoothing recipe are changed.  Each
window uses a monotone Hermite trend in log(y) for positive observables and
blends that trend with the current v2 display value by an audited strength.
Rendered first-order gaps are hard barriers; branch-local windows may operate
inside one labelled continuation branch but may not cross a rendered gap.
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
V2_SCRIPT = ROOT / "scripts" / "analysis" / "relaxtime" / "build_phase_guided_publication_clean_v2.py"
V2_SPEC = importlib.util.spec_from_file_location("phase_guided_publication_clean_v2", V2_SCRIPT)
if V2_SPEC is None or V2_SPEC.loader is None:  # pragma: no cover
    raise RuntimeError(f"unable to load v2 builder: {V2_SCRIPT}")
V2 = importlib.util.module_from_spec(V2_SPEC)
V2_SPEC.loader.exec_module(V2)


SOURCE_DIR = V2.OUT_DIR
SOURCE_TABLE_DIR = SOURCE_DIR / "tables"
SOURCE_POINTS = SOURCE_TABLE_DIR / "publication_clean_points.csv"
SOURCE_MANIFEST = SOURCE_DIR / "manifest.json"
RECIPE = (
    ROOT
    / "docs"
    / "analysis"
    / "relaxtime"
    / "phase_guided_transport"
    / "phase_guided_transport_v2_residual_smoothing"
    / "tables"
    / "residual_smoothing_windows.csv"
)
OUT_DIR = (
    ROOT
    / "docs"
    / "analysis"
    / "relaxtime"
    / "phase_guided_transport"
    / "phase_guided_transport_publication_clean_v2_residual_smoothed"
)
TABLE_DIR = OUT_DIR / "tables"
FIGURE_DIR = OUT_DIR / "figures"
OBSERVABLES = list(V2.DISPLAY_FIELDS)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relpath(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        lines = [line for line in handle if line.strip() and not line.startswith("#")]
    return list(csv.DictReader(lines))


def write_csv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def finite(row: dict[str, Any], field: str) -> float:
    value = float(row[field])
    if not math.isfinite(value):
        raise ValueError(f"non-finite {field}: {row[field]}")
    return value


def local_log_residual(rows: list[dict[str, Any]], index: int, value_field: str) -> float:
    """Return the three-point log residual, or NaN at an invalid boundary."""
    if index <= 0 or index >= len(rows) - 1:
        return math.nan
    values = [finite(rows[item], value_field) for item in (index - 1, index, index + 1)]
    if min(values) <= 0.0:
        return math.nan
    return math.log(values[1]) - 0.5 * (math.log(values[0]) + math.log(values[2]))


def point_key(row: dict[str, Any]) -> tuple[str, str, str, str, str]:
    return (
        str(row["mode_key"]),
        str(row["plot_panel"]),
        str(row["plot_series"]),
        V2.canonical_xi(row["xi"]),
        str(row["observable"]),
    )


def curve_key(row: dict[str, Any]) -> tuple[str, str, str, str]:
    return (
        str(row["mode_key"]),
        str(row["plot_panel"]),
        str(row["plot_series"]),
        str(row["observable"]),
    )


def load_source_points() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if not SOURCE_POINTS.is_file() or not SOURCE_MANIFEST.is_file():
        raise FileNotFoundError("publication_clean_v2 must exist before residual smoothing")
    source_manifest = json.loads(SOURCE_MANIFEST.read_text(encoding="utf-8"))
    expected = int(source_manifest["derived_counts"]["publication_clean_point_rows"])
    rows = read_csv(SOURCE_POINTS)
    if len(rows) != expected:
        raise ValueError(f"source point rows {len(rows)} != v2 manifest count {expected}")
    seen: set[tuple[str, str, str, str, str]] = set()
    for row in rows:
        key = point_key(row)
        if key in seen:
            raise ValueError(f"duplicate source point key: {key}")
        seen.add(key)
        raw = finite(row, "raw_value")
        clean = finite(row, "clean_value")
        row["v2_clean_value"] = clean
        row["clean_value"] = clean
        row["raw_value"] = raw
    return rows, source_manifest


def load_recipe() -> list[dict[str, Any]]:
    rows = read_csv(RECIPE)
    if not rows:
        raise ValueError("residual smoothing recipe is empty")
    out: list[dict[str, Any]] = []
    ids: set[str] = set()
    for row in rows:
        window_id = row["window_id"]
        if window_id in ids:
            raise ValueError(f"duplicate window id: {window_id}")
        ids.add(window_id)
        strength = float(row["strength"])
        if not 0.0 <= strength <= 1.0:
            raise ValueError(f"invalid smoothing strength for {window_id}: {strength}")
        observables = tuple(item for item in row["observables"].split(";") if item)
        unknown = set(observables) - set(OBSERVABLES)
        if unknown:
            raise ValueError(f"unknown observables in {window_id}: {sorted(unknown)}")
        item = dict(row)
        item["anchor_left_xi"] = float(row["anchor_left_xi"])
        item["anchor_right_xi"] = float(row["anchor_right_xi"])
        item["strength"] = strength
        item["observables_tuple"] = observables
        out.append(item)
    return out


def build_raw_phase_index(loaded: dict[str, dict[str, Any]]) -> dict[tuple[str, str, str, str], dict[str, str]]:
    out: dict[tuple[str, str, str, str], dict[str, str]] = {}
    for mode_key, payload in loaded.items():
        for row in payload["rows"]:
            key = (mode_key, row["plot_panel"], row["plot_series"], V2.canonical_xi(row["xi"]))
            if key in out:
                raise ValueError(f"duplicate raw phase key: {key}")
            out[key] = row
    return out


def gap_records(loaded: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    gaps, _ = V2.build_boundary_gap_map(loaded)
    return gaps


def gaps_for_curve(
    gaps: list[dict[str, Any]],
    mode_key: str,
    plot_panel: str,
    plot_series: str,
) -> list[dict[str, Any]]:
    """Scope phase barriers to the curve they were derived from.

    Fixture/tests may provide interval-only records; those remain applicable
    to preserve the small helper API used by the unit tests.
    """
    identity = {
        "mode_key": mode_key,
        "plot_panel": plot_panel,
        "plot_series": plot_series,
    }
    return [
        gap
        for gap in gaps
        if all(
            field not in gap or str(gap[field]) == expected
            for field, expected in identity.items()
        )
    ]


def interval_has_open_gap(x_left: float, x_right: float, gaps: list[dict[str, Any]]) -> bool:
    lo, hi = min(x_left, x_right), max(x_left, x_right)
    return any(lo < gap["gap_xi_high"] and hi > gap["gap_xi_low"] and not (
        hi <= gap["gap_xi_low"] or lo >= gap["gap_xi_high"]
    ) for gap in gaps)


def point_is_inside_gap(xi: float, gaps: list[dict[str, Any]]) -> bool:
    return any(gap["gap_xi_low"] < xi < gap["gap_xi_high"] for gap in gaps)


def _limited_slope(delta: float, slope: float) -> float:
    if delta == 0.0 or not math.isfinite(slope) or delta * slope <= 0.0:
        return 0.0
    limit = 3.0 * abs(delta)
    return max(-limit, min(limit, slope))


def _hermite(x: float, x0: float, y0: float, x1: float, y1: float, m0: float, m1: float) -> float:
    width = x1 - x0
    t = (x - x0) / width
    h00 = 2.0 * t**3 - 3.0 * t**2 + 1.0
    h10 = t**3 - 2.0 * t**2 + t
    h01 = -2.0 * t**3 + 3.0 * t**2
    h11 = t**3 - t**2
    return h00 * y0 + h10 * width * m0 + h01 * y1 + h11 * width * m1


def _slope_from_neighbor(
    rows: list[dict[str, Any]],
    index: int,
    direction: int,
    value_field: str,
    gaps: list[dict[str, Any]],
) -> float | None:
    candidates = [index - 1, index + 1] if direction < 0 else [index + 1, index - 1]
    x_here = finite(rows[index], "xi")
    for neighbor in candidates:
        if neighbor < 0 or neighbor >= len(rows):
            continue
        x_other = finite(rows[neighbor], "xi")
        if interval_has_open_gap(x_here, x_other, gaps):
            continue
        y_here = finite(rows[index], value_field)
        y_other = finite(rows[neighbor], value_field)
        if y_here <= 0.0 or y_other <= 0.0:
            return (y_here - y_other) / (x_here - x_other)
        return (math.log(y_here) - math.log(y_other)) / (x_here - x_other)
    return None


def _trend_values(
    rows: list[dict[str, Any]],
    left_index: int,
    right_index: int,
    value_field: str,
    gaps: list[dict[str, Any]],
    value_space: str,
) -> dict[str, float]:
    x0 = finite(rows[left_index], "xi")
    x1 = finite(rows[right_index], "xi")
    values = [finite(rows[index], value_field) for index in range(left_index, right_index + 1)]
    use_log = value_space == "log" and all(value > 0.0 for value in values)
    if use_log:
        z0 = math.log(values[0])
        z1 = math.log(values[-1])
        delta = (z1 - z0) / (x1 - x0)
        raw_m0 = _slope_from_neighbor(rows, left_index, -1, value_field, gaps)
        raw_m1 = _slope_from_neighbor(rows, right_index, 1, value_field, gaps)
        m0 = _limited_slope(delta, raw_m0 if raw_m0 is not None else delta)
        m1 = _limited_slope(delta, raw_m1 if raw_m1 is not None else delta)
        return {
            V2.canonical_xi(rows[index]["xi"]): math.exp(
                _hermite(finite(rows[index], "xi"), x0, z0, x1, z1, m0, m1)
            )
            for index in range(left_index, right_index + 1)
        }

    y0, y1 = values[0], values[-1]
    delta = (y1 - y0) / (x1 - x0)
    raw_m0 = _slope_from_neighbor(rows, left_index, -1, value_field, gaps)
    raw_m1 = _slope_from_neighbor(rows, right_index, 1, value_field, gaps)
    m0 = _limited_slope(delta, raw_m0 if raw_m0 is not None else delta)
    m1 = _limited_slope(delta, raw_m1 if raw_m1 is not None else delta)
    return {
        V2.canonical_xi(rows[index]["xi"]): _hermite(
            finite(rows[index], "xi"), x0, y0, x1, y1, m0, m1
        )
        for index in range(left_index, right_index + 1)
    }


def _local_roughness(
    rows: list[dict[str, Any]],
    values: dict[str, float],
    value_space: str,
) -> float:
    """Measure the largest adjacent slope change over one local window."""
    ordered = sorted(rows, key=lambda row: finite(row, "xi"))
    sampled = [values[V2.canonical_xi(row["xi"])] for row in ordered]
    use_log = value_space == "log" and all(value > 0.0 for value in sampled)
    transformed = [math.log(value) for value in sampled] if use_log else sampled
    slopes = [
        (right - left)
        / (finite(ordered[index + 1], "xi") - finite(ordered[index], "xi"))
        for index, (left, right) in enumerate(zip(transformed, transformed[1:]))
    ]
    return max((abs(right - left) for left, right in zip(slopes, slopes[1:])), default=0.0)


def _phase_gate(
    window: dict[str, Any],
    rows: list[dict[str, Any]],
    raw_phase: dict[tuple[str, str, str, str], dict[str, str]],
    gaps: list[dict[str, Any]],
) -> tuple[bool, str, bool]:
    mode_key = window["mode_key"]
    panel = window["plot_panel"]
    series = window["plot_series"]
    left = window["anchor_left_xi"]
    right = window["anchor_right_xi"]
    relevant = [row for row in rows if left <= finite(row, "xi") <= right]
    if not relevant:
        return False, "no rows in anchor interval", False
    if any(point_is_inside_gap(finite(row, "xi"), gaps) for row in relevant):
        return False, "anchor interval contains a rendered first-order gap", True
    if interval_has_open_gap(left, right, gaps):
        return False, "anchor interval crosses a rendered first-order gap", True
    try:
        phase_rows = [
            raw_phase[(mode_key, panel, series, V2.canonical_xi(row["xi"]))]
            for row in relevant
        ]
    except KeyError as error:
        return False, f"missing phase metadata: {error.args[0]}", False
    policy = window["phase_policy"]
    if policy == "crossover_only":
        passed = all(
            row["phase_reference_kind"] == "crossover" and row["phase_structure"] != "first_order"
            for row in relevant
        )
        return passed, (
            "crossover/non-first-order gate passed"
            if passed
            else "crossover-only window contains a first-order-labelled row"
        ), False
    if policy == "branch_local":
        phases = {row.get("phase_curr", "") for row in phase_rows}
        passed = len(phases) == 1 and "" not in phases
        return passed, (
            f"single continuation branch: {next(iter(phases))}"
            if passed
            else "branch-local window changes phase_curr or lacks phase_curr"
        ), False
    raise ValueError(f"unknown phase policy: {policy}")


def apply_residual_smoothing(
    points: list[dict[str, Any]],
    windows: list[dict[str, Any]],
    raw_phase: dict[tuple[str, str, str, str], dict[str, str]],
    gaps: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    grouped: dict[tuple[str, str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for point in points:
        grouped[curve_key(point)].append(point)
    for rows in grouped.values():
        rows.sort(key=lambda row: finite(row, "xi"))

    modified_keys: set[tuple[str, str, str, str, str]] = set()
    point_audit: list[dict[str, Any]] = []
    window_audit: list[dict[str, Any]] = []
    curve_audit: list[dict[str, Any]] = []
    recipe_hash = sha256_file(RECIPE)

    for window in windows:
        audit_start = len(point_audit)
        target_count = 0
        changed_count = 0
        applied_curve_count = 0
        retained_curve_count = 0
        retained_curve_ids: list[str] = []
        max_change = 0.0
        gate_errors: list[str] = []
        for observable in window["observables_tuple"]:
            key = (
                window["mode_key"],
                window["plot_panel"],
                window["plot_series"],
                observable,
            )
            rows = grouped.get(key)
            if rows is None:
                gate_errors.append(f"missing curve {key}")
                continue
            xi_to_index = {V2.canonical_xi(row["xi"]): index for index, row in enumerate(rows)}
            left_key = V2.canonical_xi(window["anchor_left_xi"])
            right_key = V2.canonical_xi(window["anchor_right_xi"])
            if left_key not in xi_to_index or right_key not in xi_to_index:
                gate_errors.append(f"missing anchor {observable}:{left_key}/{right_key}")
                continue
            left_index = xi_to_index[left_key]
            right_index = xi_to_index[right_key]
            if not left_index < right_index:
                gate_errors.append(f"unordered anchors for {observable}")
                continue
            curve_gaps = gaps_for_curve(
                gaps,
                window["mode_key"],
                window["plot_panel"],
                window["plot_series"],
            )
            passed, gate_reason, gap_guard = _phase_gate(window, rows, raw_phase, curve_gaps)
            if not passed:
                gate_errors.append(f"{observable}: {gate_reason}")
                continue
            trends = _trend_values(
                rows,
                left_index,
                right_index,
                "clean_value",
                curve_gaps,
                window["value_space"],
            )
            window_rows = rows[left_index : right_index + 1]
            source_values = {
                V2.canonical_xi(row["xi"]): finite(row, "clean_value")
                for row in window_rows
            }
            original_rows = {point_key(row): dict(row) for row in window_rows}
            pending_audits: list[tuple[int, dict[str, Any]]] = []
            curve_modified_keys: list[tuple[str, str, str, str, str]] = []
            curve_changed_count = 0
            curve_max_change = 0.0
            for index in range(left_index + 1, right_index):
                row = rows[index]
                xi = V2.canonical_xi(row["xi"])
                key_point = point_key(row)
                if key_point in modified_keys:
                    gate_errors.append(f"overlapping window at {key_point}")
                    continue
                modified_keys.add(key_point)
                curve_modified_keys.append(key_point)
                old_value = finite(row, "clean_value")
                residual_before = local_log_residual(rows, index, "clean_value")
                trend_value = trends[xi]
                strength = window["strength"]
                if window["value_space"] == "log" and old_value > 0.0 and trend_value > 0.0:
                    final_value = math.exp(
                        (1.0 - strength) * math.log(old_value) + strength * math.log(trend_value)
                    )
                else:
                    final_value = old_value + strength * (trend_value - old_value)
                if not math.isfinite(final_value) or final_value <= 0.0:
                    raise ValueError(f"invalid residual-smoothed value at {key_point}: {final_value}")
                row["clean_value"] = final_value
                row["display_status"] = "residual_smoothed"
                row["value_source"] = "v2 clean value blended with monotone local trend"
                row["canonical_data_modified"] = False
                row["smoothing_window"] = window["window_id"]
                row["smoothing_anchor_points"] = (
                    f"left={V2.canonical_xi(window['anchor_left_xi'])};"
                    f"right={V2.canonical_xi(window['anchor_right_xi'])}"
                )
                row["smoothing_fit_method"] = (
                    "monotone_hermite_log_blend"
                    if window["value_space"] == "log"
                    else "monotone_hermite_linear_blend"
                )
                row["protected_by_phase_gate"] = True
                relative_change = (final_value - old_value) / old_value if old_value else math.nan
                pending_audits.append(
                    (
                        index,
                        {
                            "window_id": window["window_id"],
                            "mode_key": row["mode_key"],
                            "plot_panel": row["plot_panel"],
                            "plot_series": row["plot_series"],
                            "observable": observable,
                            "xi": xi,
                            "v2_clean_value": old_value,
                            "trend_value": trend_value,
                            "final_clean_value": final_value,
                            "relative_change_from_v2": relative_change,
                            "strength": strength,
                            "value_space": window["value_space"],
                            "phase_policy": window["phase_policy"],
                            "anchor_points": row["smoothing_anchor_points"],
                            "fit_method": row["smoothing_fit_method"],
                            "local_residual_before": residual_before,
                            "protected_by_phase_gate": True,
                            "canonical_data_modified": False,
                        },
                    )
                )
                if abs(relative_change) > 0.0:
                    curve_changed_count += 1
                    curve_max_change = max(curve_max_change, abs(relative_change))
            candidate_values = {
                V2.canonical_xi(row["xi"]): finite(row, "clean_value")
                for row in window_rows
            }
            source_roughness = _local_roughness(window_rows, source_values, window["value_space"])
            candidate_roughness = _local_roughness(window_rows, candidate_values, window["value_space"])
            roughness_tolerance = 1.0e-12 * max(1.0, source_roughness)
            roughness_passed = candidate_roughness <= source_roughness + roughness_tolerance
            roughness_ratio = (
                candidate_roughness / source_roughness
                if source_roughness > 0.0
                else (0.0 if candidate_roughness == 0.0 else math.inf)
            )
            curve_id = f"{window['window_id']}:{observable}"
            if not roughness_passed:
                for index, _ in pending_audits:
                    row = rows[index]
                    original = original_rows[point_key(row)]
                    row.clear()
                    row.update(original)
                modified_keys.difference_update(curve_modified_keys)
                retained_curve_count += 1
                retained_curve_ids.append(curve_id)
                curve_audit.append(
                    {
                        "window_id": window["window_id"],
                        "mode_key": window["mode_key"],
                        "plot_panel": window["plot_panel"],
                        "plot_series": window["plot_series"],
                        "observable": observable,
                        "anchor_left_xi": V2.canonical_xi(window["anchor_left_xi"]),
                        "anchor_right_xi": V2.canonical_xi(window["anchor_right_xi"]),
                        "interior_point_count": right_index - left_index - 1,
                        "source_roughness": source_roughness,
                        "candidate_roughness": candidate_roughness,
                        "roughness_ratio": roughness_ratio,
                        "changed_point_count": 0,
                        "max_abs_relative_change": 0.0,
                        "action": "retained_v2",
                        "reason": "roughness_guard_candidate_not_smoother",
                        "canonical_data_modified": False,
                    }
                )
                continue
            applied_curve_count += 1
            target_count += len(pending_audits)
            changed_count += curve_changed_count
            max_change = max(max_change, curve_max_change)
            for index, audit in pending_audits:
                audit["local_residual_after"] = local_log_residual(rows, index, "clean_value")
                point_audit.append(audit)
            curve_audit.append(
                {
                    "window_id": window["window_id"],
                    "mode_key": window["mode_key"],
                    "plot_panel": window["plot_panel"],
                    "plot_series": window["plot_series"],
                    "observable": observable,
                    "anchor_left_xi": V2.canonical_xi(window["anchor_left_xi"]),
                    "anchor_right_xi": V2.canonical_xi(window["anchor_right_xi"]),
                    "interior_point_count": right_index - left_index - 1,
                    "source_roughness": source_roughness,
                    "candidate_roughness": candidate_roughness,
                    "roughness_ratio": roughness_ratio,
                    "changed_point_count": curve_changed_count,
                    "max_abs_relative_change": curve_max_change,
                    "action": "applied",
                    "reason": "roughness_guard_passed",
                    "canonical_data_modified": False,
                }
            )
        window_points = point_audit[audit_start:]
        residual_pairs = [
            (abs(float(row["local_residual_before"])), abs(float(row["local_residual_after"])))
            for row in window_points
            if math.isfinite(float(row["local_residual_before"]))
            and math.isfinite(float(row["local_residual_after"]))
        ]
        action = "applied" if not gate_errors else "applied_with_errors" if target_count else "retained_v2"
        window_audit.append(
            {
                "window_id": window["window_id"],
                "scope": window["scope"],
                "mode_key": window["mode_key"],
                "plot_panel": window["plot_panel"],
                "plot_series": window["plot_series"],
                "anchor_left_xi": V2.canonical_xi(window["anchor_left_xi"]),
                "anchor_right_xi": V2.canonical_xi(window["anchor_right_xi"]),
                "strength": window["strength"],
                "value_space": window["value_space"],
                "phase_policy": window["phase_policy"],
                "observables": ";".join(window["observables_tuple"]),
                "target_point_count": target_count,
                "changed_point_count": changed_count,
                "applied_curve_count": applied_curve_count,
                "retained_curve_count": retained_curve_count,
                "retained_curve_ids": ";".join(retained_curve_ids),
                "max_abs_relative_change": max_change,
                "mean_abs_local_residual_before": (
                    sum(before for before, _ in residual_pairs) / len(residual_pairs)
                    if residual_pairs else math.nan
                ),
                "mean_abs_local_residual_after": (
                    sum(after for _, after in residual_pairs) / len(residual_pairs)
                    if residual_pairs else math.nan
                ),
                "improved_point_count": sum(after < before for before, after in residual_pairs),
                "worsened_point_count": sum(after > before for before, after in residual_pairs),
                "action": action,
                "gate_errors": "; ".join(gate_errors),
                "recipe_source": relpath(RECIPE),
                "recipe_sha256": recipe_hash,
                "canonical_data_modified": False,
            }
        )
    if len(modified_keys) != len(point_audit):
        raise AssertionError("point audit key count does not match modified point count")
    return points, point_audit, window_audit, curve_audit


def build_curve_index(points: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str, str, str], list[dict[str, Any]]] = defaultdict(list)
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
                "v2_replacement_count": sum(row["display_status"] != "raw" for row in rows),
                "residual_smoothing_count": sum(row["display_status"] == "residual_smoothed" for row in rows),
                "xi_min": min(finite(row, "xi") for row in rows),
                "xi_max": max(finite(row, "xi") for row in rows),
                "canonical_data_modified": False,
            }
        )
    return out


def render_audit_figures(
    source_points: list[dict[str, Any]],
    final_points: list[dict[str, Any]],
    point_audit: list[dict[str, Any]],
) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    source_index = {point_key(row): row for row in source_points}
    final_index = {point_key(row): row for row in final_points}
    selected = [
        ("mode_a", "muB900.0", "alpha1.0", "tau_u"),
        ("mode_a", "muB900.0", "alpha1.0", "tau_ubar"),
        ("mode_a", "muB900.0", "alpha1.1", "tau_ubar"),
        ("mode_b", "T200.0", "muB0.0", "tau_dbar"),
        ("mode_b", "T200.0", "muB450.0", "tau_dbar"),
        ("mode_b", "T200.0", "muB900.0", "tau_sbar"),
    ]
    paths: list[Path] = []
    fig, axes = plt.subplots(3, 2, figsize=(10.5, 11.0), squeeze=False)
    for axis, key in zip(axes.flat, selected):
        rows = sorted(
            [row for row in final_points if curve_key(row) == key],
            key=lambda row: finite(row, "xi"),
        )
        source_rows = [source_index[point_key(row)] for row in rows]
        axis.plot(
            [finite(row, "xi") for row in source_rows],
            [finite(row, "clean_value") for row in source_rows],
            color="#999999",
            linestyle="--",
            linewidth=1.0,
            label="v2 display",
        )
        axis.plot(
            [finite(row, "xi") for row in rows],
            [finite(row, "clean_value") for row in rows],
            color="#2255AA",
            linewidth=1.35,
            label="residual-smoothed",
        )
        axis.set_title(f"{key[0]} / {key[1]} / {key[2]} / {key[3]}", fontsize=9)
        axis.set_xlabel(r"$\xi$")
        axis.set_ylabel(V2.OBSERVABLE_LABELS[key[3]])
        axis.grid(alpha=0.18)
        axis.legend(fontsize=7, loc="best")
    fig.suptitle("publication_clean_v2 residual smoothing: before / after", fontsize=12)
    fig.tight_layout()
    overview = FIGURE_DIR / "residual_smoothing_overview.png"
    overview.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(overview, dpi=450, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    paths.append(overview)

    fig, axes = plt.subplots(3, 2, figsize=(10.5, 11.0), squeeze=False)
    audit_by_key: dict[tuple[str, str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in point_audit:
        audit_by_key[(row["mode_key"], row["plot_panel"], row["plot_series"], row["observable"])].append(row)
    for axis, key in zip(axes.flat, selected):
        rows = sorted(audit_by_key.get(key, []), key=lambda row: float(row["xi"]))
        if rows:
            axis.axhline(0.0, color="#444444", linewidth=0.7)
            axis.plot(
                [float(row["xi"]) for row in rows],
                [100.0 * float(row["relative_change_from_v2"]) for row in rows],
                marker="o",
                color="#CC3311",
                linewidth=1.1,
            )
        axis.set_title(f"{key[0]} / {key[1]} / {key[2]} / {key[3]}", fontsize=9)
        axis.set_xlabel(r"$\xi$")
        axis.set_ylabel("change from v2 (%)")
        axis.grid(alpha=0.18)
    fig.suptitle("display-only residual smoothing adjustment", fontsize=12)
    fig.tight_layout()
    delta_path = FIGURE_DIR / "residual_smoothing_delta.png"
    fig.savefig(delta_path, dpi=450, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    paths.append(delta_path)
    return paths


def point_residual_pairs(point_audit: list[dict[str, Any]]) -> list[tuple[float, float]]:
    pairs = [
        (abs(float(row["local_residual_before"])), abs(float(row["local_residual_after"])))
        for row in point_audit
    ]
    return [(before, after) for before, after in pairs if math.isfinite(before) and math.isfinite(after)]


def render_readme(
    source_manifest: dict[str, Any],
    windows: list[dict[str, Any]],
    window_audit: list[dict[str, Any]],
    point_audit: list[dict[str, Any]],
    curve_audit: list[dict[str, Any]],
    figure_paths: list[Path],
) -> str:
    applied = sum(row["action"] in {"applied", "applied_with_errors"} for row in window_audit)
    retained = sum(row["action"] == "retained_v2" for row in window_audit)
    applied_curves = sum(row["action"] == "applied" for row in curve_audit)
    retained_curves = sum(row["action"] == "retained_v2" for row in curve_audit)
    residual_pairs = point_residual_pairs(point_audit)
    mean_before = (
        sum(before for before, _ in residual_pairs) / len(residual_pairs)
        if residual_pairs else math.nan
    )
    mean_after = (
        sum(after for _, after in residual_pairs) / len(residual_pairs)
        if residual_pairs else math.nan
    )
    improved = sum(after < before for before, after in residual_pairs)
    worsened = sum(after > before for before, after in residual_pairs)
    audit_figure_names = {"residual_smoothing_overview.png", "residual_smoothing_delta.png"}
    audit_figure_count = sum(path.name in audit_figure_names for path in figure_paths)
    publication_figure_count = len(figure_paths) - audit_figure_count
    retained_curve_rows = [
        row for row in curve_audit if row["action"] == "retained_v2"
    ]
    high_change_rows = sorted(
        [
            row
            for row in curve_audit
            if row["action"] == "applied"
            and float(row["max_abs_relative_change"]) >= 0.05
        ],
        key=lambda row: float(row["max_abs_relative_change"]),
        reverse=True,
    )
    high_change_table = (
        "| window_id | observable | max relative change | roughness ratio |\n"
        "| --- | --- | ---: | ---: |\n"
        + "\n".join(
            f"| {row['window_id']} | {row['observable']} | "
            f"{float(row['max_abs_relative_change']):.4g} | "
            f"{float(row['roughness_ratio']):.4g} |"
            for row in high_change_rows
        )
        if high_change_rows
        else "没有达到 5% 相对改动阈值的 applied 曲线。"
    )
    retained_curve_text = (
        "; ".join(
            f"{row['window_id']}:{row['observable']}"
            for row in retained_curve_rows
        )
        if retained_curve_rows
        else "无"
    )
    return f"""# Issue #130 RS `publication_clean_v2_residual_smoothed`

## 目的与边界

本包是现有 `publication_clean_v2` 的第二层 display-only 派生。它只消费 v2 的
`publication_clean_points.csv`，不调用 equilibrium/transport solver，不修改 raw CSV、
production registry、v2 原目录或 `main`。

本轮处理的是 v2 仍可见的局部显示结构：小分母敏感窗口的三点肩部、mode-A
`muB=900` 的 tau 局部斜率鼓包，以及 `alpha_T=1.1` tau 曲线的两段快速上升。
这不是传播子有限宽度、分母截断或新的物理计算。

## 显示算法

对每个显式 recipe 窗口，保留左右 anchor 的 v2 display 值；在正值曲线的
`log(y)` 空间构造单调 Hermite 趋势，再按 recipe 中的 `strength` 与原 v2 值作
几何混合。`strength=1` 用于孤立小突变，较小强度用于连续 branch response。
候选曲线还必须通过逐曲线 roughness gate：若窗口内最大的相邻 log-y 斜率跳变
没有下降，则该 observable 保留 v2 值，避免对原本平滑的曲线制造新拐点。

实际渲染的一阶 gap 是硬保护区：窗口不得进入 gap 内部或跨越 gap。branch-local
窗口只允许在同一个 `phase_curr` 分支内处理，因此不会把两个相态拼成一条曲线。

## 输入 provenance

- source artifact: `{relpath(SOURCE_DIR)}`
- source manifest SHA256: `{sha256_file(SOURCE_MANIFEST)}`
- source v2 point rows: `{source_manifest['derived_counts']['publication_clean_point_rows']}`
- recipe: `{relpath(RECIPE)}`
- recipe SHA256: `{sha256_file(RECIPE)}`
- derived solver called: `false`
- canonical/raw data modified: `false`
- manuscript eligible: `false`

## 结果摘要

- recipe windows: {len(windows)}
- applied windows: {applied}
- retained v2 windows: {retained}
- applied observable curves: {applied_curves}
- roughness-guard retained curves: {retained_curves}
- residual-smoothed point rows: {len(point_audit)}
- mean absolute local log residual: `{mean_before:.6g}` -> `{mean_after:.6g}`
- point residual comparison: {improved} improved / {worsened} worsened
- publication figures: {publication_figure_count}
- audit figures: {audit_figure_count}

`tables/publication_clean_points.csv` 同时保留 `raw_value`、`v2_clean_value` 和最终
`clean_value`。改变点还记录 `smoothing_window`、锚点、拟合方法、phase gate 状态及
`local_residual_before/after`；所有改变点都在 `tables/residual_point_audit.csv` 中逐点记录，
逐曲线 roughness 决策见 `tables/residual_curve_audit.csv`。

## 当前处理决策与作者审阅项

- roughness gate 保留 v2 的曲线数：{len(retained_curve_rows)}；这些曲线不应继续套用同一窗口。当前保留项为：`{retained_curve_text}`。
- applied 曲线中相对改动达到 5% 的项目列在下表；它们只适合作者确认，不应自动视为投稿最终值：

{high_change_table}

- `muB=900, alpha_T=1.0/1.1` 的反夸克 tau 是连续 branch response，处理目标只是降低局部斜率尖锐度；若作者不接受约 5--10% 的显示改动，应降低对应 recipe strength 或直接保留 v2。
- `mode B / T=200, muB=900 / tau_sbar` 的窗口虽可降低图面肩部，但 high-rate convergence gate 仍为空；在补 gate 前只能保留为 channel-rate candidate，不能以平滑图面替代机制证据。
- `tables/review_adjustment_map.csv` 中的 author-review 插值是输入侧审阅记录，不代表本层已自动应用；任何采用都应重新生成并核对 manifest。
- 一阶 gap、端点和跨分支连接继续禁止填补；连续宽响应不应通过扩大窗口被抹平。

## 解释边界

1. 这些值仍是 display-only 派生值，不是 solver 重算、收敛证明或物理正则化。
2. 小分母窗口的机制证据仍应引用原有 mechanism audit；本包只改变论文图的显示形状。
   特别是 `mode B / T=200 MeV / muB=900 MeV / tau_sbar` 仍是尚未完成
   high-rate convergence 的 channel-rate candidate，不能因图形变平而升级为已证实机制。
3. `muB=900, alpha_T=1` 的 `xi=-0.003/+0.003` 和 mode-B `T=120, muB=900`
   的 `xi=-0.13/-0.12` 一阶端点继续使用 v2 raw endpoint，gap 不被填充。
4. 任何定量峰值、临界行为或输运系数精确拟合仍必须使用 raw/v2 evidence，不能用
   本层的平滑值替代。

## 复现

```powershell
python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v2_residual_smoothed.py
python -m pytest tests/unit/python/test_phase_guided_publication_clean_v2_residual_smoothed.py
```
"""


def claim_ledger(
    source_manifest: dict[str, Any],
    windows: list[dict[str, Any]],
    point_audit: list[dict[str, Any]],
    curve_audit: list[dict[str, Any]],
) -> list[dict[str, str]]:
    return [
        {
            "claim_id": "PC-V2R-001",
            "status": "supported_with_scope_limit",
            "claim_zh": "本包只消费 publication_clean_v2 的 display table，raw_value 与 v2_clean_value 均保留。",
            "evidence": "manifest.json; tables/publication_clean_points.csv; tables/residual_point_audit.csv",
            "scope_limit": "不是新的 solver 结果或 production 数据。",
        },
        {
            "claim_id": "PC-V2R-002",
            "status": "supported_with_scope_limit",
            "claim_zh": f"{len(point_audit)} 个点按显式窗口作局部 display-only residual smoothing。",
            "evidence": "tables/residual_smoothing_window_audit.csv; tables/residual_curve_audit.csv; tables/residual_point_audit.csv",
            "scope_limit": "平滑参数是显示选择，不代表物理宽度或传播子正则化。",
        },
        {
            "claim_id": "PC-V2R-003",
            "status": "supported",
            "claim_zh": "实际一阶 gap 不被跨越；端点继续使用 v2 endpoint 值。",
            "evidence": "tables/boundary_gap_map.csv; tables/publication_clean_points.csv; figures/plot_manifest.json",
            "scope_limit": "gap 语义沿用 v2，不构造区间内输运值。",
        },
        {
            "claim_id": "PC-V2R-004",
            "status": "author_check",
            "claim_zh": "平滑后的曲线适合作者审阅和图形排版，不应替代 raw/v2 曲线作定量机制结论。",
            "evidence": "figures/residual_smoothing_overview.png; figures/residual_smoothing_delta.png",
            "scope_limit": "正式论文使用仍需作者确认显示强度与窗口清单。",
        },
    ]


def main() -> None:
    source_points, source_manifest = load_source_points()
    windows = load_recipe()
    loaded, inventory = V2.V1.load_inputs(OBSERVABLES)
    gaps = gap_records(loaded)
    raw_phase = build_raw_phase_index(loaded)

    expected_inputs = source_manifest.get("source_inputs", [])
    current_by_mode = {row["mode_key"]: row for row in inventory}
    for expected in expected_inputs:
        current = current_by_mode.get(expected["mode_key"])
        if current is None or current["scan_sha256"] != expected["scan_sha256"] or current["diagnostics_sha256"] != expected["diagnostics_sha256"]:
            raise ValueError(f"source input provenance changed for {expected['mode_key']}")

    final_points = [dict(row) for row in source_points]
    final_points, point_audit, window_audit, curve_audit = apply_residual_smoothing(
        final_points, windows, raw_phase, gaps
    )
    curves = build_curve_index(final_points)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    for source_table in SOURCE_TABLE_DIR.glob("*.csv"):
        shutil.copy2(source_table, TABLE_DIR / source_table.name)

    V2.FIGURE_DIR = FIGURE_DIR
    figure_paths, figure_specs = V2.render_figures(final_points, gaps, OBSERVABLES)
    audit_paths = render_audit_figures(source_points, final_points, point_audit)
    all_figure_paths = [*figure_paths, *audit_paths]

    point_fields = [
        "mode_key", "mode", "plot_panel", "plot_series", "plot_series_label", "T_MeV", "muB_MeV", "xi",
        "observable", "raw_value", "v2_clean_value", "clean_value", "display_status", "value_source",
        "phase_structure", "phase_reference_kind", "quality_flag", "quality_reason", "run_id",
        "canonical_data_modified", "smoothing_window", "smoothing_anchor_points", "smoothing_fit_method",
        "protected_by_phase_gate",
    ]
    write_csv(TABLE_DIR / "publication_clean_points.csv", final_points, point_fields)
    curve_fields = [
        "mode_key", "plot_panel", "plot_series", "observable", "point_count", "v2_replacement_count",
        "residual_smoothing_count", "xi_min", "xi_max", "canonical_data_modified",
    ]
    write_csv(TABLE_DIR / "curve_index.csv", curves, curve_fields)
    write_csv(
        TABLE_DIR / "residual_point_audit.csv",
        point_audit,
        [
            "window_id", "mode_key", "plot_panel", "plot_series", "observable", "xi", "v2_clean_value",
            "trend_value", "final_clean_value", "relative_change_from_v2", "strength", "value_space",
            "phase_policy", "anchor_points", "fit_method", "local_residual_before", "local_residual_after",
            "protected_by_phase_gate", "canonical_data_modified",
        ],
    )
    write_csv(
        TABLE_DIR / "residual_smoothing_window_audit.csv",
        window_audit,
        [
            "window_id", "scope", "mode_key", "plot_panel", "plot_series", "anchor_left_xi",
            "anchor_right_xi", "strength", "value_space", "phase_policy", "observables",
            "target_point_count", "changed_point_count", "applied_curve_count", "retained_curve_count",
            "retained_curve_ids", "max_abs_relative_change", "action", "gate_errors",
            "mean_abs_local_residual_before", "mean_abs_local_residual_after", "improved_point_count",
            "worsened_point_count", "recipe_source", "recipe_sha256", "canonical_data_modified",
        ],
    )
    write_csv(
        TABLE_DIR / "residual_curve_audit.csv",
        curve_audit,
        [
            "window_id", "mode_key", "plot_panel", "plot_series", "observable", "anchor_left_xi",
            "anchor_right_xi", "interior_point_count", "source_roughness", "candidate_roughness",
            "roughness_ratio", "changed_point_count", "max_abs_relative_change", "action", "reason",
            "canonical_data_modified",
        ],
    )
    write_csv(
        TABLE_DIR / "claim_ledger.csv",
        claim_ledger(source_manifest, windows, point_audit, curve_audit),
        ["claim_id", "status", "claim_zh", "evidence", "scope_limit"],
    )

    generator_path = Path(__file__).resolve()
    figure_assets = []
    for spec in figure_specs:
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
    for path in audit_paths:
        figure_assets.append(
            {
                "path": relpath(path),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
                "kind": "before_after_audit",
            }
        )
    axis_counts: dict[str, int] = defaultdict(int)
    for spec in figure_specs:
        axis_counts[str(spec["axis_scale"])] += 1
    plot_manifest = {
        "schema": "phase_guided_transport_publication_clean_plot_manifest_v2_residual_smoothed",
        "case": V2.CURRENT_CASE,
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "base_git_commit": V2.git_head(),
        "generator": relpath(generator_path),
        "generator_sha256": sha256_file(generator_path),
        "source_v2_manifest_sha256": sha256_file(SOURCE_MANIFEST),
        "observables": OBSERVABLES,
        "boundary_gap_count": len(gaps),
        "axis_scale_counts": dict(sorted(axis_counts.items())),
        "manuscript_eligible": False,
        "canonical_data_modified": False,
        "rendering_semantics": "v2 clean display values with explicit local residual smoothing; first-order gaps remain split and endpoint values remain v2 values; no raw mutation",
        "figures": figure_assets,
    }
    write_json(FIGURE_DIR / "plot_manifest.json", plot_manifest)

    readme_path = OUT_DIR / "README.md"
    readme_path.write_text(
        render_readme(source_manifest, windows, window_audit, point_audit, curve_audit, all_figure_paths),
        encoding="utf-8",
    )
    output_paths = [
        readme_path,
        *sorted(TABLE_DIR.glob("*.csv")),
        FIGURE_DIR / "plot_manifest.json",
        *all_figure_paths,
    ]
    residual_pairs = point_residual_pairs(point_audit)
    manifest = {
        "schema": "phase_guided_transport_publication_clean_manifest_v2_residual_smoothed",
        "case": V2.CURRENT_CASE,
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "base_git_commit": V2.git_head(),
        "generator": relpath(generator_path),
        "generator_sha256": sha256_file(generator_path),
        "source_v2_artifact": relpath(SOURCE_DIR),
        "source_v2_manifest_sha256": sha256_file(SOURCE_MANIFEST),
        "recipe_path": relpath(RECIPE),
        "recipe_sha256": sha256_file(RECIPE),
        "status": "derived_author_review_required",
        "manuscript_eligible": False,
        "canonical_data_modified": False,
        "production_write": False,
        "solver_called": False,
        "source_solver_called": True,
        "source_inputs": inventory,
        "derived_counts": {
            "recipe_window_rows": len(windows),
            "window_audit_rows": len(window_audit),
            "window_applied_rows": sum(row["action"] in {"applied", "applied_with_errors"} for row in window_audit),
            "curve_audit_rows": len(curve_audit),
            "curve_applied_rows": sum(row["action"] == "applied" for row in curve_audit),
            "curve_roughness_retained_rows": sum(row["action"] == "retained_v2" for row in curve_audit),
            "residual_point_rows": len(point_audit),
            "publication_clean_point_rows": len(final_points),
            "curve_rows": len(curves),
            "figure_count": len(all_figure_paths),
            "publication_figure_count": len(figure_paths),
            "audit_figure_count": len(audit_paths),
        },
        "residual_audit": {
            "metric": "three-point log residual",
            "finite_point_rows": len(residual_pairs),
            "mean_abs_before": (
                sum(before for before, _ in residual_pairs) / len(residual_pairs)
                if residual_pairs else None
            ),
            "mean_abs_after": (
                sum(after for _, after in residual_pairs) / len(residual_pairs)
                if residual_pairs else None
            ),
            "improved_point_rows": sum(after < before for before, after in residual_pairs),
            "worsened_point_rows": sum(after > before for before, after in residual_pairs),
        },
        "known_boundaries": [
            "existing publication_clean_v2 remains unchanged",
            "mode_a muB=900 alpha_T=1 xi=-0.003/+0.003 gap remains unfilled",
            "mode_b T=120 muB=900 xi=-0.13/-0.12 gap remains unfilled",
            "residual smoothing is display-only and not physical denominator regularization",
            "raw_value and v2_clean_value remain available for quantitative use",
        ],
        "outputs": [
            {"path": relpath(path), "sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in output_paths
        ],
    }
    write_json(OUT_DIR / "manifest.json", manifest)
    print(
        json.dumps(
            {
                "output": relpath(OUT_DIR),
                "manifest": relpath(OUT_DIR / "manifest.json"),
                "publication_figures": len(figure_paths),
                "audit_figures": len(audit_paths),
                "residual_points": len(point_audit),
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
