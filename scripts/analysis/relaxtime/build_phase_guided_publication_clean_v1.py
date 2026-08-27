#!/usr/bin/env python3
"""Build the solver-free publication-clean RS transport layer.

This command derives a display/input layer from the two author-accepted
``prod_v2`` raw result trees.  It never edits a raw CSV, the production
registry, or a canonical figure.  The existing pole-sensitive display recipe
is used only as an audited list of point keys; replacement values are
recomputed from the current raw neighbours so that the derived layer cannot
silently carry values from the superseded ``prod_v1`` tree.

The output is intentionally not manuscript eligible.  In particular, the
known mixed equilibrium/bulk branch issue and the direct-coexistence
``xi=+/-0.003`` contract remain visible in the provenance tables.
"""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import math
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[3]
RESULT_ROOT = ROOT / "data" / "outputs" / "results" / "relaxtime" / "transport" / "phase_guided"
OUT_DIR = (
    ROOT
    / "docs"
    / "analysis"
    / "relaxtime"
    / "phase_guided_transport"
    / "phase_guided_transport_publication_clean_v1"
)
TABLE_DIR = OUT_DIR / "tables"
FIGURE_DIR = OUT_DIR / "figures"

CURRENT_CASE = "first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2"
CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
WORKFLOW_HEAD_SHA = "22874505877491754eed27519ad8a7b871c82571"
RECIPE_DIR = (
    ROOT
    / "docs"
    / "analysis"
    / "relaxtime"
    / "phase_guided_transport"
    / "phase_guided_transport_v2_pole_sensitive_rendering"
    / "tables"
)
REPLACEMENT_RECIPE = RECIPE_DIR / "paper_display_replacements.csv"
MARKER_RECIPE = RECIPE_DIR / "paper_first_order_markers.csv"
BULK_AUDIT = RECIPE_DIR / "bulk_derivative_branch_audit.csv"
CEP_BOUNDARY = (
    ROOT
    / "data"
    / "reference"
    / "pnjl"
    / "issue130_phase_reference_v1"
    / "strict"
    / "tables"
    / "cep_boundary_strict_reference_v1.csv"
)

DISPLAY_FIELDS = ["eta_over_s", "zeta_over_s", "sigma_over_T"]
# These are explicit author-requested display-layer repairs for the current
# prod_v2 candidate.  They are deliberately kept outside the historical
# pole-sensitive recipe so that an audit can distinguish inherited choices
# from this review round.  The raw result tree is never changed.
REVIEW_ADJUSTMENTS = (
    {
        "window_id": "author_review_mode_b_T200p0_muB900p0_xim0p10",
        "mode_key": "mode_b",
        "plot_panel": "T200.0",
        "plot_series": "muB900.0",
        "xi": "-0.10",
        "left_xi": "-0.11",
        "right_xi": "-0.09",
        "observables": DISPLAY_FIELDS[:2],
        "reason": "author-requested local spike smoothing; current raw point is a small downward residual against adjacent xi samples",
    },
    {
        "window_id": "author_review_mode_b_T200p0_muB0p0_xip0p36",
        "mode_key": "mode_b",
        "plot_panel": "T200.0",
        "plot_series": "muB0.0",
        "xi": "0.36",
        "left_xi": "0.35",
        "right_xi": "0.37",
        "observables": DISPLAY_FIELDS[:2],
        "reason": "author-requested local spike smoothing; current raw point is a small downward residual against adjacent xi samples",
    },
)

# Publication markers are display-only midpoints of an author-confirmed
# phase-switch interval.  They are deliberately explicit rather than inferred
# from every noisy phase-label flip in the transport table: the latter contains
# diagnostic classification chatter that is not a physical CEP candidate.
# Endpoint rows remain raw evidence; the midpoint value is computed below and
# is never written into the clean point table.
PHASE_SWITCH_MARKERS = (
    {
        "window_id": "publication_mode_a_muB900p0_alpha1p0_direct_coexistence_midpoint",
        "mode_key": "mode_a",
        "plot_panel": "muB900.0",
        "plot_series": "alpha1.0",
        "left_xi": "-0.003",
        "right_xi": "0.003",
        "marker_basis": "direct_coexistence_interval",
        "reason": "direct-coexistence endpoints bracket the isotropic xi=0 phase-switch display point",
    },
    {
        "window_id": "publication_mode_b_T120p0_muB900p0_phase_switch_midpoint",
        "mode_key": "mode_b",
        "plot_panel": "T120.0",
        "plot_series": "muB900.0",
        "left_xi": "-0.14",
        "right_xi": "-0.13",
        "marker_basis": "phase_label_switch",
        "reason": "author-confirmed current prod_v2 crossover -> first_order label switch bracketing the visible transition",
    },
)
MODE_CONFIG = {
    "mode_a": {
        "mode": "mode_a_fixed_muB_phase_scaled",
        "label": "mode A: fixed muB, phase-scaled T",
    },
    "mode_b": {
        "mode": "mode_b_fixed_T_sparse_muB",
        "label": "mode B: fixed T, sparse muB",
    },
}


def relpath(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_head() -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()


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


def parse_finite(row: dict[str, str], field: str) -> float:
    value = float(row[field])
    if not math.isfinite(value):
        raise ValueError(f"non-finite {field}: {row[field]}")
    return value


def canonical_xi(value: str | float) -> str:
    return f"{float(value):.10f}"


def point_key(row: dict[str, str]) -> tuple[str, str, str]:
    return row["plot_panel"], row["plot_series"], canonical_xi(row["xi"])


def curve_key(mode_key: str, row: dict[str, str]) -> tuple[str, str, str, str]:
    panel, series, xi = point_key(row)
    return mode_key, panel, series, xi


def case_paths(mode_key: str) -> dict[str, Path]:
    mode = MODE_CONFIG[mode_key]["mode"]
    result_dir = RESULT_ROOT / mode / CURRENT_CASE
    return {
        "result_dir": result_dir,
        "scan": result_dir / "phase_guided_transport_scan.csv",
        "diagnostics": result_dir / "channel_diagnostics.csv",
        "failed": result_dir / "failed_points.csv",
        "manifest": result_dir / "manifest.json",
        "effective_config": result_dir / "effective_config.json",
    }


def validate_scan(
    mode_key: str, paths: dict[str, Path]
) -> tuple[list[dict[str, str]], dict[tuple[str, str, str], dict[str, str]], dict[tuple[str, str], list[dict[str, str]]], dict[str, Any]]:
    rows = read_csv(paths["scan"])
    if not rows:
        raise ValueError(f"{mode_key}: empty scan")
    required = {
        "T_MeV",
        "muB_MeV",
        "xi",
        "mode",
        "phase_reference_kind",
        "phase_structure",
        "plot_panel",
        "plot_series",
        "plot_series_label",
        "converged",
        "quality_flag",
        "quality_reason",
        *DISPLAY_FIELDS,
    }
    missing = sorted(required - set(rows[0]))
    if missing:
        raise ValueError(f"{mode_key}: missing scan fields {missing}")
    index: dict[tuple[str, str, str], dict[str, str]] = {}
    curves: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        for field in ["T_MeV", "muB_MeV", "xi", *DISPLAY_FIELDS]:
            parse_finite(row, field)
        if row["converged"].lower() != "true":
            raise ValueError(f"{mode_key}: non-converged scan row {point_key(row)}")
        key = point_key(row)
        if key in index:
            raise ValueError(f"{mode_key}: duplicate scan key {key}")
        index[key] = row
        curves[(row["plot_panel"], row["plot_series"])].append(row)
    for curve in curves.values():
        curve.sort(key=lambda row: parse_finite(row, "xi"))
    diagnostic_rows = read_csv(paths["diagnostics"])
    diagnostic_required = {"density", "rate", "contribution", "total", "tau_inv_species"}
    if diagnostic_rows and not diagnostic_required.issubset(diagnostic_rows[0]):
        raise ValueError(f"{mode_key}: diagnostic fields are incomplete")
    for row in diagnostic_rows:
        for field in diagnostic_required:
            value = parse_finite(row, field)
            if value < 0:
                raise ValueError(f"{mode_key}: negative diagnostic {field}: {value}")
    failed_rows = read_csv(paths["failed"])
    if failed_rows:
        raise ValueError(f"{mode_key}: failed_points.csv contains {len(failed_rows)} rows")
    manifest = json.loads(paths["manifest"].read_text(encoding="utf-8"))
    if manifest.get("calculation_sha") != CALCULATION_SHA:
        raise ValueError(f"{mode_key}: calculation SHA mismatch")
    if manifest.get("workflow_head_sha") != WORKFLOW_HEAD_SHA:
        raise ValueError(f"{mode_key}: workflow head SHA mismatch")
    if manifest.get("production_write") is not False:
        raise ValueError(f"{mode_key}: source manifest permits production write")
    if manifest.get("source_solver_called") is not True:
        raise ValueError(f"{mode_key}: source solver provenance missing")
    hashes = manifest.get("hashes", {})
    expected_scan_hash = hashes.get("phase_guided_transport_scan.csv")
    if expected_scan_hash and expected_scan_hash != sha256_file(paths["scan"]):
        raise ValueError(f"{mode_key}: scan hash disagrees with source manifest")
    expected_diag_hash = hashes.get("channel_diagnostics.csv")
    if expected_diag_hash and expected_diag_hash != sha256_file(paths["diagnostics"]):
        raise ValueError(f"{mode_key}: diagnostic hash disagrees with source manifest")
    expected_failed_hash = hashes.get("failed_points.csv")
    if expected_failed_hash and expected_failed_hash != sha256_file(paths["failed"]):
        raise ValueError(f"{mode_key}: failed-point hash disagrees with source manifest")
    scan_summary = manifest.get("scan_summary", {})
    if scan_summary.get("rows") not in (None, len(rows)):
        raise ValueError(f"{mode_key}: source scan row count disagrees with manifest")
    diagnostic_summary = manifest.get("diagnostic_summary", {})
    if diagnostic_summary.get("rows") not in (None, len(diagnostic_rows)):
        raise ValueError(f"{mode_key}: diagnostic row count disagrees with manifest")
    if diagnostic_summary.get("nonfinite", 0) or diagnostic_summary.get("duplicate_keys", 0):
        raise ValueError(f"{mode_key}: source diagnostics are not clean")
    return rows, index, dict(curves), manifest


def load_inputs() -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]]]:
    loaded: dict[str, dict[str, Any]] = {}
    inventory: list[dict[str, Any]] = []
    for mode_key in MODE_CONFIG:
        paths = case_paths(mode_key)
        for name in ("scan", "diagnostics", "failed", "manifest", "effective_config"):
            if not paths[name].exists():
                raise FileNotFoundError(f"missing {mode_key} input {paths[name]}")
        rows, index, curves, manifest = validate_scan(mode_key, paths)
        loaded[mode_key] = {
            "paths": paths,
            "rows": rows,
            "index": index,
            "curves": curves,
            "manifest": manifest,
        }
        inventory.append(
            {
                "mode_key": mode_key,
                "mode": MODE_CONFIG[mode_key]["mode"],
                "scan_rows": len(rows),
                "diagnostic_rows": manifest.get("diagnostic_summary", {}).get("rows", ""),
                "failed_rows": 0,
                "xi_count": len({canonical_xi(row["xi"]) for row in rows}),
                "scan_sha256": sha256_file(paths["scan"]),
                "diagnostics_sha256": sha256_file(paths["diagnostics"]),
                "failed_sha256": sha256_file(paths["failed"]),
                "manifest_sha256": sha256_file(paths["manifest"]),
                "effective_config_sha256": sha256_file(paths["effective_config"]),
                "calculation_sha": manifest["calculation_sha"],
                "workflow_head_sha": manifest["workflow_head_sha"],
                "source_solver_called": manifest["source_solver_called"],
                "derived_solver_called": False,
            }
        )
    return loaded, inventory


def load_recipe() -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    replacements = read_csv(REPLACEMENT_RECIPE)
    markers = read_csv(MARKER_RECIPE)
    if len(replacements) != 19:
        raise ValueError(f"display recipe rows {len(replacements)} != 19")
    if len(markers) != 6:
        raise ValueError(f"marker recipe rows {len(markers)} != 6")
    if any(row["mode_key"] not in MODE_CONFIG for row in [*replacements, *markers]):
        raise ValueError("recipe contains an unknown mode")
    return replacements, markers


def load_cep_boundary() -> list[dict[str, str]]:
    rows = read_csv(CEP_BOUNDARY)
    required = {
        "xi",
        "mu_CEP_proxy_MeV",
        "T_low_MeV",
        "T_high_MeV",
        "T_midpoint_MeV",
        "layer",
        "status",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValueError("strict CEP boundary table is empty or missing required fields")
    if any(row.get("layer") != "strict_reference_v1" for row in rows):
        raise ValueError("CEP boundary input contains a non-strict layer")
    return rows


def current_row(
    loaded: dict[str, dict[str, Any]], mode_key: str, panel: str, series: str, xi: str | float
) -> dict[str, str] | None:
    return loaded[mode_key]["index"].get((panel, series, canonical_xi(xi)))


def build_replacement_map(
    loaded: dict[str, dict[str, Any]], recipes: list[dict[str, str]]
) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for recipe in recipes:
        mode_key = recipe["mode_key"]
        target = current_row(loaded, mode_key, recipe["plot_panel"], recipe["plot_series"], recipe["xi"])
        left = current_row(loaded, mode_key, recipe["plot_panel"], recipe["plot_series"], recipe["left_xi"])
        right = current_row(loaded, mode_key, recipe["plot_panel"], recipe["plot_series"], recipe["right_xi"])
        if target is None or left is None or right is None:
            raise ValueError(
                f"recipe key missing from current {mode_key}: {recipe['plot_panel']} {recipe['plot_series']} {recipe['xi']}"
            )
        target_xi = parse_finite(target, "xi")
        left_xi = parse_finite(left, "xi")
        right_xi = parse_finite(right, "xi")
        if not left_xi < target_xi < right_xi:
            raise ValueError(f"invalid replacement bracket for {recipe['window_id']}")
        weight = (target_xi - left_xi) / (right_xi - left_xi)
        raw_value = parse_finite(target, recipe["observable"])
        left_value = parse_finite(left, recipe["observable"])
        right_value = parse_finite(right, recipe["observable"])
        derived_value = left_value + weight * (right_value - left_value)
        out.append(
            {
                "window_id": recipe["window_id"],
                "mode_key": mode_key,
                "plot_panel": recipe["plot_panel"],
                "plot_series": recipe["plot_series"],
                "observable": recipe["observable"],
                "xi": canonical_xi(target["xi"]),
                "raw_production_value_current": raw_value,
                "recipe_raw_production_value": recipe["raw_production_value"],
                "recipe_display_value": recipe["paper_display_value"],
                "left_xi": canonical_xi(left["xi"]),
                "left_value_current": left_value,
                "right_xi": canonical_xi(right["xi"]),
                "right_value_current": right_value,
                "derived_display_value": derived_value,
                "replacement_method": "linear interpolation between current prod_v2 neighbours",
                "recipe_source": relpath(REPLACEMENT_RECIPE),
                "recipe_source_sha256": sha256_file(REPLACEMENT_RECIPE),
                "canonical_data_modified": False,
                "display_status": "interpolated_noncertified",
                "adjustment_type": "inherited_recipe",
                "adjustment_reason": "historical pole-sensitive display recipe; value recomputed from current prod_v2 neighbours",
            }
        )
    return out


def build_review_adjustment_map(
    loaded: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    """Build the explicit, author-requested local smoothing candidates.

    Each target is checked against its named neighbours and the value is
    recomputed from the current raw rows.  Keeping this separate from the
    inherited recipe makes the provenance and later author acceptance
    decision unambiguous.
    """
    out: list[dict[str, Any]] = []
    for spec in REVIEW_ADJUSTMENTS:
        mode_key = spec["mode_key"]
        target = current_row(loaded, mode_key, spec["plot_panel"], spec["plot_series"], spec["xi"])
        left = current_row(loaded, mode_key, spec["plot_panel"], spec["plot_series"], spec["left_xi"])
        right = current_row(loaded, mode_key, spec["plot_panel"], spec["plot_series"], spec["right_xi"])
        if target is None or left is None or right is None:
            raise ValueError(
                f"review adjustment key missing from current {mode_key}: "
                f"{spec['plot_panel']} {spec['plot_series']} {spec['xi']}"
            )
        target_xi = parse_finite(target, "xi")
        left_xi = parse_finite(left, "xi")
        right_xi = parse_finite(right, "xi")
        if not left_xi < target_xi < right_xi:
            raise ValueError(f"invalid review adjustment bracket for {spec['window_id']}")
        weight = (target_xi - left_xi) / (right_xi - left_xi)
        for observable in spec["observables"]:
            raw_value = parse_finite(target, observable)
            left_value = parse_finite(left, observable)
            right_value = parse_finite(right, observable)
            derived_value = left_value + weight * (right_value - left_value)
            neighbour_mean = (left_value + right_value) / 2.0
            residual = raw_value - derived_value
            residual_relative = residual / neighbour_mean if neighbour_mean else math.nan
            if not all(math.isfinite(value) for value in (derived_value, residual, residual_relative)):
                raise ValueError(f"non-finite review adjustment for {spec['window_id']} {observable}")
            out.append(
                {
                    "window_id": spec["window_id"],
                    "mode_key": mode_key,
                    "plot_panel": spec["plot_panel"],
                    "plot_series": spec["plot_series"],
                    "observable": observable,
                    "xi": canonical_xi(target["xi"]),
                    "raw_production_value_current": raw_value,
                    "recipe_raw_production_value": "",
                    "recipe_display_value": "",
                    "left_xi": canonical_xi(left["xi"]),
                    "left_value_current": left_value,
                    "right_xi": canonical_xi(right["xi"]),
                    "right_value_current": right_value,
                    "derived_display_value": derived_value,
                    "local_residual": residual,
                    "local_residual_relative": residual_relative,
                    "replacement_method": "linear interpolation between current prod_v2 neighbours",
                    "recipe_source": "author_review_request",
                    "recipe_source_sha256": "",
                    "canonical_data_modified": False,
                    "display_status": "author_requested_interpolation",
                    "adjustment_type": "author_requested_smoothing",
                    "adjustment_reason": spec["reason"],
                }
            )
    return out


def build_marker_map(
    loaded: dict[str, dict[str, Any]], recipes: list[dict[str, str]]
) -> list[dict[str, Any]]:
    """Reconcile historical marker keys with current raw rows.

    Historical stars are evidence only.  The publication layer renders the
    explicit interval midpoints returned by :func:`build_publication_marker_map`
    instead, so neither the stale T=120 marker nor the direct-coexistence side
    rows are rendered as raw stars.
    """
    out: list[dict[str, Any]] = []
    for recipe in recipes:
        mode_key = recipe["mode_key"]
        direct_coexistence = (
            mode_key == "mode_a"
            and recipe["plot_panel"] == "muB900.0"
            and recipe["plot_series"] == "alpha1.0"
            and math.isclose(float(recipe["xi"]), 0.0, abs_tol=1.0e-12)
        )
        target = current_row(loaded, mode_key, recipe["plot_panel"], recipe["plot_series"], recipe["xi"])
        if target is not None:
            historical_non_cep = mode_key == "mode_b" and recipe["plot_panel"] == "T120.0"
            suppress_historical = historical_non_cep or direct_coexistence
            out.append(
                {
                    "window_id": recipe["window_id"],
                    "mode_key": mode_key,
                    "plot_panel": recipe["plot_panel"],
                    "plot_series": recipe["plot_series"],
                    "recipe_xi": recipe["xi"],
                    "render_xi": canonical_xi(target["xi"]),
                    "observable": recipe["observable"],
                    "raw_production_value_current": parse_finite(target, recipe["observable"]),
                    "recipe_raw_production_value": recipe["raw_production_value"],
                    "marker": recipe["marker"],
                    "marker_status": (
                        "suppressed_historical_non_CEP"
                        if historical_non_cep
                        else "suppressed_historical_direct_coexistence_midpoint"
                        if direct_coexistence
                        else "applied_current_raw_point"
                    ),
                    "marker_semantics": (
                        "historical first-order branch marker retained for audit only; not a CEP"
                        if historical_non_cep
                        else "historical direct-coexistence xi=0 recipe replaced by an interval-midpoint display marker"
                        if direct_coexistence
                        else "first-order/upstream branch transition point retained without smoothing"
                    ),
                    "render_marker": not suppress_historical,
                    "canonical_data_modified": False,
                }
            )
            continue
        if mode_key == "mode_a":
            minus = current_row(loaded, mode_key, recipe["plot_panel"], recipe["plot_series"], "-0.003")
            plus = current_row(loaded, mode_key, recipe["plot_panel"], recipe["plot_series"], "0.003")
            if minus is None or plus is None:
                raise ValueError("mode_a direct-coexistence side-point contract is incomplete")
            for side, side_row in (("minus", minus), ("plus", plus)):
                out.append(
                    {
                        "window_id": recipe["window_id"],
                        "mode_key": mode_key,
                        "plot_panel": recipe["plot_panel"],
                        "plot_series": recipe["plot_series"],
                        "recipe_xi": recipe["xi"],
                        "render_xi": canonical_xi(side_row["xi"]),
                        "observable": recipe["observable"],
                        "raw_production_value_current": parse_finite(side_row, recipe["observable"]),
                        "recipe_raw_production_value": recipe["raw_production_value"],
                        "marker": recipe["marker"],
                        "marker_status": "suppressed_historical_direct_coexistence_side_point",
                        "marker_semantics": "historical direct-coexistence side point retained for interval provenance; publication renders the midpoint only",
                        "render_marker": False,
                        "coexistence_side": side,
                        "canonical_data_modified": False,
                    }
                )
            continue
        raise ValueError(f"marker key missing from current {mode_key}: {recipe}")
    return out


def build_cep_slice_audit(
    loaded: dict[str, dict[str, Any]],
    markers: list[dict[str, Any]],
    cep_rows: list[dict[str, str]],
) -> list[dict[str, Any]]:
    """Prove whether a historical marker's fixed plot slice intersects CEP.

    CEP input is parameterized by ``xi`` and ``mu_q``.  The transport plots
    use fixed ``T``/``mu_B`` slices, so a marker is eligible only when one
    strict CEP row contains the plotted temperature bracket and its converted
    ``mu_B = 3 * mu_q`` agrees with the plotted series.  We intentionally do
    not place a marker at the nearest point when the two coordinates match at
    different xi values.
    """
    out: list[dict[str, Any]] = []
    seen: set[tuple[str, str, str]] = set()
    for marker in markers:
        if not (
            marker["mode_key"] == "mode_b"
            and marker["plot_panel"] == "T120.0"
        ):
            continue
        curve_key = (
            marker["mode_key"],
            marker["plot_panel"],
            marker["plot_series"],
        )
        if curve_key in seen:
            continue
        seen.add(curve_key)
        payload = loaded[marker["mode_key"]]
        curve_rows = payload["curves"].get((marker["plot_panel"], marker["plot_series"]), [])
        if not curve_rows:
            raise ValueError(f"CEP audit curve missing for {curve_key}")
        T_slice = parse_finite(curve_rows[0], "T_MeV")
        muB_slice = parse_finite(curve_rows[0], "muB_MeV")
        parsed = []
        for row in cep_rows:
            T_low = parse_finite(row, "T_low_MeV")
            T_high = parse_finite(row, "T_high_MeV")
            T_mid = parse_finite(row, "T_midpoint_MeV")
            muB_cep = 3.0 * parse_finite(row, "mu_CEP_proxy_MeV")
            parsed.append(
                {
                    "xi": canonical_xi(row["xi"]),
                    "T_low": T_low,
                    "T_high": T_high,
                    "T_mid": T_mid,
                    "muB": muB_cep,
                    "T_delta": abs(T_mid - T_slice),
                    "muB_delta": abs(muB_cep - muB_slice),
                }
            )
        nearest_T = min(parsed, key=lambda row: (row["T_delta"], row["muB_delta"], row["xi"]))
        nearest_muB = min(parsed, key=lambda row: (row["muB_delta"], row["T_delta"], row["xi"]))
        strict = [
            row
            for row in parsed
            if row["T_low"] <= T_slice <= row["T_high"]
            and math.isclose(row["muB"], muB_slice, rel_tol=0.0, abs_tol=1.0e-8)
        ]
        out.append(
            {
                "mode_key": marker["mode_key"],
                "plot_panel": marker["plot_panel"],
                "plot_series": marker["plot_series"],
                "T_slice_MeV": T_slice,
                "muB_slice_MeV": muB_slice,
                "cep_table_rows": len(parsed),
                "strict_intersection_count": len(strict),
                "strict_intersection": bool(strict),
                "nearest_T_xi": nearest_T["xi"],
                "nearest_T_midpoint_MeV": nearest_T["T_mid"],
                "nearest_T_muB_MeV": nearest_T["muB"],
                "nearest_T_delta_MeV": nearest_T["T_delta"],
                "nearest_muB_xi": nearest_muB["xi"],
                "nearest_muB_T_midpoint_MeV": nearest_muB["T_mid"],
                "nearest_muB_MeV": nearest_muB["muB"],
                "nearest_muB_delta_MeV": nearest_muB["muB_delta"],
                "marker_action": (
                    "render_strict_CEP_marker"
                    if strict
                    else "suppress_historical_non_CEP_marker"
                ),
                "reason": (
                    "strict CEP intersection exists on this fixed slice"
                    if strict
                    else "T and muB closest matches occur at different xi; nearest point is not a CEP"
                ),
                "reference_path": relpath(CEP_BOUNDARY),
                "reference_sha256": sha256_file(CEP_BOUNDARY),
                "canonical_data_modified": False,
            }
        )
    return out


def interpolate_curve_value(
    curve_rows: list[dict[str, str]], xi: float, field: str
) -> float:
    ordered = sorted(curve_rows, key=lambda row: float(row["xi"]))
    if xi < float(ordered[0]["xi"]) or xi > float(ordered[-1]["xi"]):
        raise ValueError(f"CEP xi={xi} lies outside transport curve support")
    for row in ordered:
        if math.isclose(float(row["xi"]), xi, rel_tol=0.0, abs_tol=1.0e-10):
            return parse_finite(row, field)
    for left, right in zip(ordered, ordered[1:]):
        x_left = float(left["xi"])
        x_right = float(right["xi"])
        if x_left < xi < x_right:
            y_left = parse_finite(left, field)
            y_right = parse_finite(right, field)
            weight = (xi - x_left) / (x_right - x_left)
            return y_left + weight * (y_right - y_left)
    raise ValueError(f"unable to interpolate CEP xi={xi}")


def build_publication_marker_map(
    loaded: dict[str, dict[str, Any]],
    specs: Iterable[dict[str, str]] | None = None,
) -> list[dict[str, Any]]:
    """Build display-only midpoint markers from confirmed phase-switch brackets.

    A strict CEP search is a bisection bracket, not an infinitely precise point:
    its two endpoints carry different phase labels.  The publication convention
    therefore places one star at the arithmetic midpoint of the confirmed
    minimum bracket.  The endpoint rows and their values/labels are copied into
    the audit table, while the synthetic midpoint is kept out of
    ``publication_clean_points.csv`` so no transport value is fabricated as raw
    data.
    """
    out: list[dict[str, Any]] = []
    for spec in PHASE_SWITCH_MARKERS if specs is None else specs:
        mode_key = spec["mode_key"]
        left = current_row(
            loaded, mode_key, spec["plot_panel"], spec["plot_series"], spec["left_xi"]
        )
        right = current_row(
            loaded, mode_key, spec["plot_panel"], spec["plot_series"], spec["right_xi"]
        )
        if left is None or right is None:
            raise ValueError(f"publication marker bracket missing from {mode_key}: {spec}")
        left_xi = parse_finite(left, "xi")
        right_xi = parse_finite(right, "xi")
        if not left_xi < right_xi:
            raise ValueError(f"publication marker bracket is not ordered: {spec['window_id']}")
        midpoint_xi = (left_xi + right_xi) / 2.0
        if not math.isfinite(midpoint_xi):
            raise ValueError(f"non-finite publication marker midpoint: {spec['window_id']}")
        for field in ("T_MeV", "muB_MeV"):
            if not math.isclose(
                parse_finite(left, field), parse_finite(right, field), rel_tol=0.0, abs_tol=1.0e-9
            ):
                raise ValueError(f"publication marker endpoints disagree in {field}: {spec['window_id']}")
        if spec["marker_basis"] == "phase_label_switch":
            left_kind = left["phase_reference_kind"]
            right_kind = right["phase_reference_kind"]
            if left_kind == right_kind:
                raise ValueError(f"phase-label bracket has no label switch: {spec['window_id']}")
        for observable in DISPLAY_FIELDS:
            left_value = parse_finite(left, observable)
            right_value = parse_finite(right, observable)
            display_value = (left_value + right_value) / 2.0
            if not math.isfinite(display_value):
                raise ValueError(f"non-finite publication marker value: {spec['window_id']} {observable}")
            out.append(
                {
                    "window_id": spec["window_id"],
                    "mode_key": mode_key,
                    "plot_panel": spec["plot_panel"],
                    "plot_series": spec["plot_series"],
                    "observable": observable,
                    "T_MeV": parse_finite(left, "T_MeV"),
                    "muB_MeV": parse_finite(left, "muB_MeV"),
                    "marker": "star",
                    "marker_status": "confirmed_interval_midpoint",
                    "marker_basis": spec["marker_basis"],
                    "marker_semantics": (
                        "confirmed minimum CEP/phase-switch interval midpoint; display-only proxy, "
                        "not an exact zero-width CEP coordinate"
                    ),
                    "interval_xi_low": canonical_xi(left_xi),
                    "interval_xi_high": canonical_xi(right_xi),
                    "interval_width_xi": right_xi - left_xi,
                    "render_xi": canonical_xi(midpoint_xi),
                    "left_phase_reference_kind": left["phase_reference_kind"],
                    "right_phase_reference_kind": right["phase_reference_kind"],
                    "left_phase_structure": left["phase_structure"],
                    "right_phase_structure": right["phase_structure"],
                    "left_quality_flag": left["quality_flag"],
                    "right_quality_flag": right["quality_flag"],
                    "left_value": left_value,
                    "right_value": right_value,
                    "display_value": display_value,
                    "reason": spec["reason"],
                    "render_marker": True,
                    "canonical_data_modified": False,
                }
            )
    return out


def build_cep_marker_map(
    loaded: dict[str, dict[str, Any]],
    cep_slice_audit: list[dict[str, Any]],
    cep_rows: list[dict[str, str]],
) -> list[dict[str, Any]]:
    """Build display-only stars for genuine fixed-slice CEP intersections."""
    out: list[dict[str, Any]] = []
    for audit in cep_slice_audit:
        if not audit["strict_intersection"]:
            continue
        payload = loaded[audit["mode_key"]]
        curve_rows = payload["curves"].get(
            (audit["plot_panel"], audit["plot_series"]), []
        )
        if not curve_rows:
            raise ValueError(f"CEP marker curve missing for {audit}")
        T_slice = float(audit["T_slice_MeV"])
        muB_slice = float(audit["muB_slice_MeV"])
        matches = []
        for row in cep_rows:
            T_low = parse_finite(row, "T_low_MeV")
            T_high = parse_finite(row, "T_high_MeV")
            muB_cep = 3.0 * parse_finite(row, "mu_CEP_proxy_MeV")
            if T_low <= T_slice <= T_high and math.isclose(
                muB_cep, muB_slice, rel_tol=0.0, abs_tol=1.0e-8
            ):
                matches.append(row)
        for row in matches:
            xi = float(row["xi"])
            for observable in DISPLAY_FIELDS:
                out.append(
                    {
                        "window_id": (
                            f"strict_cep_{audit['mode_key']}_{audit['plot_panel']}"
                            f"_{audit['plot_series']}_{canonical_xi(xi)}"
                        ),
                        "mode_key": audit["mode_key"],
                        "plot_panel": audit["plot_panel"],
                        "plot_series": audit["plot_series"],
                        "recipe_xi": canonical_xi(xi),
                        "render_xi": canonical_xi(xi),
                        "observable": observable,
                        "raw_production_value_current": interpolate_curve_value(
                            curve_rows, xi, observable
                        ),
                        "recipe_raw_production_value": "",
                        "marker": "star",
                        "marker_status": "strict_CEP_marker",
                        "marker_semantics": (
                            "strict CEP boundary intersection; y value linearly interpolated "
                            "from current raw transport neighbours"
                        ),
                        "render_marker": True,
                        "coexistence_side": "",
                        "canonical_data_modified": False,
                    }
                )
    return out


def build_marker_semantics_audit(
    loaded: dict[str, dict[str, Any]],
    markers: list[dict[str, Any]],
    cep_slice_audit: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Record whether each historical marker is eligible for publication."""
    out: list[dict[str, Any]] = []
    cep_by_curve = {
        (row["mode_key"], row["plot_panel"], row["plot_series"]): row
        for row in cep_slice_audit
    }
    for marker in markers:
        row = current_row(
            loaded,
            marker["mode_key"],
            marker["plot_panel"],
            marker["plot_series"],
            marker["render_xi"],
        )
        if row is None:
            raise ValueError(f"marker audit row missing for {marker}")
        direct_coexistence = (
            marker["mode_key"] == "mode_a"
            and marker["plot_panel"] == "muB900.0"
            and marker["plot_series"] == "alpha1.0"
        )
        if marker["mode_key"] == "mode_b" and marker["plot_panel"] == "T120.0":
            slice_audit = cep_by_curve[
                (marker["mode_key"], marker["plot_panel"], marker["plot_series"])
            ]
            intended = "historical_first_order_branch_marker"
            cep_semantics = "confirmed_interval_midpoint_proxy"
            verdict = "suppressed_historical_marker; publication renders the author-confirmed phase-label-switch midpoint"
            evidence = (
                "paper_first_order_markers.csv; current phase_reference_kind/phase_structure; "
                "first_order_protection.csv; tables/cep_marker_audit.csv; tables/publication_marker_map.csv"
            )
        elif direct_coexistence:
            intended = "direct_coexistence_phase_switch"
            cep_semantics = "confirmed_interval_midpoint_proxy"
            verdict = "suppressed_historical_side_marker; publication renders the xi=0 direct-coexistence interval midpoint"
            evidence = "paper_first_order_markers.csv; direct-coexistence marker contract; tables/publication_marker_map.csv"
        else:
            intended = "first_order_boundary_or_transition"
            cep_semantics = "not_a_CEP_marker"
            verdict = "consistent_non_CEP_first_order_marker"
            evidence = "paper_first_order_markers.csv; direct-coexistence marker contract"
        out.append(
            {
                "window_id": marker["window_id"],
                "mode_key": marker["mode_key"],
                "plot_panel": marker["plot_panel"],
                "plot_series": marker["plot_series"],
                "observable": marker["observable"],
                "render_xi": marker["render_xi"],
                "marker_status": marker["marker_status"],
                "intended_semantics": intended,
                "cep_semantics": cep_semantics,
                "phase_reference_kind": row["phase_reference_kind"],
                "phase_structure": row["phase_structure"],
                "quality_flag": row["quality_flag"],
                "render_marker": marker.get("render_marker", True),
                "audit_verdict": verdict,
                "evidence": evidence,
                "canonical_data_modified": False,
            }
        )
    return out


def build_clean_points(
    loaded: dict[str, dict[str, Any]],
    replacements: list[dict[str, Any]],
    markers: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    replacement_index = {
        (row["mode_key"], row["plot_panel"], row["plot_series"], row["xi"], row["observable"]): row
        for row in replacements
    }
    marker_index = {
        (row["mode_key"], row["plot_panel"], row["plot_series"], row["render_xi"], row["observable"]): row
        for row in markers
        if row.get("render_marker", True)
    }
    out: list[dict[str, Any]] = []
    for mode_key, payload in loaded.items():
        for row in payload["rows"]:
            for observable in DISPLAY_FIELDS:
                xi = canonical_xi(row["xi"])
                replacement = replacement_index.get(
                    (mode_key, row["plot_panel"], row["plot_series"], xi, observable)
                )
                marker = marker_index.get(
                    (mode_key, row["plot_panel"], row["plot_series"], xi, observable)
                )
                if replacement is not None:
                    clean_value = replacement["derived_display_value"]
                    status = replacement.get("display_status", "interpolated_noncertified")
                    source = (
                        "current prod_v2 neighbours (author-requested local smoothing)"
                        if status == "author_requested_interpolation"
                        else "current prod_v2 left/right raw neighbours"
                    )
                elif marker is not None:
                    clean_value = marker["raw_production_value_current"]
                    status = "first_order_raw_marker"
                    source = "current prod_v2 raw point"
                else:
                    clean_value = parse_finite(row, observable)
                    status = "raw"
                    source = "current prod_v2 raw point"
                out.append(
                    {
                        "mode_key": mode_key,
                        "mode": row["mode"],
                        "plot_panel": row["plot_panel"],
                        "plot_series": row["plot_series"],
                        "plot_series_label": row["plot_series_label"],
                        "T_MeV": row["T_MeV"],
                        "muB_MeV": row["muB_MeV"],
                        "xi": xi,
                        "observable": observable,
                        "raw_value": parse_finite(row, observable),
                        "clean_value": clean_value,
                        "display_status": status,
                        "value_source": source,
                        "phase_structure": row["phase_structure"],
                        "phase_reference_kind": row["phase_reference_kind"],
                        "quality_flag": row["quality_flag"],
                        "quality_reason": row["quality_reason"],
                        "run_id": row.get("run_id", ""),
                        "canonical_data_modified": False,
                    }
                )
    return out


def build_curve_index(points: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in points:
        grouped[(row["mode_key"], row["plot_panel"], row["plot_series"], row["observable"])].append(row)
    out = []
    for key, rows in sorted(grouped.items()):
        out.append(
            {
                "mode_key": key[0],
                "plot_panel": key[1],
                "plot_series": key[2],
                "observable": key[3],
                "point_count": len(rows),
                "replacement_count": sum(
                    row["display_status"]
                    in {"interpolated_noncertified", "author_requested_interpolation"}
                    for row in rows
                ),
                "marker_count": sum(row["display_status"] == "first_order_raw_marker" for row in rows),
                "xi_min": min(float(row["xi"]) for row in rows),
                "xi_max": max(float(row["xi"]) for row in rows),
                "canonical_data_modified": False,
            }
        )
    return out


def render_figures(
    points: list[dict[str, Any]],
    cep_markers: list[dict[str, Any]] | None = None,
    publication_markers: list[dict[str, Any]] | None = None,
) -> list[Path]:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - exercised only in plot environments
        raise RuntimeError("matplotlib is required to render publication-clean figures") from exc

    colors = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE"]
    labels = {"eta_over_s": r"$\eta/s$", "zeta_over_s": r"$\zeta/s$", "sigma_over_T": r"$\sigma/T$"}
    grouped: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in points:
        grouped[(row["mode_key"], row["plot_panel"], row["plot_series"])].append(row)
    cep_markers = cep_markers or []
    publication_markers = publication_markers or []
    paths: list[Path] = []
    for mode_key in MODE_CONFIG:
        panels = sorted({key[1] for key in grouped if key[0] == mode_key})
        for panel in panels:
            series_names = sorted({key[2] for key in grouped if key[:2] == (mode_key, panel)})
            for observable in DISPLAY_FIELDS:
                fig, ax = plt.subplots(figsize=(6.75, 4.6))
                marker_present = False
                publication_marker_present = False
                for series_index, series in enumerate(series_names):
                    rows = [
                        row
                        for row in points
                        if row["mode_key"] == mode_key
                        and row["plot_panel"] == panel
                        and row["plot_series"] == series
                        and row["observable"] == observable
                    ]
                    rows.sort(key=lambda row: float(row["xi"]))
                    ax.plot(
                        [float(row["xi"]) for row in rows],
                        [float(row["clean_value"]) for row in rows],
                        color=colors[series_index % len(colors)],
                        linewidth=1.5,
                        label=rows[0]["plot_series_label"],
                    )
                    marker_rows = [row for row in rows if row["display_status"] == "first_order_raw_marker"]
                    if marker_rows:
                        marker_present = True
                        ax.scatter(
                            [float(row["xi"]) for row in marker_rows],
                            [float(row["raw_value"]) for row in marker_rows],
                            marker="*",
                            s=95,
                            facecolor=colors[series_index % len(colors)],
                            edgecolor="black",
                            linewidth=0.55,
                            zorder=5,
                        )
                    cep_rows = [
                        row
                        for row in cep_markers
                        if row["mode_key"] == mode_key
                        and row["plot_panel"] == panel
                        and row["plot_series"] == series
                        and row["observable"] == observable
                    ]
                    if cep_rows:
                        ax.scatter(
                            [float(row["render_xi"]) for row in cep_rows],
                            [float(row["raw_production_value_current"]) for row in cep_rows],
                            marker="*",
                            s=120,
                            facecolor="#CCBB44",
                            edgecolor="black",
                            linewidth=0.65,
                            zorder=6,
                        )
                        ax.scatter(
                            [],
                            [],
                            marker="*",
                            s=90,
                            facecolor="#CCBB44",
                            edgecolor="black",
                            linewidth=0.65,
                            label="CEP (strict intersection)",
                        )
                    publication_rows = [
                        row
                        for row in publication_markers
                        if row["mode_key"] == mode_key
                        and row["plot_panel"] == panel
                        and row["plot_series"] == series
                        and row["observable"] == observable
                    ]
                    if publication_rows:
                        publication_marker_present = True
                        ax.scatter(
                            [float(row["render_xi"]) for row in publication_rows],
                            [float(row["display_value"]) for row in publication_rows],
                            marker="*",
                            s=120,
                            facecolor="#CCBB44",
                            edgecolor="black",
                            linewidth=0.65,
                            zorder=7,
                        )
                if publication_marker_present:
                    ax.scatter(
                        [],
                        [],
                        marker="*",
                        s=90,
                        facecolor="#CCBB44",
                        edgecolor="black",
                        linewidth=0.65,
                        label="CEP / phase-switch interval midpoint",
                    )
                if marker_present:
                    ax.scatter([], [], marker="*", s=75, facecolor="white", edgecolor="black", linewidth=0.55, label="first-order transition")
                ax.set_xlabel(r"$\xi$")
                ax.set_ylabel(labels[observable])
                ax.set_xlim(-0.52, 0.52)
                ax.legend(loc="best")
                fig.tight_layout()
                path = FIGURE_DIR / mode_key / f"plot_panel={panel}" / f"{observable}_vs_xi.png"
                path.parent.mkdir(parents=True, exist_ok=True)
                fig.savefig(path, dpi=600, bbox_inches="tight", pad_inches=0.08)
                plt.close(fig)
                paths.append(path)
    return paths


def claim_ledger(
    inherited_replacements: list[dict[str, Any]],
    review_adjustments: list[dict[str, Any]],
    markers: list[dict[str, Any]],
    publication_markers: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    midpoint_count = sum(
        row["marker_status"] == "confirmed_interval_midpoint"
        for row in publication_markers
    )
    return [
        {
            "claim_id": "PC-V1-001",
            "status": "supported",
            "claim_zh": "publication_clean_v1 只消费作者接受的两棵 prod_v2 raw result tree；raw CSV、registry 和正式 figure 均未被修改，派生过程未调用 solver。",
            "evidence": "tables/input_inventory.csv; manifest.json; figures/plot_manifest.json",
            "scope_limit": "raw numerical audit 仍为 diagnostic-only，publication_clean_v1 不改变收敛或物理结论。",
        },
        {
            "claim_id": "PC-V1-002",
            "status": "supported_with_scope_limit",
            "claim_zh": f"{len(inherited_replacements)} 个既有作者审核过的非一阶显示异常点按当前 prod_v2 左右真实邻点重新线性插值；结果只写入派生 clean_value，不覆盖 raw_value。",
            "evidence": "tables/replacement_map.csv; tables/publication_clean_points.csv",
            "scope_limit": "插值值不是新的 solver 结果，也不是有限宽度/极点正则化；只允许用于明确标注的显示层。",
        },
        {
            "claim_id": "PC-V1-003",
            "status": "supported_with_scope_limit",
            "claim_zh": f"一阶窗口保留当前 raw 值；历史 mode-A ξ=0 配方和 T=120、μB=900 历史星标均只作审计输入，论文图改为渲染 {midpoint_count} 个 observable-level 的已确认有限区间中点 marker，不伪造 raw 输运值。",
            "evidence": "tables/first_order_marker_map.csv; tables/publication_marker_map.csv; tables/publication_clean_points.csv; data/outputs/results/relaxtime/transport/phase_guided/*/manifest.json",
            "scope_limit": "星标是 CEP/phase-switch interval midpoint 的 display-only proxy；左右端点和 phase labels 必须在 publication_marker_map.csv 中保留，不能宣称零宽度精确 CEP 或新 solver 结果。",
        },
        {
            "claim_id": "PC-V1-004",
            "status": "author_check",
            "claim_zh": "派生图仍不是 manuscript eligible；mode-A μB=900、αT=1.0、ξ=−0.01 的 bulk/equilibrium 分支不一致继续作为已知排除项。",
            "evidence": "tables/bulk_derivative_branch_audit.csv; figures/plot_manifest.json",
            "scope_limit": "在稳定分支与 bulk base_state 复用修复并重跑前，不得将该曲线声明为论文定稿输入。",
        },
        {
            "claim_id": "PC-V1-005",
            "status": "not_claimed",
            "claim_zh": "本层不声称 RS numerical parity、全域 solver 收敛或旧 reference 已删除。",
            "evidence": "README.md; manifest.json",
            "scope_limit": "old-reference retirement 另立审计/PR，并在作者授权前不物理删除旧 prod_v1。",
        },
        {
            "claim_id": "PC-V1-006",
            "status": "author_check",
            "claim_zh": f"本轮作者审阅新增 {len(review_adjustments)} 个 T=200 MeV 局部显示平滑候选：mode B 的 (μB=900, ξ=−0.10) 与 (μB=0, ξ=0.36) 各覆盖 η/s、ζ/s；均由当前 raw 左右邻点线性重算。",
            "evidence": "tables/review_adjustment_map.csv; tables/publication_clean_points.csv; figures/mode_b/plot_panel=T200.0",
            "scope_limit": "这是待作者确认的显示层修正，不改变 raw result、phase 标签或收敛结论；候选数值不是新的 solver 结果。",
        },
        {
            "claim_id": "PC-V1-007",
            "status": "supported_with_scope_limit",
            "claim_zh": "历史配方中 T=120 MeV、μB=900 MeV、ξ=−0.09 的星标已从 publication-clean 图层抑制；图上改用当前 prod_v2 中作者确认的 phase-reference_kind crossover → first_order 最小切换区间中点 ξ=−0.135。",
            "evidence": "tables/marker_semantics_audit.csv; tables/first_order_marker_map.csv; tables/publication_marker_map.csv; tables/cep_marker_audit.csv; data/reference/pnjl/issue130_phase_reference_v1/*",
            "scope_limit": "strict CEP 表仍是独立坐标审计；ξ=−0.135 星标表示有限宽度 CEP/phase-switch bracket 的中点 proxy，不是旧 ξ=−0.09 分支标记，也不是零宽度精确交点。",
        },
        {
            "claim_id": "PC-V1-008",
            "status": "supported_with_scope_limit",
            "claim_zh": "mode-A μB=450 MeV、αT=1.0、ξ=−0.20 的连续区快速斜率变化与既有 simple_1m4KΠ 小分母机制窗口一致；该点 phase_reference_kind=crossover，非一阶相变。",
            "evidence": "docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_p128_xi001_analysis/tables/mechanism_window_summary.csv; tables/denominator_chain_summary.csv; tables/downstream_transport_response_summary.csv",
            "scope_limit": "机制证据是定点分解/复现支持，不把它升级为新的 solver 收敛或论文定稿结论。",
        },
    ]


def render_readme(
    inventory: list[dict[str, Any]],
    inherited_replacements: list[dict[str, Any]],
    review_adjustments: list[dict[str, Any]],
    markers: list[dict[str, Any]],
    rendered_markers: list[dict[str, Any]],
    cep_slice_audit: list[dict[str, Any]],
    cep_markers: list[dict[str, Any]],
    publication_markers: list[dict[str, Any]],
    figure_paths: list[Path],
) -> str:
    inventory_lines = "\n".join(
        f"| {row['mode_key']} | {row['scan_rows']} | {row['diagnostic_rows']} | `{row['scan_sha256']}` | `{row['diagnostics_sha256']}` |"
        for row in inventory
    )
    return f"""# Issue #130 RS `publication_clean_v1` 派生层

## 目的与边界

本包从作者接受的两套 `prod_v2` 原始 RS transport 结果生成论文显示候选层。它只做 solver-free 后处理：不修改 `data/outputs/results/**`、`production_registry.json` 或正式 figure，不调用 equilibrium/transport solver，也不把派生值写回 raw CSV。

`publication_clean` 的含义是：沿用已审核的极点/小分母显示配方，并在当前 `prod_v2` 网格上重新取左右真实邻点计算显示插值；一阶相变点继续保留 raw 值。本轮另记录作者提出的局部平滑候选，但它们仍待本轮作者确认。该层默认 `manuscript_eligible=false`，须经作者审核后才能用于论文。

## 输入 provenance

| mode | scan rows | diagnostic rows | scan SHA256 | diagnostic SHA256 |
| --- | ---: | ---: | --- | --- |
{inventory_lines}

- calculation SHA：`{CALCULATION_SHA}`
- workflow head：`{WORKFLOW_HEAD_SHA}`
- source result case：`{CURRENT_CASE}`
- source solver：已调用；本次派生：`solver_called=false`
- marker contract：`confirmed_interval_midpoint_v1`（严格 CEP/相标签切换区间的中点 display proxy）
- 旧显示配方：`{relpath(REPLACEMENT_RECIPE)}`（SHA 记录于 manifest）
- 本包生成图：{len(figure_paths)} 张 PNG，见 `figures/plot_manifest.json`
- 历史 marker：{len(markers)} 条审计记录，其中 {len(rendered_markers)} 条历史 raw marker 实际渲染；T=120、μB=900 的 3 条历史星标被抑制。
- CEP-slice 审计：{len(cep_slice_audit)} 个固定切片，严格固定切片交点星标为 {len(cep_markers)} 条。
- 论文显示 marker：{len(publication_markers)} 条（3 个 observable × 2 个已确认区间），均为区间中点的 display-only proxy。

## 派生规则

1. 继承配方中的 19 个非一阶显示替换点从当前 raw 的 `left_xi < target_xi < right_xi` 三点线性插值；`raw_value` 与左右端点同时保留。
2. 旧配方中 `mode_a, μB=900, αT=1.0, ξ=0` 当前不存在：direct-coexistence 合同保留 `ξ=−0.003/+0.003` 两侧 raw 端点，并在图上只显示其区间中点 `ξ=0` 的 display-only marker；不把中点写成输运 raw 值。
3. 一阶点不做平滑；所有派生行带 `canonical_data_modified=false`。`quality_flag`/`quality_reason` 原样带入。
4. 已知 `mode_a, μB=900, αT=1.0, ξ=−0.01` bulk 分支问题仍列为排除项；本包不会用插值掩盖它。
5. 本轮作者审阅新增 4 个 T=200 MeV 局部平滑候选：`(mode B, μB=900, ξ=−0.10)` 与 `(mode B, μB=0, ξ=0.36)` 的 `η/s`、`ζ/s`；全部按当前 raw 左右邻点线性重算，单列于 `review_adjustment_map.csv`，不改变原始结果。
6. 历史配方中 T=120 MeV、μB=900 MeV、ξ=−0.09 的星标只作为审计输入，不进入 publication-clean 图；图上改为显示作者确认的 `phase_reference_kind: crossover → first_order` 最小切换区间 `ξ∈[−0.14,−0.13]` 的中点 `ξ=−0.135`。严格 CEP 表仍作为独立坐标审计，区间中点是论文显示 proxy，不宣称零宽度精确 CEP 坐标。
7. `publication_marker_map.csv` 保存两个已确认区间的左右端点 phase label、质量标志和 observable 值；`ξ=0` 与 `ξ=−0.135` 的星标都不进入 `publication_clean_points.csv`。
8. mode-A μB=450 MeV、αT=1.0、ξ=−0.20 的非一阶斜率变化与既有 `simple_1m4KΠ` 小分母机制窗口一致；这属于机制归因证据，不是新的相变标签。

## 结果文件

- `tables/input_inventory.csv`：输入路径、行数、哈希及 solver provenance。
- `tables/replacement_map.csv`：19 个继承配方加 4 个本轮作者审阅候选的当前邻点插值映射。
- `tables/review_adjustment_map.csv`：4 个本轮作者请求的 T=200 局部平滑候选及局部残差。
- `tables/first_order_marker_map.csv`：旧标记与当前 raw/direct-coexistence 端点合同的逐项对齐，并记录是否渲染。
- `tables/marker_semantics_audit.csv`：历史 marker 语义和 publication-clean 渲染决定。
- `tables/cep_marker_audit.csv`：固定输运切片与 strict CEP 表的坐标配对；严格交点与论文区间中点是分开的语义。
- `tables/cep_marker_map.csv`：仅在存在 strict 固定切片交点时生成的 CEP 星标及其显示值；当前为空。
- `tables/publication_marker_map.csv`：两个作者确认的最小相标签切换区间中点 marker，包含左右端点 provenance 和 display 值。
- `tables/publication_clean_points.csv`：三种论文展示 observable 的长表，含 raw/clean/status。
- `tables/curve_index.csv`：18 条 panel/series/observable 曲线的覆盖和替换计数。
- `tables/claim_ledger.csv`：证据强度、范围限制和未声明事项。
- `figures/plot_manifest.json` 与 `figures/**.png`：同构的 18 张显示候选图，未声明为 manuscript-ready。

## 复现

```powershell
python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v1.py
python -m pytest tests/unit/python/test_phase_guided_publication_clean_v1.py
```

执行脚本只写入本目录；任何 raw、registry、旧 `prod_v1` 和 legacy fallback 的变更都属于越界操作。
"""


def main() -> None:
    loaded, inventory = load_inputs()
    replacement_recipe, marker_recipe = load_recipe()
    cep_rows = load_cep_boundary()
    inherited_replacements = build_replacement_map(loaded, replacement_recipe)
    review_adjustments = build_review_adjustment_map(loaded)
    replacements = [*inherited_replacements, *review_adjustments]
    markers = build_marker_map(loaded, marker_recipe)
    cep_slice_audit = build_cep_slice_audit(loaded, markers, cep_rows)
    cep_markers = build_cep_marker_map(loaded, cep_slice_audit, cep_rows)
    publication_markers = build_publication_marker_map(loaded)
    marker_semantics = build_marker_semantics_audit(loaded, markers, cep_slice_audit)
    rendered_markers = [row for row in markers if row.get("render_marker", True)]
    points = build_clean_points(loaded, replacements, rendered_markers)
    curves = build_curve_index(points)
    figure_paths = render_figures(points, cep_markers, publication_markers)

    input_fields = [
        "mode_key", "mode", "scan_rows", "diagnostic_rows", "failed_rows", "xi_count",
        "scan_sha256", "diagnostics_sha256", "failed_sha256", "manifest_sha256",
        "effective_config_sha256", "calculation_sha", "workflow_head_sha",
        "source_solver_called", "derived_solver_called",
    ]
    write_csv(TABLE_DIR / "input_inventory.csv", inventory, input_fields)
    replacement_fields = [
        "window_id", "mode_key", "plot_panel", "plot_series", "observable", "xi",
        "raw_production_value_current", "recipe_raw_production_value", "recipe_display_value",
        "left_xi", "left_value_current", "right_xi", "right_value_current", "derived_display_value",
        "local_residual", "local_residual_relative", "replacement_method", "recipe_source", "recipe_source_sha256",
        "canonical_data_modified", "display_status", "adjustment_type", "adjustment_reason",
    ]
    write_csv(TABLE_DIR / "replacement_map.csv", replacements, replacement_fields)
    write_csv(TABLE_DIR / "review_adjustment_map.csv", review_adjustments, replacement_fields)
    marker_fields = [
        "window_id", "mode_key", "plot_panel", "plot_series", "recipe_xi", "render_xi", "observable",
        "raw_production_value_current", "recipe_raw_production_value", "marker", "marker_status",
        "marker_semantics", "render_marker", "coexistence_side", "canonical_data_modified",
    ]
    write_csv(TABLE_DIR / "first_order_marker_map.csv", markers, marker_fields)
    marker_audit_fields = [
        "window_id", "mode_key", "plot_panel", "plot_series", "observable", "render_xi",
        "marker_status", "intended_semantics", "cep_semantics", "phase_reference_kind",
        "phase_structure", "quality_flag", "render_marker", "audit_verdict", "evidence", "canonical_data_modified",
    ]
    write_csv(TABLE_DIR / "marker_semantics_audit.csv", marker_semantics, marker_audit_fields)
    cep_audit_fields = [
        "mode_key", "plot_panel", "plot_series", "T_slice_MeV", "muB_slice_MeV",
        "cep_table_rows", "strict_intersection_count", "strict_intersection",
        "nearest_T_xi", "nearest_T_midpoint_MeV", "nearest_T_muB_MeV", "nearest_T_delta_MeV",
        "nearest_muB_xi", "nearest_muB_T_midpoint_MeV", "nearest_muB_MeV", "nearest_muB_delta_MeV",
        "marker_action", "reason", "reference_path", "reference_sha256", "canonical_data_modified",
    ]
    write_csv(TABLE_DIR / "cep_marker_audit.csv", cep_slice_audit, cep_audit_fields)
    write_csv(TABLE_DIR / "cep_marker_map.csv", cep_markers, marker_fields)
    publication_marker_fields = [
        "window_id", "mode_key", "plot_panel", "plot_series", "observable", "T_MeV", "muB_MeV", "marker",
        "marker_status", "marker_basis", "marker_semantics", "interval_xi_low",
        "interval_xi_high", "interval_width_xi", "render_xi", "left_phase_reference_kind",
        "right_phase_reference_kind", "left_phase_structure", "right_phase_structure",
        "left_quality_flag", "right_quality_flag", "left_value", "right_value",
        "display_value", "reason", "render_marker", "canonical_data_modified",
    ]
    write_csv(TABLE_DIR / "publication_marker_map.csv", publication_markers, publication_marker_fields)
    point_fields = [
        "mode_key", "mode", "plot_panel", "plot_series", "plot_series_label", "T_MeV", "muB_MeV", "xi",
        "observable", "raw_value", "clean_value", "display_status", "value_source", "phase_structure",
        "phase_reference_kind", "quality_flag", "quality_reason", "run_id", "canonical_data_modified",
    ]
    write_csv(TABLE_DIR / "publication_clean_points.csv", points, point_fields)
    curve_fields = [
        "mode_key", "plot_panel", "plot_series", "observable", "point_count", "replacement_count",
        "marker_count", "xi_min", "xi_max", "canonical_data_modified",
    ]
    write_csv(TABLE_DIR / "curve_index.csv", curves, curve_fields)
    write_csv(
        TABLE_DIR / "claim_ledger.csv",
        claim_ledger(inherited_replacements, review_adjustments, markers, publication_markers),
        ["claim_id", "status", "claim_zh", "evidence", "scope_limit"],
    )

    generator_path = Path(__file__).resolve()
    plot_manifest = {
        "schema": "phase_guided_transport_publication_clean_plot_manifest_v1",
        "marker_contract": "confirmed_interval_midpoint_v1",
        "case": CURRENT_CASE,
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "base_git_commit": git_head(),
        "generator": relpath(generator_path),
        "generator_sha256": sha256_file(generator_path),
        "observables": DISPLAY_FIELDS,
        "replacement_count": len(replacements),
        "inherited_replacement_count": len(inherited_replacements),
        "review_adjustment_count": len(review_adjustments),
        "marker_recipe_count": len(marker_recipe),
        "marker_audit_count": len(markers),
        "marker_render_count": len(rendered_markers) + len(cep_markers) + len(publication_markers),
        "first_order_marker_render_count": len(rendered_markers),
        "historical_marker_suppressed_count": len(markers) - len(rendered_markers),
        "cep_marker_render_count": len(cep_markers),
        "publication_marker_render_count": len(publication_markers),
        "phase_switch_midpoint_count": len(publication_markers),
        "manuscript_eligible": False,
        "canonical_data_modified": False,
        "rendering_semantics": "current prod_v2 raw curves with inherited and author-requested adjacent-neighbour display replacements; historical markers suppressed from rendering; confirmed CEP/phase-switch interval midpoints rendered as display-only stars; exact strict fixed-slice intersections remain separately audited; no raw mutation",
        "figures": [
            {"path": relpath(path), "sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in figure_paths
        ],
    }
    write_json(FIGURE_DIR / "plot_manifest.json", plot_manifest)

    readme_path = OUT_DIR / "README.md"
    readme_path.write_text(
        render_readme(
            inventory,
            inherited_replacements,
            review_adjustments,
            markers,
            rendered_markers,
            cep_slice_audit,
            cep_markers,
            publication_markers,
            figure_paths,
        ),
        encoding="utf-8",
    )
    output_paths = [
        readme_path,
        *sorted(TABLE_DIR.glob("*.csv")),
        FIGURE_DIR / "plot_manifest.json",
        *figure_paths,
    ]
    manifest = {
        "schema": "phase_guided_transport_publication_clean_manifest_v1",
        "marker_contract": "confirmed_interval_midpoint_v1",
        "case": CURRENT_CASE,
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "base_git_commit": git_head(),
        "generator": relpath(generator_path),
        "generator_sha256": sha256_file(generator_path),
        "calculation_sha": CALCULATION_SHA,
        "workflow_head_sha": WORKFLOW_HEAD_SHA,
        "status": "derived_author_review_required",
        "manuscript_eligible": False,
        "canonical_data_modified": False,
        "production_write": False,
        "solver_called": False,
        "source_solver_called": True,
        "source_case": CURRENT_CASE,
        "source_registry_status": "approved_raw_manuscript_ineligible",
        "source_inputs": inventory,
        "recipe_inputs": {
            "replacement_path": relpath(REPLACEMENT_RECIPE),
            "replacement_sha256": sha256_file(REPLACEMENT_RECIPE),
            "replacement_rows": len(replacement_recipe),
            "marker_path": relpath(MARKER_RECIPE),
            "marker_sha256": sha256_file(MARKER_RECIPE),
            "marker_rows": len(marker_recipe),
            "cep_boundary_path": relpath(CEP_BOUNDARY),
            "cep_boundary_sha256": sha256_file(CEP_BOUNDARY),
            "cep_boundary_rows": len(cep_rows),
        },
        "derived_counts": {
            "replacement_rows": len(replacements),
            "inherited_replacement_rows": len(inherited_replacements),
            "review_adjustment_rows": len(review_adjustments),
            "marker_recipe_rows": len(marker_recipe),
            "marker_render_rows": len(rendered_markers),
            "marker_audit_rows": len(markers),
            "marker_suppressed_rows": len(markers) - len(rendered_markers),
            "cep_slice_audit_rows": len(cep_slice_audit),
            "cep_marker_render_rows": len(cep_markers),
            "publication_marker_render_rows": len(publication_markers),
            "phase_switch_midpoint_rows": len(publication_markers),
            "marker_semantics_audit_rows": len(marker_semantics),
            "publication_clean_point_rows": len(points),
            "curve_rows": len(curves),
            "figure_count": len(figure_paths),
        },
        "known_boundaries": [
            "mode_a xi=0 first-order recipe retains raw xi=-0.003/+0.003 endpoints and renders their interval midpoint as a display-only marker",
            "mode_a muB=900 alpha_T=1.0 xi=-0.01 mixed bulk/equilibrium branch remains excluded",
            "T=200 mode-B xi=-0.10 (muB=900) and xi=0.36 (muB=0) have author-requested display-only smoothing candidates",
            "T=120 mode-B muB=900 xi=-0.09 historical first-order marker is suppressed; publication uses the author-confirmed phase-label-switch interval midpoint",
            "mode_a muB=450 alpha_T=1.0 xi=-0.20 retains the raw non-first-order slope change; prior mechanism evidence is simple_1m4KPi",
            "derived display values are not solver recomputations or numerical convergence evidence",
            "old prod_v1 and phase-reference legacy fallback are retained; retirement is a separate audit",
        ],
        "outputs": [
            {"path": relpath(path), "sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in output_paths
        ],
    }
    write_json(OUT_DIR / "manifest.json", manifest)
    print(json.dumps({"output": relpath(OUT_DIR), "manifest": relpath(OUT_DIR / "manifest.json"), "figures": len(figure_paths)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
