#!/usr/bin/env python3
"""Build a non-destructive pole-sensitive rendering audit for Issue #130.

The script consumes the merged v2 production artifact and the existing v1
denominator-chain analysis. It verifies that relaxation-time and channel-rate
semantics are unchanged within a strict transfer gate, then produces derived
display candidates under docs/analysis/. Production CSVs, canonical figures,
and the production registry are read-only inputs.
"""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import math
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[3]
V1_CASE = "first_canonical_v1_p128_xi001_validated_anchored_prod_v1"
V2_CASE = "first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1"
RESULT_ROOT = (
    ROOT / "data" / "outputs" / "results" / "relaxtime" / "transport" / "phase_guided"
)
FIGURE_ROOT = (
    ROOT / "data" / "outputs" / "figures" / "relaxtime" / "transport" / "phase_guided"
)
SOURCE_ANALYSIS = (
    ROOT / "docs" / "analysis" / "relaxtime" / "phase_guided_transport_p128_xi001_analysis"
)
OUT_DIR = (
    ROOT
    / "docs"
    / "analysis"
    / "relaxtime"
    / "phase_guided_transport_v2_pole_sensitive_rendering"
)
TABLE_DIR = OUT_DIR / "tables"
DERIVED_FIGURE_DIR = OUT_DIR / "figures"
PAPER_FIGURE_DIR = OUT_DIR / "paper_figures"

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

TAU_FIELDS = ["tau_u", "tau_d", "tau_s", "tau_ubar", "tau_dbar", "tau_sbar"]
TRANSPORT_FIELDS = ["eta", "eta_over_s", "zeta", "zeta_over_s", "sigma", "sigma_over_T"]
DISPLAY_FIELDS = ["eta_over_s", "zeta_over_s", "sigma_over_T"]
DISPLAY_LABELS = {
    "tau_u": r"$\tau_u$",
    "tau_d": r"$\tau_d$",
    "tau_s": r"$\tau_s$",
    "tau_ubar": r"$\tau_{\bar u}$",
    "tau_dbar": r"$\tau_{\bar d}$",
    "tau_sbar": r"$\tau_{\bar s}$",
    "eta_over_s": r"$\eta/s$",
    "zeta_over_s": r"$\zeta/s$",
    "sigma_over_T": r"$\sigma/T$",
}
PAPER_COLORS = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE"]
TRANSFER_TOL = 1.0e-3
DISPLAY_LOG_THRESHOLD = 0.03


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


def read_csv_with_comments(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        lines = [line for line in handle if line.strip() and not line.startswith("#")]
    return list(csv.DictReader(lines))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")


def f(row: dict[str, str], field: str) -> float:
    value = float(row[field])
    if not math.isfinite(value):
        raise ValueError(f"non-finite {field}: {row[field]}")
    return value


def finite_float(value: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError(f"non-finite value: {value}")
    return parsed


def signed_log_delta(baseline: float, target: float) -> float:
    if baseline <= 0 or target <= 0:
        raise ValueError(f"positive values required: baseline={baseline}, target={target}")
    return math.log(target / baseline)


def relative_delta(a: float, b: float) -> float:
    return abs(b - a) / max(abs(a), abs(b), 1.0e-300)


def row_key(row: dict[str, str]) -> tuple[str, str, str]:
    return row["plot_panel"], row["plot_series"], f"{f(row, 'xi'):.10f}"


def diag_key(row: dict[str, str]) -> tuple[str, ...]:
    return (
        row["T_MeV"],
        row["muB_MeV"],
        row["xi"],
        row["species"],
        row["channel"],
    )


def case_paths(mode: str, case: str) -> dict[str, Path]:
    result_dir = RESULT_ROOT / mode / case
    figure_dir = FIGURE_ROOT / mode / case
    return {
        "result_dir": result_dir,
        "figure_dir": figure_dir,
        "scan": result_dir / "phase_guided_transport_scan.csv",
        "diagnostics": result_dir / "channel_diagnostics.csv",
        "failed": result_dir / "failed_points.csv",
        "manifest": result_dir / "manifest.json",
        "effective_config": result_dir / "effective_config.json",
        "audit": result_dir / "PRODUCTION_AUDIT.md",
        "plot_manifest": figure_dir / "plot_manifest.json",
        "comparison": result_dir / "convergence" / "old_vs_new_full_grid_comparison.csv",
    }


def load_and_validate_inputs() -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]]]:
    loaded: dict[str, dict[str, Any]] = {}
    inventory: list[dict[str, Any]] = []
    for mode_key, cfg in MODE_CONFIG.items():
        v2_paths = case_paths(cfg["mode"], V2_CASE)
        v1_paths = case_paths(cfg["mode"], V1_CASE)
        for name, path in v2_paths.items():
            if name in {"result_dir", "figure_dir"}:
                continue
            if not path.exists():
                raise FileNotFoundError(f"missing v2 {mode_key} {name}: {path}")
        for name in ("scan", "diagnostics"):
            if not v1_paths[name].exists():
                raise FileNotFoundError(f"missing v1 {mode_key} {name}: {v1_paths[name]}")

        scan_rows = read_csv_with_comments(v2_paths["scan"])
        diag_rows = read_csv_with_comments(v2_paths["diagnostics"])
        failed_rows = read_csv_with_comments(v2_paths["failed"])
        comparison_rows = read_csv_with_comments(v2_paths["comparison"])
        old_diag_rows = read_csv_with_comments(v1_paths["diagnostics"])

        if len(scan_rows) != 909:
            raise ValueError(f"{mode_key} scan rows: {len(scan_rows)} != 909")
        if len(diag_rows) != 38178:
            raise ValueError(f"{mode_key} diagnostic rows: {len(diag_rows)} != 38178")
        if failed_rows:
            raise ValueError(f"{mode_key} has {len(failed_rows)} failed rows")
        keys = [row_key(row) for row in scan_rows]
        if len(set(keys)) != len(keys):
            raise ValueError(f"duplicate scan keys in {mode_key}")

        for row in scan_rows:
            for field in ["xi", *TAU_FIELDS, *TRANSPORT_FIELDS, "m_u", "Phi", "s_fm3inv"]:
                f(row, field)
        for row in diag_rows:
            for field in ("density", "rate", "contribution", "total", "tau_inv_species"):
                value = f(row, field)
                if value < 0:
                    raise ValueError(f"negative diagnostic {mode_key} {field}: {value}")

        tau_comparison = [
            finite_float(row["rel_delta"])
            for row in comparison_rows
            if row["field"].startswith("tau_")
        ]
        max_tau_rel_drift = max(tau_comparison)
        old_by_key = {diag_key(row): row for row in old_diag_rows}
        if len(old_by_key) != len(old_diag_rows):
            raise ValueError(f"duplicate v1 diagnostic keys in {mode_key}")
        rate_drifts: list[float] = []
        for row in diag_rows:
            key = diag_key(row)
            if key not in old_by_key:
                raise ValueError(f"v2 diagnostic key absent from v1: {mode_key} {key}")
            rate_drifts.append(relative_delta(f(old_by_key[key], "rate"), f(row, "rate")))
        max_rate_rel_drift = max(rate_drifts)
        if max_tau_rel_drift > TRANSFER_TOL or max_rate_rel_drift > TRANSFER_TOL:
            raise ValueError(
                f"mechanism transfer gate failed for {mode_key}: "
                f"tau={max_tau_rel_drift}, rate={max_rate_rel_drift}"
            )

        loaded[mode_key] = {
            "config": cfg,
            "v2_paths": v2_paths,
            "v1_paths": v1_paths,
            "scan_rows": scan_rows,
            "diag_rows": diag_rows,
            "scan_index": {row_key(row): row for row in scan_rows},
            "curves": group_curves(scan_rows),
            "max_tau_rel_drift": max_tau_rel_drift,
            "max_rate_rel_drift": max_rate_rel_drift,
        }
        inventory.append(
            {
                "mode_key": mode_key,
                "mode": cfg["mode"],
                "scan_rows": len(scan_rows),
                "diagnostic_rows": len(diag_rows),
                "failed_rows": len(failed_rows),
                "xi_count": len({round(f(row, "xi"), 10) for row in scan_rows}),
                "scan_sha256": sha256_file(v2_paths["scan"]),
                "diagnostics_sha256": sha256_file(v2_paths["diagnostics"]),
                "manifest_sha256": sha256_file(v2_paths["manifest"]),
                "plot_manifest_sha256": sha256_file(v2_paths["plot_manifest"]),
                "max_tau_v1_v2_rel_drift": max_tau_rel_drift,
                "max_rate_v1_v2_rel_drift": max_rate_rel_drift,
                "transfer_gate_tol": TRANSFER_TOL,
                "transfer_gate_passed": True,
            }
        )
    return loaded, inventory


def group_curves(rows: list[dict[str, str]]) -> dict[tuple[str, str], list[dict[str, str]]]:
    curves: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        curves[(row["plot_panel"], row["plot_series"])].append(row)
    for curve in curves.values():
        curve.sort(key=lambda row: f(row, "xi"))
    return dict(curves)


def source_window_tables() -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    window_path = SOURCE_ANALYSIS / "tables" / "tau_jump_window_summary.csv"
    mechanism_path = SOURCE_ANALYSIS / "tables" / "mechanism_window_summary.csv"
    window_rows = read_csv_with_comments(window_path)
    mechanism_rows = read_csv_with_comments(mechanism_path)
    if len(window_rows) != 10:
        raise ValueError(f"expected 10 classified source windows, got {len(window_rows)}")
    if len(mechanism_rows) != 8:
        raise ValueError(f"expected 8 denominator windows, got {len(mechanism_rows)}")
    if any(row["mechanism_verdict"] != "small_denominator_supported" for row in mechanism_rows):
        raise ValueError("source mechanism table contains a non-supported denominator window")
    return window_rows, mechanism_rows


def curve_row(
    loaded: dict[str, dict[str, Any]], mode_key: str, panel: str, series: str, xi: float
) -> dict[str, str]:
    key = (panel, series, f"{xi:.10f}")
    try:
        return loaded[mode_key]["scan_index"][key]
    except KeyError as exc:
        raise KeyError(f"missing v2 scan point {mode_key} {key}") from exc


def build_window_and_point_audits(
    loaded: dict[str, dict[str, Any]],
    source_windows: list[dict[str, str]],
    mechanism_rows: list[dict[str, str]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    source_by_id = {row["window_id"]: row for row in source_windows}
    mechanism_by_id = {row["window_id"]: row for row in mechanism_rows}
    classifications: list[dict[str, Any]] = []
    point_audit: list[dict[str, Any]] = []
    mask_rows: list[dict[str, Any]] = []
    first_order_rows: list[dict[str, Any]] = []

    ordered_ids = [row["window_id"] for row in mechanism_rows] + [
        row["window_id"]
        for row in source_windows
        if row["cause_verdict"] == "upstream_first_order_branch_jump_supported"
    ]
    if len(ordered_ids) != 10 or len(set(ordered_ids)) != 10:
        raise ValueError("window classification did not resolve to 8 pole + 2 first-order windows")

    for window_id in ordered_ids:
        source = source_by_id[window_id]
        mechanism = mechanism_by_id.get(window_id)
        mode_key = source["mode_key"]
        panel = source["plot_panel"]
        series = source["plot_series"]
        target_xi = finite_float(source["target_xi"])
        prev_xi = finite_float(source["prev_xi"])
        next_xi = finite_float(source["next_xi"])
        prev = curve_row(loaded, mode_key, panel, series, prev_xi)
        target = curve_row(loaded, mode_key, panel, series, target_xi)
        next_row = curve_row(loaded, mode_key, panel, series, next_xi)
        affected_tau = [value for value in source["affected_tau_fields"].split(";") if value]
        is_pole = mechanism is not None
        classification = "pole_sensitive_supported" if is_pole else "first_order_protected"
        primary_observable = mechanism["observable"] if is_pole else "tau_u"
        dominant_branch = (
            mechanism["dominant_denominator_branch"] if is_pole else "upstream_equilibrium_branch"
        )
        classifications.append(
            {
                "window_id": window_id,
                "mode_key": mode_key,
                "plot_panel": panel,
                "plot_series": series,
                "prev_xi": prev_xi,
                "target_xi": target_xi,
                "next_xi": next_xi,
                "classification": classification,
                "primary_observable": primary_observable,
                "affected_tau_fields": ";".join(affected_tau),
                "dominant_branch": dominant_branch,
                "upstream_branch_flag": not is_pole,
                "display_mask_eligible": is_pole,
                "evidence_source": relpath(
                    SOURCE_ANALYSIS / "tables" / (
                        "mechanism_window_summary.csv" if is_pole else "tau_jump_window_summary.csv"
                    )
                ),
            }
        )
        if not is_pole:
            first_order_rows.append(
                {
                    "window_id": window_id,
                    "mode_key": mode_key,
                    "plot_panel": panel,
                    "plot_series": series,
                    "target_xi": target_xi,
                    "protected_observables": ";".join([*TAU_FIELDS, *TRANSPORT_FIELDS]),
                    "phase_structure": target["phase_structure"],
                    "phase_reference_kind": target["phase_reference_kind"],
                    "max_background_rel_step": source["max_background_rel_step"],
                    "render_policy": "retain_raw_values_and_continuous_line_no_bridge",
                    "protection_reason": "upstream first-order/background branch jump",
                }
            )

        for observable in [*TAU_FIELDS, *TRANSPORT_FIELDS]:
            prev_value = f(prev, observable)
            target_value = f(target, observable)
            next_value = f(next_row, observable)
            interpolation_weight = (target_xi - prev_xi) / (next_xi - prev_xi)
            guide_value = prev_value + interpolation_weight * (next_value - prev_value)
            log_delta = signed_log_delta(guide_value, target_value)
            tau_candidate = observable in affected_tau
            transport_candidate = (
                observable in TRANSPORT_FIELDS and abs(log_delta) >= DISPLAY_LOG_THRESHOLD
            )
            mask_candidate = is_pole and (tau_candidate or transport_candidate)
            audit_row = {
                "window_id": window_id,
                "mode_key": mode_key,
                "plot_panel": panel,
                "plot_series": series,
                "observable": observable,
                "prev_xi": prev_xi,
                "target_xi": target_xi,
                "next_xi": next_xi,
                "prev_value": prev_value,
                "raw_target_value": target_value,
                "next_value": next_value,
                "linear_neighbor_guide_value": guide_value,
                "log_delta_vs_neighbor_mean": log_delta,
                "relative_delta_vs_neighbor_mean": relative_delta(guide_value, target_value),
                "classification": classification,
                "mask_candidate": mask_candidate,
                "mask_reason": (
                    "mechanism-supported pole-sensitive tau outlier"
                    if mask_candidate and observable in TAU_FIELDS
                    else "mechanism-supported pole window with >=3% downstream log deviation"
                    if mask_candidate
                    else "first-order protected"
                    if not is_pole
                    else "below derived-display threshold"
                ),
                "value_semantics": "raw production value; guide value is not a recomputation",
            }
            point_audit.append(audit_row)
            if mask_candidate:
                mask_rows.append(
                    {
                        **audit_row,
                        "render_action": "exclude raw point from solid guide and draw dashed neighbor bridge",
                        "policy_status": "derived_display_candidate",
                        "canonical_data_modified": False,
                    }
                )

    return classifications, point_audit, mask_rows, first_order_rows


def set_local_xlim(ax: plt.Axes, target_xi: float, half_width: float) -> None:
    ax.set_xlim(max(-0.5, target_xi - half_width), min(0.5, target_xi + half_width))


def render_pole_window(
    loaded: dict[str, dict[str, Any]],
    window: dict[str, Any],
    point_audit: list[dict[str, Any]],
) -> Path:
    mode_key = window["mode_key"]
    panel = window["plot_panel"]
    series = window["plot_series"]
    target_xi = float(window["target_xi"])
    curve = loaded[mode_key]["curves"][(panel, series)]
    audit_by_obs = {
        row["observable"]: row
        for row in point_audit
        if row["window_id"] == window["window_id"]
    }
    observables = [window["primary_observable"], *DISPLAY_FIELDS]
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 8.5))
    color = "#276FBF"
    raw_color = "#9AA0A6"
    flag_color = "#D97706"
    for ax, observable in zip(axes.flat, observables):
        xs = [f(row, "xi") for row in curve]
        ys = [f(row, observable) for row in curve]
        audit = audit_by_obs[observable]
        ax.plot(xs, ys, color=raw_color, linewidth=1.5, alpha=0.72)
        guide = list(ys)
        target_index = min(range(len(xs)), key=lambda idx: abs(xs[idx] - target_xi))
        if bool(audit["mask_candidate"]):
            guide[target_index] = math.nan
        ax.plot(xs, guide, color=color, linewidth=2.2)
        ax.axvspan(target_xi - 0.004, target_xi + 0.004, color=flag_color, alpha=0.10)
        ax.scatter(
            [target_xi],
            [float(audit["raw_target_value"])],
            marker="x",
            s=70,
            linewidths=2.0,
            color=flag_color,
            zorder=5,
        )
        if bool(audit["mask_candidate"]):
            ax.plot(
                [float(audit["prev_xi"]), float(audit["next_xi"])],
                [float(audit["prev_value"]), float(audit["next_value"])],
                color=flag_color,
                linewidth=2.0,
                linestyle="--",
                zorder=4,
            )
            ax.scatter(
                [target_xi],
                [float(audit["linear_neighbor_guide_value"])],
                facecolors="white",
                edgecolors=flag_color,
                s=45,
                zorder=6,
            )
        set_local_xlim(ax, target_xi, 0.08)
        ax.set_xlabel(r"$\xi$")
        ax.set_ylabel(DISPLAY_LABELS[observable])
        ax.grid(alpha=0.22)
        ax.set_title(
            f"{observable}: raw deviation={float(audit['relative_delta_vs_neighbor_mean']):.3g}"
        )
    fig.suptitle(
        f"Pole-sensitive derived display candidate\n"
        f"{mode_key}, {panel}, {series}, xi={target_xi:.2f}, {window['dominant_branch']}",
        fontsize=14,
        y=0.985,
    )
    handles = [
        Line2D([0], [0], color=raw_color, lw=1.5, label="raw production curve"),
        Line2D([0], [0], color=color, lw=2.2, label="solid guide excluding eligible point"),
        Line2D([0], [0], color=flag_color, marker="x", lw=0, label="raw pole-sensitive sample"),
        Line2D([0], [0], color=flag_color, lw=2.0, ls="--", label="neighbor bridge (display only)"),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.012),
        ncol=2,
        frameon=False,
    )
    fig.tight_layout(rect=(0.0, 0.09, 1.0, 0.93))
    path = DERIVED_FIGURE_DIR / f"pole_sensitive_{window['window_id']}.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def render_first_order_window(
    loaded: dict[str, dict[str, Any]], window: dict[str, Any]
) -> Path:
    mode_key = window["mode_key"]
    panel = window["plot_panel"]
    series = window["plot_series"]
    target_xi = float(window["target_xi"])
    curve = loaded[mode_key]["curves"][(panel, series)]
    observables = ["tau_u", *DISPLAY_FIELDS]
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 8.5), constrained_layout=True)
    color = "#5B3F95"
    for ax, observable in zip(axes.flat, observables):
        xs = [f(row, "xi") for row in curve]
        ys = [f(row, observable) for row in curve]
        target_index = min(range(len(xs)), key=lambda idx: abs(xs[idx] - target_xi))
        ax.plot(xs, ys, color=color, linewidth=2.2)
        ax.axvline(target_xi, color="#2E1065", linewidth=1.8, linestyle=":")
        ax.scatter(
            [target_xi], [ys[target_index]], marker="D", s=55, color="#2E1065", zorder=5
        )
        set_local_xlim(ax, target_xi, 0.12)
        ax.set_xlabel(r"$\xi$")
        ax.set_ylabel(DISPLAY_LABELS[observable])
        ax.grid(alpha=0.22)
        ax.set_title(f"{observable}: retained raw first-order response")
    fig.suptitle(
        f"Protected first-order/background branch window\n"
        f"{mode_key}, {panel}, {series}, xi={target_xi:.2f}; no display mask applied",
        fontsize=14,
    )
    path = DERIVED_FIGURE_DIR / f"first_order_protected_{window['window_id']}.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def paper_replacement_rows(mask_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for row in mask_rows:
        if row["observable"] not in DISPLAY_FIELDS:
            continue
        rows.append(
            {
                "window_id": row["window_id"],
                "mode_key": row["mode_key"],
                "plot_panel": row["plot_panel"],
                "plot_series": row["plot_series"],
                "observable": row["observable"],
                "xi": row["target_xi"],
                "raw_production_value": row["raw_target_value"],
                "paper_display_value": row["linear_neighbor_guide_value"],
                "left_xi": row["prev_xi"],
                "left_value": row["prev_value"],
                "right_xi": row["next_xi"],
                "right_value": row["next_value"],
                "replacement_method": "linear interpolation between adjacent unmasked production samples",
                "figure_semantics": "clean continuous line; no correction marker in paper candidate",
                "canonical_data_modified": False,
            }
        )
    return rows


def paper_marker_rows(
    loaded: dict[str, dict[str, Any]], first_order_rows: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for window in first_order_rows:
        target = curve_row(
            loaded,
            window["mode_key"],
            window["plot_panel"],
            window["plot_series"],
            float(window["target_xi"]),
        )
        for observable in DISPLAY_FIELDS:
            rows.append(
                {
                    "window_id": window["window_id"],
                    "mode_key": window["mode_key"],
                    "plot_panel": window["plot_panel"],
                    "plot_series": window["plot_series"],
                    "observable": observable,
                    "xi": window["target_xi"],
                    "raw_production_value": f(target, observable),
                    "marker": "star",
                    "marker_semantics": "first-order/upstream branch transition point retained without smoothing",
                }
            )
    return rows


def configure_paper_style() -> None:
    matplotlib.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "font.size": 10,
            "mathtext.fontset": "stix",
            "axes.labelsize": 10,
            "axes.linewidth": 0.6,
            "legend.fontsize": 8,
            "legend.frameon": False,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def render_paper_figures(
    loaded: dict[str, dict[str, Any]],
    replacements: list[dict[str, Any]],
    markers: list[dict[str, Any]],
) -> list[Path]:
    configure_paper_style()
    replacement_index = {
        (
            row["mode_key"],
            row["plot_panel"],
            row["plot_series"],
            row["observable"],
            f"{float(row['xi']):.10f}",
        ): float(row["paper_display_value"])
        for row in replacements
    }
    marker_index = {
        (
            row["mode_key"],
            row["plot_panel"],
            row["plot_series"],
            row["observable"],
        ): row
        for row in markers
    }
    paths: list[Path] = []
    for mode_key, payload in loaded.items():
        panels = sorted({panel for panel, _ in payload["curves"]})
        for panel in panels:
            series_names = sorted(
                series for current_panel, series in payload["curves"] if current_panel == panel
            )
            for observable in DISPLAY_FIELDS:
                fig, ax = plt.subplots(figsize=(6.75, 4.6))
                marker_present = False
                for series_index, series in enumerate(series_names):
                    curve = payload["curves"][(panel, series)]
                    color = PAPER_COLORS[series_index % len(PAPER_COLORS)]
                    xs = [f(row, "xi") for row in curve]
                    ys: list[float] = []
                    for row in curve:
                        key = (
                            mode_key,
                            panel,
                            series,
                            observable,
                            f"{f(row, 'xi'):.10f}",
                        )
                        ys.append(replacement_index.get(key, f(row, observable)))
                    ax.plot(
                        xs,
                        ys,
                        color=color,
                        linewidth=1.5,
                        label=curve[0]["plot_series_label"],
                    )
                    marker = marker_index.get((mode_key, panel, series, observable))
                    if marker is not None:
                        marker_present = True
                        ax.scatter(
                            [float(marker["xi"])],
                            [float(marker["raw_production_value"])],
                            marker="*",
                            s=95,
                            facecolor=color,
                            edgecolor="black",
                            linewidth=0.55,
                            zorder=5,
                        )
                if marker_present:
                    ax.scatter(
                        [],
                        [],
                        marker="*",
                        s=75,
                        facecolor="white",
                        edgecolor="black",
                        linewidth=0.55,
                        label="first-order transition",
                    )
                ax.set_xlabel(r"$\xi$")
                ax.set_ylabel(DISPLAY_LABELS[observable])
                ax.set_xlim(-0.52, 0.52)
                ax.legend(loc="best")
                fig.tight_layout()
                path = (
                    PAPER_FIGURE_DIR
                    / mode_key
                    / f"plot_panel={panel}"
                    / f"{observable}_vs_xi.png"
                )
                path.parent.mkdir(parents=True, exist_ok=True)
                fig.savefig(path, dpi=600, bbox_inches="tight", pad_inches=0.08)
                plt.close(fig)
                paths.append(path)
    return paths


def build_claim_ledger() -> list[dict[str, Any]]:
    return [
        {
            "claim_id": "CLAIM-V2-POLE-001",
            "status": "supported_with_scope_limit",
            "claim_zh": "v1 已定点深拆的 8 个小分母窗口可映射到 v2，因为 tau 与通道率均通过 1e-3 机制迁移门槛。",
            "evidence": "tables/input_inventory.csv; tables/window_classification.csv",
            "scope_limit": "迁移门槛不替代局部 high-rate convergence gate。",
        },
        {
            "claim_id": "CLAIM-V2-POLE-002",
            "status": "supported",
            "claim_zh": "两个一阶/上游分支窗口被硬保护，不进入 pole-sensitive display mask。",
            "evidence": "tables/first_order_protection.csv",
            "scope_limit": "仅覆盖现有标量输运产物。",
        },
        {
            "claim_id": "CLAIM-V2-POLE-003",
            "status": "supported",
            "claim_zh": "虚线邻点桥接仅用于项目内部审计，不进入论文候选图；它不是重算或物理正则化结果。",
            "evidence": "tables/pole_sensitive_mask.csv; figures/plot_manifest.json",
            "scope_limit": "内部审计图不得作为未注明处理的论文曲线。",
        },
        {
            "claim_id": "CLAIM-V2-POLE-004",
            "status": "author_check",
            "claim_zh": "若需要物理意义上的平滑曲线，应先定义有限宽度或极点正则化并用新 slug 重跑。",
            "evidence": "README.md",
            "scope_limit": "本分析 PR 不修改传播子模型。",
        },
        {
            "claim_id": "CLAIM-V2-POLE-005",
            "status": "author_directed_candidate",
            "claim_zh": "论文候选图以连续实线使用相邻真实点替换已确认的数值噪点，不显示修正标记；一阶相变点以星号标注并保留原始物理跳变。",
            "evidence": "tables/paper_display_replacements.csv; tables/paper_first_order_markers.csv; paper_figures/plot_manifest.json",
            "scope_limit": "仅适用于本表列出的 13 个派生显示替换值与 6 个一阶标记点。",
        },
    ]


def render_readme(
    inventory: list[dict[str, Any]],
    classifications: list[dict[str, Any]],
    mask_rows: list[dict[str, Any]],
    first_order_rows: list[dict[str, Any]],
    figure_paths: list[Path],
    paper_figure_paths: list[Path],
    paper_replacements: list[dict[str, Any]],
    paper_markers: list[dict[str, Any]],
) -> str:
    pole_windows = [row for row in classifications if row["classification"] == "pole_sensitive_supported"]
    mask_by_window: dict[str, list[str]] = defaultdict(list)
    for row in mask_rows:
        mask_by_window[row["window_id"]].append(row["observable"])
    inventory_lines = "\n".join(
        f"| {row['mode_key']} | {row['scan_rows']} | {row['diagnostic_rows']} | "
        f"{row['max_tau_v1_v2_rel_drift']:.8g} | {row['max_rate_v1_v2_rel_drift']:.8g} |"
        for row in inventory
    )
    pole_lines = "\n".join(
        f"| {row['window_id']} | {row['plot_panel']} | {row['plot_series']} | "
        f"{float(row['target_xi']):.2f} | {row['dominant_branch']} | "
        f"{';'.join(mask_by_window[row['window_id']])} |"
        for row in pole_windows
    )
    first_order_lines = "\n".join(
        f"| {row['window_id']} | {row['plot_panel']} | {row['plot_series']} | "
        f"{float(row['target_xi']):.2f} | {row['phase_structure']} |"
        for row in first_order_rows
    )
    figure_lines = "\n".join(f"- `{relpath(path)}`" for path in figure_paths)
    paper_figure_lines = "\n".join(
        f"- `{relpath(path)}`" for path in paper_figure_paths
    )
    return f"""# Phase-guided transport v2 极点敏感派生显示审计

## 范围与结论

本分析包消费 `{V2_CASE}`，不修改 production CSV、canonical figure 或 `production_registry.json`。它把旧 xi001 denominator-chain 机制表迁移到 v2，并生成**仅供作者审阅**的非破坏性派生显示候选。

当前结论：8 个窗口通过小分母机制迁移门槛；2 个一阶/上游分支窗口被硬保护。虚线桥接是视觉指南，不是新计算值，也不能静默替换正式图。

## 输入与迁移门槛

| mode | scan rows | diagnostic rows | max tau v1/v2 drift | max rate v1/v2 drift |
| --- | ---: | ---: | ---: | ---: |
{inventory_lines}

迁移门槛为 `{TRANSFER_TOL}`。v2 的能量语义修改不改变弛豫时间定义；实际 tau 和 channel-rate 漂移均显著低于门槛，因此允许继承已有 denominator-chain 定点证据。该门槛不替代尚未完成的局部 high-rate convergence gate。

## Pole-sensitive display mask 候选

| window | panel | series | xi | denominator branch | eligible observables |
| --- | --- | --- | ---: | --- | --- |
{pole_lines}

规则：

1. tau 只有在旧分析的 `affected_tau_fields` 中才进入 mask；
2. transport 只有在目标点相对相邻两点均值的绝对 log 偏离不小于 `{DISPLAY_LOG_THRESHOLD}` 时才进入 mask；
3. raw 点始终以橙色叉号保留；solid guide 在 eligible 点断开，并以虚线连接两侧真实邻点；
4. `linear_neighbor_guide_value` 只用于画图，不写回 production，也不声称是物理值。

## 一阶窗口硬保护

| window | panel | series | xi | phase structure |
| --- | --- | --- | ---: | --- |
{first_order_lines}

这两个窗口保留 raw 曲线和跳变，不应用桥接或 mask。保护依据是背景质量、Polyakov loop 和熵密度的同步快速变化，而不是单纯依赖 phase 字符串标签。

## 派生图

{figure_lines}

每个 pole-sensitive 图包含 primary tau 与 `eta_over_s/zeta_over_s/sigma_over_T` 四个局部面板。每个一阶图包含相同的下游比值与 `tau_u`，用于核对跳变未被删除。

## 论文候选图

{paper_figure_lines}

论文候选图沿用正式图的固定 panel、多曲线布局，只绘制 `eta_over_s`、`zeta_over_s` 和 `sigma_over_T`：

1. 共生成 `{len(paper_figure_paths)}` 张 600 DPI 图；
2. `{len(paper_replacements)}` 个已确认的小分母下游噪点在派生绘图值中由左右相邻真实样本线性插值替换，图面只显示正常连续实线，不显示叉号、空心点、虚线或修正标签；
3. `{len(paper_markers)}` 个 observable-level 一阶位置以星号标在对应曲线上，不使用竖直虚线；
4. 星号处仍使用 raw production 值，一阶/上游分支跳变没有被平滑；
5. 替换值与星号位置分别记录在 `tables/paper_display_replacements.csv` 和 `tables/paper_first_order_markers.csv`。图面不呈现内部修正痕迹，但仓库内保留完整可追溯记录。

## 证据边界与作者判断

- `supported`：v2 输入完整、无 failed/NaN/负 rate，v1→v2 tau/rate 迁移门槛通过；一阶窗口保护规则明确。
- `supported_with_scope_limit`：8 个窗口已有 denominator-chain 与生产 rate 复现证据，但没有逐窗口新的 high-rate convergence gate。
- `author_directed_candidate`：论文候选图按作者约定隐藏数值修正痕迹，并用星号标示一阶相变位置；仍需最终视觉审核后决定是否采用。
- `author_check`：是否需要计算层有限宽度/极点正则化和新 production slug。
- 本包不把小分母结构直接定性为随机数值噪声，也不把虚线桥接升级为物理预测。

## 复现

```powershell
python scripts/analysis/relaxtime/build_phase_guided_pole_sensitive_rendering.py
```

关键表：

- `tables/input_inventory.csv`
- `tables/window_classification.csv`
- `tables/rendered_point_audit.csv`
- `tables/pole_sensitive_mask.csv`
- `tables/first_order_protection.csv`
- `tables/paper_display_replacements.csv`
- `tables/paper_first_order_markers.csv`
- `tables/claim_ledger.csv`
- `figures/plot_manifest.json`
- `paper_figures/plot_manifest.json`
"""


def main() -> None:
    loaded, inventory = load_and_validate_inputs()
    source_windows, mechanism_rows = source_window_tables()
    classifications, point_audit, mask_rows, first_order_rows = build_window_and_point_audits(
        loaded, source_windows, mechanism_rows
    )

    write_csv(
        TABLE_DIR / "input_inventory.csv",
        inventory,
        [
            "mode_key",
            "mode",
            "scan_rows",
            "diagnostic_rows",
            "failed_rows",
            "xi_count",
            "scan_sha256",
            "diagnostics_sha256",
            "manifest_sha256",
            "plot_manifest_sha256",
            "max_tau_v1_v2_rel_drift",
            "max_rate_v1_v2_rel_drift",
            "transfer_gate_tol",
            "transfer_gate_passed",
        ],
    )
    write_csv(
        TABLE_DIR / "window_classification.csv",
        classifications,
        [
            "window_id",
            "mode_key",
            "plot_panel",
            "plot_series",
            "prev_xi",
            "target_xi",
            "next_xi",
            "classification",
            "primary_observable",
            "affected_tau_fields",
            "dominant_branch",
            "upstream_branch_flag",
            "display_mask_eligible",
            "evidence_source",
        ],
    )
    point_fields = [
        "window_id",
        "mode_key",
        "plot_panel",
        "plot_series",
        "observable",
        "prev_xi",
        "target_xi",
        "next_xi",
        "prev_value",
        "raw_target_value",
        "next_value",
        "linear_neighbor_guide_value",
        "log_delta_vs_neighbor_mean",
        "relative_delta_vs_neighbor_mean",
        "classification",
        "mask_candidate",
        "mask_reason",
        "value_semantics",
    ]
    write_csv(TABLE_DIR / "rendered_point_audit.csv", point_audit, point_fields)
    write_csv(
        TABLE_DIR / "pole_sensitive_mask.csv",
        mask_rows,
        [*point_fields, "render_action", "policy_status", "canonical_data_modified"],
    )
    write_csv(
        TABLE_DIR / "first_order_protection.csv",
        first_order_rows,
        [
            "window_id",
            "mode_key",
            "plot_panel",
            "plot_series",
            "target_xi",
            "protected_observables",
            "phase_structure",
            "phase_reference_kind",
            "max_background_rel_step",
            "render_policy",
            "protection_reason",
        ],
    )
    paper_replacements = paper_replacement_rows(mask_rows)
    paper_markers = paper_marker_rows(loaded, first_order_rows)
    write_csv(
        TABLE_DIR / "paper_display_replacements.csv",
        paper_replacements,
        [
            "window_id",
            "mode_key",
            "plot_panel",
            "plot_series",
            "observable",
            "xi",
            "raw_production_value",
            "paper_display_value",
            "left_xi",
            "left_value",
            "right_xi",
            "right_value",
            "replacement_method",
            "figure_semantics",
            "canonical_data_modified",
        ],
    )
    write_csv(
        TABLE_DIR / "paper_first_order_markers.csv",
        paper_markers,
        [
            "window_id",
            "mode_key",
            "plot_panel",
            "plot_series",
            "observable",
            "xi",
            "raw_production_value",
            "marker",
            "marker_semantics",
        ],
    )
    claims = build_claim_ledger()
    write_csv(
        TABLE_DIR / "claim_ledger.csv",
        claims,
        ["claim_id", "status", "claim_zh", "evidence", "scope_limit"],
    )

    figure_paths: list[Path] = []
    for window in classifications:
        if window["classification"] == "pole_sensitive_supported":
            figure_paths.append(render_pole_window(loaded, window, point_audit))
        else:
            figure_paths.append(render_first_order_window(loaded, window))

    if PAPER_FIGURE_DIR.exists():
        shutil.rmtree(PAPER_FIGURE_DIR)
    paper_figure_paths = render_paper_figures(
        loaded, paper_replacements, paper_markers
    )

    generator_path = Path(__file__).resolve()
    plot_manifest = {
        "schema": "phase_guided_v2_pole_sensitive_plot_manifest_v1",
        "case": V2_CASE,
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "base_git_commit": git_head(),
        "generator": relpath(generator_path),
        "generator_sha256": sha256_file(generator_path),
        "rendering_semantics": "raw points retained; eligible points excluded only from solid display guide; dashed neighbor bridge is not recomputed physics",
        "figures": [
            {"path": relpath(path), "sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in figure_paths
        ],
    }
    write_json(DERIVED_FIGURE_DIR / "plot_manifest.json", plot_manifest)
    paper_plot_manifest = {
        "schema": "phase_guided_v2_paper_candidate_plot_manifest_v1",
        "case": V2_CASE,
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "base_git_commit": git_head(),
        "generator": relpath(generator_path),
        "generator_sha256": sha256_file(generator_path),
        "observables": DISPLAY_FIELDS,
        "paper_display_replacement_count": len(paper_replacements),
        "first_order_marker_count": len(paper_markers),
        "rendering_semantics": "clean continuous paper-candidate lines use audited derived replacements without visible correction marks; stars retain raw first-order transition values; no vertical transition lines",
        "figures": [
            {"path": relpath(path), "sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in paper_figure_paths
        ],
    }
    write_json(PAPER_FIGURE_DIR / "plot_manifest.json", paper_plot_manifest)

    readme_path = OUT_DIR / "README.md"
    with readme_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(
            render_readme(
                inventory,
                classifications,
                mask_rows,
                first_order_rows,
                figure_paths,
                paper_figure_paths,
                paper_replacements,
                paper_markers,
            )
        )
    output_paths = [
        readme_path,
        *sorted(TABLE_DIR.glob("*.csv")),
        DERIVED_FIGURE_DIR / "plot_manifest.json",
        *figure_paths,
        PAPER_FIGURE_DIR / "plot_manifest.json",
        *paper_figure_paths,
    ]
    manifest = {
        "schema": "phase_guided_v2_pole_sensitive_analysis_manifest_v1",
        "case": V2_CASE,
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "base_git_commit": git_head(),
        "generator": relpath(generator_path),
        "generator_sha256": sha256_file(generator_path),
        "scope": "read-only production audit and explicitly labeled derived display candidates",
        "canonical_production_modified": False,
        "registry_modified": False,
        "window_counts": {
            "pole_sensitive_supported": sum(
                row["classification"] == "pole_sensitive_supported" for row in classifications
            ),
            "first_order_protected": sum(
                row["classification"] == "first_order_protected" for row in classifications
            ),
            "mask_rows": len(mask_rows),
            "paper_display_replacements": len(paper_replacements),
            "paper_first_order_markers": len(paper_markers),
            "paper_figures": len(paper_figure_paths),
        },
        "inputs": [
            {
                "mode_key": row["mode_key"],
                "scan_sha256": row["scan_sha256"],
                "diagnostics_sha256": row["diagnostics_sha256"],
                "max_tau_v1_v2_rel_drift": row["max_tau_v1_v2_rel_drift"],
                "max_rate_v1_v2_rel_drift": row["max_rate_v1_v2_rel_drift"],
            }
            for row in inventory
        ],
        "outputs": [
            {"path": relpath(path), "sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in output_paths
        ],
    }
    write_json(OUT_DIR / "manifest.json", manifest)
    print(
        f"wrote {relpath(OUT_DIR)}: {len(classifications)} windows, "
        f"{len(mask_rows)} mask rows, {len(figure_paths)} audit figures, "
        f"{len(paper_figure_paths)} paper figures"
    )


if __name__ == "__main__":
    main()
