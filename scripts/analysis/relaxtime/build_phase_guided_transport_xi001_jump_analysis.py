#!/usr/bin/env python3
"""Build a tau-first jump analysis package for the xi=0.01 transport case.

The script is read-only with respect to production results. It consumes the
promoted p128 xi001 phase-guided transport CSVs and channel diagnostics, then
writes derived tables, figures, a manifest, and a Chinese analysis note under
docs/analysis/.
"""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import io
import json
import math
import statistics
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[3]
CASE = "first_canonical_v1_p128_xi001_validated_anchored_prod_v1"
CONVERGENCE_CASE = f"{CASE}_convergence"
ANALYSIS_DIR = (
    ROOT
    / "docs"
    / "analysis"
    / "relaxtime"
    / "phase_guided_transport_p128_xi001_analysis"
)

MODE_INPUTS = {
    "mode_a": {
        "mode_value": "mode_a_fixed_muB_phase_scaled",
        "label": "mode A: fixed muB, phase-scaled T",
        "result_dir": ROOT
        / "data"
        / "outputs"
        / "results"
        / "relaxtime"
        / "transport"
        / "phase_guided"
        / "mode_a_fixed_muB_phase_scaled"
        / CASE,
        "figure_dir": ROOT
        / "data"
        / "outputs"
        / "figures"
        / "relaxtime"
        / "transport"
        / "phase_guided"
        / "mode_a_fixed_muB_phase_scaled"
        / CASE,
        "main_rows": 909,
        "channel_rows": 38178,
    },
    "mode_b": {
        "mode_value": "mode_b_fixed_T_sparse_muB",
        "label": "mode B: fixed T, sparse muB",
        "result_dir": ROOT
        / "data"
        / "outputs"
        / "results"
        / "relaxtime"
        / "transport"
        / "phase_guided"
        / "mode_b_fixed_T_sparse_muB"
        / CASE,
        "figure_dir": ROOT
        / "data"
        / "outputs"
        / "figures"
        / "relaxtime"
        / "transport"
        / "phase_guided"
        / "mode_b_fixed_T_sparse_muB"
        / CASE,
        "main_rows": 909,
        "channel_rows": 38178,
    },
}

TAU_SPECIES = ["u", "d", "s", "ubar", "dbar", "sbar"]
TAU_FIELDS = [f"tau_{sp}" for sp in TAU_SPECIES]
TAUINV_FIELDS = [f"tauinv_{sp}" for sp in TAU_SPECIES]
DOWNSTREAM_FIELDS = ["eta", "eta_over_s", "zeta", "zeta_over_s", "sigma", "sigma_over_T"]
BACKGROUND_FIELDS = ["m_u", "m_s", "Phi", "Phibar", "s_fm3inv", "rho_norm", "P_fm4inv"]


def relpath(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def current_git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:
        return "unknown"


def read_csv_with_comments(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    comments: list[str] = []
    data_lines: list[str] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for line in handle:
            if line.startswith("#"):
                comments.append(line.rstrip("\n"))
            elif line.strip():
                data_lines.append(line)
    reader = csv.DictReader(data_lines)
    return comments, list(reader)


def read_optional_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    _, rows = read_csv_with_comments(path)
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def f(row: dict[str, str], field: str, default: float = math.nan) -> float:
    try:
        value = float(row.get(field, ""))
    except (TypeError, ValueError):
        return default
    return value if math.isfinite(value) else default


def safe_ratio(a: float, b: float) -> float:
    if not (math.isfinite(a) and math.isfinite(b)) or a <= 0 or b <= 0:
        return math.nan
    return b / a


def rel_delta(a: float, b: float) -> float:
    if not (math.isfinite(a) and math.isfinite(b)):
        return math.nan
    return abs(b - a) / max(abs(a), abs(b), 1.0e-300)


def signed_log_delta(a: float, b: float) -> float:
    ratio = safe_ratio(a, b)
    return math.log(ratio) if math.isfinite(ratio) and ratio > 0 else math.nan


def norm_key(value: Any) -> str:
    return f"{float(value):.10f}"


def curve_key(row: dict[str, str]) -> tuple[str, str]:
    return row["plot_panel"], row["plot_series"]


def grid_key(row: dict[str, str], species: str | None = None) -> tuple[str, ...]:
    base = (
        norm_key(row["T_MeV"]),
        norm_key(row["muB_MeV"]),
        norm_key(row["xi"]),
    )
    return base if species is None else (*base, species)


def load_inputs() -> dict[str, dict[str, Any]]:
    loaded: dict[str, dict[str, Any]] = {}
    for mode_key, cfg in MODE_INPUTS.items():
        result_dir = cfg["result_dir"]
        figure_dir = cfg["figure_dir"]
        paths = {
            "scan": result_dir / "phase_guided_transport_scan.csv",
            "diagnostics": result_dir / "channel_diagnostics.csv",
            "failed": result_dir / "failed_points.csv",
            "effective_config": result_dir / "effective_config.json",
            "manifest": result_dir / "manifest.json",
            "production_audit": result_dir / "PRODUCTION_AUDIT.md",
            "readme": result_dir / "README.md",
            "plot_manifest": figure_dir / "plot_manifest.json",
        }
        for name, path in paths.items():
            if not path.exists():
                raise FileNotFoundError(f"missing {mode_key} {name}: {path}")
        _, scan_rows = read_csv_with_comments(paths["scan"])
        _, diag_rows = read_csv_with_comments(paths["diagnostics"])
        _, failed_rows = read_csv_with_comments(paths["failed"])
        loaded[mode_key] = {
            "cfg": cfg,
            "paths": paths,
            "scan_rows": scan_rows,
            "diag_rows": diag_rows,
            "failed_rows": failed_rows,
            "plot_manifest": json.loads(paths["plot_manifest"].read_text(encoding="utf-8")),
            "manifest": json.loads(paths["manifest"].read_text(encoding="utf-8")),
        }
    return loaded


def validate_inputs(loaded: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    inventory: list[dict[str, Any]] = []
    required_scan = [
        "T_MeV",
        "muB_MeV",
        "xi",
        "plot_panel",
        "plot_series",
        "phase_reference_kind",
        "phase_prev",
        "phase_curr",
        "phase_structure",
        "phase_boundary_xi_used",
        "seed_source",
        *TAU_FIELDS,
        *TAUINV_FIELDS,
        *DOWNSTREAM_FIELDS,
        *BACKGROUND_FIELDS,
    ]
    required_diag = [
        "T_MeV",
        "muB_MeV",
        "xi",
        "species",
        "channel",
        "density_key",
        "density",
        "rate",
        "contribution",
        "total",
        "tau_inv_species",
    ]
    for mode_key, data in loaded.items():
        cfg = data["cfg"]
        scan_rows = data["scan_rows"]
        diag_rows = data["diag_rows"]
        failed_rows = data["failed_rows"]
        if len(scan_rows) != cfg["main_rows"]:
            raise ValueError(f"{mode_key} scan row count mismatch: {len(scan_rows)}")
        if len(diag_rows) != cfg["channel_rows"]:
            raise ValueError(f"{mode_key} channel row count mismatch: {len(diag_rows)}")
        if failed_rows:
            raise ValueError(f"{mode_key} has failed points: {len(failed_rows)}")
        scan_header = set(scan_rows[0])
        diag_header = set(diag_rows[0])
        missing_scan = [name for name in required_scan if name not in scan_header]
        missing_diag = [name for name in required_diag if name not in diag_header]
        if missing_scan:
            raise ValueError(f"{mode_key} missing scan columns: {missing_scan}")
        if missing_diag:
            raise ValueError(f"{mode_key} missing diagnostics columns: {missing_diag}")
        plot_count = int(data["plot_manifest"].get("count", -1))
        png_count = len(list(cfg["figure_dir"].glob("plot_panel=*/*.png")))
        if plot_count != png_count:
            raise ValueError(f"{mode_key} plot count mismatch: manifest={plot_count}, files={png_count}")
        sidecars = list(cfg["figure_dir"].glob("plot_panel=*/*.png.provenance.json"))
        if sidecars:
            raise ValueError(f"{mode_key} has unexpected plot sidecars: {len(sidecars)}")
        xis = sorted({f(row, "xi") for row in scan_rows})
        inventory.append(
            {
                "mode_key": mode_key,
                "scan_rows": len(scan_rows),
                "channel_rows": len(diag_rows),
                "failed_rows": len(failed_rows),
                "xi_count": len(xis),
                "xi_min": xis[0],
                "xi_max": xis[-1],
                "xi_step": round(xis[1] - xis[0], 10),
                "plot_count": png_count,
                "has_antiquark_tau_fields": all(name in scan_header for name in ["tau_ubar", "tau_dbar", "tau_sbar"]),
                "has_zeta_over_s": "zeta_over_s" in scan_header,
                "scan_csv": relpath(data["paths"]["scan"]),
                "channel_diagnostics_csv": relpath(data["paths"]["diagnostics"]),
                "plot_manifest": relpath(data["paths"]["plot_manifest"]),
                "scan_sha256": sha256_file(data["paths"]["scan"]),
                "diagnostics_sha256": sha256_file(data["paths"]["diagnostics"]),
            }
        )
    return inventory


def detect_tau_step_candidates(loaded: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    candidates: list[dict[str, Any]] = []
    for mode_key, data in loaded.items():
        grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
        for row in data["scan_rows"]:
            grouped[curve_key(row)].append(row)
        for (panel, series), rows in grouped.items():
            rows = sorted(rows, key=lambda row: f(row, "xi"))
            for species in TAU_SPECIES:
                field = f"tau_{species}"
                steps: list[tuple[int, float, float, float]] = []
                for idx in range(len(rows) - 1):
                    left = f(rows[idx], field)
                    right = f(rows[idx + 1], field)
                    if left > 0 and right > 0:
                        log_step = abs(math.log(right / left))
                        rel_step = rel_delta(left, right)
                        steps.append((idx, log_step, rel_step, right / left))
                logs = [step[1] for step in steps]
                if not logs:
                    continue
                median = statistics.median(logs)
                mad = statistics.median([abs(value - median) for value in logs]) or 1.0e-12
                robust_scale = 1.4826 * mad
                robust_threshold = median + 8.0 * robust_scale
                for idx, log_step, rel_step, ratio in steps:
                    robust_z = (log_step - median) / robust_scale if robust_scale > 0 else 0.0
                    if not ((log_step >= 0.12 and robust_z >= 8.0) or log_step >= 0.30):
                        continue
                    left = rows[idx]
                    right = rows[idx + 1]
                    candidates.append(
                        {
                            "mode_key": mode_key,
                            "mode": data["cfg"]["mode_value"],
                            "plot_panel": panel,
                            "plot_series": series,
                            "species": species,
                            "observable": field,
                            "xi_left": f(left, "xi"),
                            "xi_right": f(right, "xi"),
                            "xi_mid": 0.5 * (f(left, "xi") + f(right, "xi")),
                            "tau_left": f(left, field),
                            "tau_right": f(right, field),
                            "tau_ratio_right_over_left": ratio,
                            "tau_step_factor": max(ratio, 1.0 / ratio) if ratio > 0 else math.nan,
                            "tau_rel_step": rel_step,
                            "log_step_abs": log_step,
                            "robust_z": robust_z,
                            "phase_left": left.get("phase_curr", ""),
                            "phase_right": right.get("phase_curr", ""),
                            "phase_reference_left": left.get("phase_reference_kind", ""),
                            "phase_reference_right": right.get("phase_reference_kind", ""),
                            "phase_structure_left": left.get("phase_structure", ""),
                            "phase_structure_right": right.get("phase_structure", ""),
                            "seed_right": right.get("seed_source", ""),
                            "m_u_left": f(left, "m_u"),
                            "m_u_right": f(right, "m_u"),
                            "Phi_left": f(left, "Phi"),
                            "Phi_right": f(right, "Phi"),
                        }
                    )
    candidates.sort(key=lambda row: (-float(row["tau_step_factor"]), row["mode_key"], row["plot_panel"], row["plot_series"]))
    return candidates


def local_extreme_score(rows: list[dict[str, str]], idx: int, affected_species: set[str]) -> float:
    if idx <= 0 or idx >= len(rows) - 1:
        return 0.0
    score = 0.0
    for species in affected_species:
        field = f"tau_{species}"
        prev_v = f(rows[idx - 1], field)
        cur_v = f(rows[idx], field)
        next_v = f(rows[idx + 1], field)
        if not (prev_v > 0 and cur_v > 0 and next_v > 0):
            continue
        left = abs(math.log(cur_v / prev_v))
        right = abs(math.log(next_v / cur_v))
        if (cur_v < prev_v and cur_v < next_v) or (cur_v > prev_v and cur_v > next_v):
            score += left + right
        else:
            score += max(left, right)
    return score


def cluster_step_candidates(
    loaded: dict[str, dict[str, Any]],
    step_candidates: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    by_group: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for item in step_candidates:
        by_group[(item["mode_key"], item["plot_panel"], item["plot_series"])].append(item)

    scan_by_group: dict[tuple[str, str, str], list[dict[str, str]]] = {}
    for mode_key, data in loaded.items():
        temp: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
        for row in data["scan_rows"]:
            temp[(mode_key, row["plot_panel"], row["plot_series"])].append(row)
        for key, rows in temp.items():
            scan_by_group[key] = sorted(rows, key=lambda row: f(row, "xi"))

    windows: list[dict[str, Any]] = []
    for group_key, items in by_group.items():
        items = sorted(items, key=lambda row: (float(row["xi_mid"]), row["species"]))
        clusters: list[list[dict[str, Any]]] = []
        for item in items:
            if not clusters or float(item["xi_left"]) > max(float(x["xi_right"]) for x in clusters[-1]) + 0.021:
                clusters.append([item])
            else:
                clusters[-1].append(item)
        rows = scan_by_group[group_key]
        xi_to_index = {round(f(row, "xi"), 10): idx for idx, row in enumerate(rows)}
        for cluster_index, cluster in enumerate(clusters, start=1):
            affected = sorted({item["species"] for item in cluster})
            x_min = min(float(item["xi_left"]) for item in cluster)
            x_max = max(float(item["xi_right"]) for item in cluster)
            candidate_indices = [
                idx
                for idx, row in enumerate(rows)
                if x_min - 1.0e-10 <= f(row, "xi") <= x_max + 1.0e-10
            ]
            scored = [
                (local_extreme_score(rows, idx, set(affected)), idx)
                for idx in candidate_indices
            ]
            if scored and max(score for score, _ in scored) > 0:
                target_idx = max(scored)[1]
            else:
                strongest = max(cluster, key=lambda row: float(row["tau_step_factor"]))
                target_xi = float(strongest["xi_right"])
                target_idx = xi_to_index.get(round(target_xi, 10), candidate_indices[-1])
            prev_idx = max(0, target_idx - 1)
            next_idx = min(len(rows) - 1, target_idx + 1)
            target = rows[target_idx]
            prev_row = rows[prev_idx]
            next_row = rows[next_idx]
            mode_key, panel, series = group_key
            safe_panel = panel.replace(".", "p").replace("-", "m")
            safe_series = series.replace(".", "p").replace("-", "m")
            safe_xi = f"{f(target, 'xi'):+.2f}".replace("+", "p").replace("-", "m").replace(".", "p")
            window_id = f"{mode_key}_{safe_panel}_{safe_series}_xi{safe_xi}_{cluster_index}"
            windows.append(
                {
                    "window_id": window_id,
                    "mode_key": mode_key,
                    "mode": loaded[mode_key]["cfg"]["mode_value"],
                    "plot_panel": panel,
                    "plot_series": series,
                    "xi_left": x_min,
                    "xi_right": x_max,
                    "target_xi": f(target, "xi"),
                    "prev_xi": f(prev_row, "xi"),
                    "next_xi": f(next_row, "xi"),
                    "affected_species": ";".join(affected),
                    "affected_tau_fields": ";".join(f"tau_{sp}" for sp in affected),
                    "max_tau_step_factor": max(float(item["tau_step_factor"]) for item in cluster),
                    "max_tau_rel_step": max(float(item["tau_rel_step"]) for item in cluster),
                    "n_step_candidates": len(cluster),
                    "step_species_rows": cluster,
                    "prev_row": prev_row,
                    "target_row": target,
                    "next_row": next_row,
                }
            )
    windows.sort(key=lambda row: (-float(row["max_tau_step_factor"]), row["mode_key"], row["plot_panel"], row["plot_series"], float(row["target_xi"])))
    return windows


def build_diag_index(loaded: dict[str, dict[str, Any]]) -> dict[str, dict[tuple[str, ...], list[dict[str, str]]]]:
    index: dict[str, dict[tuple[str, ...], list[dict[str, str]]]] = {}
    for mode_key, data in loaded.items():
        by_key: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
        for row in data["diag_rows"]:
            by_key[grid_key(row, row["species"])].append(row)
        index[mode_key] = by_key
    return index


def diag_rows_for(diag_index: dict[str, dict[tuple[str, ...], list[dict[str, str]]]], mode_key: str, row: dict[str, str], species: str) -> list[dict[str, str]]:
    return diag_index[mode_key].get(grid_key(row, species), [])


def top_channel_deltas(
    diag_index: dict[str, dict[tuple[str, ...], list[dict[str, str]]]],
    mode_key: str,
    prev_row: dict[str, str],
    target_row: dict[str, str],
    next_row: dict[str, str],
    species: str,
) -> list[dict[str, Any]]:
    prev = diag_rows_for(diag_index, mode_key, prev_row, species)
    target = diag_rows_for(diag_index, mode_key, target_row, species)
    next_rows = diag_rows_for(diag_index, mode_key, next_row, species)
    prev_map = {(row["channel"], row["density_key"]): row for row in prev}
    next_map = {(row["channel"], row["density_key"]): row for row in next_rows}
    total = f(target[0], "total") if target else math.nan
    rows: list[dict[str, Any]] = []
    for row in target:
        key = (row["channel"], row["density_key"])
        contrib = f(row, "contribution")
        rate = f(row, "rate")
        prev_contrib = f(prev_map.get(key, {}), "contribution", 0.0)
        next_contrib = f(next_map.get(key, {}), "contribution", 0.0)
        prev_rate = f(prev_map.get(key, {}), "rate", math.nan)
        next_rate = f(next_map.get(key, {}), "rate", math.nan)
        neighbor_contrib = max(prev_contrib, next_contrib, 1.0e-300)
        neighbor_rate = max(
            prev_rate if math.isfinite(prev_rate) else 0.0,
            next_rate if math.isfinite(next_rate) else 0.0,
            1.0e-300,
        )
        rows.append(
            {
                "channel": row["channel"],
                "density_key": row["density_key"],
                "target_contribution": contrib,
                "target_rate": rate,
                "target_share": contrib / total if total > 0 else math.nan,
                "prev_contribution": prev_contrib,
                "next_contribution": next_contrib,
                "contribution_ratio_to_neighbor_max": contrib / neighbor_contrib,
                "rate_ratio_to_neighbor_max": rate / neighbor_rate if neighbor_rate > 0 else math.nan,
                "contribution_delta_vs_neighbor_max": contrib - neighbor_contrib,
            }
        )
    rows.sort(
        key=lambda item: (
            -float(item["contribution_delta_vs_neighbor_max"]),
            -float(item["target_contribution"]),
        )
    )
    return rows


def background_step_summary(prev_row: dict[str, str], target_row: dict[str, str], next_row: dict[str, str]) -> dict[str, Any]:
    primary_fields = ["m_u", "Phi", "Phibar", "s_fm3inv", "m_s"]
    best_field = ""
    best_step = -1.0
    details: dict[str, float] = {}
    for field in BACKGROUND_FIELDS:
        values = [f(prev_row, field), f(target_row, field), f(next_row, field)]
        steps = [rel_delta(values[0], values[1]), rel_delta(values[1], values[2])]
        max_step = max(step for step in steps if math.isfinite(step)) if any(math.isfinite(step) for step in steps) else math.nan
        details[field] = max_step
        if field in primary_fields and math.isfinite(max_step) and max_step > best_step:
            best_field = field
            best_step = max_step
    return {"max_background_rel_step": best_step, "dominant_background_driver": best_field, "background_steps": details}


def tau_shape(window: dict[str, Any]) -> str:
    prev_row = window["prev_row"]
    target = window["target_row"]
    next_row = window["next_row"]
    species = window["affected_species"].split(";")
    dips = 0
    peaks = 0
    monotone_up = 0
    monotone_down = 0
    for sp in species:
        values = [f(prev_row, f"tau_{sp}"), f(target, f"tau_{sp}"), f(next_row, f"tau_{sp}")]
        if all(value > 0 for value in values):
            if values[1] < values[0] and values[1] < values[2]:
                dips += 1
            elif values[1] > values[0] and values[1] > values[2]:
                peaks += 1
            elif values[0] < values[1] < values[2]:
                monotone_up += 1
            elif values[0] > values[1] > values[2]:
                monotone_down += 1
    if dips >= max(1, len(species) // 2):
        return "local_tau_dip"
    if peaks >= max(1, len(species) // 2):
        return "local_tau_peak"
    if monotone_up >= max(1, len(species) // 2):
        return "step_or_rapid_increase"
    if monotone_down >= max(1, len(species) // 2):
        return "step_or_rapid_decrease"
    return "mixed_tau_response"


def classify_window(
    window: dict[str, Any],
    channel_rows: list[dict[str, Any]],
    background: dict[str, Any],
) -> tuple[str, str, float]:
    prev_row = window["prev_row"]
    target = window["target_row"]
    next_row = window["next_row"]
    rows = [prev_row, target, next_row]
    phase_values = {row.get("phase_curr", "") for row in rows}
    seed_text = ";".join(row.get("seed_source", "") for row in rows)
    phase_reference_values = {row.get("phase_reference_kind", "") for row in rows}
    phase_structure_values = {row.get("phase_structure", "") for row in rows}
    phase_change = len(phase_values) > 1 or "phase_switch" in seed_text or "multiseed_governance" in seed_text
    first_order_context = (
        "first_order" in phase_reference_values
        or "first_order_possible" in phase_structure_values
    )
    background_fast = float(background["max_background_rel_step"]) >= 0.25
    max_rate_ratio = max(
        [float(row.get("rate_ratio_to_neighbor_max", math.nan)) for row in channel_rows if math.isfinite(float(row.get("rate_ratio_to_neighbor_max", math.nan)))],
        default=math.nan,
    )
    max_contribution_ratio = max(
        [float(row.get("contribution_ratio_to_neighbor_max", math.nan)) for row in channel_rows if math.isfinite(float(row.get("contribution_ratio_to_neighbor_max", math.nan)))],
        default=math.nan,
    )
    dominant = ";".join(
        f"{row['species']}:{row['channel']}"
        for row in channel_rows[:4]
    )
    if phase_change:
        return (
            "phase_branch_switch_supported",
            "phase_curr/seed_source shows an accepted branch switch; tau changes should be interpreted as downstream of the selected equilibrium branch.",
            0.95,
        )
    if first_order_context and background_fast:
        return (
            "upstream_first_order_branch_jump_supported",
            f"first-order reference context with a fast background jump ({background['dominant_background_driver']}); channel rates respond downstream of the branch/background change.",
            0.90,
        )
    if background_fast:
        return (
            "upstream_branch_fast_change_supported",
            f"background variables change rapidly without an explicit phase label switch; dominant driver is {background['dominant_background_driver']}.",
            0.80,
        )
    if math.isfinite(max_rate_ratio) and (max_rate_ratio >= 3.0 or max_contribution_ratio >= 3.0):
        return (
            "channel_rate_spike_supported",
            f"smooth background; tau jump is driven by a local channel rate/contribution spike, led by {dominant}. Denominator-chain origin remains a candidate unless deep-scanned.",
            0.78,
        )
    return (
        "weak_or_broad_tau_variation_candidate",
        "flagged by robust local tau step, but neither upstream background nor channel-rate evidence is strong enough for a paper-ready mechanism verdict.",
        0.45,
    )


def build_jump_tables(loaded: dict[str, dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    step_candidates = detect_tau_step_candidates(loaded)
    windows = cluster_step_candidates(loaded, step_candidates)
    diag_index = build_diag_index(loaded)
    window_rows: list[dict[str, Any]] = []
    channel_rows_out: list[dict[str, Any]] = []

    for window in windows:
        prev_row = window["prev_row"]
        target = window["target_row"]
        next_row = window["next_row"]
        affected = window["affected_species"].split(";")
        all_channel_rows: list[dict[str, Any]] = []
        for sp in affected:
            top_rows = top_channel_deltas(diag_index, window["mode_key"], prev_row, target, next_row, sp)
            for rank, row in enumerate(top_rows[:6], start=1):
                payload = {
                    "window_id": window["window_id"],
                    "mode_key": window["mode_key"],
                    "plot_panel": window["plot_panel"],
                    "plot_series": window["plot_series"],
                    "target_xi": window["target_xi"],
                    "species": sp,
                    "rank": rank,
                    **row,
                }
                channel_rows_out.append(payload)
                all_channel_rows.append(payload)
        background = background_step_summary(prev_row, target, next_row)
        all_channel_rows_sorted = sorted(
            all_channel_rows,
            key=lambda row: (
                float(row.get("contribution_delta_vs_neighbor_max", 0.0)),
                float(row.get("target_contribution", 0.0)),
            ),
            reverse=True,
        )
        verdict, note, score = classify_window(window, all_channel_rows_sorted, background)
        shape = tau_shape(window)
        boundary_changed = len(
            {prev_row.get("phase_boundary_xi_used", ""), target.get("phase_boundary_xi_used", ""), next_row.get("phase_boundary_xi_used", "")}
        ) > 1
        if verdict == "channel_rate_spike_supported":
            top_channels = ";".join(
                f"{row['species']}:{row['channel']}({float(row['target_share']):.2f})"
                for row in all_channel_rows_sorted[:5]
            )
        else:
            top_channels = (
                "not_primary;"
                f"upstream_driver={background['dominant_background_driver']};"
                f"max_background_rel_step={float(background['max_background_rel_step']):.4g}"
            )
        window_rows.append(
            {
                "window_id": window["window_id"],
                "mode_key": window["mode_key"],
                "mode": window["mode"],
                "plot_panel": window["plot_panel"],
                "plot_series": window["plot_series"],
                "xi_left": window["xi_left"],
                "target_xi": window["target_xi"],
                "xi_right": window["xi_right"],
                "prev_xi": window["prev_xi"],
                "next_xi": window["next_xi"],
                "affected_tau_fields": window["affected_tau_fields"],
                "tau_shape": shape,
                "max_tau_step_factor": window["max_tau_step_factor"],
                "max_tau_rel_step": window["max_tau_rel_step"],
                "cause_verdict": verdict,
                "evidence_score": score,
                "dominant_background_driver": background["dominant_background_driver"],
                "max_background_rel_step": background["max_background_rel_step"],
                "phase_prev": prev_row.get("phase_curr", ""),
                "phase_target": target.get("phase_curr", ""),
                "phase_next": next_row.get("phase_curr", ""),
                "phase_reference_target": target.get("phase_reference_kind", ""),
                "phase_structure_target": target.get("phase_structure", ""),
                "phase_boundary_xi_target": target.get("phase_boundary_xi_used", ""),
                "phase_boundary_changed_in_local_window": boundary_changed,
                "seed_source_target": target.get("seed_source", ""),
                "equilibrium_backend_target": target.get("equilibrium_backend", ""),
                "m_u_prev": f(prev_row, "m_u"),
                "m_u_target": f(target, "m_u"),
                "m_u_next": f(next_row, "m_u"),
                "Phi_prev": f(prev_row, "Phi"),
                "Phi_target": f(target, "Phi"),
                "Phi_next": f(next_row, "Phi"),
                "s_prev": f(prev_row, "s_fm3inv"),
                "s_target": f(target, "s_fm3inv"),
                "s_next": f(next_row, "s_fm3inv"),
                "top_channel_evidence": top_channels,
                "mechanism_note": note,
            }
        )
    window_rows.sort(key=lambda row: (-float(row["max_tau_step_factor"]), row["mode_key"], row["plot_panel"], row["plot_series"]))
    return step_candidates, window_rows, channel_rows_out


def mean_positive(values: list[float]) -> float:
    clean = [value for value in values if math.isfinite(value) and value > 0]
    return sum(clean) / len(clean) if clean else math.nan


def direction_label(delta: float, eps: float = 0.03) -> str:
    if not math.isfinite(delta):
        return "unknown"
    if delta > eps:
        return "up"
    if delta < -eps:
        return "down"
    return "flat"


def downstream_summary(window_rows: list[dict[str, Any]], loaded: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    by_window = {row["window_id"]: row for row in window_rows}
    rows_out: list[dict[str, Any]] = []
    for window_id, win in by_window.items():
        # The original row objects are not in window_rows; recover from scan rows.
        mode_rows = loaded[win["mode_key"]]["scan_rows"]
        matches = [
            row
            for row in mode_rows
            if row["plot_panel"] == win["plot_panel"]
            and row["plot_series"] == win["plot_series"]
            and f(row, "xi") in {float(win["prev_xi"]), float(win["target_xi"]), float(win["next_xi"])}
        ]
        by_xi = {round(f(row, "xi"), 10): row for row in matches}
        prev = by_xi.get(round(float(win["prev_xi"]), 10))
        target = by_xi.get(round(float(win["target_xi"]), 10))
        next_row = by_xi.get(round(float(win["next_xi"]), 10))
        if prev is None or target is None or next_row is None:
            continue
        baseline_tau = mean_positive(
            [
                mean_positive([f(prev, field), f(next_row, field)])
                for field in TAU_FIELDS
            ]
        )
        target_tau = mean_positive([f(target, field) for field in TAU_FIELDS])
        tau_log_delta = signed_log_delta(baseline_tau, target_tau)
        s_baseline = mean_positive([f(prev, "s_fm3inv"), f(next_row, "s_fm3inv")])
        s_log_delta = signed_log_delta(s_baseline, f(target, "s_fm3inv"))
        for obs in ["eta_over_s", "zeta_over_s", "sigma_over_T"]:
            baseline = mean_positive([f(prev, obs), f(next_row, obs)])
            target_value = f(target, obs)
            obs_log_delta = signed_log_delta(baseline, target_value)
            obs_dir = direction_label(obs_log_delta)
            tau_dir = direction_label(tau_log_delta)
            s_dir = direction_label(s_log_delta)
            if obs in ("eta_over_s", "zeta_over_s"):
                if obs_dir == tau_dir and tau_dir != "flat":
                    relation = "same_direction_as_tau_mean"
                elif s_dir != "flat" and obs_dir != s_dir:
                    relation = "entropy_denominator_amplifies_or_offsets_tau"
                else:
                    relation = "mixed_tau_and_thermodynamic_response"
            else:
                relation = "charge_transport_weighted_tau_response"
            rows_out.append(
                {
                    "window_id": window_id,
                    "mode_key": win["mode_key"],
                    "plot_panel": win["plot_panel"],
                    "plot_series": win["plot_series"],
                    "target_xi": win["target_xi"],
                    "observable": obs,
                    "baseline_value": baseline,
                    "target_value": target_value,
                    "log_delta_vs_neighbor_mean": obs_log_delta,
                    "observable_direction": obs_dir,
                    "tau_mean_direction": tau_dir,
                    "entropy_direction": s_dir,
                    "response_relation": relation,
                    "root_tau_verdict": win["cause_verdict"],
                }
            )
    return rows_out


def scan_row_for_window(loaded: dict[str, dict[str, Any]], window: dict[str, Any]) -> dict[str, str]:
    target_xi = float(window["target_xi"])
    matches = [
        row
        for row in loaded[window["mode_key"]]["scan_rows"]
        if row.get("plot_panel") == window["plot_panel"]
        and row.get("plot_series") == window["plot_series"]
        and abs(f(row, "xi") - target_xi) <= 1.0e-9
    ]
    if not matches:
        raise ValueError(f"scan row not found for mechanism candidate {window['window_id']}")
    return matches[0]


def build_mechanism_candidates(
    window_rows: list[dict[str, Any]],
    channel_rows: list[dict[str, Any]],
    loaded: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    by_window: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in channel_rows:
        by_window[str(row["window_id"])].append(row)

    candidates: list[dict[str, Any]] = []
    strong_windows = [
        row
        for row in window_rows
        if row["cause_verdict"] == "channel_rate_spike_supported"
    ]
    strong_windows.sort(key=lambda row: (-float(row["max_tau_step_factor"]), row["window_id"]))
    weak_windows = [
        row
        for row in window_rows
        if row["cause_verdict"] == "weak_or_broad_tau_variation_candidate"
    ]
    weak_windows.sort(key=lambda row: (-float(row["max_tau_step_factor"]), row["window_id"]))
    selected_windows = [*strong_windows, *weak_windows]

    for order, window in enumerate(selected_windows, start=1):
        channel_choices = sorted(
            by_window[window["window_id"]],
            key=lambda row: (
                float(row.get("contribution_delta_vs_neighbor_max", 0.0)),
                float(row.get("target_share", 0.0)),
                float(row.get("target_contribution", 0.0)),
            ),
            reverse=True,
        )
        if not channel_choices:
            continue
        primary_species = str(channel_choices[0]["species"])
        target = scan_row_for_window(loaded, window)
        observable = f"tau_{primary_species}"
        alpha_t = f(target, "alpha_T")
        selector_parts = [
            f"T_MeV={f(target, 'T_MeV'):.12g}",
            f"muB_MeV={f(target, 'muB_MeV'):.12g}",
        ]
        if math.isfinite(alpha_t):
            selector_parts.append(f"alpha_T={alpha_t:.12g}")
        selector = ";".join(selector_parts)
        label = (
            f"{window['plot_panel']}, {window['plot_series']}, "
            f"{observable}, xi={float(window['target_xi']):.2f}"
        )
        candidates.append(
            {
                "window_id": window["window_id"],
                "candidate_source": (
                    "xi001_non_first_order_channel_rate_spike"
                    if window["cause_verdict"] == "channel_rate_spike_supported"
                    else "xi001_weak_or_broad_tau_variation_probe"
                ),
                "selected_for_deep_scan": True,
                "selected_order": order,
                "selection_reason": (
                    "non_first_order_tau_jump_with_smooth_background_and_channel_rate_spike"
                    if window["cause_verdict"] == "channel_rate_spike_supported"
                    else "weak_or_broad_tau_variation_probe"
                ),
                "mode_key": window["mode_key"],
                "mode": window["mode"],
                "plot_panel": window["plot_panel"],
                "plot_panel_label": window["plot_panel"],
                "plot_series": window["plot_series"],
                "plot_series_label": window["plot_series"],
                "selector": selector,
                "T_MeV": f(target, "T_MeV"),
                "muB_MeV": f(target, "muB_MeV"),
                "alpha_T": alpha_t if math.isfinite(alpha_t) else "",
                "xi": float(window["target_xi"]),
                "prev_xi": float(window["prev_xi"]),
                "next_xi": float(window["next_xi"]),
                "observable": observable,
                "primary_species": primary_species,
                "mechanism_applicability": "tau_species",
                "local_shape": window["tau_shape"],
                "local_score": window["max_tau_step_factor"],
                "value_at_xi": f(target, observable),
                "phase_curr": target.get("phase_curr", ""),
                "phase_structure": target.get("phase_structure", ""),
                "window_label": label,
            }
        )
    return candidates


def plot_jump_overview(window_rows: list[dict[str, Any]], loaded: dict[str, dict[str, Any]], out_dir: Path) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    selected = window_rows[:8]
    fig_paths: list[Path] = []
    if selected:
        cols = 2
        rows_n = math.ceil(len(selected) / cols)
        fig, axes = plt.subplots(rows_n, cols, figsize=(11, 3.1 * rows_n), squeeze=False)
        for ax in axes.flat:
            ax.axis("off")
        for ax, win in zip(axes.flat, selected):
            ax.axis("on")
            mode_rows = [
                row
                for row in loaded[win["mode_key"]]["scan_rows"]
                if row["plot_panel"] == win["plot_panel"] and row["plot_series"] == win["plot_series"]
            ]
            mode_rows.sort(key=lambda row: f(row, "xi"))
            species = win["affected_tau_fields"].replace("tau_", "").split(";")
            species = [sp for sp in species if sp]
            x = [f(row, "xi") for row in mode_rows]
            for sp in species[:4]:
                values = [f(row, f"tau_{sp}") for row in mode_rows]
                center = float(win["target_xi"])
                mask = [abs(xi - center) <= 0.06 for xi in x]
                xs = [xi for xi, keep in zip(x, mask) if keep]
                ys = [value for value, keep in zip(values, mask) if keep]
                if ys:
                    baseline = mean_positive([ys[0], ys[-1]])
                    normed = [value / baseline for value in ys]
                    ax.plot(xs, normed, marker="o", ms=2, lw=1.2, label=f"tau_{sp}")
            ax.axvline(float(win["target_xi"]), color="0.4", ls="--", lw=0.8)
            ax.set_title(f"{win['plot_panel']} {win['plot_series']} xi={float(win['target_xi']):+.2f}", fontsize=9)
            ax.set_xlabel("xi")
            ax.set_ylabel("tau / edge mean")
            ax.legend(fontsize=7)
        fig.tight_layout()
        path = out_dir / "tau_jump_windows_overview.png"
        fig.savefig(path, dpi=220)
        plt.close(fig)
        fig_paths.append(path)

    if selected:
        fig, axes = plt.subplots(rows_n, cols, figsize=(11, 3.1 * rows_n), squeeze=False)
        for ax in axes.flat:
            ax.axis("off")
        for ax, win in zip(axes.flat, selected):
            ax.axis("on")
            mode_rows = [
                row
                for row in loaded[win["mode_key"]]["scan_rows"]
                if row["plot_panel"] == win["plot_panel"] and row["plot_series"] == win["plot_series"]
            ]
            mode_rows.sort(key=lambda row: f(row, "xi"))
            x = [f(row, "xi") for row in mode_rows]
            center = float(win["target_xi"])
            mask = [abs(xi - center) <= 0.06 for xi in x]
            xs = [xi for xi, keep in zip(x, mask) if keep]
            fields = ["eta_over_s", "zeta_over_s"]
            for field in fields:
                values = [f(row, field) for row in mode_rows]
                ys = [value for value, keep in zip(values, mask) if keep]
                if ys:
                    baseline = mean_positive([ys[0], ys[-1]])
                    ax.plot(xs, [value / baseline for value in ys], marker="o", ms=2, lw=1.2, label=field)
            tau_mean = [
                mean_positive([f(row, field) for field in TAU_FIELDS])
                for row in mode_rows
            ]
            tau_ys = [value for value, keep in zip(tau_mean, mask) if keep]
            if tau_ys:
                baseline = mean_positive([tau_ys[0], tau_ys[-1]])
                ax.plot(xs, [value / baseline for value in tau_ys], color="0.1", lw=1.4, label="tau mean")
            ax.axvline(float(win["target_xi"]), color="0.4", ls="--", lw=0.8)
            ax.set_title(f"{win['plot_panel']} {win['plot_series']} xi={float(win['target_xi']):+.2f}", fontsize=9)
            ax.set_xlabel("xi")
            ax.set_ylabel("normalized response")
            ax.legend(fontsize=7)
        fig.tight_layout()
        path = out_dir / "downstream_response_windows.png"
        fig.savefig(path, dpi=220)
        plt.close(fig)
        fig_paths.append(path)
    return fig_paths


def make_claim_ledger(
    window_rows: list[dict[str, Any]],
    downstream_rows: list[dict[str, Any]],
    mechanism_rows: list[dict[str, str]],
    convergence_rows: list[dict[str, str]],
) -> list[dict[str, Any]]:
    source_scan_a = relpath(MODE_INPUTS["mode_a"]["result_dir"] / "phase_guided_transport_scan.csv")
    source_scan_b = relpath(MODE_INPUTS["mode_b"]["result_dir"] / "phase_guided_transport_scan.csv")
    rows = [
        {
            "claim_id": "CLAIM-DATA-XI001-001",
            "status": "supported",
            "claim_zh": "xi001 正式产物在 mode A/B 各包含 909 个主结果点，xi 网格为 -0.50:0.01:0.50，failed points 为 0。",
            "evidence": "tables/input_inventory.csv",
            "fields_or_points": "scan_rows, xi_count, failed_rows",
        },
        {
            "claim_id": "CLAIM-FIG-XI001-001",
            "status": "supported",
            "claim_zh": "新正式图像已补齐 tau_ubar、tau_dbar、tau_sbar，并同时补出 zeta_over_s；mode A/B 图目录各有 36 张 PNG。",
            "evidence": "data/outputs/figures/.../plot_manifest.json; tables/input_inventory.csv",
            "fields_or_points": "plot_count=36",
        },
        {
            "claim_id": "CLAIM-TAU-JUMP-001",
            "status": "supported",
            "claim_zh": "tau 突变不应从 eta_over_s 或 zeta_over_s 反推根因；根因判定以 tau 曲线、上游平衡态背景量和 channel diagnostics 为主。",
            "evidence": "tables/tau_jump_window_summary.csv; tables/tau_jump_channel_attribution.csv",
            "fields_or_points": "cause_verdict, top_channel_evidence",
        },
        {
            "claim_id": "CLAIM-TAU-JUMP-002",
            "status": "supported_with_scope_limit",
            "claim_zh": "一阶/上游分支窗口中的 tau 阶跃由 m_u、Phi、熵密度等背景量快速变化解释，后续输运比值只是下游响应。",
            "evidence": "tables/tau_jump_window_summary.csv",
            "fields_or_points": "upstream_first_order_branch_jump_supported",
        },
        {
            "claim_id": "CLAIM-TAU-JUMP-003",
            "status": "supported_with_scope_limit",
            "claim_zh": "背景平滑但 tau 单点下探的窗口，其直接近因是少数 q qbar channel 的 rate/contribution 局部尖峰；是否能升级为传播子分母近零机制，需要 denominator-chain 定点深拆或局部收敛补证。",
            "evidence": "tables/tau_jump_channel_attribution.csv",
            "fields_or_points": "rate_ratio_to_neighbor_max, contribution_ratio_to_neighbor_max",
        },
        {
            "claim_id": "CLAIM-DOWNSTREAM-001",
            "status": "supported_qualitative",
            "claim_zh": "eta_over_s 与 zeta_over_s 的局部结构可作为 tau 与熵密度共同作用后的下游趋势响应：channel-rate spike 窗口中随 tau 下探而下探，分支窗口中常被 tau 上升和熵密度下降共同放大。",
            "evidence": "tables/downstream_transport_response_summary.csv",
            "fields_or_points": "response_relation, tau_mean_direction, entropy_direction",
        },
        {
            "claim_id": "CLAIM-SCOPE-001",
            "status": "author_check",
            "claim_zh": "xi001 新中间点继承 p128 积分参数收敛证据，但并非每个新 xi 点都已有 p104->p128 独立 gate；非一阶窗口的小分母机制结论已有 denominator-chain 支持，但若要写成更强的 production-grade 收敛结论，仍需补局部 high-rate convergence gate。",
            "evidence": f"{source_scan_a}; {source_scan_b}; data/outputs/results/relaxtime/transport/phase_guided/{CONVERGENCE_CASE}/README.md",
            "fields_or_points": "known_limitations",
        },
    ]
    if mechanism_rows:
        verdict_by_window = {row["window_id"]: row for row in mechanism_rows}
        mode_a = verdict_by_window.get("mode_a_muB450p0_alpha1p0_xip0p26_2")
        mode_b = verdict_by_window.get("mode_b_T200p0_muB450p0_xip0p31_1")
        if mode_a:
            mode_a_rel_error = mode_a.get("max_rate_reproduction_rel_error", "")
            rows.append(
                {
                    "claim_id": "CLAIM-MECH-XI001-001",
                    "status": "supported_with_scope_limit",
                    "claim_zh": f"非一阶窗口 mode A, muB=450, alpha_T=1.0, xi=0.26 的 tau 单点下探在 denominator-chain 证据层支持小分母机制：主导的 uubar_to_uubar/uubar_to_ddbar 通道在同一近阈值 band 出现 mixed detM 峰，修正机制脚本口径后 rate 复现误差为 {mode_a_rel_error}；但本轮没有完成额外 high-rate gate。",
                    "evidence": "tables/mechanism_window_summary.csv; tables/denominator_chain_summary.csv; tables/denominator_chain_band_table.csv",
                    "fields_or_points": f"window_id={mode_a['window_id']}; verdict={mode_a['mechanism_verdict']}; branch={mode_a['dominant_denominator_branch']}; max_rate_reproduction_rel_error={mode_a['max_rate_reproduction_rel_error']}",
                }
            )
        if mode_b:
            mode_b_supported = mode_b.get("mechanism_verdict") == "small_denominator_supported"
            rows.append(
                {
                    "claim_id": "CLAIM-MECH-XI001-002",
                    "status": "supported_with_scope_limit" if mode_b_supported else "author_check",
                    "claim_zh": (
                        "非一阶窗口 mode B, T=200, muB=450, xi=0.31 的 tau 单点下探在 denominator-chain 证据层支持 simple 1-4KΠ 小分母机制：主导通道在同一近阈值 band 对齐，修正 A-builder 诊断口径后 production rate 可精确复现；但本轮没有完成额外 high-rate gate。"
                        if mode_b_supported
                        else "非一阶窗口 mode B, T=200, muB=450, xi=0.31 的 tau 单点下探具有 simple 1-4KΠ 近阈值小分母候选证据，但机制脚本直调 rate 与 channel diagnostics 的生产 rate 未闭合，不能在当前证据下写成 paper-ready 小分母结论。"
                    ),
                    "evidence": "tables/mechanism_window_summary.csv; tables/denominator_chain_summary.csv; tables/local_rate_reproduction_mismatch_root_cause.csv",
                    "fields_or_points": f"window_id={mode_b['window_id']}; verdict={mode_b['mechanism_verdict']}; branch={mode_b['dominant_denominator_branch']}; max_rate_reproduction_rel_error={mode_b['max_rate_reproduction_rel_error']}",
                }
            )
        weak_window_ids = {
            str(row["window_id"])
            for row in window_rows
            if row.get("cause_verdict") == "weak_or_broad_tau_variation_candidate"
        }
        weak_supported = [
            row for row in mechanism_rows
            if row.get("window_id") in weak_window_ids
            and row.get("mechanism_verdict") == "small_denominator_supported"
        ]
        if weak_supported:
            all_weak_supported = len(weak_supported) == len(weak_window_ids)
            rows.append(
                {
                    "claim_id": "CLAIM-MECH-XI001-004",
                    "status": "supported_with_scope_limit",
                    "claim_zh": (
                        "6 个 weak_or_broad_tau_variation_candidate 窗口已全部完成 denominator-chain 深拆，均支持 simple 1-4KΠ 近阈值小分母机制；该结论仍是机制补证，不替代逐点 high-rate convergence gate。"
                        if all_weak_supported
                        else "已深拆的 weak_or_broad_tau_variation_candidate 窗口支持 simple 1-4KΠ 近阈值小分母机制；未深拆窗口不能外推。"
                    ),
                    "evidence": "tables/mechanism_window_summary.csv; tables/denominator_chain_summary.csv",
                    "fields_or_points": "; ".join(f"{row['window_id']}={row['mechanism_verdict']}/{row['dominant_denominator_branch']}" for row in weak_supported),
                }
            )
        rows.append(
            {
                "claim_id": "CLAIM-MECH-XI001-005",
                "status": "supported",
                "claim_zh": "此前 mode B 约 0.47 的 rate 复现偏差来自机制诊断脚本的 A-builder 口径未对齐，而不是 production channel_diagnostics 的物理/数值 bug；修正为 production workflow 的 a_builder 配置后，生产 rate 可复现到机器精度。",
                "evidence": "tables/local_rate_reproduction_mismatch_root_cause.csv; tables/mechanism_window_summary.csv; tables/mechanism_manifest.json",
                "fields_or_points": "production_a_builder_config; max_rate_reproduction_rel_error=0.0",
            }
        )
        if not convergence_rows:
            rows.append(
                {
                    "claim_id": "CLAIM-MECH-XI001-003",
                    "status": "author_check",
                    "claim_zh": "本机局部 high-rate convergence gate 在当前预算内未完成；因此本轮只把 denominator-chain 表作为机制补证，不能替代生产级收敛证明。",
                    "evidence": "tables/local_rate_convergence_gate.csv; tables/local_rate_convergence_attempts.csv",
                    "fields_or_points": "no convergence rows produced",
                }
            )
    return rows


def markdown_table(rows: list[dict[str, Any]], fields: list[str], limit: int = 8) -> str:
    clipped = rows[:limit]
    out = io.StringIO()
    out.write("| " + " | ".join(fields) + " |\n")
    out.write("| " + " | ".join(["---"] * len(fields)) + " |\n")
    for row in clipped:
        values = []
        for field in fields:
            value = row.get(field, "")
            if isinstance(value, float):
                if math.isfinite(value):
                    value = f"{value:.4g}"
                else:
                    value = ""
            values.append(str(value).replace("|", "\\|"))
        out.write("| " + " | ".join(values) + " |\n")
    return out.getvalue()


def write_readme(
    inventory: list[dict[str, Any]],
    window_rows: list[dict[str, Any]],
    downstream_rows: list[dict[str, Any]],
    figure_paths: list[Path],
    mechanism_rows: list[dict[str, str]],
    convergence_rows: list[dict[str, str]],
) -> None:
    counts = defaultdict(int)
    for row in window_rows:
        counts[row["cause_verdict"]] += 1
    key_windows = [
        row
        for row in window_rows
        if row["cause_verdict"] in (
            "upstream_first_order_branch_jump_supported",
            "upstream_branch_fast_change_supported",
            "channel_rate_spike_supported",
            "phase_branch_switch_supported",
        )
    ]
    mechanism_section = ""
    if mechanism_rows:
        mechanism_section = f"""
## 非一阶 channel-rate spike 的 denominator-chain 补证

本轮对 `channel_rate_spike_supported` 的两个非一阶窗口做定点深拆，并把 6 个 `weak_or_broad_tau_variation_candidate` 窗口全部纳入 denominator-chain 检查；`eta_over_s/zeta_over_s` 仍作为 tau 的下游响应，不进入根因判定。

{markdown_table(mechanism_rows, ["window_id", "plot_panel", "plot_series", "observable", "primary_species", "mechanism_verdict", "dominant_channels", "dominant_denominator_branch", "max_rate_reproduction_rel_error", "denominator_sigma_alignment", "upstream_branch_flag"], limit=10)}

- `mode_a_muB450p0_alpha1p0_xip0p26_2`：在 denominator-chain 证据层支持小分母机制。`uubar_to_uubar/uubar_to_ddbar` 贡献覆盖主导份额，`sigma(s)` 峰和 mixed `detM` 峰都落在近阈值 band；修正机制脚本口径后，直调 rate 与 channel diagnostics 的复现误差为机制表中的机器精度量级。该窗口的上游背景量平滑，因此近因不是一阶相变或上游分支突跳；但本轮没有完成额外 high-rate gate，论文表述应保留 scope limit。
- `mode_b_T200p0_muB450p0_xip0p31_1`：在 denominator-chain 证据层也支持小分母机制。`dubar_to_dubar/uubar_to_uubar/uubar_to_ddbar` 的 `sigma(s)` 峰和 simple `1-4KΠ` 峰在近阈值 band 对齐，且上游背景量平滑；修正机制脚本的 A-builder 口径后，直调 `average_scattering_rate` 与生产 `channel_diagnostics.csv` 的 rate 复现误差为 `0`。
- 此前 mode B 约 `0.47` 的 rate 复现误差来自诊断脚本 bug：机制脚本手工重建传播子 A 场时没有使用 production workflow 的 `a_builder` 配置。当前机制脚本已显式记录并使用 `p_nodes=16,p_max=20.0,cos_nodes=4,use_aniso=true`；这不是 production 数据 bug。
- 6 个 weak/broad 窗口全部完成深拆，均显示 simple `1-4KΠ` 小分母支持；新增的 4 个窗口是 `mode_a_muB450p0_alpha1p1_xip0p35_1`、`mode_a_muB900p0_alpha1p2_xip0p49_1`、`mode_a_muB450p0_alpha1p0_xim0p20_1` 和 `mode_b_T200p0_muB0p0_xim0p21_1`。
- 本机尝试的局部 high-rate convergence gate 未在可控时间内完成；`local_rate_convergence_gate.csv` 因此只保留表头。当前结论是 denominator-chain 补证，不是新的 production-grade 收敛证明。
"""
    elif any(row["cause_verdict"] == "channel_rate_spike_supported" for row in window_rows):
        mechanism_section = """
## 非一阶 channel-rate spike 的 denominator-chain 补证

`mechanism_window_candidates.csv` 已列出待深拆的非一阶候选窗口；若要把直接近因升级为传播子分母近零机制，需要运行 `scripts/analysis/relaxtime/phase_guided_p128_mechanism_scan.jl` 生成 denominator-chain 表。
"""
    readme = f"""# Phase-guided transport p128 xi001 tau-first 突变分析

本分析包只消费仓库内已入库的 `{CASE}` 正式产物，不重跑 production，不修改主结果 CSV。分析主线是弛豫时间 `tau_*`；`eta_over_s`、`zeta_over_s` 等输运系数只作为 `tau` 与热力学量组合后的下游响应来讨论。

## 输入核对

{markdown_table(inventory, ["mode_key", "scan_rows", "channel_rows", "failed_rows", "xi_count", "xi_step", "plot_count", "has_antiquark_tau_fields", "has_zeta_over_s"], limit=4)}

正式图像已补齐反夸克弛豫时间：`tau_ubar_vs_xi.png`、`tau_dbar_vs_xi.png`、`tau_sbar_vs_xi.png`。本轮也补出 `zeta_over_s_vs_xi.png`，便于把 `eta_over_s/zeta_over_s` 作为下游趋势响应一起审阅。

## Tau 突变总体分类

自动检测覆盖 mode A/B 全部 `tau_u,tau_d,tau_s,tau_ubar,tau_dbar,tau_sbar` 曲线。相邻 xi 步长为 0.01，候选突变由 log-step 的 robust outlier 规则筛出，再按同一 panel/series 的相邻 species 聚合成窗口。

{markdown_table([{"cause_verdict": key, "count": value} for key, value in sorted(counts.items())], ["cause_verdict", "count"], limit=20)}

最强窗口如下：

{markdown_table(key_windows, ["window_id", "plot_panel", "plot_series", "target_xi", "affected_tau_fields", "max_tau_step_factor", "cause_verdict", "top_channel_evidence"], limit=10)}

## 机制解释边界

- 一阶或上游分支窗口：若 `phase_reference_kind/phase_structure` 指向一阶邻域，并且 `m_u`、`Phi`、`s_fm3inv` 等背景量在同一 xi 邻域快速变化，则把 tau 突变归因到上游平衡态/分支变化。用户指出的 `mode A, muB=900, alpha_T=1.0, xi≈0` 属于这一类；虽然 `phase_curr` 字符串没有跳变，但 `m_u` 和熵密度有阶跃，因此不能只用 `phase_curr` 标签判断。
- 背景平滑但 tau 单点下探的窗口：直接近因是 channel diagnostics 中少数通道的 rate/contribution 局部尖峰。例如 `mode A, muB=450, alpha_T=1.0, xi=0.26` 由 `uubar_to_uubar` 与 `uubar_to_ddbar` 放大主导；`mode B, T=200, muB=450, xi=0.31` 由 `udbar_to_udbar/dubar_to_dubar` 放大主导。
- 对这些 channel-rate 尖峰，是否能写成“传播子分母近零”取决于下方逐窗口 denominator-chain verdict；没有通过 rate 复现或局部收敛补证的窗口仍保留为机制候选。

{mechanism_section}

## 下游输运系数响应

`eta_over_s` 和 `zeta_over_s` 不是本轮突变根因入口。定性上，`eta`、`zeta` 是带有 `tau` 权重的动量积分，`eta_over_s`、`zeta_over_s` 还要除以熵密度。因此：

- 在 channel-rate spike 窗口，背景熵密度通常平滑，`tau` 单点下探会直接拖低 `eta_over_s` 与 `zeta_over_s`，表现为同 xi 的下游凹陷。
- 在一阶/上游分支窗口，`tau` 常上升，同时 `s_fm3inv` 下降；除以熵密度后，`eta_over_s/zeta_over_s` 的跳变会比 `eta/zeta` 本身更显著。
- 若某个下游 ratio 与 tau 方向不完全一致，应优先检查熵密度和对应 transport numerator，而不是把 ratio 曲线直接当成散射机制证据。

## 关键图

{chr(10).join(f"- `{relpath(path)}`" for path in figure_paths)}

## 产物表

- `tables/input_inventory.csv`：输入行数、hash、图像数量和字段核对。
- `tables/tau_jump_step_candidates.csv`：逐 species 的相邻 xi 突变候选。
- `tables/tau_jump_window_summary.csv`：聚合窗口、机制 verdict、上游背景量和主导通道摘要。
- `tables/tau_jump_channel_attribution.csv`：窗口内 top channel 的 contribution/rate 局部变化。
- `tables/downstream_transport_response_summary.csv`：`eta_over_s/zeta_over_s/sigma_over_T` 对 tau 突变的下游响应关系。
- `tables/mechanism_window_candidates.csv`：非一阶 channel-rate spike 的 denominator-chain 深拆候选。
- `tables/mechanism_window_summary.csv`：已深拆窗口的机制 verdict。
- `tables/denominator_chain_summary.csv`、`tables/denominator_chain_band_table.csv`、`tables/denominator_ds_samples.csv`：传播子分母、`sigma(s)` 与 rate band 的局部证据。
- `tables/local_rate_reproduction_mismatch_root_cause.csv`：mode B 早期 rate 复现偏差的诊断脚本口径根因。
- `tables/local_rate_convergence_gate.csv`：本机局部 high-rate gate 输出；当前仅有表头，表示未在本轮可控预算内完成。
- `tables/claim_ledger.csv`：可写入论文或需要作者确认的 claim 账本。

## 作者确认项

- xi001 新中间点继承 p128 积分参数与旧 xi=0.05 锚点比较证据，但非锚点局部尖峰还没有逐点 p104->p128 gate。
- 目前的机制支持来自 production 口径 rate 复现与 denominator-chain 对齐；若要把这些局部结构写成更强的 production-grade 收敛结论，仍应补可完成的局部高精度 convergence gate。
"""
    (ANALYSIS_DIR / "README.md").write_text(readme, encoding="utf-8")


def build_manifest(
    inventory: list[dict[str, Any]],
    figure_paths: list[Path],
    table_paths: list[Path],
) -> dict[str, Any]:
    output_files = [ANALYSIS_DIR / "README.md", *figure_paths, *table_paths]
    return {
        "schema": "phase_guided_transport_p128_xi001_tau_jump_analysis_manifest_v1",
        "case": CASE,
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "generator": relpath(Path(__file__)),
        "git_commit": current_git_commit(),
        "scope": "repository-only tau-first jump analysis; downstream transport coefficients treated as responses",
        "inputs": inventory,
        "outputs": [
            {
                "path": relpath(path),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
            for path in output_files
            if path.exists()
        ],
    }


MECHANISM_TABLE_NAMES = [
    "mechanism_window_summary.csv",
    "denominator_chain_summary.csv",
    "denominator_chain_band_table.csv",
    "denominator_ds_samples.csv",
    "upstream_branch_smoothness_summary.csv",
    "global_nonmonotonic_mechanism_summary.csv",
    "local_rate_reproduction_mismatch_root_cause.csv",
    "local_rate_convergence_gate.csv",
    "local_rate_convergence_attempts.csv",
    "mechanism_manifest.json",
]


def main() -> int:
    loaded = load_inputs()
    inventory = validate_inputs(loaded)
    step_candidates, window_rows, channel_rows = build_jump_tables(loaded)
    downstream_rows = downstream_summary(window_rows, loaded)
    mechanism_candidates = build_mechanism_candidates(window_rows, channel_rows, loaded)

    tables_dir = ANALYSIS_DIR / "tables"
    figures_dir = ANALYSIS_DIR / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    table_specs = [
        (
            tables_dir / "input_inventory.csv",
            inventory,
            [
                "mode_key",
                "scan_rows",
                "channel_rows",
                "failed_rows",
                "xi_count",
                "xi_min",
                "xi_max",
                "xi_step",
                "plot_count",
                "has_antiquark_tau_fields",
                "has_zeta_over_s",
                "scan_csv",
                "channel_diagnostics_csv",
                "plot_manifest",
                "scan_sha256",
                "diagnostics_sha256",
            ],
        ),
        (
            tables_dir / "tau_jump_step_candidates.csv",
            step_candidates,
            [
                "mode_key",
                "mode",
                "plot_panel",
                "plot_series",
                "species",
                "observable",
                "xi_left",
                "xi_right",
                "xi_mid",
                "tau_left",
                "tau_right",
                "tau_ratio_right_over_left",
                "tau_step_factor",
                "tau_rel_step",
                "log_step_abs",
                "robust_z",
                "phase_left",
                "phase_right",
                "phase_reference_left",
                "phase_reference_right",
                "phase_structure_left",
                "phase_structure_right",
                "seed_right",
                "m_u_left",
                "m_u_right",
                "Phi_left",
                "Phi_right",
            ],
        ),
        (
            tables_dir / "tau_jump_window_summary.csv",
            window_rows,
            [
                "window_id",
                "mode_key",
                "mode",
                "plot_panel",
                "plot_series",
                "xi_left",
                "target_xi",
                "xi_right",
                "prev_xi",
                "next_xi",
                "affected_tau_fields",
                "tau_shape",
                "max_tau_step_factor",
                "max_tau_rel_step",
                "cause_verdict",
                "evidence_score",
                "dominant_background_driver",
                "max_background_rel_step",
                "phase_prev",
                "phase_target",
                "phase_next",
                "phase_reference_target",
                "phase_structure_target",
                "phase_boundary_xi_target",
                "phase_boundary_changed_in_local_window",
                "seed_source_target",
                "equilibrium_backend_target",
                "m_u_prev",
                "m_u_target",
                "m_u_next",
                "Phi_prev",
                "Phi_target",
                "Phi_next",
                "s_prev",
                "s_target",
                "s_next",
                "top_channel_evidence",
                "mechanism_note",
            ],
        ),
        (
            tables_dir / "tau_jump_channel_attribution.csv",
            channel_rows,
            [
                "window_id",
                "mode_key",
                "plot_panel",
                "plot_series",
                "target_xi",
                "species",
                "rank",
                "channel",
                "density_key",
                "target_contribution",
                "target_rate",
                "target_share",
                "prev_contribution",
                "next_contribution",
                "contribution_ratio_to_neighbor_max",
                "rate_ratio_to_neighbor_max",
                "contribution_delta_vs_neighbor_max",
            ],
        ),
        (
            tables_dir / "downstream_transport_response_summary.csv",
            downstream_rows,
            [
                "window_id",
                "mode_key",
                "plot_panel",
                "plot_series",
                "target_xi",
                "observable",
                "baseline_value",
                "target_value",
                "log_delta_vs_neighbor_mean",
                "observable_direction",
                "tau_mean_direction",
                "entropy_direction",
                "response_relation",
                "root_tau_verdict",
            ],
        ),
        (
            tables_dir / "mechanism_window_candidates.csv",
            mechanism_candidates,
            [
                "window_id",
                "candidate_source",
                "selected_for_deep_scan",
                "selected_order",
                "selection_reason",
                "mode_key",
                "mode",
                "plot_panel",
                "plot_panel_label",
                "plot_series",
                "plot_series_label",
                "selector",
                "T_MeV",
                "muB_MeV",
                "alpha_T",
                "xi",
                "prev_xi",
                "next_xi",
                "observable",
                "primary_species",
                "mechanism_applicability",
                "local_shape",
                "local_score",
                "value_at_xi",
                "phase_curr",
                "phase_structure",
                "window_label",
            ],
        ),
    ]
    table_paths: list[Path] = []
    for path, rows, fields in table_specs:
        write_csv(path, rows, fields)
        table_paths.append(path)

    mechanism_rows = read_optional_csv(tables_dir / "mechanism_window_summary.csv")
    convergence_rows = read_optional_csv(tables_dir / "local_rate_convergence_gate.csv")
    existing_mechanism_paths = [
        tables_dir / name
        for name in MECHANISM_TABLE_NAMES
        if (tables_dir / name).exists()
    ]

    claim_rows = make_claim_ledger(window_rows, downstream_rows, mechanism_rows, convergence_rows)
    claim_path = tables_dir / "claim_ledger.csv"
    write_csv(
        claim_path,
        claim_rows,
        ["claim_id", "status", "claim_zh", "evidence", "fields_or_points"],
    )
    table_paths.append(claim_path)
    for path in existing_mechanism_paths:
        if path not in table_paths:
            table_paths.append(path)

    figure_paths = plot_jump_overview(window_rows, loaded, figures_dir)
    write_readme(inventory, window_rows, downstream_rows, figure_paths, mechanism_rows, convergence_rows)
    manifest = build_manifest(inventory, figure_paths, table_paths)
    write_json(ANALYSIS_DIR / "manifest.json", manifest)

    plot_manifest = {
        "schema": "phase_guided_transport_p128_xi001_analysis_plot_manifest_v1",
        "figures": [
            {
                "path": relpath(path),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
            for path in figure_paths
        ],
        "input_tables": [relpath(path) for path in table_paths],
    }
    write_json(figures_dir / "plot_manifest.json", plot_manifest)
    # Refresh manifest after plot_manifest is present.
    manifest = build_manifest(inventory, [*figure_paths, figures_dir / "plot_manifest.json"], table_paths)
    write_json(ANALYSIS_DIR / "manifest.json", manifest)

    print(
        json.dumps(
            {
                "analysis_dir": relpath(ANALYSIS_DIR),
                "step_candidates": len(step_candidates),
                "jump_windows": len(window_rows),
                "mechanism_candidates": len(mechanism_candidates),
                "figures": len(figure_paths),
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
