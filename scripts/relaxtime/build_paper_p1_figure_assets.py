#!/usr/bin/env python3
"""Build figure-ready assets for the P1 Mott/isentropic paper figures.

This script does not run the PNJL or meson solvers. It consumes existing
main-computation CSV artifacts and produces stable post-processing outputs.
CSV assets are written to the result tree; plot files and plot manifests are
written to the figure tree:

* mott_lines.csv: pion/kaon Mott lines in the T-muB plane
* isentropic_trajectories.csv: fixed-sigma path points with meson observables
* isentropic_mott_crossings.csv: path crossings with pion/kaon Mott conditions
* phase_overlay.csv: optional phase-line/CEP overlay data
* figures/p1_mott_phase_diagram.{png,pdf}
* figures/p1_isentropic_mott_paths.{png,pdf}
* plot_manifest.json
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

FM_TO_MEV = 197.327

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.pnjl.phase_reference_adapter import default_downstream_layer, load_phase_reference

CHANNELS = {
    "pi": {
        "mass": "M_pi",
        "gamma": "Gamma_pi",
        "gap": "gap_pi",
        "threshold": ["threshold_pi", "M_u_plus_M_d"],
        "mass_terms": ["m_u", "m_d"],
        "label": "pion",
    },
    "K": {
        "mass": "M_K",
        "gamma": "Gamma_K",
        "gap": "gap_K",
        "threshold": ["threshold_K", "M_u_plus_M_s"],
        "mass_terms": ["m_u", "m_s"],
        "label": "kaon",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build P1 paper figure assets from Mott, phase, and isentropic scan outputs."
    )
    parser.add_argument("--mott-grid-csv", "--mott-csv", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True, help="Figure output directory.")
    parser.add_argument(
        "--asset-dir",
        type=Path,
        default=None,
        help="CSV asset output directory. Defaults to the matching data/outputs/results/.../figure_assets path.",
    )
    parser.add_argument("--isentropic-csv", type=Path, action="append", default=[])
    parser.add_argument("--phase-dir", type=Path, action="append", default=[])
    parser.add_argument(
        "--phase-reference-root",
        type=Path,
        default=None,
        help="Directory containing reusable phase reference CSV files such as boundary_<tag>.csv.",
    )
    parser.add_argument(
        "--phase-reference-tag",
        default=None,
        help="Tag suffix for reusable phase reference CSV files.",
    )
    parser.add_argument(
        "--phase-reference-candidate-root",
        type=Path,
        default=None,
        help="Explicit Issue #130 candidate root for a diagnostic phase overlay.",
    )
    parser.add_argument(
        "--phase-reference-candidate-layer",
        choices=["strict", "derived", "render", "accepted"],
        default=None,
        help="Candidate layer; omitted selects accepted for v2 and strict for v1.",
    )
    parser.add_argument(
        "--phase-mu-scale",
        type=float,
        default=3.0,
        help="Scale phase artifact mu_q columns to the plotted mu_B axis. Default: 3.0.",
    )
    parser.add_argument("--formats", default="png,pdf", help="Comma-separated plot formats.")
    parser.add_argument("--skip-plots", action="store_true")
    return parser.parse_args()


def default_asset_dir(out_dir: Path) -> Path:
    parts = list(out_dir.parts)
    lowered = [part.lower() for part in parts]
    for idx in range(len(parts) - 2):
        if lowered[idx:idx + 3] == ["data", "outputs", "figures"]:
            parts[idx + 2] = "results"
            return Path(*parts) / "figure_assets"
    return out_dir / "figure_assets"


def _manifest_path(path: Path) -> str:
    return path.as_posix()


def read_scan_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)

    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        header: list[str] | None = None
        data_lines: list[str] = []
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if header is None:
                header = [part.strip() for part in stripped.split(",")]
            else:
                data_lines.append(line)

    if header is None:
        return rows

    reader = csv.DictReader(data_lines, fieldnames=header)
    for row in reader:
        rows.append({k: (v if v is not None else "") for k, v in row.items()})
    return rows


def read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}


def load_mott_effective_config(mott_grid_csv: Path) -> dict[str, Any]:
    candidates = [
        mott_grid_csv.parent / "mott_grid_combined_manifest.json",
        mott_grid_csv.with_suffix(".manifest.json"),
        mott_grid_csv.parent / "effective_config.json",
        mott_grid_csv.parent.parent / "effective_config.json",
    ]
    for path in candidates:
        payload = read_json(path)
        if payload:
            return payload
    return {}


def summarize_mott_source(mott_grid_csv: Path) -> dict[str, Any]:
    cfg = load_mott_effective_config(mott_grid_csv)
    mp = cfg.get("scan", {}).get("mott_phase", {}) if isinstance(cfg.get("scan"), dict) else {}
    return {
        "equilibrium_branch_mode": mp.get("equilibrium_branch_mode"),
        "equilibrium_selector_policy": mp.get("equilibrium_selector_policy"),
        "equilibrium_selector_tiebreak": mp.get("equilibrium_selector_tiebreak"),
    }


def finite_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        x = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def row_float(row: dict[str, str], key: str) -> float | None:
    return finite_float(row.get(key))


def row_float_first(row: dict[str, str], keys: Iterable[str]) -> float | None:
    for key in keys:
        val = row_float(row, key)
        if val is not None:
            return val
    return None


def fmt(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        return f"{value:.12g}"
    return str(value)


def number_tag(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return "nan"
    text = f"{value:g}"
    return text.replace(".", "p").replace("+", "")


def safe_tag(text: str) -> str:
    out = []
    for char in text:
        if char.isalnum() or char in ("-", "_"):
            out.append(char)
        elif char in ("=", ".", "+"):
            out.append(char.replace(".", "p").replace("+", ""))
        else:
            out.append("_")
    return "".join(out).strip("_") or "path"


def path_display_label(row: dict[str, Any]) -> str:
    sigma = row.get("sigma_target")
    if isinstance(sigma, (int, float)) and math.isfinite(float(sigma)):
        return f"s/nB={float(sigma):g}"
    label = row.get("path_label")
    return str(label) if label else str(row.get("path_id", "path"))


def threshold_value(row: dict[str, str], channel: str) -> tuple[float | None, str]:
    spec = CHANNELS[channel]
    for key in spec["threshold"]:
        val = row_float(row, key)
        if val is not None:
            return val, key

    terms = [row_float(row, key) for key in spec["mass_terms"]]
    if all(val is not None for val in terms):
        return float(sum(val for val in terms if val is not None)), "+".join(spec["mass_terms"])

    return None, ""


def gap_value(row: dict[str, str], channel: str) -> tuple[float | None, str]:
    spec = CHANNELS[channel]
    gap = row_float(row, spec["gap"])
    if gap is not None:
        return gap, spec["gap"]

    mass = row_float(row, spec["mass"])
    threshold, source = threshold_value(row, channel)
    if mass is not None and threshold is not None:
        return mass - threshold, f"{spec['mass']}-{source}"

    return None, ""


def interpolate(v0: float | None, v1: float | None, alpha: float) -> float | None:
    if v0 is None or v1 is None:
        return None
    return v0 + alpha * (v1 - v0)


def interpolate_row_value(row0: dict[str, str], row1: dict[str, str], key: str, alpha: float) -> float | None:
    return interpolate(row_float(row0, key), row_float(row1, key), alpha)


def first_gap_crossing(
    rows: Iterable[dict[str, str]],
    channel: str,
    order_key: str,
) -> dict[str, Any] | None:
    prepared = []
    for row in rows:
        order = row_float(row, order_key)
        gap, source = gap_value(row, channel)
        if order is None or gap is None:
            continue
        prepared.append((order, gap, source, row))

    prepared.sort(key=lambda item: item[0])
    if len(prepared) < 2:
        return None

    for idx in range(len(prepared) - 1):
        order0, gap0, source0, row0 = prepared[idx]
        order1, gap1, source1, row1 = prepared[idx + 1]

        if abs(gap0) <= 1e-14:
            alpha = 0.0
        elif gap0 * gap1 > 0.0:
            continue
        elif abs(gap1 - gap0) <= 1e-14:
            alpha = 0.0
        else:
            alpha = -gap0 / (gap1 - gap0)

        alpha = max(0.0, min(1.0, alpha))
        return {
            "alpha": alpha,
            "row0": row0,
            "row1": row1,
            "gap0": gap0,
            "gap1": gap1,
            "gap_source": source0 or source1,
            "order0": order0,
            "order1": order1,
        }

    return None


def build_mott_lines(rows: list[dict[str, str]], source_csv: Path) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, float, float], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        xi = row_float(row, "xi")
        mu_b = row_float(row, "muB_MeV")
        if xi is None or mu_b is None:
            continue
        for channel in CHANNELS:
            grouped[(channel, xi, mu_b)].append(row)

    out: list[dict[str, Any]] = []
    for (channel, xi, mu_b), group in sorted(grouped.items(), key=lambda item: (item[0][0], item[0][1], item[0][2])):
        crossing = first_gap_crossing(group, channel, "T_MeV")
        if crossing is None:
            continue

        alpha = crossing["alpha"]
        row0 = crossing["row0"]
        row1 = crossing["row1"]
        spec = CHANNELS[channel]
        t_mott = interpolate_row_value(row0, row1, "T_MeV", alpha)
        mass0 = row_float(row0, spec["mass"])
        mass1 = row_float(row1, spec["mass"])
        mass = interpolate(mass0, mass1, alpha)
        gamma = interpolate_row_value(row0, row1, spec["gamma"], alpha)
        thr0, _ = threshold_value(row0, channel)
        thr1, _ = threshold_value(row1, channel)
        threshold = interpolate(thr0, thr1, alpha)
        gap_jump = abs(crossing["gap1"] - crossing["gap0"])
        mass_jump = abs(mass1 - mass0) if mass0 is not None and mass1 is not None else None
        threshold_jump = abs(thr1 - thr0) if thr0 is not None and thr1 is not None else None
        bracket_kind = "branch_jump" if mass_jump is not None and mass_jump > 0.4 else "continuous_bracket"

        out.append(
            {
                "channel": channel,
                "xi": xi,
                "muB_MeV": mu_b,
                "T_Mott_MeV": t_mott,
                "gap_source": crossing["gap_source"],
                "bracket_T_low_MeV": row_float(row0, "T_MeV"),
                "bracket_T_high_MeV": row_float(row1, "T_MeV"),
                "gap_low_inv_fm": crossing["gap0"],
                "gap_high_inv_fm": crossing["gap1"],
                "bracket_gap_jump_inv_fm": gap_jump,
                "bracket_mass_jump_inv_fm": mass_jump,
                "bracket_threshold_jump_inv_fm": threshold_jump,
                "bracket_kind": bracket_kind,
                "M_inv_fm": mass,
                "M_MeV": mass * FM_TO_MEV if mass is not None else None,
                "Gamma_inv_fm": gamma,
                "Gamma_MeV": gamma * FM_TO_MEV if gamma is not None else None,
                "threshold_inv_fm": threshold,
                "threshold_MeV": threshold * FM_TO_MEV if threshold is not None else None,
                "source_csv": str(source_csv),
            }
        )

    return out


def path_group_key(row: dict[str, str], source_csv: Path) -> tuple[str, float, float, str]:
    xi = row_float(row, "xi")
    sigma = row_float(row, "sigma_target")
    label = row.get("path_label") or row.get("path_profile") or source_csv.stem
    path_id = f"xi{number_tag(xi)}__sigma{number_tag(sigma)}__{safe_tag(label)}"
    return path_id, xi if xi is not None else math.nan, sigma if sigma is not None else math.nan, label


def build_isentropic_assets(paths: list[Path]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    trajectories: list[dict[str, Any]] = []
    groups: dict[tuple[str, float, float, str], list[dict[str, str]]] = defaultdict(list)

    for path in paths:
        for row in read_scan_csv(path):
            group_key = path_group_key(row, path)
            groups[group_key].append(row)

            t_mev = row_float(row, "T_MeV")
            mu_b = row_float(row, "muB_MeV")
            order = row_float(row, "path_order_key")
            path_id, xi, sigma, label = group_key
            out_row: dict[str, Any] = {
                "path_id": path_id,
                "path_label": label,
                "sigma_target": sigma,
                "xi": xi,
                "path_order_key": order if order is not None else t_mev,
                "T_MeV": t_mev,
                "muB_MeV": mu_b,
                "source_csv": str(path),
            }
            for channel, spec in CHANNELS.items():
                mass = row_float(row, spec["mass"])
                gamma = row_float(row, spec["gamma"])
                gap, source = gap_value(row, channel)
                out_row[f"M_{channel}_MeV"] = mass * FM_TO_MEV if mass is not None else None
                out_row[f"Gamma_{channel}_MeV"] = gamma * FM_TO_MEV if gamma is not None else None
                out_row[f"gap_{channel}_inv_fm"] = gap
                out_row[f"gap_source_{channel}"] = source
            trajectories.append(out_row)

    crossings: list[dict[str, Any]] = []
    for (path_id, xi, sigma, label), group in sorted(groups.items(), key=lambda item: item[0]):
        for channel, spec in CHANNELS.items():
            crossing = first_gap_crossing(group, channel, "path_order_key")
            if crossing is None:
                crossing = first_gap_crossing(group, channel, "T_MeV")
            if crossing is None:
                continue

            alpha = crossing["alpha"]
            row0 = crossing["row0"]
            row1 = crossing["row1"]
            t_mev = interpolate_row_value(row0, row1, "T_MeV", alpha)
            mu_b = interpolate_row_value(row0, row1, "muB_MeV", alpha)
            order = interpolate_row_value(row0, row1, "path_order_key", alpha)
            mass = interpolate_row_value(row0, row1, spec["mass"], alpha)
            gamma = interpolate_row_value(row0, row1, spec["gamma"], alpha)

            crossings.append(
                {
                    "path_id": path_id,
                    "path_label": label,
                    "sigma_target": sigma,
                    "xi": xi,
                    "channel": channel,
                    "T_MeV": t_mev,
                    "muB_MeV": mu_b,
                    "path_order_key": order if order is not None else t_mev,
                    "M_MeV": mass * FM_TO_MEV if mass is not None else None,
                    "Gamma_MeV": gamma * FM_TO_MEV if gamma is not None else None,
                    "gap_source": crossing["gap_source"],
                    "gap_low_inv_fm": crossing["gap0"],
                    "gap_high_inv_fm": crossing["gap1"],
                }
            )

    trajectories.sort(key=lambda row: (fmt(row["path_id"]), row.get("path_order_key") or 0.0))
    crossings.sort(key=lambda row: (fmt(row["path_id"]), fmt(row["channel"])))
    return trajectories, crossings


def load_phase_overlay(phase_dirs: list[Path], phase_mu_scale: float) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for phase_dir in phase_dirs:
        summary = read_json(phase_dir / "phase_summary.json")
        xi = finite_float(summary.get("xi"))

        boundary_path = phase_dir / "first_order_boundary.csv"
        for row in read_scan_csv(boundary_path) if boundary_path.is_file() else []:
            mu_q = row_float(row, "mu_transition_MeV")
            t_mev = row_float(row, "T_MeV")
            if mu_q is None or t_mev is None:
                continue
            curve_parameter = row_float_first(row, ("curve_parameter", "plot_order_key", "T_MeV"))
            plot_order_key = row_float_first(row, ("plot_order_key", "curve_parameter", "T_MeV"))
            rows.append(
                {
                    "kind": "first_order",
                    "xi": xi,
                    "muB_MeV": mu_q * phase_mu_scale,
                    "T_MeV": t_mev,
                    "variable": "mu_transition_MeV",
                    "curve_parameter": curve_parameter,
                    "plot_order_key": plot_order_key,
                    "source_csv": str(boundary_path),
                }
            )

        spinodal_path = phase_dir / "spinodal.csv"
        for row in read_scan_csv(spinodal_path) if spinodal_path.is_file() else []:
            t_mev = row_float(row, "T_MeV")
            if t_mev is None:
                continue
            for kind, col in (
                ("spinodal_hadron", "mu_spinodal_hadron_MeV"),
                ("spinodal_quark", "mu_spinodal_quark_MeV"),
            ):
                mu_q = row_float(row, col)
                if mu_q is None:
                    continue
                curve_parameter = row_float_first(row, ("curve_parameter", "plot_order_key", "T_MeV"))
                plot_order_key = row_float_first(row, ("plot_order_key", "curve_parameter", "T_MeV"))
                rows.append(
                    {
                        "kind": kind,
                        "xi": xi,
                        "muB_MeV": mu_q * phase_mu_scale,
                        "T_MeV": t_mev,
                        "variable": col,
                        "curve_parameter": curve_parameter,
                        "plot_order_key": plot_order_key,
                        "source_csv": str(spinodal_path),
                    }
                )

        crossover_path = phase_dir / "crossover_line.csv"
        for row in read_scan_csv(crossover_path) if crossover_path.is_file() else []:
            mu_q = row_float(row, "mu_MeV")
            t_mev = row_float(row, "T_crossover_MeV")
            if mu_q is None or t_mev is None:
                continue
            curve_parameter = row_float_first(row, ("curve_parameter", "plot_order_key", "mu_MeV"))
            plot_order_key = row_float_first(row, ("plot_order_key", "curve_parameter", "mu_MeV"))
            rows.append(
                {
                    "kind": "crossover",
                    "xi": xi,
                    "muB_MeV": mu_q * phase_mu_scale,
                    "T_MeV": t_mev,
                    "variable": row.get("variable", ""),
                    "curve_parameter": curve_parameter,
                    "plot_order_key": plot_order_key,
                    "source_csv": str(crossover_path),
                }
            )

        cep = summary.get("cep", {}) if isinstance(summary.get("cep"), dict) else {}
        if cep.get("found") is True:
            t_cep = finite_float(cep.get("T_cep_MeV"))
            mu_b = finite_float(cep.get("muB_cep_MeV"))
            if mu_b is None:
                mu_q = finite_float(cep.get("muq_cep_MeV") or cep.get("mu_cep_MeV"))
                mu_b = mu_q * phase_mu_scale if mu_q is not None else None
            if t_cep is not None and mu_b is not None:
                rows.append(
                    {
                        "kind": "cep",
                        "xi": xi,
                        "muB_MeV": mu_b,
                        "T_MeV": t_cep,
                        "variable": "CEP",
                        "curve_parameter": t_cep,
                        "plot_order_key": t_cep,
                        "source_csv": str(phase_dir / "phase_summary.json"),
                    }
                )

    rows.sort(key=lambda row: (fmt(row["kind"]), row.get("xi") or 0.0, row.get("plot_order_key") or 0.0, row.get("muB_MeV") or 0.0))
    return rows


def load_phase_reference_overlay(
    phase_reference_root: Path | None,
    phase_reference_tag: str | None,
    phase_mu_scale: float,
) -> list[dict[str, Any]]:
    if phase_reference_root is None:
        return []
    if not phase_reference_tag:
        raise ValueError("--phase-reference-tag is required when --phase-reference-root is provided")

    rows: list[dict[str, Any]] = []
    found_any = False
    tag = phase_reference_tag

    boundary_path = phase_reference_root / f"boundary_{tag}.csv"
    if boundary_path.is_file():
        found_any = True
        for row in read_scan_csv(boundary_path):
            xi = row_float(row, "xi")
            mu_q = row_float(row, "mu_transition_MeV")
            t_mev = row_float(row, "T_MeV")
            if xi is None or mu_q is None or t_mev is None:
                continue
            curve_parameter = row_float_first(row, ("curve_parameter", "plot_order_key", "T_MeV"))
            plot_order_key = row_float_first(row, ("plot_order_key", "curve_parameter", "T_MeV"))
            rows.append(
                {
                    "kind": "first_order",
                    "xi": xi,
                    "muB_MeV": mu_q * phase_mu_scale,
                    "T_MeV": t_mev,
                    "variable": "mu_transition_MeV",
                    "curve_parameter": curve_parameter,
                    "plot_order_key": plot_order_key,
                    "source_csv": str(boundary_path),
                }
            )

    spinodal_path = phase_reference_root / f"spinodals_{tag}.csv"
    if spinodal_path.is_file():
        found_any = True
        for row in read_scan_csv(spinodal_path):
            xi = row_float(row, "xi")
            t_mev = row_float(row, "T_MeV")
            if xi is None or t_mev is None:
                continue
            curve_parameter = row_float_first(row, ("curve_parameter", "plot_order_key", "T_MeV"))
            plot_order_key = row_float_first(row, ("plot_order_key", "curve_parameter", "T_MeV"))
            for kind, col in (
                ("spinodal_hadron", "mu_spinodal_hadron_MeV"),
                ("spinodal_quark", "mu_spinodal_quark_MeV"),
            ):
                mu_q = row_float(row, col)
                if mu_q is None:
                    continue
                rows.append(
                    {
                        "kind": kind,
                        "xi": xi,
                        "muB_MeV": mu_q * phase_mu_scale,
                        "T_MeV": t_mev,
                        "variable": col,
                        "curve_parameter": curve_parameter,
                        "plot_order_key": plot_order_key,
                        "source_csv": str(spinodal_path),
                    }
                )

    crossover_path = phase_reference_root / f"crossover_{tag}.csv"
    if crossover_path.is_file():
        found_any = True
        for row in read_scan_csv(crossover_path):
            xi = row_float(row, "xi")
            mu_q = row_float(row, "mu_MeV")
            t_mev = row_float(row, "T_crossover_MeV")
            if xi is None or mu_q is None or t_mev is None:
                continue
            curve_parameter = row_float_first(row, ("curve_parameter", "plot_order_key", "mu_MeV"))
            plot_order_key = row_float_first(row, ("plot_order_key", "curve_parameter", "mu_MeV"))
            rows.append(
                {
                    "kind": "crossover",
                    "xi": xi,
                    "muB_MeV": mu_q * phase_mu_scale,
                    "T_MeV": t_mev,
                    "variable": row.get("variable", ""),
                    "curve_parameter": curve_parameter,
                    "plot_order_key": plot_order_key,
                    "source_csv": str(crossover_path),
                }
            )

    cep_path = phase_reference_root / f"cep_{tag}.csv"
    if cep_path.is_file():
        found_any = True
        for row in read_scan_csv(cep_path):
            xi = row_float(row, "xi")
            t_cep = row_float(row, "T_CEP_MeV")
            mu_b = row_float(row, "muB_CEP_MeV")
            if mu_b is None:
                mu_q = row_float(row, "muq_CEP_MeV")
                mu_b = mu_q * phase_mu_scale if mu_q is not None else None
            if xi is None or t_cep is None or mu_b is None:
                continue
            rows.append(
                {
                    "kind": "cep",
                    "xi": xi,
                    "muB_MeV": mu_b,
                    "T_MeV": t_cep,
                    "variable": "CEP",
                    "curve_parameter": t_cep,
                    "plot_order_key": t_cep,
                    "source_csv": str(cep_path),
                }
            )

    if not found_any:
        raise FileNotFoundError(f"no phase reference CSV files found for tag {tag!r} under {phase_reference_root}")

    rows.sort(key=lambda row: (fmt(row["kind"]), row.get("xi") or 0.0, row.get("plot_order_key") or 0.0, row.get("muB_MeV") or 0.0))
    return rows


def load_phase_candidate_overlay(
    phase_reference_root: Path | None,
    phase_reference_layer: str,
    phase_mu_scale: float,
) -> list[dict[str, Any]]:
    """Load an explicit Issue #130 candidate for diagnostic plotting only."""
    if phase_reference_root is None:
        return []
    bundle = load_phase_reference(phase_reference_root, layer=phase_reference_layer)
    layer_dir = bundle.root / phase_reference_layer / "tables"
    table_names = {
        "strict": {
            "boundary": "maxwell_surface_strict_reference_v1.csv",
            "spinodals": "spinodal_surface_strict_reference_v1.csv",
            "crossover": "crossover_surface_strict_reference_v1.csv",
            "cep": "cep_boundary_strict_reference_v1.csv",
        },
        "derived": {
            "boundary": "maxwell_surface_derived_reference_v1.csv",
            "spinodals": "spinodal_surface_derived_reference_v1.csv",
            "crossover": "crossover_surface_derived_reference_v1.csv",
            "cep": "cep_boundary_derived_reference_v1.csv",
        },
        "render": {
            "boundary": "maxwell_surface_render.csv",
            "spinodals": "spinodal_surface_render.csv",
            "crossover": "crossover_surface_render.csv",
            "cep": "cep_boundary_render.csv",
        },
        "accepted": {
            "boundary": "maxwell_surface_accepted_phase_map_v1.csv",
            "spinodals": "spinodal_surface_accepted_phase_map_v1.csv",
            "crossover": "crossover_surface_accepted_phase_map_v1.csv",
            "cep": "cep_boundary_accepted_phase_map_v1.csv",
        },
    }
    source_paths = {
        table: layer_dir / filename
        for table, filename in table_names[phase_reference_layer].items()
    }
    rows: list[dict[str, Any]] = []
    for row in bundle.tables.get("boundary", ()):
        rows.append(
            {
                "kind": "first_order",
                "xi": row["xi"],
                "muB_MeV": row["muq_MeV"] * phase_mu_scale,
                "T_MeV": row["T_MeV"],
                "variable": "mu_transition_MeV",
                "curve_parameter": row["T_MeV"],
                "plot_order_key": row["T_MeV"],
                "source_csv": str(source_paths["boundary"]),
                "candidate_status": row.get("status", ""),
                "candidate_certified": row.get("certified", False),
            }
        )
    for row in bundle.tables.get("spinodals", ()):
        for kind, key in (("spinodal_hadron", "muq_spinodal_hadron_MeV"), ("spinodal_quark", "muq_spinodal_quark_MeV")):
            rows.append(
                {
                    "kind": kind,
                    "xi": row["xi"],
                    "muB_MeV": row[key] * phase_mu_scale,
                    "T_MeV": row["T_MeV"],
                    "variable": key,
                    "curve_parameter": row["T_MeV"],
                    "plot_order_key": row["T_MeV"],
                    "source_csv": str(source_paths["spinodals"]),
                    "candidate_status": row.get("status", ""),
                    "candidate_certified": row.get("certified", False),
                }
            )
    for row in bundle.tables.get("crossover", ()):
        if row.get("physical_region", "").lower() not in ("", "crossover_below_cep"):
            continue
        rows.append(
            {
                "kind": "crossover",
                "xi": row["xi"],
                "muB_MeV": row["muq_MeV"] * phase_mu_scale,
                "T_MeV": row["T_MeV"],
                "variable": "mu_MeV",
                "curve_parameter": row["muq_MeV"],
                "plot_order_key": row["muq_MeV"],
                "source_csv": str(source_paths["crossover"]),
                "candidate_status": row.get("status", ""),
                "candidate_certified": row.get("certified", False),
            }
        )
    for row in bundle.tables.get("cep", ()):
        rows.append(
            {
                "kind": "cep",
                "xi": row["xi"],
                "muB_MeV": row["muq_CEP_proxy_MeV"] * phase_mu_scale,
                "T_MeV": row["T_midpoint_MeV"],
                "variable": "CEP",
                "curve_parameter": row["T_midpoint_MeV"],
                "plot_order_key": row["T_midpoint_MeV"],
                "source_csv": str(source_paths["cep"]),
                "candidate_status": row.get("status", ""),
                "candidate_certified": row.get("certified", False),
            }
        )
    return rows


def dedupe_phase_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    deduped: list[dict[str, Any]] = []
    seen: set[tuple[Any, ...]] = set()
    for row in rows:
        key = (
            row.get("kind"),
            None if row.get("xi") is None else round(float(row["xi"]), 10),
            None if row.get("muB_MeV") is None else round(float(row["muB_MeV"]), 8),
            None if row.get("T_MeV") is None else round(float(row["T_MeV"]), 8),
            row.get("variable"),
        )
        if key in seen:
            continue
        seen.add(key)
        deduped.append(row)
    deduped.sort(key=lambda row: (fmt(row["kind"]), row.get("xi") or 0.0, row.get("plot_order_key") or 0.0, row.get("muB_MeV") or 0.0))
    return deduped


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field)) for field in fields})


def grouped_rows(rows: list[dict[str, Any]], *keys: str) -> dict[Any, list[dict[str, Any]]]:
    groups: dict[Any, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        key = row.get(keys[0]) if len(keys) == 1 else tuple(row.get(key) for key in keys)
        groups[key].append(row)
    return groups


def xi_figure_tag(xi: float) -> str:
    return safe_tag("xi_" + number_tag(xi).replace("-", "m"))


def draw_mott_phase_panel(
    ax,
    mott_lines: list[dict[str, Any]],
    phase_overlay: list[dict[str, Any]],
    color_by_xi: dict[float, str],
    *,
    include_xi_in_label: bool,
) -> None:
    channel_style = {"pi": "-", "K": "--"}
    for (channel, xi), group in grouped_rows(mott_lines, "channel", "xi").items():
        group = sorted(group, key=lambda row: row.get("muB_MeV") or 0.0)
        label = CHANNELS[str(channel)]["label"]
        if include_xi_in_label:
            label = f"{label} xi={xi:g}"
        ax.plot(
            [row["muB_MeV"] for row in group],
            [row["T_Mott_MeV"] for row in group],
            linestyle=channel_style.get(str(channel), "-"),
            color=color_by_xi.get(xi, "#333333"),
            linewidth=1.4,
            label=label,
        )
        jump_group = [row for row in group if row.get("bracket_kind") == "branch_jump"]
        if jump_group:
            ax.scatter(
                [row["muB_MeV"] for row in jump_group],
                [row["T_Mott_MeV"] for row in jump_group],
                marker="x",
                s=24,
                color=color_by_xi.get(xi, "#333333"),
                linewidth=0.8,
                label="branch-jump bracket",
                zorder=6,
            )

    phase_style = {
        "first_order": ("-", 1.0, "first order"),
        "crossover": ("-.", 1.0, "chiral crossover"),
        "spinodal_hadron": (":", 0.9, "spinodal hadron"),
        "spinodal_quark": (":", 0.9, "spinodal quark"),
    }
    for (kind, xi), group in grouped_rows(phase_overlay, "kind", "xi").items():
        fallback_sort_field = "muB_MeV" if kind == "crossover" else "T_MeV"
        group = sorted(group, key=lambda row: row.get("plot_order_key") or row.get(fallback_sort_field) or 0.0)
        color = color_by_xi.get(xi, "#666666")
        if kind == "cep":
            label = "CEP"
            if include_xi_in_label:
                label = f"CEP xi={xi:g}" if xi is not None else "CEP"
            ax.scatter(
                [row["muB_MeV"] for row in group],
                [row["T_MeV"] for row in group],
                marker="*",
                s=80,
                color=color,
                edgecolor="#000000",
                linewidth=0.4,
                label=label,
                zorder=5,
            )
        elif kind in phase_style:
            linestyle, linewidth, label = phase_style[str(kind)]
            if include_xi_in_label:
                label = f"{label} xi={xi:g}" if xi is not None else label
            ax.plot(
                [row["muB_MeV"] for row in group],
                [row["T_MeV"] for row in group],
                color=color,
                linestyle=linestyle,
                linewidth=linewidth,
                alpha=0.75,
                label=label,
            )

    ax.set_xlabel(r"$\mu_B$ (MeV)")
    ax.set_ylabel(r"$T$ (MeV)")
    ax.minorticks_on()
    ax.tick_params(direction="in", top=True, right=True)


def add_dedup_legend(ax, *, fontsize: int = 7, ncol: int = 2) -> None:
    handles, labels = ax.get_legend_handles_labels()
    dedup: dict[str, Line2D] = {}
    for handle, label in zip(handles, labels):
        dedup.setdefault(label, handle)
    ax.legend(dedup.values(), dedup.keys(), frameon=False, fontsize=fontsize, ncol=ncol)


def plot_assets(
    out_dir: Path,
    formats: list[str],
    mott_lines: list[dict[str, Any]],
    phase_overlay: list[dict[str, Any]],
    trajectories: list[dict[str, Any]],
    crossings: list[dict[str, Any]],
) -> list[str]:
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    figure_dir = out_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    written: list[str] = []

    xi_values = sorted({row["xi"] for row in mott_lines if row.get("xi") is not None})
    palette = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377"]
    color_by_xi = {xi: palette[idx % len(palette)] for idx, xi in enumerate(xi_values)}

    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    draw_mott_phase_panel(ax, mott_lines, phase_overlay, color_by_xi, include_xi_in_label=True)
    add_dedup_legend(ax, fontsize=7, ncol=2)
    for fmt_name in formats:
        target = figure_dir / f"p1_mott_phase_diagram.{fmt_name}"
        fig.savefig(target, dpi=300, bbox_inches="tight")
        written.append(_manifest_path(target))
    plt.close(fig)

    for xi in xi_values:
        xi_mott = [row for row in mott_lines if row.get("xi") == xi]
        xi_phase = [row for row in phase_overlay if row.get("xi") == xi]
        if not xi_mott and not xi_phase:
            continue
        fig, ax = plt.subplots(figsize=(5.2, 4.2))
        draw_mott_phase_panel(ax, xi_mott, xi_phase, {xi: color_by_xi.get(xi, "#333333")}, include_xi_in_label=False)
        ax.set_title(rf"$\xi={xi:g}$")
        add_dedup_legend(ax, fontsize=7, ncol=1)
        tag = xi_figure_tag(float(xi))
        for fmt_name in formats:
            target = figure_dir / f"p1_mott_phase_diagram_{tag}.{fmt_name}"
            fig.savefig(target, dpi=300, bbox_inches="tight")
            written.append(_manifest_path(target))
        plt.close(fig)

    if not trajectories:
        return written

    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.0))
    ax0, ax1 = axes
    path_colors = {}
    for idx, path_id in enumerate(sorted({row["path_id"] for row in trajectories})):
        path_colors[path_id] = palette[idx % len(palette)]

    for path_id, group in grouped_rows(trajectories, "path_id").items():
        group = sorted(group, key=lambda row: row.get("path_order_key") or 0.0)
        ax0.plot(
            [row["muB_MeV"] for row in group],
            [row["T_MeV"] for row in group],
            color=path_colors[path_id],
            linewidth=1.3,
            label=path_display_label(group[0]),
        )

    marker_by_channel = {"pi": "o", "K": "s"}
    for row in crossings:
        ax0.scatter(
            row["muB_MeV"],
            row["T_MeV"],
            marker=marker_by_channel.get(str(row["channel"]), "o"),
            color=path_colors.get(row["path_id"], "#333333"),
            edgecolor="white",
            linewidth=0.6,
            s=38,
            zorder=5,
        )

    for path_id, group in grouped_rows(trajectories, "path_id").items():
        group = sorted(group, key=lambda row: row.get("path_order_key") or 0.0)
        x_vals = [row["T_MeV"] for row in group]
        color = path_colors[path_id]
        ax1.plot(x_vals, [row.get("M_pi_MeV") for row in group], color=color, linestyle="-", linewidth=1.1)
        ax1.plot(x_vals, [row.get("M_K_MeV") for row in group], color=color, linestyle="--", linewidth=1.1)
        ax1.plot(x_vals, [row.get("Gamma_pi_MeV") for row in group], color=color, linestyle=":", linewidth=1.0)
        ax1.plot(x_vals, [row.get("Gamma_K_MeV") for row in group], color=color, linestyle="-.", linewidth=1.0)

    ax0.set_xlabel(r"$\mu_B$ (MeV)")
    ax0.set_ylabel(r"$T$ (MeV)")
    ax1.set_xlabel(r"$T$ (MeV)")
    ax1.set_ylabel("Mass / width (MeV)")
    style_handles = [
        Line2D([], [], color="#333333", linestyle="-", label=r"$M_\pi$"),
        Line2D([], [], color="#333333", linestyle="--", label=r"$M_K$"),
        Line2D([], [], color="#333333", linestyle=":", label=r"$\Gamma_\pi$"),
        Line2D([], [], color="#333333", linestyle="-.", label=r"$\Gamma_K$"),
    ]
    ax1.legend(handles=style_handles, frameon=False, fontsize=7)
    ax0.legend(frameon=False, fontsize=6)
    for ax in axes:
        ax.minorticks_on()
        ax.tick_params(direction="in", top=True, right=True)
    for fmt_name in formats:
        target = figure_dir / f"p1_isentropic_mott_paths.{fmt_name}"
        fig.savefig(target, dpi=300, bbox_inches="tight")
        written.append(_manifest_path(target))
    plt.close(fig)

    return written


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    asset_dir = args.asset_dir if args.asset_dir is not None else default_asset_dir(args.out_dir)
    asset_dir.mkdir(parents=True, exist_ok=True)

    mott_rows = read_scan_csv(args.mott_grid_csv)
    mott_lines = build_mott_lines(mott_rows, args.mott_grid_csv)
    trajectories, isentropic_crossings = build_isentropic_assets(args.isentropic_csv)
    if args.phase_reference_root is not None and args.phase_reference_candidate_root is not None:
        raise ValueError("choose either --phase-reference-root or --phase-reference-candidate-root")
    candidate_layer = args.phase_reference_candidate_layer
    if args.phase_reference_candidate_root is not None and candidate_layer is None:
        candidate_layer = default_downstream_layer(args.phase_reference_candidate_root)
    candidate_diagnostics = None
    if args.phase_reference_candidate_root is not None:
        candidate_diagnostics = dict(
            load_phase_reference(
                args.phase_reference_candidate_root,
                layer=candidate_layer,
            ).diagnostics
        )
    phase_overlay = dedupe_phase_rows(
        load_phase_overlay(args.phase_dir, args.phase_mu_scale)
        + load_phase_reference_overlay(args.phase_reference_root, args.phase_reference_tag, args.phase_mu_scale)
        + load_phase_candidate_overlay(
            args.phase_reference_candidate_root,
            candidate_layer,
            args.phase_mu_scale,
        )
    )

    mott_fields = [
        "channel",
        "xi",
        "muB_MeV",
        "T_Mott_MeV",
        "gap_source",
        "bracket_T_low_MeV",
        "bracket_T_high_MeV",
        "gap_low_inv_fm",
        "gap_high_inv_fm",
        "bracket_gap_jump_inv_fm",
        "bracket_mass_jump_inv_fm",
        "bracket_threshold_jump_inv_fm",
        "bracket_kind",
        "M_inv_fm",
        "M_MeV",
        "Gamma_inv_fm",
        "Gamma_MeV",
        "threshold_inv_fm",
        "threshold_MeV",
        "source_csv",
    ]
    trajectory_fields = [
        "path_id",
        "path_label",
        "sigma_target",
        "xi",
        "path_order_key",
        "T_MeV",
        "muB_MeV",
        "M_pi_MeV",
        "Gamma_pi_MeV",
        "gap_pi_inv_fm",
        "gap_source_pi",
        "M_K_MeV",
        "Gamma_K_MeV",
        "gap_K_inv_fm",
        "gap_source_K",
        "source_csv",
    ]
    crossing_fields = [
        "path_id",
        "path_label",
        "sigma_target",
        "xi",
        "channel",
        "T_MeV",
        "muB_MeV",
        "path_order_key",
        "M_MeV",
        "Gamma_MeV",
        "gap_source",
        "gap_low_inv_fm",
        "gap_high_inv_fm",
    ]
    phase_fields = [
        "kind", "xi", "muB_MeV", "T_MeV", "variable", "curve_parameter", "plot_order_key",
        "source_csv", "candidate_status", "candidate_certified",
    ]

    assets = {
        "mott_lines": _manifest_path(asset_dir / "mott_lines.csv"),
        "isentropic_trajectories": _manifest_path(asset_dir / "isentropic_trajectories.csv"),
        "isentropic_mott_crossings": _manifest_path(asset_dir / "isentropic_mott_crossings.csv"),
        "phase_overlay": _manifest_path(asset_dir / "phase_overlay.csv"),
    }
    write_csv(Path(assets["mott_lines"]), mott_lines, mott_fields)
    write_csv(Path(assets["isentropic_trajectories"]), trajectories, trajectory_fields)
    write_csv(Path(assets["isentropic_mott_crossings"]), isentropic_crossings, crossing_fields)
    write_csv(Path(assets["phase_overlay"]), phase_overlay, phase_fields)

    figure_paths: list[str] = []
    formats = [part.strip().lower() for part in args.formats.split(",") if part.strip()]
    if not args.skip_plots:
        figure_paths = plot_assets(args.out_dir, formats, mott_lines, phase_overlay, trajectories, isentropic_crossings)

    manifest = {
        "schema_version": "paper_p1_assets_v1",
        "inputs": {
            "mott_grid_csv": _manifest_path(args.mott_grid_csv),
            "mott_source": summarize_mott_source(args.mott_grid_csv),
            "isentropic_csv": [_manifest_path(path) for path in args.isentropic_csv],
            "phase_dir": [_manifest_path(path) for path in args.phase_dir],
            "phase_reference_root": _manifest_path(args.phase_reference_root) if args.phase_reference_root is not None else None,
            "phase_reference_tag": args.phase_reference_tag,
            "phase_reference_candidate_root": _manifest_path(args.phase_reference_candidate_root) if args.phase_reference_candidate_root is not None else None,
            "phase_reference_candidate_layer": candidate_layer,
            "phase_reference_candidate_diagnostics": candidate_diagnostics,
            "phase_mu_scale": args.phase_mu_scale,
        },
        "counts": {
            "mott_line_points": len(mott_lines),
            "isentropic_points": len(trajectories),
            "isentropic_crossings": len(isentropic_crossings),
            "phase_overlay_points": len(phase_overlay),
            "mott_line_quality": dict(Counter(fmt(row.get("bracket_kind", "unknown")) for row in mott_lines)),
        },
        "asset_dir": _manifest_path(asset_dir),
        "assets": assets,
        "figures": figure_paths,
    }
    manifest_file = args.out_dir / "plot_manifest.json"
    manifest_file.write_text(json.dumps(manifest, indent=2, ensure_ascii=True), encoding="utf-8")

    print(f"Wrote P1 paper figure assets under: {args.out_dir}")
    print(f"  csv_asset_dir={asset_dir}")
    print(f"  mott_line_points={len(mott_lines)}")
    print(f"  isentropic_points={len(trajectories)}")
    print(f"  isentropic_crossings={len(isentropic_crossings)}")
    print(f"  phase_overlay_points={len(phase_overlay)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
