#!/usr/bin/env python3
"""Select a sparse, provenance-labelled full-gate candidate set.

The selector is solver-free.  It consumes the exact merged screening CSV and
uses only screened rows for finite-difference gradients.  Candidate points are
selected near the configured chemical-freezeout curve, large ratio gradients,
screening-mask boundaries, and the accepted historical phase-reference lines.
The output is a diagnostic dispatch manifest; it never promotes screening
values to production data and never fills failed points with zero.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[3]
PLOTTER_PATH = ROOT / "scripts" / "analysis" / "relaxtime" / "plot_charged_gbu_contour_screening.py"
SCHEMA_VERSION = "charged_gbu_contour_candidate_manifest_v1"


def _load_plotter() -> Any:
    spec = importlib.util.spec_from_file_location("charged_gbu_contour_screening_plot", PLOTTER_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot import plotter helper: {PLOTTER_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _finite(value: str, field: str, path: Path) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{path}: {field} is not numeric: {value!r}") from exc
    if not math.isfinite(result):
        raise ValueError(f"{path}: {field} is not finite: {value!r}")
    return result


def _input_csv(path: Path) -> Path:
    if path.is_dir():
        path = path / "contour_points_merged.csv"
    if not path.is_file():
        raise ValueError(f"screening CSV does not exist: {path}")
    return path


def load_rows(path: str | Path) -> tuple[list[dict[str, Any]], list[float], list[float]]:
    csv_path = _input_csv(Path(path).resolve())
    rows: list[dict[str, Any]] = []
    with csv_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"T_MeV", "muB_MeV", "status", "Kplus_over_pi_plus"}
        missing = sorted(required - set(reader.fieldnames or ()))
        if missing:
            raise ValueError(f"{csv_path}: missing columns: {', '.join(missing)}")
        for raw in reader:
            row = dict(raw)
            row["T_MeV"] = _finite(raw["T_MeV"], "T_MeV", csv_path)
            row["muB_MeV"] = _finite(raw["muB_MeV"], "muB_MeV", csv_path)
            row["status"] = str(raw["status"] or "").strip()
            value = str(raw.get("Kplus_over_pi_plus", "")).strip()
            row["Kplus_over_pi_plus"] = None if value in {"", "nan", "NaN", "null"} else _finite(value, "Kplus_over_pi_plus", csv_path)
            rows.append(row)
    if not rows:
        raise ValueError(f"{csv_path}: screening CSV is empty")
    key_map: dict[tuple[float, float], dict[str, Any]] = {}
    for row in rows:
        key = (row["T_MeV"], row["muB_MeV"])
        if key in key_map:
            raise ValueError(f"duplicate screening key: {key}")
        key_map[key] = row
    T_grid = sorted({key[0] for key in key_map})
    muB_grid = sorted({key[1] for key in key_map})
    expected = {(T, muB) for T in T_grid for muB in muB_grid}
    if set(key_map) != expected:
        raise ValueError("screening candidate selection requires a complete rectangular grid")
    return rows, T_grid, muB_grid


def _screened(row: dict[str, Any] | None) -> bool:
    return row is not None and row["status"] == "screened" and row["Kplus_over_pi_plus"] is not None


def _nearest(rows: Iterable[dict[str, Any]], T: float, muB: float) -> dict[str, Any] | None:
    available = [row for row in rows if row["status"] == "screened"]
    if not available:
        return None
    return min(available, key=lambda row: (row["T_MeV"] - T) ** 2 + (row["muB_MeV"] - muB) ** 2)


def _point_record(row: dict[str, Any], category: str, source: str, distance: float | None = None, **extra: Any) -> dict[str, Any]:
    result: dict[str, Any] = {
        "T_MeV": row["T_MeV"],
        "muB_MeV": row["muB_MeV"],
        "status": row["status"],
        "screening_ratio_Kplus_over_pi_plus": row["Kplus_over_pi_plus"],
        "category": category,
        "source": source,
    }
    if distance is not None:
        result["distance_squared_MeV2"] = distance
    result.update(extra)
    return result


def select_candidates(
    rows: list[dict[str, Any]],
    T_grid: list[float],
    muB_grid: list[float],
    reference_lines: dict[str, Any],
    *,
    max_gradient: int = 6,
    max_freezeout: int = 10,
    max_phase: int = 6,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    row_map = {(row["T_MeV"], row["muB_MeV"]): row for row in rows}
    selected: dict[tuple[float, float], dict[str, Any]] = {}
    provenance: dict[tuple[float, float], list[str]] = {}

    def add(row: dict[str, Any] | None, category: str, source: str, **extra: Any) -> None:
        if row is None or row["status"] != "screened":
            return
        key = (row["T_MeV"], row["muB_MeV"])
        if key not in selected:
            selected[key] = _point_record(row, category, source, **extra)
        provenance.setdefault(key, []).append(f"{category}:{source}")

    # Freeze-out anchors are nearest screened grid points; the requested
    # freeze-out coordinates are retained so snapping error remains auditable.
    freezeout_points = reference_lines.get("freezeout", {}).get("points", [])
    if freezeout_points:
        stride = max(1, len(freezeout_points) // max_freezeout)
        anchors = freezeout_points[::stride][:max_freezeout]
        for point in anchors:
            nearest = _nearest(rows, float(point["T_MeV"]), float(point["muB_MeV"]))
            if nearest is not None:
                distance = (nearest["T_MeV"] - float(point["T_MeV"])) ** 2 + (nearest["muB_MeV"] - float(point["muB_MeV"])) ** 2
                add(nearest, "freezeout", "default_freezeout_curve", target_T_MeV=point["T_MeV"], target_muB_MeV=point["muB_MeV"], distance_squared_MeV2=distance, sqrt_s_NN_GeV=point.get("sqrt_s_NN_GeV"))

    # Gradient ranking uses central/one-sided finite differences only where
    # neighbouring rows are screened.  No interpolation crosses a mask.
    dT = min((right - left for left, right in zip(T_grid[:-1], T_grid[1:])), default=1.0)
    dmu = min((right - left for left, right in zip(muB_grid[:-1], muB_grid[1:])), default=1.0)
    gradients: list[tuple[float, dict[str, Any], float, float]] = []
    for row in rows:
        if not _screened(row):
            continue
        T, muB = row["T_MeV"], row["muB_MeV"]
        neighbours = {
            "T_low": row_map.get((T - dT, muB)),
            "T_high": row_map.get((T + dT, muB)),
            "mu_low": row_map.get((T, muB - dmu)),
            "mu_high": row_map.get((T, muB + dmu)),
        }
        values = {name: neighbour["Kplus_over_pi_plus"] if _screened(neighbour) else None for name, neighbour in neighbours.items()}
        if values["T_low"] is not None and values["T_high"] is not None:
            grad_T = (values["T_high"] - values["T_low"]) / (2.0 * dT)
        elif values["T_high"] is not None:
            grad_T = (values["T_high"] - row["Kplus_over_pi_plus"]) / dT
        elif values["T_low"] is not None:
            grad_T = (row["Kplus_over_pi_plus"] - values["T_low"]) / dT
        else:
            continue
        if values["mu_low"] is not None and values["mu_high"] is not None:
            grad_mu = (values["mu_high"] - values["mu_low"]) / (2.0 * dmu)
        elif values["mu_high"] is not None:
            grad_mu = (values["mu_high"] - row["Kplus_over_pi_plus"]) / dmu
        elif values["mu_low"] is not None:
            grad_mu = (row["Kplus_over_pi_plus"] - values["mu_low"]) / dmu
        else:
            continue
        magnitude = math.hypot(grad_T, grad_mu)
        gradients.append((magnitude, row, grad_T, grad_mu))
    gradients.sort(key=lambda item: item[0], reverse=True)
    for magnitude, row, grad_T, grad_mu in gradients[:max_gradient]:
        add(row, "gradient_extremum", "screening_ratio_finite_difference", gradient_T_per_MeV=grad_T, gradient_muB_per_MeV=grad_mu, gradient_magnitude_per_MeV=magnitude)

    # Select screened sides adjacent to any failed/masked point.  The failed
    # point itself is retained only in the boundary summary, never dispatched.
    boundary_keys: list[tuple[float, float]] = []
    for row in rows:
        if row["status"] == "screened":
            continue
        T, muB = row["T_MeV"], row["muB_MeV"]
        for key in ((T - dT, muB), (T + dT, muB), (T, muB - dmu), (T, muB + dmu)):
            if _screened(row_map.get(key)):
                boundary_keys.append(key)
    for key in sorted(set(boundary_keys)):
        add(row_map[key], "mask_boundary", "screened_neighbor_of_failed_point")

    # Historical phase references are guides only.  Points outside the
    # screening domain are recorded as uncovered and do not create fake rows.
    phase_coverage: dict[str, int] = {}
    for name, line in reference_lines.items():
        if name == "freezeout":
            continue
        covered = 0
        points = line.get("points", [])
        if points:
            stride = max(1, len(points) // max_phase)
            for point in points[::stride][:max_phase]:
                target_T = float(point["T_MeV"])
                target_muB = float(point["muB_MeV"])
                if min(T_grid) <= target_T <= max(T_grid) and min(muB_grid) <= target_muB <= max(muB_grid):
                    covered += 1
                    nearest = _nearest(rows, target_T, target_muB)
                    if nearest is not None:
                        distance = (nearest["T_MeV"] - target_T) ** 2 + (nearest["muB_MeV"] - target_muB) ** 2
                        add(nearest, f"{name}_nearby", line.get("source", name), target_T_MeV=target_T, target_muB_MeV=target_muB, distance_squared_MeV2=distance)
        phase_coverage[name] = covered

    output = []
    for key in sorted(selected):
        row = selected[key]
        row["provenance"] = sorted(provenance[key])
        row["production_point"] = f"{row['T_MeV']:.6f}:{row['muB_MeV']:.6f}"
        output.append(row)
    summary = {
        "grid_count": len(rows),
        "screened_count": sum(row["status"] == "screened" for row in rows),
        "failed_count": sum(row["status"] != "screened" for row in rows),
        "candidate_count": len(output),
        "category_counts": {category: sum(category in row["category"] for row in output) for category in {row["category"] for row in output}},
        "phase_reference_covered_anchor_counts": phase_coverage,
        "dispatch_points": ",".join(row["production_point"] for row in output),
    }
    return output, summary


def write_outputs(output_dir: Path, csv_path: Path, rows: list[dict[str, Any]], summary: dict[str, Any], reference_lines: dict[str, Any], *, source_run_id: str | None = None, source_git_sha: str | None = None) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fields = [
        "T_MeV", "muB_MeV", "production_point", "status", "screening_ratio_Kplus_over_pi_plus",
        "category", "source", "provenance", "distance_squared_MeV2", "target_T_MeV", "target_muB_MeV",
        "sqrt_s_NN_GeV", "gradient_T_per_MeV", "gradient_muB_per_MeV", "gradient_magnitude_per_MeV",
    ]
    with (output_dir / "candidate_points.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            output = dict(row)
            output["provenance"] = ";".join(output.get("provenance", []))
            writer.writerow(output)
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "diagnostic_only": True,
        "production_default": False,
        "full_gate_required": True,
        "solver_free_selection": True,
        "source_screening_csv": str(csv_path),
        "source_screening_sha256": sha256_file(csv_path),
        "source_run_id": source_run_id,
        "source_git_sha": source_git_sha,
        "units": {"T": "MeV", "muB": "MeV", "ratio": "dimensionless"},
        "background_contract": {"scope": "fixed quark-only BQS charged GBU partial yield", "rho_Q_over_rho_B": 0.4, "rho_S_fm3": 0.0, "meson_feedback": False},
        "selection_policy": {
            "freezeout": "nearest screened grid rows to default chemical-freezeout curve",
            "gradient_extremum": "finite differences of Kplus/pi+ using screened neighbours only",
            "mask_boundary": "screened rows adjacent to failed/masked screening points",
            "phase_reference": "nearest screened rows only when reference point lies inside screening domain",
            "failed_points": "never dispatched; never zero-filled",
        },
        "reference_lines": {key: {field: value for field, value in line.items() if field != "points"} | {"point_count": len(line["points"])} for key, line in reference_lines.items()},
        "summary": summary,
        "candidates": rows,
    }
    (output_dir / "candidate_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="merged screening CSV or its containing directory")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--source-run-id", default=None)
    parser.add_argument("--git-sha", default=None)
    parser.add_argument("--freezeout-profile", type=Path, default=ROOT / "config" / "physics" / "freezeout" / "default.toml")
    parser.add_argument("--phase-reference-root", type=Path, default=ROOT / "data" / "reference" / "pnjl" / "issue130_phase_reference_v2" / "accepted")
    parser.add_argument("--reference-xi", type=float, default=0.0)
    parser.add_argument("--max-gradient", type=int, default=6)
    parser.add_argument("--max-freezeout", type=int, default=10)
    parser.add_argument("--max-phase", type=int, default=6)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.max_gradient < 1 or args.max_freezeout < 1 or args.max_phase < 1:
        raise ValueError("candidate limits must be positive")
    csv_path = _input_csv(args.input.resolve())
    rows, T_grid, muB_grid = load_rows(csv_path)
    plotter = _load_plotter()
    reference_lines = plotter.load_reference_lines(
        freezeout_profile=args.freezeout_profile.resolve(),
        phase_reference_root=args.phase_reference_root.resolve(),
        xi=args.reference_xi,
    )
    candidates, summary = select_candidates(rows, T_grid, muB_grid, reference_lines, max_gradient=args.max_gradient, max_freezeout=args.max_freezeout, max_phase=args.max_phase)
    write_outputs(args.output_dir.resolve(), csv_path, candidates, summary, reference_lines, source_run_id=args.source_run_id, source_git_sha=args.git_sha)
    print(json.dumps(summary, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
