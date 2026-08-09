#!/usr/bin/env python3
"""Render representative plots for the Stage-C density feasibility v2 replay.

The plotter is intentionally solver-free.  It reads the downloaded numerical
CSV files and the Julia replay tables, then writes only figures and plot
metadata into the evidence directory.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


CASES = (
    (-0.35, 51.0, "density_first_order"),
    (-0.5, 60.0, "first_order_control"),
    (-0.5, 160.0, "monotone_control"),
    (0.35, 101.0, "density_first_order"),
)


def _read(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _read_recursive(root: Path, name: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in sorted(root.rglob(name)):
        rows.extend(_read(path))
    return rows


def _same(value: str | float, target: float, tol: float = 2e-4) -> bool:
    try:
        return math.isclose(float(value), target, abs_tol=tol, rel_tol=0.0)
    except (TypeError, ValueError):
        return False


def _token(value: float) -> str:
    return str(value).replace("-", "m").replace(".", "p")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _finite(values: list[float]) -> list[float]:
    return [value for value in values if math.isfinite(value)]


def _local_x_window(x: list[float], y: list[float], transition: float) -> tuple[float, float]:
    if len(x) < 3:
        return min(x), max(x)
    slopes = [y[index + 1] - y[index] for index in range(len(y) - 1)]
    negative = [index for index, value in enumerate(slopes) if value < 0.0]
    if negative:
        left = max(0, negative[0] - 8)
        right = min(len(x) - 1, negative[-1] + 9)
        return x[left], x[right]
    if math.isfinite(transition):
        near = [index for index, value in enumerate(y) if abs(value - transition) <= 0.2 * max(max(y) - min(y), 1e-12)]
        if near:
            return x[max(0, near[0] - 12)], x[min(len(x) - 1, near[-1] + 12)]
    return x[max(0, len(x) // 5)], x[min(len(x) - 1, 4 * len(x) // 5)]


def _local_y_window(y: list[float], x: list[float], x_window: tuple[float, float], role: str,
                    transition: float) -> tuple[float, float]:
    values = [value for value, rho in zip(y, x) if x_window[0] <= rho <= x_window[1]]
    finite = _finite(values)
    if not finite:
        finite = _finite(y)
    low, high = min(finite), max(finite)
    span = max(high - low, 1e-8)
    if "monotone" in role:
        # Keep the middle, slowly varying part visible without pretending it
        # has a turning point that the full curve does not contain.
        center = 0.5 * (low + high)
        half = max(0.22 * span, 1e-6)
        return center - half, center + half
    if math.isfinite(transition):
        half = max(0.34 * span, 0.004)
        return transition - half, transition + half
    return low - 0.08 * span, high + 0.08 * span


def render(input_dir: Path, output_dir: Path, replay_dir: Path) -> dict[str, object]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    curves = _read_recursive(input_dir, "curve_points.csv")
    slices = _read_recursive(input_dir, "slice_metrics.csv")
    selected = _read(replay_dir / "selected_point_index.csv")
    route_rows = _read(replay_dir / "route_comparison.csv")
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    figures: list[dict[str, str]] = []

    for xi, temperature, role in CASES:
        points = [
            row for row in curves
            if row.get("method") == "independent_oracle"
            and _same(row.get("xi", "nan"), xi, 1e-8)
            and _same(row.get("T_MeV", "nan"), temperature)
            and row.get("rho_level") == "0"
        ]
        points.sort(key=lambda row: float(row["rho"]))
        if len(points) < 3:
            continue
        x = [float(row["rho"]) for row in points]
        y = [float(row["muq_MeV"]) for row in points]
        transition = math.nan
        for row in route_rows:
            if (row.get("route") == "stage_b_features_v1" and row.get("cap") == "12"
                    and _same(row.get("xi", "nan"), xi, 1e-8)
                    and _same(row.get("T_MeV", "nan"), temperature)):
                try:
                    transition = float(row.get("maxwell_area_gate", "nan"))
                except (TypeError, ValueError):
                    pass
                break
        # route_comparison stores the area gate, not mu_transition.  Prefer
        # the source slice's transition field for the plotted horizontal line.
        for row in slices:
            if (row.get("method") == "memoized_dense"
                    and _same(row.get("xi", "nan"), xi, 1e-8)
                    and _same(row.get("T_MeV", "nan"), temperature)):
                try:
                    transition = float(row.get("mu_transition_MeV", "nan"))
                except (TypeError, ValueError):
                    transition = math.nan
                break

        selected_points = [
            (float(row["rho"]), float(row["muq_MeV"]))
            for row in selected
            if row.get("route") == "stage_b_features_v1" and row.get("cap") == "12"
            and _same(row.get("xi", "nan"), xi, 1e-8)
            and _same(row.get("T_MeV", "nan"), temperature)
        ]
        zoom_x = _local_x_window(x, y, transition)
        zoom_y = _local_y_window(y, x, zoom_x, role, transition)

        fig, (full, zoom) = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
        full.plot(x, y, color="#264653", linewidth=1.15, label="independent-oracle Stage-B")
        zoom.plot(x, y, color="#264653", linewidth=1.5)
        if selected_points:
            sx, sy = zip(*selected_points)
            full.scatter(sx, sy, s=14, color="#e76f51", label="selected Stage-C points", zorder=3)
            zoom.scatter(sx, sy, s=26, color="#e76f51", zorder=3)
        if math.isfinite(transition):
            full.axhline(transition, color="#2a9d8f", linestyle="--", linewidth=0.9, label="Maxwell mu")
            zoom.axhline(transition, color="#2a9d8f", linestyle="--", linewidth=0.9)
        full.set_xlabel(r"$\rho$")
        full.set_ylabel(r"$\mu_q$ [MeV]")
        full.set_title(f"xi={xi:g}, T={temperature:g} MeV")
        full.legend(fontsize=8, loc="best")
        full.grid(alpha=0.2)
        zoom.set_xlim(*zoom_x)
        zoom.set_ylim(*zoom_y)
        zoom.set_xlabel(r"$\rho$")
        zoom.set_ylabel(r"$\mu_q$ [MeV]")
        zoom.set_title("local S/Maxwell view" if "monotone" not in role else "local monotone view")
        zoom.grid(alpha=0.2)
        name = f"rho_mu_xi_{_token(xi)}_T_{_token(temperature)}.png"
        path = figure_dir / name
        fig.savefig(path, dpi=160)
        plt.close(fig)
        figures.append({"path": f"figures/{name}", "xi": str(xi),
                        "T_MeV": str(temperature), "role": role})

    plot_manifest = {
        "schema_version": "cep_hybrid_stagec_density_certificate_feasibility_v2",
        "figures": figures,
        "input_curve_files": sorted({str(path.relative_to(input_dir).as_posix())
                                      for path in input_dir.rglob("curve_points.csv")}),
        "input_file_hashes": {
            str(path.relative_to(input_dir).as_posix()): _sha256(path)
            for path in sorted(input_dir.rglob("*.csv"))
        },
    }
    (output_dir / "plot_manifest.json").write_text(
        json.dumps(plot_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    aggregate_manifest_path = output_dir / "manifest.json"
    if aggregate_manifest_path.is_file():
        aggregate = json.loads(aggregate_manifest_path.read_text(encoding="utf-8"))
        hashes: dict[str, str] = {}
        for path in sorted(output_dir.rglob("*")):
            if path.is_file() and path.name != "manifest.json":
                hashes[path.relative_to(output_dir).as_posix()] = _sha256(path)
        aggregate["files"] = hashes
        aggregate["plot_manifest_sha256"] = _sha256(output_dir / "plot_manifest.json")
        aggregate_manifest_path.write_text(json.dumps(aggregate, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return plot_manifest


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--replay-dir", type=Path, default=None)
    args = parser.parse_args()
    result = render(args.input_dir, args.output_dir, args.replay_dir or args.output_dir)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
