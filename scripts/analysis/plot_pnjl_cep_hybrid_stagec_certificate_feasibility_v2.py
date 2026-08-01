#!/usr/bin/env python3
"""Render compact representative plots for Stage-C certificate replay v2."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


CASES = (
    (-0.5, 147.0947265625, "author_first_order"),
    (0.5, 106.9599609375, "author_first_order"),
    (-0.5, 147.2197265625, "consensus_monotone"),
    (0.5, 107.0849609375, "consensus_monotone"),
)


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _same(value: float, target: float, tol: float = 2e-4) -> bool:
    return math.isclose(value, target, abs_tol=tol, rel_tol=0.0)


def _token(value: float) -> str:
    return str(value).replace("-", "m").replace(".", "p")


def render(input_dir: Path, output_dir: Path, replay_dir: Path) -> dict[str, object]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    curves = _read(input_dir / "curve_points.csv")
    slices = _read(input_dir / "slice_metrics.csv")
    selected = _read(replay_dir / "selected_points.csv")
    figures: list[dict[str, str]] = []
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)

    for xi, temperature, label in CASES:
        points: dict[float, float] = {}
        for row in curves:
            if row.get("method") != "independent_oracle":
                continue
            if not _same(float(row["xi"]), xi, 1e-8) or not _same(float(row["T_MeV"]), temperature):
                continue
            points[float(row["rho"])] = float(row["muq_MeV"])
        if not points:
            continue
        selected_points = [
            (float(row["rho"]), float(row["muq_MeV"]))
            for row in selected
            if _same(float(row["xi"]), xi, 1e-8)
            and _same(float(row["T_MeV"]), temperature)
        ]
        transition = math.nan
        for row in slices:
            if row.get("method") == "independent_oracle" and _same(float(row["xi"]), xi, 1e-8) and _same(float(row["T_MeV"]), temperature):
                try:
                    transition = float(row.get("mu_transition_MeV", "nan"))
                except ValueError:
                    pass
                break

        x = sorted(points)
        y = [points[value] for value in x]
        fig, (ax, zoom) = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
        ax.plot(x, y, color="#264653", linewidth=1.2, label="independent oracle curve")
        zoom.plot(x, y, color="#264653", linewidth=1.5)
        if selected_points:
            sx, sy = zip(*selected_points)
            ax.scatter(sx, sy, s=12, color="#e76f51", label="selected Stage-C points", zorder=3)
            zoom.scatter(sx, sy, s=24, color="#e76f51", label="selected Stage-C points", zorder=3)
        if math.isfinite(transition):
            ax.axhline(transition, color="#2a9d8f", linestyle="--", linewidth=0.9, label="Maxwell μ")
            zoom.axhline(transition, color="#2a9d8f", linestyle="--", linewidth=0.9)
        ax.set_xlabel(r"$\rho$")
        ax.set_ylabel(r"$\mu_q$ [MeV]")
        ax.set_title(f"xi={xi:g}, T={temperature:g} MeV")
        ax.legend(fontsize=8, loc="best")
        ax.grid(alpha=0.2)

        if label == "author_first_order":
            finite_y = [value for value in y if math.isfinite(value)]
            low, high = min(finite_y), max(finite_y)
            span = max(high - low, 1e-6)
            center = (low + high) / 2
            zoom.set_ylim(center - 0.22 * span, center + 0.22 * span)
            zoom.set_title("local S/Maxwell support")
        else:
            zoom.set_title("monotone control zoom")
        zoom.set_xlabel(r"$\rho$")
        zoom.set_ylabel(r"$\mu_q$ [MeV]")
        zoom.grid(alpha=0.2)
        name = f"rho_mu_xi_{_token(xi)}_T_{_token(temperature)}.png"
        path = figure_dir / name
        fig.savefig(path, dpi=160)
        plt.close(fig)
        figures.append({"path": f"figures/{name}", "xi": str(xi), "T_MeV": str(temperature), "role": label})

    manifest = {"schema_version": "cep_hybrid_stagec_certificate_feasibility_v2", "figures": figures}
    plot_manifest_path = output_dir / "plot_manifest.json"
    plot_manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    aggregate_manifest_path = output_dir / "manifest.json"
    if aggregate_manifest_path.is_file():
        aggregate = json.loads(aggregate_manifest_path.read_text(encoding="utf-8"))
        hashes: dict[str, str] = {}
        for path in sorted(output_dir.rglob("*")):
            if path.is_file() and path.name != "manifest.json":
                hashes[path.relative_to(output_dir).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
        aggregate["files"] = hashes
        aggregate_manifest_path.write_text(json.dumps(aggregate, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    result = render(args.input_dir, args.output_dir, args.output_dir)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
