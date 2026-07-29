#!/usr/bin/env python3
"""Plot aggregated CEP narrow pilot v2 evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


COLORS = {
    "c2_dense_baseline": "#1f77b4",
    "rho_support_cascade": "#2ca02c",
    "high_resolution_oracle": "#d62728",
}


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _finite(value: str | None) -> bool:
    try:
        return value not in {None, "", "nan", "NaN", "null"} and float(value) == float(value)
    except (TypeError, ValueError):
        return False


def plot(input_dir: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    accuracy = _rows(input_dir / "cep_accuracy.csv")
    costs = _rows(input_dir / "method_costs.csv")
    slices = _rows(input_dir / "slice_metrics.csv")
    curves = _rows(input_dir / "curve_points.csv")

    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    for method in sorted({row["method"] for row in accuracy}):
        subset = sorted((row for row in accuracy if row["method"] == method), key=lambda row: float(row["xi"]))
        x = [float(row["xi"]) for row in subset]
        low = [float(row["fine_T_last_first_order_MeV"]) if _finite(row.get("fine_T_last_first_order_MeV")) else float("nan") for row in subset]
        high = [float(row["fine_T_first_monotone_MeV"]) if _finite(row.get("fine_T_first_monotone_MeV")) else float("nan") for row in subset]
        color = COLORS.get(method)
        ax.plot(x, low, marker="o", linestyle="-", label=f"{method} last FO", color=color)
        ax.plot(x, high, marker="x", linestyle="--", color=color, alpha=0.8, label=f"{method} first monotone")
        for xi, lo, hi in zip(x, low, high):
            if _finite(str(lo)) and _finite(str(hi)):
                ax.vlines(xi, lo, hi, color=color, alpha=0.25)
    ax.set_xlabel("$\\xi$")
    ax.set_ylabel("temperature [MeV]")
    ax.set_title("Three-state CEP evidence interval")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(output_dir / "cep_intervals.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for row in costs:
        ax.scatter(
            float(row.get("runner_seconds", 0.0)),
            float(row.get("unique_solves", 0.0)),
            color=COLORS.get(row["method"]),
            label=row["method"],
        )
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys(), fontsize=8)
    ax.set_xlabel("runner seconds")
    ax.set_ylabel("unique equilibrium solves")
    ax.set_title("Solver cost frontier")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_dir / "cost_frontier.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    status_y = {"confirmed_first_order": 2, "ambiguous_near_critical": 1, "confirmed_monotone": 0}
    for method in sorted({row["method"] for row in slices}):
        subset = [row for row in slices if row["method"] == method]
        grouped: dict[float, list[float]] = {}
        for row in subset:
            status = row.get("result_status", "ambiguous_near_critical")
            grouped.setdefault(float(row["T_MeV"]), []).append(float(status_y.get(status, 1)))
        temps = sorted(grouped)
        ax.plot(temps, [sum(grouped[t]) / len(grouped[t]) for t in temps], marker=".", label=method, color=COLORS.get(method))
    ax.set_yticks([0, 1, 2], ["confirmed monotone", "ambiguous", "confirmed first-order"])
    ax.set_xlabel("$T$ [MeV]")
    ax.set_title("Three-state slice classifications")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output_dir / "slice_states.png", dpi=160)
    plt.close(fig)

    # Produce one compact rho-mu panel per xi from the first/last evidence
    # temperatures.  The full raw curve table remains in the Actions artifact.
    for xi in sorted({float(row["xi"]) for row in accuracy}):
        fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.8), sharey=True)
        for axis, method in zip(axes, sorted(COLORS)):
            method_accuracy = next((row for row in accuracy if float(row["xi"]) == xi and row["method"] == method), None)
            if method_accuracy is None:
                axis.axis("off")
                continue
            temperatures = []
            for key in ("fine_T_last_first_order_MeV", "fine_T_first_monotone_MeV"):
                if _finite(method_accuracy.get(key)):
                    temperatures.append(float(method_accuracy[key]))
            subset = [row for row in curves if float(row["xi"]) == xi and row["method"] == method]
            if temperatures:
                for temp in temperatures:
                    points = [row for row in subset if abs(float(row["T_MeV"]) - temp) <= 1e-8 and row.get("rho_level") == "1"]
                    points.sort(key=lambda row: float(row["rho"]))
                    if points:
                        axis.plot([float(row["rho"]) for row in points], [float(row["muq_MeV"]) for row in points], marker=".", label=f"T={temp:.3f}")
            axis.set_title(method.replace("_", " "))
            axis.set_xlabel("rho")
            axis.grid(alpha=0.25)
            handles, labels = axis.get_legend_handles_labels()
            if handles:
                axis.legend(fontsize=7)
        axes[0].set_ylabel("$\\mu_q$ [MeV]")
        fig.suptitle(f"Representative rho-mu curves, xi={xi:g}")
        fig.tight_layout()
        fig.savefig(output_dir / f"rho_mu_xi_{xi:g}.png", dpi=160)
        plt.close(fig)

    files = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(output_dir.glob("*.png"))}
    inputs = [input_dir / name for name in ("cep_accuracy.csv", "method_costs.csv", "slice_metrics.csv", "curve_points.csv")]
    manifest = {
        "schema_version": "cep_narrow_pilot_v2",
        "generated_at_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "input_sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs if path.is_file()},
        "figures": files,
    }
    (output_dir / "plot_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    plot(args.input_dir.resolve(), args.output_dir.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
