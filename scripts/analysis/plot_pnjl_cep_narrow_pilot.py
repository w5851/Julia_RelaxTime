#!/usr/bin/env python3
"""Plot the diagnostic tables produced by collect_pnjl_cep_narrow_pilot.py."""

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


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def plot(input_dir: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    accuracy = rows(input_dir / "cep_accuracy.csv")
    costs = rows(input_dir / "method_costs.csv")
    slices = rows(input_dir / "slice_metrics.csv")

    colors = {
        "c2_dense_baseline": "#1f77b4",
        "rho_support_cascade": "#2ca02c",
        "high_resolution_oracle": "#d62728",
    }

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for method in sorted({row["method"] for row in accuracy}):
        subset = sorted((row for row in accuracy if row["method"] == method), key=lambda row: float(row["xi"]))
        ax.plot(
            [float(row["xi"]) for row in subset],
            [float(row["estimated_T_CEP_MeV"]) for row in subset],
            marker="o",
            label=method,
            color=colors.get(method),
        )
    ax.set_xlabel("$\\xi$")
    ax.set_ylabel("$T_{CEP}$ [MeV]")
    ax.set_title("CEP comparison")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output_dir / "cep_comparison.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for row in costs:
        ax.scatter(
            float(row["runner_seconds"]),
            float(row["equilibrium_requests"]),
            color=colors.get(row["method"]),
            label=row["method"],
        )
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys(), fontsize=8)
    ax.set_xlabel("runner seconds")
    ax.set_ylabel("equilibrium requests")
    ax.set_title("Cost frontier")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_dir / "cost_frontier.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    for method in sorted({row["method"] for row in slices}):
        subset = sorted((row for row in slices if row["method"] == method), key=lambda row: (float(row["xi"]), float(row["T_MeV"])))
        # Aggregate slice wall time by temperature across xi for a compact view.
        grouped: dict[float, float] = {}
        for row in subset:
            grouped[float(row["T_MeV"])] = grouped.get(float(row["T_MeV"]), 0.0) + float(row["wall_seconds"])
        ax.plot(sorted(grouped), [grouped[key] for key in sorted(grouped)], marker=".", label=method, color=colors.get(method))
    ax.set_xlabel("$T$ [MeV]")
    ax.set_ylabel("summed slice wall seconds")
    ax.set_title("Slice costs")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output_dir / "slice_costs.png", dpi=160)
    plt.close(fig)

    files = {}
    for path in sorted(output_dir.glob("*.png")):
        files[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = {
        "schema_version": "cep_narrow_pilot_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "input_sha256": {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (input_dir / "cep_accuracy.csv", input_dir / "method_costs.csv", input_dir / "slice_metrics.csv")
        },
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
