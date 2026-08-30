#!/usr/bin/env python3
"""Plot the solver-free Maxwell candidate/endpoint audit."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def pick(rows: list[dict[str, str]], xi: float, temperature: float) -> list[dict[str, str]]:
    selected = [
        row
        for row in rows
        if abs(float(row["xi"]) - xi) < 1e-9
        and abs(float(row["T_MeV"]) - temperature) < 1e-9
        and row["method"] == "independent_oracle"
    ]
    selected.sort(key=lambda row: float(row["rho"]))
    return selected


def summary_row(rows: list[dict[str, str]], xi: float, temperature: float) -> dict[str, str]:
    for row in rows:
        if abs(float(row["xi"]) - xi) < 1e-9 and abs(float(row["T_MeV"]) - temperature) < 1e-9 and row["method"] == "independent_oracle":
            return row
    raise ValueError(f"missing summary row xi={xi}, T={temperature}")


def plot_anchor(rows: list[dict[str, str]], summary: dict[str, str], output: Path, *, near_zero: bool) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    x = [float(row["rho"]) for row in rows]
    y = [float(row["muq_MeV"]) for row in rows]
    axes[0].plot(x, y, color="#1f77b4", linewidth=1.2, marker=".", markersize=2)
    axes[0].set_xlabel(r"$\rho$")
    axes[0].set_ylabel(r"$\mu_q$ [MeV]")
    axes[0].set_title("full curve")
    axes[0].grid(alpha=0.25)
    axes[0].set_xlim(left=0.0)

    if near_zero:
        local = [(rho, mu) for rho, mu in zip(x, y) if rho <= 0.02]
    else:
        local = list(zip(x, y))
    axes[1].plot([item[0] for item in local], [item[1] for item in local],
                 color="#d62728", linewidth=1.4, marker="o", markersize=3)
    candidate_mu = float(summary["candidate_mu"])
    rho_h = float(summary["rho_hadron"])
    rho_q = float(summary["rho_quark"])
    axes[1].axhline(candidate_mu, color="#2ca02c", linestyle="--", linewidth=1.0,
                    label=f"candidate $\\mu_M$={candidate_mu:.6f}")
    axes[1].axvline(rho_h, color="#9467bd", linestyle=":", linewidth=1.0,
                    label=f"$\\rho_h$={rho_h:.6g}")
    axes[1].axvline(rho_q, color="#8c564b", linestyle=":", linewidth=1.0,
                    label=f"$\\rho_q$={rho_q:.6g}")
    axes[1].set_xlabel(r"$\rho$")
    axes[1].set_ylabel(r"$\mu_q$ [MeV]")
    axes[1].set_title("near-rho=0 / candidate crossings")
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8, loc="best")
    axes[1].set_xlim(0.0, 0.02 if near_zero else max(x))
    fig.suptitle(f"independent oracle: xi={summary['xi']}, T={summary['T_MeV']} MeV")
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = read_rows(args.input_dir / "curve_points.csv")
    summaries = read_rows(args.summary)
    figure_paths = []
    for xi, temperature, slug in [
        (-0.5, 5.0, "minus05_T5"),
        (-0.5, 20.0, "minus05_T20"),
        (0.0, 5.0, "zero_T5"),
    ]:
        anchor = pick(rows, xi, temperature)
        summary = summary_row(summaries, xi, temperature)
        output = args.output_dir / f"{slug}_full_near_zero.png"
        plot_anchor(anchor, summary, output, near_zero=True)
        figure_paths.append(output)
    def sha256(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = {
        "schema_version": "pnjl_maxwell_endpoint_candidate_feasibility_v1",
        "source_curve_file": str((args.input_dir / "curve_points.csv").resolve()),
        "source_curve_sha256": sha256(args.input_dir / "curve_points.csv"),
        "summary_file": str(args.summary.resolve()),
        "summary_sha256": sha256(args.summary),
        "figures": [
            {"path": str(path.relative_to(args.output_dir)).replace("\\", "/"),
             "sha256": sha256(path),
             "scope": "full_curve_plus_rho_le_0.02_zoom"}
            for path in figure_paths
        ],
        "command": "python scripts/analysis/plot_pnjl_maxwell_endpoint_candidate.py",
    }
    (args.output_dir / "plot_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
