#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt


PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INPUT = (
    PROJECT_ROOT
    / "data"
    / "outputs"
    / "results"
    / "relaxtime"
    / "meson_density"
    / "crossover_validation"
    / "friesen2019_crossover_line"
    / "crossover_line.csv"
)


def configure_style() -> None:
    matplotlib.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times"],
            "font.size": 9,
            "mathtext.fontset": "stix",
            "axes.labelsize": 9,
            "axes.linewidth": 0.7,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "legend.frameon": False,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.06,
        }
    )


def read_rows(path: Path):
    rows = []
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            try:
                rows.append(
                    {
                        "mu_MeV": float(row["mu_MeV"]),
                        "muB_MeV": float(row["muB_MeV"]),
                        "T_crossover_MeV": float(row["T_crossover_MeV"]),
                        "T_over_muB": float(row["T_over_muB"]),
                        "converged": row["converged"].strip().lower() == "true",
                    }
                )
            except Exception:
                continue
    rows.sort(key=lambda item: item["muB_MeV"])
    return rows


def plot_crossover_line(rows, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    conv = [r for r in rows if r["converged"]]
    fail = [r for r in rows if not r["converged"]]

    if conv:
        ax.plot(
            [r["muB_MeV"] for r in conv],
            [r["T_crossover_MeV"] for r in conv],
            color="#d62728",
            linewidth=1.8,
            marker="o",
            markersize=3.2,
            label="crossover line",
        )
    if fail:
        ax.scatter(
            [r["muB_MeV"] for r in fail],
            [r["T_crossover_MeV"] for r in fail],
            color="#7f7f7f",
            s=16,
            marker="x",
            label="non-converged",
        )

    ax.set_xlabel(r"$\mu_B$ [MeV]")
    ax.set_ylabel(r"$T_c$ [MeV]")
    ax.grid(True, alpha=0.2, linewidth=0.5)
    ax.legend(loc="best")
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_scaling(rows, out_path: Path) -> None:
    conv = [r for r in rows if r["converged"] and r["muB_MeV"] > 0.0]
    fig, ax = plt.subplots(figsize=(5.6, 3.4))
    ax.plot(
        [r["muB_MeV"] for r in conv],
        [r["T_over_muB"] for r in conv],
        color="#1f77b4",
        linewidth=1.8,
        marker="o",
        markersize=3.0,
    )
    ax.set_xlabel(r"$\mu_B$ [MeV]")
    ax.set_ylabel(r"$T/\mu_B$")
    ax.grid(True, alpha=0.2, linewidth=0.5)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def write_summary(rows, out_path: Path) -> None:
    conv = [r for r in rows if r["converged"]]
    lines = [
        "# Friesen 2019 crossover line quicklook",
        "",
        f"- total points: {len(rows)}",
        f"- converged points: {len(conv)}",
    ]
    if conv:
        lines.extend(
            [
                f"- T_c(mu_B=0): {conv[0]['T_crossover_MeV']:.6f} MeV",
                f"- first point: mu_B={conv[0]['muB_MeV']:.6f} MeV, T_c={conv[0]['T_crossover_MeV']:.6f} MeV",
                f"- last point: mu_B={conv[-1]['muB_MeV']:.6f} MeV, T_c={conv[-1]['T_crossover_MeV']:.6f} MeV",
            ]
        )
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Plot temporary quicklook figures for Friesen-2019 crossover line.")
    ap.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="input crossover_line.csv")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    configure_style()
    rows = read_rows(args.input)
    out_dir = args.input.parent
    plot_crossover_line(rows, out_dir / "crossover_line_quicklook.png")
    plot_scaling(rows, out_dir / "crossover_line_scaling_quicklook.png")
    write_summary(rows, out_dir / "quicklook_summary.md")


if __name__ == "__main__":
    main()
