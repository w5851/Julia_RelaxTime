#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List

import matplotlib
import matplotlib.pyplot as plt


PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT_DIR = (
    PROJECT_ROOT
    / "data"
    / "outputs"
    / "results"
    / "relaxtime"
    / "analysis"
    / "meson_density_phase_f_minimal"
)

INPUTS = {
    "stable": PROJECT_ROOT / "data" / "outputs" / "results" / "relaxtime" / "scan" / "meson_density_scan_208_220_step2.csv",
    "strict_bw": PROJECT_ROOT / "data" / "outputs" / "results" / "relaxtime" / "scan" / "strict_bw_meson_density_scan_stage2_208_220_step2_converged.csv",
    "current_bu": PROJECT_ROOT / "data" / "outputs" / "results" / "relaxtime" / "scan" / "phase_shift_meson_density_scan_208_220_step2.csv",
    "gbu_reference": PROJECT_ROOT / "data" / "outputs" / "results" / "relaxtime" / "scan" / "phase_shift_meson_density_scan_gbu_reference_208_220_step2.csv",
}

STYLE = {
    "stable": {"label": "stable", "color": "#1f77b4", "linestyle": "-", "linewidth": 1.8},
    "strict_bw": {"label": "strict BW", "color": "#d62728", "linestyle": "--", "linewidth": 1.8},
    "current_bu": {"label": "current BU", "color": "#2ca02c", "linestyle": "-.", "linewidth": 1.8},
    "gbu_reference": {"label": "gbu reference", "color": "#9467bd", "linestyle": ":", "linewidth": 2.2},
}


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
            "legend.fontsize": 8,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.06,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def read_scan_csv(path: Path) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    with path.open("r", encoding="utf-8") as fh:
        header = None
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            header = [item.strip() for item in s.split(",")]
            break
        if header is None:
            raise ValueError(f"no header found in {path}")
        reader = csv.DictReader(fh, fieldnames=header)
        for row in reader:
            rows.append(row)
    return rows


def load_series() -> Dict[str, List[Dict[str, float]]]:
    loaded: Dict[str, List[Dict[str, float]]] = {}
    for key, path in INPUTS.items():
        rows = read_scan_csv(path)
        parsed = []
        for row in rows:
            parsed.append(
                {
                    "T_MeV": float(row["T_MeV"]),
                    "n_pi": float(row["n_pi"]),
                    "n_K": float(row["n_K"]),
                    "kpi_ratio": float(row["kpi_ratio"]),
                }
            )
        parsed.sort(key=lambda item: item["T_MeV"])
        loaded[key] = parsed
    return loaded


def write_combined_csv(series: Dict[str, List[Dict[str, float]]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["scheme", "T_MeV", "n_pi", "n_K", "kpi_ratio"])
        for scheme, rows in series.items():
            for row in rows:
                writer.writerow([scheme, row["T_MeV"], row["n_pi"], row["n_K"], row["kpi_ratio"]])


def plot_density_panels(series: Dict[str, List[Dict[str, float]]], out_path: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(6.4, 5.4), sharex=True)
    for scheme, rows in series.items():
        style = STYLE[scheme]
        T = [row["T_MeV"] for row in rows]
        n_pi = [row["n_pi"] for row in rows]
        n_k = [row["n_K"] for row in rows]
        axes[0].plot(T, n_pi, label=style["label"], color=style["color"], linestyle=style["linestyle"], linewidth=style["linewidth"])
        axes[1].plot(T, n_k, label=style["label"], color=style["color"], linestyle=style["linestyle"], linewidth=style["linewidth"])

    axes[0].set_ylabel(r"$n_{\pi}$ [fm$^{-3}$]")
    axes[1].set_ylabel(r"$n_{K}$ [fm$^{-3}$]")
    axes[1].set_xlabel(r"$T$ [MeV]")
    axes[0].legend(loc="upper left", ncol=2)
    for ax in axes:
        ax.grid(True, alpha=0.2, linewidth=0.5)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_ratio(series: Dict[str, List[Dict[str, float]]], out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    for scheme, rows in series.items():
        style = STYLE[scheme]
        T = [row["T_MeV"] for row in rows]
        ratio = [row["kpi_ratio"] for row in rows]
        ax.plot(T, ratio, label=style["label"], color=style["color"], linestyle=style["linestyle"], linewidth=style["linewidth"])

    ax.set_xlabel(r"$T$ [MeV]")
    ax.set_ylabel(r"$K/\pi$")
    ax.grid(True, alpha=0.2, linewidth=0.5)
    ax.legend(loc="best", ncol=2)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def write_summary(series: Dict[str, List[Dict[str, float]]], out_path: Path) -> None:
    current_first = series["current_bu"][0]
    current_last = series["current_bu"][-1]
    gbu_first = series["gbu_reference"][0]
    gbu_last = series["gbu_reference"][-1]
    stable_first = series["stable"][0]
    stable_last = series["stable"][-1]
    bw_first = series["strict_bw"][0]
    bw_last = series["strict_bw"][-1]

    text = f"""# Phase F minimal summary

Input window: T = {stable_first["T_MeV"]:.0f}:{series["stable"][1]["T_MeV"] - stable_first["T_MeV"]:.0f}:{stable_last["T_MeV"]:.0f} MeV

Artifacts:

- density figure: `phase_f_minimal_n_densities.png`
- ratio figure: `phase_f_minimal_kpi_ratio.png`
- merged data: `phase_f_minimal_curves.csv`

Quick read:

1. `stable` remains the lowest-complexity baseline, with `K/pi` rising from about `0.57` to `0.64`.
2. `strict BW Stage2` is not close to `stable` in the current window; it enhances `n_K` strongly and keeps `K/pi` above `1`.
3. `current BU` stays much larger than `stable` in both `n_pi` and `n_K`, with `K/pi` around `1.37-1.40`.
4. `gbu reference` pulls the BU result back toward the `stable` scale and sits well below both `current BU` and `strict BW Stage2` in `K/pi`.

Endpoint values:

- {stable_first["T_MeV"]:.0f} MeV:
  - stable: n_pi={stable_first["n_pi"]:.5f}, n_K={stable_first["n_K"]:.5f}, K/pi={stable_first["kpi_ratio"]:.5f}
  - strict BW: n_pi={bw_first["n_pi"]:.5f}, n_K={bw_first["n_K"]:.5f}, K/pi={bw_first["kpi_ratio"]:.5f}
  - current BU: n_pi={current_first["n_pi"]:.5f}, n_K={current_first["n_K"]:.5f}, K/pi={current_first["kpi_ratio"]:.5f}
  - gbu reference: n_pi={gbu_first["n_pi"]:.5f}, n_K={gbu_first["n_K"]:.5f}, K/pi={gbu_first["kpi_ratio"]:.5f}

- {stable_last["T_MeV"]:.0f} MeV:
  - stable: n_pi={stable_last["n_pi"]:.5f}, n_K={stable_last["n_K"]:.5f}, K/pi={stable_last["kpi_ratio"]:.5f}
  - strict BW: n_pi={bw_last["n_pi"]:.5f}, n_K={bw_last["n_K"]:.5f}, K/pi={bw_last["kpi_ratio"]:.5f}
  - current BU: n_pi={current_last["n_pi"]:.5f}, n_K={current_last["n_K"]:.5f}, K/pi={current_last["kpi_ratio"]:.5f}
  - gbu reference: n_pi={gbu_last["n_pi"]:.5f}, n_K={gbu_last["n_K"]:.5f}, K/pi={gbu_last["kpi_ratio"]:.5f}
"""
    out_path.write_text(text, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Generate minimal Phase F meson-density figures from four workflow scan CSVs.")
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR, help="Output directory for figures and merged csv")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    configure_style()
    series = load_series()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_combined_csv(series, args.out_dir / "phase_f_minimal_curves.csv")
    plot_density_panels(series, args.out_dir / "phase_f_minimal_n_densities.png")
    plot_ratio(series, args.out_dir / "phase_f_minimal_kpi_ratio.png")
    write_summary(series, args.out_dir / "phase_f_minimal_summary.md")


if __name__ == "__main__":
    main()
