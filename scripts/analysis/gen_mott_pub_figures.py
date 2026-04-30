#!/usr/bin/env python3
"""Generate publication-quality figures for the Mott paper.

Generates two core figures (pion and kaon Mott crossing) with:
  - 3 xi values × 3 line types = 9 curves per panel
  - M_mes (solid), M_quarksum (dashed), Gamma_mes (dotted)
  - Mott crossing: single dot only, no vertical line, no annotation text
  - Units: MeV (converted from fm^{-1} via × 197.327)
  - 3×3 legend (each of 9 lines listed individually)
  - Y-axis fixed 0–1000 MeV; minor ticks enabled
  - Output: PDF (vector, pdflatex-native) + PNG (raster)

Usage:
  python scripts/analysis/gen_mott_pub_figures.py \
      --csv data/outputs/results/relaxtime/mott_phase/reference_100_300_fine/mott_phase_derived.csv \
      --out-dir paper_lib/My_Paper/pnjl_meson_mott_xi/figures
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
from cycler import cycler

FM_TO_MEV = 197.327
XI_VALUES = [-0.3, 0.0, 0.3]
APS_COLORS = ["#4477AA", "#EE6677", "#228833"]
XI_COLOR = { -0.3: APS_COLORS[0], 0.0: APS_COLORS[1], 0.3: APS_COLORS[2] }
XI_LS    = { -0.3: "--", 0.0: "-", 0.3: "-." }
XI_LABEL = { -0.3: r"$\xi=-0.3$", 0.0: r"$\xi=0$", 0.3: r"$\xi=+0.3$" }
# Line-style descriptors for the three quantity types
TYPE_LS   = {"mes": (0, (1.5, 1.0)), "thr": "-",  "gam": (0, (3.0, 2.0))}
TYPE_LW   = {"mes": 1.0,             "thr": 1.1,  "gam": 0.9}


def configure_style():
    matplotlib.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times"],
        "font.size": 9,
        "mathtext.fontset": "stix",
        "axes.prop_cycle": cycler(color=APS_COLORS),
        "axes.labelsize": 9,
        "axes.linewidth": 0.6,
        "legend.fontsize": 6.5,
        "legend.frameon": False,
        "legend.handlelength": 1.5,
        "legend.handletextpad": 0.5,
        "legend.borderpad": 0.2,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.minor.width": 0.35,
        "ytick.minor.width": 0.35,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "lines.linewidth": 1.0,
        "lines.markersize": 3.5,
        "axes.grid": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
    })


def read_csv(path: Path) -> List[Dict[str, str]]:
    rows = []
    with path.open("r", encoding="utf-8") as f:
        header = None
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            header = [x.strip() for x in s.split(",")]
            break
        if header is None:
            return rows
        reader = csv.DictReader(f, fieldnames=header)
        for row in reader:
            rows.append(row)
    return rows


def estimate_mott_temperature(
    T_vals: np.ndarray, M_mes: np.ndarray, M_thr: np.ndarray
) -> Optional[float]:
    gap = M_mes - M_thr
    for i in range(len(gap) - 1):
        if gap[i] * gap[i + 1] <= 0:
            if abs(gap[i + 1] - gap[i]) < 1e-15:
                return float(T_vals[i])
            t = -gap[i] / (gap[i + 1] - gap[i])
            return float(T_vals[i] + t * (T_vals[i + 1] - T_vals[i]))
    return None


def plot_mott_figure(
    rows: List[Dict[str, str]],
    meson_col: str,
    thr_col: str,
    gamma_col: str,
    meson_label: str,
    thr_label: str,
    gamma_label: str,
    out_path: Path,
):
    # ── Collect data ──────────────────────────────────────────────────
    data: Dict[float, Dict[str, list]] = {}
    for r in rows:
        try:
            xi = float(r["xi"])
            T  = float(r["T_MeV"])
            M  = float(r[meson_col])  * FM_TO_MEV
            thr = float(r[thr_col])   * FM_TO_MEV
            gam = float(r[gamma_col]) * FM_TO_MEV
        except (KeyError, ValueError):
            continue
        if xi not in XI_VALUES:
            continue
        if xi not in data:
            data[xi] = {"T": [], "M": [], "thr": [], "gam": []}
        data[xi]["T"].append(T)
        data[xi]["M"].append(M)
        data[xi]["thr"].append(thr)
        data[xi]["gam"].append(gam)

    for xi in data:
        order = np.argsort(data[xi]["T"])
        for k in data[xi]:
            data[xi][k] = np.array(data[xi][k])[order]

    # ── Figure ───────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(3.375, 2.6))

    mott_points = []
    legend_handles = []
    for xi in sorted(data.keys()):
        d = data[xi]
        color = XI_COLOR[xi]

        # Meson mass (solid)
        ax.plot(d["T"], d["M"], color=color, ls=TYPE_LS["mes"], lw=TYPE_LW["mes"])
        legend_handles.append(
            mlines.Line2D([], [], color=color, ls=TYPE_LS["mes"], lw=TYPE_LW["mes"],
                          label=f"{XI_LABEL[xi]}  {meson_label}")
        )

        # Threshold (solid — same style, distinguished by color)
        ax.plot(d["T"], d["thr"], color=color, ls=TYPE_LS["thr"], lw=TYPE_LW["thr"])
        legend_handles.append(
            mlines.Line2D([], [], color=color, ls=TYPE_LS["thr"], lw=TYPE_LW["thr"],
                          label=f"{XI_LABEL[xi]}  {thr_label}")
        )

        # Width (dotted)
        ax.plot(d["T"], d["gam"], color=color, ls=TYPE_LS["gam"], lw=TYPE_LW["gam"])
        legend_handles.append(
            mlines.Line2D([], [], color=color, ls=TYPE_LS["gam"], lw=TYPE_LW["gam"],
                          label=f"{XI_LABEL[xi]}  {gamma_label}")
        )

        # Mott crossing dot
        T_mott = estimate_mott_temperature(d["T"], d["M"], d["thr"])
        if T_mott is not None:
            M_at_mott = np.interp(T_mott, d["T"], d["M"])
            ax.plot(T_mott, M_at_mott, "o", color=color, ms=4.5, zorder=10,
                    mec="white", mew=0.5)

    # Legend: 3×3 entries, upper right, 2 columns
    ax.legend(
        handles=legend_handles, loc="upper right", fontsize=6,
        ncol=2, columnspacing=0.8, handlelength=1.8, handletextpad=0.4,
    )

    # ── Axes ─────────────────────────────────────────────────────────
    ax.set_xlabel(r"$T$ (MeV)")
    ax.set_ylabel(r"Mass (MeV)")
    ax.set_xlim(100, 280)
    ax.set_ylim(0, 1000)
    ax.minorticks_on()
    ax.tick_params(which="minor", length=2, direction="in", right=True, top=True)

    # ── Save ─────────────────────────────────────────────────────────
    for fmt in ["pdf", "png", "eps"]:
        out = out_path.with_suffix(f".{fmt}")
        if fmt == "eps":
            # EPS backend needs Type 3 outlined fonts (TrueType embedding unreliable)
            matplotlib.rcParams["ps.fonttype"] = 3
            fig.savefig(out, format="eps", dpi=600)
            matplotlib.rcParams["ps.fonttype"] = 42  # restore
        elif fmt == "pdf":
            fig.savefig(out, format="pdf")
        else:
            fig.savefig(out, format="png", dpi=600)
        print(f"Saved: {out}")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Publication-quality Mott figures")
    ap.add_argument("--csv", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    if not args.csv.exists():
        print(f"CSV not found: {args.csv}")
        return 1

    configure_style()
    rows = read_csv(args.csv)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Pion
    plot_mott_figure(
        rows,
        meson_col="M_pi", thr_col="M_u_plus_M_d", gamma_col="Gamma_pi",
        meson_label=r"$M_{\pi}$",
        thr_label=r"$M_u+M_d$",
        gamma_label=r"$\Gamma_{\pi}$",
        out_path=args.out_dir / "fig_mott_pi",
    )

    # Kaon
    plot_mott_figure(
        rows,
        meson_col="M_K", thr_col="M_u_plus_M_s", gamma_col="Gamma_K",
        meson_label=r"$M_K$",
        thr_label=r"$M_u+M_s$",
        gamma_label=r"$\Gamma_K$",
        out_path=args.out_dir / "fig_mott_K",
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
