#!/usr/bin/env python3
"""Plot upstream mechanism figure for detM_im in T200 positive window."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
from cycler import cycler


APS_COLOR_CYCLE = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB"]


def configure_style(dpi: int) -> None:
    matplotlib.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times"],
            "font.size": 10,
            "mathtext.fontset": "stix",
            "axes.prop_cycle": cycler(color=APS_COLOR_CYCLE),
            "axes.labelsize": 10,
            "axes.linewidth": 0.6,
            "legend.fontsize": 8,
            "legend.frameon": False,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "axes.grid": False,
            "savefig.dpi": dpi,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.08,
        }
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path(r"D:\Desktop\Temp\relaxtime_t200_window\t200_detm_im_upstream_trace.csv"),
    )
    parser.add_argument("--out-dir", type=Path, default=Path("docs/analysis/relaxtime"))
    parser.add_argument("--dpi", type=int, default=600)
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for r in reader:
            rows.append(
                {
                    "xi": float(r["xi"]),
                    "ds": float(r["ds"]),
                    "pi_uu_im": float(r["pi_uu_im"]),
                    "pi_ss_im": float(r["pi_ss_im"]),
                    "im_M00M88": float(r["im_M00M88"]),
                    "im_M08sq": float(r["im_M08sq"]),
                    "detM_im": float(r["detM_im"]),
                }
            )
    return rows


def shade_segments(ax: plt.Axes) -> None:
    ax.axvspan(1.0e-3, 5.0, color="#DDEAF7", alpha=0.25)
    ax.axvspan(5.0, 10.0, color="#FBE7D5", alpha=0.25)
    ax.axvspan(10.0, 20.0, color="#E2F3E8", alpha=0.25)


def main() -> None:
    args = parse_args()
    configure_style(args.dpi)
    rows = read_rows(args.csv)

    xi_values = sorted({r["xi"] for r in rows})

    fig, axes = plt.subplots(3, 1, figsize=(8.4, 9.8), constrained_layout=True)

    for ax in axes:
        ax.set_xscale("log")
        shade_segments(ax)

    style_map = {
        0.34: {"ls": "-", "mk": "o", "me": 1},
        0.36: {"ls": "--", "mk": "s", "me": 2},
        0.38: {"ls": "-.", "mk": "^", "me": 3},
    }

    for xi in xi_values:
        sub = sorted([r for r in rows if abs(r["xi"] - xi) < 1.0e-12], key=lambda r: r["ds"])
        x = [r["ds"] for r in sub]
        st = style_map.get(xi, {"ls": "-", "mk": "o", "me": 1})

        axes[0].plot(
            x,
            [r["pi_uu_im"] for r in sub],
            linewidth=1.4,
            linestyle=st["ls"],
            marker=st["mk"],
            markersize=3.0,
            markevery=st["me"],
            markerfacecolor="white",
            markeredgewidth=0.9,
            label=f"xi={xi:.2f}",
        )
        axes[1].plot(
            x,
            [r["pi_ss_im"] for r in sub],
            linewidth=1.4,
            linestyle=st["ls"],
            marker=st["mk"],
            markersize=3.0,
            markevery=st["me"],
            markerfacecolor="white",
            markeredgewidth=0.9,
            label=f"xi={xi:.2f}",
        )

        axes[2].plot(x, [r["im_M00M88"] for r in sub], linewidth=1.25, marker="o", markersize=2.6, label=f"xi={xi:.2f} Im(M00*M88)")
        axes[2].plot(x, [r["im_M08sq"] for r in sub], "--", linewidth=1.2, marker="s", markersize=2.5, label=f"xi={xi:.2f} Im(M08^2)")
        axes[2].plot(x, [r["detM_im"] for r in sub], ":", linewidth=1.7, marker="^", markersize=2.6, label=f"xi={xi:.2f} detM_im")

    axes[0].set_ylabel("Im(Pi_uu^P)")
    axes[1].set_ylabel("Im(Pi_ss^P)")
    axes[2].set_ylabel("imaginary decomposition")
    axes[2].set_xlabel("Delta s")

    axes[1].set_yscale("symlog", linthresh=1.0e-8)
    axes[1].set_ylim(-1.0e-10, 2.0)
    axes[2].set_yscale("symlog", linthresh=2.0e-5)
    axes[2].axhline(0.0, color="#666666", linewidth=0.7)

    axes[0].legend(ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
    axes[1].legend(ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
    axes[2].legend(ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)

    axes[2].text(1.5e-3, 2.3e-2, "decrease", fontsize=8)
    axes[2].text(5.6, 2.3e-2, "recover & cross", fontsize=8)
    axes[2].text(11.5, 2.3e-2, "grow", fontsize=8)
    axes[0].text(1.2e-3, 2.2, "three xi almost overlap", fontsize=8)
    axes[1].text(1.2e-3, 3.0e-6, "near-zero overlap except high Delta s", fontsize=8)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.out_dir / "tauu_pos_uubaruubar_uubar_to_uubar_detM_im_upstream_mechanism.png"
    fig.savefig(out_path)
    plt.close(fig)
    print(str(out_path))


if __name__ == "__main__":
    main()
