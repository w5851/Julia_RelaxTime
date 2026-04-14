#!/usr/bin/env python3
"""Plot term-level determinant imaginary-part decomposition for T200 window."""

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
        default=Path(r"D:\Desktop\Temp\relaxtime_t200_window\t200_detm_im_upstream_term_switch_trace.csv"),
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
                    "detM_im": float(r["detM_im"]),
                    "c1": float(r["c1_ReM00_ImM88"]),
                    "c2": float(r["c2_ImM00_ReM88"]),
                    "c3": float(r["c3_minus2ReM08ImM08"]),
                    "rebuild_error": float(r["rebuild_error"]),
                    "imM00_Piuu": float(r["imM00_from_Piuu"]),
                    "imM00_Piss": float(r["imM00_from_Piss"]),
                    "imM08_Piuu": float(r["imM08_from_Piuu"]),
                    "imM08_Piss": float(r["imM08_from_Piss"]),
                    "imM88_Piuu": float(r["imM88_from_Piuu"]),
                    "imM88_Piss": float(r["imM88_from_Piss"]),
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

    style_map = {
        0.34: {"ls": "-", "mk": "o", "me": 1},
        0.36: {"ls": "--", "mk": "s", "me": 2},
        0.38: {"ls": "-.", "mk": "^", "me": 3},
    }

    fig, axes = plt.subplots(3, 1, figsize=(8.4, 10.2), constrained_layout=True)
    for ax in axes:
        ax.set_xscale("log")
        shade_segments(ax)
        ax.axhline(0.0, color="#666666", linewidth=0.7)

    for xi in xi_values:
        sub = sorted([r for r in rows if abs(r["xi"] - xi) < 1.0e-12], key=lambda r: r["ds"])
        x = [r["ds"] for r in sub]
        st = style_map.get(xi, {"ls": "-", "mk": "o", "me": 1})

        axes[0].plot(
            x,
            [r["c1"] for r in sub],
            linewidth=1.4,
            linestyle=st["ls"],
            marker=st["mk"],
            markersize=2.8,
            markevery=st["me"],
            markerfacecolor="white",
            markeredgewidth=0.8,
            label=f"xi={xi:.2f} Re(M00)Im(M88)",
        )
        axes[0].plot(
            x,
            [r["c2"] for r in sub],
            linewidth=1.4,
            linestyle=st["ls"],
            marker=st["mk"],
            markersize=2.8,
            markevery=st["me"],
            markerfacecolor="white",
            markeredgewidth=0.8,
            alpha=0.7,
            label=f"xi={xi:.2f} Im(M00)Re(M88)",
        )
        axes[0].plot(
            x,
            [r["c3"] for r in sub],
            linewidth=1.4,
            linestyle=st["ls"],
            marker=st["mk"],
            markersize=2.8,
            markevery=st["me"],
            markerfacecolor="white",
            markeredgewidth=0.8,
            alpha=0.55,
            label=f"xi={xi:.2f} -2Re(M08)Im(M08)",
        )

        axes[1].plot(
            x,
            [r["detM_im"] for r in sub],
            linewidth=1.5,
            linestyle=st["ls"],
            marker=st["mk"],
            markersize=3.0,
            markevery=st["me"],
            markerfacecolor="white",
            markeredgewidth=0.9,
            label=f"xi={xi:.2f} detM_im",
        )

        axes[2].plot(
            x,
            [r["imM00_Piuu"] for r in sub],
            linewidth=1.2,
            linestyle=st["ls"],
            marker=st["mk"],
            markersize=2.4,
            markevery=st["me"],
            markerfacecolor="white",
            markeredgewidth=0.7,
            label=f"xi={xi:.2f} Im(M00)|Piuu",
        )
        axes[2].plot(
            x,
            [r["imM08_Piuu"] for r in sub],
            linewidth=1.2,
            linestyle=st["ls"],
            marker=st["mk"],
            markersize=2.4,
            markevery=st["me"],
            markerfacecolor="white",
            markeredgewidth=0.7,
            alpha=0.75,
            label=f"xi={xi:.2f} Im(M08)|Piuu",
        )
        axes[2].plot(
            x,
            [r["imM88_Piuu"] for r in sub],
            linewidth=1.2,
            linestyle=st["ls"],
            marker=st["mk"],
            markersize=2.4,
            markevery=st["me"],
            markerfacecolor="white",
            markeredgewidth=0.7,
            alpha=0.6,
            label=f"xi={xi:.2f} Im(M88)|Piuu",
        )

    axes[0].set_ylabel("term contributions")
    axes[0].set_yscale("symlog", linthresh=1.0e-8)
    axes[1].set_ylabel("detM_im")
    axes[1].set_yscale("symlog", linthresh=1.0e-8)
    axes[2].set_ylabel("imag parts from Pi_uu")
    axes[2].set_yscale("symlog", linthresh=1.0e-8)
    axes[2].set_xlabel("Delta s")

    max_abs_err = max(abs(r["rebuild_error"]) for r in rows)
    axes[1].text(1.2e-3, 2.5e-3, f"max rebuild error={max_abs_err:.2e}", fontsize=8)
    axes[1].text(1.2e-3, 8.0e-4, "detM_im = c1 + c2 + c3", fontsize=8)
    axes[0].text(1.2e-3, 4.0e-3, "watch dominant term switch", fontsize=8)
    axes[2].text(1.2e-3, 2.0e-3, "Pi_ss-driven imag parts ~ 0 (not shown)", fontsize=8)

    axes[0].legend(ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
    axes[1].legend(ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
    axes[2].legend(ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.out_dir / "tauu_pos_uubaruubar_uubar_to_uubar_detM_im_term_switch.png"
    fig.savefig(out_path)
    plt.close(fig)
    print(str(out_path))


if __name__ == "__main__":
    main()
