#!/usr/bin/env python3
"""Plot t-sensitivity companion figure and summary for detM_im upstream decomposition."""

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
        default=Path(r"D:\Desktop\Temp\relaxtime_t200_window\t200_detm_im_t_sensitivity_trace.csv"),
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=Path(r"D:\Desktop\Temp\relaxtime_t200_window\t200_detm_im_t_sensitivity_summary.csv"),
    )
    parser.add_argument("--out-dir", type=Path, default=Path("docs/analysis/relaxtime"))
    parser.add_argument("--dpi", type=int, default=600)
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for r in reader:
            rows.append(
                {
                    "xi": float(r["xi"]),
                    "t_mode": r["t_mode"],
                    "ds": float(r["ds"]),
                    "pi_uu_im": float(r["pi_uu_im"]),
                    "pi_ss_im": float(r["pi_ss_im"]),
                    "im_M00M88": float(r["im_M00M88"]),
                    "im_M08sq": float(r["im_M08sq"]),
                    "detM_im": float(r["detM_im"]),
                }
            )
    return rows


def write_summary(rows: list[dict[str, float | str]], out_path: Path) -> None:
    grouped: dict[tuple[float, float], dict[str, float]] = {}
    for r in rows:
        key = (float(r["xi"]), float(r["ds"]))
        grouped.setdefault(key, {})[str(r["t_mode"])] = float(r["detM_im"])

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["xi", "ds", "detM_im_t_min", "detM_im_t_mid", "detM_im_t_max", "range_abs", "range_rel_to_mid"])
        for (xi, ds) in sorted(grouped.keys()):
            item = grouped[(xi, ds)]
            vmin = item["t_min"]
            vmid = item["t_mid"]
            vmax = item["t_max"]
            rabs = max(vmin, vmid, vmax) - min(vmin, vmid, vmax)
            rrel = rabs / max(abs(vmid), 1.0e-30)
            writer.writerow([f"{xi:.2f}", f"{ds:.12g}", f"{vmin:.12e}", f"{vmid:.12e}", f"{vmax:.12e}", f"{rabs:.12e}", f"{rrel:.12e}"])


def shade_segments(ax: plt.Axes) -> None:
    ax.axvspan(1.0e-3, 5.0, color="#DDEAF7", alpha=0.25)
    ax.axvspan(5.0, 10.0, color="#FBE7D5", alpha=0.25)
    ax.axvspan(10.0, 20.0, color="#E2F3E8", alpha=0.25)


def main() -> None:
    args = parse_args()
    configure_style(args.dpi)
    rows = read_rows(args.csv)
    write_summary(rows, args.summary_csv)

    xi_values = sorted({float(r["xi"]) for r in rows})
    mode_values = ["t_min", "t_mid", "t_max"]
    mode_style = {
        "t_min": {"ls": "-", "mk": "o"},
        "t_mid": {"ls": "--", "mk": "s"},
        "t_max": {"ls": "-.", "mk": "^"},
    }

    fig, axes = plt.subplots(len(xi_values), 1, figsize=(8.6, 8.8), constrained_layout=True)
    if len(xi_values) == 1:
        axes = [axes]

    for ax, xi in zip(axes, xi_values):
        ax.set_xscale("log")
        ax.set_yscale("symlog", linthresh=2.0e-5)
        shade_segments(ax)
        ax.axhline(0.0, color="#666666", linewidth=0.7)

        for mode in mode_values:
            sub = sorted(
                [r for r in rows if abs(float(r["xi"]) - xi) < 1.0e-12 and str(r["t_mode"]) == mode],
                key=lambda r: float(r["ds"]),
            )
            x = [float(r["ds"]) for r in sub]
            y = [float(r["detM_im"]) for r in sub]
            st = mode_style[mode]
            ax.plot(
                x,
                y,
                linewidth=1.4,
                linestyle=st["ls"],
                marker=st["mk"],
                markersize=2.8,
                markerfacecolor="white",
                markeredgewidth=0.8,
                markevery=2,
                label=mode,
            )

        ax.set_ylabel(f"detM_im (xi={xi:.2f})")
        ax.legend(ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)

    axes[-1].set_xlabel("Delta s")
    axes[0].text(1.3e-3, 2.3e-2, "same shape, different t slices", fontsize=8)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.out_dir / "tauu_pos_uubaruubar_uubar_to_uubar_detM_im_t_sensitivity.png"
    fig.savefig(out_path)
    plt.close(fig)
    print(str(out_path))
    print(str(args.summary_csv))


if __name__ == "__main__":
    main()
