#!/usr/bin/env python3
"""Plot Re(M08) zero-cross control equation for T200 positive window."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
from cycler import cycler


APS_COLOR_CYCLE = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB"]
PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_IO_DIR = PROJECT_ROOT / "data" / "outputs" / "tmp" / "relaxtime_t200_window"


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
            "savefig.pad_inches": 0.08,
        }
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--csv",
        type=Path,
        default=DEFAULT_IO_DIR / "t200_m08_re_zero_cross_trace.csv",
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=DEFAULT_IO_DIR / "t200_m08_re_zero_cross_summary.csv",
    )
    parser.add_argument("--out-dir", type=Path, default=Path("docs/analysis/relaxtime/transport/t200_tauu_spikes"))
    parser.add_argument("--x-scale", choices=("log", "linear"), default="log")
    parser.add_argument("--focus-ds-min", type=float, default=10.0)
    parser.add_argument("--focus-ds-max", type=float, default=15.0)
    parser.add_argument("--focus-y-pad", type=float, default=0.20)
    parser.add_argument("--dpi", type=int, default=600)
    return parser.parse_args()


def linear_fit(x: list[float], y: list[float]) -> tuple[float, float, float]:
    n = len(x)
    if n < 2:
        return 0.0, 0.0, 0.0
    mean_x = sum(x) / n
    mean_y = sum(y) / n
    sxx = sum((xi - mean_x) ** 2 for xi in x)
    sxy = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y))
    if sxx == 0.0:
        return 0.0, mean_y, 0.0
    a = sxy / sxx
    b = mean_y - a * mean_x
    y_hat = [a * xi + b for xi in x]
    sse = sum((yi - yhi) ** 2 for yi, yhi in zip(y, y_hat))
    sst = sum((yi - mean_y) ** 2 for yi in y)
    r2 = 1.0 if sst == 0.0 else (1.0 - sse / sst)
    return a, b, r2


def read_trace(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for r in reader:
            rows.append(
                {
                    "xi": float(r["xi"]),
                    "ds": float(r["ds"]),
                    "re_M08": float(r["re_M08_total"]),
                    "delta_re_pi": float(r["re_Piuu_minus_Piss"]),
                    "target_delta": float(r["re_Piuu_minus_Piss_target"]),
                    "re_m08_formula": float(r["re_M08_formula"]),
                    "residual_formula": float(r["residual_formula"]),
                }
            )
    return rows


def read_summary(path: Path) -> dict[float, float]:
    out: dict[float, float] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for r in reader:
            out[float(r["xi"])] = float(r["ds_zero_cross_linear"])
    return out


def shade_segments(ax: plt.Axes) -> None:
    ax.axvspan(1.0e-3, 5.0, color="#DDEAF7", alpha=0.25)
    ax.axvspan(5.0, 10.0, color="#FBE7D5", alpha=0.25)
    ax.axvspan(10.0, 20.0, color="#E2F3E8", alpha=0.25)


def focus_ylim(rows: list[dict[str, float]], ds_min: float, ds_max: float, pad_ratio: float) -> tuple[float, float]:
    local = [r["re_M08"] for r in rows if ds_min <= r["ds"] <= ds_max]
    if not local:
        local = [r["re_M08"] for r in rows]
    lo = min(local)
    hi = max(local)
    span = hi - lo
    if span <= 0.0:
        span = max(1.0e-4, abs(hi) * 0.2)
    pad = max(1.0e-5, span * max(0.0, pad_ratio))
    lo = min(lo - pad, -pad)
    hi = max(hi + pad, pad)
    return lo, hi


def main() -> None:
    args = parse_args()
    if not args.csv.is_file():
        raise FileNotFoundError(f"input trace csv not found: {args.csv}. Pass --csv to override.")
    if not args.summary_csv.is_file():
        raise FileNotFoundError(f"input summary csv not found: {args.summary_csv}. Pass --summary-csv to override.")
    configure_style(args.dpi)
    rows = read_trace(args.csv)
    zero_ds = read_summary(args.summary_csv)

    xi_values = sorted({r["xi"] for r in rows})
    style_map = {
        0.34: {"ls": "-", "mk": "o", "me": 6},
        0.36: {"ls": "--", "mk": "s", "me": 7},
        0.38: {"ls": "-.", "mk": "^", "me": 8},
    }

    fig, axes = plt.subplots(3, 1, figsize=(8.4, 10.2), constrained_layout=False)
    fig.subplots_adjust(left=0.10, right=0.78, top=0.97, bottom=0.08, hspace=0.22)
    for ax in axes:
        if args.x_scale == "log":
            ax.set_xscale("log")
            shade_segments(ax)
        ax.axhline(0.0, color="#666666", linewidth=0.7)

    for xi in xi_values:
        sub = sorted([r for r in rows if abs(r["xi"] - xi) < 1.0e-12], key=lambda r: r["ds"])
        x = [r["ds"] for r in sub]
        st = style_map.get(xi, {"ls": "-", "mk": "o", "me": 6})

        axes[0].plot(
            x,
            [r["re_M08"] for r in sub],
            linewidth=1.4,
            linestyle=st["ls"],
            marker=st["mk"],
            markersize=2.8,
            markevery=st["me"],
            markerfacecolor="white",
            markeredgewidth=0.8,
            label=f"xi={xi:.2f} Re(M08)",
        )

        local = [r for r in sub if 8.0 <= r["ds"] <= 16.0]
        if len(local) >= 3:
            xfit = [r["ds"] for r in local]
            yfit = [r["re_M08"] for r in local]
            a, b, r2 = linear_fit(xfit, yfit)
            yline = [a * xv + b for xv in xfit]
            axes[0].plot(xfit, yline, linewidth=1.0, linestyle=":", alpha=0.9)
            xpos = 0.04
            ypos = 0.90 - 0.08 * xi_values.index(xi)
            axes[0].text(xpos, ypos, f"xi={xi:.2f}: R^2(8<=ds<=16)={r2:.4f}", transform=axes[0].transAxes, fontsize=8)

        axes[1].plot(
            x,
            [r["delta_re_pi"] for r in sub],
            linewidth=1.4,
            linestyle=st["ls"],
            marker=st["mk"],
            markersize=2.8,
            markevery=st["me"],
            markerfacecolor="white",
            markeredgewidth=0.8,
            label=f"xi={xi:.2f} Re(Piuu-Piss)",
        )
        target = sub[0]["target_delta"]
        axes[1].axhline(target, linestyle=st["ls"], linewidth=1.0, alpha=0.75)

        axes[2].plot(
            x,
            [r["residual_formula"] for r in sub],
            linewidth=1.3,
            linestyle=st["ls"],
            marker=st["mk"],
            markersize=2.6,
            markevery=st["me"],
            markerfacecolor="white",
            markeredgewidth=0.7,
            label=f"xi={xi:.2f} residual",
        )

        if xi in zero_ds and zero_ds[xi] == zero_ds[xi]:
            axes[0].axvline(zero_ds[xi], linestyle=st["ls"], linewidth=0.9, alpha=0.7)

    axes[0].set_ylabel("Re(M08)")
    axes[1].set_ylabel("Re(Piuu-Piss)")
    axes[2].set_ylabel("M08 formula residual")
    axes[2].set_xlabel("Delta s")

    if args.x_scale == "log":
        axes[0].set_yscale("symlog", linthresh=1.0e-5)
        axes[1].set_yscale("symlog", linthresh=1.0e-5)
        axes[2].set_yscale("symlog", linthresh=1.0e-14)
        axes[0].text(0.02, 0.10, "Re(M08) crosses zero near Delta s ~ 10-15", transform=axes[0].transAxes, fontsize=8)
        axes[1].text(0.02, 0.10, "zero-cross when Re(Piuu-Piss) meets target", transform=axes[1].transAxes, fontsize=8)
        axes[2].text(0.02, 0.10, "formula residual stays near machine zero", transform=axes[2].transAxes, fontsize=8)
    else:
        for ax in axes:
            ax.set_xlim(args.focus_ds_min, args.focus_ds_max)
        y0_lo, y0_hi = focus_ylim(rows, args.focus_ds_min, args.focus_ds_max, args.focus_y_pad)
        axes[0].set_ylim(y0_lo, y0_hi)
        axes[0].text(0.02, 0.90, "linear x/y focus around zero-cross", transform=axes[0].transAxes, fontsize=8)
        axes[1].text(0.02, 0.90, "target-line intersection explains zero-cross", transform=axes[1].transAxes, fontsize=8)
        axes[2].text(0.02, 0.90, "residual near machine precision", transform=axes[2].transAxes, fontsize=8)

    axes[0].legend(ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
    axes[1].legend(ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
    axes[2].legend(ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    if args.x_scale == "linear":
        out_name = "tauu_pos_uubaruubar_uubar_to_uubar_reM08_zero_cross_linear_x.png"
    else:
        out_name = "tauu_pos_uubaruubar_uubar_to_uubar_reM08_zero_cross.png"
    out_path = args.out_dir / out_name
    fig.savefig(out_path)
    plt.close(fig)
    print(str(out_path))


if __name__ == "__main__":
    main()
