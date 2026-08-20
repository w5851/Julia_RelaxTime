#!/usr/bin/env python3
"""Compare t-scan behavior between s-channel and t-channel upstream quantities."""

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
        default=Path(r"D:\Desktop\Temp\relaxtime_t200_window\t200_detm_im_s_vs_t_channel_tscan.csv"),
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=Path(r"D:\Desktop\Temp\relaxtime_t200_window\t200_detm_im_s_vs_t_channel_tscan_summary.csv"),
    )
    parser.add_argument("--out-dir", type=Path, default=Path("docs/analysis/relaxtime/transport/t200_tauu_spikes"))
    parser.add_argument("--dpi", type=int, default=600)
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for r in reader:
            rows.append(
                {
                    "channel": r["channel"],
                    "xi": float(r["xi"]),
                    "ds": float(r["ds"]),
                    "t_alpha": float(r["t_alpha"]),
                    "k0": float(r["k0"]),
                    "k_norm": float(r["k_norm"]),
                    "detM_im": float(r["detM_im"]),
                    "invabs_detM": float(r["invabs_detM"]),
                    "abs_Dmixed_sq": float(r["abs_Dmixed_sq"]),
                }
            )
    return rows


def write_summary(rows: list[dict[str, float | str]], out_path: Path) -> None:
    keys = sorted({(str(r["channel"]), float(r["xi"]), float(r["ds"])) for r in rows})
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "channel",
                "xi",
                "ds",
                "k0_min",
                "k0_max",
                "k_norm_max",
                "detM_im_range_rel_to_mid",
                "invabs_detM_range_rel_to_mid",
                "abs_Dmixed_sq_range_rel_to_mid",
            ]
        )
        for channel, xi, ds in keys:
            sub = [r for r in rows if str(r["channel"]) == channel and float(r["xi"]) == xi and float(r["ds"]) == ds]
            sub_sorted = sorted(sub, key=lambda r: float(r["t_alpha"]))
            mid = min(sub_sorted, key=lambda r: abs(float(r["t_alpha"]) - 0.5))
            k0_vals = [float(r["k0"]) for r in sub]
            k_vals = [float(r["k_norm"]) for r in sub]
            det_vals = [float(r["detM_im"]) for r in sub]
            inv_vals = [float(r["invabs_detM"]) for r in sub]
            dm_vals = [float(r["abs_Dmixed_sq"]) for r in sub]
            det_rel = (max(det_vals) - min(det_vals)) / max(abs(float(mid["detM_im"])), 1.0e-30)
            inv_rel = (max(inv_vals) - min(inv_vals)) / max(abs(float(mid["invabs_detM"])), 1.0e-30)
            dm_rel = (max(dm_vals) - min(dm_vals)) / max(abs(float(mid["abs_Dmixed_sq"])), 1.0e-30)
            writer.writerow(
                [
                    channel,
                    f"{xi:.2f}",
                    f"{ds:.12g}",
                    f"{min(k0_vals):.12e}",
                    f"{max(k0_vals):.12e}",
                    f"{max(abs(v) for v in k_vals):.12e}",
                    f"{det_rel:.12e}",
                    f"{inv_rel:.12e}",
                    f"{dm_rel:.12e}",
                ]
            )


def main() -> None:
    args = parse_args()
    configure_style(args.dpi)
    rows = read_rows(args.csv)
    write_summary(rows, args.summary_csv)

    xi_target = 0.36
    ds_values = sorted({float(r["ds"]) for r in rows})
    channel_order = [":s", ":t"]
    channel_title = {":s": "s-channel", ":t": "t-channel"}

    fig, axes = plt.subplots(len(channel_order), 4, figsize=(13.2, 6.8), constrained_layout=True)

    for row_idx, channel in enumerate(channel_order):
        for col_idx, ds in enumerate(ds_values):
            ax = axes[row_idx, col_idx]
            sub = sorted(
                [
                    r
                    for r in rows
                    if str(r["channel"]) == channel and abs(float(r["xi"]) - xi_target) < 1.0e-12 and abs(float(r["ds"]) - ds) < 1.0e-12
                ],
                key=lambda r: float(r["t_alpha"]),
            )
            x = [float(r["t_alpha"]) for r in sub]
            ax.plot(x, [float(r["detM_im"]) for r in sub], label="detM_im", linewidth=1.35)
            ax.plot(x, [float(r["invabs_detM"]) for r in sub], "--", label="inv|detM|^2", linewidth=1.15)
            ax.plot(x, [float(r["abs_Dmixed_sq"]) for r in sub], ":", label="|D_mixed|^2", linewidth=1.25)
            ax.axhline(0.0, color="#666666", linewidth=0.6)
            ax.set_yscale("symlog", linthresh=2.0e-5)
            if row_idx == 0:
                ax.set_title(f"Delta s={ds:g}")
            if col_idx == 0:
                ax.set_ylabel(channel_title[channel])
            if row_idx == len(channel_order) - 1:
                ax.set_xlabel("t_alpha")

    axes[0, -1].legend(ncol=1, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.out_dir / "tauu_pos_uubaruubar_detM_im_s_vs_t_channel_tscan.png"
    fig.savefig(out_path)
    plt.close(fig)
    print(str(out_path))
    print(str(args.summary_csv))


if __name__ == "__main__":
    main()
