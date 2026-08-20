#!/usr/bin/env python3
"""Plot t-dependence scans for detM_im upstream and propagator quantities."""

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
        default=Path(r"D:\Desktop\Temp\relaxtime_t200_window\t200_detm_im_vs_t_scan.csv"),
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=Path(r"D:\Desktop\Temp\relaxtime_t200_window\t200_detm_im_vs_t_scan_summary.csv"),
    )
    parser.add_argument("--out-dir", type=Path, default=Path("docs/analysis/relaxtime/transport/t200_tauu_spikes"))
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
                    "t_alpha": float(r["t_alpha"]),
                    "k0_s": float(r["k0_s"]),
                    "k_s": float(r["k_s"]),
                    "detM_im": float(r["detM_im"]),
                    "invabs_detM": float(r["invabs_detM"]),
                    "abs_Dmixed_sq": float(r["abs_Dmixed_sq"]),
                }
            )
    return rows


def write_summary(rows: list[dict[str, float]], out_path: Path) -> None:
    keys = sorted({(r["xi"], r["ds"]) for r in rows})
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "xi",
                "ds",
                "k0_min",
                "k0_max",
                "k_abs_max",
                "detM_im_min",
                "detM_im_max",
                "detM_im_range_abs",
                "detM_im_range_rel_to_mid",
                "invabs_detM_range_rel_to_mid",
                "abs_Dmixed_sq_range_rel_to_mid",
            ]
        )
        for xi, ds in keys:
            sub = [r for r in rows if r["xi"] == xi and r["ds"] == ds]
            sub_sorted = sorted(sub, key=lambda r: r["t_alpha"])
            mid = min(sub_sorted, key=lambda r: abs(r["t_alpha"] - 0.5))
            k0_vals = [r["k0_s"] for r in sub]
            k_vals = [r["k_s"] for r in sub]
            det_vals = [r["detM_im"] for r in sub]
            inv_vals = [r["invabs_detM"] for r in sub]
            dm_vals = [r["abs_Dmixed_sq"] for r in sub]
            det_range = max(det_vals) - min(det_vals)
            det_rel = det_range / max(abs(mid["detM_im"]), 1.0e-30)
            inv_rel = (max(inv_vals) - min(inv_vals)) / max(abs(mid["invabs_detM"]), 1.0e-30)
            dm_rel = (max(dm_vals) - min(dm_vals)) / max(abs(mid["abs_Dmixed_sq"]), 1.0e-30)
            writer.writerow(
                [
                    f"{xi:.2f}",
                    f"{ds:.12g}",
                    f"{min(k0_vals):.12e}",
                    f"{max(k0_vals):.12e}",
                    f"{max(abs(v) for v in k_vals):.12e}",
                    f"{min(det_vals):.12e}",
                    f"{max(det_vals):.12e}",
                    f"{det_range:.12e}",
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

    xi_values = sorted({r["xi"] for r in rows})
    ds_values = sorted({r["ds"] for r in rows})

    style_map = {
        0.34: {"ls": "-", "mk": "o"},
        0.36: {"ls": "--", "mk": "s"},
        0.38: {"ls": "-.", "mk": "^"},
    }

    fig, axes = plt.subplots(4, len(ds_values), figsize=(12.6, 8.8), constrained_layout=True)

    for col, ds in enumerate(ds_values):
        ax_k0 = axes[0, col]
        ax_det = axes[1, col]
        ax_inv = axes[2, col]
        ax_dm = axes[3, col]

        for xi in xi_values:
            sub = sorted(
                [r for r in rows if abs(r["xi"] - xi) < 1.0e-12 and abs(r["ds"] - ds) < 1.0e-12],
                key=lambda r: r["t_alpha"],
            )
            x = [r["t_alpha"] for r in sub]
            st = style_map.get(xi, {"ls": "-", "mk": "o"})
            label = f"xi={xi:.2f}"
            ax_k0.plot(x, [r["k0_s"] for r in sub], linestyle=st["ls"], marker=st["mk"], markersize=2.5, markevery=6, linewidth=1.2, label=label)
            ax_det.plot(x, [r["detM_im"] for r in sub], linestyle=st["ls"], marker=st["mk"], markersize=2.5, markevery=6, linewidth=1.2, label=label)
            ax_inv.plot(x, [r["invabs_detM"] for r in sub], linestyle=st["ls"], marker=st["mk"], markersize=2.5, markevery=6, linewidth=1.2, label=label)
            ax_dm.plot(x, [r["abs_Dmixed_sq"] for r in sub], linestyle=st["ls"], marker=st["mk"], markersize=2.5, markevery=6, linewidth=1.2, label=label)

        ax_det.axhline(0.0, color="#666666", linewidth=0.6)
        ax_k0.set_title(f"Delta s={ds:g}")

    for col in range(len(ds_values)):
        axes[0, col].set_ylabel("k0 (s-channel)")
        axes[1, col].set_ylabel("detM_im")
        axes[2, col].set_ylabel("inv|detM|^2")
        axes[3, col].set_ylabel("|D_mixed|^2")
        axes[3, col].set_xlabel("t_alpha: 0=t_min, 1=t_max")

    axes[0, -1].legend(ncol=1, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.out_dir / "tauu_pos_uubaruubar_uubar_to_uubar_detM_im_vs_t_scan.png"
    fig.savefig(out_path)
    plt.close(fig)
    print(str(out_path))
    print(str(args.summary_csv))


if __name__ == "__main__":
    main()
