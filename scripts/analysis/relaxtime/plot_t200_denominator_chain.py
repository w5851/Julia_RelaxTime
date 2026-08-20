#!/usr/bin/env python3
"""Plot denominator-to-rate chain for T200 window diagnostics."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
from cycler import cycler


APS_COLOR_CYCLE = [
    "#4477AA",
    "#EE6677",
    "#228833",
    "#CCBB44",
    "#66CCEE",
    "#AA3377",
    "#BBBBBB",
]


def configure_publication_style(dpi: int) -> None:
    matplotlib.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times"],
            "font.size": 10,
            "mathtext.fontset": "stix",
            "axes.prop_cycle": cycler(color=APS_COLOR_CYCLE),
            "axes.titlesize": 10,
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
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": dpi,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.08,
        }
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sample-csv",
        type=Path,
        default=Path(r"D:\Desktop\Temp\relaxtime_t200_window\t200_denominator_ds_samples.csv"),
    )
    parser.add_argument(
        "--band-csv",
        type=Path,
        default=Path(r"D:\Desktop\Temp\relaxtime_t200_window\t200_fullchain_band_table.csv"),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("docs/analysis/relaxtime/transport/t200_tauu_spikes"),
    )
    parser.add_argument("--scenario", action="append", default=[])
    parser.add_argument("--formats", type=str, default="png")
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--show-title", action="store_true")
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def as_float(value: str) -> float:
    try:
        return float(value)
    except Exception:
        return math.nan


def finite(values: list[float]) -> bool:
    return any(math.isfinite(v) for v in values)


def robust_abs_limit(values: list[float], quantile: float, floor: float, ceil: float) -> float:
    clean = sorted(abs(v) for v in values if math.isfinite(v) and abs(v) > 0.0)
    if not clean:
        return floor
    idx = int((len(clean) - 1) * quantile)
    base = clean[max(0, min(idx, len(clean) - 1))]
    return max(floor, min(ceil, 1.35 * base))


def build_groups(rows: list[dict[str, str]], scenarios: set[str]) -> list[tuple[str, str]]:
    groups = sorted({(r["scenario"], r["process"]) for r in rows})
    if not scenarios:
        return groups
    return [g for g in groups if g[0] in scenarios]


def band_label(row: dict[str, str]) -> str:
    left = as_float(row["ds_left"])
    right = as_float(row["ds_right"])
    if math.isinf(right):
        return f"({left:g},inf]"
    return f"({left:g},{right:g}]"


def plot_group(
    scenario: str,
    process: str,
    sample_rows: list[dict[str, str]],
    band_rows: list[dict[str, str]],
    out_dir: Path,
    show_title: bool,
    formats: list[str],
) -> list[Path]:
    xi_values = sorted({as_float(r["xi"]) for r in sample_rows})
    fig, axes = plt.subplots(
        4,
        1,
        figsize=(7.0, 12.2),
        constrained_layout=True,
        gridspec_kw={"height_ratios": [1.15, 0.9, 1.15, 1.0]},
    )

    zoom_candidates: list[float] = []
    near0_candidates: list[float] = []
    near0_ds_max = 0.6

    for xi in xi_values:
        rows_xi = sorted(
            [r for r in sample_rows if as_float(r["xi"]) == xi],
            key=lambda r: as_float(r["ds"]),
        )
        x = [as_float(r["ds"]) for r in rows_xi]

        den_re = [as_float(r["den_simple_re"]) for r in rows_xi]
        den_im = [as_float(r["den_simple_im"]) for r in rows_xi]
        det_re = [as_float(r["detM_re"]) for r in rows_xi]
        det_im = [as_float(r["detM_im"]) for r in rows_xi]

        axes[0].plot(x, den_re, "--", linewidth=1.1, label=f"xi={xi:.2f} den_re")
        axes[0].plot(x, den_im, ":", linewidth=1.1, label=f"xi={xi:.2f} den_im")
        if finite(det_re):
            axes[0].plot(x, det_re, "-", linewidth=1.4, label=f"xi={xi:.2f} detM_re")
            axes[1].plot(x, det_re, "-", linewidth=1.35, marker="o", markersize=2.5, label=f"xi={xi:.2f} detM_re")
            zoom_candidates.extend([abs(v) for v in det_re if math.isfinite(v)])
            near0_candidates.extend(
                [abs(v) for xv, v in zip(x, det_re) if math.isfinite(v) and xv <= near0_ds_max]
            )
        if finite(det_im):
            axes[0].plot(x, det_im, "-.", linewidth=1.2, label=f"xi={xi:.2f} detM_im")
            axes[1].plot(x, det_im, "-.", linewidth=1.1, marker="s", markersize=2.4, label=f"xi={xi:.2f} detM_im")
            zoom_candidates.extend([abs(v) for v in det_im if math.isfinite(v)])
            near0_candidates.extend(
                [abs(v) for xv, v in zip(x, det_im) if math.isfinite(v) and xv <= near0_ds_max]
            )

        invabs_den = [max(as_float(r["invabs_den_simple"]), 1.0e-30) for r in rows_xi]
        abs_d_total = [max(as_float(r["abs_D_total_sq"]), 1.0e-30) for r in rows_xi]
        axes[2].plot(x, invabs_den, "--", linewidth=1.1, label=f"xi={xi:.2f} inv|den|^2")
        axes[2].plot(x, abs_d_total, "-", linewidth=1.4, label=f"xi={xi:.2f} |D_total|^2")

        invabs_det = [as_float(r["invabs_detM"]) for r in rows_xi]
        if finite(invabs_det):
            invabs_det_clip = [max(v, 1.0e-30) if math.isfinite(v) else math.nan for v in invabs_det]
            axes[2].plot(x, invabs_det_clip, ":", linewidth=1.1, label=f"xi={xi:.2f} inv|detM|^2")

        bands_xi = sorted(
            [r for r in band_rows if as_float(r["xi"]) == xi],
            key=lambda r: int(r["band"]),
        )
        bx = list(range(len(bands_xi)))
        omega_sigma = [as_float(r["omega_sigma_bin"]) for r in bands_xi]
        rate_bin = [as_float(r["rate_bin"]) for r in bands_xi]

        axes[3].plot(bx, omega_sigma, "-o", linewidth=1.2, markersize=3, label=f"xi={xi:.2f} omega_sigma")
        axes[3].plot(bx, rate_bin, "--s", linewidth=1.0, markersize=3, label=f"xi={xi:.2f} rate_bin")

    axes[0].set_xscale("log")
    axes[0].set_yscale("symlog", linthresh=1.0e-3)
    axes[0].axhline(0.0, color="#666666", linewidth=0.7)
    axes[0].set_ylabel("denominator kernel (full)")
    axes[0].set_xlabel("Delta s")

    zoom_max = 3.0e-4
    if zoom_candidates:
        sorted_vals = sorted(v for v in zoom_candidates if v > 0.0)
        if sorted_vals:
            p80 = sorted_vals[min(len(sorted_vals) - 1, int(0.8 * (len(sorted_vals) - 1)))]
            zoom_max = max(3.0e-4, min(2.0e-2, 1.4 * p80))
    axes[1].set_xscale("log")
    axes[1].set_yscale("symlog", linthresh=1.0e-6)
    axes[1].axhline(0.0, color="#666666", linewidth=0.7)
    near0_ylim = robust_abs_limit(near0_candidates, quantile=0.95, floor=2.0e-5, ceil=8.0e-3)
    axes[1].set_ylim(-near0_ylim, near0_ylim)
    axes[1].set_xlim(1.0e-3, near0_ds_max)
    axes[1].set_ylabel("den kernel near 0")
    axes[1].set_xlabel("Delta s")

    axes[2].set_xscale("log")
    axes[2].set_yscale("log")
    axes[2].set_ylabel("propagator magnitude")
    axes[2].set_xlabel("Delta s")

    band_ticks = sorted([r for r in band_rows if as_float(r["xi"]) == xi_values[0]], key=lambda r: int(r["band"]))
    axes[3].set_xticks(list(range(len(band_ticks))))
    axes[3].set_xticklabels([band_label(r) for r in band_ticks], rotation=30, ha="right")
    axes[3].set_ylabel("cross-section and rate bins")
    axes[3].set_xlabel("Delta s bands")

    axes[0].legend(ncol=2, loc="upper left")
    axes[1].legend(ncol=2, loc="upper left")
    axes[2].legend(ncol=2, loc="best")
    axes[3].legend(ncol=2, loc="upper left")

    if show_title:
        fig.suptitle(f"{scenario} / {process}: denominator -> propagator -> rate")

    out_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{scenario}_{process}_denominator_chain"
    saved: list[Path] = []
    for fmt in formats:
        path = out_dir / f"{stem}.{fmt}"
        fig.savefig(path)
        saved.append(path)
    plt.close(fig)
    return saved


def main() -> None:
    args = parse_args()
    configure_publication_style(args.dpi)

    sample_rows = read_csv(args.sample_csv)
    band_rows = read_csv(args.band_csv)
    groups = build_groups(sample_rows, set(args.scenario))
    if not groups:
        raise RuntimeError("No scenario/process groups selected for plotting")

    formats = [s.strip().lower() for s in args.formats.split(",") if s.strip()]
    if not formats:
        formats = ["png"]

    for scenario, process in groups:
        sample_group = [r for r in sample_rows if r["scenario"] == scenario and r["process"] == process]
        band_group = [r for r in band_rows if r["scenario"] == scenario and r["process"] == process]
        if not sample_group or not band_group:
            continue
        paths = plot_group(scenario, process, sample_group, band_group, args.out_dir, args.show_title, formats)
        for path in paths:
            print(f"Saved: {path}")


if __name__ == "__main__":
    main()
