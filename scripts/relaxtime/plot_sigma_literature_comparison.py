#!/usr/bin/env python3
"""Compare sigma(s) scan outputs against digitized literature points."""

from __future__ import annotations

import argparse
import csv
import math
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt


PROCESS_LABELS = {
    "ud_to_ud": "ud -> ud",
    "us_to_us": "us -> us",
    "udbar_to_udbar": "udbar -> udbar",
    "usbar_to_usbar": "usbar -> usbar",
}


def _find_project_root() -> Path:
    script_dir = Path(__file__).resolve().parent
    candidates = [script_dir, script_dir.parent, script_dir.parent.parent, Path.cwd()]

    for start in candidates:
        current = start
        for _ in range(5):
            if (current / "Project.toml").exists() or (current / ".git").exists():
                return current
            parent = current.parent
            if parent == current:
                break
            current = parent
    return Path.cwd()


PROJECT_ROOT = _find_project_root()
DEFAULT_SCAN_DIR = PROJECT_ROOT / "data" / "outputs" / "results" / "relaxtime" / "cross_section"
DEFAULT_LITERATURE = PROJECT_ROOT / "tests" / "validation" / "data" / "relaxtime_sigma_literature_digitized_longtable_v1.csv"
DEFAULT_OUT_DIR = PROJECT_ROOT / "data" / "outputs" / "figures" / "relaxtime" / "literature"
DEFAULT_SUMMARY = PROJECT_ROOT / "data" / "outputs" / "results" / "relaxtime" / "cross_section" / "sigma_literature_compare_summary.csv"
MB_PER_FM2 = 10.0


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan-dir", type=Path, default=DEFAULT_SCAN_DIR)
    parser.add_argument("--literature-csv", type=Path, default=DEFAULT_LITERATURE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--summary-csv", type=Path, default=DEFAULT_SUMMARY)
    return parser.parse_args()


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _read_scan_csv(path: Path) -> list[dict[str, str]]:
    data_lines: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            data_lines.append(line)
    return list(csv.DictReader(data_lines))


def _iter_scan_files(scan_dir: Path) -> list[Path]:
    files = []
    for path in sorted(scan_dir.rglob("xs_T*_muB*_xi0p0.csv")):
        if path.name.startswith("_"):
            continue
        files.append(path)
    return files


def _interpolate(xs: list[float], ys: list[float], x: float) -> float:
    if not xs:
        return math.nan
    if x < xs[0] or x > xs[-1]:
        return math.nan
    if x == xs[0]:
        return ys[0]
    if x == xs[-1]:
        return ys[-1]

    index = bisect_left(xs, x)
    if index < len(xs) and xs[index] == x:
        return ys[index]

    x0 = xs[index - 1]
    x1 = xs[index]
    y0 = ys[index - 1]
    y1 = ys[index]
    if x1 == x0:
        return y0
    weight = (x - x0) / (x1 - x0)
    return y0 + weight * (y1 - y0)


def _configure_matplotlib() -> None:
    matplotlib.rcParams.update({
        "figure.dpi": 160,
        "savefig.dpi": 300,
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "legend.fontsize": 8,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "axes.grid": True,
        "grid.alpha": 0.25,
    })


def _series_label(row: dict[str, str]) -> str:
    return f"T={int(float(row['T_MeV']))} MeV, mu_B={int(float(row['muB_MeV']))} MeV"


def _series_style(literature_rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    palette = [
        "#1f77b4",
        "#d62728",
        "#2ca02c",
        "#ff7f0e",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
        "#bcbd22",
        "#17becf",
    ]
    styles: dict[str, dict[str, str]] = {}
    for index, series in enumerate(sorted({row["series"] for row in literature_rows})):
        sample = next(row for row in literature_rows if row["series"] == series)
        styles[series] = {
            "color": palette[index % len(palette)],
            "label": _series_label(sample),
        }
    return styles


def main() -> None:
    args = _parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
    _configure_matplotlib()

    literature_rows = _read_csv(args.literature_csv)
    scan_files = _iter_scan_files(args.scan_dir)

    model_series: dict[tuple[str, float, float], dict[str, object]] = {}
    for path in scan_files:
        rows = _read_scan_csv(path)
        if not rows:
            continue
        required = {"process", "T_MeV", "muB_MeV", "sqrt_s_MeV", "sigma"}
        if not required.issubset(rows[0].keys()):
            continue
        grouped: dict[tuple[str, float, float], list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            key = (row["process"], float(row["T_MeV"]), float(row["muB_MeV"]))
            grouped[key].append(row)

        for key, items in grouped.items():
            items = sorted(items, key=lambda item: float(item["sqrt_s_MeV"]))
            xs = [float(item["sqrt_s_MeV"]) / 1000.0 for item in items]
            ys = [float(item["sigma"]) * MB_PER_FM2 for item in items]
            model_series[key] = {"xs": xs, "ys": ys, "scan_file": path.name}

    summary_rows: list[dict[str, str]] = []
    stats_by_process: dict[str, list[float]] = defaultdict(list)
    stats_by_series: dict[str, list[float]] = defaultdict(list)

    for row in literature_rows:
        process = row["process"]
        T_mev = float(row["T_MeV"])
        muB_mev = float(row["muB_MeV"])
        sqrt_s_gev = float(row["sqrt_s_GeV"])
        sigma_lit = float(row["sigma_mb"])
        model = model_series.get((process, T_mev, muB_mev))
        model_sigma = math.nan
        scan_file = ""
        if model is not None:
            model_sigma = _interpolate(model["xs"], model["ys"], sqrt_s_gev)
            scan_file = str(model["scan_file"])

        rel_error = math.nan
        abs_rel_error = math.nan
        if sigma_lit != 0.0 and math.isfinite(model_sigma):
            rel_error = (model_sigma - sigma_lit) / sigma_lit
            abs_rel_error = abs(rel_error)
            stats_by_process[process].append(abs_rel_error)
            stats_by_series[row["series"]].append(abs_rel_error)

        summary_rows.append({
            "point_id": row["point_id"],
            "series": row["series"],
            "process": process,
            "T_MeV": row["T_MeV"],
            "muB_MeV": row["muB_MeV"],
            "sqrt_s_GeV": row["sqrt_s_GeV"],
            "literature_sigma_mb": row["sigma_mb"],
            "model_sigma_mb": "" if not math.isfinite(model_sigma) else f"{model_sigma:.12g}",
            "relative_error": "" if not math.isfinite(rel_error) else f"{rel_error:.12g}",
            "abs_relative_error": "" if not math.isfinite(abs_rel_error) else f"{abs_rel_error:.12g}",
            "scan_file": scan_file,
        })

    with args.summary_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    styles = _series_style(literature_rows)
    process_order = [
        "ud_to_ud",
        "us_to_us",
        "udbar_to_udbar",
        "usbar_to_usbar",
    ]

    overlay_fig, overlay_axes = plt.subplots(2, 2, figsize=(13.5, 9.2), sharex=False, sharey=False)
    error_fig, error_axes = plt.subplots(2, 2, figsize=(13.5, 9.2), sharex=False, sharey=False)

    for process, overlay_axis, error_axis in zip(process_order, overlay_axes.ravel(), error_axes.ravel()):
        overlay_axis.set_title(PROCESS_LABELS[process])
        overlay_axis.set_xlabel("sqrt(s) [GeV]")
        overlay_axis.set_ylabel("sigma [mb]")

        error_axis.set_title(PROCESS_LABELS[process])
        error_axis.set_xlabel("sqrt(s) [GeV]")
        error_axis.set_ylabel("|relative error|")
        error_axis.set_ylim(0.0, 2.5)
        error_axis.axhline(0.2, color="#999999", linewidth=0.8, linestyle="--")
        error_axis.axhline(0.5, color="#999999", linewidth=0.8, linestyle=":")
        error_axis.axhline(1.0, color="#bbbbbb", linewidth=0.8, linestyle="-.")

        process_rows = [row for row in literature_rows if row["process"] == process]
        for series in sorted({row["series"] for row in process_rows}):
            rows = [row for row in process_rows if row["series"] == series]
            rows = sorted(rows, key=lambda item: float(item["sqrt_s_GeV"]))
            style = styles[series]
            label = style["label"]
            lit_x = [float(row["sqrt_s_GeV"]) for row in rows]
            lit_y = [float(row["sigma_mb"]) for row in rows]
            overlay_axis.scatter(
                lit_x,
                lit_y,
                s=18,
                facecolors="none",
                edgecolors=style["color"],
                linewidths=0.9,
                label=f"wpd: {label}",
            )

            key = (process, float(rows[0]["T_MeV"]), float(rows[0]["muB_MeV"]))
            model = model_series.get(key)
            if model is not None:
                overlay_axis.plot(
                    model["xs"],
                    model["ys"],
                    color=style["color"],
                    linewidth=1.4,
                    label=f"model: {label}",
                )

            err_x = []
            err_y = []
            if model is not None:
                for row in rows:
                    sigma_lit = float(row["sigma_mb"])
                    sigma_model = _interpolate(model["xs"], model["ys"], float(row["sqrt_s_GeV"]))
                    if sigma_lit == 0.0 or not math.isfinite(sigma_model):
                        continue
                    err_x.append(float(row["sqrt_s_GeV"]))
                    err_y.append(abs((sigma_model - sigma_lit) / sigma_lit))
            if err_x:
                error_axis.plot(err_x, err_y, color=style["color"], linewidth=1.1, marker="o", markersize=2.6, label=label)

        overlay_axis.legend(loc="best", frameon=False)
        error_axis.legend(loc="best", frameon=False)

    overlay_fig.suptitle("Sigma literature comparison by process")
    overlay_fig.tight_layout()
    overlay_png = args.out_dir / "sigma_literature_overlay_by_process.png"
    overlay_pdf = args.out_dir / "sigma_literature_overlay_by_process.pdf"
    overlay_fig.savefig(overlay_png, bbox_inches="tight")
    overlay_fig.savefig(overlay_pdf, bbox_inches="tight")
    plt.close(overlay_fig)

    error_fig.suptitle("Sigma literature absolute relative error by process")
    error_fig.tight_layout()
    error_png = args.out_dir / "sigma_literature_error_by_process.png"
    error_pdf = args.out_dir / "sigma_literature_error_by_process.pdf"
    error_fig.savefig(error_png, bbox_inches="tight")
    error_fig.savefig(error_pdf, bbox_inches="tight")
    plt.close(error_fig)

    print(f"Saved overlay figure: {overlay_png}")
    print(f"Saved overlay figure: {overlay_pdf}")
    print(f"Saved error figure: {error_png}")
    print(f"Saved error figure: {error_pdf}")
    print(f"Saved summary: {args.summary_csv}")

    for process in process_order:
        errors = stats_by_process.get(process, [])
        if not errors:
            continue
        mean_abs = sum(errors) / len(errors)
        max_abs = max(errors)
        print(
            f"{process}: mean_abs_rel_error={mean_abs:.4f}, "
            f"max_abs_rel_error={max_abs:.4f}, points={len(errors)}"
        )

    for series in sorted(stats_by_series):
        errors = stats_by_series[series]
        if not errors:
            continue
        mean_abs = sum(errors) / len(errors)
        max_abs = max(errors)
        print(
            f"{series}: mean_abs_rel_error={mean_abs:.4f}, "
            f"max_abs_rel_error={max_abs:.4f}, points={len(errors)}"
        )


if __name__ == "__main__":
    main()