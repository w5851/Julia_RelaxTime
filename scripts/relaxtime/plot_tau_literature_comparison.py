#!/usr/bin/env python3
"""Overlay current tau(T) scan against digitized literature points."""

from __future__ import annotations

import argparse
import csv
import math
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt


SERIES_TO_FIELD = {
    "tau_u_muB0": "tau_u",
    "tau_s_muB0": "tau_s",
    "tau_u_muB800": "tau_u",
    "tau_ubar_muB800": "tau_ubar",
    "tau_s_muB800": "tau_s",
    "tau_sbar_muB800": "tau_sbar",
}

SERIES_STYLE = {
    "tau_u_muB0": {"color": "#1f77b4", "label": "u, mu_B=0"},
    "tau_s_muB0": {"color": "#d62728", "label": "s, mu_B=0"},
    "tau_u_muB800": {"color": "#1f77b4", "label": "u, mu_B=800"},
    "tau_ubar_muB800": {"color": "#17becf", "label": "ubar, mu_B=800"},
    "tau_s_muB800": {"color": "#d62728", "label": "s, mu_B=800"},
    "tau_sbar_muB800": {"color": "#ff7f0e", "label": "sbar, mu_B=800"},
}


def _find_project_root() -> Path:
    script_dir = Path(__file__).resolve().parent
    candidates = [script_dir, script_dir.parent, script_dir.parent.parent]
    candidates.append(Path.cwd())

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
DEFAULT_SCAN = PROJECT_ROOT / "data" / "outputs" / "results" / "relaxtime" / "scan" / "relaxation_times_vs_T_literature_compare.csv"
DEFAULT_LITERATURE = PROJECT_ROOT / "tests" / "validation" / "data" / "relaxtime_tau_literature_digitized_longtable_v1.csv"
DEFAULT_OUT_DIR = PROJECT_ROOT / "data" / "outputs" / "figures" / "relaxtime" / "literature"
DEFAULT_SUMMARY = PROJECT_ROOT / "data" / "outputs" / "results" / "relaxtime" / "scan" / "relaxation_times_vs_T_literature_compare_summary.csv"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan-csv", type=Path, default=DEFAULT_SCAN)
    parser.add_argument("--literature-csv", type=Path, default=DEFAULT_LITERATURE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--summary-csv", type=Path, default=DEFAULT_SUMMARY)
    return parser.parse_args()


def _read_scan_csv(path: Path) -> list[dict[str, str]]:
    data_lines: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            data_lines.append(line)
    return list(csv.DictReader(data_lines))


def _read_scan_metadata(path: Path) -> dict[str, str]:
    metadata: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped.startswith("#"):
                continue
            payload = stripped[1:].strip()
            if ":" not in payload:
                continue
            key, value = payload.split(":", 1)
            metadata[key.strip()] = value.strip()
    return metadata


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


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


def main() -> None:
    args = _parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
    _configure_matplotlib()

    scan_metadata = _read_scan_metadata(args.scan_csv)
    scan_rows = _read_scan_csv(args.scan_csv)
    literature_rows = _read_csv(args.literature_csv)

    provenance_entrypoint = scan_metadata.get("provenance.entrypoint", "unknown")
    provenance_tau_path = scan_metadata.get("provenance.tau_path", "unknown")
    git_commit = scan_metadata.get("git_commit", "unknown")

    model_series: dict[str, tuple[list[float], list[float]]] = {}
    grouped_scan: dict[tuple[str, float], list[dict[str, str]]] = defaultdict(list)
    for row in scan_rows:
        muB = float(row["muB_MeV"])
        grouped_scan[("muB", muB)].append(row)

    for series, field in SERIES_TO_FIELD.items():
        muB = 0.0 if series.endswith("muB0") else 800.0
        rows = sorted(grouped_scan[("muB", muB)], key=lambda item: float(item["T_MeV"]))
        xs = [float(row["T_MeV"]) for row in rows if row.get(field)]
        ys = [float(row[field]) for row in rows if row.get(field)]
        model_series[series] = (xs, ys)

    summary_rows: list[dict[str, str]] = []
    stats: dict[str, list[float]] = defaultdict(list)
    literature_by_mub: dict[float, list[dict[str, str]]] = defaultdict(list)

    for row in literature_rows:
        literature_tau = float(row["tau_fm"])
        T_mev = float(row["T_MeV"])
        muB = float(row["muB_MeV"])
        series = row["series"]
        xs, ys = model_series[series]
        model_tau = _interpolate(xs, ys, T_mev)
        rel_error = math.nan
        if literature_tau != 0.0 and math.isfinite(model_tau):
            rel_error = (model_tau - literature_tau) / literature_tau
            stats[series].append(abs(rel_error))

        summary_rows.append({
            "point_id": row["point_id"],
            "series": series,
            "muB_MeV": row["muB_MeV"],
            "T_MeV": row["T_MeV"],
            "literature_tau_fm": row["tau_fm"],
            "model_tau_fm": "" if not math.isfinite(model_tau) else f"{model_tau:.12g}",
            "relative_error": "" if not math.isfinite(rel_error) else f"{rel_error:.12g}",
            "scan_entrypoint": provenance_entrypoint,
            "scan_tau_path": provenance_tau_path,
            "scan_git_commit": git_commit,
        })
        literature_by_mub[muB].append(row)

    with args.summary_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    for axis, muB in zip(axes, [0.0, 800.0]):
        axis.set_title(f"mu_B = {int(muB)} MeV")
        axis.set_xlabel("T [MeV]")
        axis.set_ylabel("tau [fm]")
        axis.set_yscale("log")
        axis.set_xlim(130.0, 405.0)
        axis.set_ylim(0.2, 30.0)

        series_names = [name for name in SERIES_TO_FIELD if (name.endswith("muB0") if muB == 0.0 else name.endswith("muB800"))]
        for series in series_names:
            xs, ys = model_series[series]
            style = SERIES_STYLE[series]
            axis.plot(xs, ys, color=style["color"], linewidth=1.6, label=f"model: {style['label']}")

            literature_points = [row for row in literature_by_mub[muB] if row["series"] == series]
            lit_x = [float(row["T_MeV"]) for row in literature_points]
            lit_y = [float(row["tau_fm"]) for row in literature_points]
            axis.scatter(
                lit_x,
                lit_y,
                s=18,
                facecolors="none",
                edgecolors=style["color"],
                linewidths=0.9,
                label=f"wpd: {style['label']}",
            )

        axis.legend(loc="upper right", ncol=1, frameon=False)

    fig.suptitle(f"Relaxation time literature comparison ({provenance_entrypoint})")
    fig.tight_layout()
    png_path = args.out_dir / "tau_literature_comparison.png"
    pdf_path = args.out_dir / "tau_literature_comparison.pdf"
    fig.savefig(png_path, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved figure: {png_path}")
    print(f"Saved figure: {pdf_path}")
    print(f"Saved summary: {args.summary_csv}")
    print(f"Scan provenance: entrypoint={provenance_entrypoint}, tau_path={provenance_tau_path}, git_commit={git_commit}")
    for series, errors in sorted(stats.items()):
        if not errors:
            continue
        mean_abs = sum(errors) / len(errors)
        max_abs = max(errors)
        print(f"{series}: mean_abs_rel_error={mean_abs:.4f}, max_abs_rel_error={max_abs:.4f}, points={len(errors)}")


if __name__ == "__main__":
    main()