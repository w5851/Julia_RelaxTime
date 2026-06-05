#!/usr/bin/env python
"""Render FIG3-like heatmaps from combined meson-density scan CSV output."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_rows(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8") as f:
        data_lines = (line for line in f if not line.lstrip().startswith("#"))
        for row in csv.DictReader(data_lines):
            rows.append(row)
    return rows


def as_float(row: dict[str, str], field: str) -> float:
    try:
        return float(row.get(field, "nan"))
    except ValueError:
        return math.nan


def matrix(rows: list[dict[str, str]], regime: str, field: str):
    selected = [
        row
        for row in rows
        if row.get("regime") == regime and row.get("status") == "ok"
    ]
    mus = np.array(sorted({as_float(row, "muq_MeV") for row in selected}), dtype=float)
    temps = np.array(sorted({as_float(row, "T_MeV") for row in selected}), dtype=float)
    mus = mus[np.isfinite(mus)]
    temps = temps[np.isfinite(temps)]
    values = np.full((temps.size, mus.size), np.nan, dtype=float)
    mu_index = {v: i for i, v in enumerate(mus)}
    temp_index = {v: i for i, v in enumerate(temps)}
    for row in selected:
        mu = as_float(row, "muq_MeV")
        temp = as_float(row, "T_MeV")
        val = as_float(row, field)
        if not (math.isfinite(mu) and math.isfinite(temp) and math.isfinite(val)):
            continue
        if val < 0.0:
            continue
        values[temp_index[temp], mu_index[mu]] = val
    return mus, temps, values


def centers_to_edges(xs: np.ndarray) -> np.ndarray:
    if xs.size == 0:
        return np.array([], dtype=float)
    if xs.size == 1:
        return np.array([xs[0] - 0.5, xs[0] + 0.5], dtype=float)
    edges = np.empty(xs.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (xs[:-1] + xs[1:])
    edges[0] = xs[0] - (edges[1] - xs[0])
    edges[-1] = xs[-1] + (xs[-1] - edges[-2])
    return edges


def finite_max(mats: list[np.ndarray]) -> float:
    vals = np.concatenate([m[np.isfinite(m)] for m in mats if np.isfinite(m).any()])
    if vals.size == 0:
        return 1.0
    return max(0.2, float(np.max(vals)))


def manifest_path(path: Path) -> str:
    root = Path.cwd().resolve()
    resolved = path.resolve()
    try:
        return resolved.relative_to(root).as_posix()
    except ValueError:
        return resolved.as_posix()


def write_plot_manifest(path: Path, csv_path: Path, out_path: Path, args: argparse.Namespace) -> None:
    payload: dict[str, object]
    if path.exists():
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            payload = {}
    else:
        payload = {}

    payload["format"] = "combined_meson_density_plot_manifest_v1"
    payload["date"] = dt.date.today().isoformat()
    payload["source_csv"] = manifest_path(csv_path)

    figures = payload.get("figures")
    if not isinstance(figures, list):
        figures = []
    out_rel = manifest_path(out_path)
    figures = [
        item for item in figures
        if not (isinstance(item, dict) and item.get("path") == out_rel)
    ]
    figures.append({
        "path": out_rel,
        "kind": "fig3_like_heatmap_png",
        "field": args.field,
        "dpi": args.dpi,
        "title": args.title,
        "generated_by": "scripts/analysis/relaxtime/render_combined_meson_density_fig3_like.py",
    })
    payload["figures"] = figures
    path.write_text(json.dumps(payload, ensure_ascii=True, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--manifest", default=None, type=Path)
    parser.add_argument("--field", default="kpi_ratio")
    parser.add_argument("--dpi", type=int, default=220)
    parser.add_argument("--vmax", type=float, default=None)
    parser.add_argument("--title", default="Combined meson-density FIG3-like scan")
    args = parser.parse_args()

    rows = read_rows(args.csv)
    regimes = sorted({row.get("regime", "") for row in rows if row.get("regime")})
    panels = []
    for regime in regimes:
        mus, temps, mat = matrix(rows, regime, args.field)
        if mus.size and temps.size:
            panels.append((regime, mus, temps, mat))
    if not panels:
        raise SystemExit("no plottable rows")

    vmax = args.vmax if args.vmax is not None else finite_max([panel[3] for panel in panels])
    ncols = 2
    nrows = int(math.ceil(len(panels) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(12.5, 4.8 * nrows), constrained_layout=True)
    axes_arr = np.atleast_1d(axes).ravel()
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("0.7")
    last_im = None
    for ax, (regime, mus, temps, mat) in zip(axes_arr, panels):
        xedges = centers_to_edges(mus)
        yedges = centers_to_edges(temps)
        masked = np.ma.masked_invalid(mat)
        last_im = ax.pcolormesh(xedges, yedges, masked, cmap=cmap, vmin=0.0, vmax=vmax, shading="auto")
        ax.set_title(regime)
        ax.set_xlabel("mu_q [MeV]")
        ax.set_ylabel("T [MeV]")
        ax.set_xlim(xedges[0], xedges[-1])
        ax.set_ylim(yedges[0], yedges[-1])
    for ax in axes_arr[len(panels):]:
        ax.axis("off")
    fig.suptitle(args.title)
    if last_im is not None:
        fig.colorbar(last_im, ax=axes_arr[:len(panels)], label=args.field)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=args.dpi)
    plt.close(fig)
    manifest = args.manifest if args.manifest is not None else args.out.parent / "plot_manifest.json"
    write_plot_manifest(manifest, args.csv, args.out, args)
    print(args.out)


if __name__ == "__main__":
    main()
