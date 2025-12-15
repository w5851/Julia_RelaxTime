#!/usr/bin/env python3
"""Plot heatmaps from gap_transport_scan.csv.

Reads the CSV produced by scripts/relaxtime/run_gap_transport_scan.jl and
renders (T, mu) heatmaps for selected fields at a fixed xi.

Example:
  python scripts/relaxtime/plot_gap_transport_scan.py \
    --csv data/outputs/results/relaxtime/gap_transport_scan.csv \
    --xi 0.0 \
    --fields eta,sigma,tau_u \
    --out-dir data/outputs/figures/relaxtime
"""

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt


def parse_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def read_rows(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"CSV not found: {path}")
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        return [row for row in reader]


def filter_by_xi(rows: List[Dict[str, str]], xi: float, tol: float = 1e-9) -> List[Dict[str, str]]:
    out = []
    for r in rows:
        xi_val = parse_float(r.get("xi", "nan"))
        if math.isnan(xi_val) or abs(xi_val - xi) > tol:
            continue
        out.append(r)
    return out


def build_grid(rows: List[Dict[str, str]], field: str) -> Tuple[List[float], List[float], List[List[float]]]:
    t_vals = sorted({parse_float(r.get("T_MeV", "nan")) for r in rows if r})
    mu_vals = sorted({parse_float(r.get("mu_MeV", "nan")) for r in rows if r})
    t_vals = [v for v in t_vals if not math.isnan(v)]
    mu_vals = [v for v in mu_vals if not math.isnan(v)]

    t_index = {t: i for i, t in enumerate(t_vals)}
    mu_index = {m: j for j, m in enumerate(mu_vals)}

    grid = [[math.nan for _ in mu_vals] for _ in t_vals]
    for r in rows:
        t = parse_float(r.get("T_MeV", "nan"))
        m = parse_float(r.get("mu_MeV", "nan"))
        if math.isnan(t) or math.isnan(m):
            continue
        i = t_index.get(t)
        j = mu_index.get(m)
        if i is None or j is None:
            continue
        grid[i][j] = parse_float(r.get(field, "nan"))
    return t_vals, mu_vals, grid


def plot_heatmap(t_vals: List[float], mu_vals: List[float], grid: List[List[float]], *, field: str, xi: float, out_path: Path) -> None:
    if not t_vals or not mu_vals:
        return

    plt.figure(figsize=(7.6, 5.2))
    x0, x1 = min(mu_vals), max(mu_vals)
    y0, y1 = min(t_vals), max(t_vals)
    if x0 == x1:
        x0 -= 0.5
        x1 += 0.5
    if y0 == y1:
        y0 -= 0.5
        y1 += 0.5

    im = plt.imshow(
        grid,
        origin="lower",
        aspect="auto",
        extent=[x0, x1, y0, y1],
        interpolation="nearest",
    )
    plt.colorbar(im, label=field)
    plt.xlabel(r"$\mu$ (MeV)")
    plt.ylabel(r"$T$ (MeV)")
    plt.title(f"{field} heatmap (xi={xi})")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=200)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot heatmaps from gap_transport_scan CSV")
    parser.add_argument("--csv", type=Path, default=Path("data/outputs/results/relaxtime/gap_transport_scan.csv"), help="Input CSV")
    parser.add_argument("--xi", type=float, default=0.0, help="Filter xi")
    parser.add_argument("--fields", type=str, default="eta,sigma,tau_u", help="Comma-separated fields to plot")
    parser.add_argument("--out-dir", type=Path, default=Path("data/outputs/figures/relaxtime"), help="Output directory")
    args = parser.parse_args()

    rows = read_rows(args.csv)
    rows = filter_by_xi(rows, args.xi)
    if not rows:
        raise RuntimeError(f"No rows for xi={args.xi} in {args.csv}")

    fields = [f.strip() for f in args.fields.split(",") if f.strip()]
    for field in fields:
        t_vals, mu_vals, grid = build_grid(rows, field)
        out_path = args.out_dir / f"gap_transport_{field}_xi{args.xi}.png"
        plot_heatmap(t_vals, mu_vals, grid, field=field, xi=args.xi, out_path=out_path)
        print(f"Saved {field} heatmap to {out_path}")


if __name__ == "__main__":
    main()
