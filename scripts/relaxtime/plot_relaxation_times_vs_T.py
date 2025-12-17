#!/usr/bin/env python3
"""Plot relaxation times vs temperature.

Reads the CSV produced by scripts/relaxtime/scan_relaxation_times_vs_T.jl and
writes line plots for u, s, ubar, sbar at selected muB values.

Example:
  python scripts/relaxtime/plot_relaxation_times_vs_T.py \
    --csv data/outputs/results/relaxtime/relaxation_times_vs_T.csv \
    --muB 0,800 \
    --out-dir data/outputs/figures/relaxtime
"""

import argparse
import csv
import io
import math
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt


def parse_float(x: str) -> float:
    try:
        return float(x)
    except Exception:
        return math.nan


def read_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8") as f:
        lines = [line for line in f if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        return []
    return list(csv.DictReader(io.StringIO("".join(lines))))


def main() -> None:
    ap = argparse.ArgumentParser(description="Plot tau(T) at fixed muB")
    ap.add_argument(
        "--csv",
        type=Path,
        default=Path("data/outputs/results/relaxtime/relaxation_times_vs_T.csv"),
        help="Input CSV",
    )
    ap.add_argument(
        "--muB",
        type=str,
        default="0,800",
        help="Comma-separated muB values (MeV) to plot",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=Path("data/outputs/figures/relaxtime"),
        help="Output directory",
    )
    args = ap.parse_args()

    rows = read_rows(args.csv)
    if not rows:
        raise RuntimeError(f"No rows in {args.csv}")

    muB_vals = [parse_float(s.strip()) for s in args.muB.split(",") if s.strip()]
    muB_vals = [m for m in muB_vals if not math.isnan(m)]
    if not muB_vals:
        raise RuntimeError("--muB is empty")

    flavors = ["u", "s", "ubar", "sbar"]

    args.out_dir.mkdir(parents=True, exist_ok=True)

    for flavor in flavors:
        plt.figure(figsize=(6.8, 4.6))
        for muB in muB_vals:
            sub = [r for r in rows if abs(parse_float(r.get("muB_MeV", "nan")) - muB) < 1e-9]
            sub.sort(key=lambda r: parse_float(r.get("T_MeV", "nan")))
            xs = [parse_float(r.get("T_MeV", "nan")) for r in sub]
            ys = [parse_float(r.get(f"tau_{flavor}", "nan")) for r in sub]
            xs2, ys2 = zip(*[(x, y) for x, y in zip(xs, ys) if not (math.isnan(x) or math.isnan(y))]) if xs else ([], [])
            if not xs2:
                continue
            plt.plot(xs2, ys2, marker="o", lw=2, label=fr"$\mu_B={muB:.0f}\,\mathrm{{MeV}}$")

        plt.xlabel(r"$T$ (MeV)")
        plt.ylabel(fr"$\tau_{{{flavor}}}$ (fm)")
        plt.title(fr"Relaxation time $\tau_{{{flavor}}}(T)$")
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        out = args.out_dir / f"tau_{flavor}_vs_T.png"
        plt.savefig(out, dpi=200)
        plt.close()
        print(f"Saved {out}")


if __name__ == "__main__":
    main()
