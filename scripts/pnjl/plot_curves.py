#!/usr/bin/env python3
"""Plot PNJL ρ-μ and P-μ curves from scan CSV files."""

import argparse
import csv
import math
from pathlib import Path
from typing import Iterable, List, Dict

import matplotlib.pyplot as plt

RESULT_DIR = Path("data/outputs/figures/pnjl")
RESULT_DIR.mkdir(parents=True, exist_ok=True)


def read_csv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"CSV not found: {path}")
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        return [row for row in reader]


def parse_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def filter_rows(rows: Iterable[Dict[str, str]], *, T: float, xi: float, tol: float = 1e-6) -> List[Dict[str, str]]:
    filtered = []
    for row in rows:
        if not row:
            continue
        try:
            T_val = parse_float(row.get("T_MeV", "nan"))
            xi_val = parse_float(row.get("xi", "nan"))
        except ValueError:
            continue
        if math.isnan(T_val) or math.isnan(xi_val):
            continue
        if abs(T_val - T) < tol and abs(xi_val - xi) < tol and row.get("converged", "").lower() == "true":
            filtered.append(row)
    filtered.sort(key=lambda r: parse_float(r.get("mu_MeV", r.get("mu_avg_MeV", "nan"))))
    return filtered


def plot_rho_mu(trho_rows: List[Dict[str, str]], temps: List[float], xi: float, output: Path) -> None:
    if not trho_rows:
        return
    plt.figure(figsize=(8, 5))
    for T in temps:
        rows = filter_rows(trho_rows, T=T, xi=xi)
        if not rows:
            continue
        mu_vals = [parse_float(r.get("mu_avg_MeV", "nan")) for r in rows]
        rho_vals = [parse_float(r.get("rho", "nan")) for r in rows]
        plt.plot(mu_vals, rho_vals, marker="o", label=f"T={T:.0f} MeV")
    plt.xlabel(r"$\mu$ (MeV)")
    plt.ylabel(r"$\rho/\rho_0$")
    plt.title(f"PNJL ρ-μ curves (ξ={xi})")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output, dpi=200)
    plt.close()


def plot_P_mu(tmu_rows: List[Dict[str, str]], temps: List[float], xi: float, output: Path) -> None:
    if not tmu_rows:
        return
    plt.figure(figsize=(8, 5))
    for T in temps:
        rows = filter_rows(tmu_rows, T=T, xi=xi)
        if not rows:
            continue
        mu_vals = [parse_float(r.get("mu_MeV", "nan")) for r in rows]
        pressure_vals = [parse_float(r.get("pressure_fm4", "nan")) for r in rows]
        plt.plot(mu_vals, pressure_vals, marker="o", label=f"T={T:.0f} MeV")
    plt.xlabel(r"$\mu$ (MeV)")
    plt.ylabel(r"Pressure (fm$^{-4}$)")
    plt.title(f"PNJL P-μ curves (ξ={xi})")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output, dpi=200)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Plot PNJL scan curves")
    parser.add_argument("--tmu-csv", type=Path, default=Path("data/outputs/results/pnjl/tmu_scan.csv"), help="Path to T-μ scan CSV")
    parser.add_argument("--trho-csv", type=Path, default=Path("data/outputs/results/pnjl/trho_scan.csv"), help="Path to T-ρ scan CSV")
    parser.add_argument("--temps", type=float, nargs="*", default=[50.0, 100.0, 150.0], help="List of temperatures to plot")
    parser.add_argument("--xi", type=float, default=0.0, help="Anisotropy parameter ξ to filter")
    parser.add_argument("--rho-output", type=Path, default=RESULT_DIR / "pnjl_rho_mu.png", help="Output path for ρ-μ plot")
    parser.add_argument("--pressure-output", type=Path, default=RESULT_DIR / "pnjl_pressure_mu.png", help="Output path for P-μ plot")
    args = parser.parse_args()

    tmu_rows = read_csv(args.tmu_csv) if args.tmu_csv.exists() else []
    trho_rows = read_csv(args.trho_csv) if args.trho_csv.exists() else []

    if trho_rows:
        plot_rho_mu(trho_rows, args.temps, args.xi, args.rho_output)
        print(f"Saved ρ-μ plot to {args.rho_output}")
    else:
        print(f"Skipping ρ-μ plot: CSV not found ({args.trho_csv})")

    if tmu_rows:
        plot_P_mu(tmu_rows, args.temps, args.xi, args.pressure_output)
        print(f"Saved P-μ plot to {args.pressure_output}")
    else:
        print(f"Skipping P-μ plot: CSV not found ({args.tmu_csv})")


if __name__ == "__main__":
    main()
