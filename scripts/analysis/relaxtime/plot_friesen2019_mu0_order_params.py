#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt


PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INPUT = (
    PROJECT_ROOT
    / "data"
    / "outputs"
    / "results"
    / "relaxtime"
    / "meson_density"
    / "crossover_validation"
    / "friesen2019_mu0_order_params"
    / "mu0_order_params.csv"
)

LITERATURE_TC_MEV = 218.0


def configure_style() -> None:
    matplotlib.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times"],
            "font.size": 9,
            "mathtext.fontset": "stix",
            "axes.labelsize": 9,
            "axes.linewidth": 0.7,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "legend.frameon": False,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.06,
        }
    )


def read_rows(path: Path):
    rows = []
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append(
                {
                    "T_MeV": float(row["T_MeV"]),
                    "phi_u_norm": float(row["phi_u_norm"]),
                    "phi_s_norm": float(row["phi_s_norm"]),
                    "Phi": float(row["Phi"]),
                    "Phibar": float(row["Phibar"]),
                    "converged": row["converged"].strip().lower() == "true",
                }
            )
    rows.sort(key=lambda item: item["T_MeV"])
    return rows


def infer_current_tc_mev(rows) -> float:
    conv = [r for r in rows if r["converged"]]
    if len(conv) < 2:
        raise ValueError("need at least two converged points to infer crossover temperature")
    deriv = []
    for left, right in zip(conv[:-1], conv[1:]):
        dT = right["T_MeV"] - left["T_MeV"]
        if dT <= 0.0:
            continue
        slope = (right["phi_u_norm"] - left["phi_u_norm"]) / dT
        deriv.append(((left["T_MeV"] + right["T_MeV"]) * 0.5, slope))
    if not deriv:
        raise ValueError("failed to build derivative samples for crossover inference")
    return min(deriv, key=lambda item: item[1])[0]


def plot_curves(rows, out_path: Path) -> None:
    conv = [r for r in rows if r["converged"]]
    current_tc_mev = infer_current_tc_mev(rows)
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    ax.plot([r["T_MeV"] for r in conv], [r["phi_u_norm"] for r in conv], color="black", linewidth=1.8, label=r"$\langle \bar q q \rangle_u / \langle \bar q q \rangle_{u,T\approx 0}$")
    ax.plot([r["T_MeV"] for r in conv], [r["phi_s_norm"] for r in conv], color="#d62728", linewidth=1.6, label=r"$\langle \bar q q \rangle_s / \langle \bar q q \rangle_{s,T\approx 0}$")
    ax.plot([r["T_MeV"] for r in conv], [r["Phi"] for r in conv], color="#1f77b4", linewidth=1.6, label=r"$\Phi$")
    ax.axvline(current_tc_mev, color="#444444", linestyle="--", linewidth=1.0, label=r"current $T_c$")
    ax.axvline(LITERATURE_TC_MEV, color="#9467bd", linestyle=":", linewidth=1.2, label=r"literature $T_c$")
    ax.set_xlim(0.0, 350.0)
    ax.set_ylim(0.0, 1.05)
    ax.set_xlabel(r"$T$ [MeV]")
    ax.set_ylabel("normalized order parameters")
    ax.grid(True, alpha=0.2, linewidth=0.5)
    ax.legend(loc="best")
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def write_summary(rows, out_path: Path) -> None:
    conv = [r for r in rows if r["converged"]]
    current_tc_mev = infer_current_tc_mev(rows)
    def nearest(target: float):
        return min(conv, key=lambda r: abs(r["T_MeV"] - target))
    a = nearest(current_tc_mev)
    b = nearest(LITERATURE_TC_MEV)
    text = "\n".join(
        [
            "# Friesen 2019 muB=0 order-parameter quicklook",
            "",
            f"- converged points: {len(conv)}/{len(rows)}",
            f"- current crossover reference: {current_tc_mev:.6f} MeV",
            f"- literature crossover reference: {LITERATURE_TC_MEV:.6f} MeV",
            f"- near current Tc: T={a['T_MeV']:.6f} MeV, phi_u_norm={a['phi_u_norm']:.6f}, phi_s_norm={a['phi_s_norm']:.6f}, Phi={a['Phi']:.6f}",
            f"- near literature Tc: T={b['T_MeV']:.6f} MeV, phi_u_norm={b['phi_u_norm']:.6f}, phi_s_norm={b['phi_s_norm']:.6f}, Phi={b['Phi']:.6f}",
        ]
    )
    out_path.write_text(text + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Plot muB=0 order-parameter curves for Friesen-2019 quicklook.")
    ap.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    configure_style()
    rows = read_rows(args.input)
    out_dir = args.input.parent
    plot_curves(rows, out_dir / "mu0_order_params_quicklook.png")
    write_summary(rows, out_dir / "quicklook_summary.md")


if __name__ == "__main__":
    main()
