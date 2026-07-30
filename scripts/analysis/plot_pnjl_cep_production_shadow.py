#!/usr/bin/env python3
"""Generate compact diagnostic figures from shadow aggregate tables."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        (args.output_dir / "plot_manifest.json").write_text(json.dumps({"schema_version": "cep_cascade_production_shadow_v1", "figures": [], "reason": "matplotlib unavailable"}) + "\n", encoding="utf-8")
        return 0
    curve_rows = rows(args.input_dir / "curve_points.csv")
    figures = []
    for xi in ("-0.5", "0.0", "0.5"):
        subset = [row for row in curve_rows if row.get("xi") == xi]
        if not subset:
            continue
        fig, ax = plt.subplots(figsize=(7, 4))
        for method in ("production_cascade", "memoized_dense", "independent_oracle"):
            points = [row for row in subset if row.get("method") == method and row.get("sampling_role") != "targeted"]
            points.sort(key=lambda row: (float(row["T_MeV"]), float(row["rho"])))
            if points:
                ax.plot([float(row["rho"]) for row in points], [float(row["muq_MeV"]) for row in points], ".", ms=2, label=method)
        ax.set_xlabel(r"$\rho$")
        ax.set_ylabel(r"$\mu_q$ [MeV]")
        ax.set_title(f"CEP production shadow, xi={xi}")
        ax.legend(fontsize=7)
        path = args.output_dir / f"rho_mu_xi_{xi.replace('.', 'p').replace('-', 'm')}.png"
        fig.tight_layout()
        fig.savefig(path, dpi=140)
        plt.close(fig)
        figures.append(str(path.name))
    (args.output_dir / "plot_manifest.json").write_text(
        json.dumps({"schema_version": "cep_cascade_production_shadow_v1", "figures": figures}, indent=2) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
