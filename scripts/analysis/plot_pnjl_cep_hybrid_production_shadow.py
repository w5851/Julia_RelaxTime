#!/usr/bin/env python3
"""Generate deterministic full-curve/support-zoom figures for hybrid shadow."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

METHODS = ("production_hybrid", "memoized_dense", "independent_oracle")
XIS = ("-0.5", "0.0", "0.5")


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _float(value: str | None) -> float:
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return math.nan


def _status(rows: list[dict[str, str]], xi: str, temperature: float, method: str) -> str:
    for row in rows:
        if row.get("method") == method and row.get("xi") == xi and math.isclose(_float(row.get("T_MeV")), temperature, abs_tol=1e-8, rel_tol=0.0):
            return row.get("result_status", "")
    return "missing"


def _select_anchors(rows: list[dict[str, str]]) -> list[tuple[str, float, str]]:
    selected: dict[tuple[str, float], str] = {}
    for xi in XIS:
        oracle = sorted([row for row in rows if row.get("xi") == xi and row.get("method") == "independent_oracle"], key=lambda row: _float(row.get("T_MeV")))
        if not oracle:
            continue
        selected[(xi, _float(oracle[0].get("T_MeV")))] = "low_temperature"
        first = [row for row in oracle if row.get("result_status") == "confirmed_first_order"]
        monotone = [row for row in oracle if row.get("result_status") == "confirmed_monotone"]
        ambiguous = [row for row in oracle if row.get("result_status") == "ambiguous_near_critical"]
        if first:
            selected[(xi, _float(first[-1].get("T_MeV")))] = "first_order"
        if monotone:
            selected[(xi, _float(monotone[0].get("T_MeV")))] = "first_monotone"
        if ambiguous:
            selected[(xi, _float(ambiguous[-1].get("T_MeV")))] = "ambiguous_near_critical"
        for oracle_row in oracle:
            temperature = _float(oracle_row.get("T_MeV"))
            hybrid = _status(rows, xi, temperature, "production_hybrid")
            oracle_status = oracle_row.get("result_status", "")
            if oracle_status == "confirmed_first_order" and hybrid == "ambiguous_near_critical":
                selected[(xi, temperature)] = "oracle_first_order_hybrid_ambiguous"
            elif oracle_status == "ambiguous_near_critical" and hybrid == "confirmed_first_order":
                selected[(xi, temperature)] = "unsupported_hybrid_confirmation"
    return sorted([(xi, temperature, reason) for (xi, temperature), reason in selected.items()], key=lambda item: (XIS.index(item[0]), item[1], item[2]))


def _bounds(curve_rows: list[dict[str, str]], slice_rows: list[dict[str, str]], xi: str, temperature: float) -> tuple[float, float]:
    values: list[float] = []
    for row in slice_rows:
        if row.get("xi") == xi and row.get("method") == "production_hybrid" and math.isclose(_float(row.get("T_MeV")), temperature, abs_tol=1e-8, rel_tol=0.0):
            for field in ("rho_hadron", "rho_quark", "rho_spinodal_hadron", "rho_spinodal_quark", "support_low", "support_high"):
                value = _float(row.get(field))
                if math.isfinite(value):
                    values.append(value)
    if len(values) < 2:
        values = [_float(row.get("rho")) for row in curve_rows if math.isfinite(_float(row.get("rho")))]
    if len(values) < 2:
        return 0.0, 4.0
    low, high = min(values), max(values)
    padding = max((high - low) * 0.15, 0.05)
    return max(0.0, low - padding), min(4.0, high + padding)


def _plot_anchor(output: Path, curve_rows: list[dict[str, str]], slice_rows: list[dict[str, str]], xi: str, temperature: float, reason: str) -> str:
    import matplotlib.pyplot as plt

    subset = [row for row in curve_rows if row.get("xi") == xi and math.isclose(_float(row.get("T_MeV")), temperature, abs_tol=1e-8, rel_tol=0.0)]
    low, high = _bounds(subset, slice_rows, xi, temperature)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    colors = {"production_hybrid": "tab:orange", "memoized_dense": "tab:blue", "independent_oracle": "tab:green"}
    for method in METHODS:
        points = sorted([row for row in subset if row.get("method") == method], key=lambda row: _float(row.get("rho")))
        finite = [(_float(row.get("rho")), _float(row.get("muq_MeV"))) for row in points]
        finite = [(x, y) for x, y in finite if math.isfinite(x) and math.isfinite(y)]
        if not finite:
            continue
        x, y = zip(*finite)
        label = f"{method} ({_status(slice_rows, xi, temperature, method)})"
        axes[0].plot(x, y, ".-", ms=2, lw=0.7, color=colors[method], label=label)
        zoom = [(xx, yy) for xx, yy in finite if low <= xx <= high]
        if zoom:
            zx, zy = zip(*zoom)
            axes[1].plot(zx, zy, ".-", ms=3, lw=0.9, color=colors[method], label=method)
    axes[0].set_xlim(0.0, 4.0)
    axes[1].set_xlim(low, high)
    for axis in axes:
        axis.set_xlabel(r"$\rho$")
        axis.grid(alpha=0.2)
    axes[0].set_ylabel(r"$\mu_q$ [MeV]")
    if axes[0].lines:
        axes[0].legend(fontsize=7, loc="best")
    if axes[1].lines:
        axes[1].legend(fontsize=7, loc="best")
    fig.suptitle(f"CEP hybrid shadow: xi={xi}, T={temperature:g} MeV ({reason})")
    fig.tight_layout()
    filename = f"rho_mu_xi_{xi.replace('.', 'p').replace('-', 'm')}_T_{temperature:g}_{reason}.png"
    fig.savefig(output / filename, dpi=150)
    plt.close(fig)
    return filename


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib

        matplotlib.use("Agg")
    except ImportError:
        (args.output_dir / "plot_manifest.json").write_text(json.dumps({"schema_version": "cep_cascade_production_shadow_v2", "figures": [], "reason": "matplotlib unavailable"}, indent=2) + "\n", encoding="utf-8")
        return 0
    curve_path = args.input_dir / "curve_points.csv"
    slice_path = args.input_dir / "slice_metrics.csv"
    curve_rows, slice_rows = _rows(curve_path), _rows(slice_path)
    figures = []
    for xi, temperature, reason in _select_anchors(slice_rows):
        filename = _plot_anchor(args.output_dir, curve_rows, slice_rows, xi, temperature, reason)
        figures.append({"file": filename, "xi": xi, "T_MeV": temperature, "reason": reason})
    manifest = {
        "schema_version": "cep_cascade_production_shadow_v2",
        "source_sha256": {"curve_points.csv": hashlib.sha256(curve_path.read_bytes()).hexdigest(), "slice_metrics.csv": hashlib.sha256(slice_path.read_bytes()).hexdigest()},
        "figures": figures,
    }
    for figure in figures:
        figure["sha256"] = hashlib.sha256((args.output_dir / figure["file"]).read_bytes()).hexdigest()
    (args.output_dir / "plot_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
