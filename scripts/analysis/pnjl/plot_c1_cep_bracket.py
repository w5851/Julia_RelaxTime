#!/usr/bin/env python3
"""Render the C1 CEP temperature brackets as a solver-free diagnostic plot.

The C1 CEP artifact contains a temperature bracket for each xi, but no
resolved ``T_CEP_MeV`` value.  This script deliberately renders the bracket
as a band and labels its midpoint as diagnostic only; it never promotes the
midpoint to a CEP result.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REQUIRED_COLUMNS = (
    "xi",
    "T_bracket_low_MeV",
    "T_bracket_high_MeV",
    "bracket_width_T_MeV",
    "result_status",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cep", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--metadata", type=Path)
    return parser.parse_args()


def fail(message: str) -> None:
    raise SystemExit(f"[c1-cep-bracket] {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def finite(row: dict[str, str], field: str, row_number: int) -> float:
    try:
        value = float(row[field])
    except (KeyError, TypeError, ValueError) as exc:
        fail(f"row {row_number} has non-numeric {field}: {row.get(field)!r}")
        raise AssertionError from exc
    if not math.isfinite(value):
        fail(f"row {row_number} has non-finite {field}: {value}")
    return value


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        fail(f"missing input: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = sorted(set(REQUIRED_COLUMNS) - set(reader.fieldnames or []))
        if missing:
            fail(f"{path.name} missing columns: {missing}")
        rows = list(reader)

    if not rows:
        fail("CEP input is empty")

    seen_xi: set[float] = set()
    parsed: list[dict[str, float | str]] = []
    for row_number, row in enumerate(rows, 2):
        xi = finite(row, "xi", row_number)
        low = finite(row, "T_bracket_low_MeV", row_number)
        high = finite(row, "T_bracket_high_MeV", row_number)
        width = finite(row, "bracket_width_T_MeV", row_number)
        if high < low:
            fail(f"row {row_number} has reversed bracket: {low}, {high}")
        if width <= 0 or not math.isclose(width, high - low, rel_tol=0.0, abs_tol=1e-10):
            fail(f"row {row_number} has inconsistent bracket width: {width} != {high - low}")
        key = round(xi, 10)
        if key in seen_xi:
            fail(f"duplicate xi key at row {row_number}: {key}")
        seen_xi.add(key)
        parsed.append(
            {
                "xi": xi,
                "low": low,
                "high": high,
                "width": width,
                "mid": 0.5 * (low + high),
                "status": row.get("result_status", "").strip() or "missing",
            }
        )
    return sorted(parsed, key=lambda row: float(row["xi"]))


def render(rows: list[dict[str, float | str]], output: Path, input_path: Path) -> dict[str, object]:
    xi = np.asarray([float(row["xi"]) for row in rows])
    low = np.asarray([float(row["low"]) for row in rows])
    high = np.asarray([float(row["high"]) for row in rows])
    midpoint = np.asarray([float(row["mid"]) for row in rows])

    output.parent.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(9.2, 5.4), constrained_layout=True)
    axis.fill_between(
        xi,
        low,
        high,
        color="#4c78a8",
        alpha=0.28,
        linewidth=0,
        label=r"CEP temperature bracket $[T_{\rm low},T_{\rm high}]$",
    )
    axis.plot(xi, low, color="#2f5597", linewidth=1.0, alpha=0.9, label=r"$T_{\rm low}$")
    axis.plot(xi, high, color="#2f5597", linewidth=1.0, alpha=0.9, label=r"$T_{\rm high}$")
    axis.plot(
        xi,
        midpoint,
        color="#222222",
        linestyle=(0, (4, 2)),
        linewidth=1.05,
        marker=".",
        markersize=3.5,
        label="bracket midpoint (diagnostic only)",
    )
    axis.set_xlabel(r"$\xi$")
    axis.set_ylabel(r"$T$ (MeV)")
    axis.set_title("C1 CEP temperature brackets versus $\\xi$")
    axis.grid(True, color="#e1e1e1", linewidth=0.65)
    axis.minorticks_on()
    axis.tick_params(direction="in", top=True, right=True, which="both")
    axis.legend(loc="best", framealpha=0.92)
    status_counts: dict[str, int] = {}
    for row in rows:
        status = str(row["status"])
        status_counts[status] = status_counts.get(status, 0) + 1
    figure.text(
        0.01,
        0.005,
        f"diagnostic-only | rows={len(rows)} | statuses={status_counts} | no resolved CEP promoted",
        fontsize=8,
        color="#444444",
    )
    figure.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(figure)
    return {
        "path": str(output),
        "sha256": sha256(output),
        "input": str(input_path),
        "input_sha256": sha256(input_path),
        "row_count": len(rows),
        "xi_min": float(xi.min()),
        "xi_max": float(xi.max()),
        "T_min_MeV": float(low.min()),
        "T_max_MeV": float(high.max()),
        "bracket_widths_MeV": sorted({float(row["width"]) for row in rows}),
        "status_counts": status_counts,
        "midpoint_is_physical_CEP": False,
    }


def main() -> None:
    args = parse_args()
    rows = read_rows(args.cep)
    figure_info = render(rows, args.output, args.cep)
    if args.metadata is not None:
        args.metadata.parent.mkdir(parents=True, exist_ok=True)
        args.metadata.write_text(
            json.dumps(
                {
                    "schema_version": "pnjl_c1_cep_bracket_plot_v1",
                    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
                    "figure": figure_info,
                    "plot_contract": {
                        "x": "xi",
                        "y": "T_bracket_low_MeV..T_bracket_high_MeV",
                        "midpoint": "diagnostic_only",
                        "solver_called": False,
                    },
                },
                indent=2,
                ensure_ascii=False,
            )
            + "\n",
            encoding="utf-8",
        )
    print(json.dumps(figure_info, ensure_ascii=False))


if __name__ == "__main__":
    main()
