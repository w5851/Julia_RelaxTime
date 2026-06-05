#!/usr/bin/env python
"""Render temperature-scan PNGs from combined meson-density scan CSV output."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt


FIELDS = ("n_pi", "n_K", "kpi_ratio")
FIELD_TITLES = {
    "n_pi": "n_pi vs T",
    "n_K": "n_K vs T",
    "kpi_ratio": "K/pi ratio vs T",
}


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


def manifest_path(path: Path) -> str:
    root = Path.cwd().resolve()
    resolved = path.resolve()
    try:
        return resolved.relative_to(root).as_posix()
    except ValueError:
        return resolved.as_posix()


def write_plot_manifest(path: Path, csv_path: Path, out_path: Path, args: argparse.Namespace) -> None:
    if path.exists():
        try:
            payload: dict[str, object] = json.loads(path.read_text(encoding="utf-8"))
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
        "kind": "temperature_scan_png",
        "fields": list(FIELDS),
        "dpi": args.dpi,
        "title": args.title,
        "generated_by": "scripts/analysis/relaxtime/render_combined_meson_density_temperature_scan.py",
    })
    payload["figures"] = figures
    path.write_text(json.dumps(payload, ensure_ascii=True, indent=2) + "\n", encoding="utf-8")


def update_readme(path: Path, out_path: Path, manifest_path: Path, args: argparse.Namespace) -> None:
    if not path.exists():
        return
    text = path.read_text(encoding="utf-8")
    png_line = f"- PNG: `{manifest_path_to_readme(out_path)}`"
    command_line = (
        "- PNG command: "
        f"`python scripts/analysis/relaxtime/render_combined_meson_density_temperature_scan.py "
        f"--csv {manifest_path_to_readme(args.csv)} "
        f"--out {manifest_path_to_readme(out_path)} "
        f"--manifest {manifest_path_to_readme(manifest_path)} "
        f"--readme {manifest_path_to_readme(path)} "
        f"--dpi {args.dpi} --title \"{args.title}\"`"
    )

    lines = text.splitlines()
    png_idx = next((idx for idx, line in enumerate(lines) if line.strip().startswith("- PNG:")), None)
    if png_idx is None:
        for idx, line in enumerate(lines):
            if line.strip().startswith("- plot manifest:"):
                lines.insert(idx, png_line)
                break
        else:
            lines.append(png_line)
    else:
        lines[png_idx] = png_line

    for idx in reversed([
        i for i, line in enumerate(lines)
        if line.strip().startswith("- PNG command:")
    ]):
        del lines[idx]

    render_idx = next((idx for idx, line in enumerate(lines) if line.strip() == "## Rendering"), None)
    if render_idx is None:
        insert_idx = next((idx for idx, line in enumerate(lines) if line.strip() == "## Regime Definitions"), len(lines))
        lines.insert(insert_idx, "")
        lines.insert(insert_idx, command_line)
        lines.insert(insert_idx, "")
        lines.insert(insert_idx, "## Rendering")
    else:
        insert_idx = render_idx + 1
        while insert_idx < len(lines) and lines[insert_idx].strip() == "":
            insert_idx += 1
        lines.insert(insert_idx, command_line)
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def manifest_path_to_readme(path: Path) -> str:
    return manifest_path(path)


def plot_temperature_scan(rows: list[dict[str, str]], out_path: Path, args: argparse.Namespace) -> None:
    regimes = sorted({row.get("regime", "") for row in rows if row.get("regime")})
    if not regimes:
        raise SystemExit("no regimes found")

    fig, axes = plt.subplots(2, 2, figsize=(13.0, 8.4), constrained_layout=True)
    axes_arr = axes.ravel()
    panel_axes = axes_arr[:3]
    legend_ax = axes_arr[3]

    for ax, field in zip(panel_axes, FIELDS):
        plotted = False
        for regime in regimes:
            selected = [
                row for row in rows
                if row.get("regime") == regime and row.get("status") == "ok"
            ]
            pairs = []
            for row in selected:
                temp = as_float(row, "T_MeV")
                value = as_float(row, field)
                if math.isfinite(temp) and math.isfinite(value) and value > 0.0:
                    pairs.append((temp, value))
            if not pairs:
                continue
            pairs.sort()
            xs = [p[0] for p in pairs]
            ys = [p[1] for p in pairs]
            ax.plot(xs, ys, marker="o", linewidth=1.8, markersize=4.2, label=regime)
            plotted = True
        ax.set_title(FIELD_TITLES[field])
        ax.set_xlabel("T [MeV]")
        ax.set_ylabel(field)
        if field in ("n_pi", "n_K"):
            ax.set_yscale("log")
        ax.grid(True, alpha=0.25)
        if not plotted:
            ax.text(0.5, 0.5, "no finite rows", ha="center", va="center", transform=ax.transAxes)

    legend_ax.axis("off")
    handles, labels = panel_axes[0].get_legend_handles_labels()
    if handles:
        legend_ax.legend(handles, labels, loc="upper left", frameon=False, title="Regime")
    legend_ax.text(
        0.0,
        0.02,
        f"source: {args.csv.name}",
        fontsize=9,
        color="0.35",
        transform=legend_ax.transAxes,
    )

    fig.suptitle(args.title)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=args.dpi)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--manifest", default=None, type=Path)
    parser.add_argument("--readme", default=None, type=Path)
    parser.add_argument("--dpi", type=int, default=260)
    parser.add_argument("--title", default="Combined meson-density temperature scan")
    args = parser.parse_args()

    rows = read_rows(args.csv)
    plot_temperature_scan(rows, args.out, args)
    manifest = args.manifest if args.manifest is not None else args.out.parent / "plot_manifest.json"
    write_plot_manifest(manifest, args.csv, args.out, args)
    if args.readme is not None:
        update_readme(args.readme, args.out, manifest, args)
    print(args.out)


if __name__ == "__main__":
    main()
