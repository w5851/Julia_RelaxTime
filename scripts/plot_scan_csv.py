#!/usr/bin/env python3
"""Generic plotter for scan CSV files (with optional metadata header).

CSV spec (scan_csv_v1):
- Lines starting with '#' are metadata/comments and ignored by the CSV parser.
- The first non-comment line is the CSV header.

Supports two modes:
- lines: y(x) curves, optionally grouped by a column
- heatmap: z(x,y) heatmaps, optionally filtered by column=value

Examples
--------
Relaxation times vs T (two curves for muB):
  python scripts/plot_scan_csv.py \
    --mode lines \
    --csv data/outputs/results/relaxtime/relaxation_times_vs_T.csv \
    --x T_MeV --ys tau_u,tau_s,tau_ubar,tau_sbar \
    --group muB_MeV \
    --out-dir data/outputs/figures/relaxtime

Gap+transport heatmaps at fixed xi:
  python scripts/plot_scan_csv.py \
    --mode heatmap \
    --csv data/outputs/results/relaxtime/gap_transport_scan.csv \
    --x muq_MeV --y T_MeV --fields eta,sigma,tau_u \
    --where xi=0.0 \
    --out-dir data/outputs/figures/relaxtime
"""

from __future__ import annotations

import argparse
import csv
import io
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt


def _find_project_root() -> Path:
    """Find project root by looking for Project.toml or .git."""
    # Start from script location
    script_dir = Path(__file__).resolve().parent
    candidates = [script_dir, script_dir.parent]
    # Also check cwd
    candidates.append(Path.cwd())
    
    for start in candidates:
        current = start
        for _ in range(5):  # max 5 levels up
            if (current / "Project.toml").exists() or (current / ".git").exists():
                return current
            parent = current.parent
            if parent == current:
                break
            current = parent
    # Fallback to cwd
    return Path.cwd()


PROJECT_ROOT = _find_project_root()
DEFAULT_OUT_DIR = PROJECT_ROOT / "data" / "outputs" / "figures"


def _parse_float(x: str) -> float:
    try:
        return float(x)
    except Exception:
        return math.nan


def _is_comment(line: str) -> bool:
    s = line.strip()
    return (not s) or s.startswith("#")


def read_scan_csv(path: Path) -> Tuple[Dict[str, str], List[Dict[str, str]]]:
    if not path.exists():
        raise FileNotFoundError(f"CSV not found: {path}")

    meta: Dict[str, str] = {}
    data_lines: List[str] = []

    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if _is_comment(line):
                s = line.strip()
                if s.startswith("#"):
                    s2 = s[1:].strip()
                    if ":" in s2:
                        k, v = s2.split(":", 1)
                        meta[k.strip()] = v.strip()
                    elif "=" in s2:
                        k, v = s2.split("=", 1)
                        meta[k.strip()] = v.strip()
                continue
            data_lines.append(line)

    if not data_lines:
        return meta, []

    reader = csv.DictReader(io.StringIO("".join(data_lines)))
    return meta, [row for row in reader]


def _meta_get(meta: Dict[str, str], key: str, default: str | None = None) -> str | None:
    v = meta.get(key)
    if v is None:
        return default
    s = str(v).strip()
    return s if s else default


def _axis_label(meta: Dict[str, str], *, axis: str, col: str) -> str:
    # Supports:
    # - x_label / y_label
    # - x_unit / y_unit
    # - y_label.<col> / y_unit.<col>
    label = _meta_get(meta, f"{axis}_label.{col}") or _meta_get(meta, f"{axis}_label") or col
    unit = _meta_get(meta, f"{axis}_unit.{col}") or _meta_get(meta, f"{axis}_unit")
    return f"{label} [{unit}]" if unit else label


def _axis_scale(meta: Dict[str, str], *, axis: str, col: str) -> str | None:
    # Supports:
    # - x_scale / y_scale
    # - y_scale.<col>
    s = _meta_get(meta, f"{axis}_scale.{col}") or _meta_get(meta, f"{axis}_scale")
    if s is None:
        return None
    s2 = s.lower()
    if s2 in {"linear", "log"}:
        return s2
    return None


def _apply_where(rows: List[Dict[str, str]], where: List[str]) -> List[Dict[str, str]]:
    if not where:
        return rows

    clauses: List[Tuple[str, str]] = []
    for w in where:
        if "=" not in w:
            raise ValueError(f"Invalid --where clause (expected col=value): {w}")
        k, v = w.split("=", 1)
        clauses.append((k.strip(), v.strip()))

    out: List[Dict[str, str]] = []
    for r in rows:
        ok = True
        for k, v in clauses:
            if str(r.get(k, "")) != v:
                ok = False
                break
        if ok:
            out.append(r)
    return out


def _sanitize_filename(s: str) -> str:
    s2 = "".join(ch if ch.isalnum() or ch in {"-", "_", "."} else "_" for ch in str(s))
    s2 = s2.strip("._")
    return s2 if s2 else "value"


def _split_rows(rows: List[Dict[str, str]], *, split: str | None) -> List[Tuple[str, List[Dict[str, str]]]]:
    if not split:
        return [("__all__", rows)]

    groups: Dict[str, List[Dict[str, str]]] = {}
    for r in rows:
        groups.setdefault(str(r.get(split, "")), []).append(r)
    return sorted(groups.items(), key=lambda kv: kv[0])


def plot_lines(
    rows: List[Dict[str, str]],
    *,
    x: str,
    ys: List[str],
    group: str | None,
    out_dir: Path,
    title_prefix: str | None,
    meta: Dict[str, str] | None,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    def group_key(r: Dict[str, str]) -> str:
        return "__all__" if not group else str(r.get(group, ""))

    groups: Dict[str, List[Dict[str, str]]] = {}
    for r in rows:
        groups.setdefault(group_key(r), []).append(r)

    for y in ys:
        plt.figure(figsize=(6.8, 4.6))
        any_plotted = False
        for gk, sub in sorted(groups.items(), key=lambda kv: kv[0]):
            sub2 = sorted(sub, key=lambda r: _parse_float(r.get(x, "nan")))
            xs = [_parse_float(r.get(x, "nan")) for r in sub2]
            ys2 = [_parse_float(r.get(y, "nan")) for r in sub2]
            pairs = [(xx, yy) for xx, yy in zip(xs, ys2) if not (math.isnan(xx) or math.isnan(yy))]
            if not pairs:
                continue
            xs3, ys3 = zip(*pairs)
            label = gk if group else y
            plt.plot(xs3, ys3, marker="o", lw=2, label=label)
            any_plotted = True

        if not any_plotted:
            plt.close()
            continue

        if meta:
            plt.xlabel(_axis_label(meta, axis="x", col=x))
            plt.ylabel(_axis_label(meta, axis="y", col=y))
            xscale = _axis_scale(meta, axis="x", col=x)
            yscale = _axis_scale(meta, axis="y", col=y)
            if xscale:
                plt.xscale(xscale)
            if yscale:
                plt.yscale(yscale)
        else:
            plt.xlabel(x)
            plt.ylabel(y)
        title = f"{title_prefix} - {y}" if title_prefix else y
        plt.title(title)
        plt.grid(True, alpha=0.3)
        if group:
            plt.legend()
        plt.tight_layout()
        out = out_dir / f"{y}_vs_{x}.png"
        plt.savefig(out, dpi=200)
        plt.close()
        print(f"Saved {out}")


def _build_grid(rows: List[Dict[str, str]], *, x: str, y: str, z: str) -> Tuple[List[float], List[float], List[List[float]]]:
    xs = sorted({_parse_float(r.get(x, "nan")) for r in rows})
    ys = sorted({_parse_float(r.get(y, "nan")) for r in rows})
    xs = [v for v in xs if not math.isnan(v)]
    ys = [v for v in ys if not math.isnan(v)]

    x_index = {v: j for j, v in enumerate(xs)}
    y_index = {v: i for i, v in enumerate(ys)}

    grid = [[math.nan for _ in xs] for _ in ys]
    for r in rows:
        xv = _parse_float(r.get(x, "nan"))
        yv = _parse_float(r.get(y, "nan"))
        if math.isnan(xv) or math.isnan(yv):
            continue
        i = y_index.get(yv)
        j = x_index.get(xv)
        if i is None or j is None:
            continue
        grid[i][j] = _parse_float(r.get(z, "nan"))

    return xs, ys, grid


def plot_heatmaps(
    rows: List[Dict[str, str]],
    *,
    x: str,
    y: str,
    fields: List[str],
    out_dir: Path,
    title_prefix: str | None,
    meta: Dict[str, str] | None,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    for field in fields:
        xs, ys, grid = _build_grid(rows, x=x, y=y, z=field)
        if not xs or not ys:
            continue

        plt.figure(figsize=(7.6, 5.2))
        x0, x1 = min(xs), max(xs)
        y0, y1 = min(ys), max(ys)
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
        if meta:
            plt.colorbar(im, label=_axis_label(meta, axis="y", col=field))
            plt.xlabel(_axis_label(meta, axis="x", col=x))
            plt.ylabel(_axis_label(meta, axis="y", col=y))
        else:
            plt.colorbar(im, label=field)
            plt.xlabel(x)
            plt.ylabel(y)
        title = f"{title_prefix} - {field}" if title_prefix else field
        plt.title(title)
        plt.tight_layout()

        out = out_dir / f"heatmap_{field}_({y}_vs_{x}).png"
        plt.savefig(out, dpi=200)
        plt.close()
        print(f"Saved {out}")


def _resolve_path(p: Path) -> Path:
    """Resolve path relative to project root if not absolute."""
    if p.is_absolute():
        return p
    return PROJECT_ROOT / p


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Plot scan CSV (lines/heatmap)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"Project root: {PROJECT_ROOT}\nDefault output: {DEFAULT_OUT_DIR}",
    )
    ap.add_argument("--csv", type=Path, required=True, help="Input scan CSV (relative to project root or absolute)")
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR, help="Output directory")
    ap.add_argument("--mode", type=str, choices=["lines", "heatmap"], required=True)
    ap.add_argument("--title", type=str, default=None, help="Optional title prefix")
    ap.add_argument("--where", action="append", default=[], help="Filter clause col=value (repeatable)")
    ap.add_argument(
        "--split",
        type=str,
        default=None,
        help="Optional split-by column (e.g. process). Produces one set of figures per distinct value.",
    )

    # lines
    ap.add_argument("--x", type=str, default=None, help="x column")
    ap.add_argument("--ys", type=str, default=None, help="Comma-separated y columns (lines mode)")
    ap.add_argument("--group", type=str, default=None, help="Group-by column (lines mode)")

    # heatmap
    ap.add_argument("--y", type=str, default=None, help="y column (heatmap mode)")
    ap.add_argument("--fields", type=str, default=None, help="Comma-separated value fields (heatmap mode)")

    args = ap.parse_args()

    # Resolve paths relative to project root
    csv_path = _resolve_path(args.csv)
    out_dir = _resolve_path(args.out_dir)

    meta, rows = read_scan_csv(csv_path)
    rows = _apply_where(rows, args.where)
    if not rows:
        raise RuntimeError(f"No data rows after filtering: {csv_path}")

    title = args.title
    if title is None:
        # If producer wrote a human-friendly title, prefer it.
        title = meta.get("title")

    splits = _split_rows(rows, split=args.split)

    for split_value, subrows in splits:
        if not subrows:
            continue

        current_out_dir = out_dir
        title2 = title
        if args.split:
            safe = _sanitize_filename(split_value)
            current_out_dir = out_dir / f"{args.split}={safe}"
            title2 = f"{title} ({args.split}={split_value})" if title else f"{args.split}={split_value}"

        if args.mode == "lines":
            if not args.x or not args.ys:
                raise ValueError("lines mode requires --x and --ys")
            ys = [s.strip() for s in args.ys.split(",") if s.strip()]
            plot_lines(subrows, x=args.x, ys=ys, group=args.group, out_dir=current_out_dir, title_prefix=title2, meta=meta)
        else:
            if not args.x or not args.y or not args.fields:
                raise ValueError("heatmap mode requires --x, --y and --fields")
            fields = [s.strip() for s in args.fields.split(",") if s.strip()]
            plot_heatmaps(subrows, x=args.x, y=args.y, fields=fields, out_dir=current_out_dir, title_prefix=title2, meta=meta)


if __name__ == "__main__":
    main()
