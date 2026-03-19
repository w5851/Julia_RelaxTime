#!/usr/bin/env python3
"""通用扫描 CSV 作图脚本（支持可选的元数据头）。

CSV 规范（scan_csv_v1）
- 以 '#' 开头的行视为元数据/注释，会被 CSV 解析器忽略。
- 第一个非注释行是 CSV 表头。

支持两种模式
- lines：绘制 y(x) 曲线，可按某一列分组（例如同一张图多条线）。
- heatmap：绘制 z(x,y) 的热力图，可用 col=value 过滤数据。

常用可选参数（新）
- --xlim xmin xmax / --ylim ymin ymax：设置坐标范围（lines 会同时裁剪超出范围的数据点）。
- --xscale / --yscale：强制指定坐标轴缩放（覆盖 CSV 元数据里的 x_scale/y_scale）。
- heatmap：--zscale / --clim / --cmap：设置颜色缩放/范围/色图。

使用示例
--------
1) 弛豫时间 vs T（按 muB 分组）
    python scripts/plot_scan_csv.py \
        --mode lines \
        --csv data/outputs/results/relaxtime/scan/relaxation_times_vs_T.csv \
        --x T_MeV --ys tau_u,tau_s,tau_ubar,tau_sbar \
        --group muB_MeV \
        --out-dir data/outputs/figures/relaxtime

2) 固定 xi 的 gap/transport 热力图
    python scripts/plot_scan_csv.py \
        --mode heatmap \
        --csv data/outputs/results/relaxtime/scan/gap_transport_scan.csv \
        --x muq_MeV --y T_MeV --fields eta,sigma,tau_u \
        --where xi=0.0 \
        --out-dir data/outputs/figures/relaxtime

3) 固定 muB（或 muq），同一张图画不同 xi 的多条线：
     以 T 为横坐标，纵坐标为 eta_over_s 和 zeta_over_s（分别输出两张图），
     并限制横坐标 100..400、纵坐标对数轴且范围 1e-3..1e2。
    python scripts/plot_scan_csv.py \
        --mode lines \
        --csv data/outputs/results/relaxtime/scan/gap_transport_scan_xi-0p6to0p6.csv \
        --where muB_MeV=800.0 \
        --x T_MeV --ys eta_over_s,zeta_over_s \
        --group xi \
        --xlim 100 400 \
        --yscale log --ylim 1e-3 1e2 \
        --out-dir data/outputs/figures/relaxtime/gap_transport_by_xi_muB800
"""

from __future__ import annotations

import argparse
import csv
import io
import math
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib
from matplotlib.ticker import LogLocator, MaxNLocator
from cycler import cycler


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


def _parse_scale_arg(s: str | None) -> str | None:
    if s is None:
        return None
    s2 = str(s).strip().lower()
    if s2 in {"linear", "log"}:
        return s2
    return None


def _apply_xrange_filter(pairs: List[Tuple[float, float]], xlim: Tuple[float, float] | None) -> List[Tuple[float, float]]:
    if not xlim:
        return pairs
    xmin, xmax = xlim
    lo, hi = (xmin, xmax) if xmin <= xmax else (xmax, xmin)
    return [(x, y) for x, y in pairs if lo <= x <= hi]


def _apply_positive_filter_for_log(pairs: List[Tuple[float, float]], *, yscale: str | None) -> List[Tuple[float, float]]:
    if yscale != "log":
        return pairs
    # 对数轴要求 y>0；将非正值丢弃，避免报错/空图。
    return [(x, y) for x, y in pairs if y > 0]


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
            actual = str(r.get(k, ""))
            actual_num = _parse_float(actual)
            expected_num = _parse_float(v)
            if not math.isnan(actual_num) and not math.isnan(expected_num):
                if not math.isclose(actual_num, expected_num, rel_tol=1e-9, abs_tol=1e-12):
                    ok = False
                    break
            elif actual != v:
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
    multi_y: bool = False,
    out_dir: Path,
    title_prefix: str | None,
    meta: Dict[str, str] | None,
    xlim: Tuple[float, float] | None,
    ylim: Tuple[float, float] | None,
    xscale_override: str | None,
    yscale_override: str | None,
    marker: str | None,
    linewidth: float,
    grid_alpha: float,
    legend_loc: str | None,
    line_style: str = "-",
    figsize: Tuple[float, float] | None = None,
    show_title: bool = False,
    formats: List[str] = None,
    dpi: int = 600,
    check: bool = False,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    _figsize = figsize or APS_FIGSIZE["default"]

    if multi_y and group:
        raise ValueError("--multi-y cannot be combined with --group")

    if multi_y:
        fig, ax = plt.subplots(figsize=_figsize)

        xscale = xscale_override or (_axis_scale(meta, axis="x", col=x) if meta else None)
        yscale = yscale_override or (_axis_scale(meta, axis="y", col=ys[0]) if meta and ys else None)

        sub2 = sorted(rows, key=lambda r: _parse_float(r.get(x, "nan")))
        xs_all = [_parse_float(r.get(x, "nan")) for r in sub2]

        any_plotted = False
        for y in ys:
            ys_all = [_parse_float(r.get(y, "nan")) for r in sub2]
            pairs = [(xx, yy) for xx, yy in zip(xs_all, ys_all) if not (math.isnan(xx) or math.isnan(yy))]
            pairs = _apply_xrange_filter(pairs, xlim)
            pairs = _apply_positive_filter_for_log(pairs, yscale=yscale)
            if not pairs:
                continue
            xs3, ys3 = zip(*pairs)
            plot_kwargs = {"lw": linewidth, "linestyle": line_style, "label": y}
            if marker:
                plot_kwargs["marker"] = marker
            ax.plot(xs3, ys3, **plot_kwargs)
            any_plotted = True

        if not any_plotted:
            plt.close(fig)
            return

        if meta:
            ax.set_xlabel(_axis_label(meta, axis="x", col=x))
        else:
            ax.set_xlabel(x)

        if all(str(y).startswith("tau_") for y in ys):
            ax.set_ylabel("τ")
        else:
            ax.set_ylabel(_axis_label(meta, axis="y", col=ys[0]) if meta and ys else "Value")

        if xscale:
            ax.set_xscale(xscale)
        if yscale:
            ax.set_yscale(yscale)

        x_values = []
        y_values = []
        for line in ax.get_lines():
            xd = line.get_xdata()
            yd = line.get_ydata()
            x_values.extend([float(v) for v in xd if not math.isnan(v)])
            y_values.extend([float(v) for v in yd if not math.isnan(v)])

        _set_axis_limits_strict(ax, axis="x", values=x_values, user_lim=tuple(xlim) if xlim else None, scale=xscale)
        _set_axis_limits_strict(ax, axis="y", values=y_values, user_lim=tuple(ylim) if ylim else None, scale=yscale)
        _apply_axis_alignment(ax)

        if show_title and title_prefix:
            ax.set_title(title_prefix)
        if grid_alpha > 0:
            ax.grid(True, alpha=grid_alpha)
        if legend_loc:
            plt.legend(loc=legend_loc, frameon=False)
        else:
            plt.legend(frameon=False)
        fmts = formats or ["pdf", "png"]
        y_tag = "_".join([str(s) for s in ys])
        # keep filename reasonably short but stable
        if len(y_tag) > 80:
            y_tag = f"{ys[0]}_to_{ys[-1]}_{len(ys)}ys"
        y_tag = _sanitize_filename(y_tag)
        saved = []
        for fmt in fmts:
            out = out_dir / f"multi_y_{y_tag}_vs_{x}.{fmt}"
            fmt_lower = fmt.lower()
            fig.savefig(out, format=fmt_lower, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
            saved.append(out)
        plt.close(fig)
        print(f"Saved {out_dir} ({', '.join(fmts)})")
        if check:
            _post_save_checks(saved, expected_dpi=dpi)
        return

    def group_key(r: Dict[str, str]) -> str:
        return "__all__" if not group else str(r.get(group, ""))

    groups: Dict[str, List[Dict[str, str]]] = {}
    for r in rows:
        groups.setdefault(group_key(r), []).append(r)

    for y in ys:
        fig, ax = plt.subplots(figsize=_figsize)
        any_plotted = False

        # 轴缩放：优先使用命令行覆盖，其次使用元数据。
        xscale = xscale_override or (_axis_scale(meta, axis="x", col=x) if meta else None)
        yscale = yscale_override or (_axis_scale(meta, axis="y", col=y) if meta else None)

        for gk, sub in sorted(groups.items(), key=lambda kv: kv[0]):
            sub2 = sorted(sub, key=lambda r: _parse_float(r.get(x, "nan")))
            xs = [_parse_float(r.get(x, "nan")) for r in sub2]
            ys2 = [_parse_float(r.get(y, "nan")) for r in sub2]
            pairs = [(xx, yy) for xx, yy in zip(xs, ys2) if not (math.isnan(xx) or math.isnan(yy))]
            pairs = _apply_xrange_filter(pairs, xlim)
            pairs = _apply_positive_filter_for_log(pairs, yscale=yscale)
            if not pairs:
                continue
            xs3, ys3 = zip(*pairs)
            label = gk if group else y
            plot_kwargs = {"lw": linewidth, "linestyle": line_style, "label": label}
            if marker:
                plot_kwargs["marker"] = marker
            ax.plot(xs3, ys3, **plot_kwargs)
            any_plotted = True

        if not any_plotted:
            plt.close()
            continue

        if meta:
            ax.set_xlabel(_axis_label(meta, axis="x", col=x))
            ax.set_ylabel(_axis_label(meta, axis="y", col=y))
        else:
            ax.set_xlabel(x)
            ax.set_ylabel(y)

        if xscale:
            ax.set_xscale(xscale)
        if yscale:
            ax.set_yscale(yscale)

        # Apply strict axis-limit rules: user-specified limits win; otherwise compute from data and align to major ticks
        # Collect all plotted data on axes for limit calculation
        # For x: use all xs across groups; for y: compute from plotted y-values
        x_values = []
        y_values = []
        for line in ax.get_lines():
            xd = line.get_xdata()
            yd = line.get_ydata()
            x_values.extend([float(v) for v in xd if not math.isnan(v)])
            y_values.extend([float(v) for v in yd if not math.isnan(v)])

        _set_axis_limits_strict(ax, axis="x", values=x_values, user_lim=tuple(xlim) if xlim else None, scale=xscale)
        _set_axis_limits_strict(ax, axis="y", values=y_values, user_lim=tuple(ylim) if ylim else None, scale=yscale)

        _apply_axis_alignment(ax)

        if show_title:
            title = f"{title_prefix} - {y}" if title_prefix else y
            ax.set_title(title)
        if grid_alpha > 0:
            ax.grid(True, alpha=grid_alpha)
        if group:
            if legend_loc:
                plt.legend(loc=legend_loc, frameon=False)
            else:
                plt.legend(frameon=False)
        # 保存为指定格式（优先支持矢量 PDF）
        fmts = formats or ["pdf", "png"]
        saved = []
        for fmt in fmts:
            out = out_dir / f"{y}_vs_{x}.{fmt}"
            fmt_lower = fmt.lower()
            fig.savefig(out, format=fmt_lower, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
            saved.append(out)
        plt.close(fig)
        print(f"Saved {out_dir} ({', '.join(fmts)})")
        if check:
            _post_save_checks(saved, expected_dpi=dpi)


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


def _pick_bounds_from_ticks(ticks: List[float], lo: float, hi: float) -> Tuple[float, float]:
    t = [float(x) for x in ticks if not math.isnan(float(x))]
    if len(t) < 2:
        return lo, hi
    t = sorted(set(t))
    low = next((x for x in t if x <= lo), t[0])
    high = next((x for x in reversed(t) if x >= hi), t[-1])
    if low == high:
        return lo, hi
    return low, high


def _set_axis_limits_strict(ax, *, axis: str, values: List[float], user_lim: Tuple[float, float] | None, scale: str | None):
    """Set axis limits according to rule:
    - If user_lim provided, use it.
    - Else compute major ticks from data range and set axis limits to first/last major tick.
    - For log scale use LogLocator; for linear use AutoLocator/MaxNLocator.
    """
    if user_lim:
        if axis == "x":
            ax.set_xlim(user_lim)
        else:
            ax.set_ylim(user_lim)
        return

    vals = [v for v in values if not math.isnan(v)]
    if not vals:
        return

    try:
        if scale == "log":
            pos_vals = [v for v in vals if v > 0]
            if not pos_vals:
                return
            vmin = min(pos_vals)
            vmax = max(pos_vals)
            if math.isclose(vmin, vmax):
                low = 10 ** math.floor(math.log10(vmin))
                high = 10 ** math.ceil(math.log10(vmax * 1.01))
            else:
                locator = LogLocator(base=10.0)
                ticks = locator.tick_values(vmin, vmax)
                ticks = [tick for tick in ticks if math.isfinite(tick) and tick > 0]
                if ticks:
                    low, high = _pick_bounds_from_ticks(list(ticks), vmin, vmax)
                else:
                    low = 10 ** math.floor(math.log10(vmin))
                    high = 10 ** math.ceil(math.log10(vmax))
        else:
            vmin = min(vals)
            vmax = max(vals)
            if vmin == vmax:
                vmin -= 0.5
                vmax += 0.5
            # prefer MaxNLocator for nice round ticks
            locator = MaxNLocator(nbins=6)
            ticks = locator.tick_values(vmin, vmax)
            low, high = _pick_bounds_from_ticks(list(ticks), vmin, vmax)
    except Exception:
        vmin = min(vals)
        vmax = max(vals)
        low, high = vmin, vmax

    if axis == "x":
        ax.set_xlim(low, high)
    else:
        ax.set_ylim(low, high)


def _apply_axis_alignment(ax) -> None:
    """Apply stable axis styling without truncating visible spines.

    Earlier versions shortened left/bottom spines to the first and last major
    ticks. Combined with tight bounding-box export this could make the axis look
    visually broken near the plot edges. Keep zero margins, but leave the spines
    at their full extent so exported figures retain continuous axes.
    """
    try:
        ax.margins(x=0.0, y=0.0)
    except Exception:
        pass


def plot_heatmaps(
    rows: List[Dict[str, str]],
    *,
    x: str,
    y: str,
    fields: List[str],
    out_dir: Path,
    title_prefix: str | None,
    meta: Dict[str, str] | None,
    xlim: Tuple[float, float] | None,
    ylim: Tuple[float, float] | None,
    xscale_override: str | None,
    yscale_override: str | None,
    zscale: str | None,
    clim: Tuple[float, float] | None,
    cmap: str | None,
    figsize: Tuple[float, float] | None = None,
    show_title: bool = False,
    formats: List[str] = None,
    dpi: int = 600,
    check: bool = False,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    _figsize_hm = figsize or APS_FIGSIZE["double-column"]

    for field in fields:
        xs, ys, grid = _build_grid(rows, x=x, y=y, z=field)
        if not xs or not ys:
            continue

        # 可选裁剪（按 x/y 轴范围截取网格）
        if xlim:
            xmin, xmax = xlim
            lo, hi = (xmin, xmax) if xmin <= xmax else (xmax, xmin)
            x_keep = [j for j, xv in enumerate(xs) if lo <= xv <= hi]
        else:
            x_keep = list(range(len(xs)))
        if ylim:
            ymin, ymax = ylim
            lo, hi = (ymin, ymax) if ymin <= ymax else (ymax, ymin)
            y_keep = [i for i, yv in enumerate(ys) if lo <= yv <= hi]
        else:
            y_keep = list(range(len(ys)))

        xs2 = [xs[j] for j in x_keep]
        ys2 = [ys[i] for i in y_keep]
        grid2 = [[grid[i][j] for j in x_keep] for i in y_keep]
        if not xs2 or not ys2:
            continue

        fig, ax = plt.subplots(figsize=_figsize_hm)

        # 轴缩放：优先使用命令行覆盖，其次使用元数据。
        xscale = xscale_override or (_axis_scale(meta, axis="x", col=x) if meta else None)
        yscale = yscale_override or (_axis_scale(meta, axis="y", col=y) if meta else None)

        x0, x1 = min(xs2), max(xs2)
        y0, y1 = min(ys2), max(ys2)
        if x0 == x1:
            x0 -= 0.5
            x1 += 0.5
        if y0 == y1:
            y0 -= 0.5
            y1 += 0.5

        norm = None
        if zscale == "log":
            # LogNorm 要求 vmin>0；若未提供 clim，则从数据里找一个正的最小值。
            vmin = None
            vmax = None
            if clim:
                vmin, vmax = clim
            else:
                positive = [v for row in grid2 for v in row if (not math.isnan(v)) and v > 0]
                if positive:
                    vmin = min(positive)
                    vmax = max(positive)
            if vmin is not None and vmax is not None and vmin > 0 and vmax > 0:
                norm = LogNorm(vmin=vmin, vmax=vmax)

        im = ax.imshow(
            grid2,
            origin="lower",
            aspect="auto",
            extent=[x0, x1, y0, y1],
            interpolation="nearest",
            cmap=cmap,
            norm=norm,
        )

        # 线性颜色范围（优先级低于 zscale=log 的 norm）
        if norm is None and clim is not None:
            im.set_clim(clim[0], clim[1])

        if meta:
            fig.colorbar(im, ax=ax, label=_axis_label(meta, axis="y", col=field))
            ax.set_xlabel(_axis_label(meta, axis="x", col=x))
            ax.set_ylabel(_axis_label(meta, axis="y", col=y))
        else:
            fig.colorbar(im, ax=ax, label=field)
            ax.set_xlabel(x)
            ax.set_ylabel(y)

        if xscale:
            ax.set_xscale(xscale)
        if yscale:
            ax.set_yscale(yscale)

        # Apply strict axis limits based on grid coordinates xs2/ys2 unless user provided limits
        _set_axis_limits_strict(ax, axis="x", values=xs2, user_lim=tuple(xlim) if xlim else None, scale=xscale)
        _set_axis_limits_strict(ax, axis="y", values=ys2, user_lim=tuple(ylim) if ylim else None, scale=yscale)

        _apply_axis_alignment(ax)
        if show_title:
            title = f"{title_prefix} - {field}" if title_prefix else field
            ax.set_title(title)

        fmts = formats or ["pdf", "png"]
        saved = []
        for fmt in fmts:
            out = out_dir / f"heatmap_{field}_({y}_vs_{x}).{fmt}"
            fmt_lower = fmt.lower()
            fig.savefig(out, format=fmt_lower, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
            saved.append(out)
        plt.close(fig)
        print(f"Saved {out_dir} ({', '.join(fmts)})")
        if check:
            _post_save_checks(saved, expected_dpi=dpi)


def _resolve_path(p: Path) -> Path:
    """Resolve path relative to project root if not absolute."""
    if p.is_absolute():
        return p
    return PROJECT_ROOT / p


# APS standard figure widths (inches)
APS_FIGSIZE = {
    "single-column": (3.375, 2.5),
    "double-column": (6.75, 4.6),
    "default": (6.8, 4.6),
}

APS_COLOR_CYCLE = [
    "#4477AA",
    "#EE6677",
    "#228833",
    "#CCBB44",
    "#66CCEE",
    "#AA3377",
    "#BBBBBB",
]


def configure_publication_style(*, font: str = "Times New Roman", font_size: int = 10, line_width: float = 1.0, dpi: int = 600) -> None:
    """Configure matplotlib rcParams for publication-quality figures.

    Key APS (Physical Review) compliance points:
    - pdf.fonttype = 42 → TrueType font embedding (required by APS)
    - ps.fonttype = 42  → same for EPS output
    - axes.linewidth ≥ 0.5 pt
    - tick direction inward
    - font size 7–10 pt
    """
    matplotlib.rcParams.update({
        "font.family": "serif",
        "font.serif": [font, "Times"],
        "font.size": font_size,
        "mathtext.fontset": "stix",
        "axes.prop_cycle": cycler(color=APS_COLOR_CYCLE),
        "axes.titlesize": font_size,
        "axes.labelsize": font_size,
        "axes.linewidth": 0.6,
        "legend.fontsize": max(8, font_size - 2),
        "legend.frameon": False,
        "xtick.labelsize": max(8, font_size - 2),
        "ytick.labelsize": max(8, font_size - 2),
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "lines.linewidth": line_width,
        "lines.markersize": 4,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "axes.grid": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.dpi": dpi,
        "savefig.format": "pdf",
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.08,
    })


def _post_save_checks(paths: List[Path], *, expected_dpi: int = 600) -> None:
    """Run simple post-save checks on saved image files.

    - For PNG: check `info['dpi']` if available and warn if below expected_dpi.
    - For PDF/EPS: print a reminder to verify vectorness and font embedding.
    """
    try:
        from PIL import Image
    except Exception:
        Image = None

    for p in paths:
        try:
            p = Path(p)
            if not p.exists():
                print(f"[check] Missing file: {p}")
                continue
            suffix = p.suffix.lower()
            if suffix in {".png", ".tiff", ".tif"}:
                if Image is None:
                    print(f"[check] PIL not available — cannot inspect {p.name}")
                    continue
                with Image.open(p) as im:
                    dpi = im.info.get("dpi")
                    if dpi:
                        # dpi may be tuple
                        if isinstance(dpi, tuple):
                            dpi_val = int(min(dpi))
                        else:
                            dpi_val = int(dpi)
                        if dpi_val < expected_dpi:
                            print(f"[check] WARNING: {p.name} DPI={dpi_val} < expected {expected_dpi}")
                        else:
                            print(f"[check] OK: {p.name} DPI={dpi_val}")
                    else:
                        print(f"[check] NOTE: {p.name} has no DPI metadata; visually verify size/resolution")
            elif suffix in {".pdf", ".eps"}:
                print(f"[check] {p.name} created (PDF/EPS) — please verify vector content and font embedding if required.")
            else:
                print(f"[check] {p.name} created (unknown ext).")
        except Exception as e:
            print(f"[check] Error inspecting {p}: {e}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Plot scan CSV (lines/heatmap)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"Project root: {PROJECT_ROOT}\nDefault output: {DEFAULT_OUT_DIR}",
    )
    ap.add_argument("--csv", type=Path, required=True, help="Input scan CSV (relative to project root or absolute)")
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR, help="Output directory")
    ap.add_argument("--mode", type=str, choices=["lines", "heatmap"], required=True)
    ap.add_argument("--title", type=str, default=None, help="Optional title prefix (hidden by default; use --show-title to render)")
    ap.add_argument("--show-title", action="store_true", help="Render in-figure titles; disabled by default for publication-style output")
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
    ap.add_argument("--multi-y", action="store_true", help="Plot all --ys on a single figure (lines mode; cannot use --group)")

    # axis / style
    ap.add_argument("--xlim", type=float, nargs=2, default=None, metavar=("XMIN", "XMAX"), help="x axis limits")
    ap.add_argument("--ylim", type=float, nargs=2, default=None, metavar=("YMIN", "YMAX"), help="y axis limits")
    ap.add_argument("--xscale", type=str, default=None, choices=["linear", "log"], help="Override x axis scale")
    ap.add_argument("--yscale", type=str, default=None, choices=["linear", "log"], help="Override y axis scale")
    ap.add_argument("--marker", type=str, default=None, help="Marker style for lines (omit to disable)")
    ap.add_argument("--linewidth", type=float, default=1.5, help="Line width")
    ap.add_argument("--grid-alpha", type=float, default=0.0, help="Grid alpha")
    ap.add_argument("--legend-loc", type=str, default=None, help="Legend location (matplotlib loc string)")
    ap.add_argument("--line-style", type=str, default="-", help="Line style for lines (e.g. '-', '--', '-.', ':')")
    ap.add_argument("--check", action="store_true", help="Run simple post-save checks (PNG DPI, basic validity)")

    # heatmap
    ap.add_argument("--y", type=str, default=None, help="y column (heatmap mode)")
    ap.add_argument("--fields", type=str, default=None, help="Comma-separated value fields (heatmap mode)")
    ap.add_argument("--zscale", type=str, default=None, choices=["linear", "log"], help="Heatmap color scale")
    ap.add_argument("--clim", type=float, nargs=2, default=None, metavar=("VMIN", "VMAX"), help="Heatmap color limits")
    ap.add_argument("--cmap", type=str, default=None, help="Heatmap colormap name")

    # figure size presets (APS journal compliance)
    ap.add_argument(
        "--figsize",
        type=str,
        default=None,
        help=(
            "Figure size preset or WxH in inches. "
            "Presets: 'single-column' (3.375x2.5), 'double-column' (6.75x4.6). "
            "Custom example: '5.0x3.5'"
        ),
    )

    # output formats and quality
    ap.add_argument("--formats", type=str, default="pdf,png", help="Comma-separated output formats (e.g. pdf,png,eps)")
    ap.add_argument("--dpi", type=int, default=600, help="Output DPI for saved figures (default 600)")

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

    # Parse figsize
    figsize: Tuple[float, float] | None = None
    if args.figsize:
        fs = args.figsize.strip().lower()
        if fs in APS_FIGSIZE:
            figsize = APS_FIGSIZE[fs]
        elif "x" in fs:
            parts = fs.split("x", 1)
            figsize = (float(parts[0]), float(parts[1]))
        else:
            raise ValueError(f"Invalid --figsize: {args.figsize}. Use preset name or WxH (e.g. '5.0x3.5')")

    # Configure publication style before plotting
    configure_publication_style(dpi=args.dpi)

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
            plot_lines(
                subrows,
                x=args.x,
                ys=ys,
                group=args.group,
                multi_y=bool(args.multi_y),
                out_dir=current_out_dir,
                title_prefix=title2,
                meta=meta,
                xlim=tuple(args.xlim) if args.xlim else None,
                ylim=tuple(args.ylim) if args.ylim else None,
                xscale_override=_parse_scale_arg(args.xscale),
                yscale_override=_parse_scale_arg(args.yscale),
                marker=args.marker,
                linewidth=float(args.linewidth),
                grid_alpha=float(args.grid_alpha),
                legend_loc=args.legend_loc,
                line_style=args.line_style,
                figsize=figsize,
                show_title=bool(args.show_title),
                formats=[s.strip() for s in args.formats.split(",") if s.strip()],
                dpi=int(args.dpi),
                check=bool(args.check),
            )
        else:
            if not args.x or not args.y or not args.fields:
                raise ValueError("heatmap mode requires --x, --y and --fields")
            fields = [s.strip() for s in args.fields.split(",") if s.strip()]
            plot_heatmaps(
                subrows,
                x=args.x,
                y=args.y,
                fields=fields,
                out_dir=current_out_dir,
                title_prefix=title2,
                meta=meta,
                xlim=tuple(args.xlim) if args.xlim else None,
                ylim=tuple(args.ylim) if args.ylim else None,
                xscale_override=_parse_scale_arg(args.xscale),
                yscale_override=_parse_scale_arg(args.yscale),
                zscale=_parse_scale_arg(args.zscale),
                clim=tuple(args.clim) if args.clim else None,
                cmap=args.cmap,
                figsize=figsize,
                show_title=bool(args.show_title),
                formats=[s.strip() for s in args.formats.split(",") if s.strip()],
                dpi=int(args.dpi),
                check=bool(args.check),
            )


if __name__ == "__main__":
    main()
