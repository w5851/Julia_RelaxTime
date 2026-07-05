#!/usr/bin/env python3
"""Build formal Figure 4 PNJL phase-diagram assets from a dense reference."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


def find_project_root() -> Path:
    script_dir = Path(__file__).resolve().parent
    candidates = [script_dir, script_dir.parent, script_dir.parent.parent, Path.cwd()]
    for start in candidates:
        current = start
        for _ in range(8):
            if (current / "Project.toml").is_file() or (current / ".git").exists():
                return current
            parent = current.parent
            if parent == current:
                break
            current = parent
    return Path.cwd()


PROJECT_ROOT = find_project_root()
CASE_SLUG = "figure4_phase_diagram_prod_v1"
DEFAULT_RESULT_DIR = PROJECT_ROOT / "data" / "outputs" / "results" / "pnjl" / "phase_diagram" / CASE_SLUG
DEFAULT_FIGURE_DIR = PROJECT_ROOT / "data" / "outputs" / "figures" / "pnjl" / "phase_diagram" / CASE_SLUG
DEFAULT_REFERENCE_ROOT = DEFAULT_RESULT_DIR / "reference"
DEFAULT_REFERENCE_TAG = "figure4_phase_diagram_prod_v1_c1_p24t8"
DEFAULT_XI_VALUES = "-0.5,-0.25,0.0,0.25,0.5"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build formal PNJL Figure 4 T-muB/T-rho assets from dense phase-reference CSVs."
    )
    parser.add_argument("--phase-reference-root", type=Path, default=DEFAULT_REFERENCE_ROOT)
    parser.add_argument("--phase-reference-tag", default=DEFAULT_REFERENCE_TAG)
    parser.add_argument("--result-dir", type=Path, default=DEFAULT_RESULT_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument(
        "--xi-values",
        default=DEFAULT_XI_VALUES,
        help="Comma-separated xi subset to plot, or 'all'. Full data are still written to figure_assets.",
    )
    parser.add_argument("--include-spinodal", action="store_true", help="Show spinodal curves in the plot.")
    parser.add_argument("--formats", default="png,pdf", help="Comma-separated output formats.")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--figsize", default="7.6,3.8", help="Matplotlib figure size in inches, width,height.")
    return parser.parse_args()


def fail(message: str) -> None:
    raise SystemExit(f"[figure4-phase-assets] {message}")


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        fail(f"missing CSV: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}


def finite_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        x = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def row_float(row: dict[str, Any], key: str) -> float | None:
    return finite_float(row.get(key))


def row_float_first(row: dict[str, Any], keys: Iterable[str]) -> float | None:
    for key in keys:
        value = row_float(row, key)
        if value is not None:
            return value
    return None


def first_item(item: tuple[Any, ...]) -> Any:
    return item[0]


def asset_sort_key(row: dict[str, Any]) -> tuple[str, float, float, str]:
    return (
        str(row["kind"]),
        float(row["xi"]),
        float(row.get("plot_order_key") or 0.0),
        str(row["branch"]),
    )


def plot_order_sort_key(row: dict[str, Any]) -> float:
    return row_float(row, "plot_order_key") or 0.0


def fmt(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        return f"{value:.12g}"
    return str(value)


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field)) for field in fields})


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=True, indent=2) + "\n", encoding="utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relpath(path: Path) -> str:
    try:
        return path.resolve().relative_to(PROJECT_ROOT.resolve()).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=str(PROJECT_ROOT),
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return "unknown"


def parse_formats(raw: str) -> list[str]:
    formats = [part.strip().lower() for part in raw.split(",") if part.strip()]
    if not formats:
        fail("--formats must contain at least one format")
    allowed = {"png", "pdf", "svg"}
    bad = [item for item in formats if item not in allowed]
    if bad:
        fail(f"unsupported formats: {', '.join(bad)}")
    return formats


def parse_figsize(raw: str) -> tuple[float, float]:
    parts = [part.strip() for part in raw.split(",")]
    if len(parts) != 2:
        fail("--figsize must be width,height")
    width = finite_float(parts[0])
    height = finite_float(parts[1])
    if width is None or height is None or width <= 0 or height <= 0:
        fail("--figsize must contain positive finite numbers")
    return width, height


def parse_xi_values(raw: str | None, source_xis: list[float]) -> list[float]:
    if raw is None or not raw.strip() or raw.strip().lower() == "all":
        return source_xis
    requested = []
    for token in raw.split(","):
        token = token.strip()
        if not token:
            continue
        value = finite_float(token)
        if value is None:
            fail(f"invalid xi value: {token}")
        requested.append(round(value, 10))
    if not requested:
        fail("--xi-values did not contain any usable values")
    source_set = {round(xi, 10) for xi in source_xis}
    missing = [xi for xi in requested if xi not in source_set]
    if missing:
        fail(f"requested xi values not present in source reference: {missing}")
    return requested


def reference_paths(root: Path, tag: str) -> dict[str, Path]:
    return {
        "boundary": root / f"boundary_{tag}.csv",
        "cep": root / f"cep_{tag}.csv",
        "crossover": root / f"crossover_{tag}.csv",
        "spinodals": root / f"spinodals_{tag}.csv",
    }


def collect_source_xis(*row_groups: list[dict[str, str]]) -> list[float]:
    values = set()
    for rows in row_groups:
        for row in rows:
            xi = row_float(row, "xi")
            if xi is not None:
                values.add(round(xi, 10))
    return sorted(values)


def plotted_flag(row: dict[str, Any], plotted_xis: set[float], include_spinodal: bool) -> bool:
    xi = row.get("xi")
    if xi is None or round(float(xi), 10) not in plotted_xis:
        return False
    if row.get("kind") == "spinodal":
        return include_spinodal
    return True


def estimate_cep_rho(boundary_rows: list[dict[str, str]], xi: float, t_cep: float) -> tuple[float | None, str]:
    candidates = []
    for row in boundary_rows:
        row_xi = row_float(row, "xi")
        t_value = row_float(row, "T_MeV")
        rho_h = row_float(row, "rho_hadron")
        rho_q = row_float(row, "rho_quark")
        if row_xi is None or t_value is None or rho_h is None or rho_q is None:
            continue
        if abs(row_xi - xi) > 1.0e-8:
            continue
        candidates.append((abs(t_value - t_cep), rho_h, rho_q, t_value))
    if not candidates:
        return None, "missing_boundary_match"
    _, rho_h, rho_q, t_match = min(candidates, key=first_item)
    return 0.5 * (rho_h + rho_q), f"mean_boundary_density_at_nearest_T={t_match:.12g}"


def build_asset_rows(
    paths: dict[str, Path],
    boundary_rows: list[dict[str, str]],
    cep_rows: list[dict[str, str]],
    crossover_rows: list[dict[str, str]],
    spinodal_rows: list[dict[str, str]],
    plotted_xis: set[float],
    include_spinodal: bool,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    tmu_rows: list[dict[str, Any]] = []
    trho_rows: list[dict[str, Any]] = []

    for idx, row in enumerate(boundary_rows, start=1):
        xi = row_float(row, "xi")
        t_mev = row_float(row, "T_MeV")
        muq = row_float(row, "mu_transition_MeV")
        rho_h = row_float(row, "rho_hadron")
        rho_q = row_float(row, "rho_quark")
        if xi is None or t_mev is None or muq is None:
            continue
        curve_parameter = row_float_first(row, ("curve_parameter", "plot_order_key", "T_MeV"))
        plot_order_key = row_float_first(row, ("plot_order_key", "curve_parameter", "T_MeV"))
        base = {
            "kind": "first_order",
            "xi": xi,
            "T_MeV": t_mev,
            "muq_MeV": muq,
            "muB_MeV": 3.0 * muq,
            "variable": "mu_transition_MeV",
            "curve_parameter": curve_parameter,
            "plot_order_key": plot_order_key,
            "source_csv": relpath(paths["boundary"]),
            "source_row": idx,
        }
        tmu_row = {**base, "branch": "coexistence", "rho_over_rho0": ""}
        tmu_row["plotted"] = plotted_flag(tmu_row, plotted_xis, include_spinodal)
        tmu_rows.append(tmu_row)
        for branch, rho in (("hadron", rho_h), ("quark", rho_q)):
            if rho is None:
                continue
            trho_row = {
                **base,
                "kind": "coexistence",
                "branch": branch,
                "rho_over_rho0": rho,
                "rho_method": "boundary_" + branch,
            }
            trho_row["plotted"] = plotted_flag(trho_row, plotted_xis, include_spinodal)
            trho_rows.append(trho_row)

    for idx, row in enumerate(crossover_rows, start=1):
        xi = row_float(row, "xi")
        muq = row_float(row, "mu_MeV")
        t_mev = row_float(row, "T_crossover_MeV")
        rho = row_float(row, "rho")
        if xi is None or muq is None or t_mev is None:
            continue
        variable = row.get("variable") or "phi_u"
        curve_parameter = row_float_first(row, ("curve_parameter", "plot_order_key", "mu_MeV"))
        plot_order_key = row_float_first(row, ("plot_order_key", "curve_parameter", "mu_MeV"))
        base = {
            "kind": "crossover",
            "branch": variable,
            "xi": xi,
            "T_MeV": t_mev,
            "muq_MeV": muq,
            "muB_MeV": 3.0 * muq,
            "rho_over_rho0": rho,
            "variable": variable,
            "curve_parameter": curve_parameter,
            "plot_order_key": plot_order_key,
            "source_csv": relpath(paths["crossover"]),
            "source_row": idx,
        }
        base["plotted"] = plotted_flag(base, plotted_xis, include_spinodal)
        tmu_rows.append(dict(base))
        trho_row = {**base, "rho_method": "crossover_reference_rho"}
        trho_rows.append(trho_row)

    for idx, row in enumerate(cep_rows, start=1):
        xi = row_float(row, "xi")
        t_cep = row_float(row, "T_CEP_MeV")
        mu_b = row_float(row, "muB_CEP_MeV")
        mu_q = row_float(row, "muq_CEP_MeV")
        if mu_b is None and mu_q is not None:
            mu_b = 3.0 * mu_q
        if mu_q is None and mu_b is not None:
            mu_q = mu_b / 3.0
        if xi is None or t_cep is None or mu_b is None:
            continue
        base = {
            "kind": "cep",
            "branch": "cep",
            "xi": xi,
            "T_MeV": t_cep,
            "muq_MeV": mu_q,
            "muB_MeV": mu_b,
            "variable": "CEP",
            "curve_parameter": t_cep,
            "plot_order_key": t_cep,
            "source_csv": relpath(paths["cep"]),
            "source_row": idx,
        }
        tmu_row = {**base, "rho_over_rho0": ""}
        tmu_row["plotted"] = plotted_flag(tmu_row, plotted_xis, include_spinodal)
        tmu_rows.append(tmu_row)
        rho_cep, rho_method = estimate_cep_rho(boundary_rows, xi, t_cep)
        trho_row = {**base, "rho_over_rho0": rho_cep, "rho_method": rho_method}
        trho_row["plotted"] = plotted_flag(trho_row, plotted_xis, include_spinodal)
        trho_rows.append(trho_row)
        if rho_cep is not None:
            connector_base = {
                **base,
                "rho_over_rho0": rho_cep,
                "variable": "CEP_visual_connector",
                "rho_method": "visual_connector_" + rho_method,
                "curve_parameter": t_cep,
                "source_csv": relpath(paths["cep"]),
                "source_row": idx,
            }
            for branch in ("hadron", "quark"):
                connector_row = {
                    **connector_base,
                    "kind": "coexistence",
                    "branch": branch,
                    "plot_order_key": t_cep,
                }
                connector_row["plotted"] = plotted_flag(connector_row, plotted_xis, include_spinodal)
                trho_rows.append(connector_row)
            crossover_connector = {
                **connector_base,
                "kind": "crossover",
                "branch": "zz_CEP_visual_connector",
                "plot_order_key": (mu_q if mu_q is not None else mu_b / 3.0) + 1.0e-9,
            }
            crossover_connector["plotted"] = plotted_flag(crossover_connector, plotted_xis, include_spinodal)
            trho_rows.append(crossover_connector)

    for idx, row in enumerate(spinodal_rows, start=1):
        xi = row_float(row, "xi")
        t_mev = row_float(row, "T_MeV")
        if xi is None or t_mev is None:
            continue
        curve_parameter = row_float_first(row, ("curve_parameter", "plot_order_key", "T_MeV"))
        plot_order_key = row_float_first(row, ("plot_order_key", "curve_parameter", "T_MeV"))
        for branch, mu_col, rho_col in (
            ("hadron", "mu_spinodal_hadron_MeV", "rho_spinodal_hadron"),
            ("quark", "mu_spinodal_quark_MeV", "rho_spinodal_quark"),
        ):
            muq = row_float(row, mu_col)
            rho = row_float(row, rho_col)
            if muq is None:
                continue
            base = {
                "kind": "spinodal",
                "branch": branch,
                "xi": xi,
                "T_MeV": t_mev,
                "muq_MeV": muq,
                "muB_MeV": 3.0 * muq,
                "rho_over_rho0": rho,
                "variable": mu_col,
                "curve_parameter": curve_parameter,
                "plot_order_key": plot_order_key,
                "source_csv": relpath(paths["spinodals"]),
                "source_row": idx,
            }
            base["plotted"] = plotted_flag(base, plotted_xis, include_spinodal)
            tmu_rows.append(dict(base))
            trho_row = {**base, "rho_method": "spinodal_" + branch}
            trho_rows.append(trho_row)

    tmu_rows.sort(key=asset_sort_key)
    trho_rows.sort(key=asset_sort_key)
    return tmu_rows, trho_rows


def grouped_rows(rows: list[dict[str, Any]], *keys: str) -> dict[Any, list[dict[str, Any]]]:
    grouped: dict[Any, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        key = row.get(keys[0]) if len(keys) == 1 else tuple(row.get(key) for key in keys)
        grouped[key].append(row)
    return grouped


def row_is_plotted(row: dict[str, Any]) -> bool:
    return bool(row.get("plotted"))


def color_map_for_xis(xis: list[float]) -> dict[float, str]:
    base_palette = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB"]
    if len(xis) <= len(base_palette):
        return {xi: base_palette[idx] for idx, xi in enumerate(xis)}
    import matplotlib.pyplot as plt

    cmap = plt.get_cmap("viridis")
    denom = max(len(xis) - 1, 1)
    return {xi: cmap(idx / denom) for idx, xi in enumerate(xis)}


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    values = []
    for row in rows:
        value = row_float(row, key)
        if value is not None:
            values.append(value)
    return values


def set_axis_limits(ax1: Any, ax2: Any, tmu_rows: list[dict[str, Any]], trho_rows: list[dict[str, Any]]) -> None:
    plotted_tmu = [row for row in tmu_rows if row_is_plotted(row)]
    plotted_trho = [row for row in trho_rows if row_is_plotted(row)]
    t_values = finite_values(plotted_tmu + plotted_trho, "T_MeV")
    mu_values = finite_values(plotted_tmu, "muB_MeV")
    rho_values = finite_values(plotted_trho, "rho_over_rho0")
    if t_values:
        t_min = max(0.0, min(t_values) - 8.0)
        t_max = max(t_values) + 8.0
        ax1.set_ylim(t_min, t_max)
        ax2.set_ylim(t_min, t_max)
    if mu_values:
        ax1.set_xlim(0.0, max(mu_values) * 1.04)
    if rho_values:
        ax2.set_xlim(0.0, max(rho_values) * 1.08)


def plot_assets(
    figure_dir: Path,
    formats: list[str],
    dpi: int,
    figsize: tuple[float, float],
    tmu_rows: list[dict[str, Any]],
    trho_rows: list[dict[str, Any]],
    plotted_xis: list[float],
    include_spinodal: bool,
) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    figure_dir.mkdir(parents=True, exist_ok=True)
    colors = color_map_for_xis(plotted_xis)

    fig, (ax_mu, ax_rho) = plt.subplots(1, 2, figsize=figsize, sharey=True)
    plotted_tmu = [row for row in tmu_rows if row_is_plotted(row)]
    plotted_trho = [row for row in trho_rows if row_is_plotted(row)]

    for xi, group in grouped_rows([row for row in plotted_tmu if row["kind"] == "first_order"], "xi").items():
        group = sorted(group, key=plot_order_sort_key)
        ax_mu.plot(
            finite_values(group, "muB_MeV"),
            finite_values(group, "T_MeV"),
            color=colors.get(round(float(xi), 10), "#333333"),
            linewidth=1.45,
            solid_capstyle="round",
        )

    for xi, group in grouped_rows([row for row in plotted_tmu if row["kind"] == "crossover"], "xi").items():
        group = sorted(group, key=plot_order_sort_key)
        ax_mu.plot(
            finite_values(group, "muB_MeV"),
            finite_values(group, "T_MeV"),
            color=colors.get(round(float(xi), 10), "#333333"),
            linestyle=(0, (4, 2)),
            linewidth=1.1,
            alpha=0.72,
        )

    if include_spinodal:
        for (xi, branch), group in grouped_rows([row for row in plotted_tmu if row["kind"] == "spinodal"], "xi", "branch").items():
            group = sorted(group, key=plot_order_sort_key)
            ax_mu.plot(
                finite_values(group, "muB_MeV"),
                finite_values(group, "T_MeV"),
                color=colors.get(round(float(xi), 10), "#333333"),
                linestyle=":",
                linewidth=0.85,
                alpha=0.42,
            )

    cep_rows = [row for row in plotted_tmu if row["kind"] == "cep"]
    if cep_rows:
        ax_mu.scatter(
            finite_values(cep_rows, "muB_MeV"),
            finite_values(cep_rows, "T_MeV"),
            c=[colors.get(round(float(row["xi"]), 10), "#333333") for row in cep_rows],
            marker="*",
            s=58,
            edgecolors="#111111",
            linewidths=0.35,
            zorder=5,
        )

    for xi in plotted_xis:
        xi_rows = [row for row in plotted_trho if row["kind"] == "coexistence" and abs(float(row["xi"]) - xi) < 1.0e-8]
        hadron = sorted([row for row in xi_rows if row["branch"] == "hadron"], key=plot_order_sort_key)
        quark = sorted([row for row in xi_rows if row["branch"] == "quark"], key=plot_order_sort_key)
        color = colors.get(round(xi, 10), "#333333")
        if hadron and quark and len(hadron) == len(quark):
            t_values = finite_values(hadron, "T_MeV")
            rho_h = finite_values(hadron, "rho_over_rho0")
            rho_q = finite_values(quark, "rho_over_rho0")
            ax_rho.fill_betweenx(t_values, rho_h, rho_q, color=color, alpha=0.10, linewidth=0)
            ax_rho.plot(rho_h, t_values, color=color, linewidth=1.15)
            ax_rho.plot(rho_q, t_values, color=color, linewidth=1.15)

    for xi, group in grouped_rows([row for row in plotted_trho if row["kind"] == "crossover"], "xi").items():
        group = sorted(group, key=plot_order_sort_key)
        ax_rho.plot(
            finite_values(group, "rho_over_rho0"),
            finite_values(group, "T_MeV"),
            color=colors.get(round(float(xi), 10), "#333333"),
            linestyle=(0, (4, 2)),
            linewidth=1.1,
            alpha=0.72,
        )

    if include_spinodal:
        for (xi, branch), group in grouped_rows([row for row in plotted_trho if row["kind"] == "spinodal"], "xi", "branch").items():
            group = sorted(group, key=plot_order_sort_key)
            ax_rho.plot(
                finite_values(group, "rho_over_rho0"),
                finite_values(group, "T_MeV"),
                color=colors.get(round(float(xi), 10), "#333333"),
                linestyle=":",
                linewidth=0.85,
                alpha=0.42,
            )

    cep_rho_rows = [row for row in plotted_trho if row["kind"] == "cep" and row_float(row, "rho_over_rho0") is not None]
    if cep_rho_rows:
        ax_rho.scatter(
            finite_values(cep_rho_rows, "rho_over_rho0"),
            finite_values(cep_rho_rows, "T_MeV"),
            c=[colors.get(round(float(row["xi"]), 10), "#333333") for row in cep_rho_rows],
            marker="*",
            s=58,
            edgecolors="#111111",
            linewidths=0.35,
            zorder=5,
        )

    ax_mu.set_title(r"$T$-$\mu_B$")
    ax_rho.set_title(r"$T$-$\rho$")
    ax_mu.set_xlabel(r"$\mu_B$ (MeV)")
    ax_rho.set_xlabel(r"$\rho/\rho_0$")
    ax_mu.set_ylabel(r"$T$ (MeV)")
    ax_mu.text(0.03, 0.96, "(a)", transform=ax_mu.transAxes, va="top", ha="left", fontweight="bold")
    ax_rho.text(0.03, 0.96, "(b)", transform=ax_rho.transAxes, va="top", ha="left", fontweight="bold")
    for ax in (ax_mu, ax_rho):
        ax.grid(True, color="#E2E2E2", linewidth=0.6)
        ax.minorticks_on()
        ax.tick_params(direction="in", top=True, right=True, which="both")
    set_axis_limits(ax_mu, ax_rho, tmu_rows, trho_rows)

    color_handles = [
        Line2D([], [], color=colors[xi], linewidth=1.8, label=rf"$\xi={xi:g}$")
        for xi in plotted_xis
    ]
    style_handles = [
        Line2D([], [], color="#333333", linewidth=1.45, linestyle="-", label="first order"),
        Line2D([], [], color="#333333", linewidth=1.1, linestyle=(0, (4, 2)), label="crossover"),
        Line2D([], [], color="#333333", marker="*", linestyle="None", markersize=7, label="CEP"),
    ]
    if include_spinodal:
        style_handles.append(Line2D([], [], color="#333333", linewidth=0.9, linestyle=":", label="spinodal"))
    fig.legend(
        handles=color_handles + style_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.02),
        ncol=min(6, len(color_handles) + len(style_handles)),
        frameon=False,
        fontsize=7.0,
        handlelength=2.0,
        columnspacing=1.2,
    )
    fig.tight_layout(rect=(0.0, 0.12, 1.0, 1.0), w_pad=1.2)

    written = []
    for fmt_name in formats:
        path = figure_dir / f"figure4_phase_diagram_TmuB_Trho.{fmt_name}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        written.append(path)
    plt.close(fig)
    return written


def summarize_counts(rows: list[dict[str, Any]]) -> dict[str, Any]:
    counts = Counter(str(row.get("kind", "")) for row in rows)
    plotted_counts = Counter(str(row.get("kind", "")) for row in rows if row_is_plotted(row))
    return {
        "rows_by_kind": dict(sorted(counts.items())),
        "plotted_rows_by_kind": dict(sorted(plotted_counts.items())),
    }


def main() -> int:
    args = parse_args()
    formats = parse_formats(args.formats)
    figsize = parse_figsize(args.figsize)
    paths = reference_paths(args.phase_reference_root, args.phase_reference_tag)
    boundary_rows = read_csv(paths["boundary"])
    cep_rows = read_csv(paths["cep"])
    crossover_rows = read_csv(paths["crossover"])
    spinodal_rows = read_csv(paths["spinodals"])
    source_xis = collect_source_xis(boundary_rows, cep_rows, crossover_rows, spinodal_rows)
    if not source_xis:
        fail("no xi values found in source reference")
    plotted_xis = parse_xi_values(args.xi_values, source_xis)
    plotted_xi_set = {round(xi, 10) for xi in plotted_xis}

    asset_dir = args.result_dir / "figure_assets"
    tmu_rows, trho_rows = build_asset_rows(
        paths,
        boundary_rows,
        cep_rows,
        crossover_rows,
        spinodal_rows,
        plotted_xi_set,
        args.include_spinodal,
    )

    tmu_asset = asset_dir / "figure4_phase_lines_TmuB.csv"
    trho_asset = asset_dir / "figure4_phase_lines_Trho.csv"
    summary_asset = asset_dir / "figure4_phase_plot_inputs_summary.json"
    tmu_fields = [
        "kind",
        "branch",
        "xi",
        "T_MeV",
        "muq_MeV",
        "muB_MeV",
        "rho_over_rho0",
        "variable",
        "curve_parameter",
        "plot_order_key",
        "plotted",
        "source_csv",
        "source_row",
    ]
    trho_fields = [
        "kind",
        "branch",
        "xi",
        "T_MeV",
        "rho_over_rho0",
        "muq_MeV",
        "muB_MeV",
        "variable",
        "rho_method",
        "curve_parameter",
        "plot_order_key",
        "plotted",
        "source_csv",
        "source_row",
    ]
    write_csv(tmu_asset, tmu_rows, tmu_fields)
    write_csv(trho_asset, trho_rows, trho_fields)

    source_manifest_path = args.result_dir / "phase_reference_source_manifest.json"
    source_manifest = read_json(source_manifest_path)
    summary_payload = {
        "schema_version": "figure4_phase_plot_inputs_summary_v1",
        "phase_reference_root": relpath(args.phase_reference_root),
        "phase_reference_tag": args.phase_reference_tag,
        "source_xi_values": source_xis,
        "plotted_xi_values": plotted_xis,
        "include_spinodal": args.include_spinodal,
        "axis_conversion": {
            "boundary": "muB_MeV = 3 * mu_transition_MeV",
            "crossover": "muB_MeV = 3 * mu_MeV",
            "spinodals": "muB_MeV = 3 * mu_spinodal_*_MeV",
            "cep": "muB_CEP_MeV is used directly when present; otherwise muB_MeV = 3 * muq_CEP_MeV",
        },
        "rho_units": "rho/rho0",
        "visual_connector_policy": {
            "Trho_CEP": (
                "Plotting-only connector rows are added to the T-rho figure asset so that "
                "the coexistence branches and crossover curve visually pass through the CEP. "
                "The formal source data remain unchanged."
            ),
            "CEP_rho_coordinate": "mean boundary density at the nearest available boundary temperature for the same xi",
        },
        "counts": {
            "source": {
                "boundary_rows": len(boundary_rows),
                "cep_rows": len(cep_rows),
                "crossover_rows": len(crossover_rows),
                "spinodal_rows": len(spinodal_rows),
            },
            "TmuB_assets": summarize_counts(tmu_rows),
            "Trho_assets": summarize_counts(trho_rows),
        },
        "source_manifest_verdict": source_manifest.get("verdict"),
        "source_manifest_residual_risks": source_manifest.get("residual_risks", []),
    }
    write_json(summary_asset, summary_payload)

    figure_paths = plot_assets(
        args.figure_dir,
        formats,
        args.dpi,
        figsize,
        tmu_rows,
        trho_rows,
        plotted_xis,
        args.include_spinodal,
    )

    asset_files = [tmu_asset, trho_asset, summary_asset]
    input_files = list(paths.values())
    if source_manifest_path.is_file():
        input_files.append(source_manifest_path)
    manifest = {
        "schema_version": "figure4_phase_diagram_plot_manifest_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "script": relpath(Path(__file__)),
        "command": " ".join([sys.executable, *sys.argv]),
        "git_commit": git_commit(),
        "phase_reference": {
            "root": relpath(args.phase_reference_root),
            "tag": args.phase_reference_tag,
            "source_manifest": relpath(source_manifest_path) if source_manifest_path.is_file() else None,
            "source_verdict": source_manifest.get("verdict"),
        },
        "source_xi_values": source_xis,
        "plotted_xi_values": plotted_xis,
        "include_spinodal": args.include_spinodal,
        "formats": formats,
        "dpi": args.dpi,
        "figsize_inches": list(figsize),
        "axis_conversion": summary_payload["axis_conversion"],
        "rho_units": summary_payload["rho_units"],
        "visual_connector_policy": summary_payload["visual_connector_policy"],
        "inputs": [
            {"path": relpath(path), "sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in input_files
            if path.is_file()
        ],
        "result_side_assets": [
            {"path": relpath(path), "sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in asset_files
        ],
        "figures": [
            {"path": relpath(path), "sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in figure_paths
        ],
        "counts": summary_payload["counts"],
        "style": {
            "first_order": "solid colored line",
            "crossover": "short dashed colored line",
            "cep": "star marker",
            "spinodal": "dotted low-alpha colored line when include_spinodal=true",
            "default_spinodal_layer": "off",
        },
        "residual_risks": source_manifest.get("residual_risks", []),
    }
    write_json(args.figure_dir / "plot_manifest.json", manifest)
    print(json.dumps({"figure_dir": relpath(args.figure_dir), "asset_dir": relpath(asset_dir), "figures": [relpath(path) for path in figure_paths]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
