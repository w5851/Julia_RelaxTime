#!/usr/bin/env python3
"""Build the semantic-correction publication-clean RS display layer.

This is a solver-free derivative of the author-accepted ``prod_v2`` result.
It intentionally leaves ``publication_clean_v1`` untouched and changes only
the display semantics around first-order intervals:

* no line segment bridges an audited first-order gap;
* historical midpoint stars are retained as provenance only and are never
  rendered as CEP markers;
* endpoint markers expose the observed quark/hadron (chiral-restored/
  chiral-broken) side when the raw scan provides that label;
* high-range first-order panels use an explicit log-y display transform and
  publication typography rounds only the rendered MeV labels;
* the large-eta/s left side of mode B remains a writing-layer caveat, not a
  numerical edit.

No raw result, registry entry, or solver output is modified.
"""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import importlib.util
import json
import math
import re
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[3]
V1_SCRIPT = ROOT / "scripts" / "analysis" / "relaxtime" / "build_phase_guided_publication_clean_v1.py"
SPEC = importlib.util.spec_from_file_location("phase_guided_publication_clean_v1", V1_SCRIPT)
if SPEC is None or SPEC.loader is None:  # pragma: no cover
    raise RuntimeError(f"unable to load v1 builder: {V1_SCRIPT}")
V1 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V1)


OUT_DIR = (
    ROOT
    / "docs"
    / "analysis"
    / "relaxtime"
    / "phase_guided_transport"
    / "phase_guided_transport_publication_clean_v2"
)
TABLE_DIR = OUT_DIR / "tables"
FIGURE_DIR = OUT_DIR / "figures"

DISPLAY_FIELDS = V1.DISPLAY_FIELDS
CORE_DISPLAY_FIELDS = V1.CORE_DISPLAY_FIELDS
OBSERVABLE_LABELS = V1.OBSERVABLE_LABELS
MODE_CONFIG = V1.MODE_CONFIG
CURRENT_CASE = V1.CURRENT_CASE
CALCULATION_SHA = V1.CALCULATION_SHA
WORKFLOW_HEAD_SHA = V1.WORKFLOW_HEAD_SHA
RECIPE_DIR = V1.RECIPE_DIR
REPLACEMENT_RECIPE = V1.REPLACEMENT_RECIPE
SMOOTHING_WINDOW_RECIPE = V1.SMOOTHING_WINDOW_RECIPE
MARKER_RECIPE = V1.MARKER_RECIPE
CEP_BOUNDARY = V1.CEP_BOUNDARY

MARKER_CONTRACT = "phase_endpoint_gap_v2"

# Log-y is a display-only remedy for panels where a first-order endpoint gap
# is visually compressed by a multi-decade dynamic range.  The threshold is
# deliberately explicit so a future renderer cannot silently change an axis
# transform based on subjective inspection.
LOG_Y_DYNAMIC_RANGE_THRESHOLD = 100.0
AXIS_SCALE_POLICY = "log_if_first_order_gap_and_positive_range_ratio_ge_100"

# The first row is the direct-coexistence two-sided contract already used by
# the accepted v1 layer.  The second row deliberately records a different
# semantic from the old phase-kind bracket: the current raw scan changes
# phase_curr between -0.13 and -0.12.  The historical [-0.14,-0.13] bracket
# remains in publication_marker_map.csv as provenance, but is not treated as
# a two-phase endpoint bracket in this display layer.
BOUNDARY_GAP_SPECS = (
    {
        "boundary_id": "mode_a_muB900p0_alpha1p0_direct_coexistence_gap",
        "source_marker_window_id": "publication_mode_a_muB900p0_alpha1p0_direct_coexistence_midpoint",
        "mode_key": "mode_a",
        "plot_panel": "muB900.0",
        "plot_series": "alpha1.0",
        "source_interval_xi_low": "-0.003",
        "source_interval_xi_high": "0.003",
        "gap_xi_low": "-0.003",
        "gap_xi_high": "0.003",
        "gap_basis": "direct_coexistence_phase_curr_switch",
        "reason": "direct-coexistence quark/hadron endpoints; no synthetic midpoint curve or marker",
    },
    {
        "boundary_id": "mode_b_T120p0_muB900p0_phase_curr_switch_gap",
        "source_marker_window_id": "publication_mode_b_T120p0_muB900p0_phase_switch_midpoint",
        "mode_key": "mode_b",
        "plot_panel": "T120.0",
        "plot_series": "muB900.0",
        "source_interval_xi_low": "-0.14",
        "source_interval_xi_high": "-0.13",
        "gap_xi_low": "-0.13",
        "gap_xi_high": "-0.12",
        "gap_basis": "phase_curr_switch_adjacent_raw_points",
        "reason": "old phase-kind bracket is audit-only; the raw scan's quark-to-hadron endpoint switch is [-0.13,-0.12]",
    },
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relpath(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def git_head() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        lines = [line for line in handle if line.strip() and not line.startswith("#")]
    return list(csv.DictReader(lines))


def write_csv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def parse_finite(row: dict[str, str], field: str) -> float:
    value = float(row[field])
    if not math.isfinite(value):
        raise ValueError(f"non-finite {field}: {row[field]}")
    return value


def canonical_xi(value: str | float) -> str:
    return f"{float(value):.10f}"


def phase_label(phase_curr: str) -> str:
    return {
        "quark": "手征恢复相 (quark)",
        "hadron": "手征破缺相 (hadron)",
    }.get(phase_curr, f"相别未判定 ({phase_curr or 'unknown'})")


def phase_legend_label(phase_curr: str) -> str:
    """Use an ASCII legend label because the CI/render font is Latin-only."""
    return {
        "quark": "chiral-restored (quark) endpoint",
        "hadron": "chiral-broken (hadron) endpoint",
    }.get(phase_curr, f"phase-unresolved ({phase_curr or 'unknown'}) endpoint")


def rounded_mev(value: str | float) -> int:
    """Round a display-only MeV label to the nearest integer."""
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"non-finite display MeV value: {value}")
    return int(round(number))


def display_series_label(mode_key: str, row: dict[str, str]) -> str:
    """Return publication typography without changing source numeric fields."""
    if mode_key == "mode_a":
        plot_series = str(row.get("plot_series", ""))
        alpha_token = plot_series.removeprefix("alpha")
        if not alpha_token:
            match = re.search(r"alpha_T=([^,]+)", str(row.get("plot_series_label", "")))
            if match is None:
                raise ValueError(f"cannot recover alpha_T label from row: {row}")
            alpha_token = match.group(1)
        alpha_value = float(alpha_token)
        if not math.isfinite(alpha_value):
            raise ValueError(f"non-finite alpha_T display value: {alpha_token}")
        return rf"$\alpha_T={alpha_value:.1f},\;T={rounded_mev(row['T_MeV'])}\,\mathrm{{MeV}}$"
    if mode_key == "mode_b":
        return rf"$\mu_B={rounded_mev(row['muB_MeV'])}\,\mathrm{{MeV}}$"
    raise ValueError(f"unsupported mode for display label: {mode_key}")


def figure_axis_spec(
    rows: list[dict[str, Any]],
    gaps: list[dict[str, Any]],
) -> dict[str, Any]:
    """Choose a reversible y-axis transform from the displayed data range.

    The log transform is used only when all displayed values are positive, an
    audited first-order gap is present on the panel, and the panel spans at
    least two orders of magnitude.  Otherwise the historical linear scale is
    retained.
    """
    values = [parse_finite(row, "clean_value") for row in rows]
    if not values:
        raise ValueError("cannot choose an axis scale for an empty figure")
    data_min = min(values)
    data_max = max(values)
    dynamic_range_ratio = data_max / data_min if data_min != 0.0 else math.inf
    common = {
        "data_min": data_min,
        "data_max": data_max,
        "dynamic_range_ratio": dynamic_range_ratio,
        "first_order_gap_present": bool(gaps),
        "threshold": LOG_Y_DYNAMIC_RANGE_THRESHOLD,
    }
    if data_min <= 0.0:
        return {
            **common,
            "axis_scale": "linear",
            "axis_scale_reason": "non_positive_display_value",
        }
    if not gaps:
        return {
            **common,
            "axis_scale": "linear",
            "axis_scale_reason": "no_rendered_first_order_gap",
        }
    if dynamic_range_ratio < LOG_Y_DYNAMIC_RANGE_THRESHOLD:
        return {
            **common,
            "axis_scale": "linear",
            "axis_scale_reason": "positive_range_ratio_below_threshold",
        }
    return {
        **common,
        "axis_scale": "log",
        "axis_scale_reason": "first_order_gap_and_positive_range_ratio_ge_100",
    }


def build_boundary_gap_map(
    loaded: dict[str, dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Validate the explicit display gaps and return render/table records."""
    records: list[dict[str, Any]] = []
    table_rows: list[dict[str, Any]] = []
    for spec in BOUNDARY_GAP_SPECS:
        left = V1.current_row(
            loaded, spec["mode_key"], spec["plot_panel"], spec["plot_series"], spec["gap_xi_low"]
        )
        right = V1.current_row(
            loaded, spec["mode_key"], spec["plot_panel"], spec["plot_series"], spec["gap_xi_high"]
        )
        source_left = V1.current_row(
            loaded, spec["mode_key"], spec["plot_panel"], spec["plot_series"], spec["source_interval_xi_low"]
        )
        source_right = V1.current_row(
            loaded, spec["mode_key"], spec["plot_panel"], spec["plot_series"], spec["source_interval_xi_high"]
        )
        if None in (left, right, source_left, source_right):
            raise ValueError(f"boundary gap input is incomplete: {spec['boundary_id']}")
        assert left is not None and right is not None and source_left is not None and source_right is not None
        left_xi = parse_finite(left, "xi")
        right_xi = parse_finite(right, "xi")
        if not left_xi < right_xi:
            raise ValueError(f"boundary gap is not ordered: {spec['boundary_id']}")
        for field in ("T_MeV", "muB_MeV"):
            if not math.isclose(parse_finite(left, field), parse_finite(right, field), rel_tol=0.0, abs_tol=1.0e-9):
                raise ValueError(f"boundary endpoints disagree in {field}: {spec['boundary_id']}")
        if left["phase_curr"] == right["phase_curr"]:
            raise ValueError(
                f"boundary endpoints do not switch phase_curr: {spec['boundary_id']} "
                f"({left['phase_curr']} -> {right['phase_curr']})"
            )
        endpoint_rows = {"left": left, "right": right}
        record = {
            **spec,
            "gap_xi_low": left_xi,
            "gap_xi_high": right_xi,
            "endpoint_rows": endpoint_rows,
            "left_phase_label": phase_label(left["phase_curr"]),
            "right_phase_label": phase_label(right["phase_curr"]),
        }
        records.append(record)
        for endpoint_role, row in endpoint_rows.items():
            table_rows.append(
                {
                    "boundary_id": spec["boundary_id"],
                    "source_marker_window_id": spec["source_marker_window_id"],
                    "mode_key": spec["mode_key"],
                    "plot_panel": spec["plot_panel"],
                    "plot_series": spec["plot_series"],
                    "gap_basis": spec["gap_basis"],
                    "source_interval_xi_low": spec["source_interval_xi_low"],
                    "source_interval_xi_high": spec["source_interval_xi_high"],
                    "gap_xi_low": canonical_xi(left_xi),
                    "gap_xi_high": canonical_xi(right_xi),
                    "endpoint_role": endpoint_role,
                    "xi": canonical_xi(row["xi"]),
                    "phase_curr": row["phase_curr"],
                    "phase_label": phase_label(row["phase_curr"]),
                    "phase_reference_kind": row["phase_reference_kind"],
                    "phase_structure": row["phase_structure"],
                    "quality_flag": row["quality_flag"],
                    "quality_reason": row["quality_reason"],
                    "run_id": row.get("run_id", ""),
                    "render_gap": True,
                    "render_endpoint_marker": True,
                    "canonical_data_modified": False,
                    "reason": spec["reason"],
                }
            )
    return records, table_rows


def build_phase_switch_inventory(
    loaded: dict[str, dict[str, Any]],
    boundary_gaps: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Inventory every raw ``phase_curr`` switch before selecting display gaps.

    ``phase_curr`` can change as a continuation bookkeeping event even when
    the phase reference is crossover.  Keeping those switches in an audit
    table prevents the renderer from silently treating every label change as
    a first-order transition.
    """
    selected = {
        (
            gap["mode_key"],
            gap["plot_panel"],
            gap["plot_series"],
            canonical_xi(gap["gap_xi_low"]),
            canonical_xi(gap["gap_xi_high"]),
        ): gap["boundary_id"]
        for gap in boundary_gaps
    }
    out: list[dict[str, Any]] = []
    for mode_key, payload in loaded.items():
        for (panel, series), curve_rows in sorted(payload["curves"].items()):
            ordered = sorted(curve_rows, key=lambda row: float(row["xi"]))
            for left, right in zip(ordered, ordered[1:]):
                if left["phase_curr"] == right["phase_curr"]:
                    continue
                left_xi = canonical_xi(left["xi"])
                right_xi = canonical_xi(right["xi"])
                key = (mode_key, panel, series, left_xi, right_xi)
                physical_candidate = (
                    left["phase_structure"] == "first_order_possible"
                    and right["phase_structure"] == "first_order_possible"
                    and "first_order" in {left["phase_reference_kind"], right["phase_reference_kind"]}
                )
                out.append(
                    {
                        "mode_key": mode_key,
                        "plot_panel": panel,
                        "plot_series": series,
                        "left_xi": left_xi,
                        "right_xi": right_xi,
                        "left_phase_curr": left["phase_curr"],
                        "right_phase_curr": right["phase_curr"],
                        "left_phase_reference_kind": left["phase_reference_kind"],
                        "right_phase_reference_kind": right["phase_reference_kind"],
                        "left_phase_structure": left["phase_structure"],
                        "right_phase_structure": right["phase_structure"],
                        "physical_first_order_candidate": physical_candidate,
                        "render_gap": selected.get(key, "") != "",
                        "boundary_id": selected.get(key, ""),
                        "verdict": (
                            "rendered_first_order_endpoint_gap"
                            if key in selected
                            else "audit_only_non_first_order_or_unselected_phase_curr_switch"
                        ),
                        "canonical_data_modified": False,
                    }
                )
    return out


def _crosses_gap(left_xi: float, right_xi: float, gaps: list[dict[str, Any]]) -> bool:
    return any(left_xi < gap["gap_xi_high"] and right_xi > gap["gap_xi_low"] for gap in gaps)


def split_curve_segments(rows: list[dict[str, Any]], gaps: list[dict[str, Any]]) -> list[list[dict[str, Any]]]:
    """Split a sorted curve wherever it enters or crosses a first-order gap."""
    ordered = sorted(rows, key=lambda row: float(row["xi"]))
    segments: list[list[dict[str, Any]]] = []
    current: list[dict[str, Any]] = []
    for row in ordered:
        xi = float(row["xi"])
        inside_gap = any(gap["gap_xi_low"] < xi < gap["gap_xi_high"] for gap in gaps)
        if inside_gap:
            if len(current) >= 2:
                segments.append(current)
            current = []
            continue
        if current and _crosses_gap(float(current[-1]["xi"]), xi, gaps):
            if len(current) >= 2:
                segments.append(current)
            current = [row]
        else:
            current.append(row)
    if len(current) >= 2:
        segments.append(current)
    return segments


def render_figures(
    points: list[dict[str, Any]],
    boundary_gaps: list[dict[str, Any]],
    observables: Iterable[str],
) -> tuple[list[Path], list[dict[str, Any]]]:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("matplotlib is required to render publication-clean figures") from exc

    colors = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE"]
    labels = OBSERVABLE_LABELS
    grouped: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in points:
        grouped[(row["mode_key"], row["plot_panel"], row["plot_series"])].append(row)
    gaps_by_curve: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for gap in boundary_gaps:
        gaps_by_curve[(gap["mode_key"], gap["plot_panel"], gap["plot_series"])].append(gap)
    paths: list[Path] = []
    figure_specs: list[dict[str, Any]] = []
    for mode_key in MODE_CONFIG:
        panels = sorted({key[1] for key in grouped if key[0] == mode_key})
        for panel in panels:
            series_names = sorted({key[2] for key in grouped if key[:2] == (mode_key, panel)})
            for observable in observables:
                fig, ax = plt.subplots(figsize=(6.75, 4.6))
                phase_labels_seen: set[str] = set()
                for series_index, series in enumerate(series_names):
                    rows = [
                        row
                        for row in points
                        if row["mode_key"] == mode_key
                        and row["plot_panel"] == panel
                        and row["plot_series"] == series
                        and row["observable"] == observable
                    ]
                    rows.sort(key=lambda row: float(row["xi"]))
                    gaps = gaps_by_curve.get((mode_key, panel, series), [])
                    curve_color = colors[series_index % len(colors)]
                    segments = split_curve_segments(rows, gaps)
                    for segment_index, segment in enumerate(segments):
                        ax.plot(
                            [float(row["xi"]) for row in segment],
                            [float(row["clean_value"]) for row in segment],
                            color=curve_color,
                            linewidth=1.5,
                            label=(
                                display_series_label(mode_key, segment[0])
                                if segment_index == 0
                                else None
                            ),
                        )
                    point_index = {
                        canonical_xi(row["xi"]): row
                        for row in rows
                    }
                    for gap in gaps:
                        for endpoint_role, phase_marker, phase_edge in (
                            ("left", "o", "#3366CC"),
                            ("right", "s", "#CC3311"),
                        ):
                            endpoint = gap["endpoint_rows"][endpoint_role]
                            xi = canonical_xi(endpoint["xi"])
                            point = point_index.get(xi)
                            if point is None:
                                raise ValueError(f"render endpoint missing from point table: {gap['boundary_id']} {xi}")
                            label = phase_label(endpoint["phase_curr"])
                            legend_label = phase_legend_label(endpoint["phase_curr"])
                            if label not in phase_labels_seen:
                                phase_labels_seen.add(label)
                            else:
                                legend_label = None
                            ax.scatter(
                                [float(point["xi"])],
                                [float(point["clean_value"])],
                                marker=phase_marker,
                                s=48,
                                facecolor="white",
                                edgecolor=phase_edge,
                                linewidth=0.9,
                                zorder=6,
                                label=legend_label,
                            )
                figure_rows = [
                    row
                    for row in points
                    if row["mode_key"] == mode_key
                    and row["plot_panel"] == panel
                    and row["observable"] == observable
                ]
                axis_spec = figure_axis_spec(figure_rows, [gap for series in series_names for gap in gaps_by_curve.get((mode_key, panel, series), [])])
                if axis_spec["axis_scale"] == "log":
                    ax.set_yscale("log")
                ax.set_xlabel(r"$\xi$")
                ax.set_ylabel(labels[observable])
                ax.set_xlim(-0.52, 0.52)
                ax.legend(loc="best")
                fig.tight_layout()
                path = FIGURE_DIR / mode_key / f"plot_panel={panel}" / f"{observable}_vs_xi.png"
                path.parent.mkdir(parents=True, exist_ok=True)
                fig.savefig(path, dpi=600, bbox_inches="tight", pad_inches=0.08)
                plt.close(fig)
                paths.append(path)
                figure_specs.append(
                    {
                        "path": path,
                        "mode_key": mode_key,
                        "plot_panel": panel,
                        "observable": observable,
                        **axis_spec,
                    }
                )
    return paths, figure_specs


def suppress_midpoint_markers(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for row in rows:
        item = dict(row)
        item["marker_status"] = "suppressed_phase_switch_midpoint"
        item["marker_semantics"] = (
            "historical interval midpoint retained for provenance only; no CEP or midpoint star is rendered"
        )
        item["render_marker"] = False
        item["reason"] = str(item.get("reason", "")) + "; semantic-v2 suppresses midpoint marker"
        out.append(item)
    return out


def suppress_marker_audit(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for row in rows:
        item = dict(row)
        item["cep_semantics"] = "not_rendered_semantic_v2"
        item["render_marker"] = False
        item["audit_verdict"] = "historical marker retained for provenance; no CEP or phase-switch midpoint rendered"
        out.append(item)
    return out


def claim_ledger_v2(
    inherited_replacements: list[dict[str, Any]],
    smoothing_replacements: list[dict[str, Any]],
    review_adjustments: list[dict[str, Any]],
    boundary_gaps: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    return [
        {
            "claim_id": "PC-V2-001",
            "status": "supported",
            "claim_zh": "publication_clean_v2 仍只消费作者接受的 prod_v2 raw result；本轮只做 solver-free 显示语义派生，不修改 raw CSV、registry 或 v1 快照。",
            "evidence": "tables/input_inventory.csv; manifest.json; figures/plot_manifest.json",
            "scope_limit": "本层不是新的数值收敛或独立 bulk 全局稳定性证明。",
        },
        {
            "claim_id": "PC-V2-002",
            "status": "supported_with_scope_limit",
            "claim_zh": f"{len(boundary_gaps)} 个已审计一阶端点区间不绘制跨区间连续线；两侧端点直接使用当前 raw 点。",
            "evidence": "tables/boundary_gap_map.csv; tables/publication_clean_points.csv; figures/plot_manifest.json",
            "scope_limit": "断线是 phase-guided display 语义，不改变 raw 数据，也不填充区间内输运量。",
        },
        {
            "claim_id": "PC-V2-003",
            "status": "supported_with_scope_limit",
            "claim_zh": "历史 CEP/phase-switch 中点只保留在 publication_marker_map.csv 作为 provenance，不再绘制星标，也不在图例中称为 CEP。",
            "evidence": "tables/publication_marker_map.csv; tables/marker_semantics_audit.csv; figures/plot_manifest.json",
            "scope_limit": "历史 bracket 仍需结合其原始 phase-kind 语义解释；它不是新的零宽度 CEP 结果。",
        },
        {
            "claim_id": "PC-V2-004",
            "status": "supported_with_scope_limit",
            "claim_zh": "端点 phase_curr=quark/hadron 时分别标为手征恢复相/手征破缺相；mode-B 使用 [-0.13,-0.12] 的实际 phase_curr 切换端点，旧 [-0.14,-0.13] 只作历史审计 bracket。",
            "evidence": "tables/boundary_gap_map.csv; data/outputs/results/relaxtime/transport/phase_guided/*/phase_guided_transport_scan.csv",
            "scope_limit": "端点标签来自当前 raw phase_curr 字段，不额外宣称独立 bulk 复核。",
        },
        {
            "claim_id": "PC-V2-005",
            "status": "author_check",
            "claim_zh": "mode-B 相变线左侧的高 eta/s 区域作为论文写作可靠度 caveat 记录；本 publication 图不修改该区域数值，也不将其自动删除。",
            "evidence": "README.md; tables/claim_ledger.csv; current prod_v2 raw scan",
            "scope_limit": "这是解释层边界，不是新的误差条、phase gate 或平滑规则。",
        },
        {
            "claim_id": "PC-V2-006",
            "status": "not_claimed",
            "claim_zh": "本轮不做独立 bulk 全局分支复核；历史 bulk_derivative_branch_audit.csv 继续保留为历史证据。",
            "evidence": "docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_v2_pole_sensitive_rendering/tables/bulk_derivative_branch_audit.csv",
            "scope_limit": "不得把本 publication 图的语义修正写成 bulk/continuation 已重新证明一致。",
        },
        {
            "claim_id": "PC-V2-007",
            "status": "supported_with_scope_limit",
            "claim_zh": "仅对存在已审计一阶断线且正值动态范围达到 100 倍的 figure 使用 log-y；其余图保持线性坐标。",
            "evidence": "figures/plot_manifest.json; README.md",
            "scope_limit": "log-y 是显示坐标变换，不代表相对跳变变大，也不是新的误差模型或数值修正；若出现非正值则强制回退线性坐标。",
        },
        {
            "claim_id": "PC-V2-008",
            "status": "supported_with_scope_limit",
            "claim_zh": "图例使用 α_T、μ_B 等出版排版符号，并将 MeV 显示值四舍五入到个位；原始数值、路径键和 provenance 精度保持不变。",
            "evidence": "figures/plot_manifest.json; tables/input_inventory.csv; tables/publication_clean_points.csv",
            "scope_limit": "仅改变图例/标注字符串，不改变 CSV 数值、phase 标签或任何计算结果。",
        },
    ]


def render_readme(
    inventory: list[dict[str, Any]],
    figure_paths: list[Path],
    figure_specs: list[dict[str, Any]],
    boundary_gaps: list[dict[str, Any]],
    publication_markers: list[dict[str, Any]],
    review_adjustments: list[dict[str, Any]],
    curve_count: int,
) -> str:
    inventory_lines = "\n".join(
        f"| {row['mode_key']} | {row['scan_rows']} | {row['diagnostic_rows']} | `{row['scan_sha256']}` | `{row['diagnostics_sha256']}` |"
        for row in inventory
    )
    axis_counts: dict[str, int] = defaultdict(int)
    for spec in figure_specs:
        axis_counts[str(spec["axis_scale"])] += 1
    axis_summary = ", ".join(f"{key}={axis_counts[key]}" for key in sorted(axis_counts))
    return f"""# Issue #130 RS `publication_clean_v2` 语义修正版

## 目的与边界

本包是 `publication_clean_v1` 的版本化、solver-free 显示语义修正。v1 保持不变；v2 不修改 `data/outputs/results/**`、production registry 或任何 raw CSV，也不调用 equilibrium/transport solver。

本轮处理三件事：

1. 一阶相变端点之间不再用实线桥接；端点使用当前 raw `phase_curr` 标注为手征恢复相（quark）或手征破缺相（hadron）。
2. 旧的 CEP/phase-switch 中点星标不再渲染；中点和旧 bracket 仅保留在 provenance 表中，图例不再出现 `CEP`。
3. 统一 figure-only 的显示规则：存在一阶断线且正值动态范围达到 100 倍的图使用 log-y；其余图保持线性坐标。图例使用 `α_T`/`μ_B` 和四舍五入到个位的 MeV 显示值；CSV 中的原始精度不变。

`mode_b, T=120 MeV, mu_B=900 MeV` 的历史 phase-kind bracket 为 `[-0.14,-0.13]`，但当前 raw 扫描的 `phase_curr` 实际由 quark 切换为 hadron 的相邻端点为 `[-0.13,-0.12]`。v2 用后者断线，并在 `boundary_gap_map.csv` 保留前者作为来源 bracket；这是语义对齐，不是新的 solver 复核。

## 输入 provenance

| mode | scan rows | diagnostic rows | scan SHA256 | diagnostic SHA256 |
| --- | ---: | ---: | --- | --- |
{inventory_lines}

- source case：`{CURRENT_CASE}`；calculation SHA：`{CALCULATION_SHA}`；workflow head：`{WORKFLOW_HEAD_SHA}`。
- source solver 已调用；本次派生 `solver_called=false`。
- v1 快照仍保留；v2 另建目录，避免破坏既有 manifest/hash。
- 本包生成图：{len(figure_paths)} 张 PNG（6 个 panel × {len(DISPLAY_FIELDS)} 个 observable）；曲线索引 {curve_count} 条。
- y 轴策略：{axis_summary}（阈值 `max/min ≥ {LOG_Y_DYNAMIC_RANGE_THRESHOLD:g}`，且必须存在已渲染一阶断线和全为正的显示值）。
- 本轮平滑候选：{len(review_adjustments)} 条；仍是 display-only，raw 值和现有 provenance 不变。

## 断线与端点合同

| boundary | gap | endpoint semantics |
| --- | --- | --- |
| mode-A μB=900, αT=1 | `[-0.003,+0.003]` | quark → 手征恢复相；hadron → 手征破缺相 |
| mode-B T=120, μB=900 | `[-0.13,-0.12]` | quark → 手征恢复相；hadron → 手征破缺相 |

`publication_marker_map.csv` 中的历史 midpoint 行全部 `render_marker=false`；它们不再是 CEP 图形标记。`boundary_gap_map.csv` 是绘图实际使用的断线/端点证据表。
`phase_switch_inventory.csv` 列出 raw `phase_curr` 的全部相邻切换，并区分 crossover/continuation bookkeeping 与实际纳入 v2 断线的两处一阶端点。

## 写作层 caveat

mode-B 相变线左侧的高 `eta/s` 区域可靠度较低。这一条是论文分析/写作时的解释边界：v2 不把它当作自动误差条，不自动删点，也不修改数值。若论文需要定量误差声明，应另立诊断任务。

## 未做事项

- 本轮不做独立 bulk 全局分支复核；历史 `bulk_derivative_branch_audit.csv` 继续作为历史证据。
- 不把断线或端点标记写回 raw/reference；不生成新的 CEP 数值。
- `manuscript_eligible=false`，待作者审核 v2 图后再决定是否作为公开候选。

## Figure-only 显示策略

- log-y 只改变坐标变换，用于避免高动态范围曲线压缩一阶跳变；它不改变数值、插值或 phase 语义。
- 当前所有 mode-B 输入观测量均为正值；若未来输入含非正值，渲染器自动回退到线性坐标并在 figure manifest 中记录原因。
- mode-A 图例格式为 `α_T=… , T=… MeV`，mode-B 图例格式为 `μ_B=… MeV`；显示温度/化学势四舍五入到个位，目录和表格仍保留原始键。

## 复现

```powershell
python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v2.py
python scripts/analysis/relaxtime/sync_phase_guided_publication_clean_figure_layer_v2.py
python -m pytest tests/unit/python/test_phase_guided_publication_clean_v2.py
```
"""


def main() -> None:
    observables = list(DISPLAY_FIELDS)
    loaded, inventory = V1.load_inputs(observables)
    replacement_recipe, marker_recipe, smoothing_windows = V1.load_recipe()
    cep_rows = V1.load_cep_boundary()
    inherited_replacements = V1.build_replacement_map(loaded, replacement_recipe)
    smoothing_replacements, smoothing_window_audit = V1.build_smoothing_window_map(
        loaded, smoothing_windows, inherited_replacements, observables
    )
    review_adjustments = V1.build_review_adjustment_map(loaded)
    replacements = [*inherited_replacements, *smoothing_replacements, *review_adjustments]
    historical_markers = V1.build_marker_map(loaded, marker_recipe)
    cep_slice_audit = V1.build_cep_slice_audit(loaded, historical_markers, cep_rows)
    cep_markers = V1.build_cep_marker_map(loaded, cep_slice_audit, cep_rows, observables)
    publication_markers = suppress_midpoint_markers(
        V1.build_publication_marker_map(loaded, observables=observables)
    )
    marker_semantics = suppress_marker_audit(
        V1.build_marker_semantics_audit(loaded, historical_markers, cep_slice_audit)
    )
    # Historical and midpoint markers are audit rows only in v2.  No stars are
    # passed into the clean point table or renderer.
    points = V1.build_clean_points(loaded, replacements, [], observables)
    curves = V1.build_curve_index(points)
    boundary_gaps, boundary_gap_rows = build_boundary_gap_map(loaded)
    phase_switch_inventory = build_phase_switch_inventory(loaded, boundary_gaps)
    figure_paths, figure_specs = render_figures(points, boundary_gaps, observables)

    input_fields = [
        "mode_key", "mode", "scan_rows", "diagnostic_rows", "failed_rows", "xi_count",
        "scan_sha256", "diagnostics_sha256", "failed_sha256", "manifest_sha256",
        "effective_config_sha256", "calculation_sha", "workflow_head_sha",
        "source_solver_called", "derived_solver_called",
    ]
    write_csv(TABLE_DIR / "input_inventory.csv", inventory, input_fields)
    replacement_fields = [
        "window_id", "mode_key", "plot_panel", "plot_series", "observable", "xi",
        "raw_production_value_current", "recipe_raw_production_value", "recipe_display_value",
        "left_xi", "left_value_current", "right_xi", "right_value_current", "derived_display_value",
        "local_residual", "local_residual_relative", "replacement_method", "recipe_source", "recipe_source_sha256",
        "canonical_data_modified", "display_status", "adjustment_type", "adjustment_reason",
    ]
    write_csv(TABLE_DIR / "replacement_map.csv", replacements, replacement_fields)
    write_csv(TABLE_DIR / "review_adjustment_map.csv", review_adjustments, replacement_fields)
    smoothing_audit_fields = [
        "window_id", "scope", "mode_key", "plot_panel", "plot_series", "xi", "left_xi", "right_xi",
        "diagnostic_observable", "cause", "source_policy", "phase_reference_kinds", "phase_structures",
        "quality_flags", "phase_gate_pass", "action", "reason", "max_abs_residual_relative",
        "recipe_source", "recipe_source_sha256", "canonical_data_modified",
    ]
    write_csv(TABLE_DIR / "display_smoothing_window_audit_v2.csv", smoothing_window_audit, smoothing_audit_fields)
    marker_fields = [
        "window_id", "mode_key", "plot_panel", "plot_series", "recipe_xi", "render_xi", "observable",
        "raw_production_value_current", "recipe_raw_production_value", "marker", "marker_status",
        "marker_semantics", "render_marker", "coexistence_side", "canonical_data_modified", "reason",
    ]
    write_csv(TABLE_DIR / "first_order_marker_map.csv", historical_markers, marker_fields)
    marker_audit_fields = [
        "window_id", "mode_key", "plot_panel", "plot_series", "observable", "render_xi",
        "marker_status", "intended_semantics", "cep_semantics", "phase_reference_kind",
        "phase_structure", "quality_flag", "render_marker", "audit_verdict", "evidence", "canonical_data_modified",
    ]
    write_csv(TABLE_DIR / "marker_semantics_audit.csv", marker_semantics, marker_audit_fields)
    cep_audit_fields = [
        "mode_key", "plot_panel", "plot_series", "T_slice_MeV", "muB_slice_MeV", "cep_table_rows",
        "strict_intersection_count", "strict_intersection", "nearest_T_xi", "nearest_T_midpoint_MeV",
        "nearest_T_muB_MeV", "nearest_T_delta_MeV", "nearest_muB_xi", "nearest_muB_T_midpoint_MeV",
        "nearest_muB_MeV", "nearest_muB_delta_MeV", "marker_action", "reason", "reference_path",
        "reference_sha256", "canonical_data_modified",
    ]
    write_csv(TABLE_DIR / "cep_marker_audit.csv", cep_slice_audit, cep_audit_fields)
    write_csv(TABLE_DIR / "cep_marker_map.csv", cep_markers, marker_fields)
    publication_marker_fields = [
        "window_id", "mode_key", "plot_panel", "plot_series", "observable", "T_MeV", "muB_MeV", "marker",
        "marker_status", "marker_basis", "marker_semantics", "interval_xi_low", "interval_xi_high",
        "interval_width_xi", "render_xi", "left_phase_reference_kind", "right_phase_reference_kind",
        "left_phase_structure", "right_phase_structure", "left_quality_flag", "right_quality_flag",
        "left_value", "right_value", "display_value", "reason", "render_marker", "canonical_data_modified",
    ]
    write_csv(TABLE_DIR / "publication_marker_map.csv", publication_markers, publication_marker_fields)
    boundary_fields = [
        "boundary_id", "source_marker_window_id", "mode_key", "plot_panel", "plot_series", "gap_basis",
        "source_interval_xi_low", "source_interval_xi_high", "gap_xi_low", "gap_xi_high", "endpoint_role",
        "xi", "phase_curr", "phase_label", "phase_reference_kind", "phase_structure", "quality_flag",
        "quality_reason", "run_id", "render_gap", "render_endpoint_marker", "canonical_data_modified", "reason",
    ]
    write_csv(TABLE_DIR / "boundary_gap_map.csv", boundary_gap_rows, boundary_fields)
    phase_switch_fields = [
        "mode_key", "plot_panel", "plot_series", "left_xi", "right_xi", "left_phase_curr",
        "right_phase_curr", "left_phase_reference_kind", "right_phase_reference_kind",
        "left_phase_structure", "right_phase_structure", "physical_first_order_candidate",
        "render_gap", "boundary_id", "verdict", "canonical_data_modified",
    ]
    write_csv(TABLE_DIR / "phase_switch_inventory.csv", phase_switch_inventory, phase_switch_fields)
    point_fields = [
        "mode_key", "mode", "plot_panel", "plot_series", "plot_series_label", "T_MeV", "muB_MeV", "xi",
        "observable", "raw_value", "clean_value", "display_status", "value_source", "phase_structure",
        "phase_reference_kind", "quality_flag", "quality_reason", "run_id", "canonical_data_modified",
    ]
    write_csv(TABLE_DIR / "publication_clean_points.csv", points, point_fields)
    curve_fields = [
        "mode_key", "plot_panel", "plot_series", "observable", "point_count", "replacement_count",
        "marker_count", "xi_min", "xi_max", "canonical_data_modified",
    ]
    write_csv(TABLE_DIR / "curve_index.csv", curves, curve_fields)
    write_csv(
        TABLE_DIR / "claim_ledger.csv",
        claim_ledger_v2(inherited_replacements, smoothing_replacements, review_adjustments, boundary_gaps),
        ["claim_id", "status", "claim_zh", "evidence", "scope_limit"],
    )

    generator_path = Path(__file__).resolve()
    figure_assets = []
    for spec in figure_specs:
        path = spec["path"]
        figure_assets.append(
            {
                "path": relpath(path),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
                "mode_key": spec["mode_key"],
                "plot_panel": spec["plot_panel"],
                "observable": spec["observable"],
                "axis_scale": spec["axis_scale"],
                "axis_scale_reason": spec["axis_scale_reason"],
                "data_min": spec["data_min"],
                "data_max": spec["data_max"],
                "dynamic_range_ratio": spec["dynamic_range_ratio"],
                "first_order_gap_present": spec["first_order_gap_present"],
            }
        )
    axis_counts = defaultdict(int)
    for spec in figure_specs:
        axis_counts[str(spec["axis_scale"])] += 1
    plot_manifest = {
        "schema": "phase_guided_transport_publication_clean_plot_manifest_v2",
        "marker_contract": MARKER_CONTRACT,
        "case": CURRENT_CASE,
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "base_git_commit": git_head(),
        "generator": relpath(generator_path),
        "generator_sha256": sha256_file(generator_path),
        "observables": observables,
        "boundary_gap_count": len(boundary_gaps),
        "boundary_endpoint_count": len(boundary_gap_rows),
        "historical_midpoint_render_count": 0,
        "cep_marker_render_count": 0,
        "axis_scale_policy": AXIS_SCALE_POLICY,
        "log_y_dynamic_range_threshold": LOG_Y_DYNAMIC_RANGE_THRESHOLD,
        "axis_scale_counts": dict(sorted(axis_counts.items())),
        "legend_format_policy": "mode_a uses alpha_T/rounded T in mathtext; mode_b uses mu_B/rounded MeV; source numeric fields unchanged",
        "manuscript_eligible": False,
        "canonical_data_modified": False,
        "rendering_semantics": "prod_v2 display curves split at audited first-order phase endpoints; historical CEP/phase-switch midpoint stars suppressed; endpoint markers use raw phase_curr; no synthetic bridge or midpoint value is rendered; high-range first-order panels use explicit log-y display only",
        "figures": figure_assets,
    }
    write_json(FIGURE_DIR / "plot_manifest.json", plot_manifest)
    readme_path = OUT_DIR / "README.md"
    readme_path.write_text(
        render_readme(
            inventory,
            figure_paths,
            figure_specs,
            boundary_gaps,
            publication_markers,
            review_adjustments,
            len(curves),
        ),
        encoding="utf-8",
    )
    output_paths = [
        readme_path,
        *sorted(TABLE_DIR.glob("*.csv")),
        FIGURE_DIR / "plot_manifest.json",
        *figure_paths,
    ]
    manifest = {
        "schema": "phase_guided_transport_publication_clean_manifest_v2",
        "marker_contract": MARKER_CONTRACT,
        "case": CURRENT_CASE,
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "base_git_commit": git_head(),
        "generator": relpath(generator_path),
        "generator_sha256": sha256_file(generator_path),
        "calculation_sha": CALCULATION_SHA,
        "workflow_head_sha": WORKFLOW_HEAD_SHA,
        "status": "derived_author_review_required",
        "manuscript_eligible": False,
        "canonical_data_modified": False,
        "production_write": False,
        "solver_called": False,
        "source_solver_called": True,
        "source_case": CURRENT_CASE,
        "source_registry_status": "approved_raw_manuscript_ineligible",
        "axis_scale_policy": AXIS_SCALE_POLICY,
        "log_y_dynamic_range_threshold": LOG_Y_DYNAMIC_RANGE_THRESHOLD,
        "axis_scale_counts": dict(sorted(axis_counts.items())),
        "legend_format_policy": "mode_a uses alpha_T/rounded T in mathtext; mode_b uses mu_B/rounded MeV; source numeric fields unchanged",
        "source_inputs": inventory,
        "recipe_inputs": {
            "replacement_path": relpath(REPLACEMENT_RECIPE),
            "replacement_sha256": sha256_file(REPLACEMENT_RECIPE),
            "replacement_rows": len(replacement_recipe),
            "smoothing_window_path": relpath(SMOOTHING_WINDOW_RECIPE),
            "smoothing_window_sha256": sha256_file(SMOOTHING_WINDOW_RECIPE),
            "smoothing_window_rows": len(smoothing_windows),
            "marker_path": relpath(MARKER_RECIPE),
            "marker_sha256": sha256_file(MARKER_RECIPE),
            "marker_rows": len(marker_recipe),
            "cep_boundary_path": relpath(CEP_BOUNDARY),
            "cep_boundary_sha256": sha256_file(CEP_BOUNDARY),
            "cep_boundary_rows": len(cep_rows),
        },
        "derived_counts": {
            "replacement_rows": len(replacements),
            "inherited_replacement_rows": len(inherited_replacements),
            "smoothing_replacement_rows": len(smoothing_replacements),
            "review_adjustment_rows": len(review_adjustments),
            "historical_marker_audit_rows": len(historical_markers),
            "historical_midpoint_render_rows": 0,
            "cep_marker_render_rows": 0,
            "boundary_gap_rows": len(boundary_gaps),
            "boundary_endpoint_rows": len(boundary_gap_rows),
            "phase_switch_inventory_rows": len(phase_switch_inventory),
            "publication_clean_point_rows": len(points),
            "curve_rows": len(curves),
            "figure_count": len(figure_paths),
            "axis_scale_counts": dict(sorted(axis_counts.items())),
        },
        "known_boundaries": [
            "v1 remains unchanged; v2 is a separate semantic-correction derivative",
            "mode_a direct-coexistence gap uses raw xi=-0.003/+0.003 quark/hadron endpoints",
            "mode_b historical [-0.14,-0.13] phase-kind bracket is audit-only; render gap uses phase_curr switch [-0.13,-0.12]",
            "historical CEP/phase-switch midpoint stars are suppressed from all v2 figures",
            "mode_b high eta/s left-of-transition reliability caveat is writing-layer only",
            "independent bulk global branch replay is intentionally out of scope",
        ],
        "outputs": [
            {"path": relpath(path), "sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in output_paths
        ],
    }
    write_json(OUT_DIR / "manifest.json", manifest)
    print(json.dumps({"output": relpath(OUT_DIR), "manifest": relpath(OUT_DIR / "manifest.json"), "figures": len(figure_paths)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
