#!/usr/bin/env python3
"""Build the review-only publication_clean_v3 display derivative.

The parent is the existing publication_clean_v2_residual_smoothed package.
This builder applies one explicit, branch-local interpolation recipe to the
four mode-B T=200, muB=0 tau curves at xi=0.36.  Raw values, the original v2
display values, production registries, and solver outputs are never changed.
"""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import importlib.util
import json
import math
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[3]
V2_SCRIPT = ROOT / "scripts" / "analysis" / "relaxtime" / "build_phase_guided_publication_clean_v2.py"
V2_SPEC = importlib.util.spec_from_file_location("phase_guided_publication_clean_v2_v3", V2_SCRIPT)
if V2_SPEC is None or V2_SPEC.loader is None:  # pragma: no cover
    raise RuntimeError(f"unable to load v2 builder: {V2_SCRIPT}")
V2 = importlib.util.module_from_spec(V2_SPEC)
V2_SPEC.loader.exec_module(V2)


PARENT_DIR = (
    ROOT
    / "docs"
    / "analysis"
    / "relaxtime"
    / "phase_guided_transport"
    / "phase_guided_transport_publication_clean_v2_residual_smoothed"
)
PARENT_TABLE_DIR = PARENT_DIR / "tables"
PARENT_POINTS = PARENT_TABLE_DIR / "publication_clean_points.csv"
PARENT_MANIFEST = PARENT_DIR / "manifest.json"
PARENT_PLOT_MANIFEST = PARENT_DIR / "figures" / "plot_manifest.json"
RECIPE = (
    ROOT
    / "docs"
    / "analysis"
    / "relaxtime"
    / "phase_guided_transport"
    / "phase_guided_transport_v3_residual_smoothing"
    / "tables"
    / "tau_xi036_display_adjustments.csv"
)
OUT_DIR = (
    ROOT
    / "docs"
    / "analysis"
    / "relaxtime"
    / "phase_guided_transport"
    / "phase_guided_transport_publication_clean_v3"
)
TABLE_DIR = OUT_DIR / "tables"
FIGURE_DIR = OUT_DIR / "figures"
AUDIT_DIR = OUT_DIR / "audit"
OBSERVABLES = list(V2.DISPLAY_FIELDS)

POINT_FIELDS = [
    "mode_key",
    "mode",
    "plot_panel",
    "plot_series",
    "plot_series_label",
    "T_MeV",
    "muB_MeV",
    "xi",
    "observable",
    "raw_value",
    "v2_clean_value",
    "clean_value",
    "display_status",
    "value_source",
    "phase_structure",
    "phase_reference_kind",
    "quality_flag",
    "quality_reason",
    "run_id",
    "canonical_data_modified",
    "smoothing_window",
    "smoothing_anchor_points",
    "smoothing_fit_method",
    "protected_by_phase_gate",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relpath(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


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


def finite(row: dict[str, Any], field: str) -> float:
    value = float(row[field])
    if not math.isfinite(value):
        raise ValueError(f"non-finite {field}: {row[field]}")
    return value


def fmt(value: float) -> str:
    return format(value, ".17g")


def point_key(row: dict[str, Any]) -> tuple[str, str, str, str, str]:
    return (
        str(row["mode_key"]),
        str(row["plot_panel"]),
        str(row["plot_series"]),
        V2.canonical_xi(row["xi"]),
        str(row["observable"]),
    )


def curve_key(row: dict[str, Any]) -> tuple[str, str, str, str]:
    return (
        str(row["mode_key"]),
        str(row["plot_panel"]),
        str(row["plot_series"]),
        str(row["observable"]),
    )


def local_log_residual(rows: list[dict[str, Any]], index: int) -> float:
    if index <= 0 or index >= len(rows) - 1:
        return math.nan
    values = [finite(rows[item], "clean_value") for item in (index - 1, index, index + 1)]
    if min(values) <= 0.0:
        return math.nan
    return math.log(values[1]) - 0.5 * (math.log(values[0]) + math.log(values[2]))


def load_parent() -> tuple[list[dict[str, str]], dict[str, Any]]:
    if not PARENT_POINTS.is_file() or not PARENT_MANIFEST.is_file():
        raise FileNotFoundError("publication_clean_v2_residual_smoothed must exist before v3")
    manifest = json.loads(PARENT_MANIFEST.read_text(encoding="utf-8"))
    expected = int(manifest["derived_counts"]["publication_clean_point_rows"])
    rows = read_csv(PARENT_POINTS)
    if len(rows) != expected:
        raise ValueError(f"parent point rows {len(rows)} != manifest count {expected}")
    keys: set[tuple[str, str, str, str, str]] = set()
    for row in rows:
        key = point_key(row)
        if key in keys:
            raise ValueError(f"duplicate parent point key: {key}")
        keys.add(key)
        finite(row, "raw_value")
        finite(row, "v2_clean_value")
        finite(row, "clean_value")
    return rows, manifest


def load_recipe() -> list[dict[str, Any]]:
    rows = read_csv(RECIPE)
    required = {
        "adjustment_id",
        "mode_key",
        "plot_panel",
        "plot_series",
        "observable",
        "target_xi",
        "anchor_left_xi",
        "anchor_right_xi",
        "phase_policy",
        "method",
        "reason",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"recipe is missing required columns: {sorted(required)}")
    out: list[dict[str, Any]] = []
    ids: set[tuple[str, str]] = set()
    for row in rows:
        key = (row["adjustment_id"], row["observable"])
        if key in ids:
            raise ValueError(f"duplicate recipe entry: {key}")
        ids.add(key)
        if row["observable"] not in OBSERVABLES:
            raise ValueError(f"unknown observable in recipe: {row['observable']}")
        item = dict(row)
        for field in ("target_xi", "anchor_left_xi", "anchor_right_xi"):
            item[field] = float(row[field])
        if not item["anchor_left_xi"] < item["target_xi"] < item["anchor_right_xi"]:
            raise ValueError(f"target is not interior to anchors: {row}")
        if item["method"] != "linear_interpolation":
            raise ValueError(f"unsupported v3 method: {item['method']}")
        out.append(item)
    return out


def load_raw_context(parent_manifest: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    loaded, inventory = V2.V1.load_inputs(OBSERVABLES)
    expected_inputs = parent_manifest.get("source_inputs", [])
    current_by_mode = {str(row["mode_key"]): row for row in inventory}
    for expected in expected_inputs:
        current = current_by_mode.get(str(expected["mode_key"]))
        if current is None:
            raise ValueError(f"missing current source inventory for {expected['mode_key']}")
        for field in ("scan_sha256", "diagnostics_sha256"):
            if str(current[field]) != str(expected[field]):
                raise ValueError(f"parent input provenance changed for {expected['mode_key']}:{field}")
    gaps, _ = V2.build_boundary_gap_map(loaded)
    return loaded, inventory, gaps


def gaps_for_curve(gaps: list[dict[str, Any]], mode_key: str, panel: str, series: str) -> list[dict[str, Any]]:
    return [
        gap
        for gap in gaps
        if all(
            field not in gap or str(gap[field]) == expected
            for field, expected in (
                ("mode_key", mode_key),
                ("plot_panel", panel),
                ("plot_series", series),
            )
        )
    ]


def interval_crosses_gap(left: float, right: float, gaps: list[dict[str, Any]]) -> bool:
    lo, hi = min(left, right), max(left, right)
    return any(lo < float(gap["gap_xi_high"]) and hi > float(gap["gap_xi_low"]) for gap in gaps)


def point_inside_gap(xi: float, gaps: list[dict[str, Any]]) -> bool:
    return any(float(gap["gap_xi_low"]) < xi < float(gap["gap_xi_high"]) for gap in gaps)


def apply_adjustments(
    points: list[dict[str, str]],
    adjustments: list[dict[str, Any]],
    gaps: list[dict[str, Any]],
) -> tuple[list[dict[str, str]], list[dict[str, Any]]]:
    grouped: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in points:
        grouped[curve_key(row)].append(row)
    for rows in grouped.values():
        rows.sort(key=lambda row: finite(row, "xi"))

    modified: set[tuple[str, str, str, str, str]] = set()
    audit: list[dict[str, Any]] = []
    for adjustment in adjustments:
        key = (
            adjustment["mode_key"],
            adjustment["plot_panel"],
            adjustment["plot_series"],
            adjustment["observable"],
        )
        rows = grouped.get(key)
        if rows is None:
            raise ValueError(f"recipe curve is missing: {key}")
        by_xi = {V2.canonical_xi(row["xi"]): index for index, row in enumerate(rows)}
        target_key = V2.canonical_xi(adjustment["target_xi"])
        left_key = V2.canonical_xi(adjustment["anchor_left_xi"])
        right_key = V2.canonical_xi(adjustment["anchor_right_xi"])
        if target_key not in by_xi or left_key not in by_xi or right_key not in by_xi:
            raise ValueError(f"recipe points are missing for {key}: {left_key}/{target_key}/{right_key}")
        target_index, left_index, right_index = by_xi[target_key], by_xi[left_key], by_xi[right_key]
        if not left_index < target_index < right_index:
            raise ValueError(f"target is not an interior point for {key}")

        curve_gaps = gaps_for_curve(gaps, *key[:3])
        left_xi = finite(rows[left_index], "xi")
        target_xi = finite(rows[target_index], "xi")
        right_xi = finite(rows[right_index], "xi")
        if interval_crosses_gap(left_xi, right_xi, curve_gaps) or point_inside_gap(target_xi, curve_gaps):
            raise ValueError(f"phase gap blocks v3 adjustment for {key}")
        window_rows = rows[left_index : right_index + 1]
        if adjustment["phase_policy"] == "crossover_only":
            if any(
                row.get("phase_reference_kind") != "crossover"
                or str(row.get("phase_structure", "")).startswith("first_order")
                for row in window_rows
            ):
                raise ValueError(f"crossover gate failed for {key}")
        else:
            raise ValueError(f"unsupported phase policy: {adjustment['phase_policy']}")

        point_key_value = point_key(rows[target_index])
        if point_key_value in modified:
            raise ValueError(f"overlapping v3 adjustment: {point_key_value}")
        modified.add(point_key_value)
        old_value = finite(rows[target_index], "clean_value")
        left_value = finite(rows[left_index], "clean_value")
        right_value = finite(rows[right_index], "clean_value")
        fraction = (target_xi - left_xi) / (right_xi - left_xi)
        derived_value = left_value + fraction * (right_value - left_value)
        if not math.isfinite(derived_value) or derived_value <= 0.0:
            raise ValueError(f"invalid v3 display value for {point_key_value}: {derived_value}")
        before_residual = local_log_residual(rows, target_index)
        rows[target_index]["clean_value"] = fmt(derived_value)
        rows[target_index]["display_status"] = "v3_residual_smoothed"
        rows[target_index]["value_source"] = "v2 residual-smoothed value with explicit local linear interpolation"
        rows[target_index]["canonical_data_modified"] = "False"
        rows[target_index]["smoothing_window"] = adjustment["adjustment_id"]
        rows[target_index]["smoothing_anchor_points"] = (
            f"left={V2.canonical_xi(adjustment['anchor_left_xi'])};"
            f"right={V2.canonical_xi(adjustment['anchor_right_xi'])}"
        )
        rows[target_index]["smoothing_fit_method"] = "linear_interpolation_display_only"
        rows[target_index]["protected_by_phase_gate"] = "True"
        after_residual = local_log_residual(rows, target_index)
        relative_change = (derived_value - old_value) / old_value if old_value else math.nan
        audit.append(
            {
                "adjustment_id": adjustment["adjustment_id"],
                "mode_key": adjustment["mode_key"],
                "plot_panel": adjustment["plot_panel"],
                "plot_series": adjustment["plot_series"],
                "observable": adjustment["observable"],
                "target_xi": V2.canonical_xi(adjustment["target_xi"]),
                "anchor_left_xi": V2.canonical_xi(adjustment["anchor_left_xi"]),
                "anchor_right_xi": V2.canonical_xi(adjustment["anchor_right_xi"]),
                "raw_value": rows[target_index]["raw_value"],
                "v2_clean_value": rows[target_index]["v2_clean_value"],
                "parent_residual_smoothed_value": fmt(old_value),
                "derived_v3_display_value": fmt(derived_value),
                "relative_change_from_parent": fmt(relative_change),
                "local_log_residual_before": fmt(before_residual),
                "local_log_residual_after": fmt(after_residual),
                "phase_policy": adjustment["phase_policy"],
                "phase_gate": "passed",
                "method": adjustment["method"],
                "canonical_data_modified": "False",
                "reason": adjustment["reason"],
            }
        )
    return points, audit


def build_curve_index(points: list[dict[str, str]], audit: list[dict[str, Any]]) -> list[dict[str, Any]]:
    changed = defaultdict(int)
    for row in audit:
        changed[(row["mode_key"], row["plot_panel"], row["plot_series"], row["observable"])] += 1
    grouped: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in points:
        grouped[curve_key(row)].append(row)
    out: list[dict[str, Any]] = []
    for key, rows in sorted(grouped.items()):
        out.append(
            {
                "mode_key": key[0],
                "plot_panel": key[1],
                "plot_series": key[2],
                "observable": key[3],
                "point_count": len(rows),
                "v2_replacement_count": sum(row.get("display_status") == "replacement" for row in rows),
                "residual_smoothing_count": sum(row.get("display_status") == "residual_smoothed" for row in rows),
                "v3_adjustment_count": changed[key],
                "xi_min": min(float(row["xi"]) for row in rows),
                "xi_max": max(float(row["xi"]) for row in rows),
                "canonical_data_modified": False,
            }
        )
    return out


def render_audit_figure(parent: list[dict[str, str]], final: list[dict[str, str]], audit: list[dict[str, Any]]) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    parent_by_curve: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    final_by_curve: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in parent:
        parent_by_curve[curve_key(row)].append(row)
    for row in final:
        final_by_curve[curve_key(row)].append(row)
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.5), squeeze=False)
    for axis, item in zip(axes.flat, audit):
        key = (
            item["mode_key"],
            item["plot_panel"],
            item["plot_series"],
            item["observable"],
        )
        before = sorted(parent_by_curve[key], key=lambda row: float(row["xi"]))
        after = sorted(final_by_curve[key], key=lambda row: float(row["xi"]))
        before = [row for row in before if 0.30 <= float(row["xi"]) <= 0.42]
        after = [row for row in after if 0.30 <= float(row["xi"]) <= 0.42]
        axis.plot([float(row["xi"]) for row in before], [float(row["clean_value"]) for row in before], label="parent residual-smoothed", color="#4477AA")
        axis.plot([float(row["xi"]) for row in after], [float(row["clean_value"]) for row in after], label="v3 display", color="#CC3311", linestyle="--")
        axis.scatter([float(item["target_xi"])], [float(item["derived_v3_display_value"])], color="#228833", zorder=5, label="adjusted point")
        axis.set_title(item["observable"])
        axis.set_xlabel(r"$\xi$")
        axis.set_ylabel(r"$\tau$")
        axis.grid(alpha=0.25)
    axes[0, 0].legend(loc="best", fontsize=8)
    fig.suptitle("publication_clean_v3 local tau display adjustment")
    fig.tight_layout()
    path = AUDIT_DIR / "tau_xi036_before_after.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=450, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return path


def claim_ledger(audit: list[dict[str, Any]]) -> list[dict[str, str]]:
    return [
        {
            "claim_id": "PC-V3-001",
            "status": "supported_with_scope_limit",
            "claim_zh": "publication_clean_v3 以 publication_clean_v2_residual_smoothed 为父级显示表，并保留 raw_value 与原始 v2_clean_value。",
            "evidence": "manifest.json; tables/publication_clean_points.csv; tables/v3_tau_adjustment_map.csv",
            "scope_limit": "不是新的 solver 结果、production 数据或传播子正则化。",
        },
        {
            "claim_id": "PC-V3-002",
            "status": "supported_with_scope_limit",
            "claim_zh": f"在同一 crossover phase gate 内，对 {len(audit)} 个 tau 点使用显式三点局部线性插值。",
            "evidence": "tables/v3_tau_adjustment_map.csv; tables/v3_curve_audit.csv; figures/plot_manifest.json",
            "scope_limit": "插值是显示层选择，不证明原始局部结构是数值伪影或传播子小分母机制。",
        },
        {
            "claim_id": "PC-V3-003",
            "status": "supported",
            "claim_zh": "v3 不修改一阶相变 gap、端点或跨分支连接；完整 publication figure 集从 v3 point table 重新渲染。",
            "evidence": "tables/boundary_gap_map.csv; tables/publication_clean_points.csv; figures/plot_manifest.json",
            "scope_limit": "phase 语义沿用父级 v2 residual-smoothed 包。",
        },
        {
            "claim_id": "PC-V3-004",
            "status": "author_check",
            "claim_zh": "v3 图形可供作者审查，但不能替代 raw/v2 数据进行精确输运量、峰值或机制结论。",
            "evidence": "README.md; audit/tau_xi036_before_after.png; tables/v3_tau_adjustment_map.csv",
            "scope_limit": "正式 manuscript eligibility 仍为 false，需作者确认。",
        },
    ]


def render_readme(
    parent_manifest: dict[str, Any],
    recipe: list[dict[str, Any]],
    audit: list[dict[str, Any]],
    figure_count: int,
) -> str:
    rows = "\n".join(
        f"| {item['observable']} | `{item['parent_residual_smoothed_value']}` | `{item['derived_v3_display_value']}` | `{float(item['relative_change_from_parent']):.6g}` | `{item['local_log_residual_before']}` → `{item['local_log_residual_after']}` |"
        for item in audit
    )
    return f"""# Issue #130 RS `publication_clean_v3` review candidate

## 目的与边界

本包是 `publication_clean_v2_residual_smoothed` 的独立 display-only 派生，供作者审查。
它只对 `mode_b / T=200 MeV / muB=0 MeV / xi=0.36` 的四条 tau 曲线做一次显式局部
线性插值；v2、v2 residual-smoothed、raw CSV、production registry 和 solver 输出均不修改。

这不是传播子有限宽度、分母截断、重新求解或新的物理正则化，也不把该局部结构自动归因于
传播子小分母。

## 变换规则

- 左右锚点：`xi=0.35` 与 `xi=0.37`，锚点值保留父级显示值；
- 目标点：`xi=0.36`；按目标点在锚点区间中的线性位置插值；
- phase gate：区间内必须为 `crossover` 且不含 rendered first-order gap；
- 受影响 observable：`tau_u`, `tau_d`, `tau_ubar`, `tau_dbar`；
- `raw_value`、`v2_clean_value` 和父级 `clean_value` 均在独立审计表中保留。

## Provenance

- parent artifact: `{relpath(PARENT_DIR)}`
- parent manifest SHA256: `{sha256_file(PARENT_MANIFEST)}`
- parent plot manifest SHA256: `{sha256_file(PARENT_PLOT_MANIFEST)}`
- adjustment recipe: `{relpath(RECIPE)}`
- adjustment recipe SHA256: `{sha256_file(RECIPE)}`
- derived solver called: `false`
- canonical/raw data modified: `false`
- production write: `false`
- manuscript eligible: `false`

## 局部变换审计

| observable | parent display value | v3 display value | relative change | log residual before → after |
| --- | ---: | ---: | ---: | ---: |
{rows}

当前四条曲线在 `muB=0` 下由对称性几乎重合；本包仍逐 observable 记录变换，避免把对称性
假设隐藏在一个未审计的共享数值中。局部前后对比图见
`audit/tau_xi036_before_after.png`。

## 审查边界

1. v3 只适合图形审查和排版候选，不能替代 raw/v2 数据作精确数值、导数、峰值或机制结论。
2. 若正式论文采用 v3 图，应在内部稿件 provenance 中保留本包及 adjustment map；不得静默覆盖 v2。
3. 一阶 gap、端点和跨分支连接沿用父级语义，未填补任何 gap。
4. `manuscript_eligible=false` 保持到作者明确审查通过为止。

## 复现

```powershell
python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v3.py
python scripts/analysis/relaxtime/sync_phase_guided_publication_clean_figure_layer_v3.py
python -m pytest tests/unit/python/test_phase_guided_publication_clean_v3.py
```

本包生成 publication figures：{figure_count} 张，完整覆盖父级的 72 张 figure 集。
"""


def main() -> None:
    if (OUT_DIR / "manifest.json").exists():
        raise FileExistsError(f"refusing to overwrite completed v3 analysis package: {OUT_DIR}")
    parent_points, parent_manifest = load_parent()
    recipe = load_recipe()
    _, inventory, gaps = load_raw_context(parent_manifest)
    final_points = [dict(row) for row in parent_points]
    final_points, audit = apply_adjustments(final_points, recipe, gaps)
    if len(audit) != len(recipe):
        raise ValueError(f"applied audit rows {len(audit)} != recipe rows {len(recipe)}")

    # An incomplete directory may remain after a rendering failure; all
    # generated tables and figures below are deterministic and are rewritten.
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    for source_table in PARENT_TABLE_DIR.glob("*.csv"):
        shutil.copy2(source_table, TABLE_DIR / source_table.name)

    V2.FIGURE_DIR = FIGURE_DIR
    figure_paths, figure_specs = V2.render_figures(final_points, gaps, OBSERVABLES)
    audit_figure = render_audit_figure(parent_points, final_points, audit)

    write_csv(TABLE_DIR / "publication_clean_points.csv", final_points, POINT_FIELDS)
    curve_fields = [
        "mode_key",
        "plot_panel",
        "plot_series",
        "observable",
        "point_count",
        "v2_replacement_count",
        "residual_smoothing_count",
        "v3_adjustment_count",
        "xi_min",
        "xi_max",
        "canonical_data_modified",
    ]
    write_csv(TABLE_DIR / "curve_index.csv", build_curve_index(final_points, audit), curve_fields)
    adjustment_fields = [
        "adjustment_id",
        "mode_key",
        "plot_panel",
        "plot_series",
        "observable",
        "target_xi",
        "anchor_left_xi",
        "anchor_right_xi",
        "raw_value",
        "v2_clean_value",
        "parent_residual_smoothed_value",
        "derived_v3_display_value",
        "relative_change_from_parent",
        "local_log_residual_before",
        "local_log_residual_after",
        "phase_policy",
        "phase_gate",
        "method",
        "canonical_data_modified",
        "reason",
    ]
    write_csv(TABLE_DIR / "v3_tau_adjustment_map.csv", audit, adjustment_fields)
    write_csv(
        TABLE_DIR / "v3_curve_audit.csv",
        audit,
        adjustment_fields,
    )
    write_csv(
        TABLE_DIR / "claim_ledger.csv",
        claim_ledger(audit),
        ["claim_id", "status", "claim_zh", "evidence", "scope_limit"],
    )

    generated_at = dt.datetime.now(dt.timezone.utc).isoformat()
    axis_counts: dict[str, int] = defaultdict(int)
    figure_assets: list[dict[str, Any]] = []
    for spec in figure_specs:
        axis_counts[str(spec["axis_scale"])] += 1
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

    plot_manifest = {
        "schema": "phase_guided_transport_publication_clean_plot_manifest_v3",
        "marker_contract": "phase_endpoint_gap_v2",
        "case": parent_manifest["case"],
        "generated_at": generated_at,
        "base_git_commit": V2.git_head(),
        "generator": relpath(Path(__file__).resolve()),
        "generator_sha256": sha256_file(Path(__file__).resolve()),
        "source_parent_artifact": relpath(PARENT_DIR),
        "source_parent_manifest_sha256": sha256_file(PARENT_MANIFEST),
        "source_parent_plot_manifest_sha256": sha256_file(PARENT_PLOT_MANIFEST),
        "adjustment_recipe": relpath(RECIPE),
        "adjustment_recipe_sha256": sha256_file(RECIPE),
        "observables": OBSERVABLES,
        "boundary_gap_count": len(gaps),
        "axis_scale_policy": "inherited from parent v2 residual-smoothed renderer",
        "axis_scale_counts": dict(sorted(axis_counts.items())),
        "manuscript_eligible": False,
        "canonical_data_modified": False,
        "production_write": False,
        "solver_called": False,
        "rendering_semantics": "parent residual-smoothed display values plus four explicit local linear tau interpolations; raw and v2 values retained; first-order gaps remain split",
        "changed_point_count": len(audit),
        "figures": figure_assets,
    }
    write_json(FIGURE_DIR / "plot_manifest.json", plot_manifest)

    readme_path = OUT_DIR / "README.md"
    readme_path.write_text(
        render_readme(parent_manifest, recipe, audit, len(figure_paths)),
        encoding="utf-8",
    )

    output_paths = [
        readme_path,
        *sorted(TABLE_DIR.glob("*.csv")),
        FIGURE_DIR / "plot_manifest.json",
        *figure_paths,
        audit_figure,
    ]
    calculation_sha = parent_manifest.get("calculation_sha")
    if calculation_sha is None:
        calculation_sha = parent_manifest["source_inputs"][0].get("calculation_sha")
    workflow_head_sha = parent_manifest.get("workflow_head_sha")
    if workflow_head_sha is None:
        workflow_head_sha = parent_manifest["source_inputs"][0].get("workflow_head_sha")
    package_manifest = {
        "schema": "phase_guided_transport_publication_clean_manifest_v3",
        "marker_contract": "phase_endpoint_gap_v2",
        "case": parent_manifest["case"],
        "generated_at": generated_at,
        "base_git_commit": V2.git_head(),
        "generator": relpath(Path(__file__).resolve()),
        "generator_sha256": sha256_file(Path(__file__).resolve()),
        "source_parent_artifact": relpath(PARENT_DIR),
        "source_parent_manifest_sha256": sha256_file(PARENT_MANIFEST),
        "source_parent_plot_manifest_sha256": sha256_file(PARENT_PLOT_MANIFEST),
        "adjustment_recipe": relpath(RECIPE),
        "adjustment_recipe_sha256": sha256_file(RECIPE),
        "status": "derived_author_review_required",
        "manuscript_eligible": False,
        "canonical_data_modified": False,
        "production_write": False,
        "solver_called": False,
        "source_solver_called": parent_manifest.get("source_solver_called", True),
        "calculation_sha": calculation_sha,
        "workflow_head_sha": workflow_head_sha,
        "source_inputs": inventory,
        "derived_counts": {
            "parent_point_rows": len(parent_points),
            "publication_clean_point_rows": len(final_points),
            "curve_rows": len(build_curve_index(final_points, audit)),
            "adjustment_recipe_rows": len(recipe),
            "adjusted_point_rows": len(audit),
            "publication_figure_count": len(figure_paths),
            "audit_figure_count": 1,
        },
        "adjustment_summary": {
            "target": "mode_b/T200.0/muB0.0/xi=0.36",
            "observables": [row["observable"] for row in audit],
            "method": "linear_interpolation_between_parent_display_anchors",
            "phase_gate": "crossover_only_passed",
        },
        "known_boundaries": [
            "publication_clean_v2 remains unchanged",
            "publication_clean_v2_residual_smoothed remains unchanged",
            "v3 is a review-only display derivative, not a solver or production rerun",
            "raw_value and v2_clean_value remain available for quantitative use",
            "first-order gaps and endpoints remain inherited from the parent package",
            "the local structure is not assigned a propagator-denominator mechanism by this package",
        ],
        "outputs": [
            {"path": relpath(path), "sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in output_paths
        ],
    }
    write_json(OUT_DIR / "manifest.json", package_manifest)
    print(
        json.dumps(
            {
                "output": relpath(OUT_DIR),
                "manifest": relpath(OUT_DIR / "manifest.json"),
                "publication_figures": len(figure_paths),
                "audit_figures": 1,
                "adjusted_points": len(audit),
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
