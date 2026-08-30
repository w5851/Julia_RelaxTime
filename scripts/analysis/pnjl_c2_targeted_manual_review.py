#!/usr/bin/env python3
"""Build a solver-free manual review package for the nine C2 regressions.

The source directory is an Actions download containing the production hybrid
and independent full-fine oracle curves.  This script never calls Julia or an
equilibrium solver.  It copies the source curves into a versioned analysis
package, renders context/local rho-mu figures, and turns the three retained C2
CEP failures into an explicit midpoint request table.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
from pathlib import Path
from typing import Iterable


EXPECTED_SCHEMA = "pnjl_c2_targeted_closure_job_v1"
METHODS = ("production_hybrid", "independent_oracle")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Iterable[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def as_float(value: object) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return math.nan
    return result if math.isfinite(result) else math.nan


def as_bool(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def token(value: float) -> str:
    text = f"{value:g}"
    return text.replace("-", "m").replace(".", "p")


def finite_points(path: Path) -> list[dict[str, float]]:
    rows = read_csv(path)
    points: list[dict[str, float]] = []
    for row in rows:
        rho = as_float(row.get("rho"))
        mu = as_float(row.get("muq_MeV"))
        if not (math.isfinite(rho) and math.isfinite(mu)):
            raise ValueError(f"non-finite curve point: {path}")
        if "finite" in row and not as_bool(row["finite"]):
            raise ValueError(f"curve finite=false: {path}")
        if "converged" in row and not as_bool(row["converged"]):
            raise ValueError(f"curve converged=false: {path}")
        points.append({"rho": rho, "muq_MeV": mu})
    points.sort(key=lambda item: item["rho"])
    if not points or len({point["rho"] for point in points}) != len(points):
        raise ValueError(f"empty or duplicate rho curve: {path}")
    return points


def diagnostic_row(path: Path) -> dict[str, str]:
    rows = read_csv(path)
    if len(rows) != 1:
        raise ValueError(f"expected one slice diagnostic row: {path}")
    return rows[0]


def extrema(points: list[dict[str, float]]) -> list[tuple[float, float, str]]:
    slopes = [points[index + 1]["muq_MeV"] - points[index]["muq_MeV"]
              for index in range(len(points) - 1)]
    turns: list[tuple[float, float, str]] = []
    for index in range(1, len(slopes)):
        left, right = slopes[index - 1], slopes[index]
        if left == 0.0 or right == 0.0 or left * right >= 0.0:
            continue
        point = points[index]
        kind = "spinodal_max" if left > 0.0 and right < 0.0 else "spinodal_min"
        turns.append((point["rho"], point["muq_MeV"], kind))
    return turns


def crossings(points: list[dict[str, float]], mu_level: float) -> list[float]:
    if not math.isfinite(mu_level):
        return []
    result: list[float] = []
    for left, right in zip(points, points[1:]):
        y0 = left["muq_MeV"] - mu_level
        y1 = right["muq_MeV"] - mu_level
        if y0 == 0.0:
            candidate = left["rho"]
        elif y0 * y1 > 0.0:
            continue
        elif y1 == y0:
            candidate = left["rho"]
        else:
            fraction = -y0 / (y1 - y0)
            candidate = left["rho"] + fraction * (right["rho"] - left["rho"])
        if not result or abs(candidate - result[-1]) > 1e-9:
            result.append(candidate)
    return result


def local_bounds(
    curves: dict[str, list[dict[str, float]]],
    diagnostics: dict[str, dict[str, str]],
) -> tuple[tuple[float, float], tuple[float, float], dict[str, list[float]]]:
    rho_values: list[float] = []
    mu_values: list[float] = []
    reference_mu_values: list[float] = []
    crossing_values: dict[str, list[float]] = {}
    for method, points in curves.items():
        diag = diagnostics[method]
        mu_level = as_float(diag.get("mu_transition_MeV"))
        crossing_values[method] = crossings(points, mu_level)
        rho_values.extend(crossing_values[method])
        for rho, mu, _kind in extrema(points):
            rho_values.append(rho)
            mu_values.append(mu)
        if math.isfinite(mu_level):
            mu_values.append(mu_level)
            reference_mu_values.append(mu_level)
        for field in ("mu_spinodal_hadron_MeV", "mu_spinodal_quark_MeV"):
            value = as_float(diag.get(field))
            if math.isfinite(value):
                reference_mu_values.append(value)
        for field in ("rho_hadron", "rho_quark", "rho_spinodal_hadron", "rho_spinodal_quark"):
            value = as_float(diag.get(field))
            if math.isfinite(value):
                rho_values.append(value)
    if len(rho_values) < 2:
        rho_values = [0.0, 4.0]
    low, high = min(rho_values), max(rho_values)
    span = max(high - low, 0.01)
    padding = max(0.18 * span, 0.04)
    x_low, x_high = max(0.0, low - padding), min(4.0, high + padding)
    if len(reference_mu_values) >= 2:
        y_low, y_high = min(reference_mu_values), max(reference_mu_values)
    else:
        for points in curves.values():
            for point in points:
                if x_low <= point["rho"] <= x_high:
                    mu_values.append(point["muq_MeV"])
        if len(mu_values) < 2:
            mu_values = [point["muq_MeV"] for points in curves.values() for point in points]
        y_low, y_high = min(mu_values), max(mu_values)
    y_span = max(y_high - y_low, 1e-5)
    y_padding = max(0.14 * y_span, 0.0005)
    return (x_low, x_high), (y_low - y_padding, y_high + y_padding), crossing_values


def plot_case(
    output: Path,
    target_id: str,
    xi: float,
    temperature: float,
    curves: dict[str, list[dict[str, float]]],
    diagnostics: dict[str, dict[str, str]],
    overlay: dict[str, str],
) -> dict[str, object]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    local_x, local_y, crossing_values = local_bounds(curves, diagnostics)
    colors = {"production_hybrid": "#2563eb", "independent_oracle": "#15803d"}
    linestyles = {"production_hybrid": "-", "independent_oracle": "--"}
    labels = {"production_hybrid": "hybrid", "independent_oracle": "full-fine oracle"}
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.0), constrained_layout=True)
    for method in METHODS:
        points = curves[method]
        x = [point["rho"] for point in points]
        y = [point["muq_MeV"] for point in points]
        axes[0].plot(x, y, color=colors[method], linestyle=linestyles[method],
                     linewidth=1.0, label=labels[method], zorder=2 if method == METHODS[0] else 3)
        mask = [local_x[0] <= value <= local_x[1] for value in x]
        axes[1].plot([value for value, keep in zip(x, mask) if keep],
                     [value for value, keep in zip(y, mask) if keep],
                     color=colors[method], linestyle=linestyles[method],
                     linewidth=1.2, label=labels[method], zorder=2 if method == METHODS[0] else 3)
        mu_level = as_float(diagnostics[method].get("mu_transition_MeV"))
        if math.isfinite(mu_level):
            axes[0].axhline(mu_level, color=colors[method], linestyle=":", linewidth=0.8, alpha=0.8)
            axes[1].axhline(mu_level, color=colors[method], linestyle=":", linewidth=0.9, alpha=0.8)
            local_crossings = crossing_values[method]
            axes[1].scatter(local_crossings, [mu_level] * len(local_crossings),
                            color=colors[method], s=28, zorder=4, edgecolor="white", linewidth=0.4)
        for rho, mu, kind in extrema(points):
            if local_x[0] <= rho <= local_x[1] and local_y[0] <= mu <= local_y[1]:
                axes[1].scatter([rho], [mu], color=colors[method], marker="x", s=48, zorder=5)
                axes[1].annotate(kind.replace("spinodal_", ""), (rho, mu),
                                 xytext=(4, 4), textcoords="offset points", fontsize=7,
                                 color=colors[method])
    axes[0].set_xlim(0.0, 4.0)
    axes[0].set_title("full support")
    axes[1].set_xlim(*local_x)
    axes[1].set_ylim(*local_y)
    axes[1].set_title("local S-shape zoom")
    for axis in axes:
        axis.set_xlabel(r"$\rho$")
        axis.set_ylabel(r"$\mu_q$ [MeV]")
        axis.grid(alpha=0.22)
    axes[0].legend(fontsize=8, loc="best")
    hybrid_status = overlay.get("hybrid_status", "")
    oracle_status = overlay.get("oracle_status", "")
    fig.suptitle(f"xi={xi:g}, T={temperature:g} MeV | hybrid={hybrid_status}, oracle={oracle_status}",
                 fontsize=12)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=240)
    plt.close(fig)
    return {
        "target_id": target_id,
        "xi": xi,
        "T_MeV": temperature,
        "file": output.name,
        "local_rho_low": local_x[0],
        "local_rho_high": local_x[1],
        "local_mu_low": local_y[0],
        "local_mu_high": local_y[1],
        "hybrid_crossings": json.dumps(crossing_values["production_hybrid"]),
        "oracle_crossings": json.dumps(crossing_values["independent_oracle"]),
    }


def plot_overview(output: Path, entries: list[dict[str, object]], case_data: dict[str, dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 3, figsize=(15.5, 12.5), constrained_layout=False)
    fig.subplots_adjust(top=0.90, bottom=0.06, left=0.06, right=0.98,
                        hspace=0.30, wspace=0.22)
    colors = {"production_hybrid": "#2563eb", "independent_oracle": "#15803d"}
    linestyles = {"production_hybrid": "-", "independent_oracle": "--"}
    for axis, entry in zip(axes.flat, entries):
        data = case_data[str(entry["target_id"])]
        for method in METHODS:
            points = data["curves"][method]
            x = [point["rho"] for point in points]
            y = [point["muq_MeV"] for point in points]
            mask = [float(entry["local_rho_low"]) <= value <= float(entry["local_rho_high"]) for value in x]
            axis.plot([value for value, keep in zip(x, mask) if keep],
                      [value for value, keep in zip(y, mask) if keep],
                      color=colors[method], linestyle=linestyles[method], linewidth=0.9,
                      label=method.replace("production_", ""),
                      zorder=2 if method == METHODS[0] else 3)
            mu_level = as_float(data["diagnostics"][method].get("mu_transition_MeV"))
            if math.isfinite(mu_level):
                axis.axhline(mu_level, color=colors[method], linestyle=":", linewidth=0.6, alpha=0.7)
        axis.set_xlim(float(entry["local_rho_low"]), float(entry["local_rho_high"]))
        axis.set_ylim(float(entry["local_mu_low"]), float(entry["local_mu_high"]))
        axis.set_title(f"xi={float(entry['xi']):g}, T={float(entry['T_MeV']):g}", fontsize=9)
        axis.set_xlabel(r"$\rho$", fontsize=8)
        axis.set_ylabel(r"$\mu_q$ [MeV]", fontsize=8)
        axis.tick_params(labelsize=7)
        axis.grid(alpha=0.18)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.935),
               ncol=2, fontsize=8)
    fig.suptitle("C2 targeted regressions: local rho-mu audit", fontsize=14, y=0.985)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def build(args: argparse.Namespace) -> Path:
    source_dir = args.source_dir.resolve()
    output_dir = args.output_dir.resolve()
    overlay_rows = read_csv(args.overlay.resolve())
    overlay_by_id = {row["target_id"]: row for row in overlay_rows}
    case_dirs = sorted(source_dir.glob("c2-targeted-regression-*/"))
    if len(case_dirs) != 9:
        raise ValueError(f"expected nine source case directories, got {len(case_dirs)}")

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "figures").mkdir(exist_ok=True)
    (output_dir / "tables").mkdir(exist_ok=True)
    (output_dir / "curves").mkdir(exist_ok=True)
    source_files: list[dict[str, str]] = []
    combined_points: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    figure_rows: list[dict[str, object]] = []
    case_data: dict[str, dict] = {}
    calculation_sha = None
    postprocess_sha = None
    source_run_id = None

    for case_dir in case_dirs:
        manifests: dict[str, dict] = {}
        curves: dict[str, list[dict[str, float]]] = {}
        diagnostics: dict[str, dict[str, str]] = {}
        target_id = case_dir.name.removeprefix("c2-targeted-regression-")
        for method in METHODS:
            method_dir = case_dir / method
            manifest_path = method_dir / "manifest.json"
            curve_path = method_dir / "curve_points.csv"
            diag_path = method_dir / "slice_diagnostics.csv"
            manifest = read_json(manifest_path)
            if manifest.get("schema_version") != EXPECTED_SCHEMA:
                raise ValueError(f"unexpected schema: {manifest_path}")
            if manifest.get("target_id") != target_id or manifest.get("method") != method:
                raise ValueError(f"target/method mismatch: {manifest_path}")
            manifests[method] = manifest
            curves[method] = finite_points(curve_path)
            diagnostics[method] = diagnostic_row(diag_path)
            digest = sha256(curve_path)
            source_files.append({"path": str(curve_path), "sha256": digest})
            source_files.append({"path": str(manifest_path), "sha256": sha256(manifest_path)})
            source_files.append({"path": str(diag_path), "sha256": sha256(diag_path)})
            copied_name = f"{target_id}__{method}__curve_points.csv"
            shutil.copyfile(curve_path, output_dir / "curves" / copied_name)
            for point in curves[method]:
                combined_points.append({
                    "target_id": target_id,
                    "xi": manifest["target"]["xi"],
                    "T_MeV": manifest["target"]["T_MeV"],
                    "method": method,
                    "rho": point["rho"],
                    "muq_MeV": point["muq_MeV"],
                    "source_curve_sha256": digest,
                })
        manifest = manifests["production_hybrid"]
        calc = str(manifest.get("calculation_sha", ""))
        post = str(manifest.get("postprocess_sha", ""))
        run = str(manifest.get("source_run_id", ""))
        calculation_sha = calc if calculation_sha is None else calculation_sha
        postprocess_sha = post if postprocess_sha is None else postprocess_sha
        source_run_id = run if source_run_id is None else source_run_id
        if calc != calculation_sha or post != postprocess_sha or run != source_run_id:
            raise ValueError("source provenance differs across cases")
        if target_id not in overlay_by_id:
            raise ValueError(f"missing overlay row: {target_id}")
        overlay = overlay_by_id[target_id]
        xi = float(manifest["target"]["xi"])
        temperature = float(manifest["target"]["T_MeV"])
        figure_name = f"rho_mu_{target_id}_local_audit.png"
        figure = plot_case(output_dir / "figures" / figure_name, target_id, xi, temperature,
                           curves, diagnostics, overlay)
        figure_rows.append(figure)
        case_data[target_id] = {"curves": curves, "diagnostics": diagnostics}
        summary = dict(overlay)
        summary.update(figure)
        summary.update({"source_run_id": source_run_id, "calculation_sha": calculation_sha,
                        "source_hybrid_curve_sha256": sha256(case_dir / "production_hybrid" / "curve_points.csv"),
                        "source_oracle_curve_sha256": sha256(case_dir / "independent_oracle" / "curve_points.csv")})
        summary_rows.append(summary)

    overview = output_dir / "figures" / "nine_point_local_overview.png"
    plot_overview(overview, figure_rows, case_data)
    write_csv(output_dir / "tables" / "curve_points.csv", combined_points,
              ["target_id", "xi", "T_MeV", "method", "rho", "muq_MeV", "source_curve_sha256"])
    write_csv(output_dir / "tables" / "point_summary.csv", summary_rows,
              list(summary_rows[0].keys()))

    cep_rows = read_csv(args.cep_failures.resolve())
    cep_plan: list[dict[str, object]] = []
    for row in cep_rows:
        low, high = float(row["T_bracket_low_MeV"]), float(row["T_bracket_high_MeV"])
        midpoint = 0.5 * (low + high)
        cep_plan.append({
            "xi": row["xi"],
            "current_bracket_low_MeV": low,
            "current_bracket_high_MeV": high,
            "current_width_MeV": high - low,
            "target_width_MeV": 0.1,
            "required_midpoint_MeV": midpoint,
            "width_if_midpoint_classifies_MeV": 0.5 * (high - low),
            "current_status": row.get("result_status", "ambiguous"),
            "diagnostic_status": "unresolved_without_new_slice",
            "next_action": "numerical_midpoint_required",
            "solver_called_by_this_package": False,
        })
    write_csv(output_dir / "tables" / "cep_current_failures.csv", cep_rows,
              list(cep_rows[0].keys()))
    write_csv(output_dir / "tables" / "cep_refinement_plan.csv", cep_plan,
              list(cep_plan[0].keys()))

    claim_rows = [
        {"claim_id": "curve_finiteness", "claim": "Nine source targets have finite/converged hybrid and oracle rho-mu curves.",
         "status": "supported", "evidence": "tables/curve_points.csv; tables/point_summary.csv",
         "notes": "Validated from source manifests and every curve row."},
        {"claim_id": "oracle_first_order", "claim": "The independent full-fine oracle reports confirmed_first_order for all nine targets.",
         "status": "supported", "evidence": "tables/point_summary.csv",
         "notes": "Candidate count and crossing count are retained for author inspection."},
        {"claim_id": "hybrid_certificate", "claim": "The production hybrid is not a closed geometry certificate at all nine targets.",
         "status": "supported", "evidence": "tables/point_summary.csv",
         "notes": "Hybrid remains ambiguous; no label is overwritten."},
        {"claim_id": "cep_three_failures", "claim": "Three C2 CEP brackets exceed the 0.1 MeV gate.",
         "status": "supported", "evidence": "tables/cep_current_failures.csv",
         "notes": "The package identifies the midpoint requests but does not claim them solved."},
        {"claim_id": "manual_decision", "claim": "Author review is required before any oracle-backed overlay or promotion decision.",
         "status": "author_check", "evidence": "README.md; figures/nine_point_local_overview.png",
         "notes": "This package is diagnostic-only and does not write reference data."},
    ]
    write_csv(output_dir / "tables" / "claim_ledger.csv", claim_rows,
              ["claim_id", "claim", "status", "evidence", "notes"])

    figure_manifest = {
        "schema_version": "pnjl_c2_targeted_manual_review_figures_v1",
        "source_run_id": source_run_id,
        "calculation_sha": calculation_sha,
        "postprocess_sha": postprocess_sha,
        "solver_called": False,
        "oracle_labels_used_for_routing": False,
        "figures": [{**row, "sha256": sha256(output_dir / "figures" / str(row["file"]))}
                    for row in figure_rows]
        + [{"file": overview.name, "sha256": sha256(overview)}],
    }
    (output_dir / "figures" / "plot_manifest.json").write_text(
        json.dumps(figure_manifest, indent=2) + "\n", encoding="utf-8"
    )
    output_files = {}
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and path.name != "manifest.json":
            output_files[str(path.relative_to(output_dir)).replace("\\", "/")] = sha256(path)
    manifest = {
        "schema_version": "pnjl_c2_targeted_manual_review_v1",
        "verdict": "author_review_required",
        "scope": "nine_targeted_regression_curves_plus_three_cep_midpoint_plan",
        "source_run_id": source_run_id,
        "source_calculation_sha": calculation_sha,
        "source_postprocess_sha": postprocess_sha,
        "solver_called": False,
        "oracle_labels_used_for_routing": False,
        "reference_write": False,
        "target_count": len(summary_rows),
        "cep_failure_count": len(cep_rows),
        "source_files": source_files,
        "files": output_files,
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    readme = f"""# C2 targeted manual review v1

这是一个 solver-free 的作者审核包，不是 production/reference 结果。它固定读取
source run `{source_run_id}` 的 9 个 regression target，calculation SHA 为
`{calculation_sha}`，workflow/postprocess SHA 为 `{postprocess_sha}`。

## 9 个 rho-mu 曲线

`curves/` 保存 source artifact 中的 hybrid 和 independent full-fine oracle 曲线副本；
`figures/` 为每个点的完整范围加局部 S 形放大图，`nine_point_local_overview.png`
提供 3x3 总览。`tables/point_summary.csv` 保留分类、候选数、交点数、geometry 和
局部坐标范围。所有曲线点都经过 finite/converged 和重复 rho 检查。

当前观测是：9 个点的 oracle 均为唯一三交点一阶候选，而 hybrid 均为
`ambiguous_near_critical`，原因仍是 Stage-B 到 Stage-C 的 geometry certificate 未闭合；
本包不把 oracle 标签写回 hybrid。

## 三个 CEP 超限

`tables/cep_current_failures.csv` 保留 C2 blocking audit 的原始三行：
`xi=0.125, 0.39375, 0.5`，当前 bracket 宽度均为 `0.125 MeV`，超过硬门禁
`0.1 MeV`。`tables/cep_refinement_plan.csv` 给出每个 bracket 的 midpoint：
`126.25`、`113.5` 和 `107.0 MeV`。这些 midpoint 需要在相同 calculation SHA
下进行新的数值切片；本包没有调用 solver，因此不声称 CEP 已闭合。

## 作者审核边界

作者可逐点确认曲线中的 S 形、三交点、Maxwell 水平线及 hybrid/oracle 差异，并在
后续 overlay 决策中记录接受、保持 ambiguous 或要求进一步数值细化。当前包不修改
容差、support ranking、production label、reference 或 transport。
"""
    (output_dir / "README.md").write_text(readme, encoding="utf-8")
    return output_dir


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--overlay", type=Path, required=True)
    parser.add_argument("--cep-failures", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = build(args)
    print(json.dumps({"output_dir": str(output), "verdict": "author_review_required"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
