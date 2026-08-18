"""Build the versioned Issue #130 current-state evidence freeze.

This is a read-only post-processing script.  It reads the fixed C2 artifacts,
existing repository evidence, and the approved manual CEP overlay, then writes
summary tables and provenance metadata under docs/analysis.  It never calls
Julia or writes production/reference data.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
HISTORICAL_C0_C1_C2_SHA = "ffa816df0a145f73d7490db1ed9ff10c92e017a4"
POSTPROCESS_SHA = "fd359e792a89beb5ab12349bba761dc58ee16761"
C2_RUN = "31862752226"
FULL_SHADOW_RUN = "31714535418"
CEP_OVERLAY_RUN = "31999149922"
TARGETED_REPLAY_RUN = "31944988128"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def finite_value(value: str | None) -> bool:
    if value is None or value == "":
        return True
    return value.lower() not in {"nan", "+nan", "-nan", "inf", "+inf", "-inf"}


def csv_quality(path: Path, key_columns: Iterable[str]) -> dict[str, Any]:
    rows = read_csv(path)
    headers = list(rows[0]) if rows else []
    nonfinite = 0
    converged_false = 0
    seen: set[tuple[str, ...]] = set()
    duplicate_count = 0
    keys = [column for column in key_columns if column in headers]
    for row in rows:
        nonfinite += sum(not finite_value(row.get(header)) for header in headers)
        if row.get("converged", "").lower() == "false":
            converged_false += 1
        if keys:
            key = tuple(row.get(column, "") for column in keys)
            if key in seen:
                duplicate_count += 1
            seen.add(key)
    return {
        "path": str(path),
        "sha256": sha256_file(path),
        "rows": len(rows),
        "headers": headers,
        "nonfinite_tokens": nonfinite,
        "converged_false": converged_false,
        "key_columns": keys,
        "duplicate_key_count": duplicate_count,
    }


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def git_info(root: Path) -> dict[str, Any]:
    def run(*args: str) -> str:
        return subprocess.check_output(
            ["git", *args], cwd=root, text=True, encoding="utf-8", errors="replace"
        ).strip()

    status = run("status", "--short").splitlines()
    return {
        "worktree": str(root),
        "branch": run("branch", "--show-current"),
        "head": run("rev-parse", "HEAD"),
        "dirty": bool(status),
        "status_lines": status,
    }


def rel(path: Path, root: Path) -> str:
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return str(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, default=Path(r"D:\Desktop\Julia_RelaxTime_issue130_artifacts"))
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--generated-at-utc", default=None)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[3]
    output = args.output or root / "docs" / "analysis" / "pnjl" / "phase_reference_current_state_freeze_v1"
    if output.exists() and any(output.iterdir()) and not args.force:
        raise SystemExit(f"refusing to overwrite non-empty output: {output}")
    output.mkdir(parents=True, exist_ok=True)

    c2_reference = args.artifact_root / "stagec_density_v2_c2_20260813_run31862752226" / "reference"
    c2_root = c2_reference.parent
    manual_root = args.artifact_root / "cep_manual_bisection_audit_v2_31999149922"
    targeted_root = args.artifact_root / "c2_targeted_closure_replay_31944988128" / (
        "c2-targeted-closure-issue130_c2_targeted_closure_v1_20260816_r2_replay4-"
        "regression_curves"
    )

    paths = {
        "boundary": c2_reference / "boundary_issue130_stagec_density_v2_c2_grid_tight_20260813.csv",
        "crossover": c2_reference / "crossover_issue130_stagec_density_v2_c2_grid_tight_20260813.csv",
        "cep": c2_reference / "cep_issue130_stagec_density_v2_c2_grid_tight_20260813.csv",
        "grid_convergence": c2_reference / "phase_grid_convergence_issue130_stagec_density_v2_c2_grid_tight_20260813.csv",
        "spinodals": c2_reference / "spinodals_issue130_stagec_density_v2_c2_grid_tight_20260813.csv",
        "phase_diagnostics": c2_reference / "phase_diagnostics_issue130_stagec_density_v2_c2_grid_tight_20260813.json",
        "c2_manifest": c2_reference / "phase_reference_issue130_stagec_density_v2_c2_grid_tight_20260813_manifest.json",
        "c2_validation": c2_root / "validation_report.json",
        "cep_overlay": manual_root / "author_decision_overlay_v2.json",
        "cep_acceptance": manual_root / "cep_bracket_acceptance_v2.csv",
        "targeted_manifest": targeted_root / "manifest.json",
        "targeted_overlay": targeted_root / "classification_overlay.csv",
    }
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise SystemExit("missing required input(s):\n" + "\n".join(missing))

    c2_manifest = read_json(paths["c2_manifest"])
    c2_validation = read_json(paths["c2_validation"])
    phase_diagnostics = read_json(paths["phase_diagnostics"])
    cep_overlay = read_json(paths["cep_overlay"])
    targeted_manifest = read_json(paths["targeted_manifest"])
    targeted_overlay = read_csv(paths["targeted_overlay"])

    quality = {
        "boundary": csv_quality(paths["boundary"], ["xi", "T_MeV"]),
        "crossover": csv_quality(paths["crossover"], ["xi", "mu_MeV"]),
        "cep": csv_quality(paths["cep"], ["xi"]),
        "spinodals": csv_quality(paths["spinodals"], ["xi", "T_MeV"]),
        "grid_convergence": csv_quality(paths["grid_convergence"], ["axis", "xi", "T_MeV", "level", "reason"]),
    }

    c2_summary = [
        {
            "layer": "historical_c0_c1_c2",
            "source": "docs/analysis/pnjl/c2_convergence_audit_v1/",
            "runs": "31149826740;31235607046;31258823755",
            "calculation_sha": HISTORICAL_C0_C1_C2_SHA,
            "verdict": "classification_regression",
            "rows_or_points": "9 regression anchors;16 CEP brackets;3 crossover risk points",
            "state": "historical_diagnostic_only",
            "note": "保留旧合同审计，不覆盖当前 3c5f6b3c C2。",
        },
        {
            "layer": "full_hybrid_shadow",
            "source": "docs/analysis/pnjl_cep_endpoint_local_production_shadow_v4_20260813/",
            "runs": "31713534102;31714535418;31710995191",
            "calculation_sha": CALCULATION_SHA,
            "verdict": "full_hybrid_candidate",
            "rows_or_points": "72 slice rows;58760 raw curve rows",
            "state": "diagnostic_candidate",
            "note": "24 anchors classification/endpoint/cost gates pass; not a reference promotion. ",
        },
        {
            "layer": "stagec_c2_grid",
            "source": str(c2_root),
            "runs": C2_RUN,
            "calculation_sha": CALCULATION_SHA,
            "verdict": "diagnostic_only_with_unresolved_grid",
            "rows_or_points": "boundary 6886; crossover 1488; CEP 93; grid 22791/5424 unresolved",
            "state": "finite_primary_curves",
            "note": "unresolved 是 refinement/geometry 记录，不是有效 Maxwell boundary 行。",
        },
        {
            "layer": "raw_unresolved_curve_audit",
            "source": "docs/analysis/pnjl/phase_reference_limited_evidence_audit_v1/",
            "runs": "32013771445;21980679",
            "calculation_sha": CALCULATION_SHA,
            "verdict": "raw_curve_coverage_complete_diagnostic_only",
            "rows_or_points": "763 unresolved rows;761 coordinates;1281 rho rows/coordinate",
            "state": "raw_topology_only",
            "note": "S-shape observation不等于 Maxwell/geometry certificate。",
        },
        {
            "layer": "targeted_nine_regression",
            "source": str(targeted_root),
            "runs": TARGETED_REPLAY_RUN,
            "calculation_sha": CALCULATION_SHA,
            "verdict": "automated_classification_inconclusive",
            "rows_or_points": "9/9 finite;hybrid ambiguous vs oracle first-order",
            "state": "author_visual_review_accepted_only",
            "note": "作者接受代表曲线，不把人工判断覆盖 automated classification_match=false。",
        },
        {
            "layer": "cep_t_manual_overlay",
            "source": str(manual_root),
            "runs": CEP_OVERLAY_RUN,
            "calculation_sha": CALCULATION_SHA,
            "verdict": "closed_by_author_review",
            "rows_or_points": "3 fixed xi slices;all brackets 0.0625 MeV",
            "state": "diagnostic_overlay",
            "note": "只关闭固定 T 方向 bracket，不写回 C2 或 reference。",
        },
    ]

    gate_rows = [
        {"gate": "primary_curve_finite_converged", "status": "pass", "evidence": "C2 boundary/spinodals/crossover", "meaning": "主曲线可用于诊断绘图", "blocking": "no"},
        {"gate": "full_hybrid_24_anchor", "status": "pass_candidate", "evidence": "full_hybrid_shadow_v4 manifest/AUDIT", "meaning": "自动 shadow 候选通过", "blocking": "promotion still blocked"},
        {"gate": "grid_refinement_certificate", "status": "unresolved", "evidence": "phase_grid_convergence: 5424 unresolved", "meaning": "不能把所有 unresolved 行当作 Maxwell boundary", "blocking": "yes for strict reference"},
        {"gate": "raw_unresolved_curve_coverage", "status": "diagnostic_complete", "evidence": "phase_reference_limited_evidence_audit_v1", "meaning": "761 条 raw curves finite/converged; topology only", "blocking": "does not close certificate"},
        {"gate": "targeted_nine_automatic_overlay", "status": "inconclusive", "evidence": "classification_overlay.csv: 9/9 match=false", "meaning": "hybrid/oracle 自动分类仍不一致", "blocking": "yes for automatic promotion"},
        {"gate": "targeted_nine_author_visual_review", "status": "accepted_diagnostic", "evidence": "author review of temporary rho-mu plots", "meaning": "可作为人工诊断事实", "blocking": "not a replacement gate"},
        {"gate": "cep_fixed_t_manual_overlay", "status": "pass_diagnostic", "evidence": "author_decision_overlay_v2.json", "meaning": "3 brackets closed at 0.0625 MeV", "blocking": "does not create a single CEP"},
        {"gate": "mu_endpoint_refinement", "status": "not_run", "evidence": "2026-08-18 CEP mu task", "meaning": "等待 solver-free target-list/cost preflight", "blocking": "next planned diagnostic"},
        {"gate": "phase_reference_write", "status": "blocked", "evidence": "task_tracks.toml + claim ledger", "meaning": "当前没有 strict production reference", "blocking": "yes"},
        {"gate": "rs_transport", "status": "blocked", "evidence": "task_tracks.toml", "meaning": "依赖 phase reference promotion", "blocking": "yes"},
    ]

    manual_rows = []
    for row in read_csv(paths["cep_acceptance"]):
        manual_rows.append({
            "xi": row["xi"],
            "first_order_temperatures_MeV": row["first_order_temperatures_MeV"],
            "crossover_side_temperatures_MeV": row["crossover_side_temperatures_MeV"],
            "bracket_low_MeV": row["bracket_low_MeV"],
            "bracket_high_MeV": row["bracket_high_MeV"],
            "bracket_width_MeV": row["bracket_width_MeV"],
            "author_decision": row["author_decision"],
            "status": row["status"],
            "source_evidence": row["source_evidence"],
        })

    targeted_rows = []
    for row in targeted_overlay:
        targeted_rows.append({
            "target_id": row["target_id"],
            "xi": row["xi"],
            "T_MeV": row["T_MeV"],
            "hybrid_status": row["hybrid_status"],
            "oracle_status": row["oracle_status"],
            "classification_match": row["classification_match"],
            "hybrid_crossing_count": row["hybrid_crossing_count"],
            "oracle_crossing_count": row["oracle_crossing_count"],
            "hybrid_geometry_converged": row["hybrid_geometry_converged"],
            "oracle_geometry_converged": row["oracle_geometry_converged"],
            "finite_and_converged": row["finite_and_converged"],
            "author_visual_review": "accepted_diagnostic_only",
            "automatic_gate_state": "inconclusive",
        })

    source_hashes = []
    for name, path in paths.items():
        source_hashes.append({"id": name, "path": str(path), "sha256": sha256_file(path), "bytes": path.stat().st_size})
    repository_sources = [
        "docs/analysis/pnjl_cep_endpoint_local_production_shadow_v4_20260813/manifest.json",
        "docs/analysis/pnjl/c2_surface_views/c2_phase_surfaces_diagnostic_v5_no_triangulation/manifest.json",
        "docs/analysis/pnjl/phase_reference_limited_evidence_audit_v1/manifest.json",
        "docs/analysis/pnjl/c2_convergence_audit_v1/manifest.json",
        "config/governance/task_tracks.toml",
    ]
    for relative in repository_sources:
        path = root / relative
        if path.exists():
            source_hashes.append({"id": relative, "path": str(path), "sha256": sha256_file(path), "bytes": path.stat().st_size})

    plot_entries = [
        {"path": "docs/analysis/pnjl/c2_surface_views/c2_phase_surfaces_diagnostic_v5_no_triangulation/figures/c2_phase_surfaces_mu_xi_T_no_triangulation.png", "role": "C2 global native-support diagnostic", "source_manifest": "docs/analysis/pnjl/c2_surface_views/c2_phase_surfaces_diagnostic_v5_no_triangulation/manifest.json"},
        {"path": "docs/analysis/pnjl_cep_endpoint_local_production_shadow_v4_20260813/figures/rho_mu_xi_0p0_T_130.962_first_order.png", "role": "full-shadow representative first-order curve", "source_manifest": "docs/analysis/pnjl_cep_endpoint_local_production_shadow_v4_20260813/manifest.json"},
        {"path": "docs/analysis/pnjl_cep_endpoint_local_production_shadow_v4_20260813/figures/rho_mu_xi_0p0_T_131.087_first_monotone.png", "role": "full-shadow representative monotone curve", "source_manifest": "docs/analysis/pnjl_cep_endpoint_local_production_shadow_v4_20260813/manifest.json"},
        {"path": "docs/analysis/pnjl/phase_reference_limited_evidence_audit_v1/figures/rho_mu_xi_0p5_T_107_raw_audit.png", "role": "unresolved raw-curve diagnostic", "source_manifest": "docs/analysis/pnjl/phase_reference_limited_evidence_audit_v1/manifest.json"},
    ]
    for entry in plot_entries:
        plot_path = root / entry["path"]
        entry["exists"] = plot_path.exists()
        entry["sha256"] = sha256_file(plot_path) if plot_path.exists() else None

    tables_dir = output / "tables"
    write_csv(tables_dir / "c2_summary.csv", c2_summary, list(c2_summary[0]))
    write_csv(tables_dir / "gate_status_matrix.csv", gate_rows, list(gate_rows[0]))
    write_csv(tables_dir / "cep_manual_overlay.csv", manual_rows, list(manual_rows[0]))
    write_csv(tables_dir / "targeted_regression_status.csv", targeted_rows, list(targeted_rows[0]))
    quality_rows = []
    for name, row in quality.items():
        quality_rows.append({
            "artifact": name,
            "path": row["path"],
            "sha256": row["sha256"],
            "rows": row["rows"],
            "headers": "|".join(row["headers"]),
            "nonfinite_tokens": row["nonfinite_tokens"],
            "converged_false": row["converged_false"],
            "key_columns": "|".join(row["key_columns"]),
            "duplicate_key_count": row["duplicate_key_count"],
            "duplicate_key_note": (
                "refinement event log; repeated (xi,T) and retained exact composite repeats are source semantics"
                if name == "grid_convergence"
                else "declared primary key must be unique"
            ),
        })
    write_csv(tables_dir / "input_validation.csv", quality_rows, list(quality_rows[0]))
    claim_ledger = [
        {"claim_id": "c2_primary_curves", "status": "supported", "claim": "C2 primary boundary, spinodal and crossover rows are finite/converged and key-unique under their declared schemas.", "evidence": "tables/input_validation.csv; external C2 validation_report.json", "boundary": "does not close grid/refinement certificates"},
        {"claim_id": "full_hybrid_candidate", "status": "candidate", "claim": "The 24-anchor endpoint-local shadow is a full_hybrid_candidate.", "evidence": "docs/analysis/pnjl_cep_endpoint_local_production_shadow_v4_20260813/manifest.json", "boundary": "diagnostic candidate, not reference promotion"},
        {"claim_id": "grid_unresolved", "status": "inconclusive", "claim": "5424 C2 grid records remain unresolved at the refinement/geometry layer.", "evidence": "tables/c2_summary.csv; v5/tables/grid_unresolved_diagnostics.csv", "boundary": "unresolved records are not valid Maxwell boundary rows"},
        {"claim_id": "raw_curve_audit", "status": "supported_diagnostic", "claim": "All 761 deduplicated unresolved coordinates have finite/converged 1281-point raw curves.", "evidence": "docs/analysis/pnjl/phase_reference_limited_evidence_audit_v1/manifest.json", "boundary": "raw +→−→+ topology is not a Maxwell/geometry certificate"},
        {"claim_id": "nine_point_overlay", "status": "inconclusive", "claim": "The nine targeted points have hybrid ambiguous versus oracle first-order and classification_match=false.", "evidence": "tables/targeted_regression_status.csv", "boundary": "author visual acceptance is diagnostic only and does not overwrite automated status"},
        {"claim_id": "cep_manual_overlay", "status": "accepted_diagnostic", "claim": "The three fixed T-direction CEP slices have author-reviewed 0.0625 MeV brackets.", "evidence": "tables/cep_manual_overlay.csv; external author_decision_overlay_v2.json", "boundary": "bracket midpoint is not a confirmed single CEP"},
        {"claim_id": "mu_refinement", "status": "not_run", "claim": "The planned CEP mu-direction endpoint refinement has not run.", "evidence": "docs/dev/active/2026-08-18_Issue130_CEP_mu端点补点与派生补全任务单.md", "boundary": "no new solver or derived completion exists"},
        {"claim_id": "phase_reference_promotion", "status": "blocked", "claim": "No current evidence authorizes phase-reference write or RS transport production.", "evidence": "tables/gate_status_matrix.csv; config/governance/task_tracks.toml", "boundary": "full_hybrid_candidate and manual overlays are not promotion"},
    ]
    write_csv(tables_dir / "claim_ledger.csv", claim_ledger, list(claim_ledger[0]))

    figures_dir = output / "figures"
    write_json(figures_dir / "plot_manifest.json", {
        "schema_version": "pnjl_phase_reference_current_state_freeze_plot_manifest_v1",
        "calculation_sha": CALCULATION_SHA,
        "copied_figures": False,
        "entries": plot_entries,
        "note": "冻结包只登记现有图的路径和 hash，不复制或改写正式/诊断图。",
    })

    generated_at = args.generated_at_utc or datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
    README = f"""# Issue #130 当前结果冻结 v1

生成时间：`{generated_at}`。本包是当前 phase-reference 讨论的状态快照，目的
是把已有 C0/C1/C2、full-shadow、C2 grid、raw-curve audit 和人工 CEP overlay
放在同一个可追溯入口下。它不重跑 solver，不修改历史 evidence，不写入
`data/reference/**`，也不把诊断候选升级为生产 reference。

## 固定身份

- 当前工作树：`{git_info(root)['worktree']}`
- 当前分支：`{git_info(root)['branch']}`
- 当前 HEAD：`{git_info(root)['head']}`；工作树在生成时为 dirty，原有改动保持不动。
- 当前数值 calculation SHA：`{CALCULATION_SHA}`。
- 最新 Stage-C C2 source run：`{C2_RUN}`；postprocess SHA：`{POSTPROCESS_SHA}`。
- full endpoint-local shadow aggregate：`{FULL_SHADOW_RUN}`；结果为 `full_hybrid_candidate`。
- 历史 C0/C1/C2 comparator 使用 SHA：`{HISTORICAL_C0_C1_C2_SHA}`，仅作历史对照。

## 已固化事实

1. 最新 C2 primary rows：boundary `6886`、spinodals `6886`、crossover `1488`、
   CEP `93`，输入文件均为 finite，主表 `converged=false` 为零且声明的主键无重复。
   C2 phase diagnostics 的状态计数为 confirmed first-order `6886`、confirmed
   monotone `2636`、ambiguous `359`。
2. `phase_grid_convergence` 共 `22791` 行，其中 `5424` 行 unresolved。它是逐层
   refinement 事件表，不是 Maxwell boundary 表；同一 `(xi,T)` 出现多层记录是设计语义，
   因此 unresolved 行只能保留为诊断，不能直接绘制成有效 Maxwell 面。
3. 独立 raw-curve audit 覆盖 `763` 个 unresolved 行、去重后 `761` 个坐标；每个坐标
   有 `1281` 个 rho 点且 finite/converged。raw 曲线中的 `+→−→+` 只是拓扑观察，
   不等于 Maxwell 面积或 geometry certificate。
4. 九个 targeted regression 点的自动 overlay 仍是 `hybrid ambiguous` 对
   `oracle confirmed_first_order`，`classification_match=false`；作者对临时曲线的
   视觉判断只记录为 `accepted_diagnostic_only`，没有覆盖自动 gate。
5. 三个固定 T 方向 CEP 切片已经人工审核：`xi=0.125`、`0.39375`、`0.5` 的
   bracket 均为 `0.0625 MeV`，状态为 `closed_by_author_review`。这关闭的是固定切片
   的诊断 bracket，不是单值 CEP 写回。
6. μ 方向 CEP endpoint refinement 尚未运行；下一步仍是 solver-free target-list/cost
   preflight。不得把已被排除的 `mu_q > mu_CEP` response peak 当作 crossover 锚点。

## 门禁结论

`full_hybrid_candidate` 说明 endpoint-local shadow 的自动候选 gate 通过，不能等同
于 phase-reference promotion。当前 strict promotion 仍被 grid unresolved 和九点
自动分类不一致阻塞；人工 overlay 只增加诊断事实，不改变 C2 artifact。RS transport
继续依赖 phase-reference promotion，尚未启动。

详见：

- `tables/c2_summary.csv`
- `tables/gate_status_matrix.csv`
- `tables/input_validation.csv`
- `tables/cep_manual_overlay.csv`
- `tables/targeted_regression_status.csv`
- `tables/claim_ledger.csv`
- `figures/plot_manifest.json`
- `checksums.sha256`（生成后写入）

## 明确不在本冻结包中的内容

- 不重跑 C0/C1/C2，不启动 μ numerical pilot、reference write、C3/O1 或 transport。
- 不覆盖 C2 的 unresolved、CEP bracket、Maxwell candidate 或已有图像。
- 不把插值/视觉闭合图当作 strict reference；后续 derived completion 必须单独版本化并
  标明 `interpolated_noncertified`。
"""
    (output / "README.md").write_text(README, encoding="utf-8")

    info = git_info(root)
    c2_artifact_summary = {
        "artifacts": c2_manifest.get("artifacts", {}),
        "config": {
            key: c2_manifest.get("config", {}).get(key)
            for key in (
                "tag", "rho_refinement_policy", "rho_hybrid_endpoint_policy",
                "rho_support_targeted_cap", "rho_support_fine_step",
                "rho_position_tol_MeV", "rho_density_tol", "rho_maxwell_area_tol",
                "temperature_max_refine_level", "xi_max_refine_level",
                "temperature_resolution_target_MeV", "cep_tol_MeV",
            )
        },
        "provenance": c2_manifest.get("provenance", {}),
        "grid_convergence": c2_manifest.get("grid_convergence", {}),
        "status_counts": phase_diagnostics.get("status_counts", {}),
        "stage_counts": phase_diagnostics.get("stage_counts", {}),
        "telemetry_status": phase_diagnostics.get("telemetry_status"),
        "unavailable_record_count": phase_diagnostics.get("unavailable_record_count"),
    }
    c2_validation_summary = {
        key: c2_validation.get(key)
        for key in (
            "tag", "git_commit", "calculation_git_commit", "postprocess_git_commit",
            "source_workflow_run_id", "generated_at_utc", "crossover_rows",
            "converged_rows", "grid_convergence_rows", "grid_unconverged_rows",
        )
    }
    manifest: dict[str, Any] = {
        "schema_version": "pnjl_phase_reference_current_state_freeze_v1",
        "generated_at_utc": generated_at,
        "solver_called_by_freeze": False,
        "reference_write": False,
        "transport_started": False,
        "calculation": {
            "current_sha": CALCULATION_SHA,
            "historical_c0_c1_c2_sha": HISTORICAL_C0_C1_C2_SHA,
            "current_c2_source_run": C2_RUN,
            "current_c2_postprocess_sha": POSTPROCESS_SHA,
            "full_shadow_aggregate_run": FULL_SHADOW_RUN,
            "cep_manual_overlay_run": CEP_OVERLAY_RUN,
            "targeted_regression_replay_run": TARGETED_REPLAY_RUN,
        },
        "repo": info,
        "source_validation": {
            "c2_artifact_summary": c2_artifact_summary,
            "c2_validation_summary": c2_validation_summary,
            "phase_diagnostics_status_counts": phase_diagnostics.get("status_counts", {}),
            "phase_diagnostics_stage_counts": phase_diagnostics.get("stage_counts", {}),
            "cep_overlay_complete": cep_overlay.get("author_review_complete_for_fixed_slices", False),
            "targeted_classification_match": targeted_manifest.get("classification_match"),
            "targeted_finite_and_converged": targeted_manifest.get("finite_and_converged"),
        },
        "sources": source_hashes,
        "output_files": {},
        "verdict": "diagnostic_state_frozen_promotion_blocked",
        "next_action": "solver-free CEP mu-direction target-list/cost preflight",
    }
    output_files: dict[str, str] = {}
    for path in sorted(output.rglob("*")):
        if path.is_file() and path.name != "manifest.json" and path.name != "checksums.sha256":
            output_files[rel(path, output)] = sha256_file(path)
    manifest["output_files"] = output_files
    write_json(output / "manifest.json", manifest)

    checksum_lines = []
    for path in sorted(output.rglob("*")):
        if path.is_file() and path.name != "checksums.sha256":
            checksum_lines.append(f"{sha256_file(path)}  {rel(path, output)}")
    (output / "checksums.sha256").write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(output), "files": len(checksum_lines), "verdict": manifest["verdict"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
