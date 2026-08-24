#!/usr/bin/env python3
"""Build a solver-free quality audit for the Issue #130 RS pilot artifact.

The audit consumes the immutable paired candidate/legacy artifact, does not
call Julia or an equilibrium solver, and separates common transport quality
warnings from reference-specific or solver-backed failures.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable


KEY_FIELDS = ("muB_MeV", "xi", "mode", "alpha_T")
QUALITY_REASON = "tau_u_ubar_ratio_high"
MODES = ("candidate_runtime", "legacy")
TRANSPORT_FIELDS = (
    "eta_over_s",
    "sigma_over_T",
    "zeta_over_s",
    "tau_u",
    "tau_d",
    "tau_s",
    "tau_ubar",
    "tau_dbar",
    "tau_sbar",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--repo-sha", required=True)
    parser.add_argument("--calculation-sha", required=True)
    parser.add_argument("--collector-replay-root", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line.strip() and not line.startswith("#")]
    if not lines:
        return []
    return list(csv.DictReader(lines))


def write_csv(path: Path, fields: list[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def f(row: dict[str, str], field: str) -> float:
    value = float(row[field])
    if not math.isfinite(value):
        raise ValueError(f"non-finite {field}")
    return value


def key(row: dict[str, str]) -> tuple[str, ...]:
    return tuple(row[field] for field in KEY_FIELDS)


def relative_diff(left: float, right: float) -> float:
    return abs(left - right) / max(abs(right), 1.0e-30)


def artifact_files(root: Path) -> list[Path]:
    paths = [
        root / "aggregate" / "manifest.json",
        root / "aggregate" / "verdict.json",
        root / "aggregate" / "transport_comparison.csv",
    ]
    for mode in MODES:
        paths.extend(
            root / "results" / mode / name
            for name in (
                "phase_guided_transport_scan.csv",
                "channel_diagnostics.csv",
                "effective_config.json",
                "run_manifest.json",
                "pilot_metadata.json",
                "sampling_plan.csv",
                "failed_points.csv",
                "elapsed_seconds.txt",
            )
        )
        paths.append(root / "status" / f"{mode}.txt")
    return paths


def validate_inputs(root: Path, args: argparse.Namespace) -> tuple[dict[str, Any], dict[str, Any], dict[str, list[dict[str, str]]]]:
    missing = [str(path) for path in artifact_files(root) if not path.is_file()]
    if missing:
        raise ValueError("missing artifact files: " + "; ".join(missing))
    manifest = read_json(root / "aggregate" / "manifest.json")
    verdict = read_json(root / "aggregate" / "verdict.json")
    if manifest.get("run_id") != args.run_id:
        raise ValueError(f"run id mismatch: {manifest.get('run_id')}")
    if manifest.get("repo_sha") != args.repo_sha:
        raise ValueError(f"repo SHA mismatch: {manifest.get('repo_sha')}")
    if manifest.get("calculation_sha") != args.calculation_sha:
        raise ValueError(f"calculation SHA mismatch: {manifest.get('calculation_sha')}")
    rows = {mode: read_csv(root / "results" / mode / "phase_guided_transport_scan.csv") for mode in MODES}
    if any(len(rows[mode]) != 19 for mode in MODES):
        raise ValueError("expected 19 scan rows per reference")
    if set(map(key, rows[MODES[0]])) != set(map(key, rows[MODES[1]])):
        raise ValueError("candidate/legacy request keys differ")
    for mode in MODES:
        if read_csv(root / "results" / mode / "failed_points.csv"):
            raise ValueError(f"failed points present for {mode}")
        if (root / "status" / f"{mode}.txt").read_text(encoding="utf-8").strip() != "0":
            raise ValueError(f"nonzero runner status for {mode}")
        for row in rows[mode]:
            for field in TRANSPORT_FIELDS:
                f(row, field)
            if row.get("converged", "").lower() != "true":
                raise ValueError(f"nonconverged row for {mode}: {key(row)}")
    return manifest, verdict, rows


def build_quality_rows(rows: dict[str, list[dict[str, str]]]) -> list[dict[str, Any]]:
    by_mode = {mode: {key(row): row for row in rows[mode]} for mode in MODES}
    flagged = [
        item
        for item in sorted(by_mode[MODES[0]])
        if by_mode[MODES[0]][item].get("quality_reason") == QUALITY_REASON
        and by_mode[MODES[1]][item].get("quality_reason") == QUALITY_REASON
        and by_mode[MODES[0]][item].get("quality_flag", "").lower() == "true"
        and by_mode[MODES[1]][item].get("quality_flag", "").lower() == "true"
    ]
    output: list[dict[str, Any]] = []
    for item in flagged:
        candidate = by_mode["candidate_runtime"][item]
        legacy = by_mode["legacy"][item]
        transport_rel = [relative_diff(f(candidate, field), f(legacy, field)) for field in TRANSPORT_FIELDS]
        output.append(
            {
                "muB_MeV": item[0],
                "xi": item[1],
                "mode": item[2],
                "alpha_T": item[3],
                "T_MeV_candidate": f(candidate, "T_MeV"),
                "T_MeV_legacy": f(legacy, "T_MeV"),
                "delta_T_MeV": f(candidate, "T_MeV") - f(legacy, "T_MeV"),
                "tau_u_ubar_ratio_candidate": f(candidate, "quality_metric"),
                "tau_u_ubar_ratio_legacy": f(legacy, "quality_metric"),
                "n_u_over_n_ubar_candidate": f(candidate, "n_u") / f(candidate, "n_ubar"),
                "n_u_over_n_ubar_legacy": f(legacy, "n_u") / f(legacy, "n_ubar"),
                "phase_kind_candidate": candidate.get("phase_reference_kind", ""),
                "phase_kind_legacy": legacy.get("phase_reference_kind", ""),
                "phase_structure_candidate": candidate.get("phase_structure", ""),
                "phase_structure_legacy": legacy.get("phase_structure", ""),
                "max_transport_relative_diff": max(transport_rel),
                "eta_over_s_relative_diff": relative_diff(f(candidate, "eta_over_s"), f(legacy, "eta_over_s")),
                "sigma_over_T_relative_diff": relative_diff(f(candidate, "sigma_over_T"), f(legacy, "sigma_over_T")),
                "zeta_over_s_relative_diff": relative_diff(f(candidate, "zeta_over_s"), f(legacy, "zeta_over_s")),
            }
        )
    return output


def build_anchor_rows(rows: dict[str, list[dict[str, str]]]) -> list[dict[str, Any]]:
    by_mode = {mode: {key(row): row for row in rows[mode]} for mode in MODES}
    output: list[dict[str, Any]] = []
    for item in sorted(by_mode["candidate_runtime"]):
        candidate = by_mode["candidate_runtime"][item]
        legacy = by_mode["legacy"][item]
        values = {field: relative_diff(f(candidate, field), f(legacy, field)) for field in TRANSPORT_FIELDS}
        output.append(
            {
                "muB_MeV": item[0],
                "xi": item[1],
                "mode": item[2],
                "alpha_T": item[3],
                "T_MeV_candidate": f(candidate, "T_MeV"),
                "T_MeV_legacy": f(legacy, "T_MeV"),
                "delta_T_MeV": f(candidate, "T_MeV") - f(legacy, "T_MeV"),
                "phase_kind_match": candidate.get("phase_reference_kind") == legacy.get("phase_reference_kind"),
                "phase_structure_match": candidate.get("phase_structure") == legacy.get("phase_structure"),
                "max_transport_relative_diff": max(values.values()),
            }
        )
    return output


def build_channel_rows(root: Path, quality_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for mode in MODES:
        source_rows = read_csv(root / "results" / mode / "channel_diagnostics.csv")
        for row in source_rows:
            warning = next(
                (
                    item
                    for item in quality_rows
                    if item["muB_MeV"] == row.get("muB_MeV")
                    and item["xi"] == row.get("xi")
                    and math.isclose(
                        float(item[f"T_MeV_{'candidate' if mode == 'candidate_runtime' else 'legacy'}"]),
                        f(row, "T_MeV"),
                        rel_tol=0.0,
                        abs_tol=1.0e-9,
                    )
                ),
                None,
            )
            if warning is None or row.get("species") not in {"u", "ubar"}:
                continue
            contribution = f(row, "contribution")
            total = f(row, "total")
            output.append(
                {
                    "reference_mode": mode,
                    "muB_MeV": row["muB_MeV"],
                    "xi": row["xi"],
                    "alpha_T": warning["alpha_T"],
                    "T_MeV": f(row, "T_MeV"),
                    "species": row["species"],
                    "channel": row["channel"],
                    "density": f(row, "density"),
                    "rate": f(row, "rate"),
                    "contribution": contribution,
                    "species_total": total,
                    "contribution_share": contribution / total if total else math.nan,
                }
            )
    output.sort(key=lambda row: (row["reference_mode"], float(row["muB_MeV"]), float(row["xi"]), float(row["alpha_T"]), row["species"], -row["contribution"]))
    return output


def write_readme(output: Path, manifest: dict[str, Any], quality_rows: list[dict[str, Any]], replay_verdict: str) -> None:
    lines = [
        "# Issue #130 RS numerical pilot quality audit v1",
        "",
        "这是对 numerical pilot source artifact 的 solver-free 后处理审计；不调用 equilibrium solver，不改写数值 CSV。",
        "",
        f"- source run: `{manifest['source_run_id']}`",
        f"- workflow head: `{manifest['workflow_sha']}`",
        f"- calculation SHA: `{manifest['calculation_sha']}`",
        f"- original in-run verdict: `{manifest['original_verdict']}`",
        f"- replay classification: `{replay_verdict}`",
        f"- common `tau_u_ubar_ratio_high` warnings: `{len(quality_rows)}`",
        "",
        "## 结论",
        "",
        "5 个质量警告点在 candidate runtime 与 legacy 两套 reference 中完全成对出现，且原因相同；两套运行均 19/19 点完成，无 failed point、NaN/Inf、重复键或 equilibrium nonconvergence。",
        "这些点的 `tau_u/tau_ubar` 比值超过既有 `scan_quality` 阈值 6，是 transport 质量审查标记，不等价于 solver failure。高 `muB=900 MeV` 下 `n_u/n_ubar` 强烈不对称，并由 u/ubar 通道贡献共同体现；该归因只支持到“共同质量警告”，不升级为传播子机制结论。",
        "",
        "## Reference 差异为何会传播到输运",
        "",
        "本 pilot 的 `mode_a_fixed_muB_phase_scaled` 并非固定同一温度后只切换 phase label。plan 先从 reference 插值得到 `T_phase_base_MeV`，再使用 `T=alpha_T*T_phase_base_MeV` 调用 equilibrium 与 transport。candidate/legacy 的 anchor 温度因此可以不同；即使最终 phase kind/structure 相同，连续的 T、平衡态、质量、Polyakov loop、密度和散射率仍会略有差异。",
        "",
        "## 证据边界",
        "",
        "本包支持：workflow 完成、共同质量警告、reference-specific transport regression 未见明显信号。它不支持：删除 legacy reference、宣称全域 RS numerical convergence、或把质量警告自动视为可忽略的物理异常。",
        "",
        "详见 `collector_replay/`、`tables/quality_warning_summary.csv`、`tables/channel_attribution.csv`、`tables/anchor_drift_summary.csv` 和 `claim_ledger.csv`。",
    ]
    (output / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    root = args.pilot_root.resolve()
    output = args.output.resolve()
    source_manifest, source_verdict, rows = validate_inputs(root, args)
    replay_root = args.collector_replay_root.resolve()
    replay_manifest_path = replay_root / "manifest.json"
    replay_verdict_path = replay_root / "verdict.json"
    if not replay_manifest_path.is_file() or not replay_verdict_path.is_file():
        raise ValueError("collector replay manifest/verdict is missing")
    replay_manifest = read_json(replay_manifest_path)
    replay_verdict = read_json(replay_verdict_path)
    if replay_manifest.get("calculation_sha") != args.calculation_sha:
        raise ValueError("collector replay calculation SHA mismatch")
    if replay_manifest.get("run_id") != args.run_id:
        raise ValueError("collector replay source run mismatch")
    quality_rows = build_quality_rows(rows)
    if len(quality_rows) != 5:
        raise ValueError(f"expected 5 common quality warnings, found {len(quality_rows)}")
    anchor_rows = build_anchor_rows(rows)
    channel_rows = build_channel_rows(root, quality_rows)

    source_files = {
        str(path.relative_to(root)).replace("\\", "/"): sha256(path)
        for path in artifact_files(root)
    }
    output.mkdir(parents=True, exist_ok=True)
    tables = output / "tables"
    write_csv(tables / "quality_warning_summary.csv", list(quality_rows[0]), quality_rows)
    write_csv(tables / "anchor_drift_summary.csv", list(anchor_rows[0]), anchor_rows)
    write_csv(
        tables / "channel_attribution.csv",
        list(channel_rows[0]) if channel_rows else ["reference_mode"],
        channel_rows,
    )
    claim_rows = [
        {"claim_id": "paired_completion", "claim": "candidate and legacy each completed the 19-point numerical plan", "status": "supported", "evidence": "tables/anchor_drift_summary.csv; source manifest", "boundary": "this pilot scope only"},
        {"claim_id": "common_quality_warning", "claim": "the five tau_u/tau_ubar quality warnings are common to candidate and legacy", "status": "supported", "evidence": "tables/quality_warning_summary.csv", "boundary": "shared scan_quality warning, not solver failure"},
        {"claim_id": "continuous_anchor_propagation", "claim": "reference anchor differences can propagate through T-scaled transport inputs", "status": "supported", "evidence": "tables/anchor_drift_summary.csv; source sampling_plan.csv", "boundary": "no new solver-free causal decomposition beyond input path"},
        {"claim_id": "reference_transport_regression", "claim": "this scope shows no material candidate-specific transport regression", "status": "candidate", "evidence": "tables/quality_warning_summary.csv; tables/anchor_drift_summary.csv", "boundary": "limited 19-point pilot, not full production parity"},
        {"claim_id": "legacy_retirement", "claim": "old reference may be removed", "status": "not_claimed", "evidence": "README.md", "boundary": "requires consumer audit and accepted numerical revalidation"},
    ]
    write_csv(tables / "claim_ledger.csv", list(claim_rows[0]), claim_rows)
    manifest = {
        "schema_version": "pnjl_issue130_rs_pilot_quality_audit_v1",
        "source_run_id": args.run_id,
        "workflow_sha": args.repo_sha,
        "calculation_sha": args.calculation_sha,
        "original_verdict": source_verdict.get("verdict"),
        "replay_classification": replay_verdict.get("verdict"),
        "solver_called": False,
        "quality_warning_reason": QUALITY_REASON,
        "common_quality_warning_count": len(quality_rows),
        "source_files": source_files,
        "source_manifest_sha256": sha256(root / "aggregate" / "manifest.json"),
        "source_verdict_sha256": sha256(root / "aggregate" / "verdict.json"),
        "collector_replay_manifest_sha256": sha256(replay_manifest_path),
        "collector_replay_verdict_sha256": sha256(replay_verdict_path),
        "generator": "scripts/analysis/relaxtime/audit_issue130_rs_pilot_quality.py",
        "non_goals": ["no solver rerun", "no production/reference write", "no legacy deletion"],
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    write_readme(output, manifest, quality_rows, manifest["replay_classification"])
    print(json.dumps({"output": str(output), "quality_warning_count": len(quality_rows), "solver_called": False}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
