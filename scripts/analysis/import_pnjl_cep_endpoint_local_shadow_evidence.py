"""Import a full endpoint-local shadow aggregate into a repository evidence package.

This is a solver-free, read-only importer.  It validates the Actions aggregate
and its raw curve index, copies only derived tables and representative figures,
and deliberately leaves ``curve_points.csv`` outside the repository package.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCHEMA = "cep_maxwell_endpoint_local_production_shadow_v4"
EVIDENCE_SCHEMA = "pnjl_cep_endpoint_local_production_shadow_v4_evidence_v1"
TABLES = (
    "slice_metrics.csv",
    "geometry_accuracy.csv",
    "method_costs.csv",
    "actions_costs.csv",
    "cep_accuracy.csv",
    "curve_index.csv",
)
REQUIRED = (*TABLES, "curve_points.csv", "manifest.json")
ENDPOINT_CERTIFICATES = {
    "endpoint_limited_first_order",
    "endpoint_local_geometry_first_order",
}
# These fields are intentionally ``NaN`` for methods/states where the value is
# not applicable (for example, a monotone slice has no Maxwell endpoints).
# They are schema sentinels, not failed numerical samples.  All other numeric
# fields must remain finite; in particular, ``Inf`` is never accepted.
OPTIONAL_NAN_COLUMNS = {
    "support_low",
    "support_high",
    "support_mu_low",
    "support_mu_high",
    "endpoint_lower_bound",
    "endpoint_upper_bound",
    "endpoint_interpolated_rho_hadron",
    "endpoint_anchor_rho",
    "endpoint_left_bracket_low",
    "endpoint_left_bracket_high",
    "endpoint_right_bracket_low",
    "endpoint_right_bracket_high",
    "area_residual",
    "rho_hadron",
    "rho_quark",
    "mu_transition_MeV",
    "mu_spinodal_hadron_MeV",
    "mu_spinodal_quark_MeV",
    "rho_spinodal_hadron",
    "rho_spinodal_quark",
    "T_last_first_order_MeV",
    "muq_last_first_order_MeV",
    "T_first_monotone_MeV",
    "ambiguity_width_T_MeV",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def _find_aggregate(input_dir: Path) -> Path:
    if (input_dir / "manifest.json").is_file():
        return input_dir
    candidates = sorted(path.parent for path in input_dir.rglob("manifest.json"))
    if len(candidates) != 1:
        raise ValueError(f"expected one aggregate manifest under {input_dir}, found {len(candidates)}")
    return candidates[0]


def _validate_derived_tables(source: Path) -> None:
    """Reject non-finite derived values except documented NaN sentinels."""

    for name in TABLES:
        for row_number, row in enumerate(_rows(source / name), start=2):
            for column, value in row.items():
                if value is None or value == "":
                    continue
                try:
                    number = float(value)
                except (TypeError, ValueError):
                    continue
                if math.isinf(number):
                    raise ValueError(f"derived table contains Inf: {name}:{row_number}:{column}")
                if math.isnan(number) and column not in OPTIONAL_NAN_COLUMNS:
                    raise ValueError(
                        f"unexpected NaN sentinel: {name}:{row_number}:{column}"
                    )


def _validate(source: Path, manifest: dict[str, Any]) -> dict[str, Any]:
    missing = [name for name in REQUIRED if not (source / name).is_file()]
    if missing:
        raise ValueError(f"aggregate is missing required files: {missing}")
    if manifest.get("schema_version") != SCHEMA:
        raise ValueError(f"unexpected source schema: {manifest.get('schema_version')}")
    if manifest.get("scope") != "full" or manifest.get("evidence_state") != "final":
        raise ValueError("full evidence import requires scope=full and evidence_state=final")
    if manifest.get("gate", {}).get("verdict") != "full_hybrid_candidate":
        raise ValueError("full evidence import requires verdict=full_hybrid_candidate")
    gate = manifest["gate"]
    for key in ("workflow_contract_errors", "oracle_errors", "classification_errors", "endpoint_errors", "coverage_errors", "performance_errors"):
        if gate.get(key):
            raise ValueError(f"full evidence gate contains {key}: {gate[key]}")

    _validate_derived_tables(source)

    slices = _rows(source / "slice_metrics.csv")
    if len(slices) != 72:
        raise ValueError(f"expected 72 full slice rows, found {len(slices)}")
    slice_keys = {(row.get("xi"), row.get("method"), row.get("T_MeV")) for row in slices}
    if len(slice_keys) != len(slices):
        raise ValueError("duplicate (xi, method, T) slice key")
    methods = {row.get("method") for row in slices}
    if methods != {"production_hybrid", "memoized_dense", "independent_oracle"}:
        raise ValueError(f"unexpected method matrix: {methods}")

    curves = _rows(source / "curve_points.csv")
    curve_keys = {(row.get("xi"), row.get("method"), row.get("T_MeV"), row.get("rho")) for row in curves}
    if len(curve_keys) != len(curves):
        raise ValueError("duplicate raw curve key")
    for row in curves:
        if row.get("converged", "").lower() != "true" or row.get("finite", "").lower() != "true":
            raise ValueError(f"raw curve contains failed point: {row}")
        if not all(math.isfinite(_float(row.get(field))) for field in ("rho", "muq_MeV")):
            raise ValueError(f"raw curve contains non-finite rho/mu: {row}")

    endpoint_rows = [
        row for row in slices
        if row.get("method") == "production_hybrid" and row.get("certificate_type") in ENDPOINT_CERTIFICATES
    ]
    if len(endpoint_rows) != 3:
        raise ValueError(f"expected three endpoint certificates, found {len(endpoint_rows)}")
    for row in endpoint_rows:
        support_low, support_high = _float(row.get("support_low")), _float(row.get("support_high"))
        anchor = _float(row.get("endpoint_anchor_rho"))
        left_low = _float(row.get("endpoint_left_bracket_low"))
        left_high = _float(row.get("endpoint_left_bracket_high"))
        if not (math.isfinite(support_low) and math.isfinite(support_high) and support_low < support_high):
            raise ValueError(f"invalid endpoint support envelope: {row}")
        if not (math.isfinite(anchor) and support_low <= anchor <= support_high):
            raise ValueError(f"endpoint support does not contain anchor: {row}")
        if math.isfinite(left_low) and math.isfinite(left_high) and (support_low > left_low or support_high < left_high):
            raise ValueError(f"endpoint support excludes initial bracket: {row}")
        if _float(row.get("endpoint_refinement_count"),) > 12:
            raise ValueError(f"endpoint refinement cap exceeded: {row}")

    source_hashes = manifest.get("file_sha256", {})
    for name, expected in source_hashes.items():
        path = source / name
        if path.is_file() and _sha256(path) != expected:
            raise ValueError(f"source hash mismatch: {name}")
    return {
        "slice_rows": len(slices),
        "curve_rows": len(curves),
        "curve_sha256": source_hashes.get("curve_points.csv", _sha256(source / "curve_points.csv")),
        "endpoint_certificates": endpoint_rows,
    }


def _write_claims(path: Path) -> None:
    rows = [
        {
            "claim_id": "full_hybrid_gate",
            "status": "pass",
            "claim_zh": "full 24-anchor shadow 的 hybrid 三态分类、数值几何和 workflow gate 全部通过",
            "evidence": "tables/slice_metrics.csv; tables/geometry_accuracy.csv; manifest.json",
            "fields_or_points": "result_status, geometry_converged, gate.*",
            "boundary": "diagnostic candidate; requires author physical review",
        },
        {
            "claim_id": "endpoint_support_contract",
            "status": "pass",
            "claim_zh": "三个 endpoint-local 证书的 support envelope 覆盖初始左 bracket 与 anchor，最终收缩界限单独记录",
            "evidence": "tables/slice_metrics.csv; tables/geometry_accuracy.csv",
            "fields_or_points": "support_low/high, endpoint_lower/upper_bound, endpoint_*_bracket",
            "boundary": "does not add a fourth physical state",
        },
        {
            "claim_id": "solver_cost",
            "status": "pass",
            "claim_zh": "hybrid unique solves、fixed-rho requests、residual/Jacobian 和 runner cost 不高于 dense",
            "evidence": "tables/method_costs.csv; tables/actions_costs.csv",
            "fields_or_points": "method_costs, performance gate",
            "boundary": "runner noise allowance and author audit still apply",
        },
        {
            "claim_id": "deep_overlay",
            "status": "pass",
            "claim_zh": "approved required-three deep overlay 与 calculation SHA、artifact provenance 一致，且未用于 route 选择",
            "evidence": "manifest.json; source_manifest.json",
            "fields_or_points": "oracle_overlay, expected_calculation_sha",
            "boundary": "approved points only",
        },
        {
            "claim_id": "promotion",
            "status": "author_check",
            "claim_zh": "full_hybrid_candidate 尚未自动晋升 production/reference，也未启动 C0/C1/C2 或 transport",
            "evidence": "README.md; AUDIT.md",
            "fields_or_points": "automatic_gate_is_not_promotion",
            "boundary": "requires author review of representative curves, endpoint certificates, and costs",
        },
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def import_evidence(input_dir: Path, output_dir: Path, aggregate_replay_run_id: str) -> dict[str, Any]:
    source = _find_aggregate(input_dir)
    source_manifest = json.loads((source / "manifest.json").read_text(encoding="utf-8"))
    audit = _validate(source, source_manifest)
    output_dir.mkdir(parents=True, exist_ok=True)
    tables_dir = output_dir / "tables"
    figures_dir = output_dir / "figures"
    tables_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)

    for name in TABLES:
        shutil.copy2(source / name, tables_dir / name)
    shutil.copy2(source / "manifest.json", output_dir / "source_manifest.json")
    source_plot_manifest = json.loads((source / "figures" / "plot_manifest.json").read_text(encoding="utf-8"))
    figure_rows = []
    for figure in source_plot_manifest.get("figures", []):
        src = source / "figures" / figure["file"]
        if not src.is_file() or src.stat().st_size == 0:
            raise ValueError(f"missing representative figure: {src}")
        shutil.copy2(src, figures_dir / src.name)
        figure_rows.append({**figure, "sha256": _sha256(figures_dir / src.name)})
    (figures_dir / "plot_manifest.json").write_text(json.dumps({
        "schema_version": "pnjl_cep_endpoint_local_shadow_plot_manifest_v1",
        "source_schema_version": source_plot_manifest.get("schema_version"),
        "source_curve_sha256": audit["curve_sha256"],
        "figures": figure_rows,
        "raw_curve_copy_in_repository": False,
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_claims(tables_dir / "claim_ledger.csv")

    overlay = source_manifest.get("oracle_overlay", {})
    calculation_sha = source_manifest.get("expected_calculation_sha", "")
    postprocess_sha = source_manifest.get("postprocess_sha", "")
    if not aggregate_replay_run_id.isdigit():
        raise ValueError("aggregate replay run id must be numeric")
    readme = f"""# PNJL CEP endpoint-local production shadow v4

本目录是 full 24-anchor endpoint-local shadow 的诊断 evidence 包，不是
reference promotion。`full_hybrid_candidate` 只表示当前 gate 通过，仍需作者审核
代表性 rho–mu 曲线、三个 endpoint 证书和成本，再决定是否进入 C0→C1→C2。

- calculation SHA: `{calculation_sha}`
- postprocess SHA: `{postprocess_sha}`
- numerical source run: `{source_manifest.get('source_run_id', source_manifest.get('run_id', ''))}`
- aggregate replay run: `{aggregate_replay_run_id}`
- approved deep run: `{overlay.get('deep_run_id', '')}`
- verdict: `{source_manifest.get('gate', {}).get('verdict')}`
- evidence state: `{source_manifest.get('evidence_state')}`

仓库只导入聚合表和代表性 PNG；完整 `curve_points.csv` 保留在 Actions/local artifact，
通过 `tables/curve_index.csv`、source manifest 和 SHA 追溯。三态物理语义、Maxwell
容差、geometry gate、equilibrium solver、旧 reference 和 transport 均未被本包改写。
"""
    (output_dir / "README.md").write_text(readme, encoding="utf-8")
    audit_text = f"""# Full endpoint-local shadow audit

输入是 Actions aggregate replay 的 final evidence。校验结果：

- slice rows: `{audit['slice_rows']}`；raw curve rows（仅用于导入前审计）: `{audit['curve_rows']}`；
  raw key 唯一且所有点 finite/converged；
- 派生表中的 `NaN` 仅出现在 schema 明确允许的“不适用字段”（例如单调切片没有
  Maxwell 端点）；所有 `Inf` 及未声明的 NaN 均拒绝；
- 三个 endpoint certificates 的 support envelope 均覆盖初始 left bracket 与 anchor；
- gate verdict: `{source_manifest.get('gate', {}).get('verdict')}`；所有 oracle/classification/
  endpoint/coverage/performance errors 为空；
- 完整 raw curve 不提交仓库，外部 SHA: `{audit['curve_sha256']}`。

边界：这是 diagnostic candidate，不自动晋升 production/reference，不启动 C0/C1/C2、
phase-reference 或 transport。作者需审阅 figures/ 和 tables/claim_ledger.csv 后再决定。
"""
    (output_dir / "AUDIT.md").write_text(audit_text, encoding="utf-8")

    repository_files = {}
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and path.name != "manifest.json":
            repository_files[str(path.relative_to(output_dir)).replace("\\", "/")] = _sha256(path)
    manifest = {
        "schema_version": EVIDENCE_SCHEMA,
        "source_schema_version": source_manifest.get("schema_version"),
        "source_run_id": source_manifest.get("source_run_id", source_manifest.get("run_id", "")),
        "aggregate_replay_run_id": aggregate_replay_run_id,
        "deep_oracle_run_id": overlay.get("deep_run_id", ""),
        "calculation_sha": calculation_sha,
        "postprocess_sha": postprocess_sha,
        "scope": source_manifest.get("scope"),
        "verdict": source_manifest.get("gate", {}).get("verdict"),
        "evidence_state": source_manifest.get("evidence_state"),
        "raw_curve_points": {
            "file": "Actions/local artifact curve_points.csv",
            "sha256": audit["curve_sha256"],
            "rows": audit["curve_rows"],
            "raw_curve_copy_in_repository": False,
        },
        "oracle_overlay": overlay,
        "repository_files_sha256": repository_files,
        "generated_at_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--aggregate-replay-run-id", required=True)
    args = parser.parse_args()
    import_evidence(args.input_dir, args.output_dir, args.aggregate_replay_run_id)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
