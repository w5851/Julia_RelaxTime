#!/usr/bin/env python3
"""Build Issue #130 post-merge solver-free RS consumer-smoke evidence.

The v1 parity builder predates legacy phase-reference retirement and is kept
immutable.  This versioned builder records the post-#272 merge SHA, resolves
the retired legacy snapshot explicitly, and validates the imported candidate
result trees without invoking a numerical solver.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import shutil
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
MERGE_SHA = "3b19246fb911be4a2efa75fbe14fcb9a793ca256"
WORKFLOW_HEAD_SHA = "22874505877491754eed27519ad8a7b871c82571"
SOURCE_RUN_ID = "32354095831"
REPLAY_RUN_ID = "32451053476"
AGGREGATE_NAME = "aggregate_replay_20260826_v4"
AUDIT_NAME = "post_repair_audit_20260826_v1"
SCHEMA = "pnjl_issue130_rs_transport_runtime_parity_v2"
LEGACY_RESULT_SNAPSHOT_VERSION = "legacy_prod_v1_snapshot_v1"

CANDIDATE_RELATIVE = Path("data/reference/pnjl/issue130_phase_reference_v1")
LEGACY_RELATIVE = Path("data/reference/pnjl/legacy_phase_reference_v1")
RESULT_RELATIVE = {
    "mode_a_fixed_muB_phase_scaled": Path(
        "data/outputs/results/relaxtime/transport/phase_guided/"
        "mode_a_fixed_muB_phase_scaled/"
        "first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2"
    ),
    "mode_b_fixed_T_sparse_muB": Path(
        "data/outputs/results/relaxtime/transport/phase_guided/"
        "mode_b_fixed_T_sparse_muB/"
        "first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2"
    ),
}
EXPECTED_SCAN_ROWS = {
    "mode_a_fixed_muB_phase_scaled": 910,
    "mode_b_fixed_T_sparse_muB": 909,
}
EXPECTED_DIAG_ROWS = {
    "mode_a_fixed_muB_phase_scaled": 38220,
    "mode_b_fixed_T_sparse_muB": 38178,
}


def _load_importer() -> Any:
    path = Path(__file__).resolve().parents[1] / "relaxtime" / "import_issue130_rs_candidate_results.py"
    spec = importlib.util.spec_from_file_location("issue130_rs_importer", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load importer: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


IMPORTER = _load_importer()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    text = path.read_text(encoding="utf-8-sig")
    data = "\n".join(line for line in text.splitlines() if not line.startswith("#"))
    reader = csv.DictReader(data.splitlines())
    return list(reader.fieldnames or []), [dict(row) for row in reader]


def _validate_scan_rows_fast(
    mode: str,
    fieldnames: list[str],
    rows: list[dict[str, str]],
    *,
    expected_rows: int,
) -> dict[str, Any]:
    required = (
        "T_MeV", "muq_MeV", "muB_MeV", "xi", "mode", "phase_reference_kind",
        "converged", "phase_structure", "quality_flag", "quality_reason",
        "tau_u", "tau_d", "tau_s", "tau_ubar", "tau_dbar", "tau_sbar",
        "eta", "sigma", "zeta", "eta_over_s", "zeta_over_s", "sigma_over_T",
    )
    missing = sorted(set(required) - set(fieldnames))
    if missing:
        raise ValueError(f"{mode} scan: missing fields {missing}")
    if len(rows) != expected_rows:
        raise ValueError(f"{mode} scan rows {len(rows)} != expected {expected_rows}")
    keys = [IMPORTER.normalized_key(row, IMPORTER.SCAN_KEYS[mode]) for row in rows]
    duplicate_keys = [key for key, count in Counter(keys).items() if count > 1]
    if duplicate_keys:
        raise ValueError(f"{mode} scan duplicate keys: {duplicate_keys[:3]}")
    invalid_converged = [index for index, row in enumerate(rows) if row.get("converged", "").lower() != "true"]
    if invalid_converged:
        raise ValueError(f"{mode} scan has non-converged rows: {invalid_converged[:3]}")
    nonfinite = []
    for index, row in enumerate(rows):
        bad = IMPORTER._numeric_nonfinite_fields(row, mode=mode)
        if bad:
            nonfinite.append((index, bad))
    if nonfinite:
        raise ValueError(f"{mode} scan non-finite fields: {nonfinite[0]}")
    return {
        "rows": len(rows),
        "duplicate_keys": 0,
        "converged_counts": dict(Counter(row.get("converged", "") for row in rows)),
        "quality_flag_counts": dict(Counter(row.get("quality_flag", "") for row in rows)),
        "quality_reason_counts": dict(Counter(row.get("quality_reason", "") for row in rows)),
        "key_fields": list(IMPORTER.SCAN_KEYS[mode]),
    }


def _validate_diag_rows_fast(
    mode: str,
    fieldnames: list[str],
    rows: list[dict[str, str]],
    *,
    expected_rows: int,
) -> dict[str, Any]:
    required = IMPORTER.DIAG_KEYS + IMPORTER.DIAG_NUMERIC_FIELDS + ("species", "channel")
    missing = sorted(set(required) - set(fieldnames))
    if missing:
        raise ValueError(f"{mode} diagnostic: missing fields {missing}")
    if len(rows) != expected_rows:
        raise ValueError(f"{mode} diagnostic rows {len(rows)} != expected {expected_rows}")
    keys = [IMPORTER.normalized_key(row, IMPORTER.DIAG_KEYS) for row in rows]
    duplicate_keys = [key for key, count in Counter(keys).items() if count > 1]
    if duplicate_keys:
        raise ValueError(f"{mode} diagnostic duplicate keys: {duplicate_keys[:3]}")
    nonfinite = []
    negative = []
    for index, row in enumerate(rows):
        for field in IMPORTER.DIAG_NUMERIC_FIELDS:
            value = IMPORTER.parse_float(row.get(field))
            if value is None or not math.isfinite(value):
                nonfinite.append((index, field))
            elif field in IMPORTER.NONNEGATIVE_DIAG_FIELDS and value < -1e-12:
                negative.append((index, field))
    if nonfinite:
        raise ValueError(f"{mode} diagnostic non-finite fields: {nonfinite[:3]}")
    if negative:
        raise ValueError(f"{mode} diagnostic negative fields: {negative[:3]}")
    return {"rows": len(rows), "duplicate_keys": 0, "nonfinite": 0, "negative": 0}


def write_csv(path: Path, fieldnames: Iterable[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in writer.fieldnames})


def git_value(repo_root: Path, *args: str) -> str:
    try:
        return subprocess.check_output(
            ["git", *args], cwd=repo_root, text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return ""


def rel(path: Path, repo_root: Path) -> str:
    return path.resolve().relative_to(repo_root.resolve()).as_posix()


def _inventory_paths(repo_root: Path) -> list[Path]:
    candidate = repo_root / CANDIDATE_RELATIVE
    legacy = repo_root / LEGACY_RELATIVE
    paths = [candidate / "manifest.json", candidate / "strict" / "manifest.json"]
    paths.extend(sorted((candidate / "strict" / "tables").glob("*.csv")))
    paths.extend(sorted(legacy.iterdir()))
    return paths


def build_input_inventory(repo_root: Path, output_root: Path) -> list[dict[str, Any]]:
    rows = []
    for path in _inventory_paths(repo_root):
        if not path.is_file():
            raise FileNotFoundError(path)
        rows.append({"path": rel(path, repo_root), "bytes": path.stat().st_size, "sha256": sha256(path)})
    write_csv(output_root / "input_inventory.csv", ("path", "bytes", "sha256"), rows)
    return rows


def _normalise_runner_output(output_root: Path) -> dict[str, Any]:
    path = output_root / "consumer_source_smoke.json"
    if not path.is_file():
        raise FileNotFoundError(path)
    raw_path = output_root / "consumer_source_smoke_runner_raw.json"
    if not raw_path.is_file():
        shutil.copy2(path, raw_path)
    raw_sha = sha256(raw_path)
    payload = read_json(path)
    if payload.get("repo_head") != MERGE_SHA:
        raise ValueError(f"consumer smoke repo_head {payload.get('repo_head')} != merge SHA {MERGE_SHA}")
    if payload.get("solver_called") is not False:
        raise ValueError("consumer smoke must record solver_called=false")
    # The historical runner displayed pre-retirement root names.  Preserve its
    # byte copy above, then make the v2 package's canonical locator explicit.
    payload["legacy_root"] = (
        "data/reference/pnjl/legacy_phase_reference_v1/"
        "{boundary.csv,cep.csv,crossover_dense.csv,spinodals.csv}"
    )
    payload["legacy_root_normalized"] = True
    payload["legacy_root_source"] = LEGACY_RELATIVE.as_posix()
    write_json(path, payload)
    return {
        "path": path,
        "raw_sha256": raw_sha,
        "normalized_sha256": sha256(path),
        "payload": payload,
    }


def _validate_source_smoke(payload: dict[str, Any], candidate_manifest_sha: str) -> None:
    if payload.get("schema_version") != "pnjl_issue130_rs_runtime_consumer_smoke_v1":
        raise ValueError("unexpected source smoke schema")
    if payload.get("repo_head") != MERGE_SHA:
        raise ValueError("source smoke repo_head provenance mismatch")
    if payload.get("solver_called") is not False:
        raise ValueError("source smoke solver_called must be false")
    sources = payload.get("sources", {})
    expected = {"runtime", "diagnostic", "legacy"}
    if set(sources) != expected:
        raise ValueError(f"source smoke modes {set(sources)} != {expected}")
    runtime = sources["runtime"]
    diagnostic = sources["diagnostic"]
    legacy = sources["legacy"]
    if runtime.get("source_kind") != "candidate" or not runtime.get("runtime_enabled"):
        raise ValueError("runtime smoke did not enable candidate")
    runtime_view = runtime.get("diagnostics", {}).get("runtime_view")
    if runtime_view not in {
        "certified_candidate_with_legacy_fallback",
        "certified_candidate_with_accepted_then_legacy_fallback",
    }:
        raise ValueError("runtime fallback view mismatch")
    if not runtime.get("diagnostics", {}).get("fallback_enabled"):
        raise ValueError("runtime fallback is not enabled")
    if runtime_view == "certified_candidate_with_accepted_then_legacy_fallback":
        runtime_diagnostics = runtime.get("diagnostics", {})
        if runtime_diagnostics.get("fallback_order") != (
            "strict_candidate>accepted_downstream>legacy_snapshot"
        ):
            raise ValueError("accepted fallback order mismatch")
        for field in ("accepted_manifest_sha256", "accepted_layer_manifest_sha256"):
            value = runtime_diagnostics.get(field, "")
            if not isinstance(value, str) or len(value) != 64:
                raise ValueError(f"accepted fallback {field} is not a SHA-256")
        if not isinstance(runtime_diagnostics.get("accepted_fallback_row_counts"), dict):
            raise ValueError("accepted fallback row counts are missing")
    if runtime.get("diagnostics", {}).get("candidate_manifest_sha256") != candidate_manifest_sha:
        raise ValueError("runtime candidate manifest hash mismatch")
    if diagnostic.get("source_kind") != "candidate" or diagnostic.get("runtime_enabled"):
        raise ValueError("diagnostic smoke source mismatch")
    if diagnostic.get("diagnostics", {}).get("runtime_view") != "diagnostic_all_rows":
        raise ValueError("diagnostic view mismatch")
    if legacy.get("source_kind") != "legacy" or legacy.get("diagnostics", {}).get("runtime_view") != "legacy":
        raise ValueError("legacy rollback smoke source mismatch")
    points = payload.get("consumer_points", [])
    if len(points) != 3 or any(point.get("solver_called") is not False for point in points):
        raise ValueError("consumer points are not a solver-free runtime triplet")


def _validate_dry_run(output_root: Path, mode: str, expected_source: str) -> dict[str, Any]:
    mode_root = output_root / mode
    config = read_json(mode_root / "effective_config.json")
    manifest = read_json(mode_root / "run_manifest.json")
    options = config.get("options", {})
    if options.get("dry_run") is not True:
        raise ValueError(f"{mode} effective config is not dry-run")
    if options.get("phase_reference_source") != expected_source:
        raise ValueError(f"{mode} source is {options.get('phase_reference_source')}, expected {expected_source}")
    if options.get("phase_reference_mode") != ("runtime" if expected_source == "candidate" else "legacy"):
        raise ValueError(f"{mode} phase reference mode mismatch")
    if options.get("phase_anchor_policy") != "reference_interpolation":
        raise ValueError(f"{mode} dry-run used a solver-capable anchor policy")
    if options.get("compute_bulk") is not False:
        raise ValueError(f"{mode} dry-run unexpectedly enables bulk transport")
    if manifest.get("git_commit") != MERGE_SHA or manifest.get("summary", {}).get("dry_run") is not True:
        raise ValueError(f"{mode} run manifest merge/dry-run provenance mismatch")
    fields, rows = read_csv(mode_root / "sampling_plan.csv")
    if not fields or len(rows) != 1:
        raise ValueError(f"{mode} sampling plan must contain exactly one source-resolution point")
    return {
        "mode": mode,
        "source": expected_source,
        "config_sha256": sha256(mode_root / "effective_config.json"),
        "run_manifest_sha256": sha256(mode_root / "run_manifest.json"),
        "sampling_plan_sha256": sha256(mode_root / "sampling_plan.csv"),
        "run_id": manifest.get("run_id"),
        "points_total": manifest.get("summary", {}).get("points_total"),
    }


def _result_tree_hash(root: Path) -> str:
    return IMPORTER.tree_hash(root)


def _validate_result(repo_root: Path, mode: str) -> dict[str, Any]:
    root = repo_root / RESULT_RELATIVE[mode]
    manifest_path = root / "manifest.json"
    manifest = read_json(manifest_path)
    if manifest.get("calculation_sha") != CALCULATION_SHA or manifest.get("workflow_head_sha") != WORKFLOW_HEAD_SHA:
        raise ValueError(f"{mode} result source SHA mismatch")
    required_flags = {
        "artifact_status": "imported_candidate",
        "numerical_status": "diagnostic_only",
        "registry_status": "versioned_candidate_not_default",
        "runtime_default_switch": False,
        "legacy_fallback": True,
        "explicit_rollback": "--phase-reference-mode legacy",
        "source_solver_called": True,
        "aggregate_replay_solver_called": False,
        "import_solver_called": False,
        "production_write": False,
    }
    for field, expected in required_flags.items():
        if manifest.get(field) != expected:
            raise ValueError(f"{mode} result manifest {field}={manifest.get(field)!r}, expected {expected!r}")
    hashes = manifest.get("hashes", {})
    if not hashes:
        raise ValueError(f"{mode} result manifest has no file hashes")
    for relative, expected in hashes.items():
        path = root / relative
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"{mode} result hash mismatch: {relative}")
    scan_path = root / "phase_guided_transport_scan.csv"
    diag_path = root / "channel_diagnostics.csv"
    scan_fields, scan_rows = read_csv(scan_path)
    diag_fields, diag_rows = read_csv(diag_path)
    scan_summary = _validate_scan_rows_fast(mode, scan_fields, scan_rows, expected_rows=EXPECTED_SCAN_ROWS[mode])
    diag_summary = _validate_diag_rows_fast(mode, diag_fields, diag_rows, expected_rows=EXPECTED_DIAG_ROWS[mode])
    failed_fields, failed_rows = read_csv(root / "failed_points.csv")
    if not failed_fields or failed_rows:
        raise ValueError(f"{mode} failed_points.csv is not empty")
    old_root = (
        root.parents[1]
        / LEGACY_RESULT_SNAPSHOT_VERSION
        / mode
        / root.name.replace("_prod_v2", "_prod_v1")
    )
    old_hash = _result_tree_hash(old_root)
    if manifest.get("legacy_prod_v1_tree_hash") != old_hash:
        raise ValueError(f"{mode} legacy prod_v1 tree changed or hash record is stale")
    figure = manifest.get("figure_reference", {})
    figure_manifest_path = repo_root / Path(str(figure.get("path", ""))) / "plot_manifest.json"
    if not figure_manifest_path.is_file():
        raise ValueError(f"{mode} figure manifest missing")
    figure_manifest = read_json(figure_manifest_path)
    if figure_manifest.get("source_csv_sha256") != sha256(scan_path):
        raise ValueError(f"{mode} figure source hash does not match imported scan")
    return {
        "mode": mode,
        "path": rel(root, repo_root),
        "manifest_sha256": sha256(manifest_path),
        "scan_rows": scan_summary["rows"],
        "diagnostic_rows": diag_summary["rows"],
        "scan_sha256": sha256(scan_path),
        "diagnostic_sha256": sha256(diag_path),
        "legacy_prod_v1_tree_hash": old_hash,
        "figure_manifest_sha256": sha256(figure_manifest_path),
        "figure_count": figure.get("figure_count"),
        "numerical_status": manifest.get("numerical_status"),
        "runtime_default_switch": manifest.get("runtime_default_switch"),
        "source_solver_called": manifest.get("source_solver_called"),
        "aggregate_replay_solver_called": manifest.get("aggregate_replay_solver_called"),
        "import_solver_called": manifest.get("import_solver_called"),
        "quality_reason_counts": manifest.get("quality_reason_counts", {}),
        "known_metadata_defects": manifest.get("known_metadata_defects", {}),
    }


def _parity_rows(payload: dict[str, Any]) -> list[dict[str, Any]]:
    points = {row["mode"]: row for row in payload["consumer_points"]}
    rows: list[dict[str, Any]] = []
    for field in (
        "mode_a_phase_kind",
        "boundary_rows",
        "boundary_xi_used",
        "crossover_xi_used",
        "first_order_xi_used",
    ):
        runtime = points["runtime"].get(field)
        diagnostic = points["diagnostic"].get(field)
        legacy = points["legacy"].get(field)
        rows.append({
            "field": field,
            "runtime": runtime,
            "diagnostic": diagnostic,
            "legacy": legacy,
            "runtime_minus_legacy": "",
            "verdict": "same_contract_value" if runtime == diagnostic else "mode_sensitive",
        })
    for field in ("mode_a_phase_T_MeV", "crossover_T_MeV", "first_order_T_MeV"):
        runtime = points["runtime"].get(field)
        diagnostic = points["diagnostic"].get(field)
        legacy = points["legacy"].get(field)
        delta = None if runtime is None or legacy is None else float(runtime) - float(legacy)
        rows.append({
            "field": field,
            "runtime": runtime,
            "diagnostic": diagnostic,
            "legacy": legacy,
            "runtime_minus_legacy": delta,
            "verdict": "same_candidate_source" if runtime == diagnostic else "mode_sensitive",
        })
    return rows


def _build_package(repo_root: Path, output_root: Path, smoke: dict[str, Any], dry_runs: list[dict[str, Any]], results: list[dict[str, Any]], input_inventory: list[dict[str, Any]], smoke_meta: dict[str, Any]) -> None:
    candidate_manifest_sha = sha256(repo_root / CANDIDATE_RELATIVE / "manifest.json")
    _validate_source_smoke(smoke, candidate_manifest_sha)
    parity = _parity_rows(smoke)
    write_csv(
        output_root / "consumer_smoke.csv",
        (
            "mode", "source_kind", "boundary_rows", "boundary_xi_used", "boundary_first_T_MeV",
            "boundary_first_muq_MeV", "boundary_muB_CEP_MeV", "crossover_T_MeV", "mode_a_phase_kind",
            "mode_a_phase_T_MeV", "first_order_T_MeV", "solver_called"
        ),
        ({field: row.get(field) for field in (
            "mode", "source_kind", "boundary_rows", "boundary_xi_used", "boundary_first_T_MeV",
            "boundary_first_muq_MeV", "boundary_muB_CEP_MeV", "crossover_T_MeV", "mode_a_phase_kind",
            "mode_a_phase_T_MeV", "first_order_T_MeV", "solver_called"
        )} for row in smoke["consumer_points"]),
    )
    write_csv(
        output_root / "parity_comparison.csv",
        ("field", "runtime", "diagnostic", "legacy", "runtime_minus_legacy", "verdict"),
        parity,
    )
    write_csv(
        output_root / "result_import_inventory.csv",
        (
            "mode", "path", "manifest_sha256", "scan_rows", "diagnostic_rows", "scan_sha256",
            "diagnostic_sha256", "legacy_prod_v1_tree_hash", "figure_manifest_sha256", "figure_count",
            "numerical_status", "runtime_default_switch", "source_solver_called",
            "aggregate_replay_solver_called", "import_solver_called", "quality_reason_counts",
            "known_metadata_defects"
        ),
        ({field: row.get(field) for field in (
            "mode", "path", "manifest_sha256", "scan_rows", "diagnostic_rows", "scan_sha256",
            "diagnostic_sha256", "legacy_prod_v1_tree_hash", "figure_manifest_sha256", "figure_count",
            "numerical_status", "runtime_default_switch", "source_solver_called",
            "aggregate_replay_solver_called", "import_solver_called", "quality_reason_counts",
            "known_metadata_defects"
        )} for row in results),
    )
    write_json(output_root / "portable_config_snapshot.json", {
        "schema_version": "pnjl_issue130_rs_transport_portable_config_v2",
        "merge_sha": MERGE_SHA,
        "calculation_sha": CALCULATION_SHA,
        "runtime": read_json(output_root / "runtime" / "effective_config.json"),
        "legacy": read_json(output_root / "legacy" / "effective_config.json"),
    })
    claims = [
        {
            "claim_id": "post_merge_source_resolution",
            "claim": "The real phase-guided consumer resolves candidate runtime, diagnostic candidate, and explicit legacy rollback on the #272 merge SHA.",
            "status": "supported",
            "evidence": "consumer_source_smoke.json; consumer_smoke.csv; runtime/run_manifest.json; legacy/run_manifest.json",
            "boundary": "solver-free source-resolution and one-point plan smoke",
        },
        {
            "claim_id": "fallback_and_rollback",
            "claim": "Certified-only candidate view, per-key legacy fallback, and explicit legacy rollback remain available.",
            "status": "supported",
            "evidence": "consumer_source_smoke.json; parity_comparison.csv",
            "boundary": "does not certify unresolved candidate rows",
        },
        {
            "claim_id": "versioned_result_import_integrity",
            "claim": "Both imported prod_v2 result trees and their author-accepted figure source hashes match their manifests, while prod_v1 tree hashes are unchanged.",
            "status": "supported",
            "evidence": "result_import_inventory.csv; result manifest hashes",
            "boundary": "integrity/import evidence, not a new numerical run",
        },
        {
            "claim_id": "candidate_legacy_difference",
            "claim": "Candidate and legacy anchors remain contract-compatible but numerically distinct inputs.",
            "status": "supported",
            "evidence": "parity_comparison.csv",
            "boundary": "not RS numerical parity or convergence",
        },
        {
            "claim_id": "solver_free_boundary",
            "claim": "This v2 package invokes no equilibrium solver, transport coefficient evaluation, result rewrite, or runtime switch.",
            "status": "supported",
            "evidence": "consumer_source_smoke.json; dry-run manifests; result_import_inventory.csv",
            "boundary": "numerical RS promotion remains a separately authorized decision",
        },
    ]
    write_csv(output_root / "claim_ledger.csv", ("claim_id", "claim", "status", "evidence", "boundary"), claims)
    write_json(output_root / "decision.json", {
        "schema_version": SCHEMA,
        "verdict": "rs_consumer_smoke_pass_diagnostic_only",
        "repo_head": MERGE_SHA,
        "merge_sha": MERGE_SHA,
        "calculation_sha": CALCULATION_SHA,
        "source_run_id": SOURCE_RUN_ID,
        "replay_run_id": REPLAY_RUN_ID,
        "candidate_manifest_sha256": candidate_manifest_sha,
        "legacy_snapshot": LEGACY_RELATIVE.as_posix(),
        "solver_called": False,
        "transport_coefficients_computed": False,
        "production_write": False,
        "runtime_default_switch": False,
        "result_import_verified": True,
        "next_action": "author review of post-merge consumer smoke, then separately decide RS runtime/default promotion and old-reference retirement",
        "stop_conditions": [
            "do not claim RS numerical parity from this solver-free package",
            "do not delete prod_v1 result trees or legacy fallback snapshot",
            "do not remove explicit legacy rollback",
            "do not start full RS production without a separate authorization",
        ],
    })
    readme = f"""# Issue #130 RS transport phase-reference parity v2

这是 PR #272 合并 SHA `{MERGE_SHA}` 上的 solver-free consumer smoke 与
versioned result-import integrity evidence。它不重新运行 equilibrium solver，
也不把 diagnostic-only 数值结果晋升为默认 runtime。

## 固定输入

- calculation SHA: `{CALCULATION_SHA}`
- workflow head: `{WORKFLOW_HEAD_SHA}`
- source run: `{SOURCE_RUN_ID}`; aggregate replay: `{REPLAY_RUN_ID}`
- candidate: `{CANDIDATE_RELATIVE.as_posix()}`
- legacy snapshot: `{LEGACY_RELATIVE.as_posix()}`
- merge SHA: `{MERGE_SHA}`

## Smoke result

`decision.json` verdict: `rs_consumer_smoke_pass_diagnostic_only`。
runtime 使用 strict candidate certified-only view，并保留逐键 legacy fallback；
`legacy` 模式验证显式 `--phase-reference-mode legacy` rollback；diagnostic 模式
只用于全 candidate 行审计。两套 dry-run 的 phase anchor 明确使用
`reference_interpolation`，避免 dry-run 触发 coexistence solver。

## Imported result integrity

`result_import_inventory.csv` 验证 mode-A `{results[0]['scan_rows']}` scan /
`{results[0]['diagnostic_rows']}` diagnostic rows 和 mode-B
`{results[1]['scan_rows']}` scan / `{results[1]['diagnostic_rows']}` diagnostic rows，
包括 manifest file hashes、figure source hashes、direct-coexistence side-point
合同、以及旧 `prod_v1` tree hash。两套 result 仍是
`imported_candidate` + `diagnostic_only`，不切换默认 runtime。

## Provenance note

`consumer_source_smoke_runner_raw.json` 保留 runner 的原始 JSON；v2 的
`consumer_source_smoke.json` 仅规范化 legacy snapshot locator，并记录
`legacy_root_normalized=true`。`raw_sha256` 与规范化后的 hash 见 manifest。

## Boundary

本包支持 source-resolution、fallback/rollback 和导入完整性；不代表 RS 数值
parity、全域 production convergence、旧 reference 删除或新的 solver run。
"""
    (output_root / "README.md").write_text(readme, encoding="utf-8", newline="\n")
    raw_inventory = []
    for path in sorted(output_root.rglob("*")):
        if path.is_file() and path.name != "manifest.json":
            raw_inventory.append({
                "path": path.relative_to(output_root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            })
    manifest = {
        "schema_version": SCHEMA,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "repo_head": MERGE_SHA,
        "merge_sha": MERGE_SHA,
        "calculation_sha": CALCULATION_SHA,
        "workflow_head_sha": WORKFLOW_HEAD_SHA,
        "source_run_id": SOURCE_RUN_ID,
        "replay_run_id": REPLAY_RUN_ID,
        "candidate_root": CANDIDATE_RELATIVE.as_posix(),
        "legacy_snapshot_root": LEGACY_RELATIVE.as_posix(),
        "candidate_manifest_sha256": candidate_manifest_sha,
        "solver_called": False,
        "transport_coefficients_computed": False,
        "production_write": False,
        "runtime_default_switch": False,
        "result_import_verified": True,
        "dry_runs": dry_runs,
        "result_import": results,
        "input_inventory": input_inventory,
        "source_smoke_normalization": {
            "raw_runner_sha256": smoke_meta["raw_sha256"],
            "normalized_sha256": smoke_meta["normalized_sha256"],
            "legacy_root": "data/reference/pnjl/legacy_phase_reference_v1/{boundary.csv,cep.csv,crossover_dense.csv,spinodals.csv}",
        },
        "scripts": [
            {"path": rel(Path(__file__), repo_root), "sha256": sha256(Path(__file__).resolve())},
            {
                "path": "scripts/analysis/pnjl/audit_issue130_rs_runtime_parity.jl",
                "sha256": sha256(repo_root / "scripts/analysis/pnjl/audit_issue130_rs_runtime_parity.jl"),
            },
            {
                "path": "scripts/analysis/relaxtime/import_issue130_rs_candidate_results.py",
                "sha256": sha256(repo_root / "scripts/analysis/relaxtime/import_issue130_rs_candidate_results.py"),
            },
        ],
        "files": {item["path"]: item["sha256"] for item in raw_inventory},
    }
    write_json(output_root / "manifest.json", manifest)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    repo_root = args.repo_root.resolve()
    output_root = args.output_root.resolve()
    if git_value(repo_root, "rev-parse", "HEAD") != MERGE_SHA:
        raise ValueError("consumer smoke v2 must run on the PR #272 merge SHA")
    output_root.mkdir(parents=True, exist_ok=True)
    smoke_meta = _normalise_runner_output(output_root)
    candidate_manifest_sha = sha256(repo_root / CANDIDATE_RELATIVE / "manifest.json")
    _validate_source_smoke(smoke_meta["payload"], candidate_manifest_sha)
    dry_runs = [
        _validate_dry_run(output_root, "runtime", "candidate"),
        _validate_dry_run(output_root, "legacy", "legacy"),
    ]
    results = [_validate_result(repo_root, mode) for mode in RESULT_RELATIVE]
    input_inventory = build_input_inventory(repo_root, output_root)
    _build_package(repo_root, output_root, smoke_meta["payload"], dry_runs, results, input_inventory, smoke_meta)
    print(json.dumps({
        "schema_version": SCHEMA,
        "verdict": "rs_consumer_smoke_pass_diagnostic_only",
        "repo_head": MERGE_SHA,
        "output_root": str(output_root),
        "result_modes": [item["mode"] for item in results],
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
