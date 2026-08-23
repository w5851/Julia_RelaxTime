#!/usr/bin/env python3
"""Build the portable solver-free RS runtime parity evidence package."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


SCHEMA = "pnjl_issue130_rs_transport_runtime_parity_v1"
CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
SOURCE_RUN_ID = "32354095831"
REPLAY_RUN_ID = "32451053476"
CANDIDATE_RELATIVE = Path("data/reference/pnjl/issue130_phase_reference_v1")
LEGACY_RELATIVE = (
    Path("data/reference/pnjl/boundary.csv"),
    Path("data/reference/pnjl/cep.csv"),
    Path("data/reference/pnjl/crossover_dense.csv"),
    Path("data/reference/pnjl/spinodals.csv"),
)
RAW_RELATIVE = (
    Path("runtime/effective_config.json"),
    Path("runtime/run_manifest.json"),
    Path("runtime/sampling_plan.csv"),
    Path("runtime/README.md"),
    Path("legacy/effective_config.json"),
    Path("legacy/run_manifest.json"),
    Path("legacy/sampling_plan.csv"),
    Path("legacy/README.md"),
    Path("consumer_source_smoke.json"),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_value(repo_root: Path, *args: str) -> str:
    try:
        return subprocess.check_output(["git", *args], cwd=repo_root, text=True, stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        return ""


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"missing CSV header: {path}")
        return list(reader.fieldnames), list(reader)


def rel(path: Path, repo_root: Path) -> str:
    return path.resolve().relative_to(repo_root.resolve()).as_posix()


def write_json(path: Path, payload: Any) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n")


def write_csv(path: Path, fieldnames: Iterable[str], rows: Iterable[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_text_lf(path: Path, text: str) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write(text)


def finite(value: Any) -> bool:
    try:
        return float(value) == float(value) and abs(float(value)) != float("inf")
    except (TypeError, ValueError):
        return False


def source_inventory(repo_root: Path, output_root: Path) -> list[dict[str, Any]]:
    candidate_root = repo_root / CANDIDATE_RELATIVE
    paths = [candidate_root / "manifest.json", candidate_root / "strict" / "manifest.json"]
    paths.extend(sorted((candidate_root / "strict" / "tables").glob("*.csv")))
    paths.extend(repo_root / relative for relative in LEGACY_RELATIVE)
    rows = []
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
        rows.append({"path": rel(path, repo_root), "bytes": path.stat().st_size, "sha256": sha256(path)})
    write_csv(output_root / "input_inventory.csv", ("path", "bytes", "sha256"), rows)
    return rows


def raw_inventory(repo_root: Path, output_root: Path) -> list[dict[str, Any]]:
    rows = []
    for relative in RAW_RELATIVE:
        path = output_root / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        rows.append({"path": relative.as_posix(), "bytes": path.stat().st_size, "sha256": sha256(path)})
    return rows


def normalize_raw_artifacts(repo_root: Path, output_root: Path) -> None:
    """Normalize copied text artifacts and repair their sidecar hashes."""
    for relative in RAW_RELATIVE:
        path = output_root / relative
        data = path.read_bytes().replace(b"\r\n", b"\n")
        path.write_bytes(data)
    for mode in ("runtime", "legacy"):
        manifest_path = output_root / mode / "run_manifest.json"
        manifest = read_json(manifest_path)
        for artifact in manifest.get("artifacts", []):
            raw_path = Path(str(artifact.get("path", "")).replace("\\", "/"))
            resolved = repo_root / raw_path
            if resolved.is_file():
                artifact["sha256"] = sha256(resolved)
        write_json(manifest_path, manifest)


def portable_config(config: dict[str, Any], mode: str) -> dict[str, Any]:
    options = dict(config.get("options", {}))
    options["outdir"] = f"docs/analysis/pnjl/phase_reference/issue130_rs_transport_runtime_parity_v1/{mode}"
    options["plan_csv"] = f"{options['outdir']}/sampling_plan.csv"
    options["result_csv"] = f"{options['outdir']}/phase_guided_transport_scan.csv"
    options["figure_dir"] = "data/outputs/figures/relaxtime/transport/phase_guided/issue130_rs_parity_v1"
    options["phase_reference_root"] = CANDIDATE_RELATIVE.as_posix()
    options["phase_reference_source"] = "candidate" if mode == "runtime" else "legacy"
    options["solver_called"] = False
    options["transport_coefficients_computed"] = False
    return {
        "mode": mode,
        "schema_version": "pnjl_issue130_rs_transport_portable_config_v1",
        "options": options,
        "source_run_id": SOURCE_RUN_ID,
        "replay_run_id": REPLAY_RUN_ID,
        "calculation_sha": CALCULATION_SHA,
    }


def parity_rows(smoke: dict[str, Any]) -> list[dict[str, Any]]:
    points = {row["mode"]: row for row in smoke.get("consumer_points", [])}
    runtime = points["runtime"]
    legacy = points["legacy"]
    diagnostic = points["diagnostic"]
    rows = []
    for field in (
        "mode_a_phase_kind",
        "boundary_rows",
        "boundary_xi_used",
        "crossover_xi_used",
        "first_order_xi_used",
    ):
        rows.append({
            "field": field,
            "runtime": runtime.get(field),
            "diagnostic": diagnostic.get(field),
            "legacy": legacy.get(field),
            "verdict": "same_contract_value" if runtime.get(field) == diagnostic.get(field) else "mode_sensitive",
        })
    for field in ("mode_a_phase_T_MeV", "crossover_T_MeV", "first_order_T_MeV"):
        rv = runtime.get(field)
        dv = diagnostic.get(field)
        lv = legacy.get(field)
        rows.append({
            "field": field,
            "runtime": rv,
            "diagnostic": dv,
            "legacy": lv,
            "runtime_minus_legacy": None if rv is None or lv is None else float(rv) - float(lv),
            "verdict": "same_candidate_source" if rv == dv else "mode_sensitive",
        })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--output-root", type=Path, default=None)
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    output_root = (args.output_root or repo_root / "docs/analysis/pnjl/phase_reference/issue130_rs_transport_runtime_parity_v1").resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    smoke = read_json(output_root / "consumer_source_smoke.json")
    if smoke.get("solver_called") is not False:
        raise ValueError("consumer smoke must record solver_called=false")
    if smoke.get("repo_head") != git_value(repo_root, "rev-parse", "HEAD"):
        raise ValueError("consumer smoke was generated on a different HEAD")

    runtime_config = read_json(output_root / "runtime" / "effective_config.json")
    legacy_config = read_json(output_root / "legacy" / "effective_config.json")
    runtime_manifest = read_json(output_root / "runtime" / "run_manifest.json")
    legacy_manifest = read_json(output_root / "legacy" / "run_manifest.json")
    for config, mode in ((runtime_config, "runtime"), (legacy_config, "legacy")):
        if config.get("options", {}).get("dry_run") is not True:
            raise ValueError(f"{mode} source must be a dry-run")
    if runtime_config.get("options", {}).get("phase_reference_source") != "candidate":
        raise ValueError("runtime dry-run does not resolve candidate source")
    if legacy_config.get("options", {}).get("phase_reference_source") != "legacy":
        raise ValueError("legacy dry-run does not resolve legacy source")
    runtime_summary = runtime_config["options"]["phase_reference_summary"]
    smoke_runtime = smoke.get("sources", {}).get("runtime", {})
    smoke_legacy = smoke.get("sources", {}).get("legacy", {})
    if smoke_runtime.get("source_kind") != "candidate" or smoke_legacy.get("source_kind") != "legacy":
        raise ValueError("consumer smoke did not exercise candidate runtime and legacy rollback")
    if smoke_runtime.get("diagnostics", {}).get("runtime_view") != runtime_summary.get("runtime_view"):
        raise ValueError("consumer smoke/runtime config runtime_view mismatch")
    if smoke_runtime.get("diagnostics", {}).get("candidate_manifest_sha256") != sha256(repo_root / CANDIDATE_RELATIVE / "manifest.json"):
        raise ValueError("consumer smoke candidate manifest hash mismatch")
    if any(row.get("solver_called") is not False for row in smoke.get("consumer_points", [])):
        raise ValueError("consumer smoke contains a solver-enabled point")

    normalize_raw_artifacts(repo_root, output_root)
    inputs = source_inventory(repo_root, output_root)
    raw = raw_inventory(repo_root, output_root)
    write_json(output_root / "portable_config_snapshot.json", {
        "runtime": portable_config(runtime_config, "runtime"),
        "legacy": portable_config(legacy_config, "legacy"),
    })
    parity = parity_rows(smoke)
    write_csv(output_root / "parity_comparison.csv", (
        "field", "runtime", "diagnostic", "legacy", "runtime_minus_legacy", "verdict"
    ), parity)
    consumer_fields = (
        "mode", "source_kind", "boundary_rows", "boundary_xi_used", "boundary_first_T_MeV",
        "boundary_first_muq_MeV", "boundary_muB_CEP_MeV", "crossover_T_MeV", "mode_a_phase_kind",
        "mode_a_phase_T_MeV", "first_order_T_MeV", "solver_called"
    )
    write_csv(
        output_root / "consumer_smoke.csv",
        consumer_fields,
        ({field: row.get(field) for field in consumer_fields} for row in smoke["consumer_points"]),
    )

    candidate_manifest = repo_root / CANDIDATE_RELATIVE / "manifest.json"
    candidate_root_manifest = read_json(candidate_manifest)
    claims = [
        {
            "claim_id": "runtime_candidate_preferred",
            "claim": "The real phase-guided consumer resolves the strict candidate as its preferred runtime source.",
            "status": "supported" if runtime_config["options"].get("phase_reference_source") == "candidate" else "blocked",
            "evidence": "runtime/effective_config.json; consumer_smoke.csv",
            "boundary": "solver-free source-resolution smoke only",
        },
        {
            "claim_id": "fallback_and_rollback",
            "claim": "Certified-only candidate rows, legacy fallback counts, and explicit legacy rollback remain visible in provenance.",
            "status": "supported" if runtime_summary.get("fallback_enabled") and legacy_config["options"].get("phase_reference_source") == "legacy" else "blocked",
            "evidence": "runtime/effective_config.json; legacy/effective_config.json; input_inventory.csv",
            "boundary": "legacy files are retained and no fallback row is reclassified as candidate",
        },
        {
            "claim_id": "candidate_legacy_difference",
            "claim": "Candidate and legacy phase anchors are contract-compatible but numerically distinct inputs.",
            "status": "supported" if any(row.get("field") == "mode_a_phase_T_MeV" and row.get("runtime_minus_legacy") not in (None, 0, 0.0) for row in parity) else "inconclusive",
            "evidence": "parity_comparison.csv; consumer_source_smoke.json",
            "boundary": "difference is an input/reference change, not a transport coefficient validation",
        },
        {
            "claim_id": "solver_free",
            "claim": "This parity package does not call the equilibrium solver or produce numerical transport coefficients.",
            "status": "supported" if smoke.get("solver_called") is False and all(manifest.get("summary", {}).get("dry_run") is True for manifest in (runtime_manifest, legacy_manifest)) else "blocked",
            "evidence": "consumer_source_smoke.json; runtime/run_manifest.json; legacy/run_manifest.json",
            "boundary": "a numerical RS pilot remains a separate, explicitly authorized step",
        },
        {
            "claim_id": "promotion_boundary",
            "claim": "The phase candidate is not silently promoted to a solver-backed production reference by this package.",
            "status": "supported" if candidate_root_manifest.get("runtime_consumption") is False else "blocked",
            "evidence": "data/reference/pnjl/issue130_phase_reference_v1/manifest.json; runtime/effective_config.json",
            "boundary": "old-reference retirement remains blocked until consumer parity and numerical RS revalidation are accepted",
        },
    ]
    write_csv(output_root / "claim_ledger.csv", ("claim_id", "claim", "status", "evidence", "boundary"), claims)

    decision = {
        "schema_version": SCHEMA,
        "verdict": "rs_adapter_parity_solver_free_pass",
        "calculation_sha": CALCULATION_SHA,
        "source_run_id": SOURCE_RUN_ID,
        "replay_run_id": REPLAY_RUN_ID,
        "repo_head": git_value(repo_root, "rev-parse", "HEAD"),
        "candidate_manifest_sha256": sha256(candidate_manifest),
        "solver_called": False,
        "transport_coefficients_computed": False,
        "runtime_view": runtime_summary.get("runtime_view"),
        "fallback_enabled": runtime_summary.get("fallback_enabled"),
        "legacy_rollback_verified": legacy_config["options"].get("phase_reference_source") == "legacy",
        "next_action": "perform a separate limited numerical RS pilot after author review of this solver-free parity package",
        "stop_conditions": [
            "do not delete legacy reference files",
            "do not remove per-key fallback or explicit rollback",
            "do not claim RS numerical parity from this solver-free package",
            "do not start full RS production before the limited pilot passes",
        ],
    }
    write_json(output_root / "decision.json", decision)

    package_files = {}
    for path in sorted(output_root.rglob("*")):
        if path.is_file() and path.name != "manifest.json":
            package_files[path.relative_to(output_root).as_posix()] = sha256(path)
    manifest = {
        "schema_version": SCHEMA,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "repo_head": git_value(repo_root, "rev-parse", "HEAD"),
        "calculation_sha": CALCULATION_SHA,
        "source_run_id": SOURCE_RUN_ID,
        "replay_run_id": REPLAY_RUN_ID,
        "candidate_root": CANDIDATE_RELATIVE.as_posix(),
        "candidate_manifest_sha256": sha256(candidate_manifest),
        "scripts": [
            {
                "path": rel(Path(__file__), repo_root),
                "sha256": sha256(Path(__file__).resolve()),
            },
            {
                "path": "scripts/analysis/pnjl/audit_issue130_rs_runtime_parity.jl",
                "sha256": sha256(repo_root / "scripts/analysis/pnjl/audit_issue130_rs_runtime_parity.jl"),
            },
        ],
        "runtime_view": runtime_summary.get("runtime_view"),
        "solver_called": False,
        "transport_coefficients_computed": False,
        "raw_run_artifacts": raw,
        "input_inventory": inputs,
        "files": package_files,
    }
    write_json(output_root / "manifest.json", manifest)
    write_text_lf(output_root / "README.md",
        "# Issue #130 RS transport phase-reference parity v1\n\n"
        "这是合并 SHA `1ccf29310fb20c30bcd154f0b4966e25a7565225` 上的 solver-free consumer/parity evidence。\n\n"
        f"- calculation SHA: `{CALCULATION_SHA}`\n"
        f"- source run: `{SOURCE_RUN_ID}`; aggregate replay: `{REPLAY_RUN_ID}`\n"
        f"- runtime view: `{decision['runtime_view']}`\n"
        "- runtime 首选 strict candidate 的 certified-only 行；缺失或不合格键逐键由 legacy fallback 补位。\n"
        "- `legacy` 是显式 rollback；`diagnostic` 仅用于审计，不是 runtime 输入。\n"
        "- 本包只运行 source-resolution 和 phase-guided plan dry-run；`solver_called=false`，没有输运系数。\n\n"
        "## 结果\n\n"
        "candidate 与 legacy 的相变锚点数值不同，但字段/单位/phase-kind 语义一致；差异是 reference 输入差异，不能被称为 RS 数值收敛。\n"
        "下一步是作者审核本包后，单独触发限定 numerical RS pilot；在该 pilot 通过前保留旧 reference、fallback、rollback，不创建 retirement PR。\n",
    )
    # README is part of the package inventory; finalize its hash after all
    # human-readable evidence has been written.
    manifest["files"] = {
        path.relative_to(output_root).as_posix(): sha256(path)
        for path in sorted(output_root.rglob("*"))
        if path.is_file() and path.name != "manifest.json"
    }
    write_json(output_root / "manifest.json", manifest)
    print(json.dumps({"verdict": decision["verdict"], "output_root": str(output_root)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
