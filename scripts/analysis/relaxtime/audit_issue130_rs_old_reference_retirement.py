#!/usr/bin/env python3
"""Build a solver-free retirement audit for the Issue #130 RS references.

The audit inventories the versioned RS transport result/figure trees, scans
tracked repository text for old ``prod_v1`` references, and verifies the
byte-preserving fallback relation recorded by the imported ``prod_v2``
manifests.  It deliberately does not run Julia, a transport solver, or a
workflow, and it never mutates or removes a result tree.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


SCHEMA_VERSION = "issue130_rs_old_reference_retirement_audit_v1"
CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
WORKFLOW_HEAD_SHA = "22874505877491754eed27519ad8a7b871c82571"
CURRENT_CASE = "first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2"
LEGACY_CASE = "first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1"
HISTORICAL_CASES = (
    "first_canonical_v1_p128_xi001_validated_anchored_prod_v1",
    "first_canonical_v1_p128_xi005_validated_anchored_prod_v1",
)
OLD_CASE_TOKENS = (
    "first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1",
    "first_canonical_v1_p128_xi001_validated_anchored_prod_v1",
    "first_canonical_v1_p128_xi005_validated_anchored_prod_v1",
)
MODES = ("mode_a_fixed_muB_phase_scaled", "mode_b_fixed_T_sparse_muB")
RESULT_ROOT = Path("data/outputs/results/relaxtime/transport/phase_guided")
FIGURE_ROOT = Path("data/outputs/figures/relaxtime/transport/phase_guided")
LEGACY_SNAPSHOT_VERSION = "legacy_prod_v1_snapshot_v1"
LEGACY_RESULT_ROOT = RESULT_ROOT / LEGACY_SNAPSHOT_VERSION
LEGACY_FIGURE_ROOT = FIGURE_ROOT / LEGACY_SNAPSHOT_VERSION
REGISTRY_PATH = RESULT_ROOT / "production_registry.json"
AUDIT_PACKAGE_ROOT = Path("docs/analysis/relaxtime/issue130_rs_old_reference_retirement_audit_v1")
SCAN_NAME = "phase_guided_transport_scan.csv"
DIAGNOSTIC_NAME = "channel_diagnostics.csv"
TEXT_SUFFIXES = {
    ".csv",
    ".jl",
    ".json",
    ".md",
    ".ps1",
    ".py",
    ".rst",
    ".sh",
    ".toml",
    ".txt",
    ".yml",
    ".yaml",
}
HISTORICAL_TOOLING = {
    "scripts/analysis/pnjl/build_figure4_phase_diagram_assets.py",
    "scripts/analysis/pnjl/build_issue130_rs_runtime_consumer_smoke_v2.py",
    "scripts/analysis/relaxtime/audit_issue130_rs_sharded_provenance.py",
    "scripts/analysis/relaxtime/build_phase_guided_publication_clean_v1.py",
    "scripts/analysis/relaxtime/build_phase_guided_pole_sensitive_rendering.py",
    "scripts/analysis/relaxtime/render_issue130_rs_candidate_legacy_figures.py",
    "scripts/plotting/build_figure_cleanup_preflight.py",
    "scripts/plotting/build_figure_retirement_manifest.py",
    "scripts/plotting/render_plotting_pilot.py",
}
ACTIVE_CONSUMER_FILES = {
    "scripts/analysis/relaxtime/build_phase_guided_transport_xi001_jump_analysis.py",
    "scripts/analysis/relaxtime/phase_guided_p128_mechanism_scan.jl",
    "scripts/relaxtime/run_phase_guided_transport_scan.jl",
    "scripts/relaxtime/phase_reference_adapter.jl",
}
NUMERIC_FIELDS = {
    "T_MeV",
    "muq_MeV",
    "muB_MeV",
    "xi",
    "alpha_T",
    "T_fm",
    "muq_fm",
    "residual_norm",
    "Phi",
    "Phibar",
    "m_u",
    "m_d",
    "m_s",
    "rho_baryon",
    "rho_norm",
    "omega_fm4inv",
    "P_fm4inv",
    "epsilon_fm4inv",
    "omega_MeV_fm3",
    "P_MeV_fm3",
    "epsilon_MeV_fm3",
    "n_u",
    "n_d",
    "n_s",
    "n_ubar",
    "n_dbar",
    "n_sbar",
    "tau_u",
    "tau_d",
    "tau_s",
    "tau_ubar",
    "tau_dbar",
    "tau_sbar",
    "tauinv_u",
    "tauinv_d",
    "tauinv_s",
    "tauinv_ubar",
    "tauinv_dbar",
    "tauinv_sbar",
    "eta",
    "sigma",
    "zeta",
    "eta_over_s",
    "zeta_over_s",
    "sigma_over_T",
    "quality_metric",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tree_hash(root: Path) -> str:
    """Use the same path+file-hash contract as the result importer."""

    digest = hashlib.sha256()
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        digest.update(path.relative_to(root).as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(sha256(path).encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    lines = [
        line
        for line in path.read_text(encoding="utf-8-sig").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not lines:
        return [], []
    reader = csv.DictReader(lines)
    if reader.fieldnames is None:
        raise ValueError(f"missing CSV header: {path}")
    return list(reader.fieldnames), [dict(row) for row in reader]


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def write_csv(path: Path, fieldnames: Iterable[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    names = list(fieldnames)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in names})


def git_head(repo_root: Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_root,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return ""


def tracked_files(repo_root: Path) -> list[Path]:
    try:
        raw = subprocess.check_output(
            ["git", "ls-files", "-z"],
            cwd=repo_root,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise RuntimeError("git ls-files failed") from exc
    return [repo_root / item.decode("utf-8") for item in raw.split(b"\0") if item]


def relative(path: Path, repo_root: Path) -> str:
    return path.resolve().relative_to(repo_root.resolve()).as_posix()


def result_case_root(repo_root: Path, mode: str, case_slug: str) -> Path:
    """Resolve a result case without treating the retired path as canonical."""

    base = LEGACY_RESULT_ROOT if case_slug == LEGACY_CASE else RESULT_ROOT
    return repo_root / base / mode / case_slug


def figure_case_root(repo_root: Path, mode: str, case_slug: str) -> Path:
    """Resolve the figure tree paired with a result case."""

    base = LEGACY_FIGURE_ROOT if case_slug == LEGACY_CASE else FIGURE_ROOT
    return repo_root / base / mode / case_slug


def classify_case(case_slug: str, layer: str) -> tuple[str, str, str]:
    if case_slug == CURRENT_CASE:
        return "current_prod_v2", "canonical_current", "retain_canonical"
    if case_slug == LEGACY_CASE:
        if layer == "result":
            return "legacy_prod_v1", "legacy_fallback_required", "retain_immutable_legacy_snapshot"
        return "legacy_prod_v1", "legacy_figure_history", "retain_immutable_legacy_figure_snapshot"
    if case_slug in HISTORICAL_CASES or case_slug.endswith("_prod_v1_convergence"):
        return "historical_prod_v1", "historical_evidence", "retain_historical_evidence"
    if "_prod_v1" in case_slug:
        return "unknown_prod_v1", "manual_review", "hold_for_manual_review"
    return "other", "out_of_scope", "do_not_touch"


def registry_index(registry: dict[str, Any]) -> dict[tuple[str, str], dict[str, Any]]:
    return {
        (str(entry.get("case_slug", "")), str(entry.get("mode", ""))): entry
        for entry in registry.get("entries", [])
        if isinstance(entry, dict)
    }


def manifest_hash_mismatches(root: Path, manifest: dict[str, Any]) -> int:
    return len(manifest_hash_mismatch_paths(root, manifest))


def manifest_hash_mismatch_paths(root: Path, manifest: dict[str, Any]) -> list[str]:
    mismatches = []
    for relative_path, expected in (manifest.get("hashes") or {}).items():
        path = root / str(relative_path)
        if not path.is_file() or sha256(path) != str(expected):
            mismatches.append(str(relative_path))
    return mismatches


def csv_stats(path: Path, mode: str, *, diagnostic: bool = False) -> dict[str, Any]:
    if not path.is_file():
        return {
            "path_exists": False,
            "rows": 0,
            "fields": 0,
            "duplicate_keys": 0,
            "nonfinite_values": 0,
            "nonconverged_rows": 0,
        }
    fields, rows = read_csv(path)
    if diagnostic:
        key_fields = ("T_MeV", "muB_MeV", "xi", "species", "channel")
    elif mode.startswith("mode_a"):
        key_fields = ("muB_MeV", "xi", "alpha_T")
    else:
        key_fields = ("T_MeV", "muB_MeV", "xi")
    keys = [tuple(row.get(field, "") for field in key_fields) for row in rows]
    duplicate_keys = sum(count - 1 for count in Counter(keys).values() if count > 1)
    nonfinite = 0
    if not diagnostic:
        for row in rows:
            for field in NUMERIC_FIELDS.intersection(fields):
                value = row.get(field, "").strip()
                if not value:
                    continue
                try:
                    number = float(value)
                except ValueError:
                    continue
                if mode.startswith("mode_b") and field == "alpha_T":
                    # Mode-B has no thermal scale factor; the producer records
                    # this column as NaN by contract.
                    continue
                if not (number == number and abs(number) != float("inf")):
                    nonfinite += 1
    nonconverged = sum(
        row.get("converged", "").strip().lower() not in {"", "true"}
        for row in rows
    )
    return {
        "path_exists": True,
        "rows": len(rows),
        "fields": len(fields),
        "duplicate_keys": duplicate_keys,
        "nonfinite_values": nonfinite,
        "nonconverged_rows": nonconverged,
    }


def expected_scan_rows(manifest: dict[str, Any]) -> int | str:
    gates = manifest.get("gates")
    if isinstance(gates, dict) and isinstance(gates.get("scan_rows"), int):
        return gates["scan_rows"]
    summary = manifest.get("scan_summary")
    if isinstance(summary, dict) and isinstance(summary.get("rows"), int):
        return summary["rows"]
    return ""


def scan_contract_ok(scan: dict[str, Any], expected_rows: int | str) -> bool:
    """Return whether a production scan has a minimally usable shape.

    Historical convergence trees intentionally have no production scan CSV and
    therefore are not evaluated by this contract.  Current and legacy result
    trees must have a non-empty scan, match the manifest row count when one is
    recorded, and contain no duplicate/non-finite values.
    """

    if not scan["path_exists"] or scan["rows"] <= 0:
        return False
    if isinstance(expected_rows, int) and scan["rows"] != expected_rows:
        return False
    return (
        scan["duplicate_keys"] == 0
        and scan["nonfinite_values"] == 0
        and scan["nonconverged_rows"] == 0
    )


def status_consistency(case_slug: str, registry_status: str, manifest_status: str) -> str:
    if case_slug == CURRENT_CASE and registry_status == "approved" and manifest_status == "versioned_candidate_not_default":
        return "semantic_split_expected"
    if case_slug == LEGACY_CASE and registry_status == "superseded_for_manuscript" and manifest_status == "current_candidate":
        return "historical_manifest_status_stale_but_registry_superseded"
    if case_slug in HISTORICAL_CASES and registry_status == "superseded_for_manuscript" and not manifest_status:
        return "historical_manifest_status_not_recorded"
    if registry_status and manifest_status and registry_status == manifest_status:
        return "match"
    if not registry_status:
        return "not_registered"
    return "metadata_review_required"


def scan_tree(repo_root: Path, root: Path, layer: str, mode: str, registry: dict[tuple[str, str], dict[str, Any]]) -> dict[str, Any]:
    manifest_path = root / "manifest.json"
    manifest = read_json(manifest_path) if manifest_path.is_file() else {}
    case_slug = str(manifest.get("case_slug") or root.name)
    family, retention_class, action = classify_case(case_slug, layer)
    registry_entry = registry.get((case_slug, mode), {})
    data_root = root
    data_manifest = manifest
    if layer == "figure":
        data_root = result_case_root(repo_root, mode, case_slug)
        data_manifest_path = data_root / "manifest.json"
        data_manifest = read_json(data_manifest_path) if data_manifest_path.is_file() else {}
    scan = csv_stats(data_root / SCAN_NAME, mode)
    diagnostic = csv_stats(data_root / DIAGNOSTIC_NAME, mode, diagnostic=True)
    expected_rows = expected_scan_rows(data_manifest)
    scan_ok = scan_contract_ok(scan, expected_rows)
    figure_root = figure_case_root(repo_root, mode, case_slug)
    figure_manifest_path = figure_root / "plot_manifest.json"
    figure_manifest = read_json(figure_manifest_path) if figure_manifest_path.is_file() else {}
    scan_path = data_root / SCAN_NAME
    figure_source_hash = str(figure_manifest.get("source_csv_sha256", ""))
    figure_source_match = bool(figure_source_hash and scan_path.is_file() and figure_source_hash == sha256(scan_path))
    files = [path for path in root.rglob("*") if path.is_file()]
    return {
        "layer": layer,
        "mode": mode,
        "case_slug": case_slug,
        "path": relative(root, repo_root),
        "family": family,
        "retention_class": retention_class,
        "recommended_action": action,
        "file_count": len(files),
        "bytes": sum(path.stat().st_size for path in files),
        "tree_sha256": tree_hash(root),
        "manifest_sha256": sha256(manifest_path) if manifest_path.is_file() else "",
        "manifest_hash_mismatch_count": manifest_hash_mismatches(root, manifest),
        "manifest_hash_mismatch_files": ",".join(manifest_hash_mismatch_paths(root, manifest)),
        "registry_status": registry_entry.get("status", ""),
        "manifest_registry_status": manifest.get("registry_status", ""),
        "status_consistency": (
            "figure_manifest_status_not_applicable"
            if layer == "figure"
            else status_consistency(case_slug, str(registry_entry.get("status", "")), str(manifest.get("registry_status", "")))
        ),
        "source_commit": data_manifest.get("source_commit", registry_entry.get("source_commit", "")),
        "calculation_sha": data_manifest.get("calculation_sha", registry_entry.get("calculation_sha", "")),
        "workflow_head_sha": data_manifest.get("workflow_head_sha", registry_entry.get("workflow_head_sha", "")),
        "verdict": data_manifest.get("verdict", ""),
        "numerical_status": data_manifest.get("numerical_status", ""),
        "manuscript_eligible": data_manifest.get("manuscript_eligible", registry_entry.get("manuscript_eligible", "")),
        "legacy_fallback": data_manifest.get("legacy_fallback", ""),
        "explicit_rollback": data_manifest.get("explicit_rollback", ""),
        "production_write": data_manifest.get("production_write", ""),
        "scan_rows": scan["rows"],
        "expected_scan_rows": expected_rows,
        "scan_path_exists": scan["path_exists"],
        "scan_contract_ok": scan_ok,
        "scan_duplicate_keys": scan["duplicate_keys"],
        "scan_nonfinite_values": scan["nonfinite_values"],
        "scan_nonconverged_rows": scan["nonconverged_rows"],
        "diagnostic_rows": diagnostic["rows"],
        "diagnostic_duplicate_keys": diagnostic["duplicate_keys"],
        "figure_manifest_path": relative(figure_manifest_path, repo_root) if figure_manifest_path.is_file() else "",
        "figure_manifest_sha256": sha256(figure_manifest_path) if figure_manifest_path.is_file() else "",
        "figure_source_csv_sha256": figure_source_hash,
        "figure_source_match": figure_source_match,
        "tree_contract_ok": bool(
            (
                (layer == "figure" and figure_manifest_path.is_file())
                or (layer != "figure" and manifest_path.is_file() and not manifest_hash_mismatches(root, manifest))
            )
            and scan_ok
            and diagnostic["duplicate_keys"] == 0
            and figure_source_match
        ),
    }


def collect_trees(repo_root: Path, registry: dict[tuple[str, str], dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    result_root = repo_root / RESULT_ROOT
    for root in sorted(result_root.iterdir() if result_root.is_dir() else []):
        if root.is_dir() and root.name.endswith("_prod_v1_convergence"):
            rows.append(scan_tree(repo_root, root, "result_convergence", "shared_convergence", registry))
    for result_base, figure_base in ((RESULT_ROOT, FIGURE_ROOT), (LEGACY_RESULT_ROOT, LEGACY_FIGURE_ROOT)):
        for mode in MODES:
            result_mode_root = repo_root / result_base / mode
            figure_mode_root = repo_root / figure_base / mode
            for root in sorted(result_mode_root.iterdir() if result_mode_root.is_dir() else []):
                if root.is_dir() and ("_prod_v1" in root.name or "_prod_v2" in root.name):
                    rows.append(scan_tree(repo_root, root, "result", mode, registry))
            for root in sorted(figure_mode_root.iterdir() if figure_mode_root.is_dir() else []):
                if root.is_dir() and ("_prod_v1" in root.name or "_prod_v2" in root.name):
                    rows.append(scan_tree(repo_root, root, "figure", mode, registry))
    return rows


def occurrence_class(path: str) -> str:
    if path in ACTIVE_CONSUMER_FILES:
        return "active_consumer_contract"
    if path in HISTORICAL_TOOLING or "/audit" in path or "cleanup" in path or "publication_clean" in path:
        return "historical_or_audit_tooling"
    if path.startswith("data/outputs/results/") or path.startswith("data/outputs/figures/"):
        return "generated_artifact_or_manifest"
    if path.startswith("docs/analysis/") or path.startswith("docs/dev/archived/"):
        return "historical_evidence"
    if path.startswith("tests/"):
        return "governance_test"
    if path.startswith("docs/dev/active/"):
        return "active_governance_document"
    if path.startswith(".github/") or path.startswith("scripts/"):
        return "code_or_workflow_review"
    return "documentation_contract"


def occurrence_role(path: str, line: str, classification: str) -> str:
    """Add a review-facing role without treating every old token as active."""

    if classification == "active_consumer_contract":
        return "active_consumer"
    context = f"{path} {line}".lower()
    if any(word in context for word in ("fallback", "rollback", "legacy")):
        return "fallback_or_rollback"
    if classification in {"historical_or_audit_tooling", "historical_evidence"}:
        return "historical"
    if classification == "generated_artifact_or_manifest":
        return "generated_artifact"
    if classification == "active_governance_document":
        return "governance"
    if classification == "governance_test":
        return "governance_test"
    return "review_required"


def scan_references(repo_root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    token_pattern = re.compile("|".join(re.escape(token) for token in OLD_CASE_TOKENS) + r"|prod_v1")
    for path in tracked_files(repo_root):
        if path.suffix.lower() not in TEXT_SUFFIXES or not path.is_file():
            continue
        rel_path = relative(path, repo_root)
        rel_lower = rel_path.lower()
        if rel_lower.startswith(f"{AUDIT_PACKAGE_ROOT.as_posix()}/"):
            # Do not let the audit's own generated CSVs become new references
            # when the committed package is replayed in CI.
            continue
        # Raw result/figure CSVs and unrelated figure registries are already
        # covered by the tree manifests.  Restrict text scanning to code,
        # governance, and RS phase-guided manifests so the audit remains a
        # bounded consumer inventory rather than a duplicate raw-data dump.
        scan_path = (
            rel_lower.startswith("scripts/")
            or rel_lower.startswith(".github/")
            or rel_lower.startswith("config/")
            or rel_lower.startswith("docs/dev/")
            or rel_lower.startswith("docs/guides/")
            or rel_lower in {"readme.md", "docs/reference/implemented_capabilities.md"}
            or rel_lower.startswith("docs/analysis/relaxtime/phase_guided_transport/")
            or rel_lower.startswith("docs/analysis/relaxtime/issue130_rs_")
            or rel_lower.startswith("data/outputs/results/relaxtime/transport/phase_guided/")
            or rel_lower.startswith("data/outputs/figures/relaxtime/transport/phase_guided/")
        )
        if not scan_path:
            continue
        if (rel_lower.startswith("data/outputs/results/") or rel_lower.startswith("data/outputs/figures/")) and path.name not in {
            "manifest.json",
            "README.md",
            "PRODUCTION_AUDIT.md",
            "effective_config.json",
            "plot_manifest.json",
            "production_registry.json",
        }:
            continue
        if path.stat().st_size > 16 * 1024 * 1024:
            continue
        try:
            lines = path.read_text(encoding="utf-8-sig").splitlines()
        except (OSError, UnicodeDecodeError):
            continue
        for line_number, line in enumerate(lines, 1):
            matches = sorted(set(match.group(0) for match in token_pattern.finditer(line)))
            # The repository also contains unrelated PNJL figure-v1 assets.
            # A bare ``prod_v1`` token is relevant here only when the line or
            # path carries the RS phase-guided context; explicit RS case
            # slugs are always retained.
            if matches and "prod_v1" in matches:
                relevant_context = (
                    "phase_guided" in line.lower()
                    or "phase_guided" in rel_path.lower()
                    or "rs_transport" in line.lower()
                    or "phase-guided" in line.lower()
                    or "relaxtime/transport" in line.lower()
                    or "issue130_rs" in rel_path.lower()
                    or "issue130-rs" in rel_path.lower()
                )
                if not relevant_context:
                    matches.remove("prod_v1")
            for token in matches:
                classification = occurrence_class(rel_path)
                rows.append(
                    {
                        "path": rel_path,
                        "line": line_number,
                        "token": token,
                        "classification": classification,
                        "role": occurrence_role(rel_path, line, classification),
                        "snippet": " ".join(line.strip().split())[:240].rstrip(),
                    }
                )
    return rows


def build_registry_rows(registry: dict[str, Any], tree_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_tree = {(row["case_slug"], row["mode"], row["layer"]): row for row in tree_rows}
    output = []
    for entry in registry.get("entries", []):
        if not isinstance(entry, dict):
            continue
        case = str(entry.get("case_slug", ""))
        mode = str(entry.get("mode", ""))
        result = by_tree.get((case, mode, "result"), {})
        output.append(
            {
                "case_slug": case,
                "mode": mode,
                "registry_status": entry.get("status", ""),
                "superseded_by": entry.get("superseded_by", ""),
                "manuscript_eligible": entry.get("manuscript_eligible", ""),
                "registry_result_path": entry.get("result_path", ""),
                "registry_figure_path": entry.get("figure_path", ""),
                "manifest_registry_status": result.get("manifest_registry_status", ""),
                "status_consistency": result.get("status_consistency", "missing_tree"),
                "result_tree_sha256": result.get("tree_sha256", ""),
                "source_commit": entry.get("source_commit", ""),
                "calculation_sha": entry.get("calculation_sha", ""),
            }
        )
    return output


def fallback_rows(repo_root: Path, tree_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_key = {(row["case_slug"], row["mode"], row["layer"]): row for row in tree_rows}
    rows = []
    for mode in MODES:
        current = by_key.get((CURRENT_CASE, mode, "result"), {})
        legacy = by_key.get((LEGACY_CASE, mode, "result"), {})
        current_manifest = result_case_root(repo_root, mode, CURRENT_CASE) / "manifest.json"
        manifest = read_json(current_manifest) if current_manifest.is_file() else {}
        expected = str(manifest.get("legacy_prod_v1_tree_hash", ""))
        actual = str(legacy.get("tree_sha256", ""))
        checks = {
            "current_manifest_exists": current_manifest.is_file(),
            "legacy_tree_exists": bool(actual),
            "legacy_hash_recorded": bool(expected),
            "legacy_hash_matches": bool(expected and expected == actual),
            "current_tree_contract": current.get("tree_contract_ok") is True,
            "legacy_tree_contract": legacy.get("tree_contract_ok") is True,
            "legacy_registry_superseded": legacy.get("registry_status") == "superseded_for_manuscript",
            "legacy_fallback_enabled": manifest.get("legacy_fallback") is True,
            "explicit_rollback_recorded": manifest.get("explicit_rollback") == "--phase-reference-mode legacy",
            "no_production_write": manifest.get("production_write") is False,
            "solver_free_audit": True,
        }
        rows.append(
            {
                "mode": mode,
                "current_case": CURRENT_CASE,
                "legacy_case": LEGACY_CASE,
                "expected_legacy_tree_sha256": expected,
                "actual_legacy_tree_sha256": actual,
                "checks_passed": sum(bool(value) for value in checks.values()),
                "checks_total": len(checks),
                "verdict": "pass" if all(checks.values()) else "fail",
                "check_details": json.dumps(checks, sort_keys=True),
            }
        )
    return rows


def default_consumer_rows(repo_root: Path) -> list[dict[str, Any]]:
    checks = [
        (
            "scripts/analysis/relaxtime/build_phase_guided_transport_xi001_jump_analysis.py",
            "current_case_default",
            CURRENT_CASE,
        ),
        (
            "scripts/analysis/relaxtime/phase_guided_p128_mechanism_scan.jl",
            "current_case_default",
            CURRENT_CASE,
        ),
        (
            "scripts/relaxtime/phase_reference_adapter.jl",
            "legacy_snapshot_and_source_switch",
            "legacy_phase_reference_v1",
        ),
        (
            "scripts/relaxtime/phase_reference_adapter.jl",
            "explicit_source_switch",
            "source in (:candidate, :legacy)",
        ),
        (
            "scripts/relaxtime/run_phase_guided_transport_scan.jl",
            "phase_reference_mode_parameter",
            "phase_reference_mode",
        ),
    ]
    rows = []
    for rel_path, check_name, needle in checks:
        path = repo_root / rel_path
        text = path.read_text(encoding="utf-8") if path.is_file() else ""
        found = needle in text
        old_case_occurrences = sum(text.count(token) for token in OLD_CASE_TOKENS)
        contract_ok = found and old_case_occurrences == 0
        rows.append(
            {
                "path": rel_path,
                "check": check_name,
                "needle": needle,
                "verdict": "pass" if contract_ok else "fail",
                "old_case_occurrences": old_case_occurrences,
                "solver_called": False,
            }
        )
    return rows


def allowlist_rows(tree_rows: list[dict[str, Any]], fallback: list[dict[str, Any]], active_count: int) -> list[dict[str, Any]]:
    fallback_by_mode = {row["mode"]: row for row in fallback}
    rows = []
    for tree in tree_rows:
        if tree["family"] == "current_prod_v2":
            action = "retain_canonical_current"
            deletion = "not_applicable"
            reason = "current approved prod_v2 raw/figure input"
        elif tree["family"] == "legacy_prod_v1":
            action = tree["recommended_action"]
            deletion = "blocked_until_explicit_deletion_gate"
            reason = "explicit fallback/rollback and byte-preserving historical evidence"
        elif tree["family"] == "historical_prod_v1":
            action = tree["recommended_action"]
            deletion = "blocked_until_historical_dependency_review"
            reason = "older convergence/publication evidence; not the current runtime fallback"
        else:
            action = "hold_for_manual_review"
            deletion = "blocked"
            reason = "unclassified prod_v1 tree"
        rows.append(
            {
                "layer": tree["layer"],
                "mode": tree["mode"],
                "case_slug": tree["case_slug"],
                "path": tree["path"],
                "active_old_reference_occurrences": active_count,
                "recommended_action": action,
                "physical_deletion": deletion,
                "reason": reason,
                "tree_sha256": tree["tree_sha256"],
                "fallback_hash_match": fallback_by_mode.get(tree["mode"], {}).get("verdict", "") if tree["family"] == "legacy_prod_v1" and tree["layer"] == "result" else "",
            }
        )
    return rows


def write_readme(output: Path, payload: dict[str, Any], tree_rows: list[dict[str, Any]], fallback: list[dict[str, Any]], consumers: list[dict[str, Any]], references: list[dict[str, Any]]) -> None:
    family_counts = Counter(row["family"] for row in tree_rows)
    active = [row for row in references if row["classification"] == "active_consumer_contract"]
    lines = [
        "# Issue #130 RS old-reference retirement audit v1",
        "",
        "这是 RS transport 旧 `prod_v1` 的 solver-free retirement audit。审计只读取仓库内已入库的 result、figure、registry、manifest、代码和文档；不调用 Julia/equilibrium/transport solver，不改写或删除任何数值产物。",
        "",
        f"- audit schema: `{SCHEMA_VERSION}`",
        f"- audit repo HEAD: `{payload['repo_head']}`",
        f"- calculation SHA: `{CALCULATION_SHA}`",
        f"- workflow head: `{WORKFLOW_HEAD_SHA}`",
        f"- verdict: `{payload['verdict']}`",
        f"- solver_called: `{payload['solver_called']}`",
        f"- production_write: `{payload['production_write']}`",
        "",
        "## 范围和结论",
        "",
        "本审计的 old reference 指 RS phase-guided transport 的 `first_canonical_v2...prod_v1` mode-A/mode-B。更早的 `first_canonical_v1`/`xi005` 树另列为历史证据，不自动纳入当前 runtime fallback 的删除范围。PNJL phase-reference 的 `legacy_phase_reference_v1` 是另一条已经完成 canonical-path retirement 的链，本包不修改它。",
        "",
        f"当前树族计数：`{dict(sorted(family_counts.items()))}`。旧 `prod_v1` 结果的 fallback tree hash 与两套当前 `prod_v2` manifest 的记录逐模式核对；详见 `fallback_smoke.csv`。",
        f"较早历史树的 manifest hash mismatch 计数为 `{payload['historical_manifest_mismatch_count']}`；这些只作为历史 provenance warning，不阻断当前 fallback tree 的完整性门禁，详见 `reference_inventory.csv`。",
        "",
        "审计建议是：旧 `prod_v1` 可以继续评估从 active/canonical 默认路径退出，但必须保留不可变、带哈希的 legacy snapshot 与显式 rollback；本包不授权物理删除。",
        "",
        "## Consumer smoke",
        "",
        f"默认/回退静态合同检查：`{sum(row['verdict'] == 'pass' for row in consumers)}/{len(consumers)}` 通过；旧 reference 的 active consumer occurrence 为 `{len(active)}`。这不是数值 parity 或 numerical convergence 证明。",
        "",
        "- 默认分析入口应解析 current `prod_v2`。",
        "- legacy 只通过显式 versioned snapshot/fallback 访问。",
        "- `production_registry.json` 的 approved/superseded 状态与 raw manifest 内的导入/runtime 状态分开记录；`metadata_review_required` 不能直接当作数据损坏。",
        "",
        "## Retention boundary",
        "",
        "| 层 | 处理 | 删除状态 |",
        "| --- | --- | --- |",
        "| current `prod_v2` | 保留 canonical | 不适用 |",
        "| current fallback `prod_v1` | 可退出 active/canonical，迁入 versioned legacy snapshot | 未授权，不删除 |",
        "| older `prod_v1`/xi005 | 保留历史和收敛证据 | 需独立历史依赖审核 |",
        "",
        "## Evidence files",
        "",
        "- `reference_inventory.csv`: 每棵 result/figure tree 的文件数、字节数、tree/manifest hash、schema/finite/duplicate 状态。",
        "- `registry_consistency.csv`: registry 与各 tree manifest 的状态配对和 superseded 关系。",
        "- `consumer_reference_map.csv`: tracked text 中旧 token 的逐行引用、分类和 review role（active consumer / fallback-or-rollback / historical 等）。",
        "- `fallback_smoke.csv`: current `prod_v2` 到旧 `prod_v1` 的逐模式 hash/fallback/rollback 检查。",
        "- `retention_deletion_allowlist.csv`: 只读的保留、路径退役和物理删除边界。",
        "- `rollback_plan.md`: 后续 implementation PR 的回退顺序；没有删除命令。",
        "- `claim_ledger.csv`: 可支持结论与未声明结论。",
        "",
        "## Stop conditions",
        "",
        "任一 active old consumer、current/fallback tree hash 不匹配、未能解析 versioned snapshot、或历史证据依赖未厘清时，停止 implementation PR。较早历史树的已知 manifest drift 必须原样记录，不能改写成当前 fallback 通过。不得删除旧树、不得重跑 solver，也不得把本包写成 RS 数值收敛证据。",
    ]
    (output / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def write_rollback_plan(output: Path) -> None:
    text = """# RS old-reference retirement rollback plan\n\n本文件只定义后续 implementation PR 的回退顺序，不执行删除。\n\n1. 保持 current `prod_v2` registry entry、分析默认值和 figure path 不变。\n2. 若要退出旧 `prod_v1` canonical path，先把 raw/figure tree 做 byte-preserving versioned snapshot，并写入 source path、tree SHA-256、manifest SHA-256。\n3. 更新 resolver/registry/documentation，使显式 legacy fallback 指向该 snapshot；默认解析仍必须指向 `prod_v2`。\n4. 在 solver-free smoke 中验证 default、explicit legacy 和 rollback 三条路径。\n5. 若任一检查失败，恢复 registry/path pointers；不改 raw bytes。\n6. 物理删除必须另开 PR，携带 deletion allowlist、历史依赖审计和作者的单独授权。\n\n回退原则：先恢复 pointer，再恢复 active path，最后才考虑文件操作；任何阶段都不重算数值。\n"""
    (output / "rollback_plan.md").write_text(text, encoding="utf-8", newline="\n")


def build_claims(payload: dict[str, Any], fallback: list[dict[str, Any]], consumers: list[dict[str, Any]], references: list[dict[str, Any]]) -> list[dict[str, Any]]:
    fallback_ok = all(row["verdict"] == "pass" for row in fallback)
    consumer_ok = all(row["verdict"] == "pass" for row in consumers)
    active_count = sum(row["classification"] == "active_consumer_contract" for row in references)
    return [
        {
            "claim_id": "prod_v2_default_contract",
            "claim": "静态默认分析入口包含 current prod_v2，且没有把旧 prod_v1 作为默认 case。",
            "status": "supported" if consumer_ok else "inconclusive",
            "evidence": "consumer_reference_map.csv; consumer_smoke.csv; production_registry.json",
            "boundary": "solver-free contract only",
        },
        {
            "claim_id": "legacy_fallback_integrity",
            "claim": "两套 current prod_v2 manifest 记录的 legacy prod_v1 tree hash 与仓库实际 tree 一致。",
            "status": "supported" if fallback_ok else "inconclusive",
            "evidence": "fallback_smoke.csv; reference_inventory.csv",
            "boundary": "does not certify numerical parity",
        },
        {
            "claim_id": "no_active_old_consumer",
            "claim": "tracked active consumer contracts 没有直接引用旧 RS prod_v1 case。",
            "status": "supported" if active_count == 0 else "author_check",
            "evidence": "consumer_reference_map.csv",
            "boundary": "historical tooling/docs are retained and separately classified",
        },
        {
            "claim_id": "canonical_path_retirement_candidate",
            "claim": "旧 RS prod_v1 可以进入 active/canonical path retirement 的 implementation review。",
            "status": "candidate" if fallback_ok and consumer_ok and active_count == 0 else "inconclusive",
            "evidence": "retention_deletion_allowlist.csv; fallback_smoke.csv; consumer_smoke.csv",
            "boundary": "requires a separate implementation PR and author review",
        },
        {
            "claim_id": "physical_deletion",
            "claim": "旧 RS prod_v1 可以物理删除。",
            "status": "not_claimed",
            "evidence": "retention_deletion_allowlist.csv; rollback_plan.md",
            "boundary": "legacy snapshot, historical evidence and diagnostic-only numerical status remain",
        },
        {
            "claim_id": "solver_free_boundary",
            "claim": "本审计没有调用 solver、写 production/reference 或改写 raw result。",
            "status": "supported",
            "evidence": "manifest.json; command_manifest.json",
            "boundary": "no numerical conclusion",
        },
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    output = args.output_root.resolve()
    registry_path = repo_root / REGISTRY_PATH
    registry = read_json(registry_path)
    tree_rows = collect_trees(repo_root, registry_index(registry))
    references = scan_references(repo_root)
    fallback = fallback_rows(repo_root, tree_rows)
    consumers = default_consumer_rows(repo_root)
    active_count = sum(row["classification"] == "active_consumer_contract" for row in references)
    fallback_ok = all(row["verdict"] == "pass" for row in fallback)
    consumer_ok = all(row["verdict"] == "pass" for row in consumers)
    integrity_blockers = [
        row
        for row in tree_rows
        if row["family"] in {"current_prod_v2", "legacy_prod_v1"}
        and not row["tree_contract_ok"]
    ]
    if active_count:
        verdict = "retirement_audit_blocked_active_consumer"
    elif not fallback_ok or not consumer_ok or integrity_blockers:
        verdict = "retirement_audit_inconclusive"
    else:
        verdict = "retirement_audit_pass_retain_legacy"
    payload: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "repo_head": git_head(repo_root),
        "calculation_sha": CALCULATION_SHA,
        "workflow_head_sha": WORKFLOW_HEAD_SHA,
        "current_case": CURRENT_CASE,
        "legacy_case": LEGACY_CASE,
        "legacy_snapshot_version": LEGACY_SNAPSHOT_VERSION,
        "legacy_snapshot_result_root": LEGACY_RESULT_ROOT.as_posix(),
        "legacy_snapshot_figure_root": LEGACY_FIGURE_ROOT.as_posix(),
        "historical_cases": list(HISTORICAL_CASES),
        "verdict": verdict,
        "solver_called": False,
        "production_write": False,
        "old_reference_deleted": False,
        "active_old_reference_occurrences": active_count,
        "fallback_smoke_pass": fallback_ok,
        "consumer_smoke_pass": consumer_ok,
        "canonical_path_retirement_candidate": bool(fallback_ok and consumer_ok and active_count == 0),
        "active_legacy_integrity_blocker_count": len(integrity_blockers),
        "historical_manifest_mismatch_count": sum(
            int(row["manifest_hash_mismatch_count"] or 0)
            for row in tree_rows
            if row["family"] == "historical_prod_v1"
        ),
        "physical_deletion_gate": "blocked_diagnostic_only_or_historical_dependency",
        "scope_boundary": "RS transport prod_v1 only; PNJL phase-reference legacy snapshot is out of scope",
    }
    output.mkdir(parents=True, exist_ok=True)
    write_csv(
        output / "reference_inventory.csv",
        tree_rows[0].keys() if tree_rows else ("path",),
        tree_rows,
    )
    registry_rows = build_registry_rows(registry, tree_rows)
    write_csv(
        output / "registry_consistency.csv",
        registry_rows[0].keys() if registry_rows else ("case_slug",),
        registry_rows,
    )
    write_csv(
        output / "consumer_reference_map.csv",
        ("path", "line", "token", "classification", "role", "snippet"),
        references,
    )
    write_csv(
        output / "fallback_smoke.csv",
        fallback[0].keys() if fallback else ("mode",),
        fallback,
    )
    write_csv(
        output / "consumer_smoke.csv",
        consumers[0].keys() if consumers else ("path",),
        consumers,
    )
    allowlist = allowlist_rows(tree_rows, fallback, active_count)
    write_csv(
        output / "retention_deletion_allowlist.csv",
        allowlist[0].keys() if allowlist else ("path",),
        allowlist,
    )
    claims = build_claims(payload, fallback, consumers, references)
    write_csv(output / "claim_ledger.csv", claims[0].keys(), claims)
    write_json(output / "registry_snapshot.json", registry)
    write_rollback_plan(output)
    payload["generator"] = "scripts/analysis/relaxtime/audit_issue130_rs_old_reference_retirement.py"
    payload["generator_sha256"] = sha256(repo_root / payload["generator"])
    payload["registry_sha256"] = sha256(registry_path)
    payload["tree_count"] = len(tree_rows)
    payload["reference_occurrence_count"] = len(references)
    payload["claims"] = {row["claim_id"]: row["status"] for row in claims}
    write_json(output / "audit_summary.json", payload)
    write_json(
        output / "command_manifest.json",
        {
            "command": "python scripts/analysis/relaxtime/audit_issue130_rs_old_reference_retirement.py",
            "repo_root": str(repo_root),
            "output_root": str(output),
            "repo_head": payload["repo_head"],
            "solver_called": False,
            "production_write": False,
        },
    )
    write_readme(output, payload, tree_rows, fallback, consumers, references)
    output_files = [
        "reference_inventory.csv",
        "registry_consistency.csv",
        "consumer_reference_map.csv",
        "fallback_smoke.csv",
        "consumer_smoke.csv",
        "retention_deletion_allowlist.csv",
        "claim_ledger.csv",
        "registry_snapshot.json",
        "rollback_plan.md",
        "audit_summary.json",
        "command_manifest.json",
        "README.md",
    ]
    write_json(
        output / "manifest.json",
        {
            "schema_version": f"{SCHEMA_VERSION}_manifest",
            "source_repo_head": payload["repo_head"],
            "calculation_sha": CALCULATION_SHA,
            "workflow_head_sha": WORKFLOW_HEAD_SHA,
            "solver_called": False,
            "production_write": False,
            "old_reference_deleted": False,
            "files": [
                {"path": name, "bytes": (output / name).stat().st_size, "sha256": sha256(output / name)}
                for name in output_files
            ],
        },
    )
    print(
        json.dumps(
            {
                "verdict": verdict,
                "tree_count": len(tree_rows),
                "active_old_reference_occurrences": active_count,
                "fallback_smoke_pass": fallback_ok,
                "consumer_smoke_pass": consumer_ok,
                "solver_called": False,
            },
            ensure_ascii=False,
        )
    )
    return 0 if verdict == "retirement_audit_pass_retain_legacy" else 2


if __name__ == "__main__":
    raise SystemExit(main())
