#!/usr/bin/env python3
"""Import author-reviewed Issue #130 RS candidate figures into formal paths.

The numerical source remains an immutable external diagnostic artifact.  This
script copies only the candidate PNG trees, verifies the source and per-image
hashes, and writes an explicit figure provenance manifest.  It never calls a
solver and refuses to overwrite an existing formal case directory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path
from typing import Any


CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
SOURCE_MANIFEST_SHA256 = "9a986e1887292309a46a963dfbc76f08cc7fb67d3fa5664f282ee7256853f01b"
CASE_NAME = "first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2"
EXPECTED_FIGURE_COUNT = 36
MODES = (
    "mode_a_fixed_muB_phase_scaled",
    "mode_b_fixed_T_sparse_muB",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"JSON object required: {path}")
    return payload


def _inside(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
    except ValueError:
        return False
    return True


def _figure_hashes(case_root: Path) -> dict[str, str]:
    return {
        path.relative_to(case_root).as_posix(): sha256(path)
        for path in sorted(case_root.rglob("*.png"))
    }


def validate_source_manifest(source_root: Path) -> tuple[dict[str, Any], str]:
    source_root = source_root.resolve()
    manifest_path = source_root / "manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(manifest_path)
    manifest_sha = sha256(manifest_path)
    if manifest_sha != SOURCE_MANIFEST_SHA256:
        raise ValueError(
            f"source manifest hash {manifest_sha} != expected {SOURCE_MANIFEST_SHA256}"
        )
    manifest = _read_json(manifest_path)
    if manifest.get("schema_version") != "issue130_rs_candidate_legacy_figures_v1":
        raise ValueError("unexpected source figure manifest schema")
    if manifest.get("calculation_sha") != CALCULATION_SHA:
        raise ValueError("source figure calculation SHA mismatch")
    if manifest.get("production_write") is not False or manifest.get("solver_called") is not False:
        raise ValueError("source figure manifest must be solver-free and non-production")
    return manifest, manifest_sha


def validate_candidate_tree(
    *,
    source_root: Path,
    mode: str,
    rendered: dict[str, Any],
) -> dict[str, Any]:
    if rendered.get("mode") != mode or rendered.get("reference") != "candidate_runtime":
        raise ValueError(f"candidate rendered record mismatch for {mode}")
    case_root = Path(rendered["figure_root"]).resolve()
    if not _inside(case_root, source_root):
        raise ValueError(f"candidate figure root escapes source root: {case_root}")
    plot_manifest_path = case_root / "plot_manifest.json"
    if not plot_manifest_path.is_file():
        raise FileNotFoundError(plot_manifest_path)
    plot_manifest_sha = sha256(plot_manifest_path)
    plot_manifest = _read_json(plot_manifest_path)
    figure_hashes = plot_manifest.get("figure_hashes")
    if not isinstance(figure_hashes, dict) or len(figure_hashes) != EXPECTED_FIGURE_COUNT:
        raise ValueError(f"{mode} candidate plot manifest must contain 36 figure hashes")
    actual_hashes = _figure_hashes(case_root)
    if actual_hashes != figure_hashes:
        raise ValueError(f"{mode} candidate figure hash map does not match its files")
    if plot_manifest.get("count") != EXPECTED_FIGURE_COUNT:
        raise ValueError(f"{mode} candidate plot count is not 36")
    if rendered.get("figure_hashes") != figure_hashes:
        raise ValueError(f"{mode} root and tree figure hash maps disagree")
    if rendered.get("figure_count") != EXPECTED_FIGURE_COUNT:
        raise ValueError(f"{mode} root figure count is not 36")
    return {
        "mode": mode,
        "source_case_root": str(case_root),
        "source_plot_manifest": str(plot_manifest_path),
        "source_plot_manifest_sha256": plot_manifest_sha,
        "source_csv": plot_manifest.get("source_csv"),
        "source_csv_sha256": plot_manifest.get("source_csv_sha256"),
        "plot_contract": {
            "entrypoint": plot_manifest.get("plot_entrypoint"),
            "x": plot_manifest.get("x"),
            "y_columns": plot_manifest.get("y_columns"),
            "format": plot_manifest.get("format"),
            "dpi": plot_manifest.get("dpi"),
        },
        "figure_hashes": figure_hashes,
    }


def import_candidate_tree(
    *,
    repo_root: Path,
    target_root: Path,
    source: dict[str, Any],
    source_manifest_sha256: str,
    workflow_head_sha: str,
    author_review_note: str,
) -> dict[str, Any]:
    mode = source["mode"]
    source_case_root = Path(source["source_case_root"])
    target_case_root = target_root / mode / CASE_NAME
    if target_case_root.exists():
        raise FileExistsError(
            f"refusing to overwrite existing formal figure case: {target_case_root}"
        )
    target_case_root.mkdir(parents=True, exist_ok=False)
    for relative_name in sorted(source["figure_hashes"]):
        source_path = source_case_root / Path(relative_name)
        target_path = target_case_root / Path(relative_name)
        target_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_path, target_path)

    copied_hashes = _figure_hashes(target_case_root)
    if copied_hashes != source["figure_hashes"]:
        raise ValueError(f"copied figure hash map does not match source for {mode}")

    manifest = {
        "schema_version": "issue130_rs_candidate_runtime_formal_figure_v1",
        "figure_status": "author_accepted_formal_layout",
        "numerical_status": "diagnostic_only",
        "calculation_sha": CALCULATION_SHA,
        "workflow_head_sha": workflow_head_sha,
        "source_audit_manifest_sha256": source_manifest_sha256,
        "source_external_root": str(source_case_root.parents[2]),
        "source_case_root": str(source_case_root),
        "source_plot_manifest": source["source_plot_manifest"],
        "source_plot_manifest_sha256": source["source_plot_manifest_sha256"],
        "source_csv": source["source_csv"],
        "source_csv_sha256": source["source_csv_sha256"],
        "target_case_root": target_case_root.relative_to(repo_root).as_posix(),
        "figure_import": True,
        "solver_called": False,
        "production_write": False,
        "author_review_note": author_review_note,
        "plot_contract": source["plot_contract"],
        "figures": sorted(copied_hashes),
        "figure_count": len(copied_hashes),
        "figure_hashes": copied_hashes,
    }
    manifest_path = target_case_root / "plot_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        "mode": mode,
        "target_case_root": str(target_case_root),
        "plot_manifest": str(manifest_path),
        "plot_manifest_sha256": sha256(manifest_path),
        "figure_count": len(copied_hashes),
        "figure_hashes": copied_hashes,
    }


def formalize(args: argparse.Namespace) -> dict[str, Any]:
    repo_root = Path(args.repo_root).resolve()
    source_root = Path(args.source_root).resolve()
    target_root = Path(args.target_root).resolve()
    official_root = (repo_root / "data" / "outputs" / "figures").resolve()
    if not _inside(target_root, official_root):
        raise ValueError("target root must be inside data/outputs/figures")
    source_manifest, source_manifest_sha256 = validate_source_manifest(source_root)
    rendered = {
        item["mode"]: item
        for item in source_manifest.get("rendered", [])
        if item.get("reference") == "candidate_runtime"
    }
    if set(rendered) != set(MODES):
        raise ValueError("source manifest must contain exactly the two candidate modes")

    sources = {
        mode: validate_candidate_tree(
            source_root=source_root,
            mode=mode,
            rendered=rendered[mode],
        )
        for mode in MODES
    }
    imported = [
        import_candidate_tree(
            repo_root=repo_root,
            target_root=target_root,
            source=sources[mode],
            source_manifest_sha256=source_manifest_sha256,
            workflow_head_sha=args.workflow_head_sha,
            author_review_note=args.author_review_note,
        )
        for mode in MODES
    ]
    return {
        "schema_version": "issue130_rs_candidate_runtime_formalization_v1",
        "calculation_sha": CALCULATION_SHA,
        "workflow_head_sha": args.workflow_head_sha,
        "source_external_root": str(source_root),
        "source_manifest_sha256": source_manifest_sha256,
        "target_root": str(target_root),
        "numerical_status": "diagnostic_only",
        "figure_status": "author_accepted_formal_layout",
        "solver_called": False,
        "production_write": False,
        "imported": imported,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--target-root", required=True, type=Path)
    parser.add_argument("--workflow-head-sha", required=True)
    parser.add_argument(
        "--author-review-note",
        default="author_accepted_full_candidate_legacy_comparison_and_review_figures",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    print(json.dumps(formalize(args), ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
