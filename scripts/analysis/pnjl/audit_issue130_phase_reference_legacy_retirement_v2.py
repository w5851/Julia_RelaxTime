#!/usr/bin/env python3
"""Build the post-acceptance, solver-free PNJL legacy coverage audit.

The key coverage calculation is shared with the immutable v1 audit.  This
version records the author-approved ``accepted`` downstream default while
keeping the runtime comparison on strict certified rows.  It never calls a
solver, rewrites either reference tree, or deletes the legacy snapshot.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.analysis.pnjl import audit_issue130_phase_reference_legacy_retirement as audit_v1


SCHEMA_VERSION = "pnjl_issue130_phase_reference_legacy_retirement_audit_v2"
OUTPUT_RELATIVE = Path(
    "docs/analysis/pnjl/phase_reference/issue130_phase_reference_legacy_audit_v2"
)
PACKAGE_RELATIVE = Path("data/reference/pnjl/issue130_phase_reference_v2")


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def sha256(path: Path) -> str:
    return audit_v1.sha256(path)


def _refresh_manifest(output_root: Path, manifest: dict[str, Any]) -> None:
    generated_files = sorted(
        path for path in output_root.rglob("*") if path.is_file() and path.name != "manifest.json"
    )
    manifest["files"] = [
        {
            "path": path.relative_to(output_root).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in generated_files
    ]
    write_json(output_root / "manifest.json", manifest)


def _refresh_repo_head(output_root: Path, repo_root: Path) -> str:
    """Replace v1 prose provenance with the head used for this v2 build."""

    head = audit_v1.git_value(repo_root, "rev-parse", "HEAD")
    for name in ("README.md", "AUDIT.md"):
        path = output_root / name
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8")
        if name == "README.md":
            text = text.replace(
                "# Issue #130 PNJL legacy phase-reference retirement audit v1",
                "# Issue #130 PNJL legacy phase-reference retirement audit v2",
                1,
            )
        text = re.sub(r"本次 repo HEAD：`[0-9a-f]+`", f"本次 repo HEAD：`{head}`", text)
        path.write_text(text, encoding="utf-8", newline="\n")
    return head


def _append_claim(output_root: Path, accepted_status: str) -> None:
    claim_path = output_root / "tables" / "claim_ledger.json"
    claims = json.loads(claim_path.read_text(encoding="utf-8"))
    if not isinstance(claims, list):
        raise ValueError("legacy audit claim ledger must be a list")
    claims.append(
        {
            "claim_id": "accepted_downstream_default",
            "claim": "The author-approved accepted layer is the downstream analysis default.",
            "status": "supported" if accepted_status == "accepted_for_downstream" else "inconclusive",
            "evidence": "data/reference/pnjl/issue130_phase_reference_v2/manifest.json; accepted/manifest.json",
            "boundary": "accepted is not runtime input and does not replace strict certified rows or legacy fallback",
        }
    )
    write_json(claim_path, claims)
    csv_path = output_root / "tables" / "claim_ledger.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(claims[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(claims)


def build_audit(
    repo_root: Path,
    output_root: Path,
    *,
    replace_existing: bool = False,
) -> dict[str, Any]:
    output_root = output_root.resolve()
    if output_root.exists() and any(output_root.iterdir()):
        if not replace_existing:
            raise FileExistsError(f"refusing to overwrite non-empty audit output: {output_root}")
        shutil.rmtree(output_root)

    # Reuse the reviewed v1 implementation for semantic key coverage and
    # consumer scanning, but write to a new immutable v2 evidence directory.
    audit_v1.build_audit(repo_root.resolve(), output_root)
    repo_head = _refresh_repo_head(output_root, repo_root.resolve())
    package_root = (repo_root / PACKAGE_RELATIVE).resolve()
    package_manifest = read_json(package_root / "manifest.json")
    accepted_manifest = read_json(package_root / "accepted" / "manifest.json")
    accepted_status = str(package_manifest.get("promotion_status", ""))
    if accepted_status != "accepted_for_downstream":
        raise ValueError("accepted package is not author-promoted for downstream use")
    if package_manifest.get("downstream_default_layer") != "accepted":
        raise ValueError("accepted package does not declare accepted as downstream default")

    decision_path = output_root / "decision.json"
    decision = read_json(decision_path)
    decision.update(
        {
            "schema_version": "pnjl_issue130_phase_reference_legacy_retirement_decision_v2",
            "audit_version": "v2",
            "accepted_downstream_default": True,
            "downstream_default_layer": "accepted",
            "accepted_promotion_status": accepted_status,
            "accepted_manifest_sha256": sha256(package_root / "accepted" / "manifest.json"),
            "accepted_interpolated_rows_remain_noncertified": True,
            "runtime_reference_layer": "strict",
            "repo_head": repo_head,
            "next_action": (
                "retain legacy snapshot; accepted is the downstream analysis default, while strict "
                "runtime still requires candidate-only migration and explicit fallback audit"
            ),
        }
    )
    write_json(decision_path, decision)
    _append_claim(output_root, accepted_status)

    readme_path = output_root / "README.md"
    readme_path.write_text(
        readme_path.read_text(encoding="utf-8")
        + "\n## Post-acceptance boundary\n\n"
        + "The v2 `accepted` layer is the author-approved default for downstream "
        + "phase-map/analysis consumers. Runtime coverage is still evaluated on "
        + "strict certified keys; accepted interpolation does not eliminate a "
        + "legacy fallback dependency.\n",
        encoding="utf-8",
        newline="\n",
    )
    audit_path = output_root / "AUDIT.md"
    audit_path.write_text(
        audit_path.read_text(encoding="utf-8")
        + "\n## Post-acceptance boundary\n\n"
        + "`accepted_for_downstream` is a downstream data decision only. It does "
        + "not change strict runtime selection, legacy rollback, or the physical "
        + "deletion gate.\n",
        encoding="utf-8",
        newline="\n",
    )

    manifest = read_json(output_root / "manifest.json")
    manifest.update(
        {
            "schema_version": SCHEMA_VERSION,
            "audit_version": "v2",
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "accepted_root": PACKAGE_RELATIVE.as_posix(),
            "accepted_manifest_sha256": sha256(package_root / "accepted" / "manifest.json"),
            "accepted_promotion_status": accepted_status,
            "downstream_default_layer": "accepted",
            "runtime_reference_layer": "strict",
            "accepted_runtime_consumption": accepted_manifest.get("runtime_consumption"),
            "accepted_reference_write": accepted_manifest.get("reference_write"),
            "accepted_interpolated_rows_noncertified": True,
            "repo_head": repo_head,
        }
    )
    _refresh_manifest(output_root, manifest)
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--output-root", type=Path, default=None)
    parser.add_argument(
        "--replace-existing",
        action="store_true",
        help="explicitly rebuild the local uncommitted v2 evidence directory",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    output_root = (args.output_root or repo_root / OUTPUT_RELATIVE).resolve()
    manifest = build_audit(repo_root, output_root, replace_existing=args.replace_existing)
    print(json.dumps({"output_root": str(output_root), "verdict": manifest["verdict"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
