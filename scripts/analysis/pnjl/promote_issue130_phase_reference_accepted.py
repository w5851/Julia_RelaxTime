#!/usr/bin/env python3
"""Record the author-approved downstream default for the Issue #130 v2 package.

This is a solver-free metadata promotion.  It updates only the v2 ``accepted``
view and its manifests/claim ledger; strict tables, render tables, the legacy
snapshot and the solver runtime contract remain unchanged.  The command is
deliberately explicit so a future rebuild cannot silently turn a pending
candidate into an accepted downstream input.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


PACKAGE_SCHEMA = "pnjl_issue130_phase_reference_v2"
LAYER_SCHEMA = "pnjl_issue130_phase_reference_layer_v2"
DEFAULT_ROOT = Path("data/reference/pnjl/issue130_phase_reference_v2")
TABLES = {
    "boundary": "maxwell_surface_accepted_phase_map_v1.csv",
    "crossover": "crossover_surface_accepted_phase_map_v1.csv",
    "cep": "cep_boundary_accepted_phase_map_v1.csv",
    "spinodals": "spinodal_surface_accepted_phase_map_v1.csv",
}
PENDING_STATUS = "candidate_pending_author_review"
ACCEPTED_STATUS = "author_accepted_for_downstream"
PENDING_SCOPE = "downstream_phase_map_candidate"
ACCEPTED_SCOPE = "downstream_phase_map_default"


class AcceptedPromotionError(ValueError):
    """Raised when the v2 package cannot satisfy the promotion contract."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AcceptedPromotionError(f"invalid JSON: {path}") from exc
    if not isinstance(value, dict):
        raise AcceptedPromotionError(f"expected JSON object: {path}")
    return value


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def file_record(path: Path, root: Path) -> dict[str, Any]:
    return {
        "path": path.relative_to(root).as_posix(),
        "sha256": sha256(path),
        "bytes": path.stat().st_size,
    }


def _accepted_table_path(root: Path, table: str, manifest: dict[str, Any]) -> Path:
    record = manifest.get("tables", {}).get(table, {}).get("output", {})
    relative = record.get("path")
    if not isinstance(relative, str) or not relative:
        relative = f"accepted/tables/{TABLES[table]}"
    path = (root / relative).resolve()
    try:
        path.relative_to(root.resolve())
    except ValueError as exc:
        raise AcceptedPromotionError(f"accepted table escapes package root: {relative}") from exc
    return path


def _promote_table(path: Path) -> tuple[int, int]:
    try:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            fields = list(reader.fieldnames or [])
            rows = [dict(row) for row in reader]
    except OSError as exc:
        raise AcceptedPromotionError(f"cannot read accepted table: {path}") from exc
    required = {"acceptance_status", "acceptance_scope", "source_status", "extrapolation"}
    missing = sorted(required - set(fields))
    if missing:
        raise AcceptedPromotionError(f"{path} missing promotion fields: {', '.join(missing)}")

    changed = 0
    for row in rows:
        status = row.get("acceptance_status", "")
        if status == PENDING_STATUS:
            row["acceptance_status"] = ACCEPTED_STATUS
            changed += 1
        elif status != ACCEPTED_STATUS:
            raise AcceptedPromotionError(
                f"unexpected acceptance_status {status!r} in {path}; refusing partial promotion"
            )
        scope = row.get("acceptance_scope", "")
        if scope == PENDING_SCOPE:
            row["acceptance_scope"] = ACCEPTED_SCOPE
        elif scope != ACCEPTED_SCOPE:
            raise AcceptedPromotionError(
                f"unexpected acceptance_scope {scope!r} in {path}; refusing partial promotion"
            )
        if row.get("extrapolation", "").strip().lower() not in {"false", "0", "no"}:
            raise AcceptedPromotionError(f"accepted table has extrapolated row: {path}")

    if changed:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    return len(rows), changed


def _update_text_files(root: Path) -> None:
    (root / "accepted" / "README.md").write_text(
        "# accepted\n\n"
        "This is the author-accepted downstream phase-map view generated from "
        "structured render tables. Native, unresolved and interpolated source "
        "states remain explicitly labelled; `author_accepted_for_downstream` "
        "does not imply strict certification. This layer is the default for "
        "downstream phase-map/analysis consumers only and keeps "
        "`runtime_consumption=false`.\n",
        encoding="utf-8",
        newline="\n",
    )
    (root / "README.md").write_text(
        "# Issue #130 phase-reference v2\n\n"
        "The public semantic layers are `strict`, `render` and `accepted`. The "
        "former v1 `derived` layer remains an internal build input and provenance "
        "record. `render` is the complete structured display package; `accepted` "
        "is the author-accepted default for downstream phase-map/analysis "
        "consumers. All layers are solver-free and `runtime_consumption=false`; "
        "the solver runtime continues to use strict certified rows with the "
        "explicit legacy fallback/rollback contract.\n",
        encoding="utf-8",
        newline="\n",
    )
    (root / "AUDIT.md").write_text(
        "# Issue #130 phase-reference v2 solver-free materialization audit\n\n"
        "The author has accepted the `accepted` layer as the downstream analysis "
        "default. This does not upgrade interpolated rows to strict certification "
        "and does not change solver runtime selection.\n\n"
        "- strict/render/accepted source values remain unchanged except for the "
        "accepted row-level decision metadata;\n"
        "- `solver_called=false`, `reference_write=false`, and "
        "`runtime_consumption=false`;\n"
        "- strict runtime and the explicit PNJL legacy fallback/rollback remain "
        "available.\n",
        encoding="utf-8",
        newline="\n",
    )


def promote_package(root: Path, *, recorded_at: str | None = None) -> dict[str, Any]:
    root = root.resolve()
    root_manifest_path = root / "manifest.json"
    accepted_manifest_path = root / "accepted" / "manifest.json"
    claim_path = root / "claim_ledger.json"
    root_manifest = read_json(root_manifest_path)
    accepted_manifest = read_json(accepted_manifest_path)
    if root_manifest.get("schema_version") != PACKAGE_SCHEMA:
        raise AcceptedPromotionError("root schema is not Issue #130 phase-reference v2")
    if accepted_manifest.get("schema_version") != LAYER_SCHEMA or accepted_manifest.get("layer") != "accepted":
        raise AcceptedPromotionError("accepted layer manifest schema/layer mismatch")
    if root_manifest.get("runtime_consumption") is not False or root_manifest.get("reference_write") is not False:
        raise AcceptedPromotionError("promotion requires runtime_consumption=false and reference_write=false")
    if accepted_manifest.get("runtime_consumption") is not False or accepted_manifest.get("reference_write") is not False:
        raise AcceptedPromotionError("accepted layer is not solver-free")
    if root_manifest.get("promotion_status") not in {"awaiting_author_review", "accepted_for_downstream"}:
        raise AcceptedPromotionError("root promotion status is not author-review pending")

    timestamp = recorded_at or datetime.now(timezone.utc).isoformat()
    row_counts: dict[str, int] = {}
    changed_rows: dict[str, int] = {}
    for table in TABLES:
        path = _accepted_table_path(root, table, accepted_manifest)
        if not path.is_file():
            raise AcceptedPromotionError(f"missing accepted table: {path}")
        row_counts[table], changed_rows[table] = _promote_table(path)

    accepted_manifest.update(
        {
            "semantics": "author-accepted downstream phase map derived from complete structured render tables",
            "promotion_status": "accepted_for_downstream",
            "acceptance_status": ACCEPTED_STATUS,
            "acceptance_scope": "downstream_phase_map_default",
            "author_decision": {
                "status": "accepted",
                "scope": "downstream_analysis_default",
                "recorded_at_utc": timestamp,
                "runtime_default_unchanged": True,
                "legacy_fallback_preserved": True,
            },
        }
    )
    accepted_manifest["tables"] = dict(accepted_manifest.get("tables", {}))
    for table in TABLES:
        path = _accepted_table_path(root, table, accepted_manifest)
        accepted_manifest["tables"].setdefault(table, {})["output"] = file_record(path, root)
        accepted_manifest["tables"][table]["rows"] = row_counts[table]
    _update_text_files(root)
    write_json(accepted_manifest_path, accepted_manifest)

    claims = read_json(claim_path)
    for claim in claims.get("claims", []):
        if claim.get("claim_id") == "accepted_status_preserved":
            claim.update(
                {
                    "status": "supported",
                    "boundary": "author accepted for downstream analysis only; interpolated rows remain non-certified",
                }
            )
        elif claim.get("claim_id") == "runtime_promotion":
            claim.update(
                {
                    "status": "blocked",
                    "boundary": "accepted is downstream-only; strict runtime and legacy fallback/rollback are unchanged",
                }
            )
    write_json(claim_path, claims)

    root_manifest.update(
        {
            "promotion_status": "accepted_for_downstream",
            "downstream_default_layer": "accepted",
            "downstream_default_scope": "analysis_and_derived_phase_maps",
            "author_decision": {
                "status": "accepted",
                "scope": "downstream_analysis_default",
                "recorded_at_utc": timestamp,
                "runtime_default_unchanged": True,
                "legacy_fallback_preserved": True,
            },
        }
    )
    constraints = dict(root_manifest.get("constraints", {}))
    constraints.update(
        {
            "accepted_for_downstream": True,
            "runtime_default_unchanged": True,
            "legacy_fallback_preserved": True,
            "interpolated_rows_noncertified": True,
        }
    )
    root_manifest["constraints"] = constraints
    root_manifest.setdefault("layers", {})["accepted"] = file_record(accepted_manifest_path, root)
    root_manifest["claim_ledger"] = file_record(claim_path, root)
    write_json(root_manifest_path, root_manifest)

    return {
        "root": str(root),
        "promotion_status": root_manifest["promotion_status"],
        "downstream_default_layer": root_manifest["downstream_default_layer"],
        "changed_rows": changed_rows,
        "row_counts": row_counts,
        "solver_called": False,
        "runtime_consumption": False,
        "reference_write": False,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--recorded-at", default=None, help="UTC timestamp for deterministic audit replay")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = promote_package(args.root, recorded_at=args.recorded_at)
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
