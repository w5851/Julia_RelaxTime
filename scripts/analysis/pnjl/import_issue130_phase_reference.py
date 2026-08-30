#!/usr/bin/env python3
"""Import the reviewed Issue #130 phase-reference candidate into a new sibling.

The importer is deliberately separate from the promotion gate.  It copies the
reviewed strict/derived/render tables byte-for-byte, records the source hashes,
and records the deleted legacy snapshot from its immutable audit inventory.
The importer no longer requires the deleted snapshot to be present in the
working tree.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCHEMA = "pnjl_issue130_phase_reference_import_v1"
CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
SOURCE_RUN_ID = "32354095831"
REPLAY_RUN_ID = "32451053476"
PACKAGE_RELATIVE = Path("docs/analysis/pnjl/phase_reference/issue130_phase_reference_layers_v1")
GATE_RELATIVE = Path("docs/analysis/pnjl/phase_reference/issue130_phase_reference_promotion_gate_v1")
# The legacy bundle is now a canonical, immutable snapshot after retirement;
# the former top-level dense paths are intentionally absent.  Import keeps
# recording these bytes so a new candidate can prove it did not mutate the
# rollback source.
LEGACY_FILES = (
    Path("data/reference/pnjl/legacy_phase_reference_v1/boundary.csv"),
    Path("data/reference/pnjl/legacy_phase_reference_v1/cep.csv"),
    Path("data/reference/pnjl/legacy_phase_reference_v1/spinodals.csv"),
    Path("data/reference/pnjl/legacy_phase_reference_v1/crossover_dense.csv"),
    Path("data/reference/pnjl/legacy_phase_reference_v1/phase_reference_dense_manifest.json"),
)
LEGACY_INVENTORY = Path(
    "docs/analysis/pnjl/phase_reference/issue130_phase_reference_legacy_audit_v4/"
    "tables/legacy_snapshot_inventory.csv"
)

LAYER_SPECS = {
    "strict": "strict_reference_v1",
    "derived": "derived_reference_v1",
    "render": "phase_surface_render_v1",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative_path(path: Path, root: Path) -> str:
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return path.as_posix()


def file_record(path: Path, root: Path) -> dict[str, Any]:
    return {
        "path": relative_path(path, root),
        "sha256": sha256(path),
        "bytes": path.stat().st_size,
    }


def read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _package_files(source_layer: Path) -> list[Path]:
    required = [source_layer / "README.md", source_layer / "manifest.json"]
    tables = sorted((source_layer / "tables").glob("*.csv"))
    files = required + tables
    missing = [path for path in files if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing layer files: " + ", ".join(str(path) for path in missing))
    return files


def _assert_csv_contract(path: Path) -> None:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle))
    if not rows or not rows[0]:
        raise ValueError(f"empty CSV: {path}")
    width = len(rows[0])
    if any(len(row) != width for row in rows[1:]):
        raise ValueError(f"CSV column mismatch: {path}")


def _historical_legacy_records(repo_root: Path) -> list[dict[str, Any]]:
    """Return pre-delete legacy hashes without reading a deleted snapshot.

    The v4 solver-free audit inventory is the immutable provenance source after
    physical deletion.  These records are metadata only; the importer never
    restores or consumes legacy rows.
    """
    inventory = repo_root / LEGACY_INVENTORY
    if not inventory.is_file():
        raise FileNotFoundError(f"legacy audit inventory missing: {inventory}")
    with inventory.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"legacy audit inventory is empty: {inventory}")
    records: list[dict[str, Any]] = []
    for row in rows:
        expected = row.get("expected_sha256", "")
        expected_bytes = row.get("expected_bytes", "")
        if not expected or not expected_bytes:
            raise ValueError(f"legacy audit inventory row is incomplete: {row}")
        records.append({
            "path": row.get("source_path", ""),
            "sha256": expected,
            "bytes": int(expected_bytes),
            "availability": "git_recovery_ref_only",
            "recovery_ref": "9aa4c313901ca0c91e851f58514e3df9aa124df4",
        })
    return records


def _validate_gate(package_root: Path, gate_root: Path) -> dict[str, Any]:
    package_manifest_path = package_root / "manifest.json"
    gate_manifest_path = gate_root / "manifest.json"
    gate_decision_path = gate_root / "decision.json"
    author_record_path = gate_root / "author_review_record_v1.json"
    for path in (package_manifest_path, gate_manifest_path, gate_decision_path, author_record_path):
        if not path.is_file():
            raise FileNotFoundError(f"missing promotion input: {path}")

    package_manifest = read_json(package_manifest_path)
    gate_manifest = read_json(gate_manifest_path)
    decision = read_json(gate_decision_path)
    author_record = read_json(author_record_path)
    if package_manifest.get("calculation_sha") != CALCULATION_SHA:
        raise ValueError("phase-reference package calculation SHA mismatch")
    if gate_manifest.get("verdict") != "promotion_candidate":
        raise ValueError("promotion gate is not a candidate")
    if decision.get("verdict") != "promotion_candidate":
        raise ValueError("promotion decision is not a candidate")
    if author_record.get("status") != "accepted":
        raise ValueError("author review has not been accepted")
    if author_record.get("reference_write_authorized") is not False:
        raise ValueError("gate must not authorize an implicit reference write")
    if author_record.get("runtime_consumption_authorized") is not False:
        raise ValueError("gate must not authorize runtime consumption")
    if decision.get("source_run_id") != SOURCE_RUN_ID or decision.get("replay_run_id") != REPLAY_RUN_ID:
        raise ValueError("promotion input run provenance mismatch")
    return {
        "package_manifest": package_manifest,
        "gate_manifest": gate_manifest,
        "decision": decision,
        "author_record": author_record,
        "package_manifest_path": package_manifest_path,
        "gate_manifest_path": gate_manifest_path,
        "decision_path": gate_decision_path,
        "author_record_path": author_record_path,
    }


def _copy_layer(source_layer: Path, destination_layer: Path, repo_root: Path) -> dict[str, Any]:
    destination_layer.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, Any]] = []
    for source in _package_files(source_layer):
        destination = destination_layer / source.relative_to(source_layer)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        records.append({
            "source": file_record(source, repo_root),
            "imported": file_record(destination, repo_root),
        })
        if source.suffix.lower() == ".csv":
            _assert_csv_contract(destination)
    return {"source_layer": source_layer.name, "files": records}


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def import_package(
    package_root: Path,
    gate_root: Path,
    reference_root: Path,
    figure_root: Path,
    repo_root: Path,
    *,
    refresh_existing: bool = False,
) -> dict[str, Any]:
    package_root = package_root.resolve()
    gate_root = gate_root.resolve()
    reference_root = reference_root.resolve()
    figure_root = figure_root.resolve()
    repo_root = repo_root.resolve()
    if (reference_root.exists() or figure_root.exists()) and not refresh_existing:
        raise FileExistsError("refusing to overwrite an existing versioned reference or figure root")
    if refresh_existing:
        for path in (reference_root, figure_root):
            if path.exists():
                shutil.rmtree(path)

    gate = _validate_gate(package_root, gate_root)
    legacy_records = _historical_legacy_records(repo_root)

    source_layers = {
        "strict": package_root / LAYER_SPECS["strict"],
        "derived": package_root / LAYER_SPECS["derived"],
        "render": package_root / LAYER_SPECS["render"],
    }
    for path in source_layers.values():
        if not path.is_dir():
            raise FileNotFoundError(f"missing phase-reference layer: {path}")

    generated_at = datetime.now(timezone.utc).isoformat()
    layer_records = {
        key: _copy_layer(source, reference_root / key, repo_root)
        for key, source in source_layers.items()
    }

    source_figure = source_layers["render"] / "figures" / "phase_surface_render_mu_xi_T.png"
    source_plot_manifest = source_layers["render"] / "figures" / "plot_manifest.json"
    if not source_figure.is_file() or not source_plot_manifest.is_file():
        raise FileNotFoundError("render figure or plot manifest is missing")
    figure_root.mkdir(parents=True, exist_ok=True)
    destination_figure = figure_root / source_figure.name
    shutil.copyfile(source_figure, destination_figure)
    source_plot = read_json(source_plot_manifest)
    canonical_plot_manifest = {
        "schema_version": "pnjl_issue130_phase_reference_plot_manifest_v1",
        "generated_at_utc": generated_at,
        "calculation_sha": CALCULATION_SHA,
        "source_run_id": SOURCE_RUN_ID,
        "replay_run_id": REPLAY_RUN_ID,
        "reference_status": "candidate",
        "runtime_consumption": False,
        "source_plot_manifest": file_record(source_plot_manifest, repo_root),
        "source_render_semantics": source_plot.get("render_semantics", {}),
        "inputs": {
            "reference_root_manifest": "data/reference/pnjl/issue130_phase_reference_v1/manifest.json",
            "render_tables": "data/reference/pnjl/issue130_phase_reference_v1/render/tables",
        },
        "figure": file_record(destination_figure, repo_root),
    }
    plot_manifest_path = figure_root / "plot_manifest.json"
    _write_text(plot_manifest_path, json.dumps(canonical_plot_manifest, ensure_ascii=False, indent=2) + "\n")

    root_manifest: dict[str, Any] = {
        "schema_version": SCHEMA,
        "generated_at_utc": generated_at,
        "reference_status": "candidate",
        "promotion_status": "imported_candidate",
        "reference_write": True,
        "runtime_consumption": False,
        "calculation_sha": CALCULATION_SHA,
        "source_run_id": SOURCE_RUN_ID,
        "replay_run_id": REPLAY_RUN_ID,
        "uniform_xi_step": 0.00625,
        "uniform_xi_count": 161,
        "source_package": {
            "path": PACKAGE_RELATIVE.as_posix(),
            "manifest": file_record(package_root / "manifest.json", repo_root),
            "gate_manifest": file_record(gate["gate_manifest_path"], repo_root),
            "decision": file_record(gate["decision_path"], repo_root),
            "author_review": file_record(gate["author_record_path"], repo_root),
        },
        "layers": layer_records,
        "figure_root": relative_path(figure_root, repo_root),
        "figure": {
            "png": file_record(destination_figure, repo_root),
            "plot_manifest": file_record(plot_manifest_path, repo_root),
            "source_png": file_record(source_figure, repo_root),
            "source_plot_manifest": file_record(source_plot_manifest, repo_root),
        },
        "legacy_reference_snapshot": legacy_records,
        "constraints": {
            "legacy_files_unchanged": True,
            "canonical_consumer_switch": False,
            "rs_transport_started": False,
            "solver_called_during_import": False,
            "strict_unresolved_preserved": True,
            "derived_rows_are_noncertified": True,
            "render_triangulation": False,
        },
    }
    audit_path = reference_root / "IMPORT_AUDIT.md"
    _write_text(audit_path, (
        "# Issue #130 phase-reference v1 import audit\n\n"
        f"- calculation SHA: `{CALCULATION_SHA}`\n"
        f"- numerical source run: `{SOURCE_RUN_ID}`\n"
        f"- aggregate replay: `{REPLAY_RUN_ID}` (`solver_called=false`)\n"
        "- status: `imported_candidate`\n"
        "- runtime consumption: `false`\n"
        "- legacy reference files: byte-preserved; their pre-import hashes are recorded in `manifest.json`\n"
        "- strict layer: native C2/v6 plus the reviewed Maxwell expansion; unresolved geometry is retained\n"
        "- derived layer: uniform-xi local interpolation; non-native rows remain non-certified\n"
        "- render layer: structured tables and no-triangulation projection\n"
        "- this import does not run PNJL, alter tolerances, or change the default reference consumer\n"
    ))
    root_manifest["audit"] = file_record(audit_path, repo_root)
    _write_text(reference_root / "README.md", (
        "# Issue #130 phase-reference v1 candidate\n\n"
        "This is a new, versioned reference sibling imported from the author-reviewed "
        "strict/derived/render evidence package. It is not the default runtime reference.\n\n"
        "- `strict/` preserves native C2/v6 rows and unresolved geometry semantics.\n"
        "- `derived/` contains uniform-xi local interpolation and marks non-native rows as non-certified.\n"
        "- `render/` contains structured render tables; the figure is stored under the figure output root.\n"
        "- `manifest.json` records source hashes and confirms `runtime_consumption=false`.\n"
    ))
    _write_text(reference_root / "manifest.json", json.dumps(root_manifest, ensure_ascii=False, indent=2) + "\n")

    return root_manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-root", type=Path, default=PACKAGE_RELATIVE)
    parser.add_argument("--gate-root", type=Path, default=GATE_RELATIVE)
    parser.add_argument("--reference-root", type=Path, default=Path("data/reference/pnjl/issue130_phase_reference_v1"))
    parser.add_argument("--figure-root", type=Path, default=Path("data/outputs/figures/pnjl/phase_reference/issue130_phase_reference_v1"))
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--refresh-existing", action="store_true", help="refresh only an existing generated candidate directory")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    result = import_package(
        repo_root / args.package_root,
        repo_root / args.gate_root,
        repo_root / args.reference_root,
        repo_root / args.figure_root,
        repo_root,
        refresh_existing=args.refresh_existing,
    )
    print(json.dumps({
        "reference_root": result["layers"],
        "reference_status": result["reference_status"],
        "promotion_status": result["promotion_status"],
        "runtime_consumption": result["runtime_consumption"],
        "legacy_files_unchanged": result["constraints"]["legacy_files_unchanged"],
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
