#!/usr/bin/env python3
"""Record author acceptance of v11 as a stage result, not manuscript promotion."""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.plotting.plot_manifest import generator_record, input_record, sha256_file, write_manifest

EXPORT_SCRIPT = ROOT / "scripts/analysis/relaxtime/export_phase_guided_publication_clean_v11_pdf_review.py"
SPEC = importlib.util.spec_from_file_location("v11_stage_frozen_exporter", EXPORT_SCRIPT)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load frozen PDF exporter: {EXPORT_SCRIPT}")
EXPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EXPORT)
V11 = EXPORT.V11
BASE = V11.V10.TRANSPORT_ANALYSIS_ROOT
ACCEPTANCE = BASE / "publication_clean_v11_stage_acceptance_v1.json"
ARCHIVE_INDEX = BASE / "publication_review_history_v6_v9_archive_v1.json"
PLACEMENT = BASE / "phase_guided_transport_publication_clean_v11_paper_size_review/placement_report.json"
HISTORICAL_CONTRACT_PATHS = {
    "docs/guides/sop/workflows/figure_production.md",
    ".agents/skills/plotting-sop/SKILL.md",
}


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def audit_v10_snapshot() -> list[dict]:
    """Check immutable v10 bytes and explicitly bound its two stale contracts."""
    v11_package = read_json(EXPORT.PNG_PACKAGE)
    frozen = [record for record in v11_package["inputs"] if record["role"] == "retained_v10_review_artifact"]
    EXPORT.verify_records(frozen)
    actual_paths = {V11.relative(path) for root in (V11.V10.FIGURE_ROOT, V11.V10.ANALYSIS_ROOT)
                    for path in root.rglob("*") if path.is_file()}
    if {record["path"] for record in frozen} != actual_paths:
        raise ValueError("retained v10 snapshot file set changed")
    package = read_json(V11.V10.ANALYSIS_ROOT / "manifest.json")
    drift = []
    for record in package["inputs"]:
        path = ROOT / record["path"]
        actual = sha256_file(path)
        if actual != record["sha256"]:
            if record["path"] not in HISTORICAL_CONTRACT_PATHS:
                raise ValueError(f"unexpected v10 historical input drift: {record['path']}")
            drift.append({"path": record["path"], "recorded_sha256": record["sha256"], "current_sha256": actual})
    if {record["path"] for record in drift} != HISTORICAL_CONTRACT_PATHS:
        raise ValueError("historical v10 contract drift no longer matches the reviewed boundary")
    return drift


def collect_payload() -> dict:
    png_package, png_index, _ = EXPORT.load_frozen_case()
    pdf_package_path = EXPORT.ANALYSIS_ROOT / "manifest.json"
    pdf_package = read_json(pdf_package_path)
    EXPORT.verify_records([*pdf_package["inputs"], *pdf_package["outputs"]])
    pdf_index = read_json(EXPORT.FIGURE_ROOT / "pdf_review_index.json")
    if len(png_index["charts"]) != 74 or len(pdf_index["charts"]) != 74:
        raise ValueError("v11 requires 72 single and two composite charts per format")
    violations = []
    for chart in pdf_index["charts"]:
        manifest_path = ROOT / chart["manifest"]
        if sha256_file(manifest_path) != chart["manifest_sha256"]:
            raise ValueError(f"PDF companion manifest changed: {manifest_path}")
        if errors := EXPORT.verify_companion(manifest_path):
            raise ValueError(f"PDF technical preflight failed: {errors}")
        if chart["submission_preflight_violations"]:
            violations.append({"manifest": chart["manifest"], "violations": chart["submission_preflight_violations"]})
    if len(violations) != 2:
        raise ValueError("expected exactly two retained composite typography blockers")
    pointer = read_json(V11.CURRENT_POINTER)
    if not pointer["current_analysis_package"].endswith("publication_clean_v5"):
        raise ValueError("current pointer must remain at v5")
    placement = read_json(PLACEMENT)
    if placement["paper_project_modified"] or placement["manuscript_eligible"]:
        raise ValueError("paper-size review changed its evidence boundary")
    closure = [ROOT / f"scripts/analysis/relaxtime/build_phase_guided_publication_clean_v{version}.py" for version in range(6, 12)]
    closure.extend([EXPORT_SCRIPT, ROOT / "scripts/analysis/relaxtime/inspect_phase_guided_publication_clean_v11_placement.py"])
    return {
        "schema": "publication_clean_v11_stage_acceptance_v1", "status": "author_accepted_stage_result",
        "acceptance_date": "2026-10-02", "task_classification": "independent", "stage_accepted": True,
        "author_authorization": "v11 passed the current review and may be retained as a stage result; manage only v10/v11 necessary artifacts",
        "accepted_scope": "frozen v11 PNG layout and same-renderer PDF companions; display-only stage result",
        "manuscript_eligible": False, "current_publication_layer": False, "final_aps_delivery_complete": False,
        "solver_called": False, "canonical_data_modified": False, "new_display_values": False,
        "raw_manuscript_eligible": False, "numerical_status": "inherited_author_accepted_display_only",
        "calculation_sha": png_package["calculation_sha"], "workflow_head_sha": png_package["workflow_head_sha"],
        "png_package": input_record(EXPORT.PNG_PACKAGE, role="accepted_v11_png_package"),
        "pdf_package": input_record(pdf_package_path, role="accepted_v11_pdf_companions"),
        "png_index": input_record(EXPORT.PNG_INDEX, role="accepted_v11_png_index"),
        "pdf_index": input_record(EXPORT.FIGURE_ROOT / "pdf_review_index.json", role="accepted_v11_pdf_index"),
        "placement_report": input_record(PLACEMENT, role="measured_paper_placement_report"),
        "historical_archive_index": input_record(ARCHIVE_INDEX, role="local_restorable_history_index"),
        "renderer_dependency_closure": [input_record(path, role="retained_frozen_renderer_dependency") for path in closure],
        "current_pointer_retained": input_record(V11.CURRENT_POINTER, role="unchanged_current_v5_pointer"),
        "retained_v10_snapshot": {"status": "immutable_historical_parent_not_current_contract_pass",
                                 "expected_contract_drift": audit_v10_snapshot()},
        "counts": {"single_charts": 72, "composite_charts": 2, "png": 74, "pdf": 74, "mode_a_singles": 36, "mode_b_singles": 36},
        "unresolved_formal_gates": {"native_pdf_submission_preflight": violations,
                                    "paper_placement": placement["unresolved_formal_gates"],
                                    "numerical_qualification": "unchanged; no new convergence gate or production promotion"},
        "generator": generator_record(Path(__file__), command="python scripts/analysis/relaxtime/formalize_phase_guided_publication_clean_v11_stage.py --apply"),
    }


def validate_acceptance() -> dict:
    retained = read_json(ACCEPTANCE)
    actual = collect_payload()
    if retained != actual:
        raise ValueError("stage acceptance no longer matches the frozen package or its reviewed boundaries")
    return retained


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        record = validate_acceptance()
        print(f"[v11-stage] verified: {record['status']}; manuscript_eligible=false; current=v5")
        return
    if ACCEPTANCE.exists():
        raise FileExistsError("refusing to overwrite a stage acceptance record")
    payload = collect_payload()
    if args.apply:
        write_manifest(ACCEPTANCE, payload)
        print(f"[v11-stage] recorded: {ACCEPTANCE}")
    else:
        print("[v11-stage] preflight passed; use --apply to record stage-only acceptance")


if __name__ == "__main__":
    main()
