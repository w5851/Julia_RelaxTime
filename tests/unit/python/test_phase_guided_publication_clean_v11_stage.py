from __future__ import annotations

import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[3] / "scripts/analysis/relaxtime/formalize_phase_guided_publication_clean_v11_stage.py"
SPEC = importlib.util.spec_from_file_location("v11_stage_tests", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_v11_stage_acceptance_is_current_complete_and_nonpromoting():
    record = MODULE.validate_acceptance()
    assert record["stage_accepted"] is True
    assert record["status"] == "author_accepted_stage_result"
    for field in ("manuscript_eligible", "current_publication_layer", "final_aps_delivery_complete",
                  "solver_called", "canonical_data_modified", "new_display_values", "raw_manuscript_eligible"):
        assert record[field] is False
    assert record["counts"]["png"] == record["counts"]["pdf"] == 74
    assert len(record["unresolved_formal_gates"]["native_pdf_submission_preflight"]) == 2
    assert len(record["renderer_dependency_closure"]) == 8
    assert {item["path"] for item in record["retained_v10_snapshot"]["expected_contract_drift"]} == MODULE.HISTORICAL_CONTRACT_PATHS
