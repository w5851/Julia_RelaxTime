from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "pnjl" / "audit_issue130_phase_reference_runtime_consumers.py"
SPEC = importlib.util.spec_from_file_location("runtime_consumer_audit", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_candidate_tables_are_not_directly_compatible_with_legacy_contracts() -> None:
    assert MODULE.TABLE_SPECS
    assert all(spec["legacy"] for spec in MODULE.TABLE_SPECS.values())


def test_consumers_do_not_implicitly_select_candidate() -> None:
    assert all(item["candidate_reachable"] != "true" for item in MODULE.CONSUMERS)
    assert all(item["compatibility"] != "direct_compatible" for item in MODULE.CONSUMERS)


def test_candidate_inventory_rejects_duplicate_and_nonfinite_rows(tmp_path: Path) -> None:
    path = tmp_path / "table.csv"
    path.write_text("xi,T_MeV\n0,1\n0,1\n", encoding="utf-8")
    _, rows = MODULE.read_csv(path)
    assert MODULE.duplicate_count(rows, ("xi", "T_MeV")) == 1
    assert not MODULE.finite("NaN")
