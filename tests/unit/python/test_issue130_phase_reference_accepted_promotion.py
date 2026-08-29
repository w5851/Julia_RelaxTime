from __future__ import annotations

import csv
import json
from pathlib import Path

from scripts.analysis.pnjl.build_issue130_phase_reference_v2 import build_package
from scripts.analysis.pnjl.promote_issue130_phase_reference_accepted import (
    ACCEPTED_SCOPE,
    ACCEPTED_STATUS,
    promote_package,
)
from tests.unit.python.test_issue130_phase_reference_v2 import _write_source


def test_author_promotion_updates_only_accepted_metadata(tmp_path: Path) -> None:
    source = tmp_path / "source"
    package = tmp_path / "package"
    _write_source(source)
    build_package(source, package)

    strict_before = (package / "strict" / "tables" / "maxwell_surface_strict_reference_v1.csv").read_bytes()
    result = promote_package(package, recorded_at="2026-08-29T00:00:00+00:00")

    assert result["promotion_status"] == "accepted_for_downstream"
    assert result["downstream_default_layer"] == "accepted"
    assert result["solver_called"] is False
    assert result["runtime_consumption"] is False
    assert result["reference_write"] is False
    assert all(value == 1 for value in result["changed_rows"].values())
    assert (package / "strict" / "tables" / "maxwell_surface_strict_reference_v1.csv").read_bytes() == strict_before

    accepted_path = package / "accepted" / "tables" / "maxwell_surface_accepted_phase_map_v1.csv"
    with accepted_path.open(encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle))
    assert row["acceptance_status"] == ACCEPTED_STATUS
    assert row["acceptance_scope"] == ACCEPTED_SCOPE
    assert row["source_status"] == "strict_certified"
    assert row["extrapolation"] == "False"

    root_manifest = json.loads((package / "manifest.json").read_text(encoding="utf-8"))
    accepted_manifest = json.loads((package / "accepted" / "manifest.json").read_text(encoding="utf-8"))
    assert root_manifest["promotion_status"] == "accepted_for_downstream"
    assert root_manifest["downstream_default_layer"] == "accepted"
    assert root_manifest["constraints"]["runtime_default_unchanged"] is True
    assert accepted_manifest["acceptance_status"] == ACCEPTED_STATUS
    assert accepted_manifest["promotion_status"] == "accepted_for_downstream"


def test_author_promotion_is_idempotent(tmp_path: Path) -> None:
    source = tmp_path / "source"
    package = tmp_path / "package"
    _write_source(source)
    build_package(source, package)
    promote_package(package, recorded_at="2026-08-29T00:00:00+00:00")
    second = promote_package(package, recorded_at="2026-08-29T00:00:00+00:00")
    assert all(value == 0 for value in second["changed_rows"].values())
