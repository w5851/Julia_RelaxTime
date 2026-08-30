from __future__ import annotations

import json
import importlib.util
from pathlib import Path

ROOT = Path(__file__).parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "relaxtime" / "audit_issue130_rs_sharded_provenance.py"
SPEC = importlib.util.spec_from_file_location("issue130_rs_sharded_provenance", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_read_csv_ignores_provenance_comment_lines(tmp_path: Path) -> None:
    path = tmp_path / "scan.csv"
    path.write_text("# producer: pre-repair\n# hash: pending\na,b\n1,2\n", encoding="utf-8")
    fields, rows = MODULE.read_csv(path)
    assert fields == ["a", "b"]
    assert rows == [{"a": "1", "b": "2"}]


def test_effective_config_and_failed_points_mismatches_are_classified(tmp_path: Path) -> None:
    root = tmp_path / "results"
    root.mkdir()
    effective = root / "effective_config.json"
    failed = root / "failed_points.csv"
    scan = root / "phase_guided_transport_scan.csv"
    effective.write_text(json.dumps({"options": {"dry_run": False}}), encoding="utf-8")
    failed.write_text("T_MeV,muB_MeV,xi\n", encoding="utf-8")
    scan.write_text("T_MeV,muB_MeV,xi,converged\n1,2,0,true\n", encoding="utf-8")
    manifest = root / "run_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "artifacts": [
                    {"path": f"data/outputs/results/{name.name}", "sha256": "0" * 64}
                    for name in (effective, failed, scan)
                ]
            }
        ),
        encoding="utf-8",
    )
    rows = MODULE.artifact_hash_rows(manifest)
    by_name = {row["basename"]: row for row in rows}
    assert by_name["effective_config.json"]["verdict"] == "historical_known_mismatch"
    assert by_name["failed_points.csv"]["verdict"] == "historical_known_mismatch"
    assert by_name["phase_guided_transport_scan.csv"]["verdict"] == "unexpected_hash_mismatch"
    assert set(MODULE.ALLOWED_HISTORICAL_MISMATCHES) == {"effective_config.json", "failed_points.csv"}


def test_candidate_legacy_review_separates_classification_and_route_differences(tmp_path: Path) -> None:
    aggregate = tmp_path / "aggregate"
    output = tmp_path / "audit"
    aggregate.mkdir()
    legacy = tmp_path / "legacy.csv"
    fields = [
        "muB_MeV", "xi", "alpha_T", "phase_reference_kind", "phase_structure",
        "quality_flag", "quality_reason", "phase_prev", "phase_curr",
        "phase_boundary_xi_used", "seed_source", "equilibrium_backend",
    ]
    candidate_row = {
        "muB_MeV": "0", "xi": "0", "alpha_T": "1", "phase_reference_kind": "crossover",
        "phase_structure": "no_transition", "quality_flag": "false", "quality_reason": "ok",
        "phase_prev": "quark", "phase_curr": "quark", "phase_boundary_xi_used": "-0.1",
        "seed_source": "phase_aware_default_quark", "equilibrium_backend": "models",
    }
    legacy_row = {**candidate_row, "phase_structure": "crossover_possible", "phase_boundary_xi_used": "0.0"}
    for path, row in ((aggregate / "mode_a_scan.csv", candidate_row), (legacy, legacy_row)):
        path.write_text(",".join(fields) + "\n" + ",".join(row[field] for field in fields) + "\n", encoding="utf-8")
    (aggregate / "candidate_legacy_comparison.csv").write_text(
        "mode,key,column,legacy,candidate,abs_delta,rel_delta\n"
        "mode_a_fixed_muB_phase_scaled,0|0|1,tau_u,1.0,1.2,0.2,0.2\n",
        encoding="utf-8",
    )
    review = MODULE.candidate_legacy_review(
        aggregate,
        {
            "candidate_legacy_summary": {
                "mode_a_fixed_muB_phase_scaled": {"legacy_path": str(legacy)}
            }
        },
        output,
    )
    summary = review["phase_summary"]["mode_a_fixed_muB_phase_scaled"]
    assert summary["phase_classification_mismatch_keys"] == 1
    assert summary["route_difference_keys"] == 1
