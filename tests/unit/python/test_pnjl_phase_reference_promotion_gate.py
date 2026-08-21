from pathlib import Path

from scripts.analysis.pnjl.audit_issue130_phase_reference_promotion import (
    audit_package,
    expected_xi_grid,
)


PACKAGE = Path(__file__).parents[3] / "docs" / "analysis" / "pnjl" / "phase_reference" / "issue130_phase_reference_layers_v1"


def test_expected_xi_grid_is_uniform_and_closed():
    values = expected_xi_grid()
    assert len(values) == 161
    assert values[0] == -0.5
    assert values[-1] == 0.5
    assert all(round(right - left, 8) == 0.00625 for left, right in zip(values, values[1:]))


def test_real_reviewed_package_is_promotion_candidate(tmp_path):
    decision = audit_package(
        PACKAGE,
        tmp_path / "gate",
        author_review_status="accepted",
        review_pr="248",
        review_merge_sha="09d1a6895100cef208c9108034d9f45b631158eb",
    )
    assert decision["verdict"] == "promotion_candidate"
    assert decision["reference_write"] is False
    assert decision["runtime_consumption"] is False
    assert decision["summary"]["maxwell_unresolved_rows"] > 0


def test_pending_author_review_blocks_gate(tmp_path):
    decision = audit_package(PACKAGE, tmp_path / "gate", author_review_status="pending")
    assert decision["verdict"] == "promotion_blocked"
