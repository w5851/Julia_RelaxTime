from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "relaxtime" / "formalize_phase_guided_publication_clean_v5.py"
SPEC = importlib.util.spec_from_file_location("phase_guided_publication_clean_v5_formalization", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_v5_acceptance_metadata_is_display_only() -> None:
    metadata = MODULE.accepted_fields("2026-09-24T00:00:00+00:00")
    assert metadata["status"] == "derived_author_accepted_display_only"
    assert metadata["figure_status"] == "author_accepted_formal_layout"
    assert metadata["numerical_status"] == "author_accepted_display_only"
    assert metadata["manuscript_eligible"] is True
    assert metadata["manuscript_eligibility_scope"] == "publication_clean_v5_display_layer_only"
    assert metadata["raw_numerical_status"] == "diagnostic_only"
    assert metadata["raw_manuscript_eligible"] is False
    assert metadata["convergence_gate_status"] == "local_high_rate_gate_not_run"
    assert metadata["author_acceptance"]["source"] == "formal_manuscript_uses_publication_clean_v5"


def test_v5_label_contract_is_exact() -> None:
    assert MODULE.EXPECTED_LABELS[("relaxation_time_y_axis", "tau_u")] == (
        r"$\tau_u$",
        r"$\tau_u\;[\mathrm{fm}]$",
    )
    assert MODULE.EXPECTED_LABELS[("first_order_endpoint_legend", "quark")][1] == (
        "chirally restored branch endpoint"
    )
    assert MODULE.EXPECTED_LABELS[("first_order_endpoint_legend", "hadron")][1] == (
        "chirally broken branch endpoint"
    )
