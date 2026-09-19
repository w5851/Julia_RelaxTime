from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "relaxtime" / "formalize_phase_guided_publication_clean_v4.py"
SPEC = importlib.util.spec_from_file_location("phase_guided_publication_clean_formalization", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _source() -> dict[str, str]:
    return {
        "v3_package_sha256": "1" * 64,
        "v3_plot_sha256": "2" * 64,
        "v3_public_plot_sha256": "3" * 64,
        "v4_package_sha256_before": "4" * 64,
        "v4_plot_sha256_before": "5" * 64,
        "v4_public_plot_sha256_before": "6" * 64,
        "v4_layer_sha256_before": "7" * 64,
    }


def test_acceptance_metadata_keeps_numerical_layer_diagnostic() -> None:
    metadata = MODULE.accepted_fields("2026-09-19T00:00:00+00:00", _source())
    assert metadata["status"] == "derived_author_accepted_display_only"
    assert metadata["figure_status"] == "author_accepted_formal_layout"
    assert metadata["numerical_status"] == "diagnostic_only"
    assert metadata["manuscript_eligible"] is False
    assert metadata["current_publication_layer"] is True
    assert "publication_clean_v3" in metadata["supersedes"][0]


def test_formal_layer_readme_states_acceptance_boundary() -> None:
    text = MODULE.formal_layer_readme()
    assert "author-accepted current publication-clean layer" in text
    assert "figure_status=author_accepted_formal_layout" in text
    assert "numerical_status=diagnostic_only" in text
