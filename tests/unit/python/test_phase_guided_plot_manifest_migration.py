from pathlib import Path

from scripts.analysis.relaxtime.migrate_phase_guided_plot_manifest_bundles import check_migration
from scripts.plotting.plot_provenance import git_source, validate_snapshot


def test_v10_v11_migration_is_lossless_and_preserves_accepted_evidence():
    report = check_migration()
    assert report["original_manifest_count"] == 225
    assert report["figure_record_count"] == 222
    assert len(report["protected"]) == 315
    assert len(report["bundles"]) == 3
    assert all(record["figure_count"] == 74 for record in report["bundles"])
    for field in ("solver_called", "figures_rendered", "canonical_data_modified",
                  "manuscript_eligibility_changed", "current_publication_layer_changed"):
        assert report[field] is False


def test_v11_snapshot_roots_do_not_require_retired_git_objects(monkeypatch):
    root = Path(__file__).resolve().parents[3]
    code_ref = "0725c0c4f87bbcbced4b7641c0a086e65c6fe42b"
    git_source.cache_clear()
    def unavailable(*args, **kwargs):
        raise FileNotFoundError("historical Git unavailable")
    monkeypatch.setattr("subprocess.check_output", unavailable)
    paths = [
        "docs/analysis/relaxtime/phase_guided_transport/publication_clean_v11_stage_acceptance_v1.json",
        "docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v11_png_review/manifest.json",
        "docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v11_pdf_review/manifest.json",
        "data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v11_png_review/plot_manifest.json",
        "data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v11_pdf_review/pdf_review_index.json",
    ]
    for path in paths:
        assert validate_snapshot(root / path, root=root, code_ref=code_ref) == []
    git_source.cache_clear()
