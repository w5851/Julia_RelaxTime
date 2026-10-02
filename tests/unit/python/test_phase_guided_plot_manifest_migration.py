from scripts.analysis.relaxtime.migrate_phase_guided_plot_manifest_bundles import check_migration


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
