from pathlib import Path

from scripts.analysis.pnjl.import_issue130_phase_reference import CALCULATION_SHA, import_package


ROOT = Path(__file__).parents[3]
PACKAGE = ROOT / "docs" / "analysis" / "pnjl" / "phase_reference" / "issue130_phase_reference_layers_v1"
GATE = ROOT / "docs" / "analysis" / "pnjl" / "phase_reference" / "issue130_phase_reference_promotion_gate_v1"


def test_import_creates_candidate_sibling_without_requiring_deleted_legacy_reference(tmp_path):
    reference_root = tmp_path / "reference" / "issue130_phase_reference_v1"
    figure_root = tmp_path / "figures" / "issue130_phase_reference_v1"

    manifest = import_package(PACKAGE, GATE, reference_root, figure_root, ROOT)

    assert manifest["schema_version"] == "pnjl_issue130_phase_reference_import_v1"
    assert manifest["calculation_sha"] == CALCULATION_SHA
    assert manifest["reference_status"] == "candidate"
    assert manifest["runtime_consumption"] is False
    assert manifest["constraints"]["strict_unresolved_preserved"] is True
    assert (reference_root / "IMPORT_AUDIT.md").is_file()
    assert (reference_root / "strict" / "tables" / "maxwell_surface_strict_reference_v1.csv").is_file()
    assert (reference_root / "derived" / "tables" / "surface_coverage_mask.csv").is_file()
    assert (reference_root / "render" / "tables" / "maxwell_surface_render.csv").is_file()
    assert (figure_root / "phase_surface_render_mu_xi_T.png").is_file()
    assert (figure_root / "plot_manifest.json").is_file()
    assert len(manifest["legacy_reference_snapshot"]) == 6
    assert all(item["availability"] == "git_recovery_ref_only" for item in manifest["legacy_reference_snapshot"])


def test_import_refuses_existing_destination(tmp_path):
    reference_root = tmp_path / "reference"
    figure_root = tmp_path / "figures"
    reference_root.mkdir()
    try:
        import_package(PACKAGE, GATE, reference_root, figure_root, ROOT)
    except FileExistsError:
        pass
    else:
        raise AssertionError("existing reference root must not be overwritten")


def test_refresh_existing_is_explicit(tmp_path):
    reference_root = tmp_path / "reference"
    figure_root = tmp_path / "figures"
    first = import_package(PACKAGE, GATE, reference_root, figure_root, ROOT)
    refreshed = import_package(PACKAGE, GATE, reference_root, figure_root, ROOT, refresh_existing=True)
    assert first["calculation_sha"] == refreshed["calculation_sha"] == CALCULATION_SHA
    assert (reference_root / "manifest.json").is_file()
