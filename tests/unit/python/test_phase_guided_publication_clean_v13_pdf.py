from __future__ import annotations

import copy
import json
import zipfile

import matplotlib.pyplot as plt
import pytest
from PIL import Image, PngImagePlugin

from scripts.analysis.relaxtime import export_phase_guided_publication_clean_v13_pdf as pdf
from scripts.plotting.plot_bundle import load_chart_records
from scripts.plotting.plot_manifest import input_record
from scripts.plotting.plot_style import configure_matplotlib, load_profile


@pytest.mark.parametrize("field,value", [
    ("author_png_accepted", False), ("pdf_export_authorized", None),
    ("manuscript_eligible", True), ("current_publication_layer", True),
])
def test_export_requires_acceptance_without_eligibility_promotion(field, value):
    acceptance = pdf.read_json(pdf.ACCEPTANCE)
    pdf.require_acceptance(acceptance)
    acceptance[field] = value
    with pytest.raises(ValueError, match="explicit PNG acceptance"):
        pdf.require_acceptance(acceptance)


def test_live_renderer_guard_rejects_drift_even_with_an_available_snapshot(tmp_path):
    source = tmp_path / "renderer.py"
    source.write_bytes(b"accepted source\n")
    frozen = input_record(source, role="renderer", root=tmp_path)
    snapshot = tmp_path / "retained.py"
    snapshot.write_bytes(source.read_bytes())
    frozen["source_snapshot"] = input_record(snapshot, role="source_snapshot", root=tmp_path)
    pdf.require_current(frozen, root=tmp_path)
    source.write_bytes(b"changed source\n")
    with pytest.raises(ValueError, match="frozen live file changed"):
        pdf.require_current(frozen, root=tmp_path)


def test_historical_snapshot_checks_do_not_authorize_a_drifted_live_renderer(monkeypatch, tmp_path):
    monkeypatch.setattr(pdf, "ROOT", tmp_path)
    require_current = pdf.require_current
    monkeypatch.setattr(pdf, "require_current", lambda record: require_current(record, root=tmp_path))
    source = tmp_path / "scripts/renderer.py"
    source.parent.mkdir()
    source.write_bytes(b"accepted source\n")
    frozen = input_record(source, role="renderer", root=tmp_path)
    archive_path = tmp_path / "snapshot.zip"
    package_hash = "a" * 64
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr("snapshot_manifest.json", json.dumps({
            "source_package_sha256": package_hash, "files": [frozen]}))
        archive.writestr(frozen["path"], source.read_bytes())
    acceptance = {"source_snapshot": input_record(archive_path, role="source_snapshot", root=tmp_path),
                  "png_package": {"sha256": package_hash}}
    package = {"inputs": [], "generator": frozen}
    source.write_bytes(b"changed source\n")
    assert pdf.verify_source_snapshot(acceptance, package, require_live=False) == ["scripts/renderer.py"]
    with pytest.raises(ValueError, match="frozen live file changed"):
        pdf.verify_source_snapshot(acceptance, package)
    archive_path.write_bytes(b"corrupt retained archive")
    with pytest.raises(ValueError, match="frozen live file changed"):
        pdf.verify_source_snapshot(acceptance, package, require_live=False)


def test_pixel_guard_ignores_png_metadata_but_rejects_changed_curve(tmp_path):
    configure_matplotlib(load_profile("candidate_aps_v2"))
    figure, axis = plt.subplots(figsize=(2, 2), dpi=72)
    try:
        line, = axis.plot([0, 1], [0, 1])
        reference = tmp_path / "accepted.png"
        figure.savefig(reference, dpi=72, bbox_inches=None)
        with Image.open(reference) as image:
            retained = image.copy()
        metadata = PngImagePlugin.PngInfo()
        metadata.add_text("comment", "metadata is not part of visual acceptance")
        retained.save(reference, pnginfo=metadata)
        assert pdf.confirm_png_pixels(figure, reference, 72)["pixel_identical"] is True
        line.set_ydata([1, 0])
        with pytest.raises(ValueError, match="no longer reproduces"):
            pdf.confirm_png_pixels(figure, reference, 72)
    finally:
        plt.close(figure)


def test_specs_or_legend_drift_fails_before_pdf_export(monkeypatch):
    _, accepted = load_chart_records(pdf.PNG_INDEX, root=pdf.ROOT)
    chart, source = accepted[0]
    figure = plt.figure()
    changed = copy.deepcopy(source["rendering"]["panel_specs"])
    changed[0]["axis_scale"] = "changed"
    monkeypatch.setattr(pdf.V13, "render_single", lambda *args: (
        figure, changed, source["rendering"]["legend_placements"]))
    with pytest.raises(ValueError, match="specifications or legend placements changed"):
        pdf.export_chart(chart, source, {}, {}, load_profile("candidate_aps_v2"), {})
    assert not plt.fignum_exists(figure.number)


def test_export_refuses_to_overwrite_existing_case(monkeypatch, tmp_path):
    monkeypatch.setattr(pdf, "FIGURE_ROOT", tmp_path)
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        pdf.build_delivery()


def test_artifact_contains_all_accepted_views_as_vector_delivery():
    pdf.check_delivery()
    _, pairs = load_chart_records(pdf.FIGURE_ROOT / "plot_manifest.json", root=pdf.ROOT)
    package = json.loads((pdf.ANALYSIS_ROOT / "manifest.json").read_text(encoding="utf-8"))
    assert package["validation_summary"]["pixel_identical_count"] == 75
    assert package["point_row_count"] == 21828
    for _, record in pairs:
        assert record["manuscript_eligible"] is False
        assert record["current_publication_layer"] is False
        assert record["solver_called"] is False
        assert record["canonical_data_modified"] is False
        assert record["new_display_values"] is False
        assert record["rendering"]["typography_exception"] is None
        assert record["rendering"]["quality"]["minimum_capital_numeral_height_mm"] >= 2
        assert record["rendering"]["pdf_placement_limits"]["maximum_width_inches"] == 7
        assert record["rendering"]["placement_limits"]["maximum_width_inches"] == 6.75
        assert record["rendering"]["pdf_placement_limits"]["single_column_reuse_qualified"] is False
        output = record["outputs"][2]
        assert output["inspection"]["page_count"] == 1
        assert output["inspection"]["raster_image_count"] == 0
        assert output["inspection"]["fonts_embedded"] is True
        assert output["inspection"]["type3_font_count"] == 0
        assert output["inspection"]["physical_size_inches"] == pytest.approx(record["rendering"]["figure_size_inches"])


def test_artifact_scientific_scope_and_png_references_cannot_be_changed():
    _, accepted = load_chart_records(pdf.PNG_INDEX, root=pdf.ROOT)
    _, delivered = load_chart_records(pdf.FIGURE_ROOT / "plot_manifest.json", root=pdf.ROOT)
    source = accepted[0][1]
    candidate = copy.deepcopy(delivered[0][1])
    candidate["axes"][0]["display_unit"] = "MeV"
    with pytest.raises(ValueError, match="accepted source field changed: axes"):
        pdf.verify_derivative(source, candidate, load_profile("candidate_aps_v2"))
    candidate = copy.deepcopy(delivered[0][1])
    candidate["outputs"][0]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="retain both accepted PNG references"):
        pdf.verify_derivative(source, candidate, load_profile("candidate_aps_v2"))
