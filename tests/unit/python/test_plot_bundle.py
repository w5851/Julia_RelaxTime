from pathlib import Path

import pytest
from PIL import Image

from scripts.plotting.plot_bundle import build_bundle, expand_bundle
from scripts.plotting.plot_manifest import build_manifest, generator_record, input_record, output_record, write_manifest
from scripts.plotting.validate_plot_artifact import validate_manifest


def fixtures(root: Path):
    source = root / "draw.py"
    source.write_text("print('fixture')\n", encoding="utf-8")
    data = root / "points.csv"
    data.write_text("x,y\n0,1\n", encoding="utf-8")
    records = []
    charts = []
    for number in range(2):
        output = root / f"figure{number}.png"
        Image.new("RGB", (80, 60), "white").save(output, dpi=(600, 600))
        records.append(build_manifest(
            asset_id=f"fixture-{number}", figure_family="fixture", case_slug=f"figure{number}",
            figure_mode="audit", semantic_status="diagnostic", style_profile="audit_v1",
            publication_scope="internal_review",
            generator=generator_record(source, command="fixture", root=root),
            inputs=[input_record(data, role="data", root=root)],
            axes=[{"field": "x", "source_unit": "dimensionless", "display_unit": "dimensionless"}],
            series=[{"series_id": "line", "state": "supported", "support_rule": "all", "mask_rule": "none"}],
            outputs=[output_record(output, fmt="png", dpi=600, vector=False, root=root)],
            selection_rule="all", interpolation_policy="none", connector_policy="forbidden",
            missing_value_policy="preserve", validation={"finite": True, "duplicate_keys": True, "support": True},
            root=root))
        charts.append({"kind": "single", "stem": f"figure{number}", "outputs": records[-1]["outputs"]})
    return {"single_figure_count": 2, "composite_figure_count": 0, "charts": charts}, records


def test_shared_provenance_is_stored_once_and_roundtrips(tmp_path):
    index, records = fixtures(tmp_path)
    bundle = build_bundle(index, records)
    assert bundle["shared"]["inputs"] == records[0]["inputs"]
    assert all("inputs" not in figure["record"] for figure in bundle["figures"])
    assert [record for _, record in expand_bundle(bundle)] == records
    path = tmp_path / "plot_manifest.json"
    write_manifest(path, bundle)
    assert validate_manifest(path, repo_root=tmp_path) == []


@pytest.mark.parametrize("damage", ["duplicate_id", "shared_override", "empty", "invalid_package"])
def test_bundle_rejects_ambiguous_or_empty_records(tmp_path, damage):
    index, records = fixtures(tmp_path)
    bundle = build_bundle(index, records)
    if damage == "duplicate_id":
        bundle["figures"][1]["figure_id"] = bundle["figures"][0]["figure_id"]
    elif damage == "shared_override":
        bundle["figures"][0]["record"]["inputs"] = []
    elif damage == "empty":
        bundle["figures"] = []
    else:
        bundle["package"] = None
    with pytest.raises(ValueError):
        expand_bundle(bundle)


def test_bundle_checks_each_output_and_declared_counts(tmp_path):
    index, records = fixtures(tmp_path)
    bundle = build_bundle(index, records)
    bundle["package"]["single_figure_count"] = 3
    path = tmp_path / "plot_manifest.json"
    write_manifest(path, bundle)
    assert "bundle single_figure_count mismatch" in validate_manifest(path, repo_root=tmp_path)
    (tmp_path / "figure1.png").write_bytes(b"tampered")
    assert any("sha256 mismatch" in error for error in validate_manifest(path, repo_root=tmp_path))


def test_missing_per_figure_evidence_is_not_replaced_by_an_index(tmp_path):
    index, records = fixtures(tmp_path)
    bundle = build_bundle(index, records)
    del bundle["shared"]["axes"]
    path = tmp_path / "plot_manifest.json"
    write_manifest(path, bundle)
    assert any("manifest missing field: axes" in error for error in validate_manifest(path, repo_root=tmp_path))
