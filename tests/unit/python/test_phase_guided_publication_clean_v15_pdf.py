import json

import pytest

from scripts.analysis.relaxtime import export_phase_guided_publication_clean_v15_pdf as export


def test_authorization_cannot_rebind_another_review_or_infer_author_approval(tmp_path):
    manifest, receipt = tmp_path / "manifest.json", tmp_path / "author.json"
    manifest.write_text("{}")
    value = dict(schema="transport_v15_author_acceptance_v1", png_accepted=True,
        pdf_export_authorized=True, figure_promotion_authorized=True,
        current_update_authorized=True, manuscript_use_authorized=True,
        author_instruction="Synthetic test fixture, not author approval.",
        review_manifest_sha256=export.sha256_file(manifest))
    receipt.write_text(json.dumps(value))
    assert export.checked_authorization(receipt, manifest) == value
    manifest.write_text("{\"changed\":true}")
    with pytest.raises(ValueError, match="different PNG bundle"):
        export.checked_authorization(receipt, manifest)
    value["png_accepted"] = False
    receipt.write_text(json.dumps(value))
    with pytest.raises(ValueError, match="explicit v15 author"):
        export.checked_authorization(receipt, manifest)


@pytest.mark.parametrize("case", ["../outside", "..", ".", "C:/outside", "parent\\child", ""])
def test_output_must_stay_in_single_new_case(case, tmp_path):
    with pytest.raises(ValueError, match="single new case"):
        export.export(case, tmp_path / "unused.json")


def test_font_relocation_requires_original_bytes(tmp_path):
    path = tmp_path / "retained.ttf"
    path.write_bytes(b"unrelated bytes")
    font = dict(path="old.ttf", bytes=10, sha256="0" * 64)
    pairs = [({}, dict(rendering=dict(font=font, quality=dict(actual_fonts=[font]))))]
    with pytest.raises(ValueError, match="accepted font bytes"):
        with export.exact_fonts(pairs, {"old.ttf": str(path)}):
            pytest.fail("unverified font used")
