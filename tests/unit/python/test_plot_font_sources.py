"""A relocated font must retain exact historical bytes and its measured family."""
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
import shutil

from scripts.plotting.plot_style import load_profile, resolve_font
from scripts.plotting.validate_plot_artifact import _check_font_evidence


def sample(tmp_path):
    profile = load_profile("candidate_aps_v2")
    data = deepcopy(profile.data)
    data.update(font_policy="error", font_candidates=["DejaVu Serif"])
    profile = replace(profile, data=data)
    font = resolve_font(profile)
    source = Path(font["path"])
    retained = tmp_path / "retained.ttf"
    shutil.copyfile(source, retained)
    missing = str(tmp_path / "old_environment" / source.name)
    font["path"] = missing
    rendering = dict(font=font, quality=dict(font=deepcopy(font), actual_fonts=[deepcopy(font)]))
    return profile, rendering, {missing: str(retained)}, retained


def test_exact_relocation_is_explicit_and_does_not_modify_record(tmp_path):
    profile, rendering, mapping, _ = sample(tmp_path)
    before = deepcopy(rendering)
    errors = []
    _check_font_evidence(rendering, profile, tmp_path, errors)
    assert errors and all("missing file" in issue for issue in errors)
    errors = []
    _check_font_evidence(rendering, profile, tmp_path, errors, font_sources=mapping)
    assert errors == []
    assert rendering == before


def test_relocation_cannot_replace_a_font_with_different_bytes(tmp_path):
    profile, rendering, mapping, retained = sample(tmp_path)
    retained.write_bytes(b"wrong font")
    errors = []
    _check_font_evidence(rendering, profile, tmp_path, errors, font_sources=mapping)
    assert any("sha256 mismatch" in issue for issue in errors)


def test_exact_bytes_do_not_override_an_incorrect_family_claim(tmp_path):
    profile, rendering, mapping, _ = sample(tmp_path)
    rendering["font"]["resolved"] = "Unrelated Family"
    rendering["quality"]["font"]["resolved"] = "Unrelated Family"
    rendering["quality"]["actual_fonts"][0]["resolved"] = "Unrelated Family"
    errors = []
    _check_font_evidence(rendering, profile, tmp_path, errors, font_sources=mapping)
    assert any("font name differs" in issue for issue in errors)
