from __future__ import annotations

import json

import pytest

from scripts.analysis.relaxtime import build_phase_guided_publication_clean_v15 as delivery


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")


def candidate(tmp_path, **overrides):
    review = dict(calculation_sha=delivery.raw.CALCULATION_SHA, points=1867,
                  full_grid_numerical_acceptance_passed=True, artifact_contract_passed=True)
    review.update(overrides)
    path = tmp_path / "returned/grid/review.json"
    write(path, review)
    write(tmp_path / "formal_promotion_v1/promotion.json", dict(formal_production_promoted=True,
          grid_review_sha256=delivery.sha256_file(path)))
    return path


@pytest.mark.parametrize("override", [dict(points=1819), dict(calculation_sha="unrelated"),
                                    dict(full_grid_numerical_acceptance_passed=False), dict(artifact_contract_passed=False)])
def test_loader_rejects_wrong_or_unaccepted_grid_before_rendering(tmp_path, override):
    candidate(tmp_path, **override)
    with pytest.raises(ValueError, match="numerical grid is not accepted"):
        delivery.load_inputs(tmp_path)


def test_promotion_must_bind_current_grid_bytes(tmp_path):
    path = candidate(tmp_path)
    path.write_text(path.read_text() + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="promotion does not bind"):
        delivery.load_inputs(tmp_path)


def test_package_must_match_verified_return(tmp_path):
    candidate(tmp_path)
    write(tmp_path / "returned/package_manifest.json", dict(files=[]))
    write(tmp_path / "local_return_verification.json", dict(status="verified_return", package_manifest_sha256="wrong"))
    with pytest.raises(ValueError, match="integrity verification"):
        delivery.load_inputs(tmp_path)


def test_existing_figure_case_is_preserved(tmp_path, monkeypatch):
    monkeypatch.setattr(delivery.raw, "FIGURE_PARENT", tmp_path / "figures")
    monkeypatch.setattr(delivery.raw, "RESULT_PARENT", tmp_path / "results")
    destination = tmp_path / "figures/existing"
    destination.mkdir(parents=True)
    marker = destination / "previous.png"
    marker.write_bytes(b"frozen previous figure")
    with pytest.raises(FileExistsError):
        delivery.build_review("existing", tmp_path / "unused_source")
    assert marker.read_bytes() == b"frozen previous figure"
    assert not (tmp_path / "results/existing").exists()
