from __future__ import annotations

import json
from pathlib import Path

from scripts.plotting.inventory_figure_assets import (
    build_registry,
    write_candidates,
    write_registry,
)


def _write_png(path: Path, *, width: int = 2, height: int = 3) -> None:
    path.write_bytes(
        b"\x89PNG\r\n\x1a\n"
        + width.to_bytes(4, "big")
        + height.to_bytes(4, "big")
        + b"\x08\x02\x00\x00\x00"
        + b"IEND"
        + b"\x00\x00\x00\x00"
    )


def _write_contract_manifest(case_dir: Path) -> None:
    (case_dir / "plot_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": "plot_manifest_v1",
                "figure_mode": "strict",
                "style_profile": "strict_origin_like_v1",
                "semantic_status": "confirmed",
                "publication_scope": "main_text_candidate",
            }
        ),
        encoding="utf-8",
    )


def test_inventory_is_read_only_and_groups_only_variants(tmp_path):
    asset_root = tmp_path / "data" / "outputs" / "figures"
    contract_case = asset_root / "domain" / "family" / "case__plotv1__strict__line_first_v2"
    contract_case.mkdir(parents=True)
    _write_contract_manifest(contract_case)
    _write_png(contract_case / "figure.png")
    (contract_case / "figure.svg").write_text(
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 2 3"></svg>\n',
        encoding="utf-8",
    )

    variant_case = asset_root / "domain" / "family" / "case__plotv1__strict__line_first_v1"
    variant_case.mkdir(parents=True)
    _write_contract_manifest(variant_case)
    (variant_case / "figure.svg").write_text(
        '<svg xmlns="http://www.w3.org/2000/svg"></svg>\n',
        encoding="utf-8",
    )

    legacy_case = asset_root / "domain" / "family" / "legacy_case"
    legacy_case.mkdir(parents=True)
    (legacy_case / "figure.pdf").write_bytes(b"%PDF-1.7\n")

    registry = build_registry(
        asset_root=asset_root,
        repo_root=tmp_path,
        tracked_only=False,
        generated_at_utc="2026-08-16T00:00:00+00:00",
    )

    assert registry["schema_version"] == "figure_asset_registry_v1"
    assert registry["protection"]["delete_performed"] is False
    assert registry["protection"]["move_performed"] is False
    assert registry["discovery"]["asset_count"] == 4
    assert any(item["metadata"].get("width_px") == 2 for item in registry["assets"])
    assert any(item["metadata"].get("vector") is True for item in registry["assets"])
    assert any(item["metadata"].get("signature_valid") is True for item in registry["assets"])

    review_groups = registry["review_groups"]
    assert any(group["proposed_classifications"] == ["legacy_unregistered"] for group in review_groups)
    variant_group = next(group for group in review_groups if "case__plotv1__variants" in group["review_group"])
    assert variant_group["proposed_action"] == "owner_review_only"
    assert all(item["proposed_action"] != "delete" for item in registry["assets"])


def test_inventory_reports_can_be_written_without_overwriting(tmp_path):
    registry = {
        "schema_version": "figure_asset_registry_v1",
        "review_groups": [
            {
                "review_group": "domain/family/case",
                "case_directories": ["domain/family/case"],
                "asset_count": 1,
                "formats": ["png"],
                "manifest_schemas": [],
                "figure_modes": [],
                "style_profiles": [],
                "proposed_classifications": ["legacy_unregistered"],
                "proposed_action": "owner_review_only",
                "manual_questions": ["confirm external use"],
            }
        ],
    }
    registry_path = tmp_path / "asset_registry.json"
    candidates_path = tmp_path / "cleanup_candidates.csv"
    write_registry(registry_path, registry, overwrite=False)
    write_candidates(candidates_path, registry, overwrite=False)
    assert json.loads(registry_path.read_text(encoding="utf-8"))["schema_version"] == "figure_asset_registry_v1"
    assert "owner_review_only" in candidates_path.read_text(encoding="utf-8")

    try:
        write_registry(registry_path, registry, overwrite=False)
    except FileExistsError:
        pass
    else:
        raise AssertionError("inventory report overwrite must be explicit")
