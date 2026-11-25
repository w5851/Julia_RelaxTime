"""Regression guard for every vector package the repository currently claims as active.

The formal two-tier architecture delivers manuscript figures as `vector_delivery`
packages.  CI otherwise never touches the delivered artifacts themselves, so a
stray edit to a manifest, a lost PDF or a dangling AUTHOR review pin would only
surface when someone reruns the manual validator.  This test pins the active set
through the same entrypoint the SOP uses.
"""

from __future__ import annotations

import json
from pathlib import Path

from scripts.plotting.validate_plot_artifact import PROJECT_ROOT, validate_manifest

POINTER = PROJECT_ROOT / "docs/analysis/relaxtime/phase_guided_transport/publication_clean_current.json"


def _active_package_manifests() -> list[Path]:
    """Return every active vector package recorded in the current-layer pointer."""
    pointer = json.loads(POINTER.read_text(encoding="utf-8"))
    figures = pointer["current_manuscript_figures"]

    manifests = [PROJECT_ROOT / figures["figure_root"] / "plot_manifest.json"]
    for package in figures.get("additional_packages", []):
        manifests.append(PROJECT_ROOT / package["figure_root"] / "plot_manifest.json")
    for package in pointer.get("independent_active_vector_packages", []):
        manifests.append(PROJECT_ROOT / package["figure_root"] / "plot_manifest.json")
    return manifests


def test_pointer_declares_a_two_tier_formal_layer() -> None:
    pointer = json.loads(POINTER.read_text(encoding="utf-8"))

    assert pointer["status"] == "author_accepted_formal_layer"
    assert pointer["layer_architecture"]["publication_clean_layer_status"] == "abolished"
    assert pointer["current_numerical_layer"]["relative_velocity_convention"] == "medium_moller"
    assert pointer["current_numerical_layer"]["raw_manuscript_eligible"] is True
    # The retired cm layer must stay traceable instead of silently disappearing.
    assert pointer["retired_publication_clean_layers"]["v15_cm_layer"]["recoverable_from_git_commit"]


def test_active_vector_manifests_exist() -> None:
    manifests = _active_package_manifests()

    assert len(manifests) == 5
    for manifest in manifests:
        assert manifest.is_file(), f"missing active vector manifest: {manifest}"


def test_every_active_vector_package_validates() -> None:
    """Run the artifact validator over the whole active set with no errors."""
    failures = []
    for manifest in _active_package_manifests():
        errors = validate_manifest(manifest)
        if errors:
            failures.append((manifest.parent.name, errors))

    assert not failures, "\n".join(f"{name}: {errs}" for name, errs in failures)


def test_formal_layer_records_match_delivered_manifests() -> None:
    """`promotion.json` and the pointer must agree with the manifests on disk."""
    import hashlib

    pointer = json.loads(POINTER.read_text(encoding="utf-8"))
    promotion = json.loads(
        (PROJECT_ROOT / pointer["current_numerical_layer"]["promotion_record"]).read_text(encoding="utf-8")
    )

    def sha256(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    recorded = {item["plot_manifest"]: item["plot_manifest_sha256"] for item in promotion["manuscript_figure_layer"]}
    for rel, digest in recorded.items():
        assert sha256(PROJECT_ROOT / rel) == digest, f"promotion.json hash drifted for {rel}"

    assert sha256(PROJECT_ROOT / pointer["current_numerical_layer"]["promotion_record"]) == pointer[
        "current_numerical_layer"
    ]["promotion_record_sha256"]
