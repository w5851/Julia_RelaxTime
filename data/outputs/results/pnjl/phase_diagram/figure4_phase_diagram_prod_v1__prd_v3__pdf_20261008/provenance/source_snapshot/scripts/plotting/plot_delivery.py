"""Verify the author-accepted PNG behind a subsequent vector delivery.

Acceptance is supplied by the caller from an actual author response. This
module verifies its bound artifacts; it never creates or infers approval.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from scripts.plotting.plot_provenance import validate_hash_record


SEMANTIC_FIELDS = ("style_profile", "axes", "series", "selection_rule", "interpolation_policy",
                   "connector_policy", "missing_value_policy", "derived_display_geometry",
                   "caption_parameters", "frozen_render_signature")
RENDER_FIELDS = ("figure_size_inches", "column", "legend_policy", "legend_placements",
                 "parameter_encoding", "coexistence_fill", "accessibility")


def _checked_json(record: Any, root: Path, label: str, errors: list[str]) -> dict[str, Any] | None:
    if not isinstance(record, dict):
        errors.append(f"{label} must be a hash-bound file record")
        return None
    issues = validate_hash_record(record, root=root, label=label, allow_historical=False)
    errors.extend(issues)
    if issues:
        return None
    path = Path(record["path"])
    try:
        value = json.loads((path if path.is_absolute() else root / path).read_text(encoding="utf-8-sig"))
        if not isinstance(value, dict):
            raise ValueError("JSON root is not an object")
        return value
    except (OSError, ValueError) as exc:
        errors.append(f"cannot read {label}: {exc}")
        return None


def validate_png_delivery(manifest: dict[str, Any], *, root: Path) -> list[str]:
    """Check acceptance, PNG identity, inputs and display semantics across stages."""
    errors: list[str] = []
    review = manifest.get("author_review")
    if not isinstance(review, dict) or review.get("status") != "author_accepted_png":
        return ["vector delivery requires explicit author_accepted_png evidence"]
    acceptance = _checked_json(review.get("record"), root, "author_review.record", errors)
    if acceptance is None:
        return errors
    if (acceptance.get("schema") != "plot_png_acceptance_v1"
            or acceptance.get("status") != "author_accepted"
            or acceptance.get("pdf_export_authorized") is not True
            or not isinstance(acceptance.get("author_instruction"), str)
            or not acceptance["author_instruction"].strip()):
        errors.append("PNG acceptance must record an actual author instruction and PDF authorization")
    accepted = _checked_json(acceptance.get("accepted_manifest"), root, "accepted_manifest", errors)
    png = acceptance.get("accepted_png")
    if not isinstance(png, dict):
        errors.append("acceptance must bind the accepted PNG hash")
        return errors
    errors.extend(validate_hash_record(png, root=root, label="accepted_png", allow_historical=False))
    if accepted is None:
        return errors
    rendering = accepted.get("rendering", {})
    if (rendering.get("delivery_stage") != "png_review"
            or rendering.get("vector_delivery_pending") is not True
            or accepted.get("manuscript_eligible") is not False
            or accepted.get("current_publication_layer") is not False):
        errors.append("accepted manifest must be the frozen non-eligible PNG review stage")
    accepted_outputs = accepted.get("outputs", [])
    if not accepted_outputs or any(item.get("format") != "png" for item in accepted_outputs):
        errors.append("accepted review outputs must contain PNG only")
    def same_png(item: dict[str, Any]) -> bool:
        return (item.get("format") == "png" and item.get("sha256") == png.get("sha256")
                and item.get("bytes") == png.get("bytes"))
    if not any(same_png(item) for item in accepted_outputs):
        errors.append("accepted PNG is not bound by the accepted manifest")
    if not any(same_png(item) for item in manifest.get("outputs", [])):
        errors.append("vector delivery PNG differs from the author-accepted PNG")
    if manifest.get("rendering", {}).get("vector_delivery_pending") is not False:
        errors.append("vector delivery must record vector_delivery_pending=false")
    for field in SEMANTIC_FIELDS:
        if field not in accepted or field not in manifest or accepted[field] != manifest[field]:
            errors.append(f"PNG/vector semantic mismatch: {field}")
    for field in RENDER_FIELDS:
        if field not in rendering or rendering[field] != manifest.get("rendering", {}).get(field):
            errors.append(f"PNG/vector rendering mismatch: {field}")
    if accepted.get("generator", {}).get("sha256") != manifest.get("generator", {}).get("sha256"):
        errors.append("PNG/vector generator hash mismatch")
    # Snapshot locations may change on delivery; content and roles may not.
    def frozen_inputs(value: dict[str, Any]) -> list[tuple[str, str, int]]:
        return sorted((item.get("role", ""), item.get("sha256", ""), item.get("bytes", -1))
                      for item in value.get("inputs", [])
                      if item.get("role") not in {"author_png_acceptance", "accepted_png_manifest"})
    if frozen_inputs(accepted) != frozen_inputs(manifest):
        errors.append("PNG/vector frozen input contents or roles changed")
    return errors
