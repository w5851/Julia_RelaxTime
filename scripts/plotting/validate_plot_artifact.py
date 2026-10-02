"""Validate the machine-readable contract of a new plot artifact."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
from typing import Any

_SCRIPT_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_SCRIPT_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_PROJECT_ROOT))

from scripts.plotting.plot_manifest import MANIFEST_SCHEMA, PROJECT_ROOT
from scripts.plotting.plot_style import ALLOWED_PROFILES, load_profile
from scripts.plotting.plot_quality import inspect_export
from scripts.plotting.plot_provenance import code_ref_for_manifest, validate_hash_record, validate_snapshot


ALLOWED_MODES = {"audit", "estimated_midpoint", "strict", "legacy"}
STRICT_FORBIDDEN_STATES = {
    "unresolved",
    "nonconverged",
    "estimated_midpoint",
    "cep_bracket",
    "plotting_connector",
}
REQUIRED_FIELDS = {
    "schema_version",
    "asset_id",
    "figure_family",
    "case_slug",
    "figure_mode",
    "semantic_status",
    "style_profile",
    "publication_scope",
    "generator",
    "inputs",
    "axes",
    "series",
    "outputs",
    "validation",
}


def _resolve_artifact_path(value: str, root: Path) -> Path:
    candidate = Path(value)
    return candidate.resolve() if candidate.is_absolute() else (root / candidate).resolve()


def _check_hash_record(record: dict[str, Any], *, root: Path, label: str, errors: list[str], code_ref: str | None = None) -> Path | None:
    issues = validate_hash_record(record, root=root, label=label, code_ref=code_ref, allow_historical=not label.startswith("outputs"))
    errors.extend(issues)
    value = record.get("path")
    return _resolve_artifact_path(value, root) if isinstance(value, str) and value and not issues else None


def _check_png_dpi(path: Path, declared_dpi: Any, *, label: str, errors: list[str]) -> None:
    if declared_dpi is None:
        errors.append(f"{label}.dpi is required for PNG")
        return
    try:
        from PIL import Image

        with Image.open(path) as image:
            actual = image.info.get("dpi")
            actual_value = min(actual) if isinstance(actual, tuple) else actual
            if actual_value is None or float(actual_value) + 1.0 < float(declared_dpi):
                errors.append(f"{label} actual PNG dpi {actual_value!r} is below declared {declared_dpi}")
    except ImportError:
        errors.append("Pillow is required to validate PNG DPI")
    except Exception as exc:
        errors.append(f"{label} PNG inspection failed: {exc}")


def validate_manifest(manifest_path: str | Path, *, repo_root: Path = PROJECT_ROOT, code_ref: str | None = None) -> list[str]:
    """Return all contract violations; an empty list means the artifact passes."""

    path = Path(manifest_path).resolve()
    errors: list[str] = []
    profile = None
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return [f"cannot read manifest {path}: {exc}"]
    if not isinstance(manifest, dict):
        return ["manifest root must be a JSON object"]

    missing = sorted(REQUIRED_FIELDS - set(manifest))
    errors.extend(f"manifest missing field: {field}" for field in missing)
    if manifest.get("schema_version") != MANIFEST_SCHEMA:
        errors.append(f"schema_version must be {MANIFEST_SCHEMA!r}")

    mode = manifest.get("figure_mode")
    if mode not in ALLOWED_MODES:
        errors.append(f"invalid figure_mode: {mode!r}")
    style = manifest.get("style_profile")
    if style not in ALLOWED_PROFILES:
        errors.append(f"invalid style_profile: {style!r}")
    else:
        try:
            profile = load_profile(style)
        except (OSError, ValueError) as exc:
            profile = None
            errors.append(str(exc))

    if mode == "estimated_midpoint":
        if manifest.get("publication_scope") != "supplement_or_internal_review":
            errors.append("estimated_midpoint must use publication_scope=supplement_or_internal_review")
    if mode == "strict":
        if manifest.get("interpolation_policy") != "none":
            errors.append("strict interpolation_policy must be none")
        if manifest.get("connector_policy") != "forbidden":
            errors.append("strict connector_policy must be forbidden")

    generator = manifest.get("generator")
    historical_code_ref = (code_ref or code_ref_for_manifest(path, repo_root)
                           or manifest.get("postprocess_sha") or manifest.get("git_commit"))
    if isinstance(generator, dict):
        _check_hash_record(generator, root=repo_root, label="generator", errors=errors, code_ref=historical_code_ref)
    else:
        errors.append("generator must be an object")

    inputs = manifest.get("inputs")
    if not isinstance(inputs, list) or not inputs:
        errors.append("inputs must be a non-empty list")
    else:
        for index, record in enumerate(inputs):
            label = f"inputs[{index}]"
            if not isinstance(record, dict):
                errors.append(f"{label} must be an object")
                continue
            if not record.get("role"):
                errors.append(f"{label}.role is required")
            _check_hash_record(record, root=repo_root, label=label, errors=errors, code_ref=historical_code_ref)

    axes = manifest.get("axes")
    if not isinstance(axes, list) or not axes:
        errors.append("axes must be a non-empty list")
    else:
        for index, axis in enumerate(axes):
            label = f"axes[{index}]"
            if not isinstance(axis, dict):
                errors.append(f"{label} must be an object")
                continue
            for field in ("field", "source_unit", "display_unit"):
                if not axis.get(field):
                    errors.append(f"{label}.{field} is required")

    series = manifest.get("series")
    if not isinstance(series, list) or not series:
        errors.append("series must be a non-empty list")
    else:
        seen_ids: set[str] = set()
        for index, item in enumerate(series):
            label = f"series[{index}]"
            if not isinstance(item, dict):
                errors.append(f"{label} must be an object")
                continue
            series_id = item.get("series_id")
            if not isinstance(series_id, str) or not series_id:
                errors.append(f"{label}.series_id is required")
            elif series_id in seen_ids:
                errors.append(f"duplicate series_id: {series_id}")
            else:
                seen_ids.add(series_id)
            for field in ("state", "support_rule", "mask_rule"):
                if not item.get(field):
                    errors.append(f"{label}.{field} is required")
            state = item.get("state")
            if mode == "strict" and state in STRICT_FORBIDDEN_STATES:
                errors.append(f"strict series[{index}] contains forbidden state: {state}")

    outputs = manifest.get("outputs")
    formats: set[str] = set()
    if not isinstance(outputs, list) or not outputs:
        errors.append("outputs must be a non-empty list")
    else:
        for index, record in enumerate(outputs):
            label = f"outputs[{index}]"
            if not isinstance(record, dict):
                errors.append(f"{label} must be an object")
                continue
            fmt = str(record.get("format", "")).lower()
            if not fmt:
                errors.append(f"{label}.format is required")
            formats.add(fmt)
            output_path = _check_hash_record(record, root=repo_root, label=label, errors=errors)
            if fmt == "svg" and record.get("vector") is not True:
                errors.append(f"{label}.vector must be true for SVG")
            if fmt == "png" and output_path is not None:
                _check_png_dpi(output_path, record.get("dpi"), label=label, errors=errors)
    if mode == "strict":
        missing_formats = set(profile.formats if profile is not None else ("png", "svg")) - formats
        if missing_formats:
            errors.append(f"strict outputs missing required formats: {sorted(missing_formats)}")
        png_records = [record for record in outputs if isinstance(record, dict) and str(record.get("format", "")).lower() == "png"]
        if not any(record.get("dpi") is not None and int(record["dpi"]) >= 600 for record in png_records):
            errors.append("strict requires a PNG output declared at >=600 dpi")

    validation = manifest.get("validation")
    if not isinstance(validation, dict):
        errors.append("validation must be an object")
    else:
        for field in ("finite", "duplicate_keys", "support"):
            if field not in validation:
                errors.append(f"validation.{field} is required")
        if mode == "strict":
            for field in ("finite", "duplicate_keys", "support", "strict_gate"):
                if validation.get(field) is not True:
                    errors.append(f"strict validation.{field} must be true")

    rendering = manifest.get("rendering")
    if mode == "strict":
        if not isinstance(rendering, dict):
            errors.append("strict rendering metadata is required")
        else:
            column = rendering.get("column")
            declared_size = rendering.get("figure_size_inches")
            if column not in {"single_column", "double_column"}:
                errors.append("strict rendering.column must be single_column or double_column")
            elif not isinstance(declared_size, list) or len(declared_size) != 2:
                errors.append("strict rendering.figure_size_inches must be [width, height]")
            else:
                expected_size = profile.data["figure_size_in"][column] if profile is not None else None
                if rendering.get("size_override_reason") and profile is not None and profile.data.get("quality", {}).get("contract") == "final_size_v2":
                    expected_size = declared_size
                if expected_size is not None and any(abs(float(a) - float(b)) > 1.0e-9 for a, b in zip(declared_size, expected_size)):
                    errors.append(f"strict rendering size {declared_size} does not match profile {expected_size}")

    if mode == "strict" and any(token in str(manifest.get("semantic_status", "")).lower() for token in ("unresolved", "estimated", "bracket")):
        errors.append("strict semantic_status cannot claim unresolved, estimated, or bracket content")

    if profile is not None and mode == "strict" and profile.data["semantics"].get("allow_connector") is not False:
        errors.append("strict profile must disallow connector rows")

    if profile is not None and profile.data.get("quality", {}).get("contract") == "final_size_v2":
        _check_final_size_contract(manifest, profile, repo_root, errors)

    return errors


def _check_final_size_contract(manifest: dict[str, Any], profile: Any, root: Path, errors: list[str]) -> None:
    """Apply v2 display checks equally to review and strict assets.

    Passing these checks never promotes the source numerical evidence.
    """
    policy = profile.data["quality"]
    rendering = manifest.get("rendering", {})
    if not isinstance(rendering, dict):
        errors.append("v2 rendering metadata must be an object")
        return
    quality = rendering.get("quality", {})
    if not isinstance(quality, dict):
        errors.append("v2 rendering.quality must be an object")
        return
    outputs = manifest.get("outputs", [])
    if not isinstance(outputs, list):
        return
    formats = {item.get("format") for item in outputs if isinstance(item, dict)}
    delivery_stage = rendering.get("delivery_stage")
    if delivery_stage not in {None, "png_review", "vector_delivery"}:
        errors.append(f"v2 delivery_stage is invalid: {delivery_stage!r}")
    if delivery_stage == "png_review":
        required_formats = {"png"}
        if formats != {"png"}:
            errors.append("v2 png_review must contain PNG only; vector delivery is a later stage")
        if manifest.get("manuscript_eligible") is not False:
            errors.append("v2 png_review must keep manuscript_eligible=false")
        if manifest.get("current_publication_layer") is not False:
            errors.append("v2 png_review must keep current_publication_layer=false")
        if rendering.get("vector_delivery_pending") is not True:
            errors.append("v2 png_review must record vector_delivery_pending=true")
    else:
        required_formats = set(profile.formats)
    if required_formats - formats:
        errors.append(f"v2 outputs missing required formats: {sorted(required_formats - formats)}")
    width = quality.get("intended_width_inches")
    if not isinstance(width, (int, float)) or not math.isfinite(width) or width <= 0:
        errors.append("v2 quality.intended_width_inches must be positive")
        return
    if width > float(policy["max_width_inches"]):
        errors.append("v2 final width exceeds journal profile maximum")
    declared_size = quality.get("figure_size_inches")
    if (not isinstance(declared_size, list) or len(declared_size) != 2
            or any(not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0 for value in declared_size)):
        errors.append("v2 measured figure_size_inches is required")
        return
    expected_size = profile.data["figure_size_in"].get(rendering.get("column"))
    if expected_size is None:
        errors.append("v2 rendering.column must be single_column or double_column")
    elif not rendering.get("size_override_reason") and declared_size != expected_size:
        errors.append("v2 figure size override requires a documented reason")
    if rendering.get("figure_size_inches") != declared_size:
        errors.append("v2 rendering size must equal measured figure size")
    scale = quality.get("placement_scale")
    if not isinstance(scale, (int, float)) or not math.isfinite(scale) or abs(scale - width / declared_size[0]) > 1e-9:
        errors.append("v2 placement_scale disagrees with intended/exported width")
    for field in ("minimum_capital_numeral_height_mm", "minimum_curve_linewidth_pt"):
        value = quality.get(field)
        if not isinstance(value, (int, float)) or not math.isfinite(value):
            errors.append(f"v2 quality.{field} must be a finite measurement")
            return
    for field in ("clipped_text", "text_overlap_pairs", "smallest_glyphs", "tick_axes"):
        if not isinstance(quality.get(field), list):
            errors.append(f"v2 quality.{field} must contain measurement evidence")
            return
    glyphs = quality["smallest_glyphs"]
    if not glyphs or any(not isinstance(glyph, dict) or not isinstance(glyph.get("final_height_mm"), (int, float)) for glyph in glyphs):
        errors.append("v2 smallest_glyphs must record actual displayed glyphs")
    elif abs(min(glyph["final_height_mm"] for glyph in glyphs) - quality["minimum_capital_numeral_height_mm"]) > 1e-9:
        errors.append("v2 minimum glyph height disagrees with glyph evidence")
    for record in outputs:
        if not isinstance(record, dict) or record.get("format") not in {"pdf", "png"}:
            continue
        output = _resolve_artifact_path(str(record.get("path", "")), root)
        if not output.is_file():
            continue
        try:
            inspection = inspect_export(output)
        except (RuntimeError, ValueError, OSError) as exc:
            errors.append(str(exc))
            continue
        if record.get("inspection") != inspection:
            errors.append(f"v2 exported inspection evidence mismatch: {record.get('path')}")
        if any(abs(float(actual) - float(expected)) > 0.005 for actual, expected in zip(inspection["physical_size_inches"], declared_size)):
            errors.append("v2 exported physical size differs from measured figure size (tight crop or resize)")
        if record["format"] == "pdf":
            if record.get("vector") is not True or inspection["raster_image_count"]:
                errors.append("v2 line chart PDF must contain vector curves, not embedded raster images")
            if not inspection["fonts_embedded"] or inspection["type3_font_count"]:
                errors.append("v2 PDF fonts must be embedded and not Type 3")
            if inspection["page_count"] != 1:
                errors.append("v2 chart PDF must have exactly one page")
        elif inspection["size_pixels"][0] / width + 1 < profile.dpi:
            errors.append("v2 PNG effective dpi at final placement is below profile dpi")
    glyph_height = quality.get("minimum_capital_numeral_height_mm", 0)
    typography_exception = rendering.get("typography_exception")
    if typography_exception in {"dense_composite_review_compact_legend", "dense_composite_review_compact_typography"}:
        if (manifest.get("figure_mode") != "audit" or manifest.get("publication_scope") != "internal_review"
                or delivery_stage != "png_review" or manifest.get("manuscript_eligible") is not False
                or manifest.get("current_publication_layer") is not False):
            errors.append("compact typography exception is restricted to non-eligible audit/internal PNG review")
        if float(glyph_height) < 1.5:
            errors.append("compact typography review glyph height is below 1.5 mm")
    elif float(glyph_height) < float(policy["min_capital_numeral_height_mm"]):
        errors.append("v2 final capital/numeral glyph height is below 2 mm")
    if float(quality.get("minimum_curve_linewidth_pt", 0)) < float(policy["min_curve_linewidth_pt"]):
        errors.append("v2 final curve linewidth is below profile minimum")
    marker = quality.get("minimum_landmark_diameter_mm")
    if marker is not None and (not isinstance(marker, (int, float)) or not math.isfinite(marker) or marker < float(policy["min_landmark_diameter_mm"])):
        errors.append("v2 final landmark diameter is below 1 mm")
    if quality.get("clipped_text") or quality.get("text_overlap_pairs"):
        errors.append("v2 contains clipped or overlapping text")
    if quality.get("legend_axes_overlap_count") != 0:
        allowed_in_axes_policy = rendering.get("legend_policy") in {
            "shared_in_first_panel_reviewed",
            "shared_in_top_row_panel_reviewed",
            "best_in_axes_reviewed",
            "shared_in_panel_reviewed_geometry_checked",
        }
        if not allowed_in_axes_policy:
            errors.append("v2 legend must be outside data axes; in-axes layout requires a separately reviewed contract")
        elif rendering.get("legend_policy") == "shared_in_panel_reviewed_geometry_checked":
            if quality.get("legend_in_axes_overflow_count") != 0:
                errors.append("v2 geometry-checked in-axes legend exceeds its host axes")
            for kind in ("curve", "landmark"):
                overlap_count = quality.get(f"legend_{kind}_overlap_count")
                records = quality.get(f"legend_{kind}_overlaps")
                if not isinstance(overlap_count, int) or not isinstance(records, list):
                    errors.append(f"v2 geometry-checked in-axes legend requires {kind} overlap evidence")
                elif overlap_count != len(records):
                    errors.append(f"v2 legend {kind} overlap count disagrees with evidence")
                elif overlap_count != 0:
                    errors.append(f"v2 geometry-checked in-axes legend intersects a plotted {kind}")
    ticks = quality.get("tick_axes", [])
    if not ticks or any(not isinstance(item, dict) or not item.get("inward") or not item.get("both_sides") or not item.get("minor_count") or not item.get("major_count") for item in ticks):
        errors.append("v2 requires inward major/minor ticks on all four sides")
    axes = manifest.get("axes", [])
    for axis in axes if isinstance(axes, list) else []:
        if not isinstance(axis, dict):
            continue
        label = axis.get("label")
        if not isinstance(label, str) or not label:
            errors.append("v2 axes.label must record the rendered label")
        elif axis.get("display_unit") != "dimensionless" and ("[" in label or "]" in label or "(" not in label or ")" not in label):
            errors.append("v2 dimensional axes must use parentheses for units")
    if rendering.get("color_route") == "color_online_grayscale_print" and not ({"eps", "ps"} & formats):
        errors.append("APS color-online/grayscale-print production route requires PS/EPS")
    if rendering.get("color_route") not in {"color_print_and_online", "color_online_grayscale_print", "undecided_review"}:
        errors.append("v2 rendering.color_route is required")
    if manifest.get("figure_mode") == "strict" and rendering.get("color_route") == "undecided_review":
        errors.append("strict v2 production must select its color/print route")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--repo-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--code-ref", help="full commit SHA for historical code/contract records")
    parser.add_argument("--snapshot", action="store_true", help="verify a retained manifest graph against its committed snapshot")
    args = parser.parse_args()
    if args.snapshot:
        if not args.code_ref:
            parser.error("--snapshot requires --code-ref")
        errors = validate_snapshot(args.manifest, root=args.repo_root.resolve(), code_ref=args.code_ref)
    else:
        errors = validate_manifest(args.manifest, repo_root=args.repo_root.resolve(), code_ref=args.code_ref)
    if errors:
        print(f"[plot-validator] FAILED: {len(errors)} violation(s)")
        for error in errors:
            print(f" - {error}")
        return 1
    scope = "historical snapshot integrity; no current-contract promotion" if args.snapshot else "artifact contract"
    print(f"[plot-validator] OK ({scope}): {args.manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
