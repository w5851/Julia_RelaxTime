"""Validate the machine-readable contract of a new plot artifact."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

_SCRIPT_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_SCRIPT_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_PROJECT_ROOT))

from scripts.plotting.plot_manifest import MANIFEST_SCHEMA, PROJECT_ROOT, sha256_file
from scripts.plotting.plot_style import ALLOWED_PROFILES, load_profile


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


def _check_hash_record(record: dict[str, Any], *, root: Path, label: str, errors: list[str]) -> Path | None:
    value = record.get("path")
    if not isinstance(value, str) or not value:
        errors.append(f"{label}.path must be a non-empty string")
        return None
    path = _resolve_artifact_path(value, root)
    if not path.is_file():
        errors.append(f"{label} missing file: {value}")
        return None
    expected_bytes = record.get("bytes")
    if expected_bytes is not None and int(expected_bytes) != path.stat().st_size:
        errors.append(f"{label}.bytes mismatch for {value}")
    expected_hash = record.get("sha256")
    if not isinstance(expected_hash, str) or len(expected_hash) != 64:
        errors.append(f"{label}.sha256 must be a 64-character hash")
    elif sha256_file(path) != expected_hash:
        errors.append(f"{label}.sha256 mismatch for {value}")
    return path


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


def validate_manifest(manifest_path: str | Path, *, repo_root: Path = PROJECT_ROOT) -> list[str]:
    """Return all contract violations; an empty list means the artifact passes."""

    path = Path(manifest_path).resolve()
    errors: list[str] = []
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
    if isinstance(generator, dict):
        _check_hash_record(generator, root=repo_root, label="generator", errors=errors)
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
            _check_hash_record(record, root=repo_root, label=label, errors=errors)

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
        missing_formats = {"png", "svg"} - formats
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
                if expected_size is not None and any(abs(float(a) - float(b)) > 1.0e-9 for a, b in zip(declared_size, expected_size)):
                    errors.append(f"strict rendering size {declared_size} does not match profile {expected_size}")

    if mode == "strict" and any(token in str(manifest.get("semantic_status", "")).lower() for token in ("unresolved", "estimated", "bracket")):
        errors.append("strict semantic_status cannot claim unresolved, estimated, or bracket content")

    if profile is not None and mode == "strict" and profile.data["semantics"].get("allow_connector") is not False:
        errors.append("strict profile must disallow connector rows")

    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--repo-root", type=Path, default=PROJECT_ROOT)
    args = parser.parse_args()
    errors = validate_manifest(args.manifest, repo_root=args.repo_root.resolve())
    if errors:
        print(f"[plot-validator] FAILED: {len(errors)} violation(s)")
        for error in errors:
            print(f" - {error}")
        return 1
    print(f"[plot-validator] OK: {args.manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
