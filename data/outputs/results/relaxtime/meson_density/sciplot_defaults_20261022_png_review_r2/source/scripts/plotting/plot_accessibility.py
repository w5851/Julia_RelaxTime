"""Measure color contrast and declare an independent parameter encoding."""

from __future__ import annotations

import math
from typing import Any


def contrast_ratio(foreground: str, background: str) -> float:
    """WCAG 2.2 sRGB relative-luminance contrast for two opaque colors."""
    from matplotlib.colors import to_rgba

    def luminance(color: str) -> float:
        red, green, blue, alpha = to_rgba(color)
        if alpha != 1:
            raise ValueError("contrast requires composited opaque colors")
        linear = [v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4
                  for v in (red, green, blue)]
        return sum(v * weight for v, weight in zip(linear, (0.2126, 0.7152, 0.0722)))

    low, high = sorted((luminance(foreground), luminance(background)))
    return (high + 0.05) / (low + 0.05)


def accessibility_record(parameters: list[dict[str, Any]], *, backgrounds: list[str],
                         background_scope: str) -> dict[str, Any]:
    """Record the actual backgrounds; callers must account for any layered fill."""
    if not parameters or not backgrounds or not background_scope.strip():
        raise ValueError("parameter encodings and explained adjacent backgrounds are required")
    records = []
    for parameter in parameters:
        records.append({**parameter, "contrasts": [contrast_ratio(parameter["color"], bg) for bg in backgrounds]})
    return {"schema": "plot_accessibility_v1", "criterion": "WCAG 2.2 1.4.1 and 1.4.11 (AA)",
            "minimum_contrast": 3.0, "backgrounds": backgrounds,
            "background_scope": background_scope, "parameters": records}


def validate_accessibility(record: Any) -> list[str]:
    """Recompute contrast; a named palette or a supplied pass flag is insufficient."""
    if not isinstance(record, dict) or record.get("schema") != "plot_accessibility_v1":
        return ["accessibility requires a plot_accessibility_v1 record"]
    errors = []
    try:
        expected = accessibility_record(record["parameters"], backgrounds=record["backgrounds"],
                                        background_scope=record["background_scope"])
        if record.get("minimum_contrast") != 3.0:
            errors.append("essential graphical contrast threshold must be 3:1")
        for actual, measured in zip(record["parameters"], expected["parameters"]):
            supplied = actual.get("contrasts")
            if (not isinstance(supplied, list) or len(supplied) != len(measured["contrasts"])
                    or any(not isinstance(a, (float, int)) or not math.isfinite(a)
                           or not math.isclose(a, b, abs_tol=1e-9, rel_tol=0)
                           for a, b in zip(supplied, measured["contrasts"]))):
                errors.append("accessibility contrast measurements disagree with recorded colors")
            if min(measured["contrasts"]) < 3.0:
                errors.append(f"parameter {actual.get('value')!r} has essential contrast below 3:1")
        encodings = [item.get("noncolor_encoding") for item in record["parameters"]]
        if len(encodings) > 1 and (any(not isinstance(value, str) or not value.strip() for value in encodings)
                                   or len(set(encodings)) != len(encodings)):
            errors.append("parameters require distinct noncolor encodings (symbols, labels, or line styles)")
    except (KeyError, TypeError, ValueError, AttributeError) as exc:
        errors.append(f"invalid accessibility evidence: {exc}")
    return errors
