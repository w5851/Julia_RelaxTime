"""Shared review styling, raster companions and SOP snapshot provenance."""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path

import matplotlib
from PIL import Image

from scripts.plotting.plot_manifest import output_record, sha256_file
from scripts.plotting.plot_quality import inspect_export
from scripts.plotting.plot_series import SeriesStyleRegistry
from scripts.plotting.plot_style import style_context


REVIEW_COLORS = ("#0072b2", "#d55e00", "#009e73")
WIDTH_IN, LABEL_PT, TICK_PT, LEGEND_PT = 6.75, 13.0, 11.0, 11.0
REVIEW_DPI = 600
SOP_SOURCE = Path(__file__).resolve().parents[2] / "docs/guides/sop/workflows/figure_production.md"


def write_sop_snapshot(results_dir: Path, *, source: Path = SOP_SOURCE) -> Path:
    """Retain the SOP text under the case results, using canonical LF bytes."""
    result = results_dir / "figure_production_snapshot.md"
    payload = source.read_bytes().replace(b"\r\n", b"\n")
    result.parent.mkdir(parents=True, exist_ok=True)
    try:
        with result.open("xb") as stream:
            stream.write(payload)
    except FileExistsError:
        if result.read_bytes() != payload:
            raise ValueError(f"existing figure SOP snapshot differs: {result}") from None
    return result


def grayscale_review_record(color: Path) -> dict:
    """Write a same-size, single-channel 600 dpi PNG and bind it to its color source."""
    gray = color.with_name(color.stem + "_gray.png")
    with Image.open(color) as image:
        image.convert("L").save(gray, dpi=(REVIEW_DPI, REVIEW_DPI))
    record = output_record(gray, fmt="png", dpi=REVIEW_DPI, vector=False)
    record.update(inspection=inspect_export(gray), role="grayscale_review",
                  source_color_sha256=sha256_file(color))
    return record


def register_shared_series_styles(series: list[dict], registry: SeriesStyleRegistry) -> dict:
    """Record one style identity for curves that share a visual encoding across panels."""
    for style in registry.styles:
        members = [row["series_id"] for row in series if style.series_id in row.get("style_ids", ())]
        if not members:
            raise ValueError(f"style has no source curves: {style.series_id}")
        series.append(dict(series_id=style.series_id, state="shared_style_identity",
                           support_rule="same encoding across the named frozen source series",
                           mask_rule="inherits the source series masks", member_series_ids=members))
    return registry.record(style.series_id for style in registry.styles)


@contextmanager
def figure_context(profile, *, compact=False):
    with style_context(profile) as context:
        with matplotlib.rc_context({"axes.labelsize": LABEL_PT, "axes.titlesize": LABEL_PT,
                                    "xtick.labelsize": TICK_PT if compact else LABEL_PT,
                                    "ytick.labelsize": TICK_PT if compact else LABEL_PT,
                                    "path.simplify": False}):
            yield context
