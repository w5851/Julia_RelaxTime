"""Measure final-size typography, layout, and exported chart contents.

The checks supplement, rather than replace, author review of scientific
meaning, grayscale legibility, and curve visibility.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
import math
import re
import subprocess
from typing import Any, Iterable

from scripts.plotting.plot_manifest import output_record
from scripts.plotting.plot_style import PlotProfile


MM_PER_POINT = 25.4 / 72


def _poppler(command: list[str]) -> str:
    try:
        return subprocess.check_output(command, text=True, encoding="utf-8", errors="replace")
    except (OSError, subprocess.CalledProcessError) as exc:
        raise RuntimeError(f"artifact inspection requires working {command[0]}") from exc


def inspect_export(path: Path) -> dict[str, Any]:
    """Inspect actual dimensions, resolution, PDF images, and embedded fonts."""
    if path.suffix.lower() == ".png":
        from PIL import Image

        with Image.open(path) as image:
            dpi = image.info.get("dpi")
            if not dpi or min(dpi) <= 0:
                raise ValueError(f"PNG has no physical resolution: {path}")
            return {
                "size_pixels": list(image.size),
                "actual_dpi": list(dpi),
                "physical_size_inches": [image.width / dpi[0], image.height / dpi[1]],
            }
    if path.suffix.lower() != ".pdf":
        raise ValueError(f"unsupported inspected format: {path.suffix}")
    info = _poppler(["pdfinfo", str(path)])
    size = re.search(r"Page size:\s*([\d.]+) x ([\d.]+) pts", info)
    pages = re.search(r"Pages:\s*(\d+)", info)
    if size is None or pages is None:
        raise ValueError(f"cannot parse PDF physical dimensions: {path}")
    images = _poppler(["pdfimages", "-list", str(path)])
    image_rows = [line for line in images.splitlines() if re.match(r"\s*\d+\s+\d+\s+", line)]
    fonts = _poppler(["pdffonts", str(path)])
    font_rows = [
        line for line in fonts.splitlines()[2:]
        if line.strip() and not set(line.strip()) <= {"-"}
    ]
    return {
        "physical_size_inches": [float(size[1]) / 72, float(size[2]) / 72],
        "page_count": int(pages[1]),
        "raster_image_count": len(image_rows),
        "font_count": len(font_rows),
        "fonts_embedded": bool(font_rows) and all(line.split()[-5] == "yes" for line in font_rows),
        "type3_font_count": sum("Type 3" in line for line in font_rows),
        "inspection_tools": ["pdfinfo", "pdfimages -list", "pdffonts"],
    }


@lru_cache(maxsize=4096)
def _glyph_height_mm(font_path: str, size: float, char: str) -> float:
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextPath

    return float(TextPath((0, 0), char, prop=FontProperties(fname=font_path, size=size)).get_extents().height) * MM_PER_POINT


def _text_glyphs(artist: Any) -> list[tuple[str, float, str]]:
    from matplotlib import font_manager
    from matplotlib.mathtext import MathTextParser

    text = artist.get_text()
    properties = artist.get_fontproperties()
    if "$" in text:
        parsed = MathTextParser("path").parse(text, dpi=72, prop=properties)
        return [(font.fname, float(size), chr(code)) for font, size, code, _, _ in parsed.glyphs]
    font_path = font_manager.findfont(properties)
    return [(font_path, float(properties.get_size_in_points()), char) for char in text]


def visible_legends(figure: Any) -> list[Any]:
    """Include every public Legend artist, also those retained with add_artist."""
    from matplotlib.legend import Legend

    legends = [*figure.legends, *(child for ax in figure.axes for child in ax.get_children()
                                 if isinstance(child, Legend))]
    return list({id(item): item for item in legends if item.get_visible()}.values())


def visible_texts(figure: Any) -> list[Any]:
    """Skip unused/out-of-range tick objects retained by Matplotlib locators."""
    items = list(figure.texts)
    for legend in visible_legends(figure):
        items.extend([legend.get_title(), *legend.get_texts()])
    for ax in figure.axes:
        items.extend([ax.xaxis.label, ax.yaxis.label, ax.title, *ax.texts])
        for axis in (ax.xaxis, ax.yaxis):
            low, high = sorted(axis.get_view_interval())
            items.append(axis.get_offset_text())
            for tick in [*axis.get_major_ticks(), *axis.get_minor_ticks()]:
                if low - 1e-10 <= tick.get_loc() <= high + 1e-10:
                    items.extend([tick.label1, tick.label2])
    return list({id(item): item for item in items if item.get_visible() and item.get_text().strip()}.values())


def _point_in_bbox(x: float, y: float, bbox: Any) -> bool:
    return bbox.x0 <= x <= bbox.x1 and bbox.y0 <= y <= bbox.y1


def _segment_intersects_bbox(
    start: tuple[float, float],
    end: tuple[float, float],
    bbox: Any,
) -> bool:
    """Return whether a display-space line segment meets a rectangle."""

    if _point_in_bbox(*start, bbox) or _point_in_bbox(*end, bbox):
        return True
    dx = end[0] - start[0]
    dy = end[1] - start[1]
    lower, upper = 0.0, 1.0
    for coefficient, constant in (
        (-dx, start[0] - bbox.x0),
        (dx, bbox.x1 - start[0]),
        (-dy, start[1] - bbox.y0),
        (dy, bbox.y1 - start[1]),
    ):
        if coefficient == 0.0:
            if constant < 0.0:
                return False
            continue
        ratio = constant / coefficient
        if coefficient < 0.0:
            if ratio > upper:
                return False
            lower = max(lower, ratio)
        else:
            if ratio < lower:
                return False
            upper = min(upper, ratio)
    return lower <= upper


def _line_intersects_bbox(line: Any, bbox: Any) -> bool:
    """Check the rendered polyline, including log/data transforms."""

    import math

    xy = line.get_xydata()
    if len(xy) < 2:
        return False
    transformed = line.get_transform().transform(xy)
    previous: tuple[float, float] | None = None
    for point in transformed:
        current = (float(point[0]), float(point[1]))
        if not all(math.isfinite(value) for value in current):
            previous = None
            continue
        if previous is not None and _segment_intersects_bbox(previous, current, bbox):
            return True
        previous = current
    return False


def measure_figure(figure: Any, *, intended_width_inches: float) -> dict[str, Any]:
    """Measure all displayed capitals/numerals, including math subscripts."""
    from matplotlib.collections import PathCollection
    from matplotlib.transforms import Bbox

    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    size = list(map(float, figure.get_size_inches()))
    scale = float(intended_width_inches) / size[0]
    artists = visible_texts(figure)
    glyphs = []
    clipped = []
    text_boxes = []
    for artist in artists:
        box = artist.get_window_extent(renderer)
        if not figure.bbox.contains(box.x0 + 0.5, box.y0 + 0.5) or not figure.bbox.contains(box.x1 - 0.5, box.y1 - 0.5):
            clipped.append(artist.get_text())
        text_boxes.append((artist.get_text(), box))
        for font_path, font_size, char in _text_glyphs(artist):
            if char.isascii() and (char.isupper() or char.isdigit()):
                glyphs.append({
                    "text": artist.get_text(), "glyph": char,
                    "font": Path(font_path).name, "source_size_pt": font_size,
                    "final_height_mm": _glyph_height_mm(font_path, font_size, char) * scale,
                    "role": "script" if font_size < 0.9 * artist.get_fontsize() else "primary",
                })
    overlaps = []
    for index, (text, box) in enumerate(text_boxes):
        for other_text, other in text_boxes[index + 1:]:
            dx = min(box.x1, other.x1) - max(box.x0, other.x0)
            dy = min(box.y1, other.y1) - max(box.y0, other.y0)
            if dx > 1 and dy > 1:
                overlaps.append([text, other_text])
    legends = visible_legends(figure)
    legend_axes_overlaps = int(sum(
        legend.get_window_extent(renderer).overlaps(ax.get_window_extent(renderer))
        for legend in legends for ax in figure.axes
    ))
    legend_curve_overlaps = []
    legend_landmark_overlaps = []
    legend_in_axes_overflows = []
    legend_layout = []
    legend_pair_overlaps = []
    for legend_index, legend in enumerate(legends):
        box = legend.get_window_extent(renderer)
        host = legend.axes
        host_index = figure.axes.index(host) if host in figure.axes else None
        overlapped_axes = [index for index, ax in enumerate(figure.axes) if box.overlaps(ax.bbox)]
        contained = bool(box.x0 >= host.bbox.x0 and box.y0 >= host.bbox.y0
                         and box.x1 <= host.bbox.x1 and box.y1 <= host.bbox.y1) if host is not None else None
        if host_index in overlapped_axes and not contained:
            legend_in_axes_overflows.append(host_index)
        legend_layout.append({
            "legend_index": legend_index, "host_axes_index": host_index,
            "overlapped_axes": overlapped_axes, "contained_in_host": contained,
            "title": legend.get_title().get_text(), "labels": [text.get_text() for text in legend.get_texts()],
            "bbox_inches": [float(value) / figure.dpi for value in box.bounds],
        })
        for other_index, other in enumerate(legends[:legend_index]):
            other_box = other.get_window_extent(renderer)
            if box.overlaps(other_box):
                legend_pair_overlaps.append([other_index, legend_index])
    for legend_index, legend in enumerate(legends):
        legend_box = legend.get_window_extent(renderer)
        for axis_index, ax in enumerate(figure.axes):
            if not legend_box.overlaps(ax.get_window_extent(renderer)):
                continue
            for line_index, line in enumerate(ax.lines):
                if not line.get_visible() or line.get_linestyle() in {"None", "", " "}:
                    continue
                stroke = renderer.points_to_pixels(line.get_linewidth() / 2)
                padded = Bbox.from_extents(legend_box.x0 - stroke, legend_box.y0 - stroke,
                                          legend_box.x1 + stroke, legend_box.y1 + stroke)
                if _line_intersects_bbox(line, padded):
                    legend_curve_overlaps.append(
                        {
                            "legend_index": legend_index,
                            "axis_index": axis_index,
                            "line_index": line_index,
                            "label": str(line.get_label()),
                        }
                    )
            for collection_index, collection in enumerate(ax.collections):
                if not isinstance(collection, PathCollection) or not collection.get_visible():
                    continue
                offsets = collection.get_offset_transform().transform(collection.get_offsets())
                sizes = collection.get_sizes()
                for point_index, (x, y) in enumerate(offsets):
                    diameter = float(sizes[min(point_index, len(sizes) - 1)]) ** 0.5 if len(sizes) else 0.0
                    radius = renderer.points_to_pixels(diameter / 2)
                    expanded = Bbox.from_extents(legend_box.x0 - radius, legend_box.y0 - radius,
                                                legend_box.x1 + radius, legend_box.y1 + radius)
                    if _point_in_bbox(float(x), float(y), expanded):
                        legend_landmark_overlaps.append({"legend_index": legend_index,
                            "axis_index": axis_index, "collection_index": collection_index,
                            "point_index": point_index})
    linewidths = [
        line.get_linewidth() * scale
        for ax in figure.axes for line in ax.lines
        if line.get_visible() and line.get_linestyle() not in {"None", "", " "}
    ]
    tick_axes = []
    marker_diameters = []
    for ax in figure.axes:
        for name, axis in (("x", ax.xaxis), ("y", ax.yaxis)):
            low, high = sorted(axis.get_view_interval())
            major = [tick for tick in axis.get_major_ticks() if low <= tick.get_loc() <= high]
            minor = [tick for tick in axis.get_minor_ticks() if low <= tick.get_loc() <= high]
            ticks = [*major, *minor]
            tick_axes.append({
                "axis": name, "scale": ax.get_xscale() if name == "x" else ax.get_yscale(),
                "minor_locator": type(axis.get_minor_locator()).__name__,
                "major_count": len(major), "minor_count": len(minor),
                "inward": bool(ticks) and all(tick._tickdir == "in" for tick in ticks),
                "both_sides": bool(ticks) and all(tick.tick1line.get_visible() and tick.tick2line.get_visible() for tick in ticks),
            })
        for collection in ax.collections:
            if hasattr(collection, "get_sizes"):
                paths = collection.get_paths()
                for index, value in enumerate(collection.get_sizes()):
                    if not paths:
                        continue
                    bounds = paths[min(index, len(paths) - 1)].get_extents()
                    marker_diameters.append(min(bounds.width, bounds.height) * float(value) ** 0.5 * MM_PER_POINT * scale)
    return {
        "figure_size_inches": size,
        "intended_width_inches": float(intended_width_inches),
        "placement_scale": scale,
        "glyph_measurement": "font outline heights of displayed ASCII capitals/numerals, including math scripts; em size is not glyph height",
        "minimum_capital_numeral_height_mm": min((item["final_height_mm"] for item in glyphs), default=0),
        "minimum_primary_capital_numeral_height_mm": min((item["final_height_mm"] for item in glyphs if item["role"] == "primary"), default=0),
        "minimum_script_capital_numeral_height_mm": min((item["final_height_mm"] for item in glyphs if item["role"] == "script"), default=None),
        "smallest_glyphs": sorted(glyphs, key=lambda item: item["final_height_mm"])[:8],
        "minimum_curve_linewidth_pt": min(linewidths, default=0),
        "minimum_landmark_diameter_mm": min(marker_diameters, default=None),
        "tick_axes": tick_axes,
        "clipped_text": clipped,
        "text_overlap_pairs": overlaps,
        "legend_axes_overlap_count": legend_axes_overlaps,
        "legend_count": len(legends),
        "legend_layout": legend_layout,
        "legend_pair_overlap_count": len(legend_pair_overlaps),
        "legend_pair_overlaps": legend_pair_overlaps,
        "legend_curve_overlap_count": len(legend_curve_overlaps),
        "legend_curve_overlaps": legend_curve_overlaps,
        "legend_landmark_overlap_count": len(legend_landmark_overlaps),
        "legend_landmark_overlaps": legend_landmark_overlaps,
        "legend_in_axes_overflow_count": len(legend_in_axes_overflows),
        "legend_in_axes_overflows": legend_in_axes_overflows,
        "human_visual_review": "required",
    }


def placement_limits(quality: dict[str, Any], profile: PlotProfile,
                     outputs: Iterable[dict[str, Any]]) -> dict[str, Any]:
    """Compute the usable insertion interval, without granting publication status.

    Glyph/line/marker limits bound reduction. Actual PNG pixels and the project
    width cap bound enlargement. The measurements must refer to one known width.
    """
    width = quality.get("intended_width_inches")
    if not isinstance(width, (int, float)) or not math.isfinite(width) or width <= 0:
        raise ValueError("intended_width_inches must be finite and positive")
    policy = profile.data["quality"]
    fields = {
        "glyph": ("minimum_capital_numeral_height_mm", "min_capital_numeral_height_mm"),
        "curve": ("minimum_curve_linewidth_pt", "min_curve_linewidth_pt"),
        "landmark": ("minimum_landmark_diameter_mm", "min_landmark_diameter_mm"),
    }
    minima = {}
    for name, (field, threshold) in fields.items():
        value = quality.get(field)
        if value is None and name == "landmark":
            continue
        if not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise ValueError(f"{field} must be finite and positive")
        minima[name] = width * float(policy[threshold]) / value
    maxima = {"project_profile": float(policy["max_width_inches"])}
    for index, output in enumerate(outputs):
        if output.get("format") == "png":
            pixels = output.get("inspection", {}).get("size_pixels", [])
            if len(pixels) != 2 or not all(isinstance(x, (int, float)) and math.isfinite(x) and x > 0 for x in pixels):
                raise ValueError("PNG placement requires inspected positive pixel dimensions")
            maxima[f"png_{index}_effective_dpi"] = pixels[0] / profile.dpi
    low, high = max(minima.values()), min(maxima.values())
    single_width = float(profile.data["figure_size_in"]["single_column"][0])
    return {
        "schema": "plot_placement_limits_v1", "measured_width_inches": width,
        "minimum_width_inches": low, "maximum_width_inches": high,
        "minimum_width_constraints": minima, "maximum_width_constraints": maxima,
        "has_usable_interval": low <= high,
        "measured_width_qualified": low <= width <= high,
        "single_column_width_inches": single_width,
        "single_column_reuse_qualified": low <= single_width <= high,
        "single_column_minimum_glyph_mm": quality["minimum_capital_numeral_height_mm"] * single_width / width,
        "scope": "size gates only; layout, grayscale, author acceptance and source qualification remain separate",
    }


def export_figure(
    figure: Any,
    stem: Path,
    profile: PlotProfile,
    *,
    intended_width_inches: float | None = None,
    formats: Iterable[str] | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Export at fixed physical dimensions and return measured quality evidence."""
    stem.parent.mkdir(parents=True, exist_ok=True)
    width = intended_width_inches or float(figure.get_size_inches()[0])
    quality = measure_figure(figure, intended_width_inches=width)
    records = []
    selected_formats = tuple(str(fmt).lower() for fmt in (formats if formats is not None else profile.formats))
    if not selected_formats or len(selected_formats) != len(set(selected_formats)):
        raise ValueError("export formats must be a non-empty unique sequence")
    for fmt in selected_formats:
        output = stem.with_suffix(f".{fmt}")
        if output.exists():
            raise FileExistsError(f"refusing to overwrite chart: {output}")
        figure.savefig(output, format=fmt, dpi=profile.dpi, bbox_inches=None)
        record = output_record(
            output,
            fmt=fmt,
            dpi=profile.dpi if fmt == "png" else None,
            vector=fmt in {"pdf", "eps", "ps", "svg"},
        )
        record["inspection"] = inspect_export(output)
        records.append(record)
    return records, quality
