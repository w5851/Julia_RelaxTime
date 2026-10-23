"""Restyle the frozen paper Fig.4 as a PNG review with legends inside panels.

The hash-checked original renderer selects every curve and branch endpoint.
Its file writes are suppressed; only layout and typography are changed here.
No solver, new interpolation, or manuscript/pointer update is performed.
"""

from __future__ import annotations

import argparse
from contextlib import redirect_stdout
import hashlib
import io
import json
from pathlib import Path
import runpy
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd

from scripts.plotting.plot_manifest import (
    build_manifest, generator_record, input_record, runtime_record,
    sha256_file, write_manifest,
)
from scripts.plotting.plot_quality import _line_intersects_bbox, export_figure, measure_figure
from scripts.plotting.plot_style import configure_axis_ticks, configure_matplotlib, load_profile
from scripts.plotting.validate_plot_artifact import validate_manifest


POINTS = ROOT / "docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v5/tables/publication_clean_points.csv"
SCAN = ROOT / "data/outputs/results/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2/phase_guided_transport_scan.csv"
PAPER_FIGURE = Path("My_Paper/pnjl_aniso_transport/writing_package/outputs/journal_en/review_revision_2026-10-12/figure_fixed_T")
FROZEN_HASHES = {
    "publication_clean_points.csv": "3783e96be9a14ac7bd2cd44f3d29526fe50de5402d49624574e655cd216b3186",
    "phase_guided_transport_scan.csv": "3b73b937f1f09cf7112b3a42219088f98b852e5ef0354174b66f71823f7bb29c",
    "plot_fixed_T_transport.py": "0e4987b45677eff3140691aa8fb05da56dee2e7ed604e337a06df40c205a9e38",
    "figure4_transport_fixed_T.png": "4e6df8201724b227872608ca5a7fdd3961ac13ad64bb01a2900b35155e455b81",
    "figure4_transport_fixed_T.pdf": "7e807fd80f003034e2afb6ddc7272244b552fa40dbc035e7aec8f29a3b510e9d",
    "first_order_endpoints.csv": "b3dabf76ef6a1ad8a5c721bd892148d6f5d7a3781dcf1484f923221f5a20d6dc",
}


def frozen_inputs(paper_root: Path) -> list[Path]:
    files = [POINTS, SCAN, *(paper_root / PAPER_FIGURE / name for name in list(FROZEN_HASHES)[2:])]
    for path in files:
        if sha256_file(path) != FROZEN_HASHES[path.name]:
            raise ValueError(f"frozen source changed: {path}")
    return files


def source_figure(script: Path) -> dict:
    # This particular script is inspected and hash-checked before runpy.
    with TemporaryDirectory(prefix="transport_fig4_source_") as scratch:
        argv = [str(script), str(POINTS), str(SCAN), scratch]
        with patch.object(sys, "argv", argv), patch.object(Figure, "savefig"), \
                patch.object(pd.DataFrame, "to_csv"), redirect_stdout(io.StringIO()):
            return runpy.run_path(str(script), run_name="__main__")


def plotted_state(fig: Figure) -> list[dict]:
    return [{
        "panel": chr(97 + i), "xlim": list(ax.get_xlim()), "ylim": list(ax.get_ylim()),
        "xscale": ax.get_xscale(), "yscale": ax.get_yscale(),
        "lines": [{"xy": line.get_xydata().tolist(), "marker": line.get_marker(),
                   "linestyle": line.get_linestyle(), "color": line.get_color()}
                  for line in ax.lines],
    } for i, ax in enumerate(fig.axes)]


def state_digest(state: list[dict]) -> str:
    return hashlib.sha256(json.dumps(state, sort_keys=True, allow_nan=False).encode()).hexdigest()


def box_is_clear(ax, box: Bbox, renderer) -> bool:
    padding = renderer.points_to_pixels(2.0)
    padded = box.padded(padding)
    if not ax.bbox.contains(padded.x0, padded.y0) or not ax.bbox.contains(padded.x1, padded.y1):
        return False
    for line in ax.lines:
        if line.get_linestyle() not in {"None", "", " "} and _line_intersects_bbox(line, padded):
            return False
        if line.get_marker() not in {None, "None", "", " "}:
            radius = renderer.points_to_pixels((line.get_markersize() + line.get_markeredgewidth()) / 2)
            expanded = padded.padded(radius)
            if any(expanded.contains(*xy) for xy in line.get_transform().transform(line.get_xydata())):
                return False
    return not any(box.overlaps(text.get_window_extent(renderer)) for text in ax.texts if text.get_visible())


def place_legend(fig, axes, handles, title=None):
    """Search empty rectangles at the unchanged ranges; never shrink text."""
    renderer = fig.canvas.get_renderer()
    for ax in axes:
        legend = ax.legend(handles=handles, title=title, loc="lower left", borderaxespad=0,
                           borderpad=0.25, labelspacing=0.20, handlelength=1.6,
                           handletextpad=0.4, frameon=False, fontsize=13, title_fontsize=13)
        box = legend.get_window_extent(renderer)
        if box.width + 6 >= ax.bbox.width or box.height + 6 >= ax.bbox.height:
            legend.remove()
            continue
        xs = np.linspace(ax.bbox.x0 + 4, ax.bbox.x1 - box.width - 4, 25)
        ys = np.linspace(ax.bbox.y0 + 4, ax.bbox.y1 - box.height - 4, 61)
        # Prefer the middle of the available whitespace over hugging a border.
        candidates = sorted(((x, y) for x in xs for y in ys), key=lambda xy:
                            ((xy[0] + box.width / 2 - ax.bbox.x0) / ax.bbox.width - 0.5) ** 2
                            + ((xy[1] + box.height / 2 - ax.bbox.y0) / ax.bbox.height - 0.5) ** 2)
        for x, y in candidates:
            candidate = Bbox.from_bounds(x, y, box.width, box.height)
            if box_is_clear(ax, candidate, renderer):
                anchor = ax.transAxes.inverted().transform((x, y))
                legend.set_bbox_to_anchor(anchor, transform=ax.transAxes)
                return {"panel": chr(97 + fig.axes.index(ax)), "anchor_axes": anchor.tolist(),
                        "labels": [item.get_label() for item in handles], "title": title}
        legend.remove()
    return None


def endpoint_quality(fig) -> dict:
    """Also check Line2D endpoint markers, which the shared checker does not cover."""
    renderer = fig.canvas.get_renderer()
    boxes = [ax.get_legend().get_window_extent(renderer) for ax in fig.axes if ax.get_legend() is not None]
    overlaps, diameters = [], []
    count = 0
    for index, ax in enumerate(fig.axes):
        for line in ax.lines:
            if line.get_marker() in {None, "None", "", " "}:
                continue
            diameter = line.get_markersize()
            diameters.append(diameter * 25.4 / 72)
            radius = renderer.points_to_pixels((diameter + line.get_markeredgewidth()) / 2)
            for xy in line.get_transform().transform(line.get_xydata()):
                count += 1
                if any(box.padded(radius).contains(*xy) for box in boxes):
                    overlaps.append({"panel": chr(97 + index), "xy_pixels": xy.tolist()})
    return {"marker_count": count, "minimum_diameter_mm": min(diameters), "legend_overlaps": overlaps}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paper-root", type=Path, default=ROOT.parent / "paper")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    if out.exists():
        raise FileExistsError(f"refusing existing case: {out}")
    sources = frozen_inputs(args.paper_root)
    original = source_figure(sources[2])
    fig, axes = original["fig"], original["axes"]
    before = plotted_state(fig)
    pts = original["pts"]
    assert np.isfinite(pts[["xi", "clean_value"]].to_numpy()).all()
    assert not pts.duplicated(["plot_panel", "plot_series", "observable", "xi"]).any()
    expected_endpoints = pd.read_csv(sources[-1])
    actual_endpoints = pd.DataFrame(original["records"], columns=expected_endpoints.columns)
    pd.testing.assert_frame_equal(actual_endpoints, expected_endpoints, check_dtype=False)
    assert len(actual_endpoints) == 3
    assert (actual_endpoints.xi_restored_end == -0.13).all()
    assert (actual_endpoints.xi_broken_end == -0.12).all()

    for legend in list(fig.legends):
        legend.remove()
    profile = load_profile("candidate_aps_v2")
    font = configure_matplotlib(profile)
    fig.set_size_inches(6.75, 7.4)
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.09, top=0.95, hspace=0.26, wspace=0.38)
    for ax in axes.flat:
        configure_axis_ticks(ax, profile)
    for text in fig.findobj(Text):
        text.set_fontsize(13)
    fig.canvas.draw()

    # First test all five original entries in a single panel, as requested.
    panel_order = [axes[0, 0], axes[1, 0], axes[2, 1], axes[2, 0], axes[2, 2],
                   axes[0, 1], axes[1, 1], axes[0, 2], axes[1, 2]]
    single = place_legend(fig, panel_order, original["handles"])
    layout = [single] if single else []
    if single is None:
        # Compact parameter notation keeps all three values together at 13 pt.
        chemical = [Line2D([], [], color=original["COLORS"][i], ls=original["STYLES"][i],
                           label=str(int(mu))) for i, mu in enumerate(original["MUB"])]
        first = place_legend(fig, panel_order, chemical, r"$\mu_B\;(\mathrm{MeV})$")
        if first is None:
            raise RuntimeError("no clear panel for the three chemical potentials")
        available = [ax for ax in panel_order if ax.get_legend() is None]
        second = place_legend(fig, available, original["handles"][3:], "First-order")
        if second is None:
            raise RuntimeError("no clear panel for the two first-order endpoints")
        layout = [first, second]
    after = plotted_state(fig)
    assert before == after, "curve samples, branch markers, or axis ranges changed"
    quality = measure_figure(fig, intended_width_inches=6.75)
    endpoints = endpoint_quality(fig)
    assert not endpoints["legend_overlaps"]
    for field in ("clipped_text", "text_overlap_pairs", "legend_curve_overlaps", "legend_in_axes_overflows"):
        if quality[field]:
            raise RuntimeError(f"layout check failed: {field}={quality[field]}")
    if quality["minimum_capital_numeral_height_mm"] < 2:
        raise RuntimeError(f"glyphs too small: {quality['smallest_glyphs']}")

    out.mkdir(parents=True, exist_ok=False)
    archive = out / "source_snapshot"
    archive.mkdir()
    generator_copy = archive / Path(__file__).name
    shutil.copy2(__file__, generator_copy)
    source_copy = archive / sources[2].name
    shutil.copy2(sources[2], source_copy)
    input_records = [input_record(path, role=role) for path, role in zip(sources, [
        "frozen_display_points", "raw_branch_support", "original_paper_renderer",
        "original_paper_png", "original_paper_pdf", "original_branch_endpoints"])]
    input_records.append(input_record(source_copy, role="original_renderer_snapshot"))
    for relative in ["scripts/plotting/plot_style.py", "scripts/plotting/plot_manifest.py",
                     "scripts/plotting/plot_quality.py", "scripts/plotting/validate_plot_artifact.py",
                     "config/plotting/candidate_aps_v2.toml"]:
        path = ROOT / relative
        snapshot = archive / path.name
        shutil.copy2(path, snapshot)
        record = input_record(path, role="plotting_contract_source")
        record["source_snapshot"] = input_record(snapshot, role="plotting_contract_snapshot")
        input_records.append(record)
    outputs, quality = export_figure(fig, out / "figure4_transport_fixed_T", profile,
                                    intended_width_inches=6.75, formats=("png",))
    frozen_inputs(args.paper_root)
    command = subprocess.list2cmdline([sys.executable, str(Path(__file__).resolve()),
                                      "--paper-root", str(args.paper_root.resolve()), "--output-dir", str(out)])
    manifest = build_manifest(
        asset_id=out.name, figure_family="transport_fixed_T_comparison", case_slug=out.name,
        figure_mode="audit", semantic_status="inherited_v5_display_values_layout_review",
        style_profile=profile.profile_id, publication_scope="internal_review",
        generator=generator_record(__file__, command=command, source_snapshot=generator_copy,
                                   runtime=runtime_record({"matplotlib": matplotlib.__version__,
                                                           "numpy": np.__version__, "pandas": pd.__version__, "font": font})),
        inputs=input_records,
        axes=[{"field": field, "source_unit": "dimensionless", "display_unit": "dimensionless",
               "label": label, "transform": "identity", "scale": scale}
              for field, label, scale in [("xi", r"$\xi$", "linear"), *original["ROWS"]]],
        series=[{"series_id": f"{obs}_T{temp:g}_muB{mu:g}", "state": "inherited_v5_display",
                 "support_rule": "all frozen mode_b clean_value rows ordered by xi",
                 "mask_rule": "split at the original mass-gap branch switch; retain both endpoints"}
                for obs, _, _ in original["ROWS"] for temp in original["TEMPS"] for mu in original["MUB"]],
        outputs=outputs,
        selection_rule="exact hash-checked paper renderer; mode_b; T=120,160,200 MeV; muB=0,450,900 MeV",
        interpolation_policy="inherited v5 display adjustments; no new interpolation or replacement",
        connector_policy="forbidden across the first-order gap",
        missing_value_policy="reject nonfinite selected display values; no gap filling",
        validation={"finite": True, "duplicate_keys": True, "support": True, "strict_gate": False,
                    "original_and_restyled_plotted_state_equal": True,
                    "plotted_state_sha256": state_digest(after),
                    "curve_and_marker_artist_count": sum(len(ax.lines) for ax in fig.axes),
                    "frozen_inputs_unchanged_after_render": True,
                    "first_order_endpoints": actual_endpoints.to_dict(orient="records"),
                    "line2d_endpoint_quality": endpoints},
        rendering={"column": "double_column", "figure_size_inches": [6.75, 7.4],
                   "size_override_reason": "nine panels with 13 pt in-panel legends at fixed width and unchanged axis ranges",
                   "delivery_stage": "png_review", "vector_delivery_pending": True,
                   "color_route": "undecided_review", "quality": quality,
                   "legend_policy": "shared_in_panel_reviewed_geometry_checked",
                   "legend_layout": layout, "all_five_entries_fit_one_panel": single is not None,
                   "in_panel_exception_basis": "author group-meeting request on 2026-10-21; geometric checks plus agent visual review; author acceptance pending"},
        calculation_sha="3c5f6b3c9bd535cff7657364dadb2efc31f2ea48",
        source_run_id="first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2",
    )
    manifest.update({"manuscript_eligible": False, "current_publication_layer": False,
                     "solver_called": False, "canonical_data_modified": False, "new_display_values": False})
    manifest_path = out / "plot_manifest.json"
    write_manifest(manifest_path, manifest)
    errors = validate_manifest(manifest_path, repo_root=ROOT)
    if errors:
        raise RuntimeError("; ".join(errors))
    print(json.dumps({"manifest": str(manifest_path), "legend_layout": layout,
                      "minimum_glyph_height_mm": quality["minimum_capital_numeral_height_mm"],
                      "endpoints": endpoints, "artifact_validation": "passed"}, indent=2))
    plt.close(fig)


if __name__ == "__main__":
    main()
