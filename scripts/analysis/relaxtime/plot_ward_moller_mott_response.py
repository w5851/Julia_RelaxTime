"""Build the Appendix A figure from the verified Ward/Møller point tables."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shlex
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FormatStrFormatter, MaxNLocator, ScalarFormatter
from PIL import Image

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from scripts.analysis.relaxtime.assemble_ward_moller_mott_response import CORE, PROFILE
from scripts.plotting.plot_accessibility import accessibility_record
from scripts.plotting.plot_manifest import (
    build_manifest, generator_record, input_record, output_record,
    runtime_record, sha256_file, write_manifest,
)
from scripts.plotting.plot_quality import export_figure, inspect_export, placement_limits
from scripts.plotting.plot_style import configure_axis_ticks, load_profile, style_context
from scripts.plotting.validate_plot_artifact import validate_manifest_record

STEM = "mott_threshold_local_response"
CURVE = "#245a9b"
REFERENCE = "#4d4d4d"
SLOPE = r"$\Delta\ln(\eta/s)/\Delta\xi$"
DETUNING = r"$\xi-\xi_{\mathrm{Mott}}$ ($10^{-3}$)"
CURVATURE = r"$\sqrt{h}\,I_S''$ ($\mathrm{fm}^2$)"
DISTANCE = r"$h$ ($10^{-6}$)"


def frozen_inputs(folder):
    record = json.loads((folder / "numerical_manifest.json").read_text(encoding="utf-8"))
    if (record["core_commit"], record["profile_sha256"], record["status"]) != (
        CORE, PROFILE, "frozen_analysis"
    ):
        raise ValueError("Figure requires a frozen Ward/Møller numerical record")
    for name, sha in record["files"].items():
        if sha256_file(folder / name) != sha:
            raise ValueError(f"Frozen input hash mismatch: {name}")
    return record


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(enumerate(csv.DictReader(stream), start=2))


def _curve(ax, indexed_rows, *, x, y, scale, panel, source, description, y_scale=1.0):
    vertices = [(scale * float(row[x]), y_scale * float(row[y])) for _, row in indexed_rows]
    if len(vertices) < 2 or not all(math.isfinite(v) for pair in vertices for v in pair):
        raise ValueError(f"Incomplete series for panel {panel}")
    if any(a[0] >= b[0] for a, b in zip(vertices, vertices[1:])):
        raise ValueError(f"Nonincreasing series for panel {panel}")
    ax.plot(*zip(*vertices), color=CURVE, linewidth=1.35)
    return {
        "series_id": panel,
        "panel": panel,
        "state": "computed_derived_series",
        "source_csv": source,
        "source_columns": [x, y],
        "source_rows": [index for index, _ in indexed_rows],
        "point_count": len(vertices),
        "vertices_sha256": hashlib.sha256(
            json.dumps(vertices, separators=(",", ":"), allow_nan=False).encode("ascii")
        ).hexdigest(),
        "support_rule": description,
        "mask_rule": "none within the stated support",
        "x_multiplier": scale,
        "y_multiplier": y_scale,
    }


def _axis(panel, field, unit, label, transform="identity"):
    return dict(panel=panel, field=field, source_unit=unit, display_unit=unit,
                label=label, transform=transform)


def _mott_axis(ax, panel, title, profile, *, show_xlabel=True):
    ax.axvline(0.0, color=REFERENCE, linestyle="--", linewidth=0.85, zorder=1)
    ax.set_xlim(-8.35, 8.35)
    ax.xaxis.set_major_locator(FixedLocator([-8, -4, 0, 4, 8]))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=3))
    ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    ax.set_title(f"({panel})  {title}", fontsize=13, loc="left", pad=7)
    if show_xlabel:
        ax.set_xlabel(DETUNING, fontsize=15)
    configure_axis_ticks(ax, profile)


def render(data, record, profile):
    slopes = read_csv(data / "local_slopes.csv")
    curvature = read_csv(data / "scalar_curvature.csv")
    checks = record["scalar_checks"]
    height = 5.2
    with style_context(profile) as context:
        fig, axs = plt.subplots(2, 2, figsize=(6.75, height))
        context.bind(fig)
        fig.subplots_adjust(left=0.122, right=0.973, bottom=0.113, top=0.929,
                            wspace=0.43, hspace=0.32)
        series, axes = [], []
        for ax, panel, title, meson in (
            (axs[0, 0], "a", r"$\pi$ threshold", "pi"),
            (axs[0, 1], "b", r"$K$ threshold", "K"),
            (axs[1, 0], "c", r"$\sigma$ threshold", "sigma"),
        ):
            rows = sorted((pair for pair in slopes if pair[1]["meson"] == meson),
                          key=lambda pair: float(pair[1]["u_midpoint"]))
            series.append(_curve(ax, rows, x="u_midpoint", y="slope", scale=1000.0,
                                 panel=panel, source="local_slopes.csv",
                                 description="adjacent native log secants on a continuous branch"))
            _mott_axis(ax, panel, title, profile, show_xlabel=panel == "c")
            if panel != "b":
                ax.set_ylabel(SLOPE, fontsize=15)
                axes.append(_axis(panel, "slope", "dimensionless", SLOPE))
            axes.append(_axis(panel, "u_midpoint", "dimensionless", DETUNING, "1000*u_midpoint"))

        parent = axs[1, 0]
        pos = parent.get_position()
        # The inset shows computed residuals only. Its footprint is in the
        # clear part of panel (c), with dimensions defined at final width.
        inset = fig.add_axes([pos.x0 + 0.54*pos.width, pos.y0 + 0.53*pos.height,
                              0.31*pos.width, 0.33*pos.height])
        half_width = record["transport_checks"]["sigma"]["inset"]["display_half_width"]
        subset = [pair for pair in rows if abs(float(pair[1]["u_midpoint"])) <= half_width]
        magnitude = max(abs(float(row["linear_background_residual"])) for _,row in subset)
        exponent = min(0, math.floor(math.log10(magnitude))) if magnitude else 0
        residual_label = rf"Residual $\times 10^{{{-exponent}}}$" if exponent else "Slope residual"
        series.append(_curve(inset, subset, x="u_midpoint", y="linear_background_residual",
                             scale=1000.0, panel="c-inset", source="local_slopes.csv", y_scale=10.0**(-exponent),
                             description="local secants minus a line fitted to the two sidebands"))
        inset.axvline(0.0, color=REFERENCE, linestyle="--", linewidth=0.65, zorder=1)
        inset.set_xlim(-1060*half_width, 1060*half_width)
        inset.xaxis.set_major_locator(FixedLocator([-800*half_width, 0, 800*half_width]))
        inset.xaxis.set_major_formatter(FormatStrFormatter("%g"))
        inset.yaxis.set_major_locator(MaxNLocator(nbins=2, min_n_ticks=2))
        inset.yaxis.set_major_formatter(FormatStrFormatter("%g"))
        inset.yaxis.tick_right()
        inset.set_title(residual_label, fontsize=12.5, loc="left", pad=3)
        configure_axis_ticks(inset, profile)
        inset.tick_params(which="major", length=3, labelsize=12, pad=1.5)
        inset.tick_params(which="minor", length=1.6)
        axes.extend([_axis("c-inset", "u_midpoint", "dimensionless", DETUNING, "1000*u_midpoint"),
                     _axis("c-inset", "linear_background_residual", "dimensionless", residual_label,
                           f"10^{(-exponent)}*linear_background_residual")])

        ax = axs[1, 1]
        usable = sorted((pair for pair in curvature
                         if pair[1]["usable"] == "True"
                         and checks["h_min"] <= float(pair[1]["h"]) <= checks["h_max"]),
                        key=lambda pair: float(pair[1]["h"]))
        series.append(_curve(ax, usable, x="h", y="scaled_curvature", scale=1e6, panel="d",
                             source="scalar_curvature.csv",
                             description="continuous scalar integral on the qualified resonance-side interval"))
        limit = checks["analytic_limit_fm2"]
        ax.axhline(limit, color=REFERENCE, linestyle="--", linewidth=0.85, zorder=1)
        lo = min(limit, *(float(row["scaled_curvature"]) for _,row in usable))
        hi = max(limit, *(float(row["scaled_curvature"]) for _,row in usable))
        pad = 0.095 * (hi-lo)
        ax.set_ylim(lo-pad, hi+pad)
        ax.set_xscale("log")
        ax.set_xlim(0.78 * checks["h_min"] * 1e6, 1.26 * checks["h_max"] * 1e6)
        ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=3))
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.4f"))
        ax.set_title("(d)  Isolated scalar integral", fontsize=13, loc="left", pad=7)
        ax.set_xlabel(DISTANCE, fontsize=16)
        ax.set_ylabel(CURVATURE, fontsize=16)
        configure_axis_ticks(ax, profile)
        axes.extend([_axis("d", "h", "dimensionless", DISTANCE, "1e6*h; logarithmic"),
                     _axis("d", "scaled_curvature", "fm2", CURVATURE)])
        return fig, axes, series


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    record = frozen_inputs(args.data)
    profile = load_profile("candidate_aps_v2")
    fig, axes, series = render(args.data, record, profile)
    args.out.mkdir(parents=True, exist_ok=False)
    outputs, quality = export_figure(fig, args.out / STEM, profile,
                                     formats=["png", "pdf"], intended_width_inches=6.75)
    png = args.out / f"{STEM}.png"
    gray = args.out / f"{STEM}_gray.png"
    with Image.open(png) as image:
        image.convert("L").save(gray, dpi=(profile.dpi, profile.dpi))
    outputs[0]["role"] = "color_review"
    outputs[1]["role"] = "vector_review"
    outputs.append({**output_record(gray, fmt="png", dpi=profile.dpi, vector=False),
                    "inspection": inspect_export(gray), "role": "grayscale_review",
                    "source_color_sha256": sha256_file(png), "conversion": "Pillow L"})
    manifest = build_manifest(
        asset_id="ward_moller_mott_threshold_local_response",
        figure_family="pnjl_transport_appendix_threshold_laws", case_slug=args.out.name,
        figure_mode="audit", semantic_status="frozen_ward_moller_local_response",
        style_profile=profile.profile_id, publication_scope="internal_review",
        generator=generator_record(__file__, command=shlex.join(sys.argv),
                                   runtime=runtime_record({"matplotlib": matplotlib.__version__})),
        inputs=[input_record(args.data / name, role="frozen_numerical_input")
                for name in ("local_slopes.csv", "scalar_curvature.csv", "numerical_manifest.json")]
               + [input_record(profile.path, role="figure_style")],
        outputs=outputs, axes=axes, series=series,
        selection_rule="native secants; resolved scalar stencil interval; no smoothing or extrapolation",
        interpolation_policy="none", connector_policy="forbidden",
        missing_value_policy="refuse missing or nonfinite series",
        validation={"finite": True, "duplicate_keys": True, "support": True, "strict_gate": False},
        rendering={
            "column": "double_column", "figure_size_inches": quality["figure_size_inches"],
            "size_override_reason": "Room for the resolved local inset and four independent axes",
            "quality": quality, "placement_limits": placement_limits(quality, profile, outputs),
            "legend_policy": "declared_geometry_checked", "legend_placements": [], "legend_outside": True,
            "color_route": "color_print_and_online", "review_contract": "prd_review_v1",
            "accessibility": accessibility_record([
                {"value": "computed response", "color": CURVE, "noncolor_encoding": "solid"},
                {"value": "threshold or analytic constant", "color": REFERENCE, "noncolor_encoding": "dashed"},
            ], backgrounds=["#ffffff"], background_scope="opaque white panel backgrounds"),
            "reference_lines": {"a": 0, "b": 0, "c": 0, "c-inset": 0,
                                "d": record["scalar_checks"]["analytic_limit_fm2"]},
            "sigma_inset_fit": record["transport_checks"]["sigma"]["inset"],
            "panel_conditions": {m: {k: c["location"][k] for k in ("T_MeV", "muB_MeV", "xi")}
                                 for m, c in record["centers"].items()},
        },
        calculation_sha=CORE,
        source_run_id=record["action_url"].rsplit("/", 1)[-1],
    )
    # The compiled review uses a vector figure, with no author acceptance asserted.
    manifest.update(status="agent_review", manuscript_eligible=False, current_publication_layer=False)
    errors = validate_manifest_record(manifest)
    if errors:
        print(json.dumps({k: quality[k] for k in (
            "smallest_glyphs", "clipped_text", "text_overlap_pairs", "tick_axes")}, indent=2))
        raise ValueError("; ".join(errors))
    write_manifest(args.out / "plot_manifest.json", manifest)
    print(json.dumps({"out": str(args.out), "minimum_glyph_mm": quality["minimum_capital_numeral_height_mm"],
                      "minimum_curve_pt": quality["minimum_curve_linewidth_pt"]}, indent=2))


if __name__ == "__main__":
    main()
