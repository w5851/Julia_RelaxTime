"""Draw the frozen magnetic PALC convergence comparison; no equation is solved."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
from PIL import Image

from scripts.plotting.plot_accessibility import accessibility_record
from scripts.plotting.plot_manifest import (
    build_manifest, generator_record, input_record, output_record,
    runtime_record, sha256_file, write_manifest,
)
from scripts.plotting.plot_quality import export_figure
from scripts.plotting.plot_series import SeriesStyle, SeriesStyleRegistry
from scripts.plotting.plot_style import configure_axis_ticks, load_profile, style_context


def read_case(path: Path, field: float) -> dict:
    cases = json.loads(path.read_text(encoding="utf-8"))["cases"]
    matching = [c for c in cases if c["eB_GeV2"] == field]
    if len(matching) != 1 or "error" in matching[0]:
        raise ValueError(f"missing or failed magnetic case: {path}, eB={field}")
    return matching[0]


def coordinates(points: list[dict]) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray([p["muB_MeV"] for p in points], dtype=float)
    y = np.asarray([p["rhoB_over_rho0"] for p in points], dtype=float)
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("nonfinite branch coordinate")
    return x, y


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results_dir", type=Path)
    parser.add_argument("figures_dir", type=Path)
    args = parser.parse_args()
    results, figures = args.results_dir.resolve(), args.figures_dir.resolve()
    figures.mkdir(parents=True, exist_ok=True)
    sources = {
        "reference": results / "resolved_mfir_o48_fine_v2.json",
        "p512_0155": results / "current_mfir_p512_e0155.json",
        "p512_04": results / "current_mfir_p512_e04.json",
        "audit_0155": results / "p1024_corrected_production_gap_e0155.json",
        "audit_04": results / "p1024_corrected_production_gap.json",
        "coverage": results / "branch_coverage_e0155.json",
    }
    refs = {b: read_case(sources["reference"], b) for b in (.155, .4)}
    production = {.155: read_case(sources["p512_0155"], .155),
                  .4: read_case(sources["p512_04"], .4)}
    audits = {b: json.loads(sources[name].read_text(encoding="utf-8"))
              for b, name in ((.155, "audit_0155"), (.4, "audit_04"))}
    coverage = json.loads(sources["coverage"].read_text(encoding="utf-8"))
    narrow = coverage["narrow_window"][0]
    if narrow["muB_MeV"] != 872.08:
        raise ValueError("unexpected narrow-window coordinate")
    profile = load_profile("candidate_aps_v2")
    registry = SeriesStyleRegistry("magnetic_palc_convergence_v1", (
        SeriesStyle("independent_mfir", "#0072B2", "-"),
        SeriesStyle("production_p512", "#D55E00", "--"),
        SeriesStyle("corrected_p1024", "#000000", "None", "o"),
    ))
    styles = {s.series_id: s for s in registry.styles}
    with style_context(profile) as context:
        fig, axes = plt.subplots(2, 2, figsize=(6.75, 7.2))
        context.bind(fig)
        fig.subplots_adjust(left=.13, right=.97, top=.85, bottom=.09,
                            hspace=.62, wspace=.40)

        def branch(ax, case, identity, label):
            ax.plot(*coordinates(case["points"]), label=label,
                    **styles[identity].plot_kwargs())

        def corrections(ax, audit):
            x, y = [], []
            for target in audit["targets"]:
                check = next(c for c in target["checks"] if c["p_num"] == 1024)
                x.append(target["muB_MeV"])
                y.append(check["rhoB_over_rho0"])
            ax.plot(x, y, label="p1024 targets", markersize=4.5,
                    markerfacecolor="none", **styles["corrected_p1024"].plot_kwargs())

        for col, b in enumerate((.155, .4)):
            ax = axes[0, col]
            branch(ax, refs[b], "independent_mfir", "Independent MFIR")
            branch(ax, production[b], "production_p512", "p512 PALC")
            corrections(ax, audits[b])
            ax.set(xlim=(820, 1220) if b == .155 else (780, 1220), ylim=(-.12, 4),
                   title=rf"$eB={b:g}\;\mathrm{{GeV}}^2$",
                   xlabel=r"$\mu_B\;(\mathrm{MeV})$")
            ax.set_xticks([900, 1050, 1200])
            ax.set_yticks([0, 2, 4])
            if col == 0:
                ax.set_ylabel(r"$\rho_B/\rho_0$")
            configure_axis_ticks(ax, profile)

        # Reuse the same p512 native vertices and independent native vertices.
        # The lower panels are viewports, not interpolated/new branch samples.
        for col, limits in enumerate(((871.5, 872.6), (871.99, 872.18))):
            ax = axes[1, col]
            branch(ax, refs[.155], "independent_mfir", "Independent MFIR")
            branch(ax, production[.155], "production_p512", "p512 PALC")
            check = next(c for c in narrow["checks"] if c["p_num"] == 1024)
            # Several original roots converge to the same state; display each
            # actual correction, without inventing a continuous p1024 branch.
            corrected = [r for r in check["roots"] if r["converged"]]
            ax.plot([narrow["muB_MeV"]] * len(corrected),
                    [r["audit"]["rhoB_over_rho0"] for r in corrected],
                    label="p1024 targets", markersize=4.5,
                    markerfacecolor="none", **styles["corrected_p1024"].plot_kwargs())
            ax.set(xlim=limits, ylim=(1.23, 1.40),
                   xlabel=r"$\mu_B\;(\mathrm{MeV})$")
            ax.set_yticks([1.25, 1.30, 1.35, 1.40])
            ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
            ax.set_xticks([871.6, 872.0, 872.4] if col == 0 else [872.00, 872.08, 872.15])
            ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f" if col == 0 else "%.2f"))
            if col == 0:
                ax.set_ylabel(r"$\rho_B/\rho_0$")
            else:
                # Stagger the dense terminal tick's label to retain readable
                # precision; tick locations remain the exact displayed values.
                ax.tick_params(axis="x", labelrotation=15)
            configure_axis_ticks(ax, profile)
            ax.set_title(r"$eB=0.155\;\mathrm{GeV}^2$ : local view", pad=8)

        labels = ["Independent MFIR", "p512 PALC", "p1024 targets"]
        fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="upper center",
                   bbox_to_anchor=(.55, .985), ncol=2, frameon=False,
                   handlelength=2.3, columnspacing=1.2)
        outputs, quality = export_figure(fig, figures / "magnetic_palc_convergence", profile,
                                         formats=["png"])
        plt.close(fig)

    color_path = ROOT / outputs[0]["path"]
    gray_path = figures / "magnetic_palc_convergence_gray.png"
    if gray_path.exists():
        raise FileExistsError(gray_path)
    with Image.open(color_path) as image:
        image.convert("L").save(gray_path, dpi=(profile.dpi, profile.dpi))
    gray = output_record(gray_path, fmt="png", dpi=profile.dpi, vector=False, root=ROOT)
    note = ("T=5 MeV; fixed n_max=32, pz_max=40 fm^-1, zeta_num=256. "
            "Top: native p512 and independent MFIR curves; bottom: two viewports of eB=0.155. "
            "At muB=872.08 the p512 five crossings reduce to three unique p1024 roots. "
            "p1024 is local correction only; no continuous p1024 branch or stable-phase selection.")
    selected = [s.series_id for s in registry.styles]
    manifest = build_manifest(
        asset_id="magnetic_palc_convergence", figure_family="magnetic_palc_validation",
        case_slug=figures.name, figure_mode="audit", semantic_status="isolated_convergence_diagnostic",
        style_profile=profile.profile_id, publication_scope="internal_review", root=ROOT,
        calculation_sha="dc183dfa561dc79ef1f8594948b1a1e705679996",
        generator=generator_record(Path(__file__), root=ROOT,
            command=subprocess.list2cmdline([sys.executable, *sys.argv]),
            runtime=runtime_record({"matplotlib": matplotlib.__version__})),
        inputs=[input_record(p, role="frozen_magnetic_validation_results", root=ROOT)
                for p in sources.values()] + [input_record(profile.path, role="plot_style", root=ROOT)],
        axes=[{"field": "muB_MeV", "source_unit": "MeV", "display_unit": "MeV",
               "label": r"$\mu_B\;(\mathrm{MeV})$"},
              {"field": "rhoB_over_rho0", "source_unit": "dimensionless",
               "display_unit": "dimensionless", "label": r"$\rho_B/\rho_0$"}],
        series=[{"series_id": s.series_id, "state": "diagnostic",
                 "support_rule": "native continuation vertices or audited target corrections",
                 "mask_rule": "none; viewport clipping only"} for s in registry.styles],
        outputs=outputs, selection_rule=note,
        interpolation_policy="none; native adjacent continuation vertices in arclength order",
        connector_policy="native branch adjacency only; p1024 markers unconnected",
        missing_value_policy="error; no fill",
        validation={"finite": True, "duplicate_keys": True, "support": True, "strict_gate": False},
        rendering={"column": "double_column", "figure_size_inches": quality["figure_size_inches"],
            "size_override_reason": "two field panels and two narrow viewports at original insertion width",
            "color_route": "undecided_review", "legend_policy": "declared_geometry_checked",
            "legend_outside": True,
            "legend_placements": [{"legend_index": 0, "placement": "outside_axes",
                "host_axes_index": None, "location": "top center", "title": "", "labels": labels,
                "scope": "methods shared by all panels"}],
            "quality": quality, "delivery_stage": "png_review", "vector_delivery_pending": True,
            "review_contract": "prd_review_v1",
            "accessibility": accessibility_record([
                {"value": s.series_id, "color": s.color,
                 "noncolor_encoding": f"{s.linestyle}/{s.marker}"} for s in registry.styles],
                backgrounds=["#ffffff"], background_scope="white axes; no filled regions"),
            "series_styles": registry.record(selected),
            "grayscale_review": {"output": gray, "conversion": "Pillow RGB to L luminance",
                                 "color_sha256": sha256_file(color_path)}},
    )
    manifest["git_commit"] = None
    manifest.update({"manuscript_eligible": False, "current_publication_layer": False,
                     "production_default_changed": False})
    write_manifest(figures / "plot_manifest.json", manifest)
    print(json.dumps({"png": str(color_path), "grayscale": str(gray_path),
                      "minimum_glyph_mm": quality["minimum_capital_numeral_height_mm"],
                      "clipped_text": quality["clipped_text"],
                      "text_overlap_pairs": quality["text_overlap_pairs"]}, indent=2))


if __name__ == "__main__":
    main()
