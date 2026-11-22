"""Render the Appendix-A local-response composite from frozen Mott data."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter, LogLocator

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from scripts.plotting.plot_manifest import (  # noqa: E402
    build_manifest,
    generator_record,
    input_record,
    output_record,
    runtime_record,
    sha256_file,
    validate_manifest_record,
    write_manifest,
)
from scripts.plotting.plot_quality import export_figure, inspect_export  # noqa: E402
from scripts.plotting.plot_style import configure_axis_ticks, load_profile, style_context  # noqa: E402


CORE = "096199b96d34045bf8e2d9611c2d786d6f92fb8e"
PROFILE = "a26d470cf09730ba79d4a440bf9e3d7bc17c202f59323591edafe4a47d249313"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def f(row: dict[str, str], key: str) -> float:
    return float(row[key])


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _series(vertices, panel, source, label, **extra):
    return {
        "series_id": label,
        "panel": panel,
        "state": "computed_native_series",
        "source_csv": str(source),
        "point_count": len(vertices),
        "vertices_sha256": hashlib.sha256(
            json.dumps(vertices, separators=(",", ":")).encode("utf-8")
        ).hexdigest(),
        "support_rule": "all frozen computed vertices in the declared local window",
        "mask_rule": "none",
        **extra,
    }


def render(data: Path, profile):
    slopes_path = data / "local_slopes.csv"
    curvature_path = data / "scalar_curvature.csv"
    slopes = read_csv(slopes_path)
    curvature = read_csv(curvature_path)

    with style_context(profile) as context:
        fig, axs = plt.subplots(2, 2, figsize=(6.75, 5.20))
        context.bind(fig)
        fig.subplots_adjust(left=0.105, right=0.976, bottom=0.115, top=0.945,
                            wspace=0.29, hspace=0.38)
        series = []

        for ax, meson, panel, title in (
            (axs[0, 0], "pi", "a", r"$\pi$ threshold, $T=138.31$ MeV"),
            (axs[0, 1], "K", "b", r"$K$ threshold, $T=200.09$ MeV"),
        ):
            selected = sorted((row for row in slopes if row["meson"] == meson),
                              key=lambda row: f(row, "u_midpoint"))
            vertices = [(f(row, "u_midpoint"), f(row, "slope")) for row in selected]
            ax.plot([x for x, _ in vertices], [y for _, y in vertices],
                    color="#245a9b", lw=1.25, marker="o", ms=2.7,
                    label=r"$\Delta\ln(\eta/s)/\Delta\xi$")
            ax.axvline(0.0, color="#555555", ls="--", lw=0.8)
            ax.set_xlim(-0.0083, 0.0083)
            ax.set_title(title, fontsize=10.7, pad=5)
            ax.text(0.97, 0.94, f"({panel})", transform=ax.transAxes,
                    ha="right", va="top", fontsize=11)
            ax.xaxis.set_major_formatter(FormatStrFormatter("%.3f"))
            configure_axis_ticks(ax, profile)
            series.append(_series(vertices, panel, slopes_path,
                                  f"{panel}.{meson}.log_eta_s_slope",
                                  observable="Delta log(eta/s) / Delta xi"))

        selected = sorted((row for row in slopes if row["meson"] == "sigma"),
                          key=lambda row: f(row, "u_midpoint"))
        vertices = [(f(row, "u_midpoint"), f(row, "slope")) for row in selected]
        residuals = [(f(row, "u_midpoint"), f(row, "linear_background_residual"))
                     for row in selected]
        ax = axs[1, 0]
        ax.plot([x for x, _ in vertices], [y for _, y in vertices],
                color="#245a9b", lw=1.25, marker="o", ms=2.7)
        ax.axvline(0.0, color="#555555", ls="--", lw=0.8)
        ax.set_xlim(-0.0083, 0.0083)
        ax.set_title(r"$\sigma$ threshold, $T=125.74$ MeV", fontsize=10.7, pad=5)
        ax.text(0.97, 0.94, "(c)", transform=ax.transAxes,
                ha="right", va="top", fontsize=11)
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.3f"))
        configure_axis_ticks(ax, profile)
        series.append(_series(vertices, "c", slopes_path, "c.sigma.log_eta_s_slope",
                              observable="Delta log(eta/s) / Delta xi"))

        inset = ax.inset_axes([0.53, 0.15, 0.43, 0.42])
        inset.plot([x for x, _ in residuals], [y for _, y in residuals],
                   color="#b54a35", lw=1.0, marker="o", ms=1.9)
        inset.axvline(0.0, color="#555555", ls="--", lw=0.65)
        inset.set_xlim(-0.0011, 0.0011)
        inset.set_title("background-subtracted", fontsize=7.3, pad=2)
        inset.tick_params(labelsize=6.5, direction="in", length=2.2)
        inset.xaxis.set_major_formatter(FormatStrFormatter("%.3f"))
        series.append(_series(residuals, "c-inset", slopes_path,
                              "c.sigma.background_residual",
                              observable="background-subtracted local slope"))

        ax = axs[1, 1]
        vertices = [(f(row, "h"), f(row, "scaled_curvature")) for row in curvature]
        ax.plot([x for x, _ in vertices], [y for _, y in vertices],
                color="#245a9b", lw=1.25, marker="o", ms=2.7,
                label=r"$\sqrt{h}\,I_S''$")
        limit = f(curvature[0], "analytic_limit")
        ax.axhline(limit, color="#555555", ls="--", lw=0.8,
                   label="analytic threshold constant")
        ax.set_xscale("log")
        ax.set_xlim(min(x for x, _ in vertices) / 0.7, max(x for x, _ in vertices) * 1.25)
        ax.set_title(r"Isolated scalar threshold integral", fontsize=10.7, pad=5)
        ax.text(0.97, 0.94, "(d)", transform=ax.transAxes,
                ha="right", va="top", fontsize=11)
        ax.set_xlabel(r"$h=|\xi-\xi_{\rm Mott}|$")
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))
        ax.legend(frameon=False, fontsize=7.4, loc="lower left",
                  handlelength=1.5, borderpad=0.1)
        configure_axis_ticks(ax, profile)
        series.append(_series(vertices, "d", curvature_path, "d.scalar.scaled_curvature",
                              observable="sqrt(h) I_S second derivative",
                              analytic_limit=limit))

        for ax in axs.flat:
            ax.set_xlabel(ax.get_xlabel() or r"$\xi-\xi_{\rm Mott}$")
            ax.set_ylabel(ax.get_ylabel() or r"$\Delta\ln(\eta/s)/\Delta\xi$")
        axs[1, 1].set_ylabel(r"$\sqrt{h}\,I_S''$")
        return fig, series


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    profile = load_profile("candidate_aps_v2")
    fig, series = render(args.data, profile)
    stem = args.out / "mott_threshold_local_response"
    export_figure(fig, stem, profile, formats=["png", "pdf"], intended_width_inches=6.75)
    png = stem.with_suffix(".png")
    pdf = stem.with_suffix(".pdf")
    manifest = build_manifest(
        asset_id="ward_moller_mott_threshold_local_response",
        figure_family="pnjl_transport_appendix_threshold_laws",
        case_slug=args.out.name,
        figure_mode="audit",
        semantic_status="new_prescription_computed_samples",
        style_profile=profile.profile_id,
        publication_scope="internal_review",
        generator=generator_record(__file__, command=" ".join(sys.argv),
                                   runtime=runtime_record({"matplotlib": matplotlib.__version__})),
        inputs=[input_record(args.data / name, role="frozen_numerical_input")
                for name in ("local_slopes.csv", "scalar_curvature.csv", "numerical_manifest.json")],
        outputs=[
            {**output_record(png, fmt="png", dpi=600, vector=False),
             "inspection": inspect_export(png), "role": "color_review"},
            {**output_record(pdf, fmt="pdf", dpi=None, vector=True),
             "inspection": inspect_export(pdf), "role": "vector_manuscript_delivery"},
        ],
        axes=[
            {"panel": panel, "observable": observable, "support": "frozen local response"}
            for panel, observable in (("a", "pi log eta/s slope"), ("b", "K log eta/s slope"),
                                      ("c", "sigma log eta/s slope and residual inset"),
                                      ("d", "isolated scalar curvature"))
        ],
        series=series,
        selection_rule="all finite frozen vertices; no smoothing or extrapolation",
        interpolation_policy="none",
        connector_policy="forbidden",
        missing_value_policy="reject missing/nonfinite rows",
        validation={"finite": True, "duplicate_keys": True, "strict_gate": False},
        rendering={"delivery_stage": "png_review", "manuscript_eligible": False,
                   "vector_delivery_pending": False, "column": "double_column",
                   "figure_size_inches": [6.75, 5.20],
                   "reference_lines": ["xi=xi_Mott", "analytic scalar threshold constant"]},
        calculation_sha={"core_commit": CORE, "profile_sha256": PROFILE,
                         "numerical_manifest_sha256": digest(args.data / "numerical_manifest.json")},
    )
    errors = validate_manifest_record(manifest)
    if errors:
        raise ValueError("; ".join(errors))
    write_manifest(args.out / "plot_manifest.json", manifest)
    print(json.dumps({"out": str(args.out), "png_sha256": sha256_file(png),
                      "pdf_sha256": sha256_file(pdf), "series": len(series)}, indent=2))


if __name__ == "__main__":
    main()
