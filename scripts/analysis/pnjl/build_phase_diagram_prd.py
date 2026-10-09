"""Render the frozen manuscript phase map; never run a numerical solver.

PNG review comes first. Vector delivery requires a hash-bound author response.
The nearest-temperature density estimate and explicit closure are display data.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator
import numpy as np
from PIL import Image, ImageOps

from scripts.plotting.plot_accessibility import accessibility_record, validate_accessibility
from scripts.plotting.plot_delivery import validate_png_delivery
from scripts.plotting.plot_manifest import (build_manifest, generator_record, input_record,
    output_record, runtime_record, sha256_file, write_manifest)
from scripts.plotting.plot_quality import export_figure, inspect_export, placement_limits
from scripts.plotting.plot_style import configure_axis_ticks, configure_matplotlib, load_profile, phase_style, resolve_font
from scripts.plotting.validate_plot_artifact import validate_manifest

TAG = "figure4_phase_diagram_prod_v1_c1_p24t8"
SOURCE = Path("data/outputs/results/pnjl/phase_diagram/figure4_phase_diagram_prod_v1")
XI_VALUES = (-0.5, -0.25, 0.0, 0.25, 0.5)
MARKERS = ("o", "s", "^", "D", "v")
TEMPERATURE_LIMITS = (60.0, 230.0)
CONTRACTS = ["config/plotting/candidate_aps_v2.toml", "docs/guides/sop/workflows/figure_production.md",
    *[f"scripts/plotting/{name}.py" for name in ("plot_style", "plot_manifest", "plot_quality",
        "plot_bundle", "plot_provenance", "validate_plot_artifact", "plot_delivery", "plot_accessibility")]]
CAPTION = r"""Chiral phase structure of anisotropic quark matter in the PNJL model in the
(a) $T$--$\mu_B$ and (b) $T$--$\rho/\rho_0$ planes, with
$\rho_0=0.16\,\mathrm{fm}^{-3}$. Colors and symbol shapes identify $\xi$.
Solid lines denote first-order chiral transitions; the two solid branches in
panel (b) bound the low- and high-density coexistence interval. Dashed lines
denote crossover transitions. Filled symbols in panel (a) give the CEP
coordinates from the source data. Open symbols in panel (b) use the same CEP
temperatures, with densities estimated as the mean of the two coexistence
densities at the nearest available temperature below each CEP. The short
joining segments close the coexistence branches and crossover curves onto
these estimates as guides to the eye. AI assistance was used to prepare the
plotting code."""


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def resolve(record: dict, root: Path) -> Path:
    path = Path(record["path"])
    return path if path.is_absolute() else root / path


def load_sources(root: Path) -> tuple[dict, dict, list[Path]]:
    metadata_path = root / SOURCE / "phase_reference_source_manifest.json"
    metadata = read_json(metadata_path)
    expected = {Path(item["path"]).name: item for item in metadata["artifacts"]}
    paths = [root / SOURCE / "reference" / f"{name}_{TAG}.csv" for name in ("boundary", "crossover", "cep")]
    rows = {}
    for name, path in zip(("boundary", "crossover", "cep"), paths):
        if path.stat().st_size != expected[path.name]["bytes"] or sha256_file(path) != expected[path.name]["sha256"]:
            raise ValueError(f"frozen numerical input drift: {path}")
        with path.open(encoding="utf-8-sig", newline="") as stream:
            rows[name] = list(csv.DictReader(line for line in stream if line.strip() and not line.startswith("#")))
    return rows, metadata, [*paths, metadata_path]


def select(rows: list[dict], xi: float, key: str) -> list[dict]:
    selected = sorted((r for r in rows if abs(float(r["xi"]) - xi) < 1e-10), key=lambda r: float(r[key]))
    keys = [float(r[key]) for r in selected]
    if not keys or len(keys) != len(set(keys)) or not np.isfinite(keys).all():
        raise ValueError(f"missing, repeated or nonfinite support for xi={xi}: {key}")
    return selected


def values(rows: list[dict], field: str) -> np.ndarray:
    result = np.asarray([float(r[field]) for r in rows])
    if not np.isfinite(result).all():
        raise ValueError(f"nonfinite source field: {field}")
    return result


def estimate_density(temperature: np.ndarray, low: np.ndarray, high: np.ndarray, cep_temperature: float) -> dict:
    if not np.isfinite([*temperature, *low, *high, cep_temperature]).all() or not np.all(low < high):
        raise ValueError("density estimation requires finite, ordered coexistence pairs")
    nearest = int(np.argmin(np.abs(temperature - cep_temperature)))
    if nearest != len(temperature) - 1 or not 0 < cep_temperature - temperature[nearest] <= 5.0:
        raise ValueError("nearest native coexistence row must be subcritical and within the retained 5 MeV support")
    return {"state": "estimated_density", "estimated_field": "rho_over_rho0",
            "method": "nearest_subcritical_coexistence_mean", "native_density_pair": [float(low[nearest]), float(high[nearest])],
            "rho_over_rho0": float((low[nearest] + high[nearest]) / 2), "T_MeV": float(cep_temperature),
            "boundary_temperature_MeV": float(temperature[nearest]),
            "interval_semantics": "coexistence_pair_not_cep_bracket", "creates_numerical_source_row": False}


def make_figure(rows: dict, profile):
    configure_matplotlib(profile)
    figure, axes = plt.subplots(1, 2, figsize=(6.75, 4.6), sharey=True)
    figure.subplots_adjust(left=0.10, right=0.98, bottom=0.32, top=0.92, wspace=0.14)
    series, closures, encodings = [], [], []
    for index, xi in enumerate(XI_VALUES):
        color, marker = profile.colors[index], MARKERS[index]
        encodings.append({"value": xi, "color": color, "noncolor_encoding": f"CEP marker {marker}"})
        boundary = select(rows["boundary"], xi, "plot_order_key")
        crossover = select(rows["crossover"], xi, "plot_order_key")
        cep = select(rows["cep"], xi, "T_CEP_MeV")
        if len(cep) != 1 or any(r["converged"].lower() != "true" for r in crossover):
            raise ValueError(f"invalid retained CEP/crossover state at xi={xi}")
        temp, low, high = (values(boundary, key) for key in ("T_MeV", "rho_hadron", "rho_quark"))
        mu = 3 * values(boundary, "mu_transition_MeV")
        mu_c, temp_c, rho_c = 3 * values(crossover, "mu_MeV"), values(crossover, "T_crossover_MeV"), values(crossover, "rho")
        cep_mu, cep_temp = float(cep[0]["muB_CEP_MeV"]), float(cep[0]["T_CEP_MeV"])
        estimate = estimate_density(temp, low, high, cep_temp)
        if not np.isfinite(cep_mu) or np.any(mu_c > cep_mu + 1e-6):
            raise ValueError("crossover support exceeds source CEP chemical potential")
        if min(temp.min(), temp_c.min(), cep_temp) < 60 or max(temp.max(), temp_c.max(), cep_temp) > 230:
            raise ValueError("temperature limits hide a source point")
        curves = [(0, mu, temp, temp, "first_order", "first_order", 5.0),
                  (1, low, temp, temp, "low_density", "first_order", 5.0),
                  (1, high, temp, temp, "high_density", "first_order", 5.0),
                  (0, mu_c, temp_c, mu_c, "crossover_a", "crossover", float(np.median(np.diff(mu_c))) * 1.5),
                  (1, rho_c, temp_c, mu_c, "crossover_b", "crossover", float(np.median(np.diff(mu_c))) * 1.5)]
        for panel, x, y, support, name, state, max_step in curves:
            segments = np.split(np.arange(len(x)), np.flatnonzero(np.diff(support) > max_step + 1e-7) + 1)
            for indices in segments:
                line, = axes[panel].plot(x[indices], y[indices], color=color, linewidth=1.2,
                    linestyle=phase_style(profile, state)["linestyle"], zorder=3)
                assert np.array_equal(line.get_xdata(), x[indices]) and np.array_equal(line.get_ydata(), y[indices])
            series.append({"series_id": f"{name}_xi_{xi:g}", "state": state, "panel": panel, "xi": xi,
                "row_count": len(x), "segment_count": len(segments), "x": x.tolist(), "y": y.tolist(),
                "support_rule": "native frozen rows at exact xi; source plotting order retained",
                "mask_rule": f"reject nonfinite and duplicate support; break support increments above {max_step:.12g}"})
        endpoint = [estimate["rho_over_rho0"], cep_temp]
        for branch, start, state in (("low_density", [float(low[-1]), float(temp[-1])], "first_order"),
                                     ("high_density", [float(high[-1]), float(temp[-1])], "first_order"),
                                     ("crossover", [float(rho_c[-1]), float(temp_c[-1])], "crossover")):
            x, y = np.asarray([start, endpoint]).T
            axes[1].plot(x, y, color=color, linewidth=1.2, linestyle=phase_style(profile, state)["linestyle"], zorder=3)
            record = {"series_id": f"closure_{branch}_xi_{xi:g}", "state": "connector", "parent_phase": state,
                "panel": 1, "xi": xi, "branch": branch, "row_count": 2, "start": start, "end": endpoint,
                "support_rule": "explicit display closure to the estimated density; no numerical row",
                "mask_rule": "other native gaps retained", "creates_numerical_source_row": False}
            series.append(record)
            closures.append(record)
        axes[0].scatter([cep_mu], [cep_temp], s=36, marker=marker, color=color, edgecolor=color, linewidth=0.7, zorder=5)
        axes[1].scatter([estimate["rho_over_rho0"]], [cep_temp], s=36, marker=marker,
                        facecolor="white", edgecolor=color, linewidth=1.1, zorder=5)
        series.append({"series_id": f"cep_xi_{xi:g}", "state": "cep_confirmed", "xi": xi, "panel": 0,
                       "muB_MeV": cep_mu, "T_MeV": cep_temp, "support_rule": "direct retained CEP table coordinates",
                       "mask_rule": "reject nonfinite"})
        series.append({**estimate, "series_id": f"estimated_density_xi_{xi:g}", "xi": xi, "panel": 1,
                       "support_rule": "nearest native subcritical coexistence pair",
                       "mask_rule": "not a certified CEP-density bracket"})
    for index, ax in enumerate(axes):
        ax.set_ylim(*TEMPERATURE_LIMITS)
        ax.yaxis.set_major_locator(FixedLocator([60, 100, 140, 180, 220]))
        ax.set_xlim(0, 1250 if index == 0 else 3.15)
        ax.xaxis.set_major_locator(FixedLocator([0, 400, 800, 1200] if index == 0 else [0, 1, 2, 3]))
        ax.set_xlabel(r"$\mu_B\;(\mathrm{MeV})$" if index == 0 else r"$\rho/\rho_0$")
        ax.text(0.02, 1.025, f"({chr(97 + index)})", transform=ax.transAxes, fontweight="bold", va="bottom")
        configure_axis_ticks(ax, profile)
    axes[0].set_ylabel(r"$T\;(\mathrm{MeV})$")
    figure.legend(handles=[Line2D([], [], color=profile.colors[i], marker=MARKERS[i], linestyle="None",
        markersize=6, label=rf"$\xi={xi:g}$") for i, xi in enumerate(XI_VALUES)], loc="lower center",
        bbox_to_anchor=(0.53, 0.11), ncol=5, handlelength=1.2, columnspacing=0.9, handletextpad=0.35)
    handles = [Line2D([], [], color="#252525", linewidth=1.2, label="First order"),
        Line2D([], [], color="#252525", linewidth=1.2, linestyle="--", label="Crossover"),
        Line2D([], [], color="#252525", marker="o", linestyle="None", markersize=6, label="CEP (a)"),
        Line2D([], [], color="#252525", marker="o", markerfacecolor="white", markeredgewidth=1.1,
               linestyle="None", markersize=6, label=r"$\rho_{\mathrm{CEP}}$ est. (b)")]
    figure.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.53, 0.02), ncol=4,
                  handlelength=1.25, columnspacing=0.75, handletextpad=0.45)
    access = accessibility_record(encodings, backgrounds=["#FFFFFF"],
        background_scope="opaque curves on white axes; no region fills or overlapping translucent backgrounds")
    if validate_accessibility(access):
        raise ValueError(validate_accessibility(access))
    return figure, series, closures, access


def metadata_directory(out: Path, *, root: Path = ROOT) -> Path:
    """Keep captions and provenance outside the image-only figures tree."""
    root, out = root.resolve(), out.resolve()
    try:
        relative = out.relative_to(root / "data/outputs/figures")
    except ValueError:
        if out.is_relative_to(root / "data/outputs"):
            raise ValueError("repository image output must be under data/outputs/figures")
        return out.with_name(out.name + "__metadata")
    return root / "data/outputs/results" / relative


def build_case(out: Path, *, root: Path = ROOT, acceptance_path: Path | None = None) -> dict:
    root, out = root.resolve(), out.resolve()
    metadata = metadata_directory(out, root=root)
    for directory in (out, metadata):
        if directory.exists():
            raise FileExistsError(f"refusing to overwrite case: {directory}")
    accepted = acceptance = None
    if acceptance_path is not None:
        acceptance = read_json(acceptance_path)
        accepted = read_json(resolve(acceptance["accepted_manifest"], root))
        probe = copy.deepcopy(accepted)
        probe["rendering"].update(delivery_stage="vector_delivery", vector_delivery_pending=False)
        probe["author_review"] = {"status": "author_accepted_png", "record": input_record(acceptance_path, role="author_png_acceptance", root=root)}
        issues = validate_png_delivery(probe, root=root)
        if issues:
            raise ValueError(issues)
    rows, source, input_paths = load_sources(root)
    signature = {"generator_sha256": sha256_file(Path(__file__)),
                 "input_sha256": {p.relative_to(root).as_posix(): sha256_file(p) for p in input_paths},
                 "contract_sha256": {p: sha256_file(root / p) for p in CONTRACTS}}
    if accepted is not None and signature != accepted.get("frozen_render_signature"):
        raise ValueError("frozen rendering inputs/code differ from the accepted PNG")
    profile = load_profile(root / "config/plotting/candidate_aps_v2.toml")
    with matplotlib.rc_context():
        figure, series, closures, access = make_figure(rows, profile)
        out.mkdir(parents=True)
        outputs, quality = export_figure(figure, out / "phase_diagram_TmuB_Trho", profile, formats=("png",))
        outputs[0]["role"] = "color_figure"
        if accepted is not None:
            if outputs[0]["sha256"] != acceptance["accepted_png"]["sha256"]:
                raise ValueError("rendered PNG differs from the accepted PNG; PDF was not generated")
            pdf, _ = export_figure(figure, out / "phase_diagram_TmuB_Trho", profile, formats=("pdf",))
            outputs.extend(pdf)
        plt.close(figure)
    gray_path = out / "phase_diagram_TmuB_Trho_gray.png"
    with Image.open(out / "phase_diagram_TmuB_Trho.png") as image:
        ImageOps.grayscale(image.convert("RGB")).save(gray_path, dpi=(profile.dpi, profile.dpi))
    gray = output_record(gray_path, fmt="png", dpi=profile.dpi, vector=False)
    gray.update(role="grayscale_review", inspection=inspect_export(gray_path),
                source_color_sha256=outputs[0]["sha256"], conversion="Pillow ImageOps.grayscale, ITU-R 601 luminance")
    outputs.append(gray)
    snapshot = metadata / "provenance/source_snapshot"
    inputs = [input_record(path, role="calculation_result" if path.suffix == ".csv" else "calculation_source_manifest", root=root) for path in input_paths]
    for relative in CONTRACTS:
        target = snapshot / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(root / relative, target)
        inputs.append(input_record(target, role="plotting_contract", root=root))
    generator_copy = snapshot / "build_phase_diagram_prd.py"
    shutil.copyfile(Path(__file__), generator_copy)
    placements = [{"legend_index": item["legend_index"], "host_axes_index": item["host_axes_index"],
        "title": item["title"], "labels": item["labels"], "placement": "outside_axes", "location": "below panels",
        "scope": "xi color/shape mapping in both panels" if i == 0 else "phase lines and direct/estimated CEP coordinates"}
        for i, item in enumerate(quality["legend_layout"])]
    stage = "vector_delivery" if accepted is not None else "png_review"
    manifest = build_manifest(asset_id=f"pnjl.phase.prd.{stage}", figure_family="pnjl_phase_diagram", case_slug=out.name,
        figure_mode="audit", semantic_status="phase_structure_with_explicit_estimated_density_and_display_closure",
        style_profile=profile.profile_id, publication_scope="internal_review",
        generator=generator_record(generator_copy, command=" ".join(sys.argv), root=root,
            runtime=runtime_record({"matplotlib": matplotlib.__version__, "font": resolve_font(profile)})), inputs=inputs,
        axes=[{"field": "muB_MeV", "source_unit": "MeV", "display_unit": "MeV", "label": r"$\mu_B\;(\mathrm{MeV})$", "transform": "muB=3*mu_q; source muB_CEP is direct", "limits": [0, 1250]},
              {"field": "rho_over_rho0", "source_unit": "dimensionless", "display_unit": "dimensionless", "label": r"$\rho/\rho_0$", "transform": "identity; input densities already normalized by rho0=0.16 fm^-3", "limits": [0, 3.15]},
              {"field": "T_MeV", "source_unit": "MeV", "display_unit": "MeV", "label": r"$T\;(\mathrm{MeV})$", "transform": "identity", "limits": list(TEMPERATURE_LIMITS)}],
        series=series, outputs=outputs, selection_rule="five exact xi slices of the retained frozen phase tables",
        interpolation_policy="native piecewise-linear vertices; no resampling or smoothing",
        connector_policy="three explicit density closure segments per xi, authorized display only",
        missing_value_policy="native gaps retained outside the declared endpoint closure; no new source row",
        validation={"finite": True, "duplicate_keys": True, "support": True, "strict_gate": False,
                    "frozen_source_hashes_checked": True,
                    "native_series_checked": sum("x" in item and "y" in item for item in series)},
        rendering={"column": "double_column", "figure_size_inches": [6.75, 4.6], "bbox_inches": "fixed",
            "delivery_stage": stage, "vector_delivery_pending": stage == "png_review", "review_contract": "prd_review_v1",
            "legend_policy": "declared_geometry_checked", "legend_placements": placements, "legend_outside": True,
            "color_route": "undecided_review", "quality": quality, "accessibility": access,
            "parameter_encoding": "color and distinct CEP shape; filled (a), open density estimate (b); phase uses solid/dashed",
            "coexistence_fill": "none; all native branches and explicitly declared CEP closure are retained"},
        calculation_sha=source["source_head_sha"], source_run_id=source["source_run_url"], root=root)
    metadata_path = metadata.relative_to(root).as_posix() if metadata.is_relative_to(root) else str(metadata)
    manifest.update(manuscript_eligible=False, current_publication_layer=False, frozen_render_signature=signature,
        metadata_directory=metadata_path,
        derived_display_geometry=closures, caption_parameters={"caption_latex": CAPTION, "rho0_fm3inv": 0.16,
            "xi_values": list(XI_VALUES), "xi_markers": list(MARKERS), "estimate": "density only; nearest subcritical coexistence mean",
            "closure": "explicit display guide", "filled_region": "none"},
        placement_limits=placement_limits(quality, profile, outputs), source_qualification={"status": "inherited_without_promotion",
            "verdict": source["verdict"], "residual_risks": source["residual_risks"], "new_numerical_computation": False})
    if accepted is not None:
        accepted_copy = metadata / "provenance/png_review_manifest.json"
        shutil.copyfile(resolve(acceptance["accepted_manifest"], root), accepted_copy)
        acceptance = {**acceptance, "accepted_manifest": input_record(accepted_copy, role="accepted_png_manifest", root=root),
                      "accepted_png": input_record(out / "phase_diagram_TmuB_Trho.png", role="accepted_png", root=root)}
        receipt = metadata / "provenance/png_acceptance.json"
        receipt.write_text(json.dumps(acceptance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        manifest["author_review"] = {"status": "author_accepted_png", "record": input_record(receipt, role="author_png_acceptance", root=root)}
        manifest["inputs"].extend([manifest["author_review"]["record"], acceptance["accepted_manifest"]])
    (metadata / "caption.tex").write_text(CAPTION + "\n", encoding="utf-8")
    limits = manifest["placement_limits"]
    (metadata / "caption_handoff.md").write_text(
        "# PNJL 相图交接\n\n" + f"阶段：{stage}。画布 171.45 × 116.84 mm；最低插入宽度 {limits['minimum_width_inches'] * 25.4:.3f} mm。\n\n"
        "图注见 caption.tex。实心符号为源表 CEP 坐标；空心符号只估计 CEP 密度。两支共存密度来自更低温度，不是 CEP 密度 bracket。闭合线段仅用于显示；无共存区叠加填色。\n\n"
        "相图来源资格保持不变；正文采用依据作者另行授权。数值表、输运图和 current 指针不由本生成器修改。\n", encoding="utf-8")
    (metadata / "README.md").write_text(f"# PNJL 相结构图 PRD 修订\n\n阶段：{stage}。\n\n单图合同见对应 figures case 的 plot_manifest.json；图注和插入限制见 caption_handoff.md。"
        "彩色和灰度图保持相同尺寸。15 条显示闭合段与经验密度估计单独记录，原生顶点与 gap 保留。"
        "provenance 保存冻结绘图源码；矢量阶段还保存已接受的 PNG manifest 和接受记录。\n", encoding="utf-8")
    for path in input_paths:
        if sha256_file(path) != signature["input_sha256"][path.relative_to(root).as_posix()]:
            raise ValueError("numerical source changed during plotting")
    write_manifest(out / "plot_manifest.json", manifest)
    errors = validate_manifest(out / "plot_manifest.json", repo_root=root)
    (metadata / "validation.json").write_text(json.dumps({"passed": not errors, "errors": errors,
        "author_visual_review": "accepted_png" if accepted else "pending", "stage": stage}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if errors:
        raise ValueError(errors)
    return manifest


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--acceptance", type=Path, help="explicit author acceptance of a frozen PNG manifest")
    args = parser.parse_args(argv)
    result = build_case(args.output_dir, root=args.repo_root, acceptance_path=args.acceptance)
    print(json.dumps({"case": str(args.output_dir), "stage": result["rendering"]["delivery_stage"],
                      "metadata_directory": result["metadata_directory"],
                      "minimum_glyph_mm": result["rendering"]["quality"]["minimum_capital_numeral_height_mm"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
