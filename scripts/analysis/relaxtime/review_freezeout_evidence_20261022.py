"""Recheck this frozen research batch after restoring its raw evidence."""
import csv
from datetime import datetime, timezone
from decimal import Decimal, localcontext
import json
import math
from pathlib import Path
import sys
import tomllib

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from scripts.analysis.relaxtime.freezeout_parameter_envelope import (
    DEFAULT_CATALOG, POLYNOMIAL_FAMILY, assess_models, complete_curve_band,
    load_designs, trajectory, verify_sources,
)
from scripts.analysis.relaxtime.plot_freezeout_comparison import compare_resolutions, load_experiment
from scripts.analysis.relaxtime.plot_freezeout_thermodynamic_envelope import projected_temperature
from scripts.analysis.relaxtime.restore_freezeout_review_evidence import verify
from scripts.analysis.relaxtime.summarize_freezeout_envelope import checked_runs, compatible_physics
from scripts.plotting.plot_manifest import sha256_file
from scripts.plotting.validate_plot_artifact import validate_manifest

BASE = ROOT / "data/outputs/results/relaxtime/meson_density"
RELEASE = BASE / "freezeout_review_20261022"


def current_path(path):
    parts = str(path).replace("\\", "/").split("/")
    if "data" not in parts:
        raise ValueError(f"expected retained data path: {path}")
    return ROOT / Path(*parts[parts.index("data"):])


def table(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def same(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys(), (a.keys(), b.keys())
        for key in a:
            same(a[key], b[key])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b):
            same(x, y)
    elif isinstance(a, float):
        assert math.isclose(a, b, rel_tol=2e-13, abs_tol=1e-14), (a, b)
    else:
        assert a == b, (a, b)


def main():
    _, archive_checks = verify(RELEASE / "archive_manifest.json", ROOT)
    runs, hash_checks, quotients = [], 0, 0
    prefixes = ("freezeout_defaults_20261017", "freezeout_literature_candidates_20261017",
                "freezeout_envelope_20261018", "andronic_envelope_20261018")
    for path in sorted(BASE.glob("*/run.json")):
        if not path.parent.name.startswith(prefixes):
            continue
        models, states, _ = checked_runs([path.parent])
        for points in models.values():
            for point in points.values():
                if point["passed"]:
                    ratio = float(point["K_plus_density_inv_fm3"]) / float(point["pi_plus_density_inv_fm3"])
                    assert math.isclose(point["ratio"], ratio, rel_tol=1e-13, abs_tol=1e-15)
                    assert float(point["background_residual"]) <= 1e-7
                    quotients += 1
                else:
                    assert not math.isfinite(point["ratio"])
        manifest = json.loads((path.parent / "manifest.json").read_text(encoding="utf-8"))
        hash_checks += len(states[0]["source_hashes"]) + len(manifest["outputs"]) + 2
        runs.append(dict(case=path.parent.name, scan_identity=manifest["scan_identity"],
                         point_count=manifest["point_count"], successful_points=manifest["successful_points"],
                         snapshot_hashes_checked=len(states[0]["source_hashes"]),
                         manifest_sha256=sha256_file(path.parent / "manifest.json")))

    experiments = load_experiment(ROOT / "data/inputs/meson_freezeout_experiment_20261018/experimental_ratios.csv")
    families, decimal_scores = [], []
    for prefix in ("freezeout_envelope_20261018", "andronic_envelope_20261018"):
        folder = BASE / (prefix + "_summary")
        summary = json.loads((folder / "summary.json").read_text(encoding="utf-8"))
        assessment_path = current_path(summary["screen_assessment"]["path"])
        assert sha256_file(assessment_path) == summary["screen_assessment"]["sha256"]
        assessment = json.loads(assessment_path.read_text(encoding="utf-8"))
        _, candidates = load_designs([current_path(d["path"]) for d in assessment["design_inputs"]])
        family = summary.get("family", POLYNOMIAL_FAMILY)
        screening, screen_states, _ = checked_runs([current_path(d["directory"]) for d in assessment["run_inputs"]])
        assert len(candidates) == len(screening) == 59
        energies, factor = assessment["primary_energies_GeV"], assessment["near_best_factor"]
        screen_scores = assess_models(screening, experiments, candidates, energies, factor)
        same(screen_scores["profiles"], assessment["profiles"])
        assert screen_scores["ranking"] == assessment["ranking"]
        _, members = complete_curve_band(screening, list(screening), assessment["band_energies_GeV"])
        _, near = complete_curve_band(screening, screen_scores["near_best_profiles"], assessment["band_energies_GeV"])
        assert members == assessment["band_complete_profiles"]
        assert near == assessment["near_best_band_complete_profiles"]
        fine, fine_states, _ = checked_runs([current_path(d["directory"]) for d in summary["fine_inputs"]])
        standard, standard_states, _ = checked_runs([current_path(d["directory"]) for d in summary["standard_inputs"]])
        baselines, baseline_states, _ = checked_runs([current_path(d["directory"]) for d in summary["baseline_inputs"]])
        for state in fine_states + baseline_states:
            compatible_physics(fine_states[0], state, same_resolution=True)
        for state in standard_states + screen_states:
            compatible_physics(fine_states[0], state, same_resolution=False)
        coordinate_error = 0.
        for models in (screening, standard, fine):
            for name, points in models.items():
                for energy, point in points.items():
                    t, mu = trajectory(candidates[name]["values"], energy, family)
                    error = max(abs(t * 1000 - point["T_MeV"]), abs(mu * 1000 - point["muB_MeV"]))
                    assert error < 1e-9
                    coordinate_error = max(coordinate_error, error)
        combined_meta = dict(candidates, **{n: dict(role="published_default") for n in baselines})
        scores = assess_models(dict(baselines, **fine), experiments, combined_meta, energies, factor)
        same(scores["profiles"], summary["metrics"]["profiles"])
        std_scores = assess_models(standard, experiments, candidates, energies, factor)
        same(std_scores["profiles"], summary["standard_metrics"]["profiles"])
        same(compare_resolutions(fine, standard),
             json.loads((folder / "resolution_comparison.json").read_text(encoding="utf-8")))
        assert sorted(fine, key=lambda n: scores["profiles"][n]["rmse"]) == summary["selected_fine_ranking"]
        assert sorted(standard, key=lambda n: std_scores["profiles"][n]["rmse"]) == summary["selected_standard_ranking"]
        projection_checks = 0
        for row in table(BASE / (prefix + "_Tmu_prepared") / "thermodynamic_band.csv"):
            mu = float(row["muB_MeV"]) / 1000
            values = [1000 * projected_temperature(candidates[n]["values"], mu, family) for n in members]
            near_values = [1000 * projected_temperature(candidates[n]["values"], mu, family) for n in near]
            for key, value in (("minimum_T_MeV", min(values)), ("maximum_T_MeV", max(values)),
                               ("near_minimum_T_MeV", min(near_values)), ("near_maximum_T_MeV", max(near_values))):
                assert math.isclose(float(row[key]), value, rel_tol=2e-13, abs_tol=1e-11)
                projection_checks += 1
            assert int(row["member_count"]) == len(members) and int(row["near_member_count"]) == len(near)
        best = summary["best_selected"]
        peak = max(fine[best], key=lambda e: fine[best][e]["ratio"])
        families.append(dict(family=family, candidate_count=len(screening),
            screening_points=sum(len(v) for v in screening.values()),
            screening_successful_points=sum(p["passed"] for v in screening.values() for p in v.values()),
            complete_members=len(members), near_best_complete_members=len(near), best_selected=best,
            fine_rmse=scores["profiles"][best]["rmse"], maximum_coordinate_difference_MeV=coordinate_error,
            standard_points=sum(len(v) for v in standard.values()), fine_points=sum(len(v) for v in fine.values()),
            maximum_standard_fine_relative_difference=summary["maximum_standard_fine_relative_difference"],
            ranking_stable=True, sampled_all_fine_peak_GeV=peak, analytic_projection_values_checked=projection_checks))
        actual = {Decimal(str(e["energy"])): Decimal(str(e["ratio"])) for e in experiments if e["primary_auau"]}
        for profile in (best, "andronic2018"):
            filename = "fine_curves.csv" if profile == best else "baseline_curves.csv"
            values = {Decimal(row["energy_GeV"]): Decimal(row["ratio"]) for row in table(folder / filename)
                      if row["profile"] == profile and row["ratio"]}
            with localcontext() as context:
                context.prec = 50
                rmse = (sum((values[Decimal(str(e))] - actual[Decimal(str(e))]) ** 2 for e in energies)
                        / Decimal(len(energies))).sqrt()
            assert abs(float(rmse) - scores["profiles"][profile]["rmse"]) < 1e-14
            decimal_scores.append(dict(prefix=prefix, profile=profile, decimal_rmse=str(rmse)))

    old = table(ROOT / "data/inputs/meson_freezeout_experiment_20261017/experimental_ratios.csv")
    new = table(ROOT / "data/inputs/meson_freezeout_experiment_20261018/experimental_ratios.csv")
    changes = [dict(id=a["id"], field=k, old=a[k], new=b[k]) for a, b in zip(old, new) for k in a if a[k] != b[k]]
    assert changes == [dict(id="NA49_17.27", field="centrality", old="0-7.2%", new="0-5%")]
    catalog = tomllib.loads(DEFAULT_CATALOG.read_text(encoding="utf-8"))
    frozen_catalog = tomllib.loads((BASE / "freezeout_envelope_20261018_design/literature_envelope_v1.toml").read_text(encoding="utf-8"))
    assert {k: v for k, v in catalog.items() if k != "source_root"} == {k: v for k, v in frozen_catalog.items() if k != "source_root"}
    source_hashes = verify_sources(catalog)
    cases = ("freezeout_defaults_20261017_fine_png_review_r2", "freezeout_envelope_20261018_png_review_r2",
             "freezeout_envelope_20261018_Tmu_png_review_r2", "andronic_envelope_20261018_png_review_r1",
             "andronic_envelope_20261018_Tmu_png_review_r1")
    figures = []
    for case in cases:
        path = ROOT / "data/outputs/figures/relaxtime/meson_density" / case / "plot_manifest.json"
        errors = validate_manifest(path, repo_root=ROOT)
        assert not errors, errors
        manifest = json.loads(path.read_text(encoding="utf-8"))
        assert manifest["delivery_stage"] == "png_review" and not manifest["manuscript_eligible"]
        figures.append(dict(case=case, manifest_sha256=sha256_file(path), outputs_checked=len(manifest["outputs"]),
                            manifest_contract="passed", author_review="pending"))
    report = dict(schema="freezeout_pr_evidence_review_v1", recorded_at_utc=datetime.now(timezone.utc).isoformat(),
                  archive=archive_checks, numerical_run_count=len(runs), run_hash_records_checked=hash_checks,
                  successful_density_quotients_checked=quotients, numerical_runs=runs, families=families,
                  decimal_rmse_checks=decimal_scores, experiment_metadata_changes=changes,
                  portable_source_catalog=dict(source_pdf_hashes_checked=len(source_hashes), semantics_identical_except_source_root=True),
                  figures=figures, scientific_status="diagnostic_only; finite sampled exploration and finite-window resolution review; author PNG acceptance pending")
    (RELEASE / "evidence_review.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: report[k] for k in ("numerical_run_count", "run_hash_records_checked", "successful_density_quotients_checked", "families")}))


if __name__ == "__main__":
    main()
