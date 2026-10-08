from __future__ import annotations

import copy
import json
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np
import pytest

from scripts.analysis.pnjl.build_phase_diagram_prd import ROOT, MARKERS, build_case, estimate_density, load_sources, make_figure
from scripts.plotting.plot_manifest import input_record
from scripts.plotting.plot_style import load_profile
from scripts.plotting.validate_plot_artifact import validate_manifest_record


def test_density_estimate_uses_nearest_native_pair_without_claiming_a_bracket():
    result = estimate_density(np.array([60., 65.]), np.array([1., 2.]), np.array([4., 3.]), 67.)
    assert result["rho_over_rho0"] == 2.5
    assert result["boundary_temperature_MeV"] == 65
    assert result["interval_semantics"] == "coexistence_pair_not_cep_bracket"
    assert result["creates_numerical_source_row"] is False
    with pytest.raises(ValueError, match="within"):
        estimate_density(np.array([60., 65.]), np.array([1., 2.]), np.array([4., 3.]), 76.)


def test_native_curves_and_accepted_r2_estimates_are_preserved():
    rows, _, _ = load_sources(ROOT)
    with matplotlib.rc_context():
        figure, series, closures, access = make_figure(rows, load_profile("candidate_aps_v2"))
        estimates = [s for s in series if s["state"] == "estimated_density"]
        assert [s["rho_over_rho0"] for s in estimates] == pytest.approx(
            [2.1681788874013033, 2.0026268738857538, 1.9252909326618552, 1.7404335535729727, 1.5206102118261988], abs=1e-12)
        assert len(closures) == 15
        assert all(c["creates_numerical_source_row"] is False for c in closures)
        assert len({p["noncolor_encoding"] for p in access["parameters"]}) == len(MARKERS) == 5
        assert all(ax.get_ylim() == (60., 230.) for ax in figure.axes)
        assert not any(isinstance(c, PolyCollection) for ax in figure.axes for c in ax.collections)
        low = next(s for s in series if s["series_id"] == "low_density_xi_-0.5")
        native = sorted((r for r in rows["boundary"] if float(r["xi"]) == -0.5), key=lambda r: float(r["plot_order_key"]))
        assert low["x"] == [float(r["rho_hadron"]) for r in native]
        assert low["y"] == [float(r["T_MeV"]) for r in native]
        plt.close(figure)


def test_full_two_stage_delivery_and_negative_public_validator(tmp_path):
    if not all(shutil.which(tool) for tool in ("pdfinfo", "pdffonts", "pdfimages")):
        pytest.skip("Poppler is required for actual vector delivery checks")
    png_dir = tmp_path / "review"
    png = build_case(png_dir)
    assert {r["format"] for r in png["outputs"]} == {"png"}
    assert not list(png_dir.glob("*.pdf"))
    receipt = tmp_path / "acceptance.json"
    receipt.write_text(json.dumps({"schema": "plot_png_acceptance_v1", "status": "author_accepted",
        "pdf_export_authorized": True, "author_instruction": "Synthetic integration fixture; not a real author approval.",
        "accepted_manifest": input_record(png_dir / "plot_manifest.json", role="accepted_png_manifest"),
        "accepted_png": input_record(png_dir / "phase_diagram_TmuB_Trho.png", role="accepted_png")}), encoding="utf-8")
    vector = build_case(tmp_path / "vector", acceptance_path=receipt)
    assert vector["outputs"][0]["sha256"] == png["outputs"][0]["sha256"]
    assert vector["manuscript_eligible"] is False
    assert validate_manifest_record(vector) == []
    missing_acceptance = copy.deepcopy(vector)
    missing_acceptance.pop("author_review")
    assert any("author_accepted_png" in error for error in validate_manifest_record(missing_acceptance))
    bad_density = copy.deepcopy(vector)
    estimate = next(s for s in bad_density["series"] if s["state"] == "estimated_density")
    estimate["rho_over_rho0"] += 0.1
    assert any("invalid estimated_density" in error for error in validate_manifest_record(bad_density))
    missing_accessibility = copy.deepcopy(png)
    missing_accessibility["rendering"].pop("accessibility")
    assert any("plot_accessibility_v1" in error for error in validate_manifest_record(missing_accessibility))
