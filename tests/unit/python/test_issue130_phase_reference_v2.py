from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from scripts.analysis.pnjl.build_issue130_phase_reference_v2 import (
    CALCULATION_SHA,
    PhaseReferenceV2Error,
    build_package,
    sha256,
)


TABLE_ROWS = {
    "boundary": (
        "surface,xi,T_MeV,mu_MeV,rho_hadron,rho_quark,area_residual,"
        "grid_status,grid_unresolved,layer,status,geometry_converged,finite_and_converged\n"
        "maxwell,0.0,100.0,300.0,1.0,2.0,1e-6,ok,false,strict_reference_v1,native,true,true\n"
    ),
    "crossover": (
        "surface,xi,mu_MeV,T_MeV,rho,mu_CEP_proxy_MeV,layer,status,physical_region\n"
        "crossover,0.0,100.0,160.0,1.0,300.0,strict_reference_v1,native,crossover_below_cep\n"
    ),
    "cep": (
        "surface,xi,mu_CEP_proxy_MeV,T_low_MeV,T_high_MeV,T_midpoint_MeV,"
        "layer,status,boundary_mode\n"
        "cep_boundary,0.0,300.0,120.0,120.0625,120.03125,strict_reference_v1,"
        "bracket_preserved,estimated_midpoint\n"
    ),
    "spinodals": (
        "surface,xi,T_MeV,mu_spinodal_hadron_MeV,mu_spinodal_quark_MeV,"
        "rho_spinodal_hadron,rho_spinodal_quark,layer,status\n"
        "spinodal,0.0,100.0,280.0,320.0,1.2,1.8,strict_reference_v1,native\n"
    ),
}


FILE_NAMES = {
    "strict": {
        "boundary": "maxwell_surface_strict_reference_v1.csv",
        "crossover": "crossover_surface_strict_reference_v1.csv",
        "cep": "cep_boundary_strict_reference_v1.csv",
        "spinodals": "spinodal_surface_strict_reference_v1.csv",
    },
    "render": {
        "boundary": "maxwell_surface_render.csv",
        "crossover": "crossover_surface_render.csv",
        "cep": "cep_boundary_render.csv",
    },
    "derived": {
        "boundary": "maxwell_surface_derived_reference_v1.csv",
        "crossover": "crossover_surface_derived_reference_v1.csv",
        "cep": "cep_boundary_derived_reference_v1.csv",
        "spinodals": "spinodal_surface_derived_reference_v1.csv",
    },
}


def _write_source(root: Path, *, bad: str | None = None) -> None:
    root.mkdir(parents=True)
    (root / "manifest.json").write_text(
        json.dumps(
            {
                "schema_version": "pnjl_issue130_phase_reference_import_v1",
                "calculation_sha": CALCULATION_SHA,
                "source_run_id": "source",
                "replay_run_id": "replay",
                "uniform_xi_step": 0.00625,
                "uniform_xi_count": 161,
                "runtime_consumption": False,
            }
        ),
        encoding="utf-8",
    )
    for layer, names in FILE_NAMES.items():
        (root / layer).mkdir(parents=True, exist_ok=True)
        (root / layer / "manifest.json").write_text(
            json.dumps({"layer": layer, "runtime_consumption": False}), encoding="utf-8"
        )
        for table, filename in names.items():
            text = TABLE_ROWS[table]
            if bad == "nonfinite" and table == "boundary":
                text = text.replace("1e-6", "NaN")
            (root / layer / "tables").mkdir(parents=True, exist_ok=True)
            (root / layer / "tables" / filename).write_text(text, encoding="utf-8")
    # v1 render deliberately has no spinodal coordinate table.
    coverage = "surface,xi,coverage_status\nspinodal,0.0,native_support\n"
    if bad == "coverage_nonfinite":
        coverage = coverage.replace("0.0", "NaN")
    (root / "render" / "tables" / "mesh_coverage.csv").write_text(coverage, encoding="utf-8")
    # Give the derived spinodal a non-native row so accepted status is tested.
    derived_spinodal = root / "derived" / "tables" / FILE_NAMES["derived"]["spinodals"]
    derived_spinodal.write_text(
        "surface,xi,T_MeV,mu_spinodal_hadron_MeV,mu_spinodal_quark_MeV,"
        "rho_spinodal_hadron,rho_spinodal_quark,layer,status,interpolation_method,"
        "source_xi_left,source_xi_right,source_axis_low,source_axis_high,xi_gap,reference_write\n"
        "spinodal,0.0,100.0,280.0,320.0,1.2,1.8,interpolated_noncertified,"
        "piecewise_linear_common_support,bilinear_axis_xi_common_support,0.0,0.1,"
        "100.0,100.0,0.1,False\n",
        encoding="utf-8",
    )


def test_v2_materializer_completes_render_and_preserves_status(tmp_path: Path) -> None:
    source = tmp_path / "source"
    output = tmp_path / "output"
    _write_source(source)
    before = sha256(source / "render" / "tables" / "maxwell_surface_render.csv")

    result = build_package(source, output)

    assert result["schema_version"] == "pnjl_issue130_phase_reference_v2"
    assert result["solver_called"] is False
    assert result["runtime_consumption"] is False
    assert sha256(source / "render" / "tables" / "maxwell_surface_render.csv") == before
    assert (output / "render" / "tables" / "spinodal_surface_render.csv").is_file()
    assert (output / "accepted" / "tables" / "spinodal_surface_accepted_phase_map_v1.csv").is_file()
    assert result["constraints"]["derived_public_layer"] is False
    assert result["public_layers"] == ["strict", "render", "accepted"]
    assert result["internal_build_inputs"] == ["data/reference/pnjl/issue130_phase_reference_v1/derived"]

    with (output / "accepted" / "tables" / "spinodal_surface_accepted_phase_map_v1.csv").open(
        encoding="utf-8-sig", newline=""
    ) as handle:
        row = next(csv.DictReader(handle))
    assert row["source_status"] == "interpolated_noncertified"
    assert row["acceptance_status"] == "candidate_pending_author_review"
    assert row["extrapolation"] == "False"


def test_v2_materializer_rejects_nonfinite_source(tmp_path: Path) -> None:
    source = tmp_path / "source"
    _write_source(source, bad="nonfinite")
    with pytest.raises(PhaseReferenceV2Error, match="non-finite"):
        build_package(source, tmp_path / "output")


def test_v2_materializer_rejects_nonfinite_coverage(tmp_path: Path) -> None:
    source = tmp_path / "source"
    _write_source(source, bad="coverage_nonfinite")
    with pytest.raises(PhaseReferenceV2Error, match="non-finite"):
        build_package(source, tmp_path / "output")


def test_v2_materializer_refuses_overwrite(tmp_path: Path) -> None:
    source = tmp_path / "source"
    output = tmp_path / "output"
    _write_source(source)
    output.mkdir()
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        build_package(source, output)
