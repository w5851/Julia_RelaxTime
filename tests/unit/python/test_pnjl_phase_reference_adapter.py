from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.pnjl.phase_reference_adapter import (
    PhaseReferenceContractError,
    build_runtime_view,
    load_phase_reference,
    to_legacy_views,
)
from scripts.pnjl.plot_phase_diagram import DEFAULT_BOUNDARY_PATH, load_candidate_phase_data
from scripts.relaxtime.build_paper_p1_figure_assets import load_phase_candidate_overlay


ROOT = Path(__file__).parents[3]
CANDIDATE = ROOT / "data" / "reference" / "pnjl" / "issue130_phase_reference_v1"


TABLES = {
    "boundary": "maxwell_surface_strict_reference_v1.csv",
    "crossover": "crossover_surface_strict_reference_v1.csv",
    "cep": "cep_boundary_strict_reference_v1.csv",
    "spinodals": "spinodal_surface_strict_reference_v1.csv",
}


def _write_candidate(root: Path, *, duplicate: bool = False, nonfinite: bool = False) -> None:
    layer = root / "strict"
    tables = layer / "tables"
    tables.mkdir(parents=True)
    (root / "manifest.json").write_text(
        json.dumps(
            {
                "schema_version": "pnjl_issue130_phase_reference_import_v1",
                "reference_status": "candidate",
                "runtime_consumption": False,
            }
        ),
        encoding="utf-8",
    )
    (layer / "manifest.json").write_text(json.dumps({"layer": "strict_reference_v1"}), encoding="utf-8")

    area = "NaN" if nonfinite else "1e-6"
    boundary = (
        "surface,xi,T_MeV,mu_MeV,rho_hadron,rho_quark,area_residual,"
        "grid_status,grid_unresolved,layer,status,geometry_converged,finite_and_converged\n"
        f"maxwell,0.0,100.0,300.0,1.0,2.0,{area},ok,false,strict_reference_v1,native,true,true\n"
    )
    if duplicate:
        boundary += boundary.split("\n", 1)[1]
    (tables / TABLES["boundary"]).write_text(boundary, encoding="utf-8")
    (tables / TABLES["crossover"]).write_text(
        "surface,xi,mu_MeV,T_MeV,rho,mu_CEP_proxy_MeV,layer,status,physical_region\n"
        "crossover,0.0,100.0,160.0,1.0,300.0,strict_reference_v1,native,crossover_below_CEP\n",
        encoding="utf-8",
    )
    (tables / TABLES["cep"]).write_text(
        "surface,xi,mu_CEP_proxy_MeV,T_low_MeV,T_high_MeV,T_midpoint_MeV,layer,status,boundary_mode\n"
        "cep_boundary,0.0,300.0,120.0,120.0625,120.03125,strict_reference_v1,bracket_preserved,estimated_midpoint\n",
        encoding="utf-8",
    )
    (tables / TABLES["spinodals"]).write_text(
        "surface,xi,T_MeV,mu_spinodal_hadron_MeV,mu_spinodal_quark_MeV,"
        "rho_spinodal_hadron,rho_spinodal_quark,layer,status\n"
        "spinodal,0.0,100.0,280.0,320.0,1.2,1.8,strict_reference_v1,native\n",
        encoding="utf-8",
    )


def test_candidate_layer_maps_mu_q_and_preserves_cep_bracket(tmp_path: Path) -> None:
    root = tmp_path / "candidate"
    _write_candidate(root)

    bundle = load_phase_reference(root, layer="strict")
    assert bundle.runtime_consumption is False
    assert bundle.tables["boundary"][0]["muq_MeV"] == 300.0
    assert bundle.tables["boundary"][0]["muB_MeV"] == 900.0
    assert bundle.tables["cep"][0]["bracket_width_MeV"] == pytest.approx(0.0625)
    assert bundle.tables["cep"][0]["certified"] is True

    legacy = to_legacy_views(bundle)
    assert legacy["boundary"][0]["mu_transition_MeV"] == 300.0
    assert legacy["cep"][0]["T_CEP_MeV"] == pytest.approx(120.03125)
    assert legacy["cep"][0]["muB_CEP_MeV"] == 900.0

    bracket = to_legacy_views(bundle, cep_mode="bracket")
    assert bracket["cep"][0]["T_CEP_MeV"] is None
    assert bracket["cep"][0]["T_low_MeV"] == 120.0


def test_runtime_view_rejects_unresolved_candidate(tmp_path: Path) -> None:
    root = tmp_path / "candidate"
    _write_candidate(root)
    boundary = root / "strict" / "tables" / TABLES["boundary"]
    boundary.write_text(boundary.read_text(encoding="utf-8").replace("native,true,true", "rho_geometry_not_converged,false,true"), encoding="utf-8")

    bundle = load_phase_reference(root, layer="strict")
    assert bundle.diagnostics["uncertified_rows"] >= 1
    with pytest.raises(PhaseReferenceContractError, match="unresolved/non-certified"):
        to_legacy_views(bundle)
    with pytest.raises(PhaseReferenceContractError, match="runtime view rejected"):
        load_phase_reference(root, layer="strict", allow_runtime=True)


def test_duplicate_and_nonfinite_keys_are_rejected(tmp_path: Path) -> None:
    duplicate_root = tmp_path / "duplicate"
    _write_candidate(duplicate_root, duplicate=True)
    with pytest.raises(PhaseReferenceContractError, match="duplicate"):
        load_phase_reference(duplicate_root)

    nonfinite_root = tmp_path / "nonfinite"
    _write_candidate(nonfinite_root, nonfinite=True)
    with pytest.raises(PhaseReferenceContractError, match="non-finite"):
        load_phase_reference(nonfinite_root)


def test_imported_candidate_is_solver_free_and_has_real_manifest() -> None:
    bundle = load_phase_reference(CANDIDATE, layer="strict")
    assert bundle.manifest["runtime_consumption"] is False
    assert bundle.diagnostics["manifest_sha256"]
    assert bundle.diagnostics["row_counts"]["boundary"] > 0


def test_plot_consumer_requires_explicit_candidate_root(tmp_path: Path) -> None:
    root = tmp_path / "candidate"
    _write_candidate(root)
    boundary, ceps, spinodals, crossovers, diagnostics = load_candidate_phase_data(root)
    assert len(boundary) == len(ceps) == len(spinodals) == len(crossovers) == 1
    assert crossovers[0].T_deconf_MeV is None
    assert diagnostics["runtime_consumption"] is False
    assert DEFAULT_BOUNDARY_PATH == ROOT / "data" / "reference" / "pnjl" / "legacy_phase_reference_v1" / "boundary.csv"


def test_paper_p1_candidate_overlay_preserves_mu_scale_and_status(tmp_path: Path) -> None:
    root = tmp_path / "candidate"
    _write_candidate(root)
    rows = load_phase_candidate_overlay(root, "strict", 3.0)
    first_order = next(row for row in rows if row["kind"] == "first_order")
    crossover = next(row for row in rows if row["kind"] == "crossover")
    cep = next(row for row in rows if row["kind"] == "cep")
    assert first_order["muB_MeV"] == 900.0
    assert crossover["muB_MeV"] == 300.0
    assert cep["muB_MeV"] == 900.0
    assert first_order["candidate_certified"] is True
    assert cep["candidate_status"] == "bracket_preserved"


def test_runtime_view_uses_certified_candidate_then_legacy_fallback(tmp_path: Path) -> None:
    root = tmp_path / "candidate"
    _write_candidate(root)
    boundary = root / "strict" / "tables" / TABLES["boundary"]
    boundary.write_text(
        boundary.read_text(encoding="utf-8").replace("native,true,true", "rho_geometry_not_converged,false,true"),
        encoding="utf-8",
    )
    candidate = load_phase_reference(root, layer="strict")
    legacy_row = dict(candidate.tables["boundary"][0])
    legacy_row["muq_MeV"] = 301.0
    runtime = build_runtime_view(candidate, legacy_tables={"boundary": [legacy_row]})
    assert runtime.source == "candidate"
    assert runtime.diagnostics["runtime_view"] == "certified_candidate_with_legacy_fallback"
    assert runtime.diagnostics["fallback_row_counts"]["boundary"] == 1
    assert runtime.tables["boundary"][0]["source_layer"] == "legacy_fallback"
