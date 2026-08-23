from scripts.analysis.pnjl.render_issue130_phase_surface_v8 import (
    append_cep_boundary_display,
    build_surface_quads,
    group_by_xi,
)


def _crossover(xi, mu, temperature):
    return {
        "surface": "crossover",
        "xi": str(xi),
        "mu_MeV": str(mu),
        "T_MeV": str(temperature),
        "rho": "1.0",
    }


def _cep(xi, mu, temperature):
    return {
        "surface": "cep_boundary",
        "xi": str(xi),
        "mu_CEP_proxy_MeV": str(mu),
        "T_midpoint_MeV": str(temperature),
    }


def test_cep_endpoint_is_display_only_and_does_not_change_source_rows():
    source = [_crossover(0.0, 0.0, 180.0), _crossover(0.0, 10.0, 170.0)]
    groups, closure = append_cep_boundary_display(
        group_by_xi(source),
        [_cep(0.0, 12.0, 168.0)],
    )

    assert len(source) == 2
    assert len(groups[0.0]) == 3
    assert groups[0.0][-1]["status"] == "boundary_constrained_cep_display_endpoint"
    assert closure[0]["closure_status"] == "added_boundary_constrained_display_endpoint"


def test_surface_quads_use_only_adjacent_common_support():
    rows = [
        _crossover(0.0, 0.0, 180.0),
        _crossover(0.0, 1.0, 179.0),
        _crossover(0.0, 2.0, 178.0),
        _crossover(1.0, 1.0, 170.0),
        _crossover(1.0, 2.0, 169.0),
        _crossover(1.0, 3.0, 168.0),
    ]
    quads, summary = build_surface_quads(
        group_by_xi(rows),
        axis="mu_MeV",
        value_field="T_MeV",
        grid_points=8,
        max_gap=2.0,
    )

    assert quads
    assert summary["xi_pair_count_with_common_support"] == 1
    assert all(1.0 - 1.0e-12 <= vertex[0] <= 2.0 + 1.0e-12 for quad in quads for vertex in quad)


def test_large_local_gap_is_masked_instead_of_connected():
    rows = [
        _crossover(0.0, 0.0, 180.0),
        _crossover(0.0, 10.0, 170.0),
        _crossover(1.0, 0.0, 175.0),
        _crossover(1.0, 10.0, 165.0),
    ]
    quads, summary = build_surface_quads(
        group_by_xi(rows),
        axis="mu_MeV",
        value_field="T_MeV",
        grid_points=8,
        max_gap=5.0,
    )

    assert quads == []
    assert summary["quad_count"] == 0
