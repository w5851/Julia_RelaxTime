from scripts.analysis.pnjl.render_issue130_phase_surface_v2 import render_style_metadata


def test_v2_uses_inner_origin_axes_and_marks_mu_zero_intersection_outside_legend():
    style = render_style_metadata()

    assert style["style_profile"] == "inner_origin_cartesian_3d_mu0_intersection_v2"
    assert style["coordinate_grid"] == "no_pane_or_grid_inner_axis_triad"
    assert style["coordinate_axes"] == "inner_origin_cartesian_axes_with_manual_ticks_and_labels"
    assert style["projection"] == "orthographic"
    assert style["mu_zero_marker"]["enabled"] is True
    assert style["mu_zero_marker"]["legend"] is False
    assert style["mu_zero_marker"]["style"] == "thin_dashed_crossover_plane_intersection"
    assert style["cartesian_scaffold"]["origin"] == {"mu_q": 0.0, "xi": 0.0, "T": 0.0}
    assert style["cartesian_scaffold"]["planes"] == []
    assert style["cartesian_scaffold"]["pane_grid"] is False
    assert style["cep_boundary"]["linewidth"] > 2.1
    assert style["cep_boundary"]["halo_linewidth"] > style["cep_boundary"]["linewidth"]


def test_v2_formal_output_is_a_versioned_figure_layer_sibling():
    from scripts.analysis.pnjl.render_issue130_phase_surface_v2 import DEFAULT_OUTPUT_ROOT

    assert DEFAULT_OUTPUT_ROOT.as_posix().endswith("phase_reference/issue130_phase_reference_v2")


def test_v2_uses_the_phase_surface_render_basename():
    assert "phase_surface_render_mu_xi_T" in "phase_surface_render_mu_xi_T.png"
