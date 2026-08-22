from scripts.analysis.pnjl.render_issue130_phase_surface_v4 import (
    DEFAULT_OUTPUT_ROOT,
    render_style_metadata,
)


def test_v4_preserves_cartesian_semantics_and_display_boundaries():
    style = render_style_metadata()

    assert style["style_profile"] == "balanced_cartesian_inner_axes_contrast_v4"
    assert style["projection"] == "orthographic"
    assert style["coordinate_grid"] == "no_pane_or_grid_inner_axis_triad"
    assert style["cartesian_scaffold"]["origin"] == {"mu_q": 0.0, "xi": 0.0, "T": 0.0}
    assert style["cartesian_scaffold"]["pane_grid"] is False
    assert style["surface_projection_walls"] is False
    assert style["mu_zero_marker"]["legend"] is False
    assert style["cep_boundary"]["color"] == "#c2185b"
    assert style["tick_labels"]["grid"] is False


def test_v4_emits_a_sibling_figure_layer_root_without_pdf_contract():
    assert DEFAULT_OUTPUT_ROOT.as_posix().endswith("phase_reference/issue130_phase_reference_v4")
    assert "phase_surface_render_mu_xi_T" in "phase_surface_render_mu_xi_T.png"
