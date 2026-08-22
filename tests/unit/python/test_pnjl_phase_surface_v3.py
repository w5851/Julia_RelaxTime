from scripts.analysis.pnjl.render_issue130_phase_surface_v3 import (
    DEFAULT_OUTPUT_ROOT,
    render_style_metadata,
)


def test_v3_uses_explicit_physical_axes_and_no_panes():
    style = render_style_metadata()

    assert style["style_profile"] == "balanced_orthographic_inner_axes_neutral_light_v3"
    assert style["coordinate_grid"] == "no_pane_or_grid_inner_axis_triad"
    assert style["projection"] == "orthographic"
    assert style["axis_labels"]["xi"] == r"$\xi$ (dimensionless)"
    assert style["cartesian_scaffold"]["origin"] == {"mu_q": 0.0, "xi": 0.0, "T": 0.0}
    assert style["cartesian_scaffold"]["pane_grid"] is False
    assert style["surface_shading"] == "neutral_directional_light_geometry_only"
    assert style["mu_zero_marker"]["legend"] is False
    assert style["legend"]["placement"] == "compact_figure_top_center_outside_data"


def test_v3_is_a_sibling_of_v2():
    assert DEFAULT_OUTPUT_ROOT.as_posix().endswith("phase_reference/issue130_phase_reference_v3")


def test_v3_uses_the_phase_surface_render_basename_and_no_pdf_contract():
    assert "phase_surface_render_mu_xi_T" in "phase_surface_render_mu_xi_T.png"
    assert "pdf" not in render_style_metadata().get("formats", [])
