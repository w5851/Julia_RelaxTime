from scripts.analysis.pnjl.render_issue130_phase_surface_v5 import (
    DEFAULT_OUTPUT_ROOT,
    render_style_metadata,
)


def test_v5_equal_aspect_and_depth_contract():
    style = render_style_metadata()

    assert style["style_profile"] == "balanced_cartesian_equal_aspect_lit_v5"
    assert style["projection"] == "orthographic"
    assert style["box_aspect"] == [1.0, 1.0, 1.0]
    assert style["view"] == {"elevation_deg": 30.0, "azimuth_deg": -45.0}
    assert style["surface_projection_walls"] is False
    assert style["cartesian_scaffold"]["pane_grid"] is False
    assert style["mu_zero_marker"]["legend"] is False
    assert style["cep_boundary"]["color"] == "#c2185b"
    assert style["formats"] == ["png", "svg"]
    assert style["pdf_emitted"] is False


def test_v5_uses_the_uniform_phase_surface_basename():
    assert DEFAULT_OUTPUT_ROOT.as_posix().endswith("phase_reference/issue130_phase_reference_v5")
    assert "phase_surface_render_mu_xi_T" in "phase_surface_render_mu_xi_T.svg"
