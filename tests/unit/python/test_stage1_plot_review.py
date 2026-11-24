from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from PIL import Image
import pytest

from scripts.analysis.relaxtime.plot_ward_moller_manuscript import transport_inputs
from scripts.pnjl.phase_reference_adapter import (
    PhaseReferenceContractError, estimate_c2_endpoint_density, load_c2_review_sources,
)
from scripts.plotting.plot_manifest import sha256_file
from scripts.plotting.plot_style import load_profile
from scripts.plotting.publication_review import (
    figure_context, grayscale_review_record, write_sop_snapshot,
)
from scripts.plotting.transport_panels import render_panel

ROOT = Path(__file__).resolve().parents[3]


def test_review_companions_keep_dimensions_and_bind_source_bytes(tmp_path):
    color = tmp_path / "phase.png"
    image = Image.new("RGBA", (2, 1))
    image.putdata([(255, 0, 0, 255), (0, 0, 255, 255)])
    image.save(color, dpi=(600, 600))
    record = grayscale_review_record(color)
    with Image.open(tmp_path / "phase_gray.png") as gray:
        assert gray.mode == "L"
        assert gray.size == (2, 1)
        assert tuple(gray.getdata()) == (76, 29)
        assert gray.info["dpi"] == pytest.approx((600, 600), abs=0.02)
    assert record["source_color_sha256"] == sha256_file(color)
    assert record["sha256"] == sha256_file(tmp_path / "phase_gray.png")

    sop = tmp_path / "SOP.md"
    sop.write_bytes(b"# SOP\r\n\r\nPNG review\r\n")
    case = tmp_path / "results"
    snapshot = write_sop_snapshot(case, source=sop)
    assert snapshot.read_bytes() == b"# SOP\n\nPNG review\n"
    assert write_sop_snapshot(case, source=sop) == snapshot
    sop.write_text("another revision", encoding="utf-8")
    with pytest.raises(ValueError, match="existing figure SOP snapshot differs"):
        write_sop_snapshot(case, source=sop)


def test_phase_c2_uses_native_rows_and_estimates_only_from_supported_pair():
    root = ROOT / "data/reference/pnjl/issue130_phase_reference_v2"
    rows, sources, _ = load_c2_review_sources(root, (-0.5, 0.5))
    assert all({float(row["xi"]) for row in group} == {-0.5, 0.5} for group in rows.values())
    assert len(rows["cep"]) == 2
    assert len(sources) == 4
    record = estimate_c2_endpoint_density([100, 102], [0.3, 0.4], [1.3, 1.2], 103)
    assert record["native_density_pair"] == [0.4, 1.2]
    assert record["rho_over_rho0"] == pytest.approx(0.8)
    assert record["creates_numerical_source_row"] is False
    with pytest.raises(PhaseReferenceContractError, match="subcritical"):
        estimate_c2_endpoint_density([100, 102], [0.3, 0.4], [1.3, 1.2], 102)


def test_transport_native_vertices_do_not_join_first_order_gap():
    inputs = ROOT / "data/outputs/results/relaxtime/transport/ward_moller/ward_moller_sphere_20261121_v1"
    data = transport_inputs(inputs / "mode_a_grid.csv", inputs / "mode_b_grid.csv")
    profile = load_profile("candidate_aps_v2")
    with figure_context(profile) as context:
        figure, axis = plt.subplots()
        context.bind(figure)
        try:
            spec = render_panel(axis, "mode_a", "muB900.0", "tau_u", data, profile)
            branch = [line for line in axis.lines if line.get_gid() == "mode_a.alpha1.0"]
            assert [endpoint["xi"] for endpoint in spec["endpoints"]] == [-0.003, 0.003]
            assert len(branch) == 2
            assert branch[0].get_xdata()[-1] == pytest.approx(-0.003)
            assert branch[1].get_xdata()[0] == pytest.approx(0.003)
            assert spec["axis_scale"] == "log"
        finally:
            plt.close(figure)
