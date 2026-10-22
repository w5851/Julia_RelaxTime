"""Behavior of explicit identities, scoped styles, and declared font policies."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import replace
from pathlib import Path
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from PIL import Image
import pytest

from scripts.plotting.plot_manifest import sha256_file
from scripts.plotting.plot_quality import export_figure, measure_figure
from scripts.plotting.plot_series import SeriesStyle, SeriesStyleRegistry, validate_series_styles
from scripts.plotting.plot_style import configure_matplotlib, load_profile, resolve_font, style_context
from scripts.analysis.relaxtime.plot_freezeout_comparison import CURVE_STYLES
from scripts.analysis.relaxtime.plot_freezeout_thermodynamic_envelope import curve_styles


def theme(font="DejaVu Serif", math="stix", **changes):
    base = load_profile("candidate_aps_v2")
    data = deepcopy(base.data)
    data.update(font_policy="error", font_candidates=[font], mathtext_fontset=math, **changes)
    return replace(base, data=data)


def test_registered_styles_survive_subsets_and_reordering_of_actual_lines():
    profile = theme()
    registry = SeriesStyleRegistry("test", tuple(SeriesStyle.from_profile(profile, name, index)
        for index, name in enumerate(("xi=-0.2", "xi=0", "xi=0.2"))))
    expected = {style.series_id: style.record() for style in registry.styles}
    for selected in (("xi=-0.2", "xi=0", "xi=0.2"), ("xi=0.2", "xi=-0.2"), ("xi=0",)):
        with style_context(profile) as context:
            figure, ax = plt.subplots()
            context.bind(figure)
            try:
                for style in registry.select(selected):
                    ax.plot([1, 2, 3], [1, float("nan"), 2], **style.plot_kwargs())
                evidence = measure_figure(figure, intended_width_inches=6.75)["series_styles"]
                assert evidence == [expected[name] for name in selected]
                assert validate_series_styles(registry.record(selected),
                    [{"series_id": name} for name in selected], evidence) == []
                evidence[0]["color"] = "#000000"
                assert validate_series_styles(registry.record(selected),
                    [{"series_id": name} for name in selected], evidence)
            finally:
                plt.close(figure)


@pytest.mark.parametrize("bad", ["", " name", "name ", 0, None])
def test_invalid_identities_fail_before_rendering(bad):
    with pytest.raises(ValueError, match="series_id"):
        SeriesStyle(bad, "black")


def test_duplicate_unknown_and_exhausted_style_slots_require_explicit_handling():
    style = SeriesStyle("first", "black")
    with pytest.raises(ValueError, match="duplicate"):
        SeriesStyleRegistry("test", (style, style))
    with pytest.raises(ValueError, match="distinct"):
        SeriesStyleRegistry("test", (style, SeriesStyle("second", "black")))
    registry = SeriesStyleRegistry("test", (style,))
    for selected, message in ((["first", "first"], "duplicate"), (["unknown"], "unregistered"), ("first", "sequence")):
        with pytest.raises(ValueError, match=message):
            registry.select(selected)
    profile = theme()
    for index in (-1, True, 100):
        with pytest.raises(ValueError, match="style_index"):
            SeriesStyle.from_profile(profile, "overflow", index)
    with pytest.raises(ValueError, match="linestyle"):
        SeriesStyle.from_profile(profile, "requires-line", len(profile.data["palette"]["parameter_linestyles"]))
    assert SeriesStyle.from_profile(profile, "declared", 5, linestyle=":").linestyle == ":"
    assert SeriesStyleRegistry("explicit", (SeriesStyle("beyond-palette", "#333333", "--", "D"),)).select(["beyond-palette"])


def test_adopted_freezeout_families_keep_their_scientific_encodings():
    assert CURVE_STYLES.select(["andronic2018"])[0] == CURVE_STYLES.select(["cleymans2006", "andronic2018"])[1]
    assert CURVE_STYLES.select(["andronic2018"])[0].plot_kwargs() == {
        "color": "#d55e00", "linestyle": "--", "marker": "", "gid": "andronic2018"}
    registry = curve_styles({"best_selected": "frozen_candidate_42"})
    assert registry.select(["andronic2018"])[0].marker == "s"
    assert registry.select(["frozen_candidate_42"])[0].marker == "o"
    with pytest.raises(ValueError, match="duplicate"):
        curve_styles({"best_selected": "andronic2018"})


def create_figure(profile):
    with style_context(profile) as context:
        figure, ax = plt.subplots(figsize=(6.75, 4.6))
        context.bind(figure)
        ax.plot([0, 1, 2], [1, 2, 1])
        ax.set_xlabel(r"$T_2\;(\mathrm{MeV})$")
        ax.set_ylabel("Value 2")
    return figure


def test_nested_themes_and_failed_exit_restore_all_caller_rcparams():
    first, second = theme(), theme("DejaVu Sans", "cm", font_size_pt=18)
    with matplotlib.rc_context({"font.size": 7, "font.weight": "bold", "text.usetex": True,
                                "mathtext.fontset": "cm", "savefig.transparent": True}):
        before = dict(matplotlib.rcParams)
        with pytest.raises(RuntimeError, match="failure"):
            with style_context(first):
                first_settings = dict(matplotlib.rcParams)
                assert matplotlib.rcParams["text.usetex"] is False
                assert matplotlib.rcParams["font.weight"] == "normal"
                with style_context(second):
                    assert matplotlib.rcParams["font.family"] == ["DejaVu Sans"]
                assert dict(matplotlib.rcParams) == first_settings
                raise RuntimeError("failure")
        assert dict(matplotlib.rcParams) == before


def test_delayed_measurement_and_png_save_keep_the_creation_theme(tmp_path):
    first, second = theme(), theme("DejaVu Sans", "cm", font_size_pt=18)
    figure = create_figure(first)
    other = create_figure(second)
    try:
        original, original_quality = export_figure(figure, tmp_path / "before", first, formats=["png"])
        with matplotlib.rc_context({"mathtext.default": "regular"}), style_context(second):
            export_figure(other, tmp_path / "other", second, formats=["png"])
            foreign = dict(matplotlib.rcParams)
            repeated, repeated_quality = export_figure(figure, tmp_path / "after", first, formats=["png"])
            assert dict(matplotlib.rcParams) == foreign
            assert repeated_quality == original_quality
            assert repeated[0]["sha256"] == original[0]["sha256"]
            assert sha256_file(tmp_path / "before.png") == sha256_file(tmp_path / "after.png")
            with Image.open(tmp_path / "after.png") as image:
                assert image.size == (4050, 2760)
                assert image.getpixel((0, 0)) == (255, 255, 255, 255)
        actual = original_quality["actual_fonts"]
        assert all(Path(item["path"]).is_file() and sha256_file(item["path"]) == item["sha256"] for item in actual)
        assert original_quality["font"]["resolved"] == "DejaVu Serif"
    finally:
        plt.close(figure)
        plt.close(other)


def test_failed_save_and_mismatched_profile_restore_external_style(tmp_path, monkeypatch):
    profile = theme()
    figure = create_figure(profile)
    before = dict(matplotlib.rcParams)
    try:
        with pytest.raises(ValueError, match="profile differs"):
            export_figure(figure, tmp_path / "bad", theme("DejaVu Sans"), formats=["png"])
        assert not (tmp_path / "bad.png").exists()
        def failed_save(*args, **kwargs):
            assert matplotlib.rcParams["font.family"] == ["DejaVu Serif"]
            raise OSError("failed save")
        monkeypatch.setattr(figure, "savefig", failed_save)
        with pytest.raises(OSError, match="failed save"):
            export_figure(figure, tmp_path / "failure", profile, formats=["png"])
        assert dict(matplotlib.rcParams) == before
    finally:
        plt.close(figure)


def test_caller_cannot_change_a_bound_profile_after_creation(tmp_path):
    profile = theme()
    figure = create_figure(profile)
    profile.data["font_size_pt"] = 6
    try:
        with pytest.raises(ValueError, match="profile differs"):
            export_figure(figure, tmp_path / "changed", profile, formats=["png"])
        assert measure_figure(figure, intended_width_inches=6.75)["font"]["resolved"] == "DejaVu Serif"
    finally:
        plt.close(figure)


def test_line_artist_landmark_sizes_are_measured_without_counting_error_caps():
    with style_context(theme()) as context:
        figure, ax = plt.subplots()
        context.bind(figure)
        try:
            ax.plot([0, 1], [0, 1], marker="s", markersize=2.5)
            ax.plot([0, 1], [1, 0], marker="_", markersize=1)
            quality = measure_figure(figure, intended_width_inches=figure.get_size_inches()[0])
            assert quality["minimum_landmark_diameter_mm"] == pytest.approx(2.5 * 25.4 / 72)
            assert quality["minimum_landmark_diameter_mm"] < 1
        finally:
            plt.close(figure)


def test_font_error_declared_substitution_and_no_implicit_default():
    base = theme()
    missing = replace(base, data={**base.data, "font_candidates": ["Missing_SciPlot_Font_7c0d"]})
    before = dict(matplotlib.rcParams)
    with pytest.raises(ValueError, match="required font not found"):
        configure_matplotlib(missing)
    assert dict(matplotlib.rcParams) == before
    substitute = replace(missing, data={**missing.data, "font_policy": "candidate_fallback",
        "font_candidates": ["Missing_SciPlot_Font_7c0d", "DejaVu Serif"]})
    with pytest.warns(UserWarning, match="declared font substitute"):
        record = resolve_font(substitute)
    assert record["fallback_used"] is True
    assert record["missing_candidates"] == ["Missing_SciPlot_Font_7c0d"]
    assert record["selected"] == "DejaVu Serif"
    assert record["path"] == str(Path(font_manager.findfont("DejaVu Serif", fallback_to_default=False)).resolve())
    with pytest.raises(ValueError, match="required font not found"):
        resolve_font(replace(missing, data={**missing.data, "font_policy": "candidate_fallback"}))


def test_legacy_cli_keeps_its_existing_global_defaults_and_font_record():
    profile = load_profile("audit_v1")
    with matplotlib.rc_context({"text.usetex": True, "savefig.transparent": True, "mathtext.default": "regular"}):
        record = configure_matplotlib(profile)
        assert set(record) == {"requested", "resolved", "path"}
        assert matplotlib.rcParams["font.family"] == ["serif"]
        assert matplotlib.rcParams["font.serif"] == profile.data["font_candidates"]
        assert matplotlib.rcParams["text.usetex"] is True
        assert matplotlib.rcParams["savefig.transparent"] is True
        assert matplotlib.rcParams["mathtext.default"] == "regular"


def test_old_frozen_v2_profile_can_be_read_but_new_render_requires_a_policy(tmp_path):
    current = load_profile("candidate_aps_v2")
    snapshot = tmp_path / "candidate_aps_v2.toml"
    text = current.path.read_text(encoding="utf-8")
    snapshot.write_text("\n".join(line for line in text.splitlines() if not line.startswith("font_policy =")), encoding="utf-8")
    historical = load_profile(snapshot)
    assert historical.profile_id == current.profile_id
    assert historical.formats == current.formats
    with pytest.raises(ValueError, match="explicit missing-font policy"):
        with style_context(historical):
            pass


def test_delayed_pdf_uses_embedded_non_type3_vector_fonts(tmp_path):
    if not all(shutil.which(tool) for tool in ("pdfinfo", "pdffonts", "pdfimages")):
        pytest.skip("Poppler is required for the PDF behavior check")
    profile = theme()
    figure = create_figure(profile)
    try:
        with matplotlib.rc_context({"pdf.fonttype": 3, "mathtext.fontset": "cm", "savefig.transparent": True}):
            before = dict(matplotlib.rcParams)
            outputs, _ = export_figure(figure, tmp_path / "vector", profile, formats=["pdf"])
            assert dict(matplotlib.rcParams) == before
        assert outputs[0]["inspection"]["fonts_embedded"] is True
        assert outputs[0]["inspection"]["type3_font_count"] == 0
        assert outputs[0]["inspection"]["raster_image_count"] == 0
    finally:
        plt.close(figure)
