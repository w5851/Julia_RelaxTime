"""Exercise all Stage 1 renderers through the hash-bound two-stage CLI."""

def _times_new_roman_available() -> bool:
    """The rendering profile needs the manuscript serif; Linux CI lacks it."""
    from matplotlib import font_manager

    try:
        font_manager.findfont("Times New Roman", fallback_to_default=False)
        return True
    except Exception:
        return False


import pytest as _pytest  # noqa: E402  (module-level skip marker)

_pytestmark = _pytest.mark.skipif(
    not _times_new_roman_available(),
    reason="module-level skip: Times New Roman is unavailable, so the rendering profile cannot run",
)


import json
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

from PIL import Image
import pytest

from scripts.plotting.plot_bundle import BUNDLE_SCHEMA, expand_bundle
from scripts.plotting.plot_manifest import input_record, sha256_file
from scripts.plotting.validate_plot_artifact import validate_manifest


ROOT = Path(__file__).resolve().parents[3]
MESON = "data/outputs/results/relaxtime/meson_density"
WARD = "data/outputs/results/relaxtime/transport/ward_moller/ward_moller_sphere_20261121_v1"
C2 = "scripts/pnjl/plot_phase_diagram.py"
POPPLER = ("pdfinfo", "pdffonts", "pdfimages")


def _run(script, args, *, check=True):
    process = subprocess.run([sys.executable, "-X", "utf8", script, *map(str, args)], cwd=ROOT,
                             capture_output=True, text=True, encoding="utf-8", errors="replace")
    if check:
        assert process.returncode == 0, f"{process.stdout}\n{process.stderr}"
    return process


def _records(path):
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") == BUNDLE_SCHEMA:
        assert manifest["shared"]["generator"]["path"].endswith("plot_ward_moller_manuscript.py")
        return [record for _, record in expand_bundle(manifest)]
    return [manifest]


def _synthetic_acceptance(directory, review_manifest, records):
    """Build ephemeral test fixtures; no approval or PDF ever enters the repository."""
    flags = []
    for i, record in enumerate(records):
        acceptance = directory / f"synthetic-acceptance-{i}.json"
        color = next(out for out in record["outputs"] if out["role"] == "color_review")
        payload = dict(schema="plot_png_acceptance_v1", status="author_accepted",
                       author_instruction="Synthetic pytest fixture, not an author response.",
                       pdf_export_authorized=True,
                       accepted_manifest=input_record(review_manifest, role="accepted_png_manifest"),
                       accepted_png=color)
        if len(records) > 1:
            payload["accepted_figure_id"] = record["asset_id"]
        acceptance.write_text(json.dumps(payload), encoding="utf-8")
        flags.extend(("--acceptance", acceptance))
    return flags


@pytest.mark.skipif(not all(shutil.which(name) for name in POPPLER), reason="requires Poppler PDF tools")
@pytest.mark.parametrize(
    "script,base_args,out_opt,results_opt,count",
    [
        pytest.param(C2, ["--review-c2"], "--output-dir", "--results-dir", 1, id="C2"),
        pytest.param("scripts/analysis/relaxtime/plot_ward_moller_manuscript.py",
                     ["--kind", "all", "--results", WARD], "--out", "--review-results-dir", 5, id="Ward"),
        *[
            pytest.param("scripts/analysis/relaxtime/plot_freezeout_thermodynamic_envelope.py",
                         ["render", "--prepared-dir", f"{MESON}/{family}_20261018_Tmu_prepared"],
                         "--output-dir", "--results-dir", 1, id=family)
            for family in ("freezeout_envelope", "andronic_envelope")
        ],
    ],
)
def test_stage1_acceptance_to_vector_delivery(tmp_path, script, base_args, out_opt, results_opt, count):
    with TemporaryDirectory(dir=tmp_path, prefix="stage1-") as scratch:
        root = Path(scratch)
        review_dir, vector_dir = root / "review", root / "vector"
        _run(script, [*base_args, out_opt, review_dir, results_opt, root / "review-results"])
        review_manifest = review_dir / "plot_manifest.json"
        assert b"\r" not in review_manifest.read_bytes()
        if script == C2:
            assert b"\r" not in (root / "review-results" / "selected_slices.json").read_bytes()
        assert validate_manifest(review_manifest, repo_root=ROOT) == []
        reviewed = _records(review_manifest)
        assert len(reviewed) == count

        flags = _synthetic_acceptance(root, review_manifest, reviewed)
        if script == C2:
            permit_file = flags[1]
            permit = json.loads(permit_file.read_text(encoding="utf-8"))
            unsigned = root / "not-authorized.json"
            unsigned.write_text(json.dumps(dict(permit, pdf_export_authorized=False)), encoding="utf-8")
            rejected = _run(script, [*base_args, out_opt, root / "rejected", results_opt, root / "rejected-results",
                                     "--deliver-vector", "--acceptance", unsigned], check=False)
            assert rejected.returncode != 0 and "PDF authorization" in rejected.stderr
            assert not (root / "rejected").exists()

        _run(script, [*base_args, out_opt, vector_dir, results_opt, root / "vector-results",
                      "--deliver-vector", *flags])
        vector_manifest = vector_dir / "plot_manifest.json"
        assert b"\r" not in vector_manifest.read_bytes()
        assert validate_manifest(vector_manifest, repo_root=ROOT) == []
        delivered = _records(vector_manifest)
        assert len(delivered) == count
        for review, vector in zip(reviewed, delivered, strict=True):
            assert review["asset_id"] == vector["asset_id"]
            assert review["generator"]["sha256"] == vector["generator"]["sha256"] == sha256_file(ROOT / script)
            for record, stage, pending in ((review, "png_review", True), (vector, "vector_delivery", False)):
                assert record["delivery_stage"] == record["rendering"]["delivery_stage"] == stage
                assert record["vector_delivery_pending"] is record["rendering"]["vector_delivery_pending"] is pending
                assert record["rendering"]["color_route"] == "color_print_and_online"
                assert record["manuscript_eligible"] is record["current_publication_layer"] is False
            assert vector["rendering"]["quality"]["actual_fonts"] == review["rendering"]["quality"]["actual_fonts"]
            inputs = lambda record: sorted((item["role"], item["bytes"], item["sha256"]) for item in record["inputs"])
            assert inputs(vector) == inputs(review)
            for role in ("color_review", "grayscale_review"):
                source, exported = (next(item for item in record["outputs"] if item.get("role") == role)
                                    for record in (review, vector))
                assert source["sha256"] == exported["sha256"]
                assert (ROOT / source["path"]).read_bytes() == (ROOT / exported["path"]).read_bytes()
                if role == "grayscale_review":
                    with Image.open(ROOT / exported["path"]) as gray:
                        assert gray.mode == "L"
                        assert gray.info["dpi"] == pytest.approx((600, 600), abs=.02)
            pdfs = [out for out in vector["outputs"] if out["format"] == "pdf"]
            assert len(pdfs) == 1 and pdfs[0]["vector"]
            pdf = pdfs[0]["inspection"]
            assert pdf["page_count"] == 1 and pdf["fonts_embedded"]
            assert pdf["type3_font_count"] == pdf["raster_image_count"] == 0
            assert {out["format"] for out in vector["outputs"]} == {"png", "pdf"}

        if count == 5:
            mechanism = next(item for item in reviewed if item["asset_id"] == "ward_moller_mott_transport_mechanism")
            windows = {p["panel"]: p["window_display"] for p in mechanism["rendering"]["panel_specs"]
                       if "window_display" in p}
            # The current frozen source contains the strange-channel correction in #365.
            assert windows == {
                "g": {"delta_tau_percent": "+1.62%", "delta_J_percent": "\u22122.38%", "sigma_change_percent": "\u22120.76%"},
                "h": {"delta_tau_percent": "+1.40%", "delta_J_percent": "\u22121.35%", "sigma_change_percent": "+0.05%"},
            }
