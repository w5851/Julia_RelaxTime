from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
import unittest

from scripts.plotting.plot_provenance import validate_hash_record

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts/analysis/relaxtime/inspect_phase_guided_publication_clean_v11_placement.py"
SPEC = importlib.util.spec_from_file_location("v11_placement_tests", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
try:
    SPEC.loader.exec_module(MODULE)
    DOCUMENT_DEPENDENCIES_AVAILABLE = True
except ModuleNotFoundError:
    DOCUMENT_DEPENDENCIES_AVAILABLE = False


@unittest.skipUnless(DOCUMENT_DEPENDENCIES_AVAILABLE, "run with bundled document Python")
class PlacementTests(unittest.TestCase):
    def test_a4_fit_uses_page_box_and_declared_margins(self):
        self.assertAlmostEqual(MODULE.a4_scale(612, 792, 0), 210 / 215.9)
        self.assertAlmostEqual(MODULE.a4_scale(612, 792, 5), 200 / 215.9)

    def test_recorded_placement_report_is_portable(self):
        path = MODULE.DEFAULT_REPORT_ROOT / "placement_report.json"
        report = json.loads(path.read_text(encoding="utf-8"))
        self.assertFalse(report["paper_project_modified"])
        self.assertFalse(report["manuscript_eligible"])
        self.assertEqual(len(report["placement"]), 2)
        for record in report["placement"]:
            self.assertAlmostEqual(record["placed_width_mm"], 179.25, places=2)
            self.assertGreater(record["minimum_primary_glyph_height_mm"], 2)
            self.assertLess(record["minimum_math_script_height_mm"], 2)
            self.assertGreater(record["caption_gap_mm"], 0)
            self.assertLess(record["image_bbox_pdf_pt"][3], record["page_size_pdf_pt"][1])
        for record in report["placement"]:
            path = ROOT / "data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v11_png_review/composites" / (record["figure"] + ".plot_manifest.json")
            historical = {**record["native_png_manifest"], "path": path.relative_to(ROOT).as_posix()}
            self.assertEqual(validate_hash_record(historical, root=ROOT, label="native_png_manifest"), [])
            source = MODULE.figure_manifest(record["figure"])["rendering"]["quality"]
            self.assertAlmostEqual(source["minimum_primary_capital_numeral_height_mm"] * record["placement_scale"], record["minimum_primary_glyph_height_mm"])
            for scenario in record["print_scenarios"]:
                self.assertAlmostEqual(scenario["figure_width_mm"], record["placed_width_mm"] * scenario["print_scale"])

    def test_current_composites_resolve_from_the_bundle(self):
        for figure in MODULE.FIGURES:
            manifest = MODULE.figure_manifest(figure)
            self.assertTrue(manifest["outputs"][0]["path"].endswith(f"/composites/{figure}.png"))
            self.assertFalse((MODULE.PNG_ROOT / f"{figure}.plot_manifest.json").exists())

    def test_missing_composite_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "exactly one composite record"):
            MODULE.figure_manifest("missing_composite")

    @unittest.skipUnless(os.environ.get("JRT_LOCAL_PAPER_PREVIEW_CHECK") == "1", "external manuscript recheck requires explicit local opt-in")
    def test_local_paper_inputs_and_vector_print_previews(self):
        report = json.loads((MODULE.DEFAULT_REPORT_ROOT / "placement_report.json").read_text(encoding="utf-8"))
        measured = MODULE.measure_placement(Path(report["inputs"][0]["path"]))
        self.assertEqual(len(measured), len(report["placement"]))
        for current, original in zip(measured, report["placement"]):
            self.assertEqual({key: value for key, value in current.items() if key != "native_png_manifest"},
                             {key: value for key, value in original.items() if key != "native_png_manifest"})
            self.assertEqual(current["native_png_manifest"]["path"], str(MODULE.PNG_MANIFEST.resolve()))
        for record in [*report["inputs"], *report["outputs"]]:
            if Path(record["path"]).resolve() == SCRIPT.resolve():
                self.assertEqual(validate_hash_record(record, root=ROOT, label="generator",
                                 code_ref="0725c0c4f87bbcbced4b7641c0a086e65c6fe42b"), [])
            else:
                self.assertEqual(MODULE.file_record(Path(record["path"])), record)
        for record in report["outputs"]:
            reader = MODULE.PdfReader(record["path"])
            self.assertEqual(len(reader.pages), 2)
            for page in reader.pages:
                # PDF decimal serialization rounds page-box coordinates; this
                # checks physical dimensions to 0.00001 mm, not bit identity.
                self.assertAlmostEqual(float(page.mediabox.width) * MODULE.MM_PER_PDF_POINT, 210, places=5)
                self.assertAlmostEqual(float(page.mediabox.height) * MODULE.MM_PER_PDF_POINT, 297, places=5)
                resources = page["/Resources"]
                images = resources.get("/XObject", {})
                self.assertFalse(any(item.get_object().get("/Subtype") == "/Image" for item in images.values()))


if __name__ == "__main__":
    unittest.main()
