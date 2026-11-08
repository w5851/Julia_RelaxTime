"""Keep displayed observations separate from frozen parameter selection."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from scripts.analysis.relaxtime.plot_freezeout_envelope import display_experiments
from scripts.analysis.relaxtime.plot_freezeout_comparison import (
    load_experiment, experiment_groups, prepare_review_directories,
)
from scripts.plotting.plot_manifest import sha256_file


class DisplayUpdateTests(unittest.TestCase):
    def setUp(self):
        self.old = ROOT / "data/inputs/meson_freezeout_experiment_20261018/experimental_ratios.csv"
        self.new = ROOT / "data/inputs/meson_freezeout_experiment_20261107/experimental_ratios.csv"
        self.summary = {"experiment": {"path": str(self.old), "sha256": sha256_file(self.old)}}

    def test_new_observation_only_changes_display(self):
        original = copy.deepcopy(self.summary)
        path, scoring_path, displayed = display_experiments(self.summary, self.new)
        old_rows = {r["id"]: r for r in load_experiment(self.old)}
        all_new = {r["id"]: r for r in load_experiment(self.new)}
        self.assertEqual(self.summary, original)
        self.assertEqual(scoring_path, self.old)
        self.assertEqual(path, self.new)
        self.assertEqual({k: all_new[k] for k in old_rows}, old_rows)
        added = {r["id"]: r for r in displayed if r["id"] not in old_rows}
        self.assertTrue({"STAR_54.4", "NA61_PbPb_7.62_4pi", "BRAHMS_62.4_y0",
            "STAR_BESII_14.6", "STAR_BESII_19.6", "STAR_BESII_27"}.issubset(added))
        self.assertEqual(added["STAR_54.4"]["energy"], 54.4)
        self.assertAlmostEqual(added["STAR_54.4"]["ratio"], .15879)
        na61 = added["NA61_PbPb_7.62_4pi"]
        self.assertEqual(na61["energy"], 7.62)
        self.assertAlmostEqual(na61["ratio"], .23625)
        self.assertIn("full phase space", na61["rapidity"])
        self.assertEqual(na61["method"], "vector_figure_extraction")
        for row in added.values():
            self.assertFalse(row["primary_auau"])
            self.assertFalse(row["auxiliary_horn"])

    def test_full_phase_space_point_is_not_labeled_as_na49_midrapidity(self):
        _, _, rows = display_experiments(self.summary, self.new)
        groups = {g["name"]: g for g in experiment_groups(rows, na49_filled=False)}
        self.assertEqual([r["id"] for r in groups["NA61"]["rows"]], ["NA61_PbPb_7.62_4pi"])
        self.assertTrue(all(r["experiment"] == "NA49" for r in groups["NA49"]["rows"]))
        self.assertNotEqual(groups["NA61"]["marker"], groups["NA49"]["marker"])
        self.assertTrue(groups["NA49"]["filled"])
        self.assertFalse(groups["NA494pi"]["filled"])
        self.assertTrue(all("full phase space" not in r["rapidity"] for r in groups["NA49"]["rows"]))
        self.assertTrue(all("full phase space" in r["rapidity"] for r in groups["NA494pi"]["rows"]))
        self.assertEqual([r["id"] for r in groups["BRAHMS"]["rows"]], ["BRAHMS_62.4_y0"])
        self.assertTrue(all(r["experiment"] == "STAR BES-II" for r in groups["STARBESII"]["rows"]))
        self.assertEqual(sum(len(g["rows"]) for g in groups.values()), len(rows))

    def test_full_phase_space_scan_and_derived_ratios_keep_source_definitions(self):
        rows = {r["id"]: r for r in load_experiment(self.new)}
        raw = self.new.parent / "raw"
        source = json.loads((raw / "na49_4pi_energy_hepdata.json").read_text(encoding="utf-8"))
        _, _, displayed = display_experiments(self.summary, self.new)
        displayed_ids = {r["id"] for r in displayed}
        for item in source["values"]:
            energy = float(item["x"][0]["value"])
            experiment = "AGS" if energy < 6 else "NA49" if energy < 100 else "BRAHMS"
            key = f"{experiment}_{energy:g}_4pi"
            self.assertIn(key, displayed_ids)
            self.assertEqual(rows[key]["ratio"], float(item["y"][0]["value"]))
        derived = json.loads((raw / "star_besii2024_midrapidity.json").read_text(encoding="utf-8"))
        for point in derived["points"]:
            row = rows[f"STAR_BESII_{point['energy_GeV']:g}"]
            self.assertAlmostEqual(row["ratio"], point["K_dNdy"] / point["pi_dNdy"], places=5)
            self.assertIn("zero_covariance", row["error_method"])
            self.assertEqual(row["method"], "ratio_of_digitized_yields")

    def test_override_does_not_hide_changed_scoring_source(self):
        self.summary["experiment"]["sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "scoring input changed"):
            display_experiments(self.summary, self.new)

    def test_explicit_draft_replacement_preserves_other_files_and_rejects_acceptance(self):
        with tempfile.TemporaryDirectory() as temporary:
            output, results = Path(temporary) / "figures", Path(temporary) / "results"
            prepare_review_directories(output, results)
            extra = results / "author-note.md"
            extra.write_text("keep", encoding="utf-8")
            manifest = output / "plot_manifest.json"
            manifest.write_text(json.dumps({"delivery_stage": "png_review",
                "validation": {"author_review": "pending"}}), encoding="utf-8")
            with self.assertRaises(FileExistsError):
                prepare_review_directories(output, results)
            prepare_review_directories(output, results, replace_draft=True)
            self.assertEqual(extra.read_text(encoding="utf-8"), "keep")
            manifest.write_text(json.dumps({"delivery_stage": "png_review",
                "validation": {"author_review": "accepted"}}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unaccepted PNG"):
                prepare_review_directories(output, results, replace_draft=True)

    def test_retired_worktree_reuses_identical_retained_table(self):
        self.summary["experiment"]["path"] = (
            "D:/retired-worktree/" + self.old.relative_to(ROOT).as_posix())
        _, scoring_path, _ = display_experiments(self.summary, self.new)
        self.assertEqual(scoring_path, self.old)


if __name__ == "__main__":
    unittest.main()
