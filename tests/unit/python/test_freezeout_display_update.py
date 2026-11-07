"""Keep displayed observations separate from frozen parameter selection."""
import copy
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from scripts.analysis.relaxtime.plot_freezeout_envelope import display_experiments
from scripts.analysis.relaxtime.plot_freezeout_comparison import load_experiment
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
        added = [r for r in displayed if r["id"] not in old_rows]
        self.assertEqual(len(added), 1)
        self.assertEqual(added[0]["energy"], 54.4)
        self.assertAlmostEqual(added[0]["ratio"], .15879)
        self.assertFalse(added[0]["primary_auau"])
        self.assertFalse(added[0]["auxiliary_horn"])

    def test_override_does_not_hide_changed_scoring_source(self):
        self.summary["experiment"]["sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "scoring input changed"):
            display_experiments(self.summary, self.new)

    def test_retired_worktree_reuses_identical_retained_table(self):
        self.summary["experiment"]["path"] = (
            "D:/retired-worktree/" + self.old.relative_to(ROOT).as_posix())
        _, scoring_path, _ = display_experiments(self.summary, self.new)
        self.assertEqual(scoring_path, self.old)


if __name__ == "__main__":
    unittest.main()
