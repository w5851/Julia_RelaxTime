"""A production baseline must retain its independent qualification targets."""
import importlib.util
import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
import zipfile

ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location(
    "transport_baseline", ROOT / "scripts/dev/export_thermal_weighted_transport_baseline.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class ProductionBaselineTests(unittest.TestCase):
    def setUp(self):
        self.target = dict(point_id="p", T_MeV="200", muB_MeV="0", xi="0.37",
                           observable="tau_u", candidate="2", reference="2.001")
        self.actual = dict(point_id="p", observable="tau_u", qualified_candidate="2",
                           qualified_reference="2.001", production="2.000003", passed="true")

    def test_qualified_production_value_keeps_both_targets(self):
        row = MODULE.production_values([self.target], [self.actual])[0]
        self.assertEqual(row["candidate"], "2.000003")
        self.assertEqual(row["qualified_candidate"], "2")
        self.assertEqual(row["reference"], "2.001")
        self.assertEqual(row["T_MeV"], self.target["T_MeV"])

    def test_missing_duplicate_unknown_or_changed_targets_fail(self):
        for rows in ([], [self.actual, self.actual],
                     [{**self.actual, "point_id": "other"}],
                     [{**self.actual, "qualified_candidate": "2.000003"}],
                     [{**self.actual, "qualified_reference": "2"}]):
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                MODULE.production_values([self.target], rows)

    def test_failed_nonfinite_or_unqualified_values_fail(self):
        for row in ({**self.actual, "passed": "false"},
                    *({**self.actual, "production": value} for value in
                      ("nan", "inf", "0", "-1", "2.001"))):
            with self.subTest(row=row), self.assertRaises(ValueError):
                MODULE.production_values([self.target], [row])

    def test_independent_limit_is_separate_from_candidate_replay(self):
        target = {**self.target, "reference": "1.9"}
        actual = {**self.actual, "qualified_reference": "1.9"}
        with self.assertRaises(ValueError):
            MODULE.production_values([target], [actual])


class ArchivedProductionBaselineTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.archive = self.root / "source.zip"
        self.contents = {"src/kernel.jl": b"kernel(x) = x * x\n",
                         "scripts/relaxtime/run_scan.jl": b"include(\"kernel.jl\")\n"}
        with zipfile.ZipFile(self.archive, "w") as frozen:
            for relative, data in self.contents.items():
                frozen.writestr(relative, data)
                path = self.root / "source" / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(data)
        files = {name: hashlib.sha256(data).hexdigest() for name, data in self.contents.items()}
        self.package = dict(head_base="a" * 40, source_branch="test", calculation_commit=None,
                            uncommitted_source=True, files=len(files), source_files=files,
                            zip_sha256=MODULE.digest(self.archive))
        self.manifest = dict(head=None, source_archive=dict(
            path=str(self.archive), sha256=self.package["zip_sha256"],
            head_base="a" * 40, source_branch="test", calculation_commit=None,
            uncommitted_source=True, verified_files=len(files),
            package_manifest_snapshot="source_package_manifest.json"),
            source_files=[dict(path=name, sha256=value) for name, value in files.items()])
        self.save_package()

    def save_package(self):
        snapshot = self.root / "source_package_manifest.json"
        snapshot.write_text(json.dumps(self.package), encoding="utf-8")
        self.manifest["source_archive"]["package_manifest_sha256"] = MODULE.digest(snapshot)

    def test_verified_archive_keeps_commit_unknown_and_base_separate(self):
        commit, identity = MODULE.replay_artifact_commit(self.root, self.manifest)
        self.assertEqual(commit, "unknown")
        self.assertIsNone(identity["calculation_commit"])
        self.assertEqual(identity["head_base"], "a" * 40)
        self.assertEqual(identity["sha256"], MODULE.digest(self.archive))
        self.assertEqual(identity["source_files_verified"], 2)

    def test_git_replay_retains_its_recorded_head(self):
        self.assertEqual(MODULE.replay_artifact_commit(self.root, {"head": "b" * 40}),
                         ("b" * 40, None))

    def test_unknown_without_archive_is_rejected(self):
        for head in (None, "unknown", "abc"):
            with self.subTest(head=head), self.assertRaises(ValueError):
                MODULE.replay_artifact_commit(self.root, {"head": head})

    def test_relocated_archive_must_keep_the_recorded_bytes(self):
        self.manifest["source_archive"]["path"] = str(self.root / "unavailable.zip")
        self.assertEqual(MODULE.replay_artifact_commit(self.root, self.manifest, self.archive)[0], "unknown")
        changed = self.root / "changed.zip"
        changed.write_bytes(self.archive.read_bytes() + b"changed")
        with self.assertRaises(ValueError):
            MODULE.replay_artifact_commit(self.root, self.manifest, changed)

    def test_worktree_base_cannot_be_used_as_calculation_commit(self):
        for key, value in (("head", "a" * 40), ("calculation_commit", "a" * 40),
                           ("uncommitted_source", False), ("head_base", "b" * 40)):
            manifest = copy.deepcopy(self.manifest)
            (manifest if key == "head" else manifest["source_archive"])[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                MODULE.replay_artifact_commit(self.root, manifest)

    def test_changed_package_snapshot_is_rejected(self):
        snapshot = self.root / "source_package_manifest.json"
        snapshot.write_bytes(snapshot.read_bytes() + b"\n")
        with self.assertRaises(ValueError):
            MODULE.replay_artifact_commit(self.root, self.manifest)

    def test_changed_source_snapshot_is_rejected(self):
        (self.root / "source/src/kernel.jl").write_bytes(b"changed source\n")
        with self.assertRaises(ValueError):
            MODULE.replay_artifact_commit(self.root, self.manifest)

    def test_archive_contents_must_match_the_file_hashes(self):
        with zipfile.ZipFile(self.archive, "w") as frozen:
            for relative, data in self.contents.items():
                frozen.writestr(relative, data + b"changed")
        self.package["zip_sha256"] = MODULE.digest(self.archive)
        self.manifest["source_archive"]["sha256"] = self.package["zip_sha256"]
        self.save_package()
        with self.assertRaises(ValueError):
            MODULE.replay_artifact_commit(self.root, self.manifest)

    def test_missing_or_duplicate_replay_snapshots_are_rejected(self):
        for entries in ([], self.manifest["source_files"] * 2):
            manifest = {**self.manifest, "source_files": entries}
            with self.subTest(entries=entries), self.assertRaises(ValueError):
                MODULE.replay_artifact_commit(self.root, manifest)


if __name__ == "__main__":
    unittest.main()
