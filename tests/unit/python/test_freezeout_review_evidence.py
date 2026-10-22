"""Archive restoration preserves byte identity and rejects unsafe replacement."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
import warnings
import zipfile

from scripts.analysis.relaxtime.restore_freezeout_review_evidence import restore, verify


class FreezeoutEvidenceTests(unittest.TestCase):
    def fixture(self, root, members):
        archive = root / "raw.zip"
        with zipfile.ZipFile(archive, "w") as handle:
            for name, data in members:
                handle.writestr(name, data)
        data = archive.read_bytes()
        manifest = root / "manifest.json"
        manifest.write_text(json.dumps(dict(schema="freezeout_review_evidence_v1", retained_files=[],
            archives=[dict(path="raw.zip", bytes=len(data), sha256=hashlib.sha256(data).hexdigest(),
                           files=len(members), uncompressed_bytes=sum(len(v) for _, v in members))])),
            encoding="utf-8")
        return manifest

    def test_restore_is_byte_exact_and_idempotent(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = self.fixture(root, [("points/one.json", b"frozen\r\n")])
            self.assertEqual(restore(manifest, root)["restored_files"], 1)
            self.assertEqual((root / "points/one.json").read_bytes(), b"frozen\r\n")
            self.assertEqual(restore(manifest, root)["restored_files"], 0)

    def test_conflict_is_rejected_before_any_member_is_written(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = self.fixture(root, [("new.json", b"new"), ("existing.json", b"old")])
            (root / "existing.json").write_bytes(b"user changes")
            with self.assertRaisesRegex(ValueError, "replace existing evidence"):
                restore(manifest, root)
            self.assertFalse((root / "new.json").exists())
            self.assertEqual((root / "existing.json").read_bytes(), b"user changes")

    def test_escaping_and_duplicate_paths_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = self.fixture(root, [("../outside.json", b"outside")])
            with self.assertRaisesRegex(ValueError, "repository-relative"):
                verify(manifest, root)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                manifest = self.fixture(root, [("one.json", b"one"), ("one.json", b"two")])
            with self.assertRaisesRegex(ValueError, "duplicate"):
                verify(manifest, root)
            manifest = self.fixture(root, [("one.json", b"one"), ("./one.json", b"two")])
            with self.assertRaisesRegex(ValueError, "duplicate"):
                verify(manifest, root)

    def test_modified_archive_bytes_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = self.fixture(root, [("one.json", b"one")])
            with (root / "raw.zip").open("ab") as handle:
                handle.write(b"changed")
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                restore(manifest, root)
            self.assertFalse((root / "one.json").exists())


if __name__ == "__main__":
    unittest.main()
