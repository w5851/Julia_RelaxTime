"""Evidence round trips and rejection of overwrite, traversal and damaged archives."""
import hashlib
import json
from pathlib import Path
import zipfile

import pytest

from scripts.analysis.relaxtime import transport_delivery_archive as archive


def example(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    (source / "table.csv").write_bytes(b"x,y\r\n1,2\r\n")
    (source / "original.bin").write_bytes(bytes(range(256)))
    package, receipt_file = tmp_path / "evidence.zip", tmp_path / "receipt.json"
    receipt = archive.pack(source, ["table.csv", "original.bin"], package, receipt_file)
    return source, package, receipt


def test_roundtrip_preserves_raw_bytes_and_requires_new_destination(tmp_path):
    source, package, receipt = example(tmp_path)
    target = tmp_path / "restored"
    result = archive.restore(package, receipt, target)
    assert result["files"] == 2
    for name in ("table.csv", "original.bin"):
        assert (source / name).read_bytes() == (target / name).read_bytes()
    with pytest.raises(FileExistsError):
        archive.restore(package, receipt, target)


@pytest.mark.parametrize("name", ["../outside", "C:/outside", "/outside", "a/../b", "a\\b", "a/./b", "a /b"])
def test_invalid_names_are_rejected_before_writing(tmp_path, name):
    with pytest.raises(ValueError):
        archive.pack(tmp_path, [name], tmp_path / "evidence.zip", tmp_path / "receipt.json")
    assert not (tmp_path / "evidence.zip").exists()


def test_archive_tampering_does_not_create_restore_destination(tmp_path):
    _, package, receipt = example(tmp_path)
    with package.open("ab") as stream:
        stream.write(b"changed")
    with pytest.raises(ValueError, match="identity mismatch"):
        archive.restore(package, receipt, tmp_path / "restored")
    assert not (tmp_path / "restored").exists()


def test_manifest_cannot_hide_corrupted_payload(tmp_path):
    package = tmp_path / "bad.zip"
    raw = json.dumps(dict(files=[dict(path="table.csv", bytes=3, sha256=hashlib.sha256(b"old").hexdigest())])).encode()
    with zipfile.ZipFile(package, "w") as bundle:
        bundle.writestr(archive.MANIFEST, raw)
        bundle.writestr("table.csv", b"new")
    receipt = dict(archive_bytes=package.stat().st_size, archive_sha256=archive.digest(package),
                   manifest_sha256=hashlib.sha256(raw).hexdigest())
    with pytest.raises(ValueError, match="member identity mismatch"):
        archive.restore(package, receipt, tmp_path / "restored")
    assert not (tmp_path / "restored").exists()


def test_case_collision_and_overwrite_are_rejected(tmp_path):
    with pytest.raises(ValueError, match="case-colliding"):
        archive.pack(tmp_path, ["a", "A"], tmp_path / "evidence.zip", tmp_path / "receipt.json")
    source, package, _ = example(tmp_path)
    before = package.read_bytes()
    with pytest.raises(FileExistsError):
        archive.pack(source, ["table.csv"], package, tmp_path / "another_receipt.json")
    assert package.read_bytes() == before
