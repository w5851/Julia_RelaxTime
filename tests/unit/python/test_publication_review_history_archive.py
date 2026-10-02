from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path
import zipfile

import pytest

SCRIPT = Path(__file__).resolve().parents[3] / "scripts/analysis/relaxtime/archive_phase_guided_publication_review_history.py"
SPEC = importlib.util.spec_from_file_location("review_archive_tests", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_review_archive_verifies_exact_bytes_and_members(tmp_path):
    archive = tmp_path / "review.zip"
    payload = b"a,b\r\n1,2\r\n"
    with zipfile.ZipFile(archive, "w") as bundle:
        bundle.writestr("tables/points.csv", payload)
    records = [{"path": "tables/points.csv", "bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()}]
    MODULE.verify_archive(archive, records)
    records[0]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="byte/hash mismatch"):
        MODULE.verify_archive(archive, records)
    with pytest.raises(ValueError, match="exact index"):
        MODULE.verify_archive(archive, [])


@pytest.mark.parametrize("name", ["../x", "/absolute", "D:/outside", "a\\b"])
def test_review_archive_rejects_unsafe_member_names(name):
    assert not MODULE.safe_member(name)


def test_review_archive_leaves_v10_v11_outside_historical_allowlist():
    assert all("v10" not in name and "v11" not in name for name in MODULE.HISTORY_DIRECTORIES)
    assert len(MODULE.HISTORY_DIRECTORIES) == 11
