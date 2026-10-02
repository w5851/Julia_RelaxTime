import json
import hashlib
from pathlib import Path
import shutil
import subprocess
import zipfile

import pytest

from scripts.plotting.plot_manifest import generator_record, input_record
from scripts.plotting.plot_provenance import code_ref_for_manifest, git_source, validate_hash_record, validate_snapshot


@pytest.fixture
def frozen(tmp_path):
    if shutil.which("git") is None:
        pytest.skip("Git unavailable")
    code = tmp_path / "scripts" / "draw.py"
    code.parent.mkdir()
    code.write_bytes(b"print('original renderer')\n")
    data = tmp_path / "input.csv"
    data.write_bytes(b"x,y\n1,2\n")
    output = tmp_path / "figure.png"
    output.write_bytes(b"fixture image bytes")
    manifest = tmp_path / "manifest.json"
    record = generator_record(code, command="fixture", root=tmp_path)
    manifest.write_text(json.dumps({"generator": record, "inputs": [input_record(data, role="data", root=tmp_path)],
                                    "outputs": [input_record(output, role="output", root=tmp_path)]}), encoding="utf-8")
    def git(*args):
        return subprocess.check_output(["git", "-C", str(tmp_path), *args], stderr=subprocess.PIPE).decode().strip()
    git("init")
    git("add", ".")
    git("-c", "user.name=Fixture", "-c", "user.email=fixture@example.invalid", "commit", "-m", "fixture")
    return tmp_path, code, data, output, manifest, record, git("rev-parse", "HEAD")


def test_historical_code_survives_current_edits_and_removal(frozen):
    root, code, _, _, manifest, record, commit = frozen
    code.write_bytes(b"print('new renderer')\n")
    assert validate_hash_record(record, root=root, label="generator")
    assert validate_hash_record(record, root=root, label="generator", code_ref=commit) == []
    assert validate_snapshot(manifest, root=root, code_ref=commit) == []
    code.unlink()
    assert validate_hash_record(record, root=root, label="generator", code_ref=commit) == []


@pytest.mark.parametrize("which", ["input", "output", "manifest"])
def test_snapshot_still_rejects_artifact_tampering(frozen, which):
    root, _, data, output, manifest, _, commit = frozen
    {"input": data, "output": output, "manifest": manifest}[which].write_bytes(b"modified")
    assert validate_snapshot(manifest, root=root, code_ref=commit)


def test_generator_snapshot_preserves_exact_bytes(frozen):
    root, code, _, _, _, _, _ = frozen
    snapshot = root / "original.py"
    snapshot.write_bytes(code.read_bytes())
    record = generator_record(code, command="fixture", root=root, source_snapshot=snapshot)
    code.unlink()
    assert validate_hash_record(record, root=root, label="generator") == []
    snapshot.write_bytes(b"tampered source")
    assert validate_hash_record(record, root=root, label="generator")


def test_contract_history_does_not_substitute_historical_data(frozen):
    root, code, data, _, _, record, commit = frozen
    code.write_bytes(b"changed")
    assert validate_hash_record(record, root=root, label="generator", code_ref="0" * 40)
    original = input_record(data, role="data", root=root)
    data.write_bytes(b"changed")
    assert validate_hash_record(original, root=root, label="inputs[0]", code_ref=commit)


def test_registered_legacy_case_selects_its_retained_code(frozen):
    root, _, _, _, _, _, commit = frozen
    registry = root / "config/plotting/historical_snapshots.toml"
    registry.parent.mkdir(parents=True)
    registry.write_text(f'[fixture]\ncode_commit = "{commit}"\nmanifest_prefix = "retained_case"\n', encoding="utf-8")
    assert code_ref_for_manifest(root / "retained_case/figure.plot_manifest.json", root) == commit
    assert code_ref_for_manifest(root / "new_case/figure.plot_manifest.json", root) is None


def test_retained_source_archive_works_without_git_history(frozen, monkeypatch):
    root, code, _, _, _, record, commit = frozen
    archive = root / "code.zip"
    with zipfile.ZipFile(archive, "w") as snapshot:
        snapshot.write(code, "scripts/draw.py")
    registry = root / "config/plotting/historical_snapshots.toml"
    registry.parent.mkdir(parents=True)
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    registry.write_text(f'[[code_snapshots]]\ncode_commit = "{commit}"\narchive = "code.zip"\nsha256 = "{digest}"\n', encoding="utf-8")
    code.unlink()
    git_source.cache_clear()
    def unavailable(*args, **kwargs):
        raise AssertionError("historical Git access should not be required")
    monkeypatch.setattr(subprocess, "check_output", unavailable)
    assert validate_hash_record(record, root=root, label="generator", code_ref=commit) == []
    archive.write_bytes(b"damaged archive")
    git_source.cache_clear()
    assert validate_hash_record(record, root=root, label="generator", code_ref=commit)
