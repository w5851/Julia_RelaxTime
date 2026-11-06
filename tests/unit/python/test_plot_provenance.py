import json
import hashlib
from pathlib import Path
import shutil
import subprocess
import zipfile

import pytest

from scripts.plotting.plot_manifest import generator_record, input_record
from scripts.plotting.plot_provenance import code_ref_for_manifest, git_source, read_manifest_record, resolve_artifact_path, validate_hash_record, validate_snapshot


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


def test_moved_artifact_keeps_frozen_hash_and_rejects_stale_or_changed_bytes(tmp_path):
    old = tmp_path / "data/outputs/figures/old/plot_manifest.json"
    old.parent.mkdir(parents=True)
    old.write_bytes(b'{"status":"diagnostic"}\n')
    record = input_record(old, role="manifest", root=tmp_path)
    new = tmp_path / "data/outputs/figures/diagnostics/example/plot_manifest.json"
    new.parent.mkdir(parents=True)
    old.rename(new)
    registry = tmp_path / "config/plotting/artifact_locations.toml"
    registry.parent.mkdir(parents=True)
    registry.write_text('[[relocations]]\nfrom="data/outputs/figures/old"\n'
                        'to="data/outputs/figures/diagnostics/example"\n', encoding="utf-8")
    assert resolve_artifact_path(record["path"], root=tmp_path) == new
    assert validate_hash_record(record, root=tmp_path, label="outputs", allow_historical=False) == []
    assert read_manifest_record(record, root=tmp_path) == {"status": "diagnostic"}
    old.write_bytes(b"a stale duplicate must not be silently skipped")
    assert validate_hash_record(record, root=tmp_path, label="outputs", allow_historical=False)
    old.unlink()
    new.write_bytes(b"changed moved artifact")
    assert validate_hash_record(record, root=tmp_path, label="outputs", allow_historical=False)
    assert resolve_artifact_path("data/outputs/figures/old_extra/item.png", root=tmp_path) == (
        tmp_path / "data/outputs/figures/old_extra/item.png")


def test_relocation_resolves_source_snapshot_and_original_git_manifest(frozen):
    root, code, _, _, _, _, _ = frozen
    original = root / "data/outputs/results/old"
    original.mkdir(parents=True)
    snapshot = original / "renderer.py"
    snapshot.write_bytes(code.read_bytes())
    record = generator_record(code, command="fixture", root=root, source_snapshot=snapshot)
    manifest = original / "plot_manifest.json"
    manifest.write_text(json.dumps({"generator": record}), encoding="utf-8")
    subprocess.run(["git", "-C", str(root), "add", "."], check=True, capture_output=True)
    subprocess.run(["git", "-C", str(root), "-c", "user.name=Fixture", "-c",
                    "user.email=fixture@example.invalid", "commit", "-m", "snapshot"], check=True, capture_output=True)
    commit = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"]).decode().strip()
    moved = original.with_name("diagnostics")
    original.rename(moved)
    registry = root / "config/plotting/artifact_locations.toml"
    registry.parent.mkdir(parents=True)
    registry.write_text('[[relocations]]\nfrom="data/outputs/results/old"\n'
                        'to="data/outputs/results/diagnostics"\n', encoding="utf-8")
    code.write_bytes(b"new renderer")
    assert validate_hash_record(record, root=root, label="generator") == []
    assert validate_snapshot(moved / manifest.name, root=root, code_ref=commit) == []
    subprocess.run(["git", "-C", str(root), "add", "."], check=True, capture_output=True)
    subprocess.run(["git", "-C", str(root), "-c", "user.name=Fixture", "-c",
                    "user.email=fixture@example.invalid", "commit", "-m", "relocated"], check=True, capture_output=True)
    moved_commit = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"]).decode().strip()
    assert validate_snapshot(moved / manifest.name, root=root, code_ref=moved_commit) == []
    (moved / snapshot.name).write_bytes(b"damaged source")
    assert validate_hash_record(record, root=root, label="generator")


def test_relocation_cannot_escape_output_tree(tmp_path):
    registry = tmp_path / "config/plotting/artifact_locations.toml"
    registry.parent.mkdir(parents=True)
    registry.write_text('[[relocations]]\nfrom="data/outputs/old"\nto="../outside"\n', encoding="utf-8")
    with pytest.raises(ValueError, match="data/outputs"):
        resolve_artifact_path("data/outputs/old/figure.png", root=tmp_path)


def test_contract_history_does_not_substitute_historical_data(frozen):
    root, code, data, _, _, record, commit = frozen
    code.write_bytes(b"changed")
    assert validate_hash_record(record, root=root, label="generator", code_ref="0" * 40)
    original = input_record(data, role="data", root=root)
    data.write_bytes(b"changed")
    assert validate_hash_record(original, root=root, label="inputs[0]", code_ref=commit)


@pytest.mark.parametrize("relative,is_selector", [
    ("docs/analysis/relaxtime/phase_guided_transport/publication_clean_current.json", True),
    ("docs/analysis/relaxtime/phase_guided_transport/numerical_result.json", False),
    ("data/outputs/results/publication_clean_current.json", False),
])
def test_archived_publication_selector_requires_exact_path_and_hash(tmp_path, relative, is_selector):
    path = tmp_path / relative
    path.parent.mkdir(parents=True)
    path.write_bytes(b'{"selected":"old"}\n')
    record = input_record(path, role="historical_input", root=tmp_path)
    archive = tmp_path / "metadata.zip"
    with zipfile.ZipFile(archive, "w") as snapshot:
        snapshot.write(path, relative)
    registry = tmp_path / "config/plotting/historical_snapshots.toml"
    registry.parent.mkdir(parents=True)
    commit = "a" * 40
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    registry.write_text(f'[[code_snapshots]]\ncode_commit="{commit}"\narchive="metadata.zip"\nsha256="{digest}"\n', encoding="utf-8")
    path.write_bytes(b'{"selected":"new"}\n')
    assert bool(validate_hash_record(record, root=tmp_path, label="inputs", code_ref=commit)) is not is_selector
    assert validate_hash_record(record, root=tmp_path, label="outputs", code_ref=commit, allow_historical=False)
    assert validate_hash_record({**record, "sha256": "0" * 64}, root=tmp_path, label="inputs", code_ref=commit)
    archive.write_bytes(b"damaged metadata archive")
    git_source.cache_clear()
    assert validate_hash_record(record, root=tmp_path, label="inputs", code_ref=commit)


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


def test_retired_manifest_archive_preserves_only_explicit_manifest_references(tmp_path):
    manifest = tmp_path / "data/outputs/figures/case/figure.plot_manifest.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text('{"status":"frozen"}\n', encoding="utf-8")
    record = input_record(manifest, role="frozen_manifest", root=tmp_path)
    data = tmp_path / "data/outputs/results/raw.json"
    data.parent.mkdir(parents=True)
    data.write_text('{"value":1}', encoding="utf-8")
    data_record = input_record(data, role="numerical_data", root=tmp_path)
    archive = tmp_path / "manifests.zip"
    with zipfile.ZipFile(archive, "w") as snapshot:
        snapshot.write(manifest, record["path"])
        snapshot.write(data, data_record["path"])
    registry = tmp_path / "config/plotting/historical_snapshots.toml"
    registry.parent.mkdir(parents=True)
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    registry.write_text(f'[[manifest_archives]]\narchive="manifests.zip"\nsha256="{digest}"\n', encoding="utf-8")
    manifest.unlink()
    assert validate_hash_record(record, root=tmp_path, label="outputs", allow_historical=False) == []
    assert read_manifest_record(record, root=tmp_path) == {"status": "frozen"}
    assert validate_hash_record({**record, "sha256": "0" * 64}, root=tmp_path, label="outputs")
    data.write_text('{"value":2}', encoding="utf-8")
    assert validate_hash_record(data_record, root=tmp_path, label="inputs")
    archive.write_bytes(b"tampered archive")
    assert validate_hash_record(record, root=tmp_path, label="outputs")


def test_uncommitted_source_archive_matches_exact_hash_and_preserves_artifact_checks(tmp_path):
    originals = {
        "scripts/draw.py": b"# uncommitted v12 renderer\n",
        "docs/analysis/example/plotting_case_contract.md": b"v12 geometry contract\n",
        "docs/analysis/example/input.csv": b"x,y\n1,2\n",
        "data/outputs/figures/example/plot_manifest.json": b'{"status":"review"}',
    }
    records = []
    archive = tmp_path / "source.zip"
    with zipfile.ZipFile(archive, "w") as snapshot:
        for relative, payload in originals.items():
            path = tmp_path / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(payload)
            records.append(input_record(path, role="fixture", root=tmp_path))
            snapshot.writestr(relative, payload)
            path.write_bytes(b"changed")
    registry = tmp_path / "config/plotting/historical_snapshots.toml"
    registry.parent.mkdir(parents=True)
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    registry.write_text(f'[[working_tree_source_archives]]\narchive="source.zip"\nsha256="{digest}"\n', encoding="utf-8")
    for record in records[:2]:
        assert validate_hash_record(record, root=tmp_path, label="inputs") == []
        assert validate_hash_record(record, root=tmp_path, label="outputs", allow_historical=False)
    for record in records[2:]:
        assert validate_hash_record(record, root=tmp_path, label="inputs")
    assert validate_hash_record({**records[0], "sha256": "0" * 64}, root=tmp_path, label="inputs")
    archive.write_bytes(b"tampered source archive")
    assert validate_hash_record(records[0], root=tmp_path, label="inputs")


@pytest.mark.parametrize("relative,is_code", [("tests/unit/test_plot.py", True),
                                            ("tests/baselines/result.csv", False)])
def test_test_source_snapshot_does_not_relax_baseline_identity(tmp_path, relative, is_code):
    path = tmp_path / relative
    path.parent.mkdir(parents=True)
    path.write_bytes(b"original\n")
    record = input_record(path, role="test_source", root=tmp_path)
    snapshot = tmp_path / "frozen" / Path(relative).name
    snapshot.parent.mkdir()
    snapshot.write_bytes(path.read_bytes())
    record["source_snapshot"] = input_record(snapshot, role="immutable_source_copy", root=tmp_path)
    path.write_bytes(b"updated\n")
    assert bool(validate_hash_record(record, root=tmp_path, label="inputs")) is not is_code
    assert validate_hash_record(record, root=tmp_path, label="outputs", allow_historical=False)
    snapshot.write_bytes(b"damaged\n")
    assert validate_hash_record(record, root=tmp_path, label="inputs")
