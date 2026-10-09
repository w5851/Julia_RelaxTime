"""Verify artifact bytes and historical code without rewriting frozen records."""

from functools import lru_cache
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tomllib
import zipfile
from typing import Any


@lru_cache(maxsize=8)
def _verified_archive(path: str, expected: str, mtime_ns: int, size: int) -> dict[str, bytes]:
    payload = Path(path).read_bytes()
    if hashlib.sha256(payload).hexdigest() != expected:
        raise ValueError(f"historical archive hash mismatch: {path}")
    with zipfile.ZipFile(Path(path)) as snapshot:
        return {name: snapshot.read(name) for name in snapshot.namelist()}


def record_bytes(record: dict[str, Any], *, root: Path) -> bytes:
    """Resolve current bytes, or an explicitly archived pre-migration manifest."""
    value = record["path"]
    path = Path(value)
    path = path if path.is_absolute() else root / path
    expected = record.get("sha256")
    current = path.read_bytes() if path.is_file() else None
    if current is not None and (not expected or hashlib.sha256(current).hexdigest() == expected):
        return current
    registry = root / "config/plotting/historical_snapshots.toml"
    if expected and registry.is_file() and path.suffix == ".json":
        try:
            relative = path.resolve().relative_to(root.resolve()).as_posix()
        except ValueError:
            relative = ""
        manifest_name = (path.name in {"plot_manifest.json", "pdf_review_index.json"}
                         or path.name.endswith((".plot_manifest.json", ".pdf_review_manifest.json")))
        if not relative.startswith("data/outputs/figures/") or not manifest_name:
            relative = ""
        for entry in tomllib.loads(registry.read_text(encoding="utf-8")).get("manifest_archives", []):
            if not relative:
                break
            archive = root / entry["archive"]
            stat = archive.stat()
            contents = _verified_archive(str(archive), entry["sha256"], stat.st_mtime_ns, stat.st_size)
            payload = contents.get(relative)
            if payload is not None and hashlib.sha256(payload).hexdigest() == expected:
                return payload
    if current is not None:
        return current
    raise FileNotFoundError(value)


def read_manifest_record(record: dict[str, Any], *, root: Path) -> dict:
    payload = record_bytes(record, root=root)
    if record.get("sha256") and hashlib.sha256(payload).hexdigest() != record["sha256"]:
        raise ValueError(f"manifest reference hash mismatch: {record['path']}")
    manifest = json.loads(payload)
    if record.get("figure_id"):
        from scripts.plotting.plot_bundle import expand_bundle
        matches = [item for entry, (_, item) in zip(manifest["figures"], expand_bundle(manifest))
                   if entry["figure_id"] == record["figure_id"]]
        if len(matches) != 1:
            raise ValueError("source figure ID is not unique in its bundle")
        return matches[0]
    return manifest


@lru_cache(maxsize=512)
def git_source(root: str, commit: str, relative: str) -> bytes:
    if not re.fullmatch(r"[0-9a-fA-F]{40}", commit):
        raise ValueError("code_ref must be a full Git commit SHA")
    registry = Path(root) / "config/plotting/historical_snapshots.toml"
    if registry.is_file():
        registered = tomllib.loads(registry.read_text(encoding="utf-8"))
        entries = registered.get("code_snapshots", []) + [
            {**entry, "code_commit": entry.get("source_commit")}
            for entry in registered.get("manifest_archives", [])
        ]
        for entry in entries:
            if entry.get("code_commit") != commit:
                continue
            archive = Path(root) / entry["archive"]
            try:
                if hashlib.sha256(archive.read_bytes()).hexdigest() != entry["sha256"]:
                    raise ValueError(f"historical code archive hash mismatch: {archive}")
                with zipfile.ZipFile(archive) as snapshot:
                    return snapshot.read(relative)
            except KeyError:
                # Several archives can jointly retain one historical dependency closure.
                continue
            except (OSError, zipfile.BadZipFile) as exc:
                raise ValueError(f"historical code archive unavailable: {archive}") from exc
    try:
        return subprocess.check_output(["git", "show", f"{commit}:{relative}"], cwd=root, stderr=subprocess.PIPE)
    except (OSError, subprocess.CalledProcessError) as exc:
        raise ValueError(f"historical source unavailable: {commit}:{relative}") from exc


def is_code_record(path: Path, root: Path) -> bool:
    try:
        rel = path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return False
    return ((rel.startswith("scripts/") and path.suffix in {".py", ".jl", ".sh", ".ps1"})
            or (rel.startswith("config/plotting/") and path.suffix == ".toml")
            or (rel.startswith("docs/guides/sop/") and path.suffix == ".md")
            or (rel.startswith(".agents/skills/") and path.name == "SKILL.md")
            or (rel.startswith("docs/analysis/") and path.name == "plotting_case_contract.md"))


def working_tree_source(path: Path, expected: str, root: Path) -> bytes | None:
    """Resolve exact source bytes without treating a dirty base commit as authority."""
    registry = root / "config/plotting/historical_snapshots.toml"
    if not registry.is_file() or not is_code_record(path, root):
        return None
    relative = path.resolve().relative_to(root.resolve()).as_posix()
    for entry in tomllib.loads(registry.read_text(encoding="utf-8")).get("working_tree_source_archives", []):
        archive = root / entry["archive"]
        stat = archive.stat()
        contents = _verified_archive(str(archive), entry["sha256"], stat.st_mtime_ns, stat.st_size)
        payload = contents.get(relative)
        if payload is not None and hashlib.sha256(payload).hexdigest() == expected.lower():
            return payload
    return None


def code_ref_for_manifest(manifest_path: Path, root: Path) -> str | None:
    registry = root / "config/plotting/historical_snapshots.toml"
    if not registry.is_file():
        return None
    try:
        entries = tomllib.loads(registry.read_text(encoding="utf-8"))
    except (OSError, tomllib.TOMLDecodeError):
        return None
    name = manifest_path.resolve().as_posix()
    for entry in entries.values():
        prefix = entry.get("manifest_prefix") if isinstance(entry, dict) else None
        commit = entry.get("code_commit") if isinstance(entry, dict) else None
        if prefix and commit and prefix in name:
            return str(commit)
    return None


def validate_hash_record(record: dict[str, Any], *, root: Path, label: str,
                         code_ref: str | None = None, allow_historical: bool = True) -> list[str]:
    value = record.get("path")
    if not isinstance(value, str) or not value:
        return [f"{label}.path must be a non-empty string"]
    path = Path(value)
    path = path.resolve() if path.is_absolute() else (root / path).resolve()
    expected = record.get("sha256")
    if not isinstance(expected, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", expected):
        return [f"{label}.sha256 must be a 64-character hash"]
    try:
        payload = record_bytes(record, root=root)
    except FileNotFoundError:
        payload = None
    except (OSError, ValueError, zipfile.BadZipFile) as exc:
        return [f"{label} cannot read {value}: {exc}"]

    def matches(data: bytes | None) -> bool:
        return data is not None and hashlib.sha256(data).hexdigest() == expected.lower()
    if not matches(payload) and allow_historical and (label.startswith("generator") or is_code_record(path, root)):
        snapshot = record.get("source_snapshot")
        if isinstance(snapshot, dict):
            errors = validate_hash_record(snapshot, root=root, label=f"{label}.source_snapshot", allow_historical=False)
            if errors:
                return errors
            snapshot_path = Path(snapshot["path"])
            payload = (snapshot_path if snapshot_path.is_absolute() else root / snapshot_path).read_bytes()
        elif is_code_record(path, root):
            try:
                retained = working_tree_source(path, expected, root)
                if retained is not None:
                    payload = retained
                elif record.get("git_commit") or code_ref:
                    commit = record.get("git_commit") or code_ref
                    payload = git_source(str(root.resolve()), commit, path.relative_to(root.resolve()).as_posix())
            except (OSError, ValueError, zipfile.BadZipFile) as exc:
                return [f"{label}: {exc}"]
    if payload is None:
        return [f"{label} missing file: {value}"]
    errors = []
    if not matches(payload):
        errors.append(f"{label}.sha256 mismatch for {value}")
    if "bytes" in record and record["bytes"] != len(payload):
        errors.append(f"{label}.bytes mismatch for {value}")
    return errors


def validate_snapshot(manifest_path: Path, *, root: Path, code_ref: str) -> list[str]:
    """Check a committed manifest graph; only code/contracts may use Git bytes.

    This verifies retained evidence. It does not execute a historical renderer
    or promote an artifact under the current style/numerical contract.
    """
    path = manifest_path.resolve()
    try:
        relative = path.relative_to(root.resolve()).as_posix()
        frozen = git_source(str(root.resolve()), code_ref, relative)
        if record_bytes({"path": relative, "sha256": hashlib.sha256(frozen).hexdigest()}, root=root) != frozen:
            return [f"snapshot manifest differs from {code_ref}:{relative}"]
    except (ValueError, OSError, zipfile.BadZipFile) as exc:
        return [str(exc)]
    errors: list[str] = []
    visited: set[tuple[Path, str | None]] = set()
    checked: dict[tuple, bool] = {}

    def read_manifest(target: Path, expected: str | None = None) -> None:
        target = target.resolve()
        key = (target, expected)
        if key in visited:
            return
        visited.add(key)
        try:
            walk(read_manifest_record({"path": str(target), "sha256": expected}, root=root), target.name)
        except (OSError, ValueError) as exc:
            errors.append(f"cannot read snapshot manifest {target}: {exc}")

    def check(record: dict, label: str, *, historical: bool = True) -> bool:
        key = (record.get("path"), record.get("sha256"), record.get("bytes"), historical)
        if key not in checked:
            issues = validate_hash_record(record, root=root, label=label, code_ref=code_ref, allow_historical=historical)
            errors.extend(issues)
            checked[key] = not issues
        value = record.get("path", "")
        follow = "outputs" in label or label.endswith(("_package", ".source_png_manifest"))
        if checked[key] and follow and value.endswith(".json") and "manifest" in Path(value).name:
            read_manifest(root / value, record.get("sha256"))
        return checked[key]

    def walk(node: Any, label: str) -> None:
        if isinstance(node, list):
            for index, item in enumerate(node):
                walk(item, f"{label}[{index}]")
        elif isinstance(node, dict):
            if "path" in node and "sha256" in node:
                check(node, label, historical="outputs" not in label)
                return
            if "manifest" in node and "manifest_sha256" in node:
                chart = {"path": node["manifest"], "sha256": node["manifest_sha256"]}
                if check(chart, label, historical=False):
                    read_manifest(root / chart["path"], chart["sha256"])
            if isinstance(node.get("generator"), str) and "generator_sha256" in node:
                check({"path": node["generator"], "sha256": node["generator_sha256"]}, f"{label}.generator")
            for key, value in node.items():
                walk(value, f"{label}.{key}")

    read_manifest(path, hashlib.sha256(frozen).hexdigest())
    return errors
