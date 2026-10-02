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


@lru_cache(maxsize=512)
def git_source(root: str, commit: str, relative: str) -> bytes:
    if not re.fullmatch(r"[0-9a-fA-F]{40}", commit):
        raise ValueError("code_ref must be a full Git commit SHA")
    registry = Path(root) / "config/plotting/historical_snapshots.toml"
    if registry.is_file():
        entries = tomllib.loads(registry.read_text(encoding="utf-8")).get("code_snapshots", [])
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
                # A registered snapshot contains only its retained dependency closure.
                break
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
            or (rel.startswith(".agents/skills/") and path.name == "SKILL.md"))


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
        payload = path.read_bytes() if path.is_file() else None
    except OSError as exc:
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
        elif (record.get("git_commit") or code_ref) and is_code_record(path, root):
            commit = record.get("git_commit") or code_ref
            try:
                payload = git_source(str(root.resolve()), commit, path.relative_to(root.resolve()).as_posix())
            except ValueError as exc:
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
        if not path.is_file() or path.read_bytes() != frozen:
            return [f"snapshot manifest differs from {code_ref}:{relative}"]
    except ValueError as exc:
        return [str(exc)]
    errors: list[str] = []
    visited: set[Path] = set()

    def read_manifest(target: Path) -> None:
        target = target.resolve()
        if target in visited:
            return
        visited.add(target)
        try:
            walk(json.loads(target.read_text(encoding="utf-8")), target.name)
        except (OSError, json.JSONDecodeError) as exc:
            errors.append(f"cannot read snapshot manifest {target}: {exc}")

    def check(record: dict, label: str, *, historical: bool = True) -> bool:
        issues = validate_hash_record(record, root=root, label=label, code_ref=code_ref, allow_historical=historical)
        errors.extend(issues)
        value = record.get("path", "")
        follow = "outputs" in label or label.endswith(("_package", ".source_png_manifest"))
        if not issues and follow and value.endswith(".json") and "manifest" in Path(value).name:
            read_manifest(root / value)
        return not issues

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
                    read_manifest(root / chart["path"])
            if isinstance(node.get("generator"), str) and "generator_sha256" in node:
                check({"path": node["generator"], "sha256": node["generator_sha256"]}, f"{label}.generator")
            for key, value in node.items():
                walk(value, f"{label}.{key}")

    read_manifest(path)
    return errors
