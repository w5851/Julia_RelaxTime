#!/usr/bin/env python3
"""Archive superseded v6-v9 review files without deleting or rewriting them."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import zipfile

ROOT = Path(__file__).resolve().parents[3]
INDEX = ROOT / "docs/analysis/relaxtime/phase_guided_transport/publication_review_history_v6_v9_archive_v1.json"
FIGURES = "data/outputs/figures/relaxtime/transport/phase_guided"
ANALYSIS = "docs/analysis/relaxtime/phase_guided_transport"
HISTORY_DIRECTORIES = (
    *(f"{FIGURES}/publication_clean_{version}" for version in ("v6", "v7", "v7_png_review", "v8_png_review", "v9_png_review")),
    *(f"{ANALYSIS}/phase_guided_transport_publication_clean_{version}" for version in
      ("figure_layer_v6", "v6", "v7", "v7_png_review", "v8_png_review", "v9_png_review")),
)
RETAINED_DEPENDENCIES = tuple(f"scripts/analysis/relaxtime/build_phase_guided_publication_clean_v{version}.py" for version in range(6, 10))
HISTORY_ONLY_FILES = (
    "scripts/analysis/relaxtime/sync_phase_guided_publication_clean_figure_layer_v6.py",
    *(f"tests/unit/python/test_phase_guided_publication_clean_v{version}.py" for version in range(6, 10)),
)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def safe_member(name: str) -> bool:
    path = PurePosixPath(name)
    return bool(name) and not path.is_absolute() and ".." not in path.parts and "\\" not in name and ":" not in name


def verify_archive(archive: Path, records: list[dict]) -> None:
    expected = {record["path"]: record for record in records}
    if len(expected) != len(records):
        raise ValueError("archive index contains duplicate paths")
    with zipfile.ZipFile(archive) as bundle:
        names = bundle.namelist()
        if len(names) != len(set(names)) or set(names) != set(expected):
            raise ValueError("archive members do not match the exact index")
        for name in names:
            if not safe_member(name):
                raise ValueError(f"unsafe archive member: {name}")
            data = bundle.read(name)
            if len(data) != expected[name]["bytes"] or digest(data) != expected[name]["sha256"]:
                raise ValueError(f"archive byte/hash mismatch: {name}")


def create_archive(destination: Path) -> None:
    if INDEX.exists() or destination.exists():
        raise FileExistsError("refusing to overwrite an archive or its published index")
    paths = []
    for relative in HISTORY_DIRECTORIES:
        directory = ROOT / relative
        if not directory.is_dir():
            raise FileNotFoundError(directory)
        paths.extend(path for path in directory.rglob("*") if path.is_file())
    paths.extend(ROOT / relative for relative in (*RETAINED_DEPENDENCIES, *HISTORY_ONLY_FILES))
    paths = sorted(set(paths))
    records = [{"path": path.relative_to(ROOT).as_posix(), "bytes": path.stat().st_size,
                "sha256": digest(path.read_bytes())} for path in paths]
    drift = []
    for relative in HISTORY_DIRECTORIES:
        manifest = ROOT / relative / "manifest.json"
        if not manifest.is_file():
            continue
        payload = json.loads(manifest.read_text(encoding="utf-8"))
        for record in payload.get("inputs", []):
            path = ROOT / record["path"]
            if path.is_file() and digest(path.read_bytes()) != record["sha256"]:
                drift.append({"manifest": manifest.relative_to(ROOT).as_posix(), "input": record["path"],
                              "recorded_sha256": record["sha256"], "current_sha256": digest(path.read_bytes())})
    destination.mkdir(parents=True, exist_ok=False)
    archive = destination / "publication_review_history_v6_v9.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as bundle:
        for path in paths:
            bundle.write(path, path.relative_to(ROOT).as_posix())
    verify_archive(archive, records)
    for record in records:
        if digest((ROOT / record["path"]).read_bytes()) != record["sha256"]:
            raise ValueError(f"source changed while archiving: {record['path']}")
    index = {
        "schema": "publication_review_history_v6_v9_archive_v1", "status": "locally_verified_restorable",
        "scope": "superseded review snapshots; not current plotting or numerical gate evidence",
        "archive": {"path": str(archive.resolve()), "bytes": archive.stat().st_size, "sha256": digest(archive.read_bytes())},
        "remote_archive_available": False, "member_count": len(records), "uncompressed_bytes": sum(record["bytes"] for record in records),
        "members": records, "source_directories": list(HISTORY_DIRECTORIES),
        "retained_renderer_dependencies": list(RETAINED_DEPENDENCIES), "history_only_files": list(HISTORY_ONLY_FILES),
        "observed_input_contract_drift": drift, "sources_deleted": False,
        "restore_command": "python scripts/analysis/relaxtime/archive_phase_guided_publication_review_history.py restore --archive <zip> --destination <new-empty-directory>",
    }
    text = json.dumps(index, indent=2) + "\n"
    INDEX.write_text(text, encoding="utf-8")
    (destination / "archive_index.json").write_text(text, encoding="utf-8")
    print(json.dumps({"archive": str(archive), "member_count": len(records), "sha256": index["archive"]["sha256"]}))


def restore_archive(archive: Path, destination: Path, index_path: Path) -> None:
    index = json.loads(index_path.read_text(encoding="utf-8"))
    if digest(archive.read_bytes()) != index["archive"]["sha256"]:
        raise ValueError("whole archive SHA-256 mismatch")
    verify_archive(archive, index["members"])
    if destination.exists():
        raise FileExistsError("restore requires a new absent destination")
    destination.mkdir(parents=True, exist_ok=False)
    root = destination.resolve()
    with zipfile.ZipFile(archive) as bundle:
        for record in index["members"]:
            target = (root / record["path"]).resolve()
            if not target.is_relative_to(root):
                raise ValueError(f"restore target escapes destination: {target}")
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(bundle.read(record["path"]))
    for record in index["members"]:
        if digest((root / record["path"]).read_bytes()) != record["sha256"]:
            raise ValueError(f"restored hash mismatch: {record['path']}")
    print(f"[review-history] restored {len(index['members'])} exact files to {root}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    create = commands.add_parser("create")
    create.add_argument("--destination", type=Path, required=True)
    restore = commands.add_parser("restore")
    restore.add_argument("--archive", type=Path, required=True)
    restore.add_argument("--destination", type=Path, required=True)
    restore.add_argument("--index", type=Path, default=INDEX)
    args = parser.parse_args()
    if args.command == "create":
        create_archive(args.destination)
    else:
        restore_archive(args.archive, args.destination, args.index)


if __name__ == "__main__":
    main()
