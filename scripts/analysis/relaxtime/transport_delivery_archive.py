"""Pack or restore an explicit transport evidence file list without overwriting files."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import stat
import zipfile

MANIFEST = "__transport_delivery__/manifest.json"
CHUNK = 1024 * 1024


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(CHUNK), b""):
            value.update(block)
    return value.hexdigest()


def safe_name(name: str) -> str:
    path = PurePosixPath(name)
    if (not name or "\\" in name or ":" in name or path.is_absolute()
            or any(part in {"", ".", ".."} for part in name.split("/"))
            or any(part.endswith((".", " ")) for part in path.parts)):
        raise ValueError(f"unsafe archive path: {name!r}")
    return name


def validate_entries(entries: list[dict]) -> None:
    names = [safe_name(row["path"]) for row in entries]
    if len({name.casefold() for name in names}) != len(names):
        raise ValueError("duplicate or case-colliding archive path")
    if MANIFEST.casefold() in {name.casefold() for name in names}:
        raise ValueError("reserved manifest path")


def verify(archive: Path, receipt: dict) -> dict:
    if archive.stat().st_size != receipt["archive_bytes"] or digest(archive) != receipt["archive_sha256"]:
        raise ValueError("archive identity mismatch")
    with zipfile.ZipFile(archive) as bundle:
        infos = bundle.infolist()
        if len({info.filename.casefold() for info in infos}) != len(infos):
            raise ValueError("duplicate zip member")
        raw = bundle.read(MANIFEST)
        if hashlib.sha256(raw).hexdigest() != receipt["manifest_sha256"]:
            raise ValueError("manifest identity mismatch")
        manifest = json.loads(raw)
        validate_entries(manifest["files"])
        expected = {row["path"] for row in manifest["files"]} | {MANIFEST}
        if {info.filename for info in infos} != expected:
            raise ValueError("archive members differ from manifest")
        for info in infos:
            if stat.S_ISLNK(info.external_attr >> 16) or info.is_dir():
                raise ValueError("archive must contain only regular files")
        for row in manifest["files"]:
            value, size = hashlib.sha256(), 0
            with bundle.open(row["path"]) as stream:
                for block in iter(lambda: stream.read(CHUNK), b""):
                    size += len(block)
                    value.update(block)
            if size != row["bytes"] or value.hexdigest() != row["sha256"]:
                raise ValueError(f"member identity mismatch: {row['path']}")
    return manifest


def pack(root: Path, names: list[str], archive: Path, receipt_path: Path) -> dict:
    root = root.resolve()
    validate_entries([{"path": name} for name in names])
    if archive.exists() or receipt_path.exists():
        raise FileExistsError("archive and receipt must be new files")
    sources = []
    for name in sorted(names):
        source = root / name
        if source.is_symlink() or not source.resolve().is_relative_to(root) or not source.is_file():
            raise ValueError(f"not a regular file inside source root: {name}")
        sources.append((name, source))
    if archive.resolve() in {source.resolve() for _, source in sources}:
        raise ValueError("archive cannot contain itself")
    archive.parent.mkdir(parents=True, exist_ok=True)
    records = []
    with zipfile.ZipFile(archive, "x", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as bundle:
        for index, (name, source) in enumerate(sources, 1):
            value, size = hashlib.sha256(), 0
            with source.open("rb") as stream, bundle.open(name, "w", force_zip64=True) as target:
                for block in iter(lambda: stream.read(CHUNK), b""):
                    target.write(block)
                    size += len(block)
                    value.update(block)
            if source.stat().st_size != size or digest(source) != value.hexdigest():
                raise ValueError(f"source changed while packing: {name}")
            records.append(dict(path=name, bytes=size, sha256=value.hexdigest()))
            if index % 5000 == 0:
                print(f"packed {index}/{len(sources)} files", flush=True)
        manifest = dict(schema="transport_delivery_archive_v1", created_utc=datetime.now(timezone.utc).isoformat(),
                        source_root=root.as_posix(), files=records)
        raw = (json.dumps(manifest, ensure_ascii=False, indent=2) + "\n").encode("utf-8")
        bundle.writestr(MANIFEST, raw)
    receipt = dict(schema="transport_delivery_archive_receipt_v1", archive=archive.resolve().as_posix(),
                   archive_bytes=archive.stat().st_size, archive_sha256=digest(archive),
                   manifest_sha256=hashlib.sha256(raw).hexdigest(), files=len(records),
                   uncompressed_bytes=sum(row["bytes"] for row in records))
    verify(archive, receipt)
    receipt["archive_verified"] = True
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    with receipt_path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(receipt, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    return receipt


def restore(archive: Path, receipt: dict, destination: Path) -> dict:
    if destination.exists():
        raise FileExistsError("restore destination must not exist")
    manifest = verify(archive, receipt)
    destination.mkdir(parents=True)
    resolved = destination.resolve()
    with zipfile.ZipFile(archive) as bundle:
        for index, row in enumerate(manifest["files"], 1):
            target = destination / safe_name(row["path"])
            if not target.resolve().is_relative_to(resolved):
                raise ValueError("restore path escapes destination")
            target.parent.mkdir(parents=True, exist_ok=True)
            with bundle.open(row["path"]) as stream, target.open("xb") as output:
                shutil.copyfileobj(stream, output, length=CHUNK)
            if target.stat().st_size != row["bytes"] or digest(target) != row["sha256"]:
                raise ValueError(f"restored file mismatch: {row['path']}")
            if index % 5000 == 0:
                print(f"restored {index}/{len(manifest['files'])} files", flush=True)
    return dict(status="restored_and_verified", files=len(manifest["files"]),
                destination=resolved.as_posix(), archive_sha256=receipt["archive_sha256"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("pack", "verify", "restore"))
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--paths", type=Path, help="UTF-8 file with one repository-relative file per line")
    parser.add_argument("--destination", type=Path)
    args = parser.parse_args()
    if args.action == "pack":
        if args.paths is None:
            parser.error("pack requires --paths")
        names = args.paths.read_text(encoding="utf-8").splitlines()
        result = pack(args.root, names, args.archive, args.receipt)
    else:
        receipt = json.loads(args.receipt.read_text(encoding="utf-8"))
        if args.action == "restore":
            if args.destination is None:
                parser.error("restore requires --destination")
            result = restore(args.archive, receipt, args.destination)
        else:
            manifest = verify(args.archive, receipt)
            result = dict(status="verified", files=len(manifest["files"]))
    print(json.dumps(result, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
