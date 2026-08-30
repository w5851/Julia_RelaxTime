"""Build a read-only inventory for historical figure asset retirement.

The inventory records existing PNG/PDF/SVG files and proposes review groups.
It deliberately has no delete, move, or overwrite operation for figure assets.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import struct
import subprocess
import sys
from typing import Any, Iterable
import xml.etree.ElementTree as ET


_SCRIPT_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_SCRIPT_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_PROJECT_ROOT))

from scripts.plotting.plot_manifest import PROJECT_ROOT, git_commit, sha256_file


REGISTRY_SCHEMA = "figure_asset_registry_v1"
ASSET_SUFFIXES = frozenset({".png", ".pdf", ".svg"})
REFERENCE_TEXT_SUFFIXES = frozenset({".csv", ".jl", ".json", ".md", ".py", ".toml", ".txt", ".yaml", ".yml"})
DEFAULT_REGISTRY = Path("docs/analysis/governance/figure_asset_registry_v1/asset_registry.json")
DEFAULT_CANDIDATES = Path("docs/analysis/governance/figure_asset_registry_v1/cleanup_candidates.csv")
VARIANT_SUFFIX_RE = re.compile(r"(?:__line_first_v\d+|__pilot_v\d+|__new_review.*)$")


def _relative_path(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def _resolve_path(path: Path, root: Path) -> Path:
    return path if path.is_absolute() else (root / path)


def _git_tracked_paths(repo_root: Path, asset_root: Path) -> set[Path]:
    root_relative = _relative_path(asset_root, repo_root)
    try:
        raw = subprocess.check_output(
            ["git", "ls-files", "-z", "--", root_relative],
            cwd=repo_root,
        )
    except (OSError, subprocess.CalledProcessError):
        return set()
    return {
        (repo_root / Path(item)).resolve()
        for item in raw.decode("utf-8").split("\0")
        if item
    }


def discover_assets(
    asset_root: Path,
    *,
    repo_root: Path,
    tracked_only: bool,
) -> tuple[list[Path], int, int]:
    """Return selected files, visible untracked exclusions, and missing tracked files."""

    on_disk = sorted(
        (
            path
            for path in asset_root.rglob("*")
            if path.is_file() and path.suffix.lower() in ASSET_SUFFIXES
        ),
        key=lambda path: path.as_posix().lower(),
    )
    if not tracked_only:
        return on_disk, 0, 0

    tracked = _git_tracked_paths(repo_root, asset_root)
    selected = [path for path in on_disk if path.resolve() in tracked]
    excluded = len([path for path in on_disk if path.resolve() not in tracked])
    missing = len([path for path in tracked if path.suffix.lower() in ASSET_SUFFIXES and not path.is_file()])
    return selected, excluded, missing


def _nearest_manifest(asset_path: Path, asset_root: Path) -> Path | None:
    current = asset_path.parent.resolve()
    root = asset_root.resolve()
    while True:
        candidate = current / "plot_manifest.json"
        if candidate.is_file():
            return candidate
        if current == root or root not in current.parents:
            return None
        current = current.parent


def _load_manifest(path: Path | None) -> tuple[dict[str, Any] | None, str | None]:
    if path is None:
        return None, None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return None, str(exc)
    if not isinstance(payload, dict):
        return None, "manifest root is not an object"
    return payload, None


def _parse_png(path: Path) -> dict[str, Any]:
    signature = b"\x89PNG\r\n\x1a\n"
    metadata: dict[str, Any] = {"vector": False}
    try:
        with path.open("rb") as handle:
            if handle.read(8) != signature:
                return {**metadata, "signature_valid": False}
            metadata["signature_valid"] = True
            width, height = struct.unpack(">II", handle.read(8))
            metadata.update({"width_px": width, "height_px": height})
            handle.seek(8)
            while True:
                header = handle.read(8)
                if len(header) != 8:
                    break
                length, chunk_type = struct.unpack(">I4s", header)
                chunk = handle.read(length)
                handle.seek(4, 1)
                if chunk_type == b"pHYs" and len(chunk) >= 9 and chunk[8] == 1:
                    x_ppu = struct.unpack(">I", chunk[0:4])[0]
                    y_ppu = struct.unpack(">I", chunk[4:8])[0]
                    metadata["dpi"] = {
                        "x": round(x_ppu * 0.0254, 3),
                        "y": round(y_ppu * 0.0254, 3),
                    }
                    break
                if chunk_type == b"IEND":
                    break
    except (OSError, struct.error) as exc:
        metadata["inspection_error"] = str(exc)
    return metadata


def _parse_svg(path: Path) -> dict[str, Any]:
    metadata: dict[str, Any] = {"vector": True}
    try:
        root = ET.fromstring(path.read_text(encoding="utf-8"))
        tag = root.tag.rsplit("}", 1)[-1]
        metadata["signature_valid"] = tag == "svg"
        metadata["width"] = root.attrib.get("width")
        metadata["height"] = root.attrib.get("height")
        metadata["viewBox"] = root.attrib.get("viewBox")
        if tag != "svg":
            metadata["inspection_error"] = "root element is not svg"
    except (OSError, UnicodeError, ET.ParseError) as exc:
        metadata["signature_valid"] = False
        metadata["inspection_error"] = str(exc)
    return metadata


def _parse_pdf(path: Path) -> dict[str, Any]:
    try:
        with path.open("rb") as handle:
            signature_valid = handle.read(5) == b"%PDF-"
    except OSError as exc:
        return {"vector": None, "signature_valid": False, "inspection_error": str(exc)}
    return {"vector": None, "signature_valid": signature_valid}


def inspect_file(path: Path) -> dict[str, Any]:
    suffix = path.suffix.lower()
    if suffix == ".png":
        metadata = _parse_png(path)
    elif suffix == ".svg":
        metadata = _parse_svg(path)
    else:
        metadata = _parse_pdf(path)
    return {
        "format": suffix.removeprefix("."),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "metadata": metadata,
    }


def _case_group(case_directory: Path, asset_root: Path, repo_root: Path) -> str:
    case_name = case_directory.name
    if "__plotv1__" in case_name:
        case_name = case_name.split("__plotv1__", 1)[0] + "__plotv1__variants"
    case_name = VARIANT_SUFFIX_RE.sub("", case_name)
    parent = case_directory.parent
    return f"{_relative_path(parent, repo_root)}/{case_name}".replace("//", "/")


def _manifest_record(manifest_path: Path | None, *, repo_root: Path) -> dict[str, Any] | None:
    if manifest_path is None:
        return None
    payload, parse_error = _load_manifest(manifest_path)
    record: dict[str, Any] = {
        "path": _relative_path(manifest_path, repo_root),
        "bytes": manifest_path.stat().st_size,
        "sha256": sha256_file(manifest_path),
    }
    if parse_error is not None:
        record["parse_error"] = parse_error
        return record
    assert payload is not None
    for field in (
        "schema_version",
        "figure_mode",
        "style_profile",
        "semantic_status",
        "publication_scope",
    ):
        if field in payload:
            record[field] = payload[field]
    return record


def _classify(manifest: dict[str, Any] | None, manifest_record: dict[str, Any] | None) -> str:
    if manifest_record is None:
        return "legacy_unregistered"
    if manifest_record.get("parse_error"):
        return "invalid_manifest"
    if manifest_record.get("schema_version") == "plot_manifest_v1":
        return "contract_case"
    return "legacy_manifest"


def _sidecar_record(asset_path: Path, *, repo_root: Path) -> dict[str, Any] | None:
    sidecar = Path(str(asset_path) + ".provenance.json")
    if not sidecar.is_file():
        return None
    return {
        "path": _relative_path(sidecar, repo_root),
        "bytes": sidecar.stat().st_size,
        "sha256": sha256_file(sidecar),
    }


def _repository_reference_scan(
    asset_paths: Iterable[Path],
    *,
    repo_root: Path,
) -> dict[str, list[str]]:
    """Find best-effort references in tracked text files, excluding figure binaries."""

    try:
        raw = subprocess.check_output(["git", "ls-files", "-z"], cwd=repo_root)
    except (OSError, subprocess.CalledProcessError):
        return {path.resolve().as_posix(): [] for path in asset_paths}

    repository_files = [
        (repo_root / Path(item)).resolve()
        for item in raw.decode("utf-8").split("\0")
        if item and Path(item).suffix.lower() in REFERENCE_TEXT_SUFFIXES
    ]
    asset_list = list(asset_paths)
    relative_paths = {_relative_path(path, repo_root): path.resolve().as_posix() for path in asset_list}
    basename_counts: dict[str, int] = defaultdict(int)
    for path in asset_list:
        basename_counts[path.name] += 1
    token_to_assets: dict[str, set[str]] = defaultdict(set)
    for path in asset_list:
        relative = _relative_path(path, repo_root)
        token_to_assets[relative].add(relative)
        token_to_assets[relative.replace("/", "\\")].add(relative)
        if basename_counts[path.name] == 1 and len(path.name) >= 5:
            token_to_assets[path.name].add(relative)

    references: dict[str, set[str]] = {relative: set() for relative in relative_paths}
    for source_path in repository_files:
        source_relative = _relative_path(source_path, repo_root)
        if source_relative.startswith("data/outputs/figures/"):
            continue
        if source_relative.startswith("docs/analysis/governance/figure_asset_registry_v1/"):
            continue
        try:
            content = source_path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for token, assets in token_to_assets.items():
            if token in content:
                for relative in assets:
                    references[relative].add(source_relative)
    return {relative_paths[relative]: sorted(files) for relative, files in references.items()}


def _manual_reasons(
    *,
    classification: str,
    case_group: str,
    group_case_counts: dict[str, int],
) -> list[str]:
    reasons: list[str] = []
    if classification != "contract_case":
        reasons.append("legacy_or_unregistered_manifest")
    if group_case_counts[case_group] > 1:
        reasons.append("multiple_case_variants")
    reasons.append("external_references_unknown")
    return reasons


def build_registry(
    *,
    asset_root: Path,
    repo_root: Path = PROJECT_ROOT,
    tracked_only: bool = True,
    generated_at_utc: str | None = None,
) -> dict[str, Any]:
    """Build a registry and review-only classifications without mutating assets."""

    asset_root = asset_root.resolve()
    repo_root = repo_root.resolve()
    paths, excluded_count, missing_count = discover_assets(
        asset_root,
        repo_root=repo_root,
        tracked_only=tracked_only,
    )
    tracked_paths = _git_tracked_paths(repo_root, asset_root)
    preliminary: list[dict[str, Any]] = []
    for path in paths:
        manifest_path = _nearest_manifest(path, asset_root)
        manifest, _ = _load_manifest(manifest_path)
        case_directory = manifest_path.parent if manifest_path is not None else path.parent
        group = _case_group(case_directory, asset_root, repo_root)
        manifest_record = _manifest_record(manifest_path, repo_root=repo_root)
        classification = _classify(manifest, manifest_record)
        record = {
            "path": _relative_path(path, repo_root),
            "case_directory": _relative_path(case_directory, repo_root),
            "review_group": group,
            "tracked": path.resolve() in tracked_paths,
            **inspect_file(path),
            "manifest": manifest_record,
            "provenance_sidecar": _sidecar_record(path, repo_root=repo_root),
            "proposed_classification": classification,
        }
        preliminary.append(record)

    repository_references = _repository_reference_scan(
        [path for path in paths],
        repo_root=repo_root,
    )
    for record in preliminary:
        references = repository_references.get(
            (repo_root / Path(record["path"])).resolve().as_posix(),
            [],
        )
        record["repository_references"] = references
        record["repository_reference_count"] = len(references)

    group_case_directories: dict[str, set[str]] = defaultdict(set)
    for record in preliminary:
        group_case_directories[record["review_group"]].add(record["case_directory"])
    group_case_counts = {group: len(case_directories) for group, case_directories in group_case_directories.items()}

    assets: list[dict[str, Any]] = []
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in preliminary:
        reasons = _manual_reasons(
            classification=record["proposed_classification"],
            case_group=record["review_group"],
            group_case_counts=group_case_counts,
        )
        record["review_required"] = bool(reasons)
        record["review_reasons"] = reasons
        record["proposed_action"] = "owner_review_only" if reasons else "keep_contract_case"
        assets.append(record)
        groups[record["review_group"]].append(record)

    case_groups: list[dict[str, Any]] = []
    review_groups: list[dict[str, Any]] = []
    for group_name in sorted(groups):
        group_assets = groups[group_name]
        classifications = sorted({item["proposed_classification"] for item in group_assets})
        manifest_records = [item["manifest"] for item in group_assets if item["manifest"] is not None]
        modes = sorted({str(item["figure_mode"]) for item in manifest_records if "figure_mode" in item})
        styles = sorted({str(item["style_profile"]) for item in manifest_records if "style_profile" in item})
        group_record = {
            "review_group": group_name,
            "case_directories": sorted({item["case_directory"] for item in group_assets}),
            "asset_count": len(group_assets),
            "formats": sorted({item["format"] for item in group_assets}),
            "manifest_schemas": sorted({str(item["schema_version"]) for item in manifest_records if "schema_version" in item}),
            "figure_modes": modes,
            "style_profiles": styles,
            "proposed_classifications": classifications,
            "assets_with_repository_references": sum(
                item["repository_reference_count"] > 0 for item in group_assets
            ),
            "repository_reference_count": sum(
                item["repository_reference_count"] for item in group_assets
            ),
            "review_required": any(item["review_required"] for item in group_assets),
            "proposed_action": "owner_review_only" if any(item["review_required"] for item in group_assets) else "keep_contract_case",
            "manual_questions": [
                "Is this group referenced by the current paper, supplement, report, or slides outside the repository?",
                "Which case or variant is canonical if multiple renders exist?",
                "Should any asset remain legacy, be quarantined, or be deleted in a later approved cleanup PR?",
            ],
        }
        case_groups.append(group_record)
        if group_record["review_required"]:
            review_groups.append(group_record)

    return {
        "schema_version": REGISTRY_SCHEMA,
        "generated_at_utc": generated_at_utc or datetime.now(timezone.utc).isoformat(),
        "source_git_commit": git_commit(repo_root),
        "generator": {
            "path": _relative_path(Path(__file__).resolve(), repo_root),
            "sha256": sha256_file(Path(__file__).resolve()),
        },
        "asset_root": _relative_path(asset_root, repo_root),
        "selection": {
            "tracked_only": tracked_only,
            "asset_extensions": sorted(ASSET_SUFFIXES),
            "untracked_policy": "excluded_and_not_modified" if tracked_only else "included_for_local_review_only",
            "diagnostic_roots": ["docs/analysis"],
            "external_reference_status": "unknown_outside_repository",
            "repository_reference_scan": "git_tracked_text_files_excluding_figure_binaries",
        },
        "protection": {
            "delete_performed": False,
            "move_performed": False,
            "formal_assets_overwritten": False,
            "cleanup_authorization": "not_granted",
        },
        "discovery": {
            "asset_count": len(assets),
            "case_group_count": len(case_groups),
            "review_group_count": len(review_groups),
            "assets_with_repository_references": sum(
                item["repository_reference_count"] > 0 for item in assets
            ),
            "repository_reference_count": sum(
                item["repository_reference_count"] for item in assets
            ),
            "visible_untracked_excluded": excluded_count,
            "missing_tracked_assets": missing_count,
        },
        "assets": assets,
        "case_groups": case_groups,
        "review_groups": review_groups,
    }


def _write_text(path: Path, content: str, *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"refusing to overwrite report: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(content, encoding="utf-8", newline="\n")
    temporary.replace(path)


def write_registry(path: Path, registry: dict[str, Any], *, overwrite: bool) -> None:
    _write_text(path, json.dumps(registry, ensure_ascii=False, indent=2) + "\n", overwrite=overwrite)


def write_candidates(path: Path, registry: dict[str, Any], *, overwrite: bool) -> None:
    import csv
    import io

    fields = [
        "review_group",
        "case_directories",
        "asset_count",
        "formats",
        "manifest_schemas",
        "figure_modes",
        "style_profiles",
        "proposed_classifications",
        "assets_with_repository_references",
        "repository_reference_count",
        "proposed_action",
        "manual_questions",
    ]
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for group in registry["review_groups"]:
        row = dict(group)
        for field in fields:
            if isinstance(row.get(field), list):
                row[field] = "; ".join(str(item) for item in row[field])
        writer.writerow({field: row.get(field, "") for field in fields})
    _write_text(path, buffer.getvalue(), overwrite=overwrite)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=PROJECT_ROOT / "data/outputs/figures")
    parser.add_argument("--repo-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--candidates", type=Path, default=DEFAULT_CANDIDATES)
    parser.add_argument(
        "--include-untracked",
        action="store_true",
        help="include untracked assets for a local review only; never performs cleanup",
    )
    parser.add_argument("--generated-at-utc", default=None)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    asset_root = _resolve_path(args.root, repo_root).resolve()
    registry_path = _resolve_path(args.registry, repo_root).resolve()
    candidates_path = _resolve_path(args.candidates, repo_root).resolve()
    if not asset_root.is_dir():
        print(f"[figure-inventory] FAILED: asset root does not exist: {asset_root}")
        return 1

    registry = build_registry(
        asset_root=asset_root,
        repo_root=repo_root,
        tracked_only=not args.include_untracked,
        generated_at_utc=args.generated_at_utc,
    )
    try:
        write_registry(registry_path, registry, overwrite=args.overwrite)
        write_candidates(candidates_path, registry, overwrite=args.overwrite)
    except (OSError, ValueError) as exc:
        print(f"[figure-inventory] FAILED: {exc}")
        return 1

    discovery = registry["discovery"]
    print("[figure-inventory] OK")
    print(f"  assets={discovery['asset_count']}")
    print(f"  case_groups={discovery['case_group_count']}")
    print(f"  review_groups={discovery['review_group_count']}")
    print(f"  visible_untracked_excluded={discovery['visible_untracked_excluded']}")
    print(f"  registry={registry_path}")
    print(f"  candidates={candidates_path}")
    print("  delete_performed=false move_performed=false")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
