"""Relocate the accepted PRD case metadata without rendering or losing old bytes."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

root = Path(sys.argv[1]).resolve()
case = "pnjl/phase_diagram/figure4_phase_diagram_prod_v1__prd_v3__pdf_20261008"
figures = (root / "data/outputs/figures" / case).resolve()
metadata = (root / "data/outputs/results" / case).resolve()
assert figures.is_relative_to(root / "data/outputs/figures")
assert metadata.is_relative_to(root / "data/outputs/results")
assert not metadata.exists(), "refusing to overwrite relocated metadata"
assert figures.is_dir()

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

def file_record(path, role):
    raw = path.read_bytes()
    return {"path": path.relative_to(root).as_posix(), "path_kind": "repository_relative",
            "role": role, "bytes": len(raw), "sha256": sha(raw)}

source_commit = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
assert source_commit == "afcb9194970f59483cfc7db7068b187b25beaf33"
original = {p.relative_to(root).as_posix(): p.read_bytes()
            for p in figures.rglob("*") if p.is_file()}
assert len(original) == 24
for relative, raw in original.items():
    assert (root / relative).resolve().is_relative_to(figures)
    assert subprocess.check_output(["git", "-C", str(root), "show", "HEAD:" + relative]) == raw

archive = metadata / "provenance/accepted_delivery_before_relocation.zip"
archive.parent.mkdir(parents=True)
with zipfile.ZipFile(archive, "x", compression=zipfile.ZIP_DEFLATED) as z:
    for relative, raw in sorted(original.items()):
        z.writestr(relative, raw)
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    assert set(z.namelist()) == set(original)
    assert all(z.read(relative) == raw for relative, raw in original.items())

mappings = {}
for relative in original:
    source = root / relative
    if source.name == "plot_manifest.json" or source.suffix.lower() in {".png", ".svg", ".pdf"}:
        continue
    destination = metadata / source.relative_to(figures)
    assert source.resolve().is_relative_to(figures)
    assert destination.resolve().is_relative_to(metadata)
    assert not destination.exists()
    destination.parent.mkdir(parents=True, exist_ok=True)
    source.rename(destination)
    mappings[relative] = destination.relative_to(root).as_posix()
for directory in sorted((p for p in figures.rglob("*") if p.is_dir()),
                        key=lambda p: len(p.parts), reverse=True):
    assert directory.resolve().is_relative_to(figures)
    directory.rmdir()  # Empty directories only; no recursive removal.

def relocate_paths(value):
    if isinstance(value, list):
        for child in value:
            relocate_paths(child)
    elif isinstance(value, dict):
        if value.get("path") in mappings:
            value["path"] = mappings[value["path"]]
        for child in value.values():
            relocate_paths(child)

receipt_path = metadata / "provenance/png_acceptance.json"
receipt = json.loads(receipt_path.read_bytes())
relocate_paths(receipt)
write_json(receipt_path, receipt)
receipt_record = file_record(receipt_path, "author_png_acceptance")

retention_path = metadata / "retention.json"
retention = json.loads(retention_path.read_bytes())
retention.update(
    policy="retain final images and plot_manifest in figures; retain captions, acceptance, source snapshots and archives in the same-named results case",
    original_delivery_archive={"path": "provenance/accepted_delivery_before_relocation.zip",
                               "bytes": archive.stat().st_size, "sha256": sha(archive.read_bytes())},
    figure_directory=figures.relative_to(root).as_posix(),
    metadata_directory=metadata.relative_to(root).as_posix(),
    storage_migration="storage_relocation.json records the byte-preserving mapping; original delivery and accepted PNG graphs can be restored from their archives",
)
write_json(retention_path, retention)
(metadata / "README.md").write_text(
    "# PNJL 相结构图 PRD 修订\n\n"
    "阶段：vector_delivery；作者已接受。\n\n"
    f"图像与唯一 plot_manifest.json 位于 `data/outputs/figures/{case}/`。"
    "本目录保存 caption、插入限制、接受记录、冻结源码及验收证据。"
    "图像字节和显示语义保持不变；15 条显示闭合段与经验密度估计仍单独记录。\n\n"
    "`storage_relocation.json` 记录 2026-10-09 的存储迁移；"
    "`provenance/accepted_delivery_before_relocation.zip` 保存迁移前完整 24 文件包，"
    "`provenance/accepted_png_case.zip` 保存作者接受的完整 PNG 阶段。"
    "两个归档均以仓库相对路径为成员名，可在单独目录恢复原始图谱。\n",
    encoding="utf-8",
)
script_copy = metadata / "provenance/relocate_storage.py"
shutil.copyfile(Path(__file__), script_copy)

manifest_path = figures / "plot_manifest.json"
before_manifest = json.loads(original[manifest_path.relative_to(root).as_posix()])
manifest = copy.deepcopy(before_manifest)
relocate_paths(manifest)
manifest["author_review"]["record"] = receipt_record
for record in manifest["inputs"]:
    if record["role"] == "author_png_acceptance":
        record.update(receipt_record)
manifest["metadata_directory"] = metadata.relative_to(root).as_posix()
write_json(manifest_path, manifest)
sys.path.insert(0, str(root))
from scripts.plotting.plot_delivery import SEMANTIC_FIELDS, RENDER_FIELDS
from scripts.plotting.validate_plot_artifact import validate_manifest
assert all(before_manifest[key] == manifest[key] for key in SEMANTIC_FIELDS)
assert all(before_manifest["rendering"][key] == manifest["rendering"][key] for key in RENDER_FIELDS)
assert before_manifest["outputs"] == manifest["outputs"]
assert validate_manifest(manifest_path, repo_root=root) == []

relocation = {
    "schema": "plot_storage_relocation_v1", "date": "2026-10-09", "timezone": "Asia/Shanghai",
    "reason": "Conform to the figures/results storage contract found by PR #326 CI",
    "source_commit": source_commit, "source_pr": "https://github.com/w5851/Julia_RelaxTime/pull/326",
    "original_delivery_archive": file_record(archive, "original_delivery_before_storage_relocation"),
    "migration_script": file_record(script_copy, "storage_migration_script"),
    "original_manifest": {"archive_member": manifest_path.relative_to(root).as_posix(),
                          "sha256": sha(original[manifest_path.relative_to(root).as_posix()])},
    "moved_files": [{"from": old, "to": new, "before_sha256": sha(original[old]),
                     "after_sha256": sha((root / new).read_bytes())}
                    for old, new in sorted(mappings.items())],
    "unchanged_render_outputs": before_manifest["outputs"],
    "original_archive_files_verified": len(original), "image_bytes_unchanged": True,
    "rendered_again": False, "display_semantics_unchanged": True, "new_author_acceptance": False,
    "author_acceptance_note": "Only artifact locations in the existing receipt changed; the complete original receipt remains archived",
    "artifact_contract_errors": [],
}
relocation_path = metadata / "storage_relocation.json"
write_json(relocation_path, relocation)
manifest["storage_relocation"] = file_record(relocation_path, "storage_relocation_proof")
write_json(manifest_path, manifest)
assert validate_manifest(manifest_path, repo_root=root) == []
for record in before_manifest["outputs"]:
    assert sha((root / record["path"]).read_bytes()) == record["sha256"]
assert {p.name for p in figures.iterdir()} == {
    "phase_diagram_TmuB_Trho.png", "phase_diagram_TmuB_Trho_gray.png",
    "phase_diagram_TmuB_Trho.pdf", "plot_manifest.json"}
print(json.dumps({"metadata_directory": str(metadata), "moved_files": len(mappings),
                  "original_archive_files_verified": len(original), "validation": "passed"}, ensure_ascii=False))
