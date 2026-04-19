from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any
import hashlib
import json
import sys


SCHEMA_VERSION_V1 = "v1"

REQUIRED_V1_FIELDS = (
    "schema_version",
    "generated_at_utc",
    "image_path",
    "image_sha256",
    "script_path",
    "command",
    "git_commit",
    "julia_version",
    "input_data_hashes",
)


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def compute_sha256(path: Path | str) -> str:
    p = Path(path)
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _normalize_input_hashes(input_data_hashes: list[dict[str, Any]] | None, image_path: Path) -> list[dict[str, str]]:
    if input_data_hashes:
        out: list[dict[str, str]] = []
        for entry in input_data_hashes:
            if not isinstance(entry, dict):
                raise ValueError("input_data_hashes entries must be dict")
            path = str(entry.get("path", "")).strip()
            sha = str(entry.get("sha256", "")).strip()
            if not path or not sha:
                raise ValueError("input_data_hashes entry must contain non-empty path and sha256")
            out.append({"path": path, "sha256": sha})
        return out

    # keep non-empty default for minimal contract tests
    return [{"path": str(image_path), "sha256": compute_sha256(image_path)}]


def build_image_provenance(
    *,
    image_path: Path | str,
    script_path: str,
    command: str,
    git_commit: str = "unknown",
    julia_version: str | None = None,
    input_data_hashes: list[dict[str, Any]] | None = None,
    generated_at_utc: str | None = None,
    config_path: str | None = None,
    config_hash: str | None = None,
    seed: int | None = None,
    artifact_paths: list[str] | None = None,
    notes: str | None = None,
) -> dict[str, Any]:
    img = Path(image_path)
    payload: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION_V1,
        "generated_at_utc": generated_at_utc or _utc_now_iso(),
        "image_path": str(img),
        "image_sha256": compute_sha256(img),
        "script_path": script_path,
        "command": command,
        "git_commit": git_commit,
        "julia_version": julia_version or f"python-{sys.version.split()[0]}",
        "input_data_hashes": _normalize_input_hashes(input_data_hashes, img),
    }

    if config_path is not None:
        payload["config_path"] = config_path
    if config_hash is not None:
        payload["config_hash"] = config_hash
    if seed is not None:
        payload["seed"] = int(seed)
    if artifact_paths:
        payload["artifact_paths"] = [str(p) for p in artifact_paths]
    if notes is not None:
        payload["notes"] = notes

    return payload


def write_image_provenance_sidecar(
    *,
    image_path: Path | str,
    sidecar_path: Path | str | None = None,
    script_path: str,
    command: str,
    git_commit: str = "unknown",
    julia_version: str | None = None,
    input_data_hashes: list[dict[str, Any]] | None = None,
    generated_at_utc: str | None = None,
    config_path: str | None = None,
    config_hash: str | None = None,
    seed: int | None = None,
    artifact_paths: list[str] | None = None,
    notes: str | None = None,
) -> Path:
    img = Path(image_path)
    sidecar = Path(sidecar_path) if sidecar_path is not None else img.with_name(f"{img.name}.provenance.json")
    payload = build_image_provenance(
        image_path=img,
        script_path=script_path,
        command=command,
        git_commit=git_commit,
        julia_version=julia_version,
        input_data_hashes=input_data_hashes,
        generated_at_utc=generated_at_utc,
        config_path=config_path,
        config_hash=config_hash,
        seed=seed,
        artifact_paths=artifact_paths,
        notes=notes,
    )
    sidecar.parent.mkdir(parents=True, exist_ok=True)
    sidecar.write_text(json.dumps(payload, ensure_ascii=False, sort_keys=True, indent=2), encoding="utf-8")
    return sidecar


def verify_image_provenance(
    *, image_path: Path | str, sidecar_path: Path | str | None = None
) -> tuple[bool, str]:
    img = Path(image_path)
    sidecar = Path(sidecar_path) if sidecar_path is not None else img.with_name(f"{img.name}.provenance.json")
    if not sidecar.exists():
        return False, "sidecar not found"

    try:
        payload = json.loads(sidecar.read_text(encoding="utf-8"))
    except Exception as e:
        return False, f"invalid sidecar json: {e}"

    for key in REQUIRED_V1_FIELDS:
        if key not in payload:
            return False, f"missing required field: {key}"

    input_hashes = payload.get("input_data_hashes")
    if not isinstance(input_hashes, list):
        return False, "input_data_hashes must be list"
    for entry in input_hashes:
        if not isinstance(entry, dict) or "path" not in entry or "sha256" not in entry:
            return False, "invalid input_data_hashes entry"

    actual = compute_sha256(img)
    expected = str(payload.get("image_sha256", ""))
    if actual != expected:
        return False, "image sha256 mismatch"
    return True, "ok"


def verify_image_provenance_sidecar(*, image_path: Path | str, sidecar_path: Path | str | None = None) -> bool:
    ok, _ = verify_image_provenance(image_path=image_path, sidecar_path=sidecar_path)
    return ok
