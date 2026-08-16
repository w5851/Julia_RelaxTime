"""Build traceable plot manifests without touching numerical inputs."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
from typing import Any, Iterable, Mapping


PROJECT_ROOT = Path(__file__).resolve().parents[2]
MANIFEST_SCHEMA = "plot_manifest_v1"


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _repo_relative(path: Path, root: Path) -> str | None:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return None


def _display_path(path: Path, root: Path) -> tuple[str, str]:
    relative = _repo_relative(path, root)
    if relative is not None:
        return relative, "repository_relative"
    return str(path.resolve()), "external_absolute"


def input_record(
    path: str | Path,
    *,
    role: str,
    root: Path = PROJECT_ROOT,
    schema: str | None = None,
    units: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Describe a read-only CSV/JSON input, including its current hash."""

    resolved = Path(path).resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"plot input not found: {resolved}")
    display, path_kind = _display_path(resolved, root)
    record: dict[str, Any] = {
        "path": display,
        "path_kind": path_kind,
        "role": role,
        "bytes": resolved.stat().st_size,
        "sha256": sha256_file(resolved),
    }
    if schema is not None:
        record["schema"] = schema
    if units is not None:
        record["units"] = dict(units)
    return record


def generator_record(
    script_path: str | Path,
    *,
    command: str,
    root: Path = PROJECT_ROOT,
    runtime: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Describe the exact plotting generator used for a case."""

    resolved = Path(script_path).resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"plot generator not found: {resolved}")
    display, path_kind = _display_path(resolved, root)
    result: dict[str, Any] = {
        "path": display,
        "path_kind": path_kind,
        "sha256": sha256_file(resolved),
        "command": command,
        "runtime": dict(runtime or {}),
    }
    return result


def output_record(
    path: str | Path,
    *,
    fmt: str,
    dpi: int | None,
    vector: bool,
    root: Path = PROJECT_ROOT,
) -> dict[str, Any]:
    """Describe one generated output after it has been written."""

    resolved = Path(path).resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"plot output not found: {resolved}")
    display, path_kind = _display_path(resolved, root)
    return {
        "path": display,
        "path_kind": path_kind,
        "format": str(fmt).lower(),
        "dpi": dpi,
        "vector": bool(vector),
        "bytes": resolved.stat().st_size,
        "sha256": sha256_file(resolved),
    }


def git_commit(root: Path = PROJECT_ROOT) -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, text=True, encoding="utf-8"
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def runtime_record(extra: Mapping[str, Any] | None = None) -> dict[str, Any]:
    record: dict[str, Any] = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "executable": sys.executable,
    }
    if extra:
        record.update(dict(extra))
    return record


def build_manifest(
    *,
    asset_id: str,
    figure_family: str,
    case_slug: str,
    figure_mode: str,
    semantic_status: str,
    style_profile: str,
    publication_scope: str,
    generator: Mapping[str, Any],
    inputs: Iterable[Mapping[str, Any]],
    axes: Iterable[Mapping[str, Any]],
    series: Iterable[Mapping[str, Any]],
    outputs: Iterable[Mapping[str, Any]],
    selection_rule: str,
    interpolation_policy: str,
    connector_policy: str,
    missing_value_policy: str,
    validation: Mapping[str, Any],
    rendering: Mapping[str, Any] | None = None,
    calculation_sha: str | None = None,
    postprocess_sha: str | None = None,
    source_run_id: str | None = None,
    root: Path = PROJECT_ROOT,
) -> dict[str, Any]:
    """Create the stable top-level manifest shape used by new plot cases."""

    return {
        "schema_version": MANIFEST_SCHEMA,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "asset_id": asset_id,
        "figure_family": figure_family,
        "case_slug": case_slug,
        "asset_kind": "figure",
        "figure_mode": figure_mode,
        "semantic_status": semantic_status,
        "style_profile": style_profile,
        "publication_scope": publication_scope,
        "generator": dict(generator),
        "git_commit": git_commit(root),
        "calculation_sha": calculation_sha,
        "postprocess_sha": postprocess_sha,
        "source_run_id": source_run_id,
        "inputs": [dict(item) for item in inputs],
        "axes": [dict(item) for item in axes],
        "series": [dict(item) for item in series],
        "selection_rule": selection_rule,
        "interpolation_policy": interpolation_policy,
        "connector_policy": connector_policy,
        "missing_value_policy": missing_value_policy,
        "outputs": [dict(item) for item in outputs],
        "rendering": dict(rendering or {}),
        "validation": dict(validation),
        "external_layout": None,
    }


def write_manifest(path: str | Path, manifest: Mapping[str, Any], *, overwrite: bool = False) -> None:
    """Write a manifest atomically and refuse accidental historical overwrite."""

    target = Path(path).resolve()
    if target.exists() and not overwrite:
        raise FileExistsError(f"refusing to overwrite existing manifest: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.name + ".tmp")
    temporary.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    temporary.replace(target)
