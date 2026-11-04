"""Freeze mode-B qualification orchestration without changing production numerics."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import zipfile

ROOT = Path(r'D:\w\meson-freezeout-defaults')
OUT = Path(__file__).resolve().parent
PREVIOUS = OUT.parent / 'full_rs_baseline_acceptance_setup_v1'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


previous = json.loads((PREVIOUS / 'package_manifest.json').read_text(encoding='utf-8-sig'))
extras = {'scripts/analysis/relaxtime/TransportValidationSelection.jl',
          'scripts/analysis/relaxtime/prepare_transport_mode_b_qualification.py',
          'scripts/analysis/relaxtime/run_transport_mode_b_qualification.jl',
          'tests/unit/relaxtime/test_transport_validation_selection.jl',
          'tests/unit/python/test_transport_mode_b_qualification.py'}
paths = sorted(set(previous['source_files']) | extras)
source_files = {name: digest(ROOT / name) for name in paths}
changed = [name for name in paths if source_files[name] != previous['source_files'].get(name)]
assert not any(name.startswith(('src/', 'config/models/', 'config/physics/', 'scripts/relaxtime/'))
               or name in ('Project.toml', 'Manifest.toml') for name in changed), changed
assert not (OUT / 'source.zip').exists() and not (OUT / 'package_manifest.json').exists()
with zipfile.ZipFile(OUT / 'source.zip', 'x', zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
    for name in paths:
        archive.write(ROOT / name, name)
manifest = dict(purpose='mode-B background-derived selection, independent layers, joint convergence and actual entrypoint checks',
                created_utc=datetime.now(timezone.utc).isoformat(), source_root=str(ROOT),
                source_branch=previous['source_branch'], head_base=previous['head_base'],
                calculation_commit=None, uncommitted_source=True,
                remote_workspace=r'D:\CodexCompute\Julia_RelaxTime\full-rs-mode-b-qualification-20261025-v1',
                runtime_version='1.12.5', source_files=source_files, files=len(source_files),
                zip_sha256=digest(OUT / 'source.zip'), zip_bytes=(OUT / 'source.zip').stat().st_size,
                previous_source_archive_sha256=previous['zip_sha256'], changed_source_files=changed,
                production_numerical_sources_unchanged=True, packager_sha256=digest(Path(__file__)))
(OUT / 'package_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
print(json.dumps({key: manifest[key] for key in ('files', 'zip_bytes', 'zip_sha256', 'changed_source_files')}, indent=2))
