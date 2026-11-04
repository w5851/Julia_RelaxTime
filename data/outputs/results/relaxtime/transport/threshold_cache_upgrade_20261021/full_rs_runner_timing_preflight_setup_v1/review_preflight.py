"""Review timing metadata and the actual dry-run artifact on quote-remote."""
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import sys
import tomllib
import zipfile

sys.dont_write_bytecode = True
WORK = Path(r'D:\CodexCompute\Julia_RelaxTime\full-rs-runner-timing-preflight-20261026-v1')
PROJECT = WORK / 'project'


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


package = read(WORK / 'package_manifest.json')
assert digest(Path(package['base_source_zip'])) == package['base_source_archive_sha256']
assert digest(WORK / 'overlay.zip') == package['overlay_archive_sha256']
for relative, checksum in package['source_files'].items():
    assert digest(PROJECT / relative) == checksum, relative
for stage in ('setup', 'unit', 'dryrun'):
    status = read(WORK / f'{stage}.status.json')
    assert status['finished'] and status['exit_code'] == 0 and status['source_verified_before_and_after'], stage
unit = (WORK / 'unit.stdout.log').read_text(encoding='utf-8-sig')
match = re.search(r'Runner timing preserves CLI execution and failure semantics\s*\|\s+(\d+)\s+(\d+)\s', unit)
assert match and int(match[1]) == int(match[2]) and int(match[1]) >= 25, unit[-5000:]
timing = tomllib.loads((WORK / 'runtime_timing.toml').read_text(encoding='utf-8-sig'))
assert timing['entrypoint_returned'] and timing['julia_version'] == '1.12.5'
assert timing['julia_threads'] == timing['blas_threads_at_start'] == timing['blas_threads_after_load'] == timing['blas_threads_at_end'] == 1
assert timing['argv'] == read(WORK / 'dryrun.status.json')['arguments'][5:]
assert timing['entrypoint_sha256'] == digest(PROJECT / 'scripts/relaxtime/run_phase_guided_transport_scan.jl')
assert timing['timing_script_sha256'] == digest(PROJECT / 'benchmark/relaxtime/transport_runner_timing.jl')
assert timing['total_memory_bytes'] > 0 and timing['cpu_threads'] > 0
assert timing['compile_time_included_in_elapsed']
for stage in ('load', 'scan'):
    phase = timing[stage]
    assert phase['compile_timing_available']
    for field in ('elapsed_seconds', 'gc_seconds', 'compile_seconds', 'recompile_seconds'):
        assert math.isfinite(phase[field]) and phase[field] >= 0
    assert phase['allocated_bytes'] >= 0
assert timing['measured_wrapper_seconds'] >= timing['load']['elapsed_seconds'] + timing['scan']['elapsed_seconds']
spec = importlib.util.spec_from_file_location('artifact', PROJECT / 'scripts/relaxtime/check_transport_production_artifact.py')
artifact = importlib.util.module_from_spec(spec)
spec.loader.exec_module(artifact)
checked = artifact.check_artifact(WORK / 'dry_run', PROJECT / 'config/models/pnjl/transport_sampling/thermal_weighted_v1.toml', 'unknown', 'dry-run', PROJECT)
assert checked == dict(dry_run=True, planned_points=2)
review = dict(host='W', timing_contract_passed=True, unit_assertions=int(match[1]),
              actual_dry_run_planned_points=2, source_files_verified=package['files'],
              base_source_archive_sha256=package['base_source_archive_sha256'],
              overlay_archive_sha256=package['overlay_archive_sha256'],
              timing_file_sha256=digest(WORK / 'runtime_timing.toml'),
              reviewer_sha256=digest(Path(__file__)),
              numerical_sources_unchanged=True, runner_cost_qualified=False,
              scope='CLI timing and failure semantics only; W dry-run timings are not GitHub runner costs')
(WORK / 'reviewed_result.json').write_text(json.dumps(review, indent=2) + '\n', encoding='utf-8')
files = [p for p in WORK.iterdir() if p.is_file() and p.suffix in ('.json', '.toml', '.log', '.ps1', '.py')]
files += [WORK / 'overlay.zip']
files += [p for p in (WORK / 'dry_run').rglob('*') if p.is_file()]
archive = WORK / 'returned_evidence.zip'
with zipfile.ZipFile(archive, 'x', zipfile.ZIP_DEFLATED) as output:
    for path in sorted(set(files)):
        output.write(path, path.relative_to(WORK).as_posix())
info = dict(sha256=digest(archive), bytes=archive.stat().st_size, files=len(set(files)))
(WORK / 'returned_archive_manifest.json').write_text(json.dumps(info, indent=2) + '\n', encoding='utf-8')
print(json.dumps(review, indent=2))
print(json.dumps(info, indent=2))
