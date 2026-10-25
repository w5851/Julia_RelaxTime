"""Independently review mode-B qualification and prepare provisional runner inputs."""
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import zipfile

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
WORK = Path(r'D:\CodexCompute\Julia_RelaxTime\full-rs-mode-b-qualification-20261025-v1')
PROJECT = WORK / 'project'
BACKGROUND = Path(r'D:\CodexCompute\Julia_RelaxTime\full-rs-mode-b-background-20261025-v1')


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(line for line in stream if line.strip() and not line.startswith('#')))


def helper(relative, name):
    spec = importlib.util.spec_from_file_location(name, PROJECT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def same(a, b):
    return math.isclose(float(a), float(b), rel_tol=1e-12, abs_tol=1e-14)


package = read(WORK / 'package_manifest.json')
assert digest(WORK / 'source.zip') == package['zip_sha256']
for relative, checksum in package['source_files'].items():
    assert digest(PROJECT / relative) == checksum, relative
selection = read(WORK / 'inputs/selection.json')
input_manifest = read(WORK / 'input_manifest.json')
for name, checksum in input_manifest['files'].items():
    assert digest(WORK / 'inputs' / name) == checksum
assert digest(BACKGROUND / 'reviewed_result.json') == input_manifest['background_review_sha256']
assert read(BACKGROUND / 'reviewed_result.json')['mode_b_background_qualified']
for filename, field in (('manifest.json', 'background_manifest_sha256'), ('status.json', 'background_status_sha256'),
                        ('sampling_plan.csv', 'sampling_plan_sha256'), ('background.csv', 'background_csv_sha256')):
    assert digest(BACKGROUND/'background'/filename) == selection[field]
for path, checksum in selection['trace_files'].items():
    assert digest(Path(path)) == checksum
if selection['boundary_sha256'] is not None:
    assert digest(BACKGROUND/'background/phase_boundaries.csv') == selection['boundary_sha256']
points = set(selection['points'])
assert len(points) == len(selection['points']) and len(points) >= 15
snapshots = 0
for stage in ('setup', 'selection', 'unit', 'cheap', 'kernel', 'outer', 'joint', 'production'):
    status = read(WORK / f'{stage}.status.json')
    assert status['finished'] and status['exit_code'] == 0, stage
    assert status['runs'] and all(r['process_exit_code'] == 0 and r['source_verified_before_and_after'] for r in status['runs'])
    if stage in ('setup', 'selection', 'unit'):
        continue
    run = WORK / stage
    manifest = read(run / 'manifest.json')
    assert manifest['head'] is None and manifest['source_archive']['sha256'] == package['zip_sha256']
    assert digest(run / 'source_package_manifest.json') == manifest['source_archive']['package_manifest_sha256']
    for item in manifest['source_files']:
        relative = item['path'].replace('\\', '/')
        assert digest(run / 'source' / relative) == item['sha256'] == package['source_files'][relative]
        snapshots += 1
    for item in manifest['input_files']:
        assert digest(run / item['snapshot'].replace('\\', '/')) == item['sha256']
    assert digest(run / 'inputs/mode_b_selection.json') == digest(WORK / 'inputs/selection.json')
    status = read(run / 'status.json')
    assert status['complete'] and status['failures'] == 0
    if stage in ('cheap', 'kernel', 'outer'):
        assert status['all_cases_complete'] and status['layer_target_passed']
assert re.search(r'\|\s+7\s+7\s', (WORK / 'unit.stdout.log').read_text(encoding='utf-8-sig'))
assert 'Ran 6 tests' in (WORK / 'python_unit.stderr.log').read_text(encoding='utf-8-sig')

artifacts = helper('scripts/relaxtime/check_transport_production_artifact.py', 'artifacts')
inputs = helper('scripts/analysis/relaxtime/prepare_transport_production_inputs.py', 'inputs')
density = helper('scripts/analysis/relaxtime/verify_transport_joint_density.py', 'density')
baseline = helper('scripts/dev/export_thermal_weighted_transport_baseline.py', 'baseline')
background_rows = rows(WORK / 'cheap/background.csv')
backgrounds = {(r['point_id'], r['configuration']): r for r in background_rows}
expected_backgrounds = {(p, c) for p in points for c in ('candidate', 'reference')}
assert len(background_rows) == len(backgrounds) == 2*len(points) and set(backgrounds) == expected_backgrounds
joint_backgrounds = {(r['point_id'], r['configuration']): r for r in rows(WORK / 'joint/background.csv')}
assert set(joint_backgrounds) == expected_backgrounds
background_changes = []
for key, state in backgrounds.items():
    for name in ('m_u', 'm_d', 'm_s', 'Phi', 'Phibar', 'A_u', 'A_d', 'A_s'):
        relative = abs(float(state[name])/float(joint_backgrounds[key][name])-1)
        assert math.isfinite(relative) and relative <= 1e-10
        background_changes.append(relative)
identities = rows(WORK / 'cheap/identity.csv')
assert len(identities) == len(expected_backgrounds)
for row in identities:
    assert row['A_exact'] == row['passed'] == 'true'
    assert 0 <= float(row['coupling_max_relative']) <= float(row['coupling_roundoff_limit']) <= 2e-15

density_rows = rows(WORK / 'cheap/densities.csv')
assert len(density_rows) == 24*len(expected_backgrounds)
python_changes = []
for row in density_rows:
    point = row['point_id'], row['background_configuration']
    state = backgrounds[point]
    config = row['configuration']
    if config == 'workflow_p64_c8':
        value = density.density(state, row['species'], 64, 8, 25 if point[1] == 'candidate' else 30, semi=False)
    else:
        p, c = {'thermal_p32_c36': (32, 36), 'thermal_p32_c48': (32, 48), 'identity_p128_c64': (128, 64)}[config]
        value = density.density(state, row['species'], p, c, 10.0)
    ref = density.density(state, row['species'], 256, 1, 5.0, xi=0.0)/math.sqrt(1+float(state['xi']))
    finer = density.density(state, row['species'], 128, 1, 5.0, xi=0.0)/math.sqrt(1+float(state['xi']))
    assert all(math.isfinite(v) and v > 0 for v in (value, ref, finer))
    assert abs(finer/ref-1) <= 1e-9
    python_changes.extend((abs(value/float(row['density'])-1), abs(ref/float(row['reference'])-1)))
    relative = float(row['density'])/float(row['reference'])-1
    target = 1e-9 if config == 'identity_p128_c64' else 1e-3
    assert same(relative, row['relative_change']) and float(row['target']) == target and abs(relative) <= target
    assert row['passed'] == 'true'
assert max(python_changes) <= 1e-10
pairs = rows(WORK / 'cheap/pair_bounds.csv')
assert len(pairs) == 2*len(expected_backgrounds)
for row in pairs:
    selected = [r for r in density_rows if all(r[k] == row[k] for k in ('point_id', 'background_configuration', 'configuration'))]
    assert len(selected) == 6
    relative = [float(r['relative_change']) for r in selected]
    bound = max(abs(1/((1+a)*(1+b))-1) for a in relative for b in relative)
    assert same(bound, row['maximum_pair_rate_error']) and bound <= 1e-3 and row['passed'] == 'true'
moments, a_rows = rows(WORK / 'cheap/moments.csv'), rows(WORK / 'cheap/A.csv')
assert len(moments) == 6*len(expected_backgrounds) and len(a_rows) == 9*len(expected_backgrounds)
for row in moments:
    assert row['passed'] == 'true'
    assert all(0 <= float(row[k]) <= .01 for k in ('eta_any_positive_tau_bound', 'sigma_any_positive_tau_bound', 'bulk_any_bracket_bound'))
for row in a_rows:
    change = float(row['A'])/float(row['reference'])-1
    assert same(change, row['relative_change']) and abs(change) <= .01 and row['passed'] == 'true'

layer_maxima = {}
for layer, nconfig, ncompare in (('kernel', 6, 5), ('outer', 12, 11)):
    cases = {tuple(row) for row in selection[layer+'_cases']}
    rates, compared = rows(WORK/layer/'rates.csv'), rows(WORK/layer/'comparisons.csv')
    assert len(rates) == len(cases)*nconfig and len(compared) == len(cases)*ncompare
    values = {(r['point_id'], r['channel'], r['configuration']): float(r['rate']) for r in rates}
    assert len(values) == len(rates)
    assert {(r['point_id'], r['channel']) for r in rates} == cases
    for row in rates:
        assert math.isfinite(float(row['rate'])) and float(row['rate']) > 0 and 0 <= float(row['reference_l1']) <= 2e-4
    for row in compared:
        pair = row['point_id'], row['channel']
        a, b = values[*pair, row['candidate_configuration']], values[*pair, row['reference_configuration']]
        assert same(a, row['candidate']) and same(b, row['reference'])
        assert same(a/b-1, row['relative_change']) and abs(a/b-1) <= .01 and row['passed'] == 'true'
    layer_maxima[layer] = max(abs(float(r['relative_change'])) for r in compared)

joint = WORK/'joint'
joint_status, joint_manifest = read(joint/'status.json'), read(joint/'manifest.json')
assert joint_status['joint_target_passed'] and joint_status['points'] == len(points)
config = joint_manifest['parameters']['configurations']
assert (config['candidate']['p'], config['candidate']['c'], config['candidate']['phi'], config['candidate']['a_p'], config['candidate']['a_c'], config['candidate']['t']) == (48, 36, 36, 32, 16, 24)
rates = rows(joint/'rates.csv')
expected = {(p, c, channel) for p in points for c in ('candidate', 'reference') for channel in artifacts.RATE_KEYS}
assert len(rates) == len(expected) and {(r['point_id'], r['configuration'], r['channel']) for r in rates} == expected
for row in rates:
    rate, error = float(row['rate']), float(row['estimated_or_reference_l1'])
    assert math.isfinite(rate) and rate > 0 and math.isfinite(error) and error >= 0
    assert error <= (2e-4 if row['configuration'] == 'reference' else 1e-3+1e-12/rate)
compared = rows(joint/'comparisons.csv')
expected = {(p, obs) for p in points for obs in artifacts.OBSERVABLES}
targets = {(r['point_id'], r['observable']): r for r in compared}
assert len(compared) == len(targets) and set(targets) == expected
for row in compared:
    a, b = float(row['candidate']), float(row['reference'])
    assert a > 0 and b > 0 and math.isfinite(a) and math.isfinite(b)
    assert same(a/b-1, row['relative_change']) and abs(a/b-1) <= .01 and row['passed'] == 'true'
for row in rows(joint/'observables.csv'):
    assert 0 <= float(row['density_relative_error']) <= 1e-3

production = WORK/'production'
prod_manifest, prod_status = read(production/'manifest.json'), read(production/'status.json')
source_hashes = inputs.verified_sources(production, prod_manifest, PROJECT.resolve())
assert source_hashes == selection['numerical_sources']
commit, identity = baseline.replay_artifact_commit(production, prod_manifest)
assert commit == 'unknown' and identity['sha256'] == package['zip_sha256']
assert prod_status['production_replay_passed'] and prod_status['completed_points'] == len(points)
production_rows = rows(production/'comparisons.csv')
assert len(production_rows) == len(expected) and {(r['point_id'], r['observable']) for r in production_rows} == expected
input_rows = {r['point_id']: r for r in rows(WORK/'inputs/review_points.csv')}
actual = {}
profile = PROJECT/'config/models/pnjl/transport_sampling/thermal_weighted_v1.toml'
for point in points:
    folder = production/'points'/point
    checked = artifacts.check_artifact(folder, profile, commit, 'production', PROJECT)
    assert checked['completed_points'] == 1 and checked['qualified_channels'] == 17
    result = rows(folder/'phase_guided_transport_scan.csv')
    assert len(result) == 1 and inputs.point_key(result[0]) == inputs.point_key(input_rows[point])
    actual[point] = result[0]
for row in production_rows:
    point, observable = row['point_id'], row['observable']
    target = targets[point, observable]
    assert float(row['qualified_candidate']) == float(target['candidate'])
    assert float(row['qualified_reference']) == float(target['reference'])
    assert float(row['production']) == float(actual[point][observable])
    inputs.verify_comparison(row)

plan = rows(BACKGROUND/'background/sampling_plan.csv')
assert len(plan) == 909
groups = {}
for row in plan:
    groups.setdefault((float(row['T_MeV']), float(row['muB_MeV'])), []).append(row)
assert len(groups) == 9 and all(len(group) == 101 for group in groups.values())


def jobs(selected, stage):
    selected_groups = {}
    for row in selected:
        selected_groups.setdefault((float(row['T_MeV']), float(row['muB_MeV'])), []).append(row)
    result = []
    for (t, mu), group in sorted(selected_groups.items()):
        group.sort(key=lambda r: float(r['xi']))
        for start in range(0, len(group), 20):
            chunk = group[start:start+20]
            label = f'{stage}_modeB_T{t:g}_muB{mu:g}_{start//20:02d}'
            result.append(dict(label=label, expected_points=len(chunk), expected_T_MeV=t,
                               inputs=dict(case_name='thermal_weighted_full_rs_mode_b_20261025', shard_label=label,
                                           sampling_profile='thermal_weighted_v1', mode='fixed-T-sparse-muB',
                                           t_list=f'{t:g}', muB_list=f'{mu:g}', alpha_t_list='1.0',
                                           xi_list=','.join(str(r['xi']) for r in chunk), run_tier='production', render_plots=False)))
    return result


pilot_ids = {point for point, _ in selection['kernel_cases']}
pilot_keys = {inputs.point_key(input_rows[point]) for point in pilot_ids}
pilot = [row for row in plan if inputs.point_key(row) in pilot_keys]
assert len(pilot) == len(pilot_keys)
payload = dict(schema_version=1, mode='B', case_name='thermal_weighted_full_rs_mode_b_20261025',
               sampling_profile='thermal_weighted_v1', sampling_profile_sha256=digest(profile), calculation_sha=None,
               numerical_source_sha256=hashlib.sha256(json.dumps(source_hashes, sort_keys=True).encode()).hexdigest(),
               source_archive_sha256=package['zip_sha256'], phase_source=str((BACKGROUND/'background').resolve()),
               sampling_plan_sha256=digest(BACKGROUND/'background/sampling_plan.csv'),
               representative_points=len(points), qualified_observables=len(production_rows),
               pilot_jobs=jobs(pilot, 'pilot'), provisional_full_jobs=jobs(plan, 'full'), full_points=909,
               pilot_input_ready=False, full_dispatch_ready=False, runner_cost_qualified=False,
               remaining_before_dispatch=['select the commit containing the reviewed implementation'],
               full_shard_policy='Provisional 20-point shards; calibrate actual runner startup, point times and peak memory before full dispatch')
assert sum(job['expected_points'] for job in payload['provisional_full_jobs']) == 909
runner_dir = WORK/'runner_inputs'
runner_dir.mkdir()
(runner_dir/'runner_inputs.json').write_text(json.dumps(payload, indent=2)+'\n', encoding='utf-8')

review = dict(host='W', mode='B', mode_b_background_qualified=True, background_grid_points=909,
              mode_b_transport_representatives_qualified=True, representative_points=len(points),
              observable_comparisons=len(production_rows), qualified_production_channels=17*len(points),
              joint_channel_rates=len(rates), independent_reference_rates=17*len(points), layer_targets_passed=True,
              maximum_numpy_julia_density_relative=max(python_changes),
              maximum_background_replay_relative=max(background_changes), maximum_layer_relative_changes=layer_maxima,
              maximum_joint_relative_change=max(abs(float(r['relative_change'])) for r in compared),
              maximum_production_replay_relative_change=max(abs(float(r['replay_relative_change'])) for r in production_rows),
              maximum_production_reference_relative_change=max(abs(float(r['reference_relative_change'])) for r in production_rows),
              source_archive_sha256=package['zip_sha256'], source_files_verified=package['files'], source_snapshots_verified=snapshots,
              python_tests_passed=6, selection_assertions_passed=7,
              runner_inputs_sha256=digest(runner_dir/'runner_inputs.json'), reviewer_sha256=digest(Path(__file__)),
              calculation_commit=None, whole_grid_transport_calculated=False, runner_cost_qualified=False, formal_production_promoted=False)
(WORK/'reviewed_result.json').write_text(json.dumps(review, indent=2)+'\n', encoding='utf-8')
archive = WORK/'returned_evidence.zip'
assert not archive.exists()
files = []
for stage in ('inputs', 'cheap', 'kernel', 'outer', 'joint', 'joint_summary', 'production', 'runner_inputs'):
    files.extend(p for p in (WORK/stage).rglob('*') if p.is_file())
files.extend(p for p in WORK.iterdir() if p.is_file() and p.suffix in ('.json', '.log', '.ps1', '.py'))
with zipfile.ZipFile(archive, 'x', zipfile.ZIP_DEFLATED, compresslevel=6) as output:
    for path in sorted(set(files)):
        output.write(path, path.relative_to(WORK).as_posix())
info = dict(sha256=digest(archive), bytes=archive.stat().st_size, files=len(set(files)))
(WORK/'returned_archive_manifest.json').write_text(json.dumps(info, indent=2)+'\n', encoding='utf-8')
print(json.dumps(review, indent=2))
print(json.dumps(info, indent=2))
