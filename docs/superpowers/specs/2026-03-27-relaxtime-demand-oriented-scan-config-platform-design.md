# Relaxtime Demand-Oriented Scan Config Platform Design

## Context

Current `scripts/relaxtime/run_*.jl` scripts already cover most computation paths, but behavior contracts and parameter sources are still fragmented across scripts.

The user request is to freeze behavior for demand-oriented scanning workflows and decide whether optional parameters should be externalized to TOML for easier tuning and reproducibility.

## Confirmed User Requirements (Anti-Drift Contract)

The following three requirements are explicitly frozen in this design and become acceptance targets.

1. Transport-vs-temperature scans at `muB=0`
   - Scan multiple `xi` values.
   - Frozen default `xi_list` (schema v1): `[-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3]`.
   - Optional alternate list such as `[-0.3, -0.15, 0.0, 0.15, 0.3]` is allowed only via TOML/CLI override.
   - Produce transport observables vs `T`: three base + two ratio metrics.
   - Produce quark/antiquark relaxation times vs `T`.

2. Cross-section-vs-`sqrt(s)` scans at `muB=0`
   - Scan multiple `xi` values with same range policy as requirement 1.
   - Include temperature dimension.
   - Frozen default `T_list_MeV` (schema v1): `[150.0, 250.0]`.
   - Alternate temperature lists are allowed only via TOML/CLI override.
   - Default processes are fixed in TOML as:
     - `ud_to_ud`
     - `us_to_us`
     - `udbar_to_udbar`
     - `usbar_to_usbar`
   - CLI process override must remain supported.

3. Plot organization contract
   - For requirement 1: one figure per physical quantity, with multiple `xi` curves in the same figure.
   - For requirement 2: split by temperature and process (different figures), with multiple `xi` curves in each figure.

## Feasibility From Existing Scripts

### Already achievable directly

- Requirement 1 is directly supported by `scripts/relaxtime/run_gap_transport_scan.jl` + `scripts/plot_scan_csv.py`.
- Output columns already include the required metrics:
  - `eta_over_s`
  - `zeta_over_s`
  - `sigma_over_T`
  - `sigma_over_T_over_eta_over_s`
  - `zeta_over_s_over_eta_over_s`
  - `tau_u,tau_d,tau_s,tau_ubar,tau_dbar,tau_sbar`

### Partially achievable, needs orchestration

- Requirement 2 calculation is supported by `scripts/relaxtime/scan_cross_section_vs_s_by_process.jl`.
- But it currently runs one `(T, xi)` at a time, so multi-`T` x multi-`xi` and grouped plotting need orchestration.

## Decision

Adopt a platformized configuration architecture with schema governance.

Core policy is fixed as:

- Parameter precedence: `built-in defaults < TOML < CLI`
- Effective configuration snapshot must be persisted in output metadata.
- Governance strength: **strong governance** (single authoritative workflow path + static gate + consistency regression).

## Architecture

### Layering

1. Core execution layer
   - `scripts/relaxtime/run_gap_transport_scan.jl`
   - single source of truth for gap -> tau -> transport scan chain

2. Orchestration layer
   - `scripts/relaxtime/run_manual_relaxation_scan_workflow.jl`
   - `scripts/relaxtime/run_outputs_2026_02_05.jl` (legacy scenario runner)
   - future cross-section multi-`T`/multi-`xi` runner

3. Patch layer
   - `scripts/relaxtime/run_offline_transport_patch.jl`

4. Regression gate layer
   - `scripts/relaxtime/run_transport_fixedpoint_regression.jl`
   - `scripts/relaxtime/run_total_cross_section_fixedpoint_regression.jl`

## Drift Prevention Governance (Strong)

To prevent script-level process drift, this spec freezes a 4-part governance policy.

1. Script classification (mandatory tags)
   - `public-authoritative`: the only user-facing official entrypoint category.
   - `internal-authoritative`: internal execution nodes callable by orchestrator, not official entrypoints.
   - `regression-gate`: verification-only scripts used by CI/regression gates.
   - `analysis-experimental`: research exploration only, not official output baselines.

2. Single official entrypoint
   - Official scan/report generation must go through `public-authoritative` orchestrator entrypoints only.
   - `internal-authoritative` scripts may still exist for phased migration but are not official user entrypoints.
   - Direct `src/relaxtime/*` script paths are kept for analysis only.

3. Static governance gate
   - Add a governance checker script that fails if any `public-authoritative` script directly includes `src/relaxtime/*`.
   - `internal-authoritative` direct includes are temporarily allowed only through explicit allowlist during migration.
   - Allowlist is explicit and versioned.

4. Workflow-vs-direct consistency regression
   - For retained direct analysis scripts, add a small fixed-point comparison suite against authoritative workflow outputs.
   - Breach beyond configured tolerance fails the gate.

## Parameter Consumption Audit (No-Intercept Contract)

Current codebase may still have "config parsed but not truly consumed" risk.
This spec adds a mandatory no-intercept audit contract.

1. Strict config mode
   - Unknown keys fail fast.
   - Unconsumed keys fail fast.
   - Enabled by default in CI and official profiles.

2. Consumed-key audit trail
   - Each layer returns `consumed_keys` and `overridden_keys`.
   - Run artifacts must include:
     - `effective_config.json`
     - `consumption_report.json`
   - `consumption_report.json` must include:
     - `consumed_keys`
     - `unused_keys`
     - `overridden_keys`
     - `fallback_used`

3. Observable fallback policy
   - Any runtime fallback (example: `compute_bulk=true` but bulk path fails and degrades) must be logged in `fallback_events`.
   - Add switch `--fail-on-fallback` to convert fallback into hard failure when required.

4. End-to-end pass-through tests
   - Add probe-value tests for key parameters to prove deep consumption, including at least:
     - `tr_p_nodes`
     - `tau_n_sigma`
     - cross-section `processes`
   - Acceptance requires changed probe values to be reflected in downstream effective execution evidence.

5. Single normalized config object
   - Entry layer must normalize TOML/CLI into one schema-validated config object before dispatch.
   - Avoid multi-hop ad-hoc key propagation with silent drops.

### Initial classification in this repository

- Scope note:
  - Classification below covers all scripts in scope of this spec (requirements 1-3 and strong-governance anti-drift paths).

- `public-authoritative`:
  - `scripts/relaxtime/run_relaxtime_orchestrator.jl` (to be introduced in this plan)
- `internal-authoritative`:
  - `scripts/relaxtime/run_gap_transport_scan.jl`
  - `scripts/relaxtime/scan_relaxation_times_vs_T.jl`
  - `scripts/relaxtime/run_manual_relaxation_scan_workflow.jl` (legacy orchestration bridge during migration)
  - `scripts/relaxtime/run_outputs_2026_02_05.jl` (legacy scenario runner during migration)
  - `scripts/relaxtime/run_offline_transport_patch.jl`
- `regression-gate`:
  - `scripts/relaxtime/run_transport_fixedpoint_regression.jl`
  - `scripts/relaxtime/run_total_cross_section_fixedpoint_regression.jl`
- `analysis-experimental`:
  - `scripts/relaxtime/scan_cross_section_vs_s_by_process.jl`
  - `scripts/relaxtime/scan_total_cross_section.jl`
  - `scripts/relaxtime/scan_scattering_amplitude.jl`
  - `scripts/relaxtime/scan_scattering_amplitude_st_grid.jl`
  - `scripts/relaxtime/scan_scattering_amplitude_vs_s_mid_t.jl`

## Configuration Platform (Schema v1)

### New configuration tree

- `config/workflows/relaxtime/default.toml`
- `config/workflows/relaxtime/profiles/*.toml`
- `config/workflows/relaxtime/schema/relaxtime_workflow_schema_v1.toml`
- `config/workflows/relaxtime/schema/aliases_v1.toml`

### Canonical entrypoints (for acceptance execution)

- Transport orchestration command:
  - `julia --project=. scripts/relaxtime/run_relaxtime_orchestrator.jl transport --config <path> --outdir <path>`
- Cross-section orchestration command:
  - `julia --project=. scripts/relaxtime/run_relaxtime_orchestrator.jl cross-section --config <path> --outdir <path>`
- Optional CLI options for both subcommands:
  - `--config` (default `config/workflows/relaxtime/default.toml`)
  - `--outdir` (default `data/outputs/results/relaxtime/orchestrated/default`)
  - `--resume`
  - `--overwrite`
  - `--fail-on-fallback`
  - if omitted, values are taken from effective config defaults
- Additional optional CLI option for cross-section (override only):
  - `--processes`
  - if omitted, default process set is read from effective TOML config

### Suggested section groups

- `[scan.transport]`
- `[scan.transport.solver]`
- `[scan.transport.tau]`
- `[scan.transport.transport]`
- `[scan.cross_section]`
- `[scan.cross_section.energy]`
- `[plot.transport]`
- `[plot.cross_section]`

### Cross-section energy-grid contract (schema v1)

- `scan.cross_section.energy.mode`: `list | linspace`
- if `mode=list`: use `scan.cross_section.energy.sqrt_s_list_MeV`
- if `mode=linspace`: use `scan.cross_section.energy.sqrt_s_min_MeV`, `sqrt_s_max_MeV`, `sqrt_s_num`
- normalization rules:
  - remove non-finite points
  - sort ascending
  - deduplicate with exact numeric equality after normalization
- normalized energy grid participates in both `config_hash` and cross-section point-key generation.

### Process naming and alias normalization

Canonical process names are underscore style (e.g. `ud_to_ud`).

Alias compatibility is supported in v1, including:

- `udtoud -> ud_to_ud`
- `ustous -> us_to_us`
- `udbartoudbar -> udbar_to_udbar`
- `usbartousbar -> usbar_to_usbar`

After normalization, metadata persists canonical `effective.processes`.

## Behavior Contracts

0. Config processing pipeline (fixed order)
   - `raw(TOML+CLI) -> alias normalization -> source-level schema validation -> precedence merge -> effective schema validation -> consumption audit -> dispatch`
   - source-level validation checks parseability/type/shape per source.
   - effective validation checks required keys, ranges, and cross-field constraints on merged config.
   - strict-mode unknown-key check runs after alias normalization.
   - strict-mode unconsumed-key check runs after consumption audit.

1. Precedence
   - No conflict error by default; higher-precedence source wins.

2. Traceability metadata
   - Required persisted files in run directory:
     - `run_manifest.json`
     - `effective_config.json`
   - Required fields in `run_manifest.json`:
     - `schema_version`
     - `config_path`
     - `config_hash`
     - `run_id`
   - CSV outputs must include `run_id` column linking each row to `run_manifest.json`.
   - `config_hash` definition:
     - SHA256 over normalized `effective_config.json` (UTF-8, sorted keys, canonical JSON serialization).
   - `timestamp_utc` format is fixed to RFC3339 UTC (example: `2026-03-27T12:34:56Z`).

3. Failure semantics
   - Scan workflows: point-level fault tolerance, continue run, record status fields.
   - Regression workflows: any threshold breach fails the whole run with non-zero exit.
   - Required row-level status fields for scan outputs:
     - `status` in `{ok,error,skipped}`
     - `error_code`
     - `error_message`
     - `timestamp_utc`

4. File semantics
   - `resume=true && overwrite=false`: append/skip-compatible resume.
   - `overwrite=true`: reset output before run.
   - both true: `overwrite` wins, emit warning.
   - Unique point-key contracts:
     - transport scan key: `(T_MeV, muB_MeV, xi)`
     - cross-section scan key: `(T_MeV, muB_MeV, xi, process, sqrt_s_MeV)`
   - `resume=true` must skip existing `status=ok` keys and must not produce duplicate keys.

5. Compatibility semantics
   - Deprecated keys and aliases remain accepted in v1 with warning.
   - Planned removal can happen only in v2 with explicit migration note.
   - `aliases_v1.toml` must list all deprecated key names and process aliases.
   - Deprecation warning format is fixed:
     - `DEPRECATION[v1]: <old_key_or_alias> -> <canonical_key>`

## Acceptance Criteria

### AC-1 (Requirement 1)

- One command can generate a single transport CSV for multi-`xi` at `muB=0`.
- Required columns exist for transport metrics and six `tau_*` fields.
- AC must verify all produced rows satisfy `muB_MeV == 0.0`.
- Plot output provides one figure per physical quantity with `xi` as line grouping:
  - 5 transport metric figures:
    - `eta_over_s`
    - `zeta_over_s`
    - `sigma_over_T`
    - `sigma_over_T_over_eta_over_s`
    - `zeta_over_s_over_eta_over_s`
  - 6 relaxation-time figures:
    - `tau_u`
    - `tau_d`
    - `tau_s`
    - `tau_ubar`
    - `tau_dbar`
    - `tau_sbar`
- AC must verify metadata has frozen default `xi_list` when not overridden.
- AC expected figure count for requirement 1 is fixed to `11` under schema v1 defaults.

### AC-2 (Requirement 2)

- One orchestration command can run `T={150,250}` x `xi_list` x process-set.
- Default process set is the 4-process TOML set above.
- CLI process override works and has higher precedence than TOML.
- AC must verify all produced rows satisfy `muB_MeV == 0.0`.
- Plot output is split by `(T, process)`, each figure grouped by `xi`.
- Default `T_list_MeV=[150.0,250.0]` must appear in effective config metadata when not overridden.

### AC-3 (Requirement 3)

- Directory and filename layout are deterministic by profile + `T` + process.
- Re-running same profile under `resume/overwrite` contracts remains reproducible.
- Required filename pattern:
  - transport figures: `transport__{metric}__muB{muB_tag}.png`
  - cross-section figures: `xsec__T{T_tag}__{process}.png`
- Required numeric tag formatting:
  - integer-like values use integer string (example: `150`)
  - otherwise replace decimal dot with `p` (example: `150.5 -> 150p5`)
- Required figure-count checks:
  - transport: `11`
  - cross-section: `N_T * N_process`
  - for schema v1 default (`T_list=[150,250]`, 4 default processes): `2 * 4 = 8`

### AC-4 (Strong governance anti-drift)

- Classification manifest exists and lists all `scripts/relaxtime/run_*.jl` and key scan scripts with tags: `public-authoritative` / `internal-authoritative` / `regression-gate` / `analysis-experimental`.
- Static governance checker fails when a `public-authoritative` script directly includes `src/relaxtime/*`.
- Orchestrator path is the only path used by official output jobs.
- Consistency regression (workflow vs direct) runs on fixed mini-set and passes configured tolerances.

### AC-5 (No parameter intercept)

- Strict mode rejects unknown/unconsumed keys in official profile.
- `consumption_report.json` is produced and `unused_keys` is empty for official runs.
- Fallback events are persisted and can be promoted to failure with `--fail-on-fallback`.
- Probe-value pass-through tests prove deep parameter consumption for required key set.
- For cross-section orchestration, produced row count matches:
  - `N_rows = N_T * N_xi * N_process * N_sqrt_s`
  - where all `N_*` are taken from normalized effective config.

## Preparation Checklist Before Implementation

1. Define canonical key map and alias map across all `run_*.jl` scripts.
2. Lock schema v1 required/optional/default fields.
3. Add dedicated cross-section orchestrator for multi-`T` x multi-`xi`.
4. Add smoke tests for:
   - config load/merge precedence
   - process alias normalization
   - output metadata required keys
   - orchestration path generation and plotting command assembly
5. Add governance assets for strong mode:
   - script classification manifest
   - static governance checker
   - workflow-vs-direct consistency mini-regression
6. Add no-intercept audit assets:
   - strict config mode enforcement
   - consumption report emitter
   - fallback event recorder + `--fail-on-fallback`
   - probe-value end-to-end pass-through tests

## Recommended Next Step

After user review of this spec, switch to implementation planning and break work into phases:

- Phase 1: config loader + schema/alias foundation
- Phase 2: transport/cross-section orchestration unification
- Phase 3: plotting contract automation + smoke coverage
