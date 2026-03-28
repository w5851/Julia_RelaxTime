# Mott Transition Temperature-xi Scan Script Design

## Context

Current repository already has:

- a meson scan script that can output `M_pi`, `M_K`, `Gamma_pi`, `Gamma_K`, and quark masses (`m_u`, `m_d`, `m_s`):
  - `scripts/relaxtime/run_gap_meson_mass_scan.jl`
- a generic plotting script that supports grouping and splitting dimensions for line plots:
  - `scripts/plot_scan_csv.py`

User request is to freeze one demand-oriented Mott scan workflow and avoid future drift.

## Confirmed User Requirements (Anti-Drift Contract)

The following requirements are frozen from user input and treated as acceptance targets.

1. Physics scan scope
   - Fix baryon chemical potential at `mu_B = 0 MeV`.
   - Scan temperature as the continuous independent variable (`T`).
   - Scan anisotropy factor `xi` as discrete independent variable.
   - Initial frozen `xi` policy:
     - `xi_min = -0.3`
     - `xi_max = 0.3`
     - `xi_step = 0.15`
     - equivalent list: `[-0.3, -0.15, 0.0, 0.15, 0.3]`

2. Required observables (dependent variables)
   - `M_pi`
   - `M_K`
   - `Gamma_pi`
   - `Gamma_K`
   - `M_u + M_d`
   - `M_u + M_s`

3. Required plotting organizations
   - Mode A: different observables in one figure, different `xi` in different figures.
   - Mode B: different `xi` in one figure, different observables in different figures.

## Requirement Clarity Assessment

Current requirement is clear enough for direct implementation. Key contracts are now frozen as follows.

1. Temperature continuity contract (frozen)
   - Keep physical meaning as continuous variable, but execution uses finite grid.
   - Executable contract is fixed to:
     - `T_min_MeV`, `T_max_MeV`, `T_step_MeV`
   - v1 recommended defaults:
     - `T_min_MeV = 120.0`
     - `T_max_MeV = 260.0`
     - `T_step_MeV = 2.0`

2. Output unit and naming contract (partially frozen)
   - Existing script stores masses in fm^-1 (`m_u`, `m_d`, `m_s`) while meson fields are named `M_*`/`Gamma_*`.
   - Need to freeze unit metadata in output CSV for all six observables.
   - Recommended v1: keep native script values and write explicit metadata keys (`unit.<field>`).

3. Plot file layout/naming contract (frozen in this spec)
   - User gave plotting semantics, but deterministic folder/file schema is not yet defined.
   - This spec now freezes deterministic folder and filename conventions.

Conclusion:
- Requirement is clear enough to produce a formal design and implement a v1 workflow.
- Remaining work is implementation and validation, not requirement disambiguation.

## Feasibility From Existing Scripts

### Already achievable directly

- Computation path exists in `scripts/relaxtime/run_gap_meson_mass_scan.jl`.
  - It already supports `muB` scan control and `xi` list control.
  - It already outputs `M_pi`, `M_K`, `Gamma_pi`, `Gamma_K`, `m_u`, `m_d`, `m_s`.
- Plotting path exists in `scripts/plot_scan_csv.py`.
  - It supports `--group` for multi-curve in one figure.
  - It supports `--split` for separating figures by a dimension (usable for `xi`).

### Needs orchestration and contract freezing

- Need one dedicated script/orchestrator preset for this exact Mott demand to avoid ad hoc command drift.
- Need a post-processing step for derived fields:
  - `M_u_plus_M_d = m_u + m_d`
  - `M_u_plus_M_s = m_u + m_s`
- Need two fixed plotting pipelines (Mode A, Mode B) with deterministic output directories.

## Decision

Adopt a dedicated demand script contract:

- one authoritative scan command for `mu_B=0` and frozen default `xi` policy;
- one derived-observable postprocess step;
- two authoritative plotting commands matching Mode A and Mode B;
- persist effective config snapshot for reproducibility.

## Proposed v1 Entrypoints

1. Compute scan
   - `julia --project=. scripts/relaxtime/run_mott_phase_scan.jl --config <path> --outdir <path>`

2. Build derived fields
   - `julia --project=. scripts/relaxtime/run_mott_phase_derived_csv.jl --in <scan_csv> --out <derived_csv>`

3. Plot Mode A / Mode B (authoritative wrapper)
   - `julia --project=. scripts/relaxtime/run_mott_phase_plot_modes.jl --in <derived_csv> --out-dir <fig_dir>`
   - wrapper internally calls `scripts/plot_scan_csv.py` with frozen Mode A / Mode B parameters.

## Configuration Contract (v1)

Recommended config section:

- `[scan.mott_phase]`
  - `muB_MeV = 0.0`
  - `xi_list = [-0.3, -0.15, 0.0, 0.15, 0.3]`
  - `T_min_MeV = 120.0`
  - `T_max_MeV = 260.0`
  - `T_step_MeV = 2.0`
  - `mesons = ["pi", "K"]`

Precedence policy:

- `built-in defaults < TOML < CLI`

Required metadata artifacts:

- `effective_config.json`
- `run_manifest.json`

## Output Contract (v1)

### Required CSV columns

- `T_MeV`, `muB_MeV`, `xi`
- `M_pi`, `M_K`, `Gamma_pi`, `Gamma_K`
- `m_u`, `m_d`, `m_s`
- `M_u_plus_M_d`, `M_u_plus_M_s`
- `status`, `error_code`, `error_message`, `timestamp_utc`

### Deterministic figure layout

- Mode A:
  - directory: `.../mode_a/`
  - one figure per `xi`
  - filename: `mott_mode_a__xi{xi_tag}.png`

- Mode B:
  - directory: `.../mode_b/`
  - one figure per observable
  - filename: `mott_mode_b__{observable}.png`

Numeric tag format:

- integer-like value -> integer string
- otherwise decimal point replaced by `p`
  - example: `-0.15 -> -0p15`

## Acceptance Criteria

### AC-1 Physics scan

- Single command generates scan rows for all default `xi` values at `muB_MeV = 0.0`.
- All rows satisfy `muB_MeV == 0.0`.
- Required base fields and derived fields exist in output CSV.

### AC-2 Plot Mode A

- One figure per `xi` is produced.
- Each figure contains exactly six curves:
  - `M_pi`, `M_K`, `Gamma_pi`, `Gamma_K`, `M_u_plus_M_d`, `M_u_plus_M_s`.

### AC-3 Plot Mode B

- One figure per observable is produced.
- Each figure contains all default `xi` curves.

### AC-4 Reproducibility

- `effective_config.json` and `run_manifest.json` are written.
- Re-run with same config under resume policy does not duplicate point keys.

## Preparation Checklist Before Implementation

- [x] Freeze temperature executable contract as `T_min/T_max/T_step`.
- [x] Add dedicated `muB=0` Mott scan entrypoint with frozen default `xi` list.
- [x] Add derived-field builder (`M_u_plus_M_d`, `M_u_plus_M_s`).
- [x] Add Mode A and Mode B plotting wrappers.
- [x] Add smoke tests for row count and figure count contracts.
- [x] Document commands and metadata schema in user-facing docs.

## Recommended Next Step

Core v1 implementation is in place with smoke-level acceptance checks passing. Next step can move to stronger governance/consistency gates if needed.
