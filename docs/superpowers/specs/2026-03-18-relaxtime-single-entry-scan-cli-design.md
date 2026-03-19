# Relaxtime Single-Entry Scan CLI Design

## Context

After PR #12, the relaxtime core kernel is already separated and stabilized. The next line is script/workflow consolidation.

Current execution paths are split across:

- `scripts/relaxtime/run_gap_transport_scan.jl`
- `scripts/relaxtime/scan_relaxation_times_vs_T.jl`
- `scripts/relaxtime/run_manual_relaxation_scan_workflow.jl`
- `.github/workflows/relaxtime-manual-scans.yml`

They contain duplicated parameter contracts (`mode`, tau node controls, sigma-grid controls, overwrite/resume/bulk toggles) and partially inconsistent defaults/validation styles.

At the same time, some parameter divergence is intentional and linked to numerical-error mitigation. This means consolidation must not accidentally regress to low-precision defaults.

## Goal

Introduce a single scan entrypoint with subcommands while preserving numerical behavior and precision safeguards.

This round should:

- make `scripts/relaxtime/run_scan.jl` the unified dispatch surface
- centralize shared CLI parameter contracts
- keep old script entrypoints functional through thin compatibility wrappers
- switch the manual GitHub workflow to the new entrypoint
- prevent precision-regression caused by over-simplified default unification

## Non-Goals

This design does not include:

- changing physics algorithms or kernel equations
- deleting old script files in this round
- regenerating large output baselines under `data/outputs/`
- full workflow matrix redesign beyond switching to the new entrypoint

## Scope

In scope:

- create `scripts/relaxtime/run_scan.jl` (single entrypoint + subcommands)
- create `scripts/relaxtime/common_scan_options.jl` (shared contract/validation/normalization)
- refactor the three existing scripts into callable execution-layer functions
- keep legacy script CLI compatibility through forwarding wrappers with deprecation notice
- update `.github/workflows/relaxtime-manual-scans.yml` to call the new entrypoint

Out of scope:

- deleting legacy wrappers
- broad output governance or artifact cleanup
- additional CI workflows unrelated to manual scans

## Design Summary

Adopt a three-layer CLI architecture.

1. Dispatch layer (`run_scan.jl`)
- Owns command selection and help.
- Forwards to execution functions only.

2. Shared contract layer (`common_scan_options.jl`)
- Owns common option schema for scanning-related controls.
- Owns normalization (including old aliases such as `--mumin`/`--mumax`/`--mustep` -> canonical `muB` contract where relevant).
- Owns validation and unified `ArgumentError` style that names offending parameter.

3. Execution layer (existing scripts, function-oriented)
- Keeps actual scan/workflow operations.
- Receives already-normalized options.
- Does not duplicate common option parsing logic.

## Command Surface

New unified command surface:

- `run_scan.jl gap-transport ...`
  - maps to current `run_gap_transport_scan.jl`
- `run_scan.jl tau-vs-t ...`
  - maps to current `scan_relaxation_times_vs_T.jl`
- `run_scan.jl manual-workflow ...`
  - maps to current `run_manual_relaxation_scan_workflow.jl`

Legacy wrappers remain available in this round and forward into the new entrypoint.

## Subcommand Option Matrix (Minimum Auditable Contract)

The implementation must include and preserve the following option/default profiles.

### Common precision-sensitive options

- `--mode`: `semi_infinite | finite_15 | finite_lambda`
- `--tau-p-nodes`
- `--tau-angle-nodes`
- `--tau-phi-nodes`
- `--tau-n-sigma`
- `--sigma-grid-n`

### `gap-transport` (from `run_gap_transport_scan.jl`)

- default `mode`: `semi_infinite`
- default `tau-p-nodes`: module default (`RT_ASR.DEFAULT_P_NODES`)
- default `tau-angle-nodes`: module default (`RT_ASR.DEFAULT_ANGLE_NODES`)
- default `tau-phi-nodes`: module default (`RT_ASR.DEFAULT_PHI_NODES`)
- default `tau-n-sigma`: module default (`RT_TCS.DEFAULT_T_INTEGRAL_POINTS`)
- default `sigma-grid-n`: module default (`RT_ASR.DEFAULT_SIGMA_GRID_N`)

### `tau-vs-t` (from `scan_relaxation_times_vs_T.jl`)

- default `mode`: `semi_infinite`
- default `tau-p-nodes`: module default (`RT_ASR.DEFAULT_P_NODES`)
- default `tau-angle-nodes`: module default (`RT_ASR.DEFAULT_ANGLE_NODES`)
- default `tau-phi-nodes`: module default (`RT_ASR.DEFAULT_PHI_NODES`)
- default `tau-n-sigma`: module default (`RT_TCS.DEFAULT_T_INTEGRAL_POINTS`)
- default `sigma-grid-n`: module default (`RT_ASR.DEFAULT_SIGMA_GRID_N`)

### `manual-workflow` (from `run_manual_relaxation_scan_workflow.jl`)

- default `mode`: `finite_15`
- default `tau-p-nodes`: `20`
- default `tau-angle-nodes`: `4`
- default `tau-phi-nodes`: `8`
- default `tau-n-sigma`: `6`
- default `sigma-grid-n`: `60`

Note: differences between `manual-workflow` and other subcommands are intentional and treated as precision-sensitive behavior contract, not accidental inconsistency.

## Precision-Safety Contract (Critical)

The repository currently includes high-precision parameter choices introduced to mitigate numerical error in outputs.

This round must preserve those effective defaults and should not collapse everything into lower precision in the name of “uniformity”.

Required policy:

- treat precision-relevant defaults as part of behavior contract
- when defaults differ for justified numerical reasons, keep per-command defaults explicit
- shared contract layer may centralize representation, but must preserve per-command precision profile
- any default change affecting numerical precision requires explicit rationale and dedicated review note

## Validation Contract

Validation behavior should be uniform and explicit:

- illegal `mode` rejected with named-parameter `ArgumentError`
- non-positive nodes/steps rejected with named-parameter `ArgumentError`
- malformed `xi` list/range rejected with named-parameter `ArgumentError`
- deprecated aliases accepted only via explicit normalization path (with clear summary output)

Alias/canonical conflict policy:

- alias-only: normalize to canonical and continue
- canonical-only: continue
- both-same: continue, canonical wins in summary output
- both-conflict: throw `ArgumentError` naming both conflicting parameters

## Compatibility Strategy

For legacy script files:

- keep the file paths and top-level CLI usability
- print one deprecation line indicating the new `run_scan.jl` equivalent
- forward arguments to the corresponding subcommand path

Wrapper compatibility contract:

- preserve `ARGS` forwarding semantics (no silent dropping of user flags)
- preserve non-zero exit codes from downstream execution
- preserve stderr visibility for downstream failures

This allows existing users and automation to continue working while migrating to the single entrypoint.

## Workflow Strategy

For `.github/workflows/relaxtime-manual-scans.yml`:

- keep existing workflow-dispatch input surface in this round (avoid UI churn)
- only switch internal invocation target to `scripts/relaxtime/run_scan.jl manual-workflow ...`
- preserve current parameter intent and defaults

Required mapping table in implementation notes:

- `integration_mode` -> `--mode`
- `tau_p_nodes` -> `--tau-p-nodes`
- `tau_angle_nodes` -> `--tau-angle-nodes`
- `tau_phi_nodes` -> `--tau-phi-nodes`
- `tau_n_sigma` -> `--tau-n-sigma`
- `sigma_grid_n` -> `--sigma-grid-n`
- `compute_bulk` false -> `--no-bulk`
- `overwrite` true -> `--overwrite`

## Testing Strategy

### Required tests

- command dispatch smoke for all three subcommands
- shared-option validation tests for invalid values
- alias-normalization tests for old parameter names
- legacy-wrapper forwarding tests (entrypoint compatibility)
- precision-default non-regression tests:
  - verify command-specific precision defaults remain at intended levels
  - verify no accidental fallback to lower precision defaults

### Minimum coverage matrix

- dispatch smoke:
  - `run_scan.jl gap-transport --help`
  - `run_scan.jl tau-vs-t --help`
  - `run_scan.jl manual-workflow --help`
- alias normalization:
  - alias-only accepted and normalized
  - canonical-only accepted
  - alias+canonical conflict rejected with named-parameter `ArgumentError`
- wrapper equivalence:
  - one representative command per legacy script reaches same downstream execution path and exit semantics
- precision profile guard:
  - assert each subcommand default profile matches the option matrix above
  - compare against pre-refactor defaults, not newly inferred values

### Verification commands (expected)

- targeted script-smoke commands for each subcommand
- targeted tests for option parsing/validation helpers
- workflow dry-run argument assembly check (where feasible in local harness)

## Risks

- accidental precision regression while centralizing defaults
- compatibility break if wrapper forwarding misses edge arguments
- partial migration where workflow and local CLI diverge

## Mitigations

- explicit per-command precision profile table in implementation
- test coverage for wrapper forwarding and alias normalization
- keep workflow input keys stable in this round
- require PR note when any precision default is intentionally changed

## Definition of Done

- `scripts/relaxtime/run_scan.jl` exists and supports three subcommands
- shared scan options are centralized in `common_scan_options.jl`
- legacy script entrypoints still run via forwarding wrappers
- manual workflow uses the new entrypoint
- invalid shared-parameter inputs fail with unified named-parameter `ArgumentError`
- precision-sensitive defaults do not regress to lower-accuracy profiles
- no physics algorithm changes and no large output regeneration are introduced

## Decision

Proceed with single-entry scan CLI consolidation (`run_scan.jl` + subcommands) using compatibility wrappers and precision-preserving option contracts as the first scripts/workflow follow-up line after core-kernel stabilization.
