# P-mu Branch Competition Diagnostic Design

## Context

Current CEP and first-order boundary detection in the phase pipeline relies on `mu(rho)` curves, `detect_s_shape`, and `maxwell_construction`. This path now gives a stable high-accuracy reference near `xi=0.0`, `legacy+p24/t8`, with non-interpolating production CEP converging to about `T_CEP = 131.028125 MeV`, `mu_CEP = 291.154919 MeV`.

The next question is not to replace Maxwell immediately, but to add a compare-only diagnostic based on fixed-`T` branch competition in `P(mu)` or equivalently `Omega(mu)`. The purpose is to compare two physically meaningful first-order criteria in the CEP neighborhood and explain where they agree or diverge.

## Goal

Add a compare-only `P-mu` / `Omega-mu` phase diagnostic that:

- tracks hadron-like and quark-like branches at fixed `T`
- identifies stable transition points from branch competition
- identifies metastable endpoints from branch existence limits
- writes branch-level artifacts that can be compared directly with Maxwell results

## Scope

First phase only:

- `xi=0.0`
- `solver_backend=:legacy`
- `p_num=24`, `t_num=8`
- CEP neighborhood only, initially around `T = 130.8:0.1:131.1 MeV`
- compare-only diagnostic output; no replacement of production CEP result

## Non-goals

- no replacement of `maxwell_construction` in production mode
- no immediate support for all `xi` values or all model backends
- no Hessian or curvature-based stability criterion in v1
- no direct coupling to adaptive `rho` optimization in v1

## High-Level Design

Introduce a new diagnostic layer in the phase topic for fixed-`T`, fixed-`mu` branch competition.

Recommended public-facing entrypoint shape:

```julia
analyze_pm_branch_competition(; T_values, mu_grid, xi, solver_backend, p_num, t_num, seed_pair, output_dir)
```

This entrypoint should:

1. iterate over target temperatures
2. scan a fixed `mu` grid for each temperature
3. solve each `mu` point twice, once from a hadron-like seed and once from a quark-like seed
4. preserve branch continuity instead of treating each point independently
5. summarize stable transition and metastable endpoints
6. compare the result with the Maxwell output at the same temperature

The module should remain independent from `CEPDetector.jl` and `ProductionPhasePipeline.jl`, with only optional integration hooks for artifact generation.

## Operational Definitions

To keep the first implementation reproducible, v1 uses the following explicit rules.

- `branch point accepted`
  - solve returns `converged=true`
  - `omega`, `pressure`, `rho_norm`, and `residual_norm` are finite
- `branch exists`
  - the branch point at a given `(T, mu, branch)` is accepted
- `branch break`
  - the point fails acceptance, or
  - the accepted solution jumps outside the allowed continuity band relative to the previous point on the same logical branch
- `valid competing branch`
  - an accepted branch point that also passes branch continuity against its own previous accepted point
- `metastable endpoint`
  - the last accepted point before a branch break inside a region where the competing branch also exists nearby
- `no_bistable_window`
  - no finite `mu` interval exists where both hadron-like and quark-like branches are accepted simultaneously

For v1, endpoint cause must be recorded explicitly as one of:

- `physical_loss_candidate`
- `nonconvergence`
- `branch_jump`
- `out_of_grid`

These labels remain diagnostic, not final physical truth claims.

## Core Algorithm

### 1. Fixed-T selection

Run only in the CEP neighborhood for the first iteration. Initial reference set:

- `130.8 MeV`
- `130.9 MeV`
- `131.0 MeV`
- `131.1 MeV`

### 2. Fixed-mu branch scan

For each temperature and each `mu` in `mu_grid`:

- solve once from a hadron-like seed
- solve once from a quark-like seed
- record for each branch:
  - `x_state`
  - `pressure`
  - `omega`
  - `rho_norm`
  - `residual_norm`
  - `converged`

### 2.1 Mu-grid contract

The first version must freeze a reproducible scan contract.

Recommended v1 contract:

- initial `mu` range: `289.5:0.1:292.5 MeV`
- optional local refinement around a sign change in `Delta Omega` or `Delta P`
- local refinement step target: `0.01 MeV`
- interpolation for a transition between adjacent accepted points:
  - linear interpolation in `Delta Omega(mu)` for `mu_transition_pm`

This keeps the first comparison deterministic and directly comparable with the current CEP neighborhood results.

### 2.2 Seed-pair contract

The first version must make seed handling explicit.

Recommended v1 `seed_pair` schema:

- `hadron_seed0`: initial hadron-like seed vector
- `quark_seed0`: initial quark-like seed vector
- `continuity_mode`: `:branch_local`
- `fallback_mode`: `:none` for v1 compare-only diagnostics

Branch propagation rule:

- each accepted hadron-branch solution seeds the next `mu` point on the hadron branch
- each accepted quark-branch solution seeds the next `mu` point on the quark branch
- failed points do not overwrite the branch seed cache

The first implementation may obtain `hadron_seed0` and `quark_seed0` from nearby known low-`mu` and high-`mu` equilibria at the same `T`.

### 3. Branch continuity tracking

Branch identity must be preserved across neighboring `mu` points.

Rules:

- hadron-like and quark-like tracks are separate logical branches
- a sudden jump to the opposite phase is a branch break, not a silent overwrite
- failed or unstable solves are preserved as branch termination evidence

Recommended continuity observables for v1:

- `norm(x_state_i - x_state_{i-1})`
- `abs(rho_norm_i - rho_norm_{i-1})`
- optional phase-label proxy from condensate / Polyakov-loop ordering

The exact thresholds should be configurable and written into the artifact metadata.

### 4. Phase criteria

For temperatures where both branches coexist over a finite `mu` interval:

- stable transition point: `Delta Omega = Omega_h - Omega_q` crosses zero
- equivalent check: `Delta P = P_h - P_q` crosses zero
- hadron metastable endpoint: largest `mu` where the hadron-like branch still exists as a valid competing branch
- quark metastable endpoint: smallest `mu` where the quark-like branch exists as a valid competing branch
- bistable window: overlap interval where both branches exist simultaneously

If there is no overlap interval, report `no_bistable_window`.

### 5. Maxwell comparison

For each temperature, compare:

- `mu_transition_pm`
- `mu_transition_maxwell`
- `Delta mu`
- bistable-window width
- which branch disappears first
- whether the point is `pm_only`, `maxwell_only`, `both`, or `neither`

For v1, Maxwell comparison must use the same physical configuration:

- same `xi`
- same backend
- same `p_num`, `t_num`
- same temperature set

Agreement should be reported using a configurable tolerance, with recommended default:

- `abs(mu_transition_pm - mu_transition_maxwell) <= 0.05 MeV`

## Output Artifacts

Recommended first-phase artifacts:

- `pm_branch_scan.csv`
  - one row per `(T, mu, branch)`
- `pm_phase_summary.json`
  - one summary object per `T`
- `pm_vs_maxwell.csv`
  - direct side-by-side comparison table

Suggested summary fields per temperature:

- `T_MeV`
- `mu_transition_pm_MeV`
- `hadron_endpoint_mu_MeV`
- `quark_endpoint_mu_MeV`
- `bistable_window_width_MeV`
- `mu_transition_maxwell_MeV`
- `delta_mu_pm_minus_maxwell_MeV`
- `comparison_status`
- `notes`

Required branch-scan fields:

- `branch`
- `branch_status`
- `endpoint_cause`
- `seed_source`
- `continuity_ok`

## Integration Plan

Minimal integration for v1:

- keep it as a standalone diagnostic path
- allow `run_phase_pipeline` or a dedicated script to call it optionally
- do not alter current CEP result selection logic
- keep current `phase_summary.json` unchanged unless the diagnostic is explicitly requested

## Verification Plan

### Functional

- verify that branch scan runs for `xi=0.0`, `legacy+p24/t8`, CEP-neighborhood temperatures
- verify that branch artifacts are written and readable
- verify that branch coexistence and transition detection behave sensibly at:
  - `130.9 MeV`
  - `131.0 MeV`
  - `131.1 MeV`

### Comparative

- compare `P-mu` transition points against current Maxwell results
- identify temperatures where both methods agree within tolerance
- identify temperatures where one method fails and the other still resolves a signal

### Acceptance

Phase 1 is acceptable when:

- the diagnostic runs end-to-end on the agreed CEP neighborhood
- it outputs stable transition and metastable endpoint summaries for the tested temperatures
- it produces a direct comparison table against Maxwell
- it gives enough evidence to discuss method strengths, failure modes, and disagreement sources

## Risks

- branch identity may drift if seed strategy is not stable enough
- branch disappearance may mix physical loss of metastability with numerical solve failure
- `P-mu` stability may depend strongly on seed construction even if it is less sensitive to `rho` discretization
- compare-only output could still be misread as a replacement unless clearly labeled

## Decision

Proceed with a compare-only `P-mu` / `Omega-mu` diagnostic based on dual-branch competition, scoped to the CEP neighborhood and benchmarked directly against Maxwell before any promotion to a primary phase criterion.
