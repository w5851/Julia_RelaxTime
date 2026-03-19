# Relaxtime Core Kernel Refinement Design

## Context

The current `relaxtime` work in the main workspace mixes three concerns:

- core numerical kernel changes in `src/relaxtime/`
- script and workflow evolution under `scripts/relaxtime/` and `scripts/analysis/`
- generated study outputs under `data/outputs/`

The largest active diff sits in the kernel itself, especially `src/relaxtime/AverageScatteringRate.jl`, with supporting changes in `src/relaxtime/RelaxationTime.jl` and the related unit tests. The numerical changes already encode meaningful behavior around `sigma(s)` caching, w0cdf sampling, threshold and resonance layering, adaptive refinement, and tau aggregation. The immediate need is to extract and stabilize that kernel as a reviewable line of work without expanding the scope to scripts, workflow automation, or regenerated artifacts.

## Goal

Create a first `relaxtime` follow-up PR that refines and stabilizes the core kernel while preserving current behavior.

This PR should:

- isolate the kernel changes into a clean review surface
- make internal responsibilities clearer inside the existing files
- convert implicit numerical assumptions into explicit test-covered contracts
- avoid introducing new user-facing APIs unless strictly necessary

## Non-Goals

This design explicitly does not include:

- reorganizing `scripts/relaxtime/`
- reorganizing `scripts/analysis/`
- rerunning full scans or regenerating all `data/outputs/`
- absorbing manual scan workflow or GitHub workflow changes
- pursuing a new round of physics tuning for the T190/T200 neighborhoods

## Scope

In scope for the first PR:

- `src/relaxtime/AverageScatteringRate.jl`
- `src/relaxtime/RelaxationTime.jl`
- `tests/unit/relaxtime/test_average_scattering_rate.jl`
- `tests/unit/relaxtime/test_relaxation_time.jl`
- optionally `tests/regression/relaxtime/test_tau_xi_probe_regression.jl` if the refinement changes fixture-facing semantics or reveals that an existing regression guard must move with the kernel contract

Out of scope for the first PR:

- `scripts/relaxtime/*`
- `scripts/analysis/*`
- `data/outputs/*`
- workflow YAML changes tied to manual scans or batch runs

## Design Summary

The first PR will use a conservative design: keep the current module entrypoints and avoid a multi-file split, but impose clearer internal layering inside the two kernel files.

This is intentionally a behavior-preserving refinement rather than an algorithmic rewrite. The main risk in the current workspace is reviewability and accidental behavior drift, not the absence of new public abstractions.

## Internal Structure

### `AverageScatteringRate`

Refine the file into five explicit logical layers.

1. Parameter and particle helpers
   - flavor mass and chemical-potential lookup
   - isotropic vs anisotropic distribution dispatch
   - low-level kinematic helpers

2. `sigma(s)` cache and interpolation behavior
   - `CrossSectionCache` storage model
   - PCHIP slope maintenance and interpolation
   - local-window averaging and out-of-range behavior

3. Sampling-grid design for cross-section evaluation
   - w0cdf-based `s`-grid design
   - threshold-layer point injection
   - resonance-layer sample injection
   - adaptive sigma refinement passes

4. Quadrature-grid reuse
   - interval grid caching
   - semi-infinite momentum grid caching
   - helper boundaries for default vs explicit grids

5. Public computational entrypoints
   - `number_density`
   - `average_scattering_rate`
   - `build_w0cdf_pchip_cache`

The intended outcome is that each helper cluster has one numerical responsibility, and the public entrypoints mainly orchestrate validated helper calls.

### `RelaxationTime`

Refine the file into four explicit logical layers.

1. Key normalization and alias resolution
   - mapping between canonical rate keys and compatibility aliases
   - lookup helpers for `NamedTuple` and `AbstractDict` inputs

2. Average-rate generation and reuse
   - reuse of existing rates when already available
   - cache handoff to `average_scattering_rate`
   - normalization of struct and `NamedTuple` inputs

3. Tau aggregation core
   - deterministic tau-inverse algebra
   - contribution row derivation from the shared layout table
   - negative-rate clamping policy as an explicit contract

4. Public entrypoints and bounded file leftovers
   - `relaxation_rates`
   - `relaxation_times`

Auxiliary CSV helpers may remain in the file, but they are outside the active refinement scope for the first PR and should be left untouched unless a minimal mechanical edit is unavoidable.

The intended outcome is that tau computation becomes easy to audit as a pure aggregation layer on top of already-defined rate contracts.

## Public API Contract

The first PR should preserve the current public surface as much as possible.

### In-scope callable surface classification

#### `AverageScatteringRate`

- `average_scattering_rate` - preserve unchanged as a stable public entrypoint
- `build_w0cdf_pchip_cache` - preserve unchanged as a stable public entrypoint
- `precompute_cross_section!` - preserve unchanged as an exported supporting entrypoint even if the first PR does not actively refactor it
- `number_density` - preserve unchanged as a callable entrypoint used by kernel logic and tests

#### `RelaxationTime`

- `compute_average_rates` - preserve unchanged as a stable public entrypoint
- `relaxation_rates` - preserve unchanged as a stable public entrypoint
- `relaxation_times` - preserve unchanged as a stable public entrypoint
- `relaxation_rate_contribution_rows` - preserve unchanged as a public diagnostic entrypoint

#### Explicitly not part of the first active refactoring surface

- CSV loader helpers in `RelaxationTime.jl` such as `_read_sigma_table_csv` and `load_cross_section_caches_from_dir`
- any script-facing wrapper outside the two kernel files

If one of these excluded helpers needs a purely mechanical adjustment to keep the file valid after local cleanup, that change should remain minimal and must not be used to broaden the PR scope.

Stable entrypoints to preserve:

- `average_scattering_rate`
- `build_w0cdf_pchip_cache`
- `precompute_cross_section!`
- `number_density`
- `compute_average_rates`
- `relaxation_rates`
- `relaxation_times`
- `relaxation_rate_contribution_rows`

Rules:

- do not rename these entrypoints in the first PR
- do not move users onto new public wrappers unless required by correctness
- prefer internal helper extraction over outward-facing API expansion

## Behavior Contracts To Make Explicit

The refinement should keep existing numerical strategy but make the following contracts more explicit in code and tests.

### Cross-section and interpolation contracts

- `CrossSectionCache` remains the storage boundary for sampled `sigma(s)` tables
- PCHIP interpolation remains the default interpolation mechanism
- queries outside the cache support range continue to return zero rather than clamp to edge values
- local-window behavior, if present, must remain consistent with the current interpolation contract

### Sampling and refinement contracts

- w0cdf grid design remains the base strategy for `sigma(s)` sampling
- threshold layering remains responsible for adding support close to threshold when the first gap indicates under-resolution
- resonance layering remains responsible for adding samples around sufficiently sharp far peaks
- adaptive refinement remains a correction mechanism for under-resolved interpolation error, not a new optimization pass with changed semantics

### Tau aggregation contracts

- alias-backed rate lookup continues to support both canonical keys and compatibility keys
- `relaxation_rates` remains deterministic from `(densities, rates)` alone
- `relaxation_rate_contribution_rows` remains explainable from the shared contribution layout
- negative tau-inverse totals remain warned-and-clamped rather than silently propagated

## Testing Strategy

The first PR should strengthen tests only where they reinforce the current kernel contract.

Required:

- `tests/unit/relaxtime/test_average_scattering_rate.jl`
- `tests/unit/relaxtime/test_relaxation_time.jl`

Required contract matrix:

- struct and `NamedTuple` equivalence for the preserved public entrypoints touched by the PR
- interpolation and cache invariants, including out-of-range `sigma(s) -> 0`
- threshold-layer and resonance-layer injection invariants
- adaptive refinement improving a deliberately under-resolved interpolation case without changing the high-level contract
- default-grid and explicit-grid equivalence where the code already intends that equivalence
- tau aggregation algebra and alias resolution
- negative-rate warn-and-clamp behavior for tau-inverse totals

Conditional regression guard:

- run `tests/regression/relaxtime/test_tau_xi_probe_regression.jl` if the PR changes any of these fixture-facing semantics:
  - interpolation semantics
  - cache support-range semantics
  - rate-key resolution semantics
  - tau aggregation semantics
- otherwise explicitly skip that regression in the PR notes with the statement that no fixture-facing contract changed
- do not expand scope to full scan regeneration just to satisfy the first kernel PR

## Delivery Plan For The First PR

The PR description should explicitly say that this round only stabilizes the core kernel.

Expected reviewer-facing framing:

- what moved into clearer helper boundaries
- which behavior contracts are now explicit
- which tests guard those contracts
- which areas were intentionally deferred to later follow-up PRs

Deferred follow-up lines after this PR:

1. script and workflow alignment under `scripts/relaxtime/`
2. study-output and artifact reconciliation under `data/outputs/`
3. any later targeted numerical investigation beyond behavior-preserving kernel refinement

## Risks

- file-internal cleanup can still shift behavior if helper extraction changes evaluation order or default keyword propagation
- over-eager cleanup could accidentally pull script assumptions into the kernel PR
- regression fixtures can create hidden coupling to generated outputs if used too broadly in the first round

## Mitigations

- keep public entrypoints stable
- prefer file-internal layering over cross-file relocation in the first PR
- use unit tests as the primary guardrail
- only pull in the existing probe regression if kernel-facing semantics truly intersect it
- document non-goals in both the spec and eventual PR message

## Definition of Done

- the first PR is limited to the agreed kernel scope
- `AverageScatteringRate` and `RelaxationTime` are clearer internally without changing intended public behavior
- the two unit test files pass
- the probe regression is run only if the kernel refinement materially touches its contract
- no script/workflow or artifact-regeneration work is smuggled into the PR
- the PR description clearly states that this round stabilizes kernel behavior only

## Decision

Proceed with a first `relaxtime` follow-up PR that extracts the core kernel into a cleaner, behavior-preserving review surface centered on `AverageScatteringRate`, `RelaxationTime`, and their tests, while deferring scripts, workflows, and generated outputs to later follow-up work.
