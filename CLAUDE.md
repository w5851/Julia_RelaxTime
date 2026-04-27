# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Julia-based QCD physics computing platform: PNJL/NJL models, phase diagrams, transport coefficients (Relaxation Time Approximation), scattering amplitudes, meson masses, and magnetized/rotating QCD. The repo is **include-driven** (not a standard registered Julia package) — modules are assembled via explicit `include()` calls into `Main`.

Default communication language with the user: **Chinese**.

## Essential Commands

### Environment setup
```sh
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

### Tests (always prefer smoke profile for local edit-run loops)

**Unit tests** (fastest, deterministic):
```sh
julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'
```

**Single unit test file:**
```sh
julia --project=. -e 'ENV["UNIT_FILES"]="relaxtime/test_transport_coefficients.jl"; include("tests/unit/runtests.jl")'
```

**Integration tests:**
```sh
julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'
```

**Regression tests:**
```sh
julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'
```

**Validation tests:**
```sh
julia --project=. -e 'ENV["VALIDATION_PROFILE"]="smoke"; include("tests/validation/runtests.jl")'
```

All layers support three profiles: `smoke` (fast, <1min), `core` (CI gate), `full` (everything). File selectors use env vars `UNIT_FILES`, `INTEGRATION_FILES`, `REGRESSION_FILES`, `VALIDATION_FILES`.

### Governance / lint-equivalent checks
```sh
julia --project=. scripts/dev/check_docs_consistency.jl
julia --project=. scripts/dev/check_models_entry_contract.jl
julia --project=. scripts/dev/check_solver_contract_leakage.jl
julia --project=. scripts/dev/check_script_entrypoints.jl
```

### Key scripts (stable entrypoints)
```sh
# Phase diagram
julia --project=. scripts/pnjl/calculate_phase_structure.jl --preset=smoke

# T-mu scan
julia --project=. scripts/pnjl/run_tmu_scan.jl --help

# Transport coefficient scan
julia --project=. scripts/relaxtime/run_gap_transport_scan.jl --help

# HTTP server + web frontend
julia --project=. scripts/server/server_full.jl
```

Full command reference and all profiles are in `AGENTS.md`.

## Architecture

### Include-driven module system

There is no `src/PackageName.jl` entry point. Modules live in `Main` namespace and guard against re-inclusion:
```julia
if !isdefined(Main, :ModuleName)
    include("src/path/to/module.jl")
end
```
`using` and `import` go at the top of each module. Relative imports (`using ..GaussLegendre`) are normal. Prefer explicit imports (`using Main.ParameterTypes: QuarkParams`).

### Unified model system (`src/models/`)

All user-facing operations go through **`src/models/entrypoints.jl`** — the `Models` module. Never route through legacy `src/pnjl/PNJL.jl` (removed; historical references are compatibility notes only).

**Abstract type hierarchy:**
```
AbstractQCDModel → AbstractNJLModel → AbstractPNJLModel (5D: phi_u, phi_d, phi_s, Phi, Phibar)
                                     → NJLModel, NJL2Model (3D)
                  → RotationModel
                  → GasLiquidModel
```

Every concrete model follows a canonical per-model contract of 5 files: `capabilities.jl`, `api.jl`, `adapters/entrypoint_adapter.jl`, `workflows/noop.jl`, and `core/` implementation.

### Solver subsystem (`src/models/solver/`)

Multi-mode constraint solver: `FixedMu`, `FixedRho`, `FixedEntropy`, `FixedSigma`, `FixedAsymmetricRho`. Solving flows through:
1. `solve()` / `solve_constraint()` — top-level API
2. `build_problem_spec()` → `ProblemSpec`
3. `ProblemSpec.forward_solve()` → `ConstraintSolver` (backed by NLsolve) or `GapSolver`

Multi-seed governance (`CandidateGovernance`) evaluates solutions from different seed strategies and picks the physical branch. The differentiation engine (`src/models/solver/diff/`) uses ForwardDiff + ImplicitDifferentiation for thermodynamic derivatives.

**Critical invariant (PR76):** non-fixedmu modes use flattened joint-solve semantics. Mixed-meson governance must remain unchanged unless explicitly scoped in the task.

### Phase diagram pipeline (`src/models/phase/`)

Multi-stage pipeline: phase boundary detection → CEP location → crossover line → spinodal decomposition → adaptive rho refinement. The production pipeline orchestrates these stages and emits `phase_summary.json` + `phase_report.md`.

### RelaxTime / transport (`src/relaxtime/`)

Computational pipeline in dependency order:
One-loop integrals → polarization (anisotropic) → A-field builder → effective couplings → meson propagator → total propagator → meson masses → Mott transition → differential cross sections → scattering amplitudes → total cross sections → average scattering rates → relaxation times → transport coefficients (η, ζ, σ).

### Workflow engine (`src/models/workflow_engine/`)

Generic pipeline runner with three pipeline kinds:
- **Workflow**: `solve_gap_and_transport`, `solve_gap_and_meson_point`
- **Scan**: `run_tmu_scan`, `run_trho_scan`
- **RelaxTime orchestrator**: unified `transport` / `cross-section` orchestration

### Configuration system (`config/`)

TOML-based with inheritance. Model params: `config/models/<model>/default.toml` (with `unittest.toml` overrides). Shared physics constants: `config/physics/default.toml`. Workflow profiles: `config/workflows/`. Loaded by `src/config/ConfigLoader.jl`.

## Key Conventions

- **Naming**: modules/types in `PascalCase`, functions in `snake_case`, test files `test_*.jl`, benchmark files `bench_*.jl`
- **Units**: natural units internally; `_inv_fm` suffix for fm⁻¹ values; `T_MeV`, `mu_MeV` for external MeV-facing names; higher-order units like `sigma_fm4`, `coupling_inv_fm4`
- **Tests**: four layers (`unit/`, `integration/`, `regression/`, `validation/`). Smoke tests must be deterministic, fast, no external deps. New core modules require unit tests. Stable public entrypoints must update `docs/api/`.
- **Scripts**: analysis in `scripts/analysis/`, perf probes in `scripts/perf/`, never put non-test scripts in `tests/`
- **Performance**: only `@inline`/`@inbounds` when justified by hot-path profiling. Benchmarks in `benchmark/`, not in smoke tests.
- **Error handling**: validate public inputs early with `ArgumentError` naming the offending parameter. Use `@warn` for recoverable fallback. Never silently ignore invalid physics inputs.
- **Commit messages**: must match historical style. Sample with `git log -10 --oneline` and match the dominant prefix pattern (see AGENTS.md for full governance rules).
