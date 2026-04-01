# Models AD/Implicit Contract

## Scope

This document freezes the Wave-A contract for the implicit differentiation boundary used by `src/models`.

Primary code paths:

- `src/models/implicit_gap.jl`
- `src/models/solver/ImplicitSolver.jl`

## Contract

The implicit boundary is standardized as two functions:

1. `forward_solve(theta, cfg, mode) -> (x, meta)`
2. `conditions(theta, x, meta, cfg, mode) -> residual_vector`

In current Wave-A implementation, `cfg` and `mode` may be closed over by adapter builders, but semantic roles stay the same.

## Boundary Rules

### Rule 1: Primal solve may normalize to Float64

`forward_solve` is primal-only and may use `Float64(...)` conversion for solver robustness.

### Rule 2: conditions must be Dual-safe

`conditions` must not force `theta` or `x` back to `Float64`.

Required pattern:

- keep `T_fm = theta[1]`
- keep chemical potential components from `theta` in their incoming numeric type
- call `gap_residual`/`gap_conditions` with these values directly

### Rule 3: fixed-mu and flavor-mu share adapter shape

Both fixed-mu and flavor-mu routes expose the same pair of adapters:

- `forward_solve`
- `conditions`

This avoids signature drift between scalar-mu and flavor-mu differentiation paths.

## Verification

Wave-A smoke tests that validate this contract:

- `tests/integration/models/test_ad_implicit_contract_smoke.jl`
- `tests/integration/models/test_models_implicitdiff_flavor_mu_smoke.jl`
- `tests/unit/models/test_implicit_gap_flavor_mu.jl`

These tests ensure Jacobian evaluation works through both fixed-mu and flavor-mu adapter paths without type-forcing regressions.
