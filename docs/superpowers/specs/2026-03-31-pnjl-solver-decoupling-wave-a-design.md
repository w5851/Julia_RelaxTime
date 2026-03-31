# PNJL Solver Decoupling Wave-A Design

## Context

The repository already has a functional Models-centered solver path, but Wave-A must prevent architecture drift before deeper implementation.

Current tension points:

1. Responsibility overlap between `constraint_solver.jl` and `solver/ImplicitSolver.jl`.
2. Residual construction uses useful closure patterns, but explicit parameter contracts are not yet first-class.
3. AD/Implicit differentiation works, but global mutable config hotspots threaten reentrancy and reproducibility.
4. Candidate governance exists, but lacks a formal context contract across scan/single-point paths.

## Goals

1. Freeze solver decoupling contracts before broad refactoring.
2. Preserve existing behavior while creating one extensible execution spine.
3. Keep AD/ID chain stable and Dual-safe.
4. Define governance interfaces reusable by both solver and scan workflows.

## Non-Goals

1. No one-shot rewrite of solver internals.
2. No backend replacement (NLsolve remains valid).
3. No immediate removal of compatibility layers in Wave-A.

## Design Decisions

### D1: Introduce ProblemSpec as the primary problem contract

ProblemSpec is the canonical object describing a constraint problem instance:

- mode
- x_dim
- theta_dim
- residual! (explicit-parameter form)
- forward_solve
- conditions
- unpack_solution
- postprocess
- hard_rules
- selector

Closure builders are retained only as adapters over explicit contracts.

### D2: Separate primal and differentiation boundaries

- `forward_solve(theta, cfg, mode)` is primal-only and may normalize to concrete numeric runtime types.
- `conditions(theta, x, meta, cfg, mode)` must be Dual-safe and avoid forced `Float64` conversion on AD path inputs.

### D3: Governance interface becomes explicit and contextual

Candidate acceptance and selection use standardized signatures:

- `hard_rule(candidate, params, context) -> (ok, reason)`
- `selector(candidates, params, context) -> selected_index`

`context` minimally includes mode and continuity/phase hints so single-point and scan selection remain coherent.

### D4: Compatibility layers are managed by mapping table, not ad-hoc retention

A migration mapping table records:

- old entry
- new entry
- compatibility period
- removal threshold

This enables incremental convergence without losing traceability.

## Architecture Sketch

1. L0 Physics: model equations and observables.
2. L1 Equation Builder: mode-to-residual assembly via ProblemSpec.
3. L2 Solve Policy: seed continuation, root policy, fallback.
4. L3 Candidate Governance: hard rules and selector.
5. L4 Workflow/Scan: orchestration only, no model-specific solver forks.

## Verification Strategy

Wave-A minimum matrix:

1. FixedMu Jacobian smoke.
2. flavor-mu Jacobian smoke.
3. multi-mode dimension consistency.
4. candidate selection stability regression.
5. config isolation (no cross-call contamination).

## Risks and Mitigations

1. Risk: Contract introduction increases interim boilerplate.
   - Mitigation: start with minimal ProblemSpec fields and adapter bridge.
2. Risk: AD path regressions from boundary edits.
   - Mitigation: lock smoke tests before large refactor steps.
3. Risk: Behavior drift between scan and solver paths.
   - Mitigation: enforce shared governance interfaces and context payload.

## Deliverables

1. ProblemSpec contract artifact in code.
2. AD/ID contract documentation and boundary rules.
3. Candidate governance interface module.
4. Compatibility mapping table document.
5. Wave-A validation test additions.
