# PNJL Solver Decoupling Wave-B Design

## Context

Wave-A has delivered core contracts and governance foundations:

1. `ProblemSpec` and explicit residual interfaces are in place.
2. AD/implicit boundary adapters are unified and documented.
3. Candidate governance context and compatibility-safe interfaces are established.
4. Compatibility mapping and migration markers are frozen.

Wave-B should move from "contract hardening" to "execution-path convergence" while preserving Wave-A safety guarantees.

## Goals

1. Converge scan and single-point selection behavior onto shared governance execution paths.
2. Start deprecation-ready routing for legacy `solve_fixed*` compatibility helpers (without forced removal yet).
3. Remove remaining cross-module behavior drift hotspots between `solver/*`, `constraint_solver.jl`, and scan orchestrators.
4. Keep regression behavior stable with explicit fallback controls and traceable diagnostics.

## Non-Goals

1. No removal of legacy entrypoints in Wave-B.
2. No solver backend replacement.
3. No broad rearchitecture outside solver/governance convergence scope.

## Design Decisions

### D1: Shared governance execution contract for scan/orchestrator paths

Wave-B requires scan-side candidate handling to call the same governance contract used by solver-side selection:

- shared hard-rule application
- shared selector semantics
- shared context payload fields

This removes hidden policy forks in scan helpers.

### D2: Deprecation-ready compatibility routing

Legacy helpers remain callable but move to explicit compatibility tier semantics:

- default new call sites use `solve_constraint`
- legacy helpers are compatibility shims with migration markers
- migration status is queryable from code and docs

### D3: Deterministic fallback policy surface

Fallback behavior (continuity/multiseed/weighted-block) should be represented in diagnostics with stable reason tags.

This allows Wave-B regression tests to validate selection/fallback stability in representative points.

### D4: Governance-doc sync as release gate

Any Wave-B governance change must update:

- active task sheet evidence
- migration mapping table
- architecture contract notes if boundary semantics change

## Target Artifacts

1. Shared governance routing updates in scan-side modules.
2. Extended migration/deprecation markers for compatibility helpers.
3. Regression tests for fallback/selection stability under representative scan points.
4. Wave-B active task sheet with checkboxes and verification matrix.

## Verification Strategy

Wave-B minimum matrix:

1. scan-side governance routing regression
2. legacy-vs-unified entry behavior parity at representative points
3. fallback reason stability regression
4. unit/integration/regression smoke re-run
5. docs governance checks

## Risks and Mitigations

1. Risk: scan routing convergence changes candidate selection unexpectedly.
   - Mitigation: lock representative regression points before routing edits.
2. Risk: compatibility markers drift from actual call graph.
   - Mitigation: maintain code-level map + task-sheet evidence links.
3. Risk: fallback diagnostics become noisy/inconsistent.
   - Mitigation: normalize reason tags and add deterministic assertions.
