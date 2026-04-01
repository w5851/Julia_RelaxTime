# PNJL Solver Decoupling Wave-D Design

## Context

Wave-A/B/C established the decoupling baseline:

1. Wave-A froze contract boundaries (`ProblemSpec`, constraint mode, governance hooks).
2. Wave-B converged solver governance and compatibility routing metadata.
3. Wave-C converged scan/workflow script routes to model-driven unified entrypoints, and kept compat adapters deprecation-ready.

Wave-D should finish the migration by converting "deprecation-ready" into "deprecation-executed" under controlled governance.

## Goals

1. Remove redundant compatibility wrappers and duplicate entrypoints that are no longer needed after Wave-C.
2. Keep a single model-driven primary path for solver + scan + workflow orchestration.
3. Preserve auditability: every removal must map to migration evidence, tests, and docs updates.
4. Freeze post-cleanup stability with deterministic regression/smoke coverage.

## Non-Goals

1. No solver backend replacement.
2. No new model feature scope beyond migration cleanup.
3. No broad refactor unrelated to decoupling closure and governance hardening.

## Design Decisions

### D1: Compatibility removal is threshold-gated, not time-gated

Wave-D only removes adapters that satisfy explicit thresholds: call-site migration complete, parity coverage green, docs mapping synced.

### D2: Primary path uniqueness is enforced at entrypoint level

`Models.solve_constraint` + model-driven scan/workflow entrypoints remain the only stable public path. For `solve_fixed*` compat helpers, Wave-D default policy is hard-deprecate first, then remove only after thresholds are satisfied.

### D3: Governance metadata must remain queryable after cleanup

Even after deleting compat code, migration status query surfaces should continue to report historical mapping and removal status for auditability.

### D4: Cleanup must be regression-first

Every adapter removal/change is preceded by failing tests that capture expected post-removal behavior (or explicit deprecation behavior), then minimal implementation, then full smoke/governance verification.

## Target Artifacts

1. Wave-D active task sheet with ordered checkpoints.
2. Wave-D implementation plan with TDD-first chunks.
3. Updated migration mapping (state from `active/deprecation-ready` -> `removed` or `hard-deprecated`).
4. Targeted integration/regression tests for post-cleanup stability.
5. API docs sync notes for stable public entrypoint changes (`docs/api/`).

## Verification Strategy

Wave-D minimum matrix:

1. targeted adapter-removal/deprecation tests,
2. scan/workflow and solver routing regression checks,
3. unit/integration/regression smoke reruns,
4. docs governance checks.
5. API docs consistency checks for exported stable entrypoints.

## Risks and Mitigations

1. Risk: hidden external dependency still calls compat path.
   - Mitigation: keep explicit migration map + failure messaging + staged removal order.
2. Risk: cleanup changes diagnostic outputs and breaks baselines.
   - Mitigation: lock deterministic message normalization and update baselines only with evidence.
3. Risk: over-cleanup crosses Wave-D scope.
   - Mitigation: restrict file list and enforce non-goal checks in task sheet.
