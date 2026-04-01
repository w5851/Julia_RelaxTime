# PNJL Solver Decoupling Wave-C Design

## Context

Wave-A froze contracts and governance interfaces.
Wave-B converged scan/single-point governance behavior and added deprecation-ready compatibility routing metadata.

Wave-C should shift the decoupling focus from "solver-path convergence" to "scan/workflow convergence":

1. make scan/workflow entrypoints model-driven by default,
2. eliminate model-specific SOP duplication in scripts/workflows,
3. keep migration traceability and non-breaking rollout discipline.

## Goals

1. Converge script and workflow scan orchestration on unified `Models` entrypoints (`run_tmu_scan`, `run_trho_scan`, and workflow-facing model-driven APIs).
2. Replace model-specialized scan SOP forks with parameterized model dispatch.
3. Keep backward compatibility during Wave-C via migration markers and deprecation-ready routing.
4. Lock deterministic regression coverage for representative model-driven scan outputs.

## Non-Goals

1. No solver backend replacement.
2. No aggressive hard deletion of compatibility wrappers in Wave-C.
3. No broad refactor outside scan/workflow convergence and migration governance scope.

## Design Decisions

### D1: Model-driven scan orchestration as the only primary path

Wave-C treats `model_kind` + unified mode interfaces as the source of truth for scan/workflow routing. New call sites must go through the unified model-driven path, not model-specialized script forks.

### D2: SOP fork contraction through adapter normalization

Legacy SOP scripts stay callable short-term, but become thin adapters over unified scan/workflow entrypoints. Their behavior should be externally equivalent at representative points.

### D3: Migration governance mirrors solver compatibility governance

Wave-C introduces script/workflow migration status tracking with explicit state, target wave, and deletion thresholds to prevent untracked drift.

### D4: Regression-first stability for convergence points

Representative points across scan/workflow routes must have deterministic baseline assertions, including parity between legacy SOP adapters and unified model-driven paths.

## Target Artifacts

1. Wave-C active task sheet with ordered execution checkpoints (C1 -> C2 -> C3).
2. Wave-C implementation plan with TDD steps and verification matrix.
3. Migration status notes for script/workflow SOP compatibility-to-unified mapping.
4. Targeted integration/regression tests for model-driven scan/workflow parity.

## Verification Strategy

Wave-C minimum matrix:

1. script/workflow parity smoke for representative model-driven scan points,
2. regression stability for unified scan diagnostics/output fields,
3. unit/integration/regression smoke reruns,
4. docs governance checks.

## Risks and Mitigations

1. Risk: script path contraction changes implicit defaults.
   - Mitigation: freeze representative parity points before routing edits.
2. Risk: migration markers drift from real call graph.
   - Mitigation: keep code-queryable migration status and doc mapping synchronized.
3. Risk: model-driven parameterization introduces hidden coupling.
   - Mitigation: isolate adapter layers and enforce single responsibility for routing.
