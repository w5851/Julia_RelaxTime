# PNJL Solver Decoupling Wave-E Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Complete compatibility retirement by removing distributed scan script entrypoints, enforcing a single governance-layer path, and extending unified scanning to all currently implemented model families.

**Architecture:** Preserve Wave-A/B/C/D contracts; perform Wave-E as an entrypoint governance consolidation and CLI convergence wave.

**Tech Stack:** Julia (`src/models/entrypoints.jl`, `src/models/scans/*`, `src/models/solver/*`, `scripts/models/*`, `scripts/pnjl/*`), tests in `tests/unit|integration|regression`, docs in `docs/dev|docs/superpowers|docs/api`.

---

## File Structure (Planned)

- Modify: `src/models/entrypoints.jl`
- Modify: `src/models/scans/TmuScan.jl`
- Modify: `src/models/scans/TrhoScan.jl`
- Modify: `src/models/solver/Solver.jl`
- Modify: `src/models/solver/ImplicitSolver.jl`
- Create: `scripts/models/run_unified_scan.jl`
- Remove/Archive (threshold-gated): `scripts/pnjl/run_tmu_scan.jl`
- Remove/Archive (threshold-gated): `scripts/pnjl/run_dense_trho_scan.jl`
- Remove/Archive (threshold-gated): `scripts/pnjl/run_adaptive_trho_scan.jl`
- Create: `tests/integration/models/test_wavee_unified_scan_cli_smoke.jl`
- Create: `tests/regression/pnjl/test_wavee_unified_scan_stability.jl`
- Modify: `docs/dev/active/2026-04-01_PNJL求解器解耦Wave-E任务单.md`
- Create/Modify: `docs/dev/active/2026-04-01_统一扫描入口兼容退场映射表_Wave-E.md`
- Modify (if behavior changes): `docs/api/models/scans/*.md`, `docs/api/models/solver/*.md`

## Chunk 0: PR #45 Closure Gate

### Task E0: Complete Wave-D PR closure as Wave-E precondition

**Files:**
- Modify: `docs/dev/active/2026-04-01_PNJL求解器解耦Wave-E任务单.md`

- [ ] **Step 1: Check PR #45 CI/checks and review status**
- [ ] **Step 2: Resolve review feedback (if any) without scope drift**
- [ ] **Step 3: Merge PR #45 and record merge evidence in Wave-E task sheet**

## Chunk 1: Unified CLI Contract (single entrypoint)

### Task E1: Introduce single unified scan/workflow CLI

**Files:**
- Create: `scripts/models/run_unified_scan.jl`
- Modify: `src/models/entrypoints.jl`

- [ ] **Step 1: Add failing integration tests for single CLI entry contract**
- [ ] **Step 2: Implement subcommands (`scan tmu`, `scan trho`, `workflow phase`) on unified entrypoints**
- [ ] **Step 3: Validate deterministic argument parsing and error messaging**
- [ ] **Step 4: Re-run targeted tests and verify green**

## Chunk 2: Model Coverage Matrix Expansion

### Task E2: Extend unified scanning to all implemented model families

**Target model set:** `PNJL`, `NJL`, `RPNJL`, `PNJLMagnetic`, `Rotation`, `GasLiquid`

**Note:** `pnjl_aniso` is represented as profile/parameter mode (e.g., `xi`) under existing model kinds; not a standalone model_kind.

**Files:**
- Modify: `src/models/scans/TmuScan.jl`
- Modify: `src/models/scans/TrhoScan.jl`
- Modify: `src/models/solver/ImplicitSolver.jl`
- Modify: `src/models/solver/Solver.jl`

- [ ] **Step 1: Add failing matrix tests for model_kind × scan_mode combinations**
- [ ] **Step 2: Implement capability-aware routing/guards and normalize unsupported-mode errors**
- [ ] **Step 3: Ensure deterministic outputs for representative points per supported model**
- [ ] **Step 4: Re-run matrix tests and verify green**

## Chunk 3: Compatibility Retirement (hard-deprecated -> removed/archived)

### Task E3: Retire distributed legacy scan scripts and keep only single path files

**Files:**
- Remove/Archive: `scripts/pnjl/run_tmu_scan.jl`
- Remove/Archive: `scripts/pnjl/run_dense_trho_scan.jl`
- Remove/Archive: `scripts/pnjl/run_adaptive_trho_scan.jl`
- Create/Modify: `docs/dev/active/2026-04-01_统一扫描入口兼容退场映射表_Wave-E.md`

- [ ] **Step 1: Add failing tests asserting legacy script retirement behavior**
- [ ] **Step 2: Remove/archive compat scripts after thresholds pass**
- [ ] **Step 3: Keep migration status queryable with explicit `removed/archived` metadata**
- [ ] **Step 4: Re-run parity/regression checks and verify no routing drift**

## Chunk 4: Docs/API and Governance Closure

### Task E4: Sync docs, task evidence, and governance status

**Files:**
- Modify: `docs/dev/active/2026-04-01_PNJL求解器解耦Wave-E任务单.md`
- Create/Modify: `docs/dev/active/2026-04-01_统一扫描入口兼容退场映射表_Wave-E.md`
- Modify (if needed): `docs/api/models/scans/*.md`, `docs/api/models/solver/*.md`

- [ ] **Step 1: Update migration states and retirement thresholds evidence**
- [ ] **Step 2: Backfill commands/results and DoD evidence into task sheet**
- [ ] **Step 3: Ensure docs governance checks remain green**

## Chunk 5: Verification Matrix and Final Closure

### Task E5: Execute Wave-E verification matrix

**Files:**
- Create/Modify tests listed above

- [ ] **Step 1: Run targeted Wave-E tests**
  - `julia --project=. -e 'include("tests/integration/models/test_wavee_unified_scan_cli_smoke.jl")'`
  - `julia --project=. -e 'include("tests/regression/pnjl/test_wavee_unified_scan_stability.jl")'`
- [ ] **Step 2: Run unit smoke profile**
  - `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
- [ ] **Step 3: Run integration smoke profile**
  - `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
- [ ] **Step 4: Run regression smoke profile**
  - `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
- [ ] **Step 5: Run governance checks**
  - `julia --project=. scripts/dev/check_unit_skip_policy.jl`
  - `julia --project=. scripts/dev/check_docs_consistency.jl`
  - `julia --project=. scripts/dev/check_active_docs_governance.jl`
- [ ] **Step 6: Update Wave-E task sheet completion evidence**
