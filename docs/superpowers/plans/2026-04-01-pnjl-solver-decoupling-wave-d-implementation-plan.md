# PNJL Solver Decoupling Wave-D Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close PNJL decoupling migration by removing or hard-deprecating redundant compatibility entrypoints, preserving unified model-driven routing and governance auditability.

**Architecture:** Keep Wave-A/B/C contracts intact and perform threshold-gated cleanup across solver/scan/workflow compatibility layers.

**Tech Stack:** Julia (`src/models/solver/*`, `src/models/scans/*`, `src/models/entrypoints.jl`, `scripts/pnjl/*`), tests in `tests/unit|integration|regression`, docs in `docs/dev|docs/superpowers`.

---

## File Structure (Planned)

- Modify: `src/models/solver/Solver.jl`
- Modify: `src/models/entrypoints.jl`
- Modify: `src/models/scans/TmuScan.jl`
- Modify: `src/models/scans/TrhoScan.jl`
- Modify: `src/models/Models.jl`
- Modify (threshold-gated): `scripts/pnjl/run_tmu_scan.jl`
- Modify (threshold-gated): `scripts/pnjl/run_dense_trho_scan.jl`
- Modify (threshold-gated): `scripts/pnjl/run_adaptive_trho_scan.jl`
- Create/Modify: `tests/integration/models/test_waved_compat_cleanup_smoke.jl`
- Create/Modify: `tests/regression/pnjl/test_waved_cleanup_stability.jl`
- Modify: `docs/dev/active/2026-04-01_PNJL求解器解耦Wave-D任务单.md`
- Create/Modify: `docs/dev/active/2026-04-01_扫描与工作流兼容层迁移映射表_Wave-D.md`

## Chunk 1: Solver/Scan Compatibility Cleanup Gate

### Task D1: Add failing tests for post-Wave-C compatibility removal contract

**Files:**
- Create/Modify: `tests/integration/models/test_waved_compat_cleanup_smoke.jl`
- Create/Modify: `tests/regression/pnjl/test_waved_cleanup_stability.jl`

- [ ] **Step 1: Add failing integration tests for legacy entry behavior after cleanup threshold**
- [ ] **Step 2: Add failing regression for diagnostics/output stability under cleaned routes**
- [ ] **Step 3: Verify failures are expected (missing/old compat behavior)**

## Chunk 2: Threshold-Gated Removal / Hard Deprecation

### Task D2: Execute compat cleanup in solver + scan/workflow entrypoints

**Files:**
- Modify: `src/models/solver/Solver.jl`
- Modify: `src/models/entrypoints.jl`
- Modify: `src/models/scans/TmuScan.jl`
- Modify: `src/models/scans/TrhoScan.jl`
- Modify: `src/models/Models.jl`
- Modify (if thresholds pass): `scripts/pnjl/run_tmu_scan.jl`
- Modify (if thresholds pass): `scripts/pnjl/run_dense_trho_scan.jl`
- Modify (if thresholds pass): `scripts/pnjl/run_adaptive_trho_scan.jl`

- [ ] **Step 1: Apply hard-deprecate-first policy for `solve_fixed*` and high-risk compat paths; remove only after thresholds pass**
- [ ] **Step 2: Keep migration metadata queryable with explicit `removed/hard_deprecated` state**
- [ ] **Step 3: Ensure primary unified model-driven path remains unchanged for callers**
- [ ] **Step 4: Re-run targeted tests and verify green**

## Chunk 3: Migration Mapping and Governance Sync

### Task D3: Sync migration docs and deletion thresholds evidence

**Files:**
- Modify: `docs/dev/active/2026-04-01_扫描与工作流兼容层迁移映射表_Wave-D.md`
- Modify: `docs/dev/active/2026-04-01_PNJL求解器解耦Wave-D任务单.md`

- [ ] **Step 1: Update mapping states (`active/deprecation-ready` -> `removed/hard_deprecated`)**
- [ ] **Step 2: Backfill evidence for each removed/deprecated compat entry**
- [ ] **Step 3: Record rollback note for high-risk cleanup points**
- [ ] **Step 4: Sync `docs/api/` for any stable public entrypoint behavior changes**

## Chunk 4: Verification Matrix and Closure

### Task D4: Execute Wave-D verification matrix

**Files:**
- Create/Modify tests listed above

- [ ] **Step 1: Run targeted Wave-D tests**
  - `julia --project=. -e 'include("tests/integration/models/test_waved_compat_cleanup_smoke.jl")'`
  - `julia --project=. -e 'include("tests/regression/pnjl/test_waved_cleanup_stability.jl")'`
- [ ] **Step 2: Run unit smoke profile**
  - `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
- [ ] **Step 3: Run integration smoke profile**
  - `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
- [ ] **Step 4: Run regression smoke profile**
  - `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
- [ ] **Step 5: Run docs/governance checks**
  - `julia --project=. scripts/dev/check_unit_skip_policy.jl`
  - `julia --project=. scripts/dev/check_docs_consistency.jl`
  - `julia --project=. scripts/dev/check_active_docs_governance.jl`
- [ ] **Step 6: Update Wave-D task sheet completion evidence**
