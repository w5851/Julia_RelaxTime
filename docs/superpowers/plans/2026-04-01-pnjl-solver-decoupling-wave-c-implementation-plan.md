# PNJL Solver Decoupling Wave-C Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Converge scan/workflow execution onto model-driven unified entrypoints, contract legacy SOP forks into compatibility adapters, and freeze migration+regression stability for Wave-D cleanup.

**Architecture:** Keep solver contracts from Wave-A/B unchanged, and focus on L4 workflow/scan convergence: model-parameterized orchestration, explicit compatibility mapping, deterministic parity diagnostics.

**Tech Stack:** Julia (`src/models/scans/*`, `src/models/workflows/*`, `scripts/pnjl/*`), tests in `tests/unit|integration|regression`, docs in `docs/dev|docs/superpowers`.

---

## File Structure (Planned)

- Modify: `scripts/pnjl/run_tmu_scan.jl`
- Modify: `scripts/pnjl/run_dense_trho_scan.jl`
- Modify: `scripts/pnjl/run_adaptive_trho_scan.jl`
- Modify: `src/models/entrypoints.jl`
- Modify: `src/models/scans/TmuScan.jl`
- Modify: `src/models/scans/TrhoScan.jl`
- Modify: `src/models/Models.jl`
- Modify: `docs/dev/active/2026-04-01_PNJL求解器解耦Wave-C任务单.md`
- Create/Modify: `docs/dev/active/2026-04-01_扫描与工作流兼容层迁移映射表_Wave-C.md`
- Create: `tests/integration/models/test_wavec_scan_workflow_parity_smoke.jl`
- Create: `tests/regression/pnjl/test_wavec_model_driven_scan_stability.jl`

## Chunk 1: Script/Workflow Unified Routing

### Task C1: Converge scan scripts onto model-driven unified entrypoints

**Files:**
- Modify: `scripts/pnjl/run_tmu_scan.jl`
- Modify: `scripts/pnjl/run_dense_trho_scan.jl`
- Modify: `scripts/pnjl/run_adaptive_trho_scan.jl`
- Modify: `src/models/entrypoints.jl`

- [ ] **Step 1: Add failing integration smoke for script path vs unified model-driven path parity**
- [ ] **Step 2: Route script scan entrypoints through unified model-driven APIs**
- [ ] **Step 3: Keep CLI signatures and behavior backward-compatible**
- [ ] **Step 4: Re-run parity smoke and verify no unintended drift**

## Chunk 2: SOP Fork Contraction and Migration Governance

### Task C2: Convert model-specialized SOP forks into compatibility adapters

**Files:**
- Modify: `src/models/scans/TmuScan.jl`
- Modify: `src/models/scans/TrhoScan.jl`
- Modify: `src/models/Models.jl`
- Create/Modify: `docs/dev/active/2026-04-01_扫描与工作流兼容层迁移映射表_Wave-C.md`

- [ ] **Step 1: Add failing tests for compatibility adapter routing metadata/queryability**
- [ ] **Step 2: Implement/extend migration status query surface for scan/workflow adapters**
- [ ] **Step 3: Ensure old SOP paths remain non-breaking but unified-first**
- [ ] **Step 4: Re-run tests and confirm migration status traceability**

## Chunk 3: Model-Driven Scan Stability Baselines

### Task C3: Freeze deterministic regression around unified model-driven scan outputs

**Files:**
- Create: `tests/regression/pnjl/test_wavec_model_driven_scan_stability.jl`
- Modify: `src/models/scans/TmuScan.jl`
- Modify: `src/models/scans/TrhoScan.jl`

- [ ] **Step 1: Add failing regression for representative model-driven scan points**
- [ ] **Step 2: Normalize required diagnostics/output ordering for stable assertions**
- [ ] **Step 3: Re-run regression and verify deterministic behavior**

## Chunk 4: Documentation and Governance Sync

### Task C4: Update Wave-C active task/evidence and migration mapping

**Files:**
- Modify: `docs/dev/active/2026-04-01_PNJL求解器解耦Wave-C任务单.md`
- Create/Modify: `docs/dev/active/2026-04-01_扫描与工作流兼容层迁移映射表_Wave-C.md`

- [ ] **Step 1: Update migration states and deletion thresholds for Wave-C progress**
- [ ] **Step 2: Backfill evidence and command outcomes to task sheet**
- [ ] **Step 3: Verify docs governance checks remain green**

## Chunk 5: Verification Matrix and Closure

### Task C5: Execute Wave-C verification matrix

**Files:**
- Create/Modify tests listed above

- [ ] **Step 1: Run targeted new integration/regression tests**
  - `julia --project=. -e 'include("tests/integration/models/test_wavec_scan_workflow_parity_smoke.jl")'`
  - `julia --project=. -e 'include("tests/regression/pnjl/test_wavec_model_driven_scan_stability.jl")'`
- [ ] **Step 2: Run unit smoke profile**
  - `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
- [ ] **Step 3: Run integration smoke profile**
  - `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
- [ ] **Step 4: Run regression smoke profile**
  - `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
- [ ] **Step 5: Run docs/governance checks**
  - `julia --project=. scripts/dev/check_docs_consistency.jl`
  - `julia --project=. scripts/dev/check_active_docs_governance.jl`
- [ ] **Step 6: Update Wave-C task sheet completion evidence**
