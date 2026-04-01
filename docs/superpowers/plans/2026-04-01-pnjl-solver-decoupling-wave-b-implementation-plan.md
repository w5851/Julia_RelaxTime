# PNJL Solver Decoupling Wave-B Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Converge Wave-A contracts into shared execution paths (solver + scan), prepare compatibility deprecation rhythm, and freeze fallback/selection stability for next-wave cleanup.

**Architecture:** Reuse Wave-A `ProblemSpec`/governance contracts, route scan-side logic through shared selectors/rules, and maintain compatibility wrappers with explicit migration markers.

**Tech Stack:** Julia (`src/models/*`), NLsolve-based solvers, ForwardDiff/ImplicitDifferentiation chain, tests in `tests/unit|integration|regression`, docs in `docs/dev|docs/architecture|docs/superpowers`.

---

## File Structure (Planned)

- Modify: `src/models/scans/TmuScan.jl`
- Modify: `src/models/scans/TrhoScan.jl`
- Modify: `src/models/solver/Solver.jl`
- Modify: `src/models/constraint_solver.jl`
- Modify: `src/models/solver/CandidateGovernance.jl`
- Modify: `src/models/Models.jl`
- Modify: `docs/dev/active/2026-04-01_PNJL求解器解耦Wave-B任务单.md`
- Create: `docs/dev/active/2026-04-01_求解器兼容层迁移映射表_Wave-B.md`
- Create: `tests/regression/pnjl/test_waveb_scan_governance_stability.jl`
- Create: `tests/integration/models/test_waveb_compat_routing_smoke.jl`

## Chunk 1: Shared Governance Routing in Scan Paths

### Task 1: Unify scan-side candidate selection with governance interfaces

**Files:**
- Modify: `src/models/scans/TmuScan.jl`
- Modify: `src/models/scans/TrhoScan.jl`
- Modify: `src/models/solver/CandidateGovernance.jl`

- [ ] **Step 1: Add failing regression for scan-side governance routing parity**
- [ ] **Step 2: Route scan candidate filtering/selection through shared governance interfaces**
- [ ] **Step 3: Ensure context payload fields are consistently propagated**
- [ ] **Step 4: Re-run scan regression and verify no unintended selection drift**

## Chunk 2: Compatibility Deprecation-Ready Routing

### Task 2: Harden unified-entry first policy while keeping compatibility non-breaking

**Files:**
- Modify: `src/models/solver/Solver.jl`
- Modify: `src/models/constraint_solver.jl`
- Modify: `src/models/Models.jl`

- [ ] **Step 1: Add failing integration smoke for unified-entry and compatibility parity**
- [ ] **Step 2: Strengthen compatibility shims and migration map query surface**
- [ ] **Step 3: Keep external signatures backward-compatible**
- [ ] **Step 4: Re-run parity smoke and confirm consistent outputs**

## Chunk 3: Fallback Diagnostics Stability

### Task 3: Normalize fallback reason semantics for regression tracking

**Files:**
- Modify: `src/models/scans/TrhoScan.jl`
- Modify: `src/models/scans/TmuScan.jl`
- Create: `tests/regression/pnjl/test_waveb_scan_governance_stability.jl`

- [ ] **Step 1: Add failing regression for fallback reason tags at representative points**
- [ ] **Step 2: Normalize reason tags and ordering in scan diagnostics**
- [ ] **Step 3: Re-run regression and verify deterministic behavior**

## Chunk 4: Governance Documentation Synchronization

### Task 4: Update migration/architecture docs for Wave-B state

**Files:**
- Modify: `docs/dev/active/2026-04-01_PNJL求解器解耦Wave-B任务单.md`
- Create/Modify: `docs/dev/active/2026-04-01_求解器兼容层迁移映射表_Wave-B.md`
- Modify: `docs/architecture/models_solver_contract.md` (if boundary semantics changed)

- [ ] **Step 1: Update migration status and removal thresholds for Wave-B progress**
- [ ] **Step 2: Link updated evidence and verification commands**
- [ ] **Step 3: Verify docs governance checks remain green**

## Chunk 5: Verification Matrix and Closure

### Task 5: Execute Wave-B verification matrix

**Files:**
- Create/Modify tests listed above

- [ ] **Step 1: Run targeted new unit/integration/regression tests**
  - `julia --project=. -e 'include("tests/integration/models/test_waveb_compat_routing_smoke.jl")'`
  - `julia --project=. -e 'include("tests/regression/pnjl/test_waveb_scan_governance_stability.jl")'`
- [ ] **Step 2: Run unit smoke profile**
  - `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
- [ ] **Step 3: Run integration smoke profile**
  - `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
- [ ] **Step 4: Run regression smoke profile**
  - `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
- [ ] **Step 5: Run docs/governance checks**
  - `julia --project=. scripts/dev/check_docs_consistency.jl`
  - `julia --project=. scripts/dev/check_active_docs_governance.jl`
- [ ] **Step 6: Update Wave-B task sheet completion evidence**
