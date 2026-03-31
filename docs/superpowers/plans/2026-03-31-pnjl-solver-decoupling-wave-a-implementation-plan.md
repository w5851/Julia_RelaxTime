# PNJL Solver Decoupling Wave-A Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Freeze and implement Wave-A contracts for PNJL solver decoupling so AD/implicit differentiation stays stable while solver/governance paths converge incrementally.

**Architecture:** Introduce a first-class ProblemSpec contract, unify AD boundary semantics, standardize candidate governance interfaces, and keep compatibility paths explicit through a migration table. Prioritize contract-and-test-first migration over broad rewrites.

**Tech Stack:** Julia (`src/models/*`), existing NLsolve-based root solving, ImplicitDifferentiation + ForwardDiff paths, tests in `tests/unit|integration|regression`, docs under `docs/dev|docs/decisions|docs/superpowers`.

---

## File Structure (Planned)

- Create: `src/models/solver/ProblemSpec.jl`
- Create: `src/models/solver/CandidateGovernance.jl`
- Modify: `src/models/solver/Conditions.jl`
- Modify: `src/models/solver/ImplicitSolver.jl`
- Modify: `src/models/solver/Solver.jl`
- Modify: `src/models/constraint_solver.jl`
- Modify: `src/models/Models.jl`
- Modify: `src/models/implicit_gap.jl`
- Modify: `src/models/derivatives/ThermoDerivatives.jl`
- Modify: `src/models/derivatives/ConservedChargeSusceptibilities.jl`
- Modify: `src/models/scans/TmuScan.jl`
- Modify: `src/models/scans/TrhoScan.jl`
- Create: `docs/architecture/models_ad_implicit_contract.md`
- Create: `docs/dev/active/2026-03-31_求解器兼容层迁移映射表.md`
- Create: `tests/unit/models/test_problem_spec_contract.jl`
- Create: `tests/integration/models/test_ad_implicit_contract_smoke.jl`
- Create: `tests/integration/models/test_solver_config_isolation_smoke.jl`
- Create: `tests/regression/pnjl/test_candidate_governance_stability.jl`

## Chunk 1: ProblemSpec Contract Introduction

### Task 1: Add ProblemSpec core type and mode builders

**Files:**
- Create: `src/models/solver/ProblemSpec.jl`
- Modify: `src/models/solver/ConstraintModes.jl` (if lightweight metadata helpers needed)
- Modify: `src/models/Models.jl`

- [x] **Step 1: Write failing unit tests for ProblemSpec shape and mode dimensions**
  - Add assertions for `FixedMu`, `FixedRho`, `FixedAsymmetricRho`, `FixedEntropy`, `FixedSigma`.

- [x] **Step 2: Run tests to verify failures**
  - Run only the new contract test file.

- [x] **Step 3: Implement ProblemSpec struct and constructors**
  - Keep fields minimal and aligned with design doc.

- [x] **Step 4: Wire module include/export in `Models.jl`**

- [x] **Step 5: Re-run unit tests and confirm pass**

### Task 2: Bridge Conditions residual builders to explicit-parameter API

**Files:**
- Modify: `src/models/solver/Conditions.jl`
- Modify: `src/models/solver/ProblemSpec.jl`

- [x] **Step 1: Add failing tests for explicit residual interface parity with existing closure path**

- [x] **Step 2: Implement explicit residual entrypoints per mode**

- [x] **Step 3: Keep existing closure builders as adapters**

- [x] **Step 4: Re-run tests and verify parity**

## Chunk 2: AD/Implicit Differentiation Contract Hardening

### Task 3: Normalize forward/conditions contract

**Files:**
- Modify: `src/models/implicit_gap.jl`
- Modify: `src/models/solver/ImplicitSolver.jl`
- Create: `docs/architecture/models_ad_implicit_contract.md`

- [ ] **Step 1: Write failing integration smoke for FixedMu and flavor-mu Jacobians**

- [ ] **Step 2: Refactor adapter layer so forward/conditions signatures are unified**

- [ ] **Step 3: Enforce Dual-safe conditions boundary (no forced Float64 on AD path)**

- [ ] **Step 4: Document boundary rules in architecture doc**

- [ ] **Step 5: Re-run integration smoke and confirm pass**

### Task 4: Remove global mutable configuration hotspots from derivative call paths

**Files:**
- Modify: `src/models/solver/ImplicitSolver.jl`
- Modify: `src/models/derivatives/ThermoDerivatives.jl`
- Modify: `src/models/derivatives/ConservedChargeSusceptibilities.jl`
- Create: `tests/integration/models/test_solver_config_isolation_smoke.jl`

- [ ] **Step 1: Add failing isolation test for back-to-back distinct config invocations**

- [ ] **Step 2: Introduce solver-scoped config payloads instead of global mutable refs where feasible**

- [ ] **Step 3: Keep backward-compatible wrappers if public signatures are impacted**

- [ ] **Step 4: Re-run isolation test and key derivative smoke tests**

## Chunk 3: Candidate Governance Unification

### Task 5: Introduce CandidateGovernance interfaces and context payload

**Files:**
- Create: `src/models/solver/CandidateGovernance.jl`
- Modify: `src/models/constraint_solver.jl`
- Modify: `src/models/solver/ImplicitSolver.jl`
- Modify: `src/models/Models.jl`

- [ ] **Step 1: Add failing tests for governance interface signatures and default behavior**

- [ ] **Step 2: Implement `hard_rule` and `selector` interface adapters**

- [ ] **Step 3: Add context payload plumbing (mode/continuity/phase/residual policy)**

- [ ] **Step 4: Re-run tests and confirm legacy selector parity**

### Task 6: Align scan and single-point selection behavior through shared governance

**Files:**
- Modify: `src/models/scans/TmuScan.jl`
- Modify: `src/models/scans/TrhoScan.jl`
- Create: `tests/regression/pnjl/test_candidate_governance_stability.jl`

- [ ] **Step 1: Add failing regression for selection stability at representative points**

- [ ] **Step 2: Route scan-side selection through shared governance interfaces**

- [ ] **Step 3: Re-run regression and verify no unintended drift**

## Chunk 4: Compatibility Layer Mapping and Documentation

### Task 7: Build migration mapping table and mark compatibility boundaries

**Files:**
- Create: `docs/dev/active/2026-03-31_求解器兼容层迁移映射表.md`
- Modify: `src/models/solver/Solver.jl`
- Modify: `src/models/constraint_solver.jl`

- [ ] **Step 1: Enumerate old/new entry mappings with removal thresholds**

- [ ] **Step 2: Add non-breaking migration markers in compatibility paths**

- [ ] **Step 3: Verify call sites still route through intended primary entrypoints**

## Chunk 5: Verification and Governance Closure

### Task 8: Run Wave-A verification matrix and governance checks

**Files:**
- Create/Modify tests listed above
- Modify docs index links if required

- [ ] **Step 1: Run unit smoke profile**
  - `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`

- [ ] **Step 2: Run integration smoke profile**
  - `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`

- [ ] **Step 3: Run regression smoke profile**
  - `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`

- [ ] **Step 4: Run docs/governance checks**
  - `julia --project=. scripts/dev/check_docs_consistency.jl`
  - `julia --project=. scripts/dev/check_active_docs_governance.jl`

- [ ] **Step 5: Update active task sheet evidence links and completion states**
