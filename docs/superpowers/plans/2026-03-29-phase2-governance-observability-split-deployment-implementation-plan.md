# Phase-2 Governance/Observability/Split Deployment Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement Phase-2 capabilities for job governance, observability baseline, and split deployment readiness while preserving localhost workflow.

**Architecture:** Extend existing fullserver job subsystem incrementally: first enforce deterministic state/policy contracts, then add observability hooks, then finalize split deployment configuration and runbooks. All changes remain API-contract-first and backward-compatible for localhost usage.

**Tech Stack:** Julia HTTP backend (`src/simulation/fullserver/*.jl`), existing Models scan execution path, frontend runtime config (`web/js`), docs under `docs/dev` and `docs/superpowers`.

---

## File Structure (Planned)

- Modify: `src/simulation/fullserver/pnjl_scan_jobs.jl`
- Modify: `src/simulation/fullserver/routing.jl`
- Modify: `src/simulation/fullserver/http_utils.jl` (if shared response helpers needed)
- Modify: `src/simulation/fullserver/shared.jl` (if shared error or metadata helpers needed)
- Modify: `src/simulation/ServerLauncher.jl` (if runtime policy env display needed)
- Modify: `web/js/runtime_config.js` and/or `web/js/api.js` (only for Phase-2-compatible config surface)
- Create: `docs/dev/active/2026-03-29_phase2_governance_observability_split_deployment_tasklist.md`
- Create: `docs/reference/ops/` runbook markdown (exact file chosen during implementation)
- Test/Check: targeted unit/integration tests under `tests/unit/` and `tests/integration/`

## Chunk 1: Governance Contract Hardening

### Task 1: Define and enforce state transition table

**Files:**
- Modify: `src/simulation/fullserver/pnjl_scan_jobs.jl`
- Create: `tests/unit/simulation/test_scan_job_state_machine.jl`

- [x] **Step 1: Write failing unit tests for legal/illegal transitions**
  - Include transitions for terminal states and cancel flow.

- [x] **Step 2: Run tests to confirm failures**
  - Run only new test file through unit selector.

- [x] **Step 3: Implement state transition validator**
  - Centralize transition logic in one helper.

- [x] **Step 4: Route all status writes through validator**
  - Prevent ad-hoc state mutation bypass.

- [x] **Step 5: Re-run tests and ensure pass**

- [ ] **Step 6: Commit**
  - `git add src/simulation/fullserver/pnjl_scan_jobs.jl tests/unit/simulation/test_scan_job_state_machine.jl`
  - `git commit -m "refactor: enforce scan job state transition contract"`

### Task 2: Add idempotency for job creation

**Files:**
- Modify: `src/simulation/fullserver/pnjl_scan_jobs.jl`
- Modify: `src/simulation/fullserver/routing.jl` (if header/body passthrough adjustments needed)
- Create: `tests/integration/relaxtime/test_pnjl_scan_idempotency.jl`

- [x] **Step 1: Write failing integration tests for idempotent replay and conflict**
  - same key + same payload -> same job reference
  - same key + different payload -> conflict

- [x] **Step 2: Implement idempotency key extraction and payload fingerprint**

- [x] **Step 3: Add lookup/store path in create-job handler**

- [x] **Step 4: Return deterministic response payload for replayed requests**

- [x] **Step 5: Re-run integration tests and ensure pass**

- [ ] **Step 6: Commit**
  - `git add src/simulation/fullserver/pnjl_scan_jobs.jl src/simulation/fullserver/routing.jl tests/integration/relaxtime/test_pnjl_scan_idempotency.jl`
  - `git commit -m "feat: add idempotent create behavior for PNJL scan jobs"`

### Task 3: Add cancel and timeout policy

**Files:**
- Modify: `src/simulation/fullserver/pnjl_scan_jobs.jl`
- Modify: `src/simulation/fullserver/routing.jl`
- Create: `tests/integration/relaxtime/test_pnjl_scan_cancel_timeout.jl`

- [x] **Step 1: Write failing tests for cancel endpoint behavior**
  - cancel queued/running allowed
  - cancel terminal state rejected (or no-op with explicit semantics)

- [x] **Step 2: Write failing tests for timeout policy**
  - forced timeout leads to deterministic terminal code/state.

- [x] **Step 3: Implement cancel route and handler**
  - define endpoint path and response schema.

- [x] **Step 4: Implement timeout enforcement hook in runner loop**
  - use configurable env policy.

- [x] **Step 5: Re-run tests and ensure pass**

- [ ] **Step 6: Commit**
  - `git add src/simulation/fullserver/pnjl_scan_jobs.jl src/simulation/fullserver/routing.jl tests/integration/relaxtime/test_pnjl_scan_cancel_timeout.jl`
  - `git commit -m "feat: add PNJL scan cancel and timeout governance"`

## Chunk 2: Observability Baseline

### Task 4: Structured logging with job lifecycle events

**Files:**
- Modify: `src/simulation/fullserver/pnjl_scan_jobs.jl`
- Modify: `src/simulation/fullserver/shared.jl` (if shared formatter/helpers needed)
- Create: `tests/unit/simulation/test_scan_job_logging_contract.jl`

- [x] **Step 1: Write failing tests for lifecycle log payload shape**
  - required fields: `job_id`, `kind`, `state`, `timestamp`, `event_code`.

- [x] **Step 2: Implement structured log helper**

- [x] **Step 3: Emit logs at create/start/progress/end transitions**

- [x] **Step 4: Re-run tests and ensure pass**

- [ ] **Step 5: Commit**
  - `git add src/simulation/fullserver/pnjl_scan_jobs.jl src/simulation/fullserver/shared.jl tests/unit/simulation/test_scan_job_logging_contract.jl`
  - `git commit -m "feat: add structured lifecycle logging for scan jobs"`

### Task 5: Metrics and diagnostics payload enrichment

**Files:**
- Modify: `src/simulation/fullserver/pnjl_scan_jobs.jl`
- Modify: `src/simulation/fullserver/pnjl_handlers.jl` (if diagnostics shared path needed)
- Create: `tests/integration/relaxtime/test_pnjl_scan_metrics_diagnostics.jl`

- [ ] **Step 1: Write failing tests for queue/metrics fields in status payload**

- [ ] **Step 2: Add counters/gauges for success/failure/queue depth/duration buckets (minimal)**

- [ ] **Step 3: Expose governance + diagnostics metadata in status responses**

- [ ] **Step 4: Re-run tests and ensure pass**

- [ ] **Step 5: Commit**
  - `git add src/simulation/fullserver/pnjl_scan_jobs.jl src/simulation/fullserver/pnjl_handlers.jl tests/integration/relaxtime/test_pnjl_scan_metrics_diagnostics.jl`
  - `git commit -m "feat: enrich scan status diagnostics and metrics visibility"`

## Chunk 3: Split Deployment Readiness and Docs

### Task 6: Split deployment config contract

**Files:**
- Modify: `web/js/runtime_config.js` and/or `web/js/api.js`
- Modify: `src/simulation/ServerLauncher.jl`
- Create: `tests/integration/relaxtime/test_frontend_backend_config_contract.jl` (if practical)

- [ ] **Step 1: Define failing checklist for environment matrix**
  - localhost profile and staged/remote profile both represented.

- [ ] **Step 2: Implement backend env policy exposure (CORS/queue/timeouts) with safe defaults**

- [ ] **Step 3: Ensure frontend config supports split deployment target cleanly**

- [ ] **Step 4: Verify localhost still works as default mode**

- [ ] **Step 5: Commit**
  - `git add web/js/runtime_config.js web/js/api.js src/simulation/ServerLauncher.jl tests/integration/relaxtime/test_frontend_backend_config_contract.jl`
  - `git commit -m "feat: add split deployment runtime config contract"`

### Task 7: Runbook and active task alignment

**Files:**
- Create: `docs/reference/ops/phase2_split_deployment_runbook.md`
- Create: `docs/dev/active/2026-03-29_phase2_governance_observability_split_deployment_tasklist.md`
- Modify: `docs/dev/active/2026-03-29_前后端分离中长期目标与前端短期待开发任务单.md`

- [ ] **Step 1: Write runbook for localhost/staging split modes**
  - include startup, env vars, health checks, rollback steps.

- [ ] **Step 2: Create Phase-2 active tasklist mapped to Stage C/E/F**

- [ ] **Step 3: Sync original active doc with Phase-2 cross-reference links**

- [ ] **Step 4: Run docs consistency governance check**
  - `julia --project=. scripts/dev/check_docs_consistency.jl`

- [ ] **Step 5: Commit**
  - `git add docs/reference/ops/phase2_split_deployment_runbook.md docs/dev/active/2026-03-29_phase2_governance_observability_split_deployment_tasklist.md docs/dev/active/2026-03-29_前后端分离中长期目标与前端短期待开发任务单.md`
  - `git commit -m "docs: add phase-2 governance and split deployment runbook"`

## Definition of Done

- Stage C governance capabilities have explicit API and tests.
- Stage E observability baseline is implemented and visible in payload/log outputs.
- Stage F split deployment contract is documented and locally verifiable.
- Localhost default workflow remains intact after all Phase-2 changes.

## Execution Status (2026-03-29)

- This Phase-2 plan has not entered implementation yet.
- Current session completed only planning/setup artifacts:
  - `docs/superpowers/specs/2026-03-29-phase2-governance-observability-split-deployment-design.md`
  - `docs/dev/active/2026-03-29_phase2_governance_observability_split_deployment_tasklist.md`
- All task checkboxes in this plan remain pending by design.
