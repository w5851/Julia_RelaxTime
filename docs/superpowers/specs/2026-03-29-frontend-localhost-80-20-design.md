# Frontend Localhost 80/20 Design

## Context

Current repository is in a same-machine integrated mode for development: frontend static pages and Julia backend API run locally, with no remote compute/server support in this phase. Existing backend already exposes API routes for synchronous compute and asynchronous PNJL scan jobs.

The active development task list is documented at:

- `docs/dev/active/2026-03-29_前后端分离中长期目标与前端短期待开发任务单.md`

User decision has been fixed to an 80/20 strategy:

- 80% effort: local usability and stable localhost workflow.
- 20% effort: minimum remote compatibility (config decoupling only), without production-grade infra complexity.

## Goals

1. Deliver a usable frontend task workflow for PNJL scan jobs on localhost.
2. Keep frontend dependent only on HTTP API contract, not on `run_*.jl` scripts.
3. Add minimum future-proofing so API endpoint switching does not require business-flow refactor.
4. Preserve path to future frontend/backend split deployment after remote resources become available.

## Non-Goals (Current Phase)

1. No production authentication or multi-tenant security design.
2. No cloud deployment, gateway orchestration, or remote cluster scheduling.
3. No broad frontend framework migration unless needed for P0 workflow.
4. No full observability platform rollout in this phase.

## Constraints

1. Development and computation are local-machine first (`localhost`) for current stage.
2. Existing backend contract is source of truth (`/api/modules/pnjl-scan/jobs`, status/result routes, and current error payload shape).
3. Frontend changes should prioritize low-risk incremental delivery.

## Selected Architecture

### A. Contract-first frontend layer

- Frontend API module becomes single contract adapter for:
  - `createScanJob`
  - `getJobStatus`
  - `getJobResult`
- Uniform timeout, retry policy, and error mapping is centralized in API module.
- UI components do not handle raw transport semantics directly.

### B. Localhost-first runtime mode

- Default base URL is local (`http://localhost:8080`).
- Job polling and UI behavior optimized for single-machine resource constraints.
- Concurrency defaults conservative (1-2 recommended visible limit in UI guidance).

### C. Minimum remote compatibility layer (20%)

- API base URL must be configurable without editing feature logic.
- Endpoint switch should be configuration-only change.
- Error model retains `code/message/diagnostics` path for future remote debugging parity.

### D. Deferred split deployment stage

- Full frontend/backend split remains explicit follow-up milestone.
- Current design ensures no coupling debt that blocks later split.

## Functional Scope (P0/P1 Mapping)

### P0

1. API client refactor with jobs API methods.
2. Job submission and polling UI closed loop.
3. Form validation for `kind=tmu|trho` and `xi` strategy exclusivity.
4. Result panel with stats and output path.
5. Unified error UX with user-friendly + technical details.
6. Localhost-oriented smoke scenarios for integration acceptance.
7. Base URL configuration support and local resource hints.

### P1

1. Parameter templates.
2. Recent history list by `job_id`.
3. Basic result plotting.

## Data and Flow Contract

1. User submits scan request from frontend form.
2. Frontend validates request locally.
3. Frontend calls `createScanJob` and gets `job_id`.
4. Frontend polls `getJobStatus` until terminal state.
5. On success, frontend calls `getJobResult` and renders output metadata.
6. On failure, frontend shows mapped error with diagnostics details.

## Error Handling Contract

1. Parameter validation errors are blocked in frontend before request.
2. Backend `4xx` maps to user-actionable messages.
3. Backend `5xx` maps to system failure messages with diagnostics panel.
4. Poll timeout and network errors map to retryable state in UI.
5. Queue-full (`429`) has explicit guidance message.

## Localhost Resource Governance

1. UI indicates long-running behavior and polling state.
2. UI guidance constrains parallel submissions to avoid local saturation.
3. Retry is explicit user action; no uncontrolled infinite retry loops.

## Verification Plan

1. Localhost happy path:
   - create job -> poll -> get result succeeds.
2. Queue full path:
   - API returns `429`, frontend shows clear constrained-capacity message.
3. Job failed path:
   - frontend displays user + technical error views.
4. Result-not-ready path:
   - frontend handles non-succeeded result status safely.
5. Base URL switch path:
   - changing config value updates API target without feature-code edits.

## Risks and Mitigations

1. Risk: frontend overfits localhost and accumulates remote migration debt.
   - Mitigation: keep base URL and transport behavior centralized in API adapter.
2. Risk: job polling causes noisy UI or resource overhead.
   - Mitigation: controlled polling interval + terminal-state stop condition.
3. Risk: backend error payload variance causes fragile UI parsing.
   - Mitigation: defensive parser with fallback message strategy.

## Acceptance Criteria

1. Frontend P0 workflow runs end-to-end on localhost without direct script coupling.
2. Base URL can be switched by configuration; business workflow code remains unchanged.
3. Frontend validation enforces `kind` and `xi` strategy constraints prior to submission.
4. Error UI covers invalid request, queue full, failed job, timeout/network cases.
5. One reproducible localhost integration record exists for the P0 workflow.
