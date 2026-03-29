# Phase-2 Governance/Observability/Split Deployment Design

## Context

Current project has completed (or is planning) Phase-1 direction: localhost-first frontend workflow with API-contract-first integration. Active task sheet still contains broader mid/long-term targets that are not fully covered by the current Phase-1 spec/plan, especially for:

1. Task system governance (active doc Stage C).
2. Observability and operations readiness (active doc Stage E).
3. Split frontend/backend deployment readiness (active doc Stage F).

This Phase-2 design defines the next bounded scope so the repository can move from localhost MVP toward controlled production-readiness without introducing unnecessary platform complexity in one jump.

## Goals

1. Standardize asynchronous scan job governance policy and state model.
2. Introduce practical observability baseline for job lifecycle and failure diagnosis.
3. Establish deployable split architecture contract (frontend static host + Julia API service).
4. Preserve localhost developer workflow while enabling low-friction migration to remote resources.

## Non-Goals

1. No enterprise IAM, tenant isolation, or full security hardening suite.
2. No Kubernetes-first deployment requirement.
3. No rewrite of backend compute kernels or physics model internals.
4. No full BI/reporting platform for metrics analytics.

## Scope

### In Scope

1. Job lifecycle governance additions:
   - cancel API
   - timeout policy
   - idempotency key handling for create requests
   - queue/resource policy visibility
2. Observability baseline:
   - structured logs with `job_id`
   - core runtime metrics (queue depth, duration, success/failure counts)
   - error categorization map
3. Split deployment contract:
   - API base URL and CORS policy strategy
   - reverse proxy baseline (path forwarding + health checks)
   - environment configuration matrix for localhost/staging/remote
4. Documentation and verification assets for governance + deployment runbook.

### Out of Scope

1. Multi-cluster scheduling and autoscaling.
2. Fine-grained RBAC and token lifecycle management.
3. Distributed tracing platform integration (optional later enhancement).

## Architecture Overview

### A. Job Governance Layer

- Keep existing queued/running/succeeded/failed flow and extend with `cancelled` terminal state.
- Enforce create-job idempotency via request key + normalized payload hash.
- Add timeout enforcement for long-running jobs with explicit final reason code.
- Maintain bounded retention and pruning policy with clear metadata in status API.

### B. Observability Layer

- Structured logging contract for every job transition (`created`, `started`, `progress`, `ended`).
- Metric counters and gauges exposed in lightweight form (log-derived or endpoint-derived).
- Error classification map aligned with frontend error rendering (`INVALID_REQUEST`, `QUEUE_FULL`, `JOB_NOT_SUCCEEDED`, `TIMEOUT`, `INTERNAL_ERROR`).

### C. Split Deployment Layer

- Frontend remains static artifact; backend remains Julia HTTP service.
- Reverse proxy routes:
  - `/` -> static frontend
  - `/api/*` -> Julia backend
- CORS strategy:
  - localhost dev: permissive for productivity
  - staged/remote: explicit allowlist origins
- Runtime config strategy:
  - frontend reads environment/boot config for API base URL
  - backend reads env for CORS/origin policy and queue constraints

## Data and Flow Contract (Phase-2)

1. Client creates job with optional idempotency key.
2. Server returns either new job or existing equivalent job (idempotent replay).
3. Job status transitions are logged and measurable.
4. Client can request cancel; server transitions job to cancelled when valid.
5. Terminal state always carries normalized outcome code and timestamps.

## Policy Contracts

### Job State Contract

- Allowed states: `queued`, `running`, `succeeded`, `failed`, `cancelled`.
- Forbidden transitions (examples):
  - `succeeded -> running`
  - `failed -> running`
  - `cancelled -> running`

### Timeout Contract

- Each job kind has a configurable max duration.
- Timeout produces terminal `failed` or `cancelled` (single consistent project choice required in implementation plan) with explicit timeout code.

### Idempotency Contract

- Request key in header or body (single normalized implementation path).
- Same key + equivalent payload returns existing job reference.
- Same key + different payload returns conflict error.

## Verification Strategy

1. Governance tests:
   - create duplicate requests with idempotency key.
   - cancel queued/running job.
   - verify timeout behavior and final state.
2. Observability tests:
   - logs contain job_id and terminal outcome.
   - queue depth and success/failure metrics observable.
3. Deployment contract tests:
   - localhost mode works unchanged.
   - split proxy mode routes frontend and API correctly.
   - CORS allowlist blocks disallowed origins in staged mode.

## Risks and Mitigations

1. Risk: governance features increase complexity and unstable transitions.
   - Mitigation: explicit state transition table + transition tests first.
2. Risk: observability adds overhead/noise.
   - Mitigation: compact structured schema and level-based logging.
3. Risk: split deployment breaks localhost development speed.
   - Mitigation: keep localhost mode as default profile and test gate.

## Acceptance Criteria

1. Job API supports cancel and idempotency with deterministic behavior.
2. Job lifecycle includes observable timestamps and normalized outcome codes.
3. Timeout and queue policy are configurable and externally visible via status/governance payload.
4. Minimal split deployment runbook is executable on one staging-like environment.
5. Localhost workflow remains functional after Phase-2 changes.
