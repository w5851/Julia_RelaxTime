# Frontend Localhost 80/20 Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a localhost-first PNJL scan frontend workflow with minimum remote compatibility via configuration-only API endpoint switching.

**Architecture:** Keep frontend business flow API-contract-first by centralizing transport in `web/js/api.js`, adding configurable base URL support, and implementing a minimal jobs workflow UI in existing pages. Prioritize localhost reliability and explicit error handling while deferring production deployment concerns.

**Tech Stack:** Julia HTTP backend (`src/simulation/fullserver`), browser ES modules (`web/js/*.js`), static HTML/CSS (`web/index.html`, `web/css/style.css`), repository docs under `docs/dev` and `docs/superpowers`.

---

## File Structure (Planned)

- Modify: `web/js/api.js`
- Modify: `web/js/ui.js`
- Modify: `web/index.html`
- Modify: `web/css/style.css` (only for minimal new UI states)
- Create: `web/js/runtime_config.js` (or equivalent config module for API base URL)
- Modify: `docs/dev/active/2026-03-29_前后端分离中长期目标与前端短期待开发任务单.md` (status checkboxes during execution)
- Test/Check: manual localhost flow (`scripts/server/server_full.jl` + browser)

## Chunk 1: API Contract Adapter and Config Decoupling

### Task 1: Add runtime base URL configuration (P0-1a)

**Files:**
- Create: `web/js/runtime_config.js`
- Modify: `web/js/api.js`

- [x] **Step 1: Write failing expectation checklist**
  - Define checks: default base URL resolves to `http://localhost:8080`; override mechanism exists; API call paths build from one source.

- [x] **Step 2: Run baseline check against current code**
  - Confirm current hardcoded `API_BASE_URL` in `web/js/api.js` fails override requirement.

- [x] **Step 3: Implement minimal runtime config module**
  - Export `get_api_base_url()` with default localhost and optional override from query/localStorage/global runtime variable (choose one minimal path and document it).

- [x] **Step 4: Wire API module to config module**
  - Replace hardcoded constant usage with config getter.

- [x] **Step 5: Verify URL switching**
  - Manual test: default localhost works; override value changes request target without modifying business methods.

- [ ] **Step 6: Commit**
  - `git add web/js/runtime_config.js web/js/api.js`
  - `git commit -m "feat: make frontend API base URL runtime-configurable"`

### Task 2: Extend jobs API methods and unified error mapping (P0-1)

**Files:**
- Modify: `web/js/api.js`

- [x] **Step 1: Write failing method checklist**
  - Define missing methods: `createScanJob`, `getJobStatus`, `getJobResult`.
  - Define expected normalized error shape for UI consumption.

- [x] **Step 2: Add method signatures and endpoint mapping**
  - `POST /api/modules/pnjl-scan/jobs`
  - `GET /api/modules/pnjl-scan/jobs/{job_id}`
  - `GET /api/modules/pnjl-scan/jobs/{job_id}/result`

- [x] **Step 3: Implement timeout and retry policy**
  - Keep conservative retry for transient network errors only.
  - Ensure no infinite retry loop.

- [x] **Step 4: Implement error normalization**
  - Map backend `code/message/diagnostics` with fallback text for unknown shapes.

- [x] **Step 5: Verify API adapter behavior manually**
  - Simulate invalid params and unreachable server paths.

- [ ] **Step 6: Commit**
  - `git add web/js/api.js`
  - `git commit -m "feat: add PNJL scan jobs API adapter and error normalization"`

## Chunk 2: Localhost Jobs Workflow UI

### Task 3: Build minimal jobs form and validation flow (P0-2 + P0-3)

**Files:**
- Modify: `web/index.html`
- Modify: `web/js/ui.js`
- Modify: `web/css/style.css` (if needed)

- [x] **Step 1: Define UI acceptance checklist**
  - Form supports `kind=tmu|trho`.
  - `xi` strategy is mutually exclusive (`xi` vs `xi_values` vs `xi_grid`).
  - Submit blocked on invalid state.

- [x] **Step 2: Add minimal form fields in HTML**
  - Keep existing page structure; append scan panel instead of redesigning full page.

- [x] **Step 3: Implement frontend validation helpers in UI module**
  - Validate required fields and strategy exclusivity before API call.

- [x] **Step 4: Implement submit handler for job creation**
  - Call `createScanJob`, store returned `job_id`, enter polling mode.

- [x] **Step 5: Verify validation and submission manually**
  - Invalid form blocked locally; valid form reaches backend.

- [ ] **Step 6: Commit**
  - `git add web/index.html web/js/ui.js web/css/style.css`
  - `git commit -m "feat: add PNJL scan submission form with local validation"`

### Task 4: Add polling/status/result rendering and retry UX (P0-2 + P0-4 + P0-5)

**Files:**
- Modify: `web/js/ui.js`
- Modify: `web/index.html`
- Modify: `web/css/style.css` (if needed)

- [x] **Step 1: Define polling state machine checklist**
  - States: idle, submitting, polling, succeeded, failed.
  - Poll stop on terminal states.

- [x] **Step 2: Implement polling loop**
  - Use controlled interval, track attempts, stop on terminal status.

- [x] **Step 3: Implement success rendering**
  - Show `stats(total/success/failure/skipped)` and `output_path`.

- [x] **Step 4: Implement error UX mapping**
  - User-facing concise message + collapsible technical details (`code/diagnostics`).

- [x] **Step 5: Add retry control**
  - Explicit user-triggered retry from failed state.

- [x] **Step 6: Verify end-to-end loop manually**
  - Create job, observe progress, complete and render result.

- [ ] **Step 7: Commit**
  - `git add web/js/ui.js web/index.html web/css/style.css`
  - `git commit -m "feat: add scan job polling, result panel, and retry UX"`

## Chunk 3: Localhost Governance and Verification

### Task 5: Add localhost mode guidance and resource governance (P0-7)

**Files:**
- Modify: `web/index.html`
- Modify: `web/js/ui.js`

- [x] **Step 1: Add localhost mode guidance copy**
  - Explain local compute behavior and expected long-running tasks.

- [x] **Step 2: Add soft concurrency guidance in UI**
  - Recommend 1-2 concurrent jobs and show warning if user exceeds intended pattern.

- [x] **Step 3: Verify guidance clarity**
  - Confirm messages are visible and non-blocking.

- [ ] **Step 4: Commit**
  - `git add web/index.html web/js/ui.js`
  - `git commit -m "docs: add localhost execution guidance for scan workflow"`

### Task 6: Execute smoke acceptance scenarios and update active task status (P0-6)

**Files:**
- Modify: `docs/dev/active/2026-03-29_前后端分离中长期目标与前端短期待开发任务单.md`

- [x] **Step 1: Run local server for validation**
  - Run: `julia --project=. scripts/server/server_full.jl`
  - Expected: server starts and exposes API/static routes.

- [x] **Step 2: Run localhost scenario A (happy path)**
  - Submit one scan job and verify terminal success + result rendering.

- [x] **Step 3: Run localhost scenario B (invalid request)**
  - Trigger frontend validation and ensure request blocked pre-submit.

- [x] **Step 4: Run localhost scenario C (failure/diagnostics path)**
  - Trigger a controlled backend failure path and verify UI diagnostics rendering.

- [x] **Step 5: Run localhost scenario D (result not ready / polling state)**
  - Verify UI behavior before job reaches success.

- [x] **Step 6: Update active doc checkboxes based on evidence**
  - Mark completed P0 items and add short execution evidence notes.

- [ ] **Step 7: Commit**
  - `git add docs/dev/active/2026-03-29_前后端分离中长期目标与前端短期待开发任务单.md`
  - `git commit -m "docs: record localhost 80/20 frontend workflow verification"`

## Definition of Done

- Frontend can submit and track PNJL scan jobs on localhost end-to-end.
- API base URL is configurable without business-flow code changes.
- Frontend enforces `kind` and `xi` strategy constraints before submission.
- Errors are mapped into user-readable and technical-detail layers.
- Active development task sheet reflects current completion status with evidence.

## Execution Status (2026-03-29)

- Completed in working tree: Task 1-6 implementation and verification steps (except commit steps).
- Evidence commands run:
  - `node web/js/api.scan.test.mjs` (PASS)
  - `node web/js/api.scan.integration.test.mjs` with local server warm-up (PASS)
- Backend unblock included during execution:
  - `src/simulation/FullServerApp.jl` loading path fixed for `Main.Models`.
- Remaining action in this plan: commit steps are pending.
