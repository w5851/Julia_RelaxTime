# Relaxtime Single-Entry Scan CLI Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Introduce a single relaxtime scan CLI entrypoint with subcommands while preserving precision-sensitive defaults, legacy script compatibility, and workflow behavior.

**Architecture:** Build a thin command dispatcher (`run_scan.jl`) on top of a shared option-contract layer (`common_scan_options.jl`), then adapt existing scripts into execution functions and compatibility wrappers. Keep behavior-preserving precision profiles per subcommand, with explicit guard tests to prevent regression to lower-precision defaults.

**Tech Stack:** Julia 1.10+, existing `scripts/relaxtime/*.jl`, GitHub Actions workflow YAML, `Test` via script-smoke and targeted CLI-contract tests.

---

## File Map

- Create: `scripts/relaxtime/run_scan.jl`
  - single entrypoint + subcommand dispatch (`gap-transport`, `tau-vs-t`, `manual-workflow`)
- Create: `scripts/relaxtime/common_scan_options.jl`
  - shared option schema, alias normalization, validation, precision-profile defaults
- Modify: `scripts/relaxtime/run_gap_transport_scan.jl`
  - expose callable execution function + keep legacy wrapper forwarding
- Modify: `scripts/relaxtime/scan_relaxation_times_vs_T.jl`
  - expose callable execution function + keep legacy wrapper forwarding
- Modify: `scripts/relaxtime/run_manual_relaxation_scan_workflow.jl`
  - expose callable execution function + keep legacy wrapper forwarding
- Modify: `.github/workflows/relaxtime-manual-scans.yml`
  - switch invocation to `run_scan.jl manual-workflow ...` with unchanged dispatch inputs
- Create: `tests/integration/relaxtime/test_single_entry_scan_cli_smoke.jl`
  - dispatch smoke + wrapper compatibility smoke + alias conflict checks

## Precision Contract Table (Must Preserve)

- `gap-transport`
  - mode: `semi_infinite`
  - tau nodes / sigma settings: module defaults (`RT_ASR` / `RT_TCS`)
- `tau-vs-t`
  - mode: `semi_infinite`
  - tau nodes / sigma settings: module defaults (`RT_ASR` / `RT_TCS`)
- `manual-workflow`
  - mode: `finite_15`
  - tau nodes: `20/4/8`
  - tau-n-sigma: `6`
  - sigma-grid-n: `60`

Any drift from the above is a behavior change and must be explicitly justified.

## Chunk 1: Shared Option Contract Layer

### Task 1: Create `common_scan_options.jl` with failing contract tests first

**Files:**
- Create: `scripts/relaxtime/common_scan_options.jl`
- Create: `tests/integration/relaxtime/test_single_entry_scan_cli_smoke.jl`

- [ ] **Step 1: Write failing test for subcommand default profiles**

Add assertions for the precision contract table per subcommand.

- [ ] **Step 2: Capture pre-refactor default snapshots from current scripts**

Before implementation changes, extract old-script effective defaults for:
- mode
- tau node parameters
- sigma-grid parameters

Use these snapshots as comparison baseline in tests to avoid self-referential validation.

- [ ] **Step 3: Write failing test for alias normalization and conflict policy**

Cover `alias-only`, `canonical-only`, `both-same`, `both-conflict` behavior.

- [ ] **Step 4: Write failing test for normalization summary output**

Assert normalized summary exposes canonical parameter names when aliases are used.

- [ ] **Step 5: Write failing test for validation error style**

Assert invalid values throw `ArgumentError` with offending parameter name in message.

Explicitly include malformed `xi` list/range cases.

- [ ] **Step 6: Run tests to confirm expected failures**

Run: `julia --project=. -e 'include("tests/integration/relaxtime/test_single_entry_scan_cli_smoke.jl")'`
Expected: FAIL due to missing shared option layer.

- [ ] **Step 7: Implement `common_scan_options.jl` minimal schema**

Implement:
- common option structs
- subcommand-specific default profile selectors
- alias normalization helpers
- validation helpers with named-parameter `ArgumentError`

- [ ] **Step 8: Re-run tests and ensure pass for implemented cases**

Run same command as Step 6.
Expected: PASS for option-contract tests.

- [ ] **Step 9: Sample commit-message style before commit**

Run: `git log -10 --oneline`
Use an existing recent prefix pattern for the commit message.

- [ ] **Step 10: Commit shared contract layer**

```bash
git add scripts/relaxtime/common_scan_options.jl tests/integration/relaxtime/test_single_entry_scan_cli_smoke.jl
git commit -m "refactor(relaxtime): centralize scan option contracts and precision profiles"
```

## Chunk 2: Single Entry Dispatcher

### Task 2: Add `run_scan.jl` subcommand dispatch with TDD

**Files:**
- Create: `scripts/relaxtime/run_scan.jl`
- Modify: `tests/integration/relaxtime/test_single_entry_scan_cli_smoke.jl`

- [ ] **Step 1: Write failing dispatch smoke tests for three subcommands**

Add smoke checks for:
- `run_scan.jl gap-transport --help`
- `run_scan.jl tau-vs-t --help`
- `run_scan.jl manual-workflow --help`

- [ ] **Step 2: Run smoke tests and confirm dispatch failures**

Run: `julia --project=. -e 'include("tests/integration/relaxtime/test_single_entry_scan_cli_smoke.jl")'`
Expected: FAIL because dispatcher is missing.

- [ ] **Step 3: Implement `run_scan.jl` dispatcher**

Implement:
- subcommand parse + help
- parameter handoff into `common_scan_options.jl`
- execution-call hooks (initially stubs if needed)

- [ ] **Step 4: Re-run dispatch smoke tests**

Run same command.
Expected: PASS for dispatch/help paths.

- [ ] **Step 5: Sample commit-message style before commit**

Run: `git log -10 --oneline`
Use an existing recent prefix pattern for the commit message.

- [ ] **Step 6: Commit dispatcher**

```bash
git add scripts/relaxtime/run_scan.jl tests/integration/relaxtime/test_single_entry_scan_cli_smoke.jl
git commit -m "refactor(relaxtime): add single-entry scan dispatcher"
```

## Chunk 3: Adapt Existing Scripts to Execution Functions + Legacy Wrappers

### Task 3: Adapt `run_gap_transport_scan.jl`

**Files:**
- Modify: `scripts/relaxtime/run_gap_transport_scan.jl`
- Modify: `tests/integration/relaxtime/test_single_entry_scan_cli_smoke.jl`

- [ ] **Step 1: Write failing wrapper-equivalence smoke test**

Test that legacy script invocation and new subcommand route reach equivalent execution path (at least help/arg parse/exit behavior).

- [ ] **Step 2: Run tests to confirm failure**

Run targeted integration test file.
Expected: FAIL before adaptation.

- [ ] **Step 3: Refactor script into callable execution function**

Expose a callable `run_gap_transport_scan_main(args...)`-style function while preserving existing behavior.

- [ ] **Step 4: Keep legacy wrapper with deprecation line**

Ensure wrapper:
- forwards all args
- preserves downstream non-zero exit
- preserves stderr output

- [ ] **Step 5: Re-run tests**

Expected: PASS for this route.

- [ ] **Step 6: Sample commit-message style before commit**

Run: `git log -10 --oneline`
Use an existing recent prefix pattern for the commit message.

- [ ] **Step 7: Commit adaptation**

```bash
git add scripts/relaxtime/run_gap_transport_scan.jl tests/integration/relaxtime/test_single_entry_scan_cli_smoke.jl
git commit -m "refactor(relaxtime): route gap-transport script through unified entrypoint"
```

### Task 4: Adapt `scan_relaxation_times_vs_T.jl`

**Files:**
- Modify: `scripts/relaxtime/scan_relaxation_times_vs_T.jl`
- Modify: `tests/integration/relaxtime/test_single_entry_scan_cli_smoke.jl`

- [ ] **Step 1: Write failing wrapper-equivalence smoke test for `tau-vs-t`**
- [ ] **Step 2: Run test to confirm failure**
- [ ] **Step 3: Refactor to callable execution function**
- [ ] **Step 4: Keep legacy forwarding wrapper with deprecation line**

Ensure and test:
- full `ARGS` forwarding
- non-zero exit code passthrough
- stderr passthrough
- [ ] **Step 5: Re-run tests to pass**
- [ ] **Step 6: Sample commit-message style before commit**

Run: `git log -10 --oneline`
Use an existing recent prefix pattern for the commit message.

- [ ] **Step 7: Commit adaptation**

```bash
git add scripts/relaxtime/scan_relaxation_times_vs_T.jl tests/integration/relaxtime/test_single_entry_scan_cli_smoke.jl
git commit -m "refactor(relaxtime): route tau-vs-t script through unified entrypoint"
```

### Task 5: Adapt `run_manual_relaxation_scan_workflow.jl`

**Files:**
- Modify: `scripts/relaxtime/run_manual_relaxation_scan_workflow.jl`
- Modify: `tests/integration/relaxtime/test_single_entry_scan_cli_smoke.jl`

- [ ] **Step 1: Write failing wrapper-equivalence smoke test for `manual-workflow`**
- [ ] **Step 2: Run test to confirm failure**
- [ ] **Step 3: Refactor to callable execution function**
- [ ] **Step 4: Keep legacy forwarding wrapper with deprecation line**

Ensure and test:
- full `ARGS` forwarding
- non-zero exit code passthrough
- stderr passthrough
- [ ] **Step 5: Re-run tests to pass**
- [ ] **Step 6: Sample commit-message style before commit**

Run: `git log -10 --oneline`
Use an existing recent prefix pattern for the commit message.

- [ ] **Step 7: Commit adaptation**

```bash
git add scripts/relaxtime/run_manual_relaxation_scan_workflow.jl tests/integration/relaxtime/test_single_entry_scan_cli_smoke.jl
git commit -m "refactor(relaxtime): route manual workflow through unified entrypoint"
```

## Chunk 4: Workflow Migration + Final Verification

### Task 6: Switch GitHub manual workflow invocation to new entrypoint

**Files:**
- Modify: `.github/workflows/relaxtime-manual-scans.yml`

- [ ] **Step 1: Write failing mapping test/checklist in integration test file**

Assert that workflow input semantics map to new CLI flags without changing defaults.

Include explicit parse-based checks for:
- `integration_mode` -> `--mode`
- `tau_*` inputs -> matching CLI flags
- `compute_bulk` boolean semantics (`false` implies `--no-bulk`)
- `overwrite` boolean semantics (`true` implies `--overwrite`)

- [ ] **Step 2: Run tests/checks to capture failure pre-change**

Run targeted test file (and any lightweight YAML validation available in repo flow).

Run: `julia --project=. -e 'include("tests/integration/relaxtime/test_single_entry_scan_cli_smoke.jl")'`
Expected: FAIL on workflow-mapping assertions before migration.

- [ ] **Step 3: Update workflow run step to `run_scan.jl manual-workflow`**

Keep input keys stable. Update only invocation target and flag wiring.

- [ ] **Step 4: Re-run tests/checks**

Expected: PASS.

Run: `julia --project=. -e 'include("tests/integration/relaxtime/test_single_entry_scan_cli_smoke.jl")'`
Expected: PASS including workflow-mapping assertions.

- [ ] **Step 5: Sample commit-message style before commit**

Run: `git log -10 --oneline`
Use an existing recent prefix pattern for the commit message.

- [ ] **Step 6: Commit workflow migration**

```bash
git add .github/workflows/relaxtime-manual-scans.yml tests/integration/relaxtime/test_single_entry_scan_cli_smoke.jl
git commit -m "ci(relaxtime): invoke manual scans via unified scan entrypoint"
```

### Task 7: Final precision and compatibility verification sweep

**Files:**
- No mandatory new file edits

- [ ] **Step 1: Verify scope boundary**

Run:
- `git status --short`
- `git diff --name-only --cached`
- `git diff --name-only`

Expected: only planned scripts/workflow/tests files changed.

- [ ] **Step 2: Run target integration test file**

Run: `julia --project=. -e 'include("tests/integration/relaxtime/test_single_entry_scan_cli_smoke.jl")'`
Expected: PASS.

- [ ] **Step 3: Run relevant smoke selectors**

Run repo-supported minimal selectors that cover modified pathways.

- [ ] **Step 4: Build DoD evidence checklist for PR**

Must include evidence for:
- single entrypoint with three subcommands
- alias conflict policy enforcement
- wrapper compatibility semantics
- workflow mapping migration
- no precision-default regression

## Notes For Implementers

- Do not reduce precision-sensitive defaults during unification.
- Treat per-subcommand default differences as contractual unless explicitly changed with rationale.
- Keep wrappers thin; avoid introducing logic forks between legacy and new path.
- Avoid broad refactors outside scan CLI/workflow scope.
