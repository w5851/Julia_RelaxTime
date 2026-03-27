# Relaxtime Config Platform Strong Governance Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a single authoritative relaxtime orchestration path with TOML+CLI config governance, anti-drift gates, and parameter no-intercept guarantees.

**Architecture:** Introduce one public orchestrator entrypoint that normalizes config once, dispatches to internal authoritative scripts, and emits audit artifacts (`effective_config.json`, `consumption_report.json`, `run_manifest.json`). Enforce strong governance via script classification + static gate + workflow-vs-direct consistency regression. Add schema/alias and pass-through test coverage to ensure parameters are truly consumed.

**Tech Stack:** Julia 1.10+, include-driven scripts, existing `scripts/relaxtime/*`, `scripts/dev/*` governance checks, `tests/unit` + `tests/integration` + `tests/regression`.

---

## File Structure (Planned)

- Create: `scripts/relaxtime/run_relaxtime_orchestrator.jl`
- Create: `scripts/relaxtime/config/WorkflowConfig.jl`
- Create: `scripts/relaxtime/config/WorkflowConfigSchema.jl`
- Create: `scripts/relaxtime/config/WorkflowConfigAudit.jl`
- Create: `scripts/relaxtime/workflow/classification_manifest.toml`
- Create: `scripts/dev/check_relaxtime_script_governance.jl`
- Create: `scripts/relaxtime/run_cross_section_orchestrated_scan.jl`
- Create: `tests/unit/relaxtime/test_workflow_config_schema.jl`
- Create: `tests/unit/relaxtime/test_workflow_config_audit.jl`
- Create: `tests/unit/relaxtime/test_script_governance_checker.jl`
- Create: `tests/integration/relaxtime/test_orchestrator_smoke.jl`
- Create: `tests/integration/relaxtime/test_cross_section_orchestrated_smoke.jl`
- Create: `tests/regression/relaxtime/test_workflow_vs_direct_consistency.jl`
- Create: `config/workflows/relaxtime/default.toml`
- Create: `config/workflows/relaxtime/profiles/muB0_transport_xi.toml`
- Create: `config/workflows/relaxtime/profiles/muB0_cross_section_xi.toml`
- Create: `config/workflows/relaxtime/schema/relaxtime_workflow_schema_v1.toml`
- Create: `config/workflows/relaxtime/schema/aliases_v1.toml`
- Modify: `scripts/relaxtime/run_gap_transport_scan.jl`
- Modify: `scripts/relaxtime/scan_cross_section_vs_s_by_process.jl`
- Modify: `scripts/relaxtime/run_manual_relaxation_scan_workflow.jl`
- Modify: `scripts/relaxtime/run_outputs_2026_02_05.jl`

## Chunk 1: Config Foundation + Public Entrypoint

### Task 1: Build schema/alias/default config assets

**Files:**
- Create: `config/workflows/relaxtime/schema/relaxtime_workflow_schema_v1.toml`
- Create: `config/workflows/relaxtime/schema/aliases_v1.toml`
- Create: `config/workflows/relaxtime/default.toml`
- Create: `config/workflows/relaxtime/profiles/muB0_transport_xi.toml`
- Create: `config/workflows/relaxtime/profiles/muB0_cross_section_xi.toml`
- Test: `tests/unit/relaxtime/test_workflow_config_schema.jl`

- [ ] **Step 1: Write failing schema test cases**
  - Cover required keys, default xi list, default T list, default 4-process set, alias normalization table.

- [ ] **Step 2: Run failing tests**
  - Run: `julia --project=. -e 'include("tests/unit/relaxtime/test_workflow_config_schema.jl")'`
  - Expected: FAIL for missing schema/config artifacts.

- [ ] **Step 3: Implement schema/config artifacts minimally**
  - Add schema v1 and profile TOML files with frozen defaults from spec.

- [ ] **Step 4: Re-run test to pass**
  - Run: `julia --project=. -e 'include("tests/unit/relaxtime/test_workflow_config_schema.jl")'`
  - Expected: PASS.

- [ ] **Step 5: Commit**
  - `git add config/workflows/relaxtime tests/unit/relaxtime/test_workflow_config_schema.jl`
  - `git commit -m "feat: add relaxtime workflow schema and default profiles"`

### Task 2: Implement normalized config loader pipeline

**Files:**
- Create: `scripts/relaxtime/config/WorkflowConfig.jl`
- Create: `scripts/relaxtime/config/WorkflowConfigSchema.jl`
- Modify: `scripts/relaxtime/run_relaxtime_orchestrator.jl`
- Test: `tests/unit/relaxtime/test_workflow_config_schema.jl`

- [ ] **Step 1: Add failing tests for pipeline order**
  - Assert fixed order: raw -> alias normalize -> source validate -> merge -> effective validate.

- [ ] **Step 2: Run tests to confirm failure**
  - Run same unit command for schema test.
  - Expected: FAIL on missing loader and pipeline behavior.

- [ ] **Step 3: Implement minimal loader + normalization**
  - Support `--config`, optional CLI overrides, and precedence `default < TOML < CLI`.

- [ ] **Step 4: Re-run tests**
  - Expected: PASS for precedence and alias normalization.

- [ ] **Step 5: Commit**
  - `git add scripts/relaxtime/config scripts/relaxtime/run_relaxtime_orchestrator.jl tests/unit/relaxtime/test_workflow_config_schema.jl`
  - `git commit -m "feat: add relaxtime config normalization pipeline"`

### Task 3: Add public orchestrator command surface

**Files:**
- Create: `scripts/relaxtime/run_relaxtime_orchestrator.jl`
- Test: `tests/integration/relaxtime/test_orchestrator_smoke.jl`

- [ ] **Step 1: Write failing integration smoke for subcommands**
  - `transport` and `cross-section` command parse + dispatch.

- [ ] **Step 2: Run failing smoke**
  - Run: `julia --project=. -e 'include("tests/integration/relaxtime/test_orchestrator_smoke.jl")'`
  - Expected: FAIL (entrypoint missing/incomplete).

- [ ] **Step 3: Implement minimal orchestrator dispatch**
  - Support options: `--config`, `--outdir`, `--resume`, `--overwrite`, `--fail-on-fallback`, cross-section `--processes`.

- [ ] **Step 4: Re-run smoke**
  - Expected: PASS.

- [ ] **Step 5: Commit**
  - `git add scripts/relaxtime/run_relaxtime_orchestrator.jl tests/integration/relaxtime/test_orchestrator_smoke.jl`
  - `git commit -m "feat: add public relaxtime orchestrator entrypoint"`

## Chunk 2: Deep Consumption + Cross-Section Orchestration

### Task 4: Implement consumption audit and strict mode

**Files:**
- Create: `scripts/relaxtime/config/WorkflowConfigAudit.jl`
- Modify: `scripts/relaxtime/run_relaxtime_orchestrator.jl`
- Test: `tests/unit/relaxtime/test_workflow_config_audit.jl`

- [ ] **Step 1: Write failing audit tests**
  - unknown key rejection, unconsumed key rejection, consumed/unused/overridden report structure.

- [ ] **Step 2: Run failing tests**
  - Run: `julia --project=. -e 'include("tests/unit/relaxtime/test_workflow_config_audit.jl")'`
  - Expected: FAIL.

- [ ] **Step 3: Implement audit report emission**
  - Emit `effective_config.json`, `consumption_report.json`, `run_manifest.json`.

- [ ] **Step 4: Re-run tests**
  - Expected: PASS.

- [ ] **Step 5: Commit**
  - `git add scripts/relaxtime/config/WorkflowConfigAudit.jl scripts/relaxtime/run_relaxtime_orchestrator.jl tests/unit/relaxtime/test_workflow_config_audit.jl`
  - `git commit -m "feat: add strict config audit and consumption reports"`

### Task 5: Wire fallback observability and fail-on-fallback

**Files:**
- Modify: `scripts/relaxtime/run_gap_transport_scan.jl`
- Modify: `scripts/relaxtime/run_relaxtime_orchestrator.jl`
- Test: `tests/integration/relaxtime/test_orchestrator_smoke.jl`

- [ ] **Step 1: Add failing fallback test**
  - simulate/force fallback event and assert manifest `fallback_events` record.

- [ ] **Step 2: Run failing test**
  - Run integration smoke target for fallback scenario.

- [ ] **Step 3: Implement fallback event collection**
  - Capture bulk fallback and expose `--fail-on-fallback` behavior.

- [ ] **Step 4: Re-run tests**
  - Expected: PASS for both fallback-tolerant and fail-on-fallback mode.

- [ ] **Step 5: Commit**
  - `git add scripts/relaxtime/run_gap_transport_scan.jl scripts/relaxtime/run_relaxtime_orchestrator.jl tests/integration/relaxtime/test_orchestrator_smoke.jl`
  - `git commit -m "feat: add fallback event governance and fail-on-fallback"`

### Task 6: Add cross-section multi-T multi-xi orchestrated runner

**Files:**
- Create: `scripts/relaxtime/run_cross_section_orchestrated_scan.jl`
- Modify: `scripts/relaxtime/scan_cross_section_vs_s_by_process.jl`
- Test: `tests/integration/relaxtime/test_cross_section_orchestrated_smoke.jl`

- [ ] **Step 1: Write failing integration test for row-count formula**
  - Assert `N_rows = N_T * N_xi * N_process * N_sqrt_s` for normalized effective config.

- [ ] **Step 2: Run failing integration test**
  - Run: `julia --project=. -e 'include("tests/integration/relaxtime/test_cross_section_orchestrated_smoke.jl")'`
  - Expected: FAIL.

- [ ] **Step 3: Implement orchestrated batch scan**
  - Loop `T_list x xi_list x processes`, normalize process aliases, merge outputs, emit deterministic columns and run_id.

- [ ] **Step 4: Re-run integration test**
  - Expected: PASS.

- [ ] **Step 5: Commit**
  - `git add scripts/relaxtime/run_cross_section_orchestrated_scan.jl scripts/relaxtime/scan_cross_section_vs_s_by_process.jl tests/integration/relaxtime/test_cross_section_orchestrated_smoke.jl`
  - `git commit -m "feat: add orchestrated cross-section multi-T multi-xi pipeline"`

### Task 7: Add plotting contract automation

**Files:**
- Modify: `scripts/relaxtime/run_manual_relaxation_scan_workflow.jl`
- Modify: `scripts/relaxtime/run_outputs_2026_02_05.jl`
- Modify: `scripts/relaxtime/run_relaxtime_orchestrator.jl`
- Test: `tests/integration/relaxtime/test_cross_section_orchestrated_smoke.jl`

- [ ] **Step 1: Add failing tests for deterministic filenames/counts**
  - Validate patterns:
    - `transport__{metric}__muB{muB_tag}.png`
    - `xsec__T{T_tag}__{process}.png`

- [ ] **Step 2: Run failing tests**
  - Expected: FAIL on current naming/paths.

- [ ] **Step 3: Implement deterministic naming + figure-count checks**
  - Enforce 11 transport figures and cross-section `N_T * N_process` under defaults.

- [ ] **Step 4: Re-run tests**
  - Expected: PASS.

- [ ] **Step 5: Commit**
  - `git add scripts/relaxtime/run_manual_relaxation_scan_workflow.jl scripts/relaxtime/run_outputs_2026_02_05.jl scripts/relaxtime/run_relaxtime_orchestrator.jl tests/integration/relaxtime/test_cross_section_orchestrated_smoke.jl`
  - `git commit -m "feat: enforce deterministic plotting contracts in orchestrated workflows"`

## Chunk 3: Governance Gates + Consistency Regression + Closure

### Task 8: Add script classification manifest and static governance checker

**Files:**
- Create: `scripts/relaxtime/workflow/classification_manifest.toml`
- Create: `scripts/dev/check_relaxtime_script_governance.jl`
- Test: `tests/unit/relaxtime/test_script_governance_checker.jl`

- [ ] **Step 1: Write failing governance tests**
  - Ensure `public-authoritative` scripts cannot directly include `src/relaxtime/*`.

- [ ] **Step 2: Run failing tests**
  - Run: `julia --project=. -e 'include("tests/unit/relaxtime/test_script_governance_checker.jl")'`
  - Expected: FAIL.

- [ ] **Step 3: Implement checker + manifest parsing**
  - Fail on policy violation; support explicit allowlist for migration-only internal scripts.

- [ ] **Step 4: Re-run tests**
  - Expected: PASS.

- [ ] **Step 5: Commit**
  - `git add scripts/relaxtime/workflow/classification_manifest.toml scripts/dev/check_relaxtime_script_governance.jl tests/unit/relaxtime/test_script_governance_checker.jl`
  - `git commit -m "ci: add relaxtime script governance checker and classification manifest"`

### Task 9: Add workflow-vs-direct mini consistency regression

**Files:**
- Create: `tests/regression/relaxtime/test_workflow_vs_direct_consistency.jl`
- Modify: `scripts/relaxtime/run_transport_fixedpoint_regression.jl`
- Modify: `scripts/relaxtime/run_total_cross_section_fixedpoint_regression.jl`

- [ ] **Step 1: Write failing regression test on fixed mini-set**
  - Compare authoritative path and selected direct path outputs under tolerances.

- [ ] **Step 2: Run failing regression test**
  - Run: `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
  - Expected: FAIL initially.

- [ ] **Step 3: Implement minimal harmonization hooks**
  - Ensure both paths use aligned config profiles and observable provenance.

- [ ] **Step 4: Re-run regression**
  - Expected: PASS on mini consistency set.

- [ ] **Step 5: Commit**
  - `git add tests/regression/relaxtime/test_workflow_vs_direct_consistency.jl scripts/relaxtime/run_transport_fixedpoint_regression.jl scripts/relaxtime/run_total_cross_section_fixedpoint_regression.jl`
  - `git commit -m "test: add workflow-versus-direct consistency regression guard"`

### Task 10: Final verification + docs sync

**Files:**
- Modify: `docs/superpowers/specs/2026-03-27-relaxtime-demand-oriented-scan-config-platform-design.md` (if needed for implementation deltas)
- Optional Modify: `docs/dev/active/2026-03-01_PNJL可选功能盘点与优先级任务单.md` (only if governance status is updated)

- [ ] **Step 1: Run full targeted verification set**
  - `julia --project=. -e 'include("tests/unit/relaxtime/test_workflow_config_schema.jl")'`
  - `julia --project=. -e 'include("tests/unit/relaxtime/test_workflow_config_audit.jl")'`
  - `julia --project=. -e 'include("tests/unit/relaxtime/test_script_governance_checker.jl")'`
  - `julia --project=. -e 'include("tests/integration/relaxtime/test_orchestrator_smoke.jl")'`
  - `julia --project=. -e 'include("tests/integration/relaxtime/test_cross_section_orchestrated_smoke.jl")'`
  - `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`

- [ ] **Step 2: Run governance script directly**
  - `julia --project=. scripts/dev/check_relaxtime_script_governance.jl`
  - Expected: PASS.

- [ ] **Step 3: Verify official profile artifact completeness**
  - Confirm each official run directory includes `run_manifest.json`, `effective_config.json`, `consumption_report.json`.

- [ ] **Step 4: Update docs deltas if needed**
  - Keep spec and implementation behavior aligned.

- [ ] **Step 5: Commit**
  - `git add docs/superpowers/specs/2026-03-27-relaxtime-demand-oriented-scan-config-platform-design.md docs/dev/active/2026-03-01_PNJL可选功能盘点与优先级任务单.md`
  - `git commit -m "docs: align relaxtime governance spec with implementation evidence"`

## Execution Notes

- Prefer one task per branch commit to keep rollback/review cheap.
- For heavy scans, keep smoke data size minimal and deterministic.
- Any new default that affects numerics must include rationale in commit body and regression evidence.

## Definition of Done

- Official output path is only `run_relaxtime_orchestrator.jl`.
- Schema v1 + alias map + profile defaults exist and are validated.
- Strong governance gate is active and tested.
- AC-1/2/3/4/5 from spec are all test-covered and passing.
- Parameter no-intercept guarantees are enforced through strict mode + consumption report + probe tests.
