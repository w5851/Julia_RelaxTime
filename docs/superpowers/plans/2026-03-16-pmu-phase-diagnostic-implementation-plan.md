# P-mu Phase Diagnostic Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a compare-only `P-mu / Omega-mu` branch-competition diagnostic for the CEP neighborhood and produce side-by-side artifacts against Maxwell results.

**Architecture:** Implement the diagnostic as a small group of focused phase-topic units: one file for seed and branch-solve contracts, one for result summarization and Maxwell comparison, and one for artifact writing/public orchestration. The implementation must reuse the existing solver/thermo path that already returns `omega`, `pressure`, `rho_norm`, and `residual_norm`, while leaving the current Maxwell-based production CEP path unchanged.

**Tech Stack:** Julia 1.10+, existing `Models` solver/thermo stack, CSV/JSON artifact generation, `Test`.

---

## File Map

- Create: `src/models/phase/PMPhaseTypes.jl`
  - result structs, status enums encoded as symbols/strings, and schema constants
- Create: `src/models/phase/PMPhaseSeeds.jl`
  - `seed_pair` schema, default seed derivation, branch-seed propagation helpers
- Create: `src/models/phase/PMPhaseDiagnostic.jl`
  - fixed-`T` / fixed-`mu` solve loop, continuity checks, temperature-level summaries
- Create: `src/models/phase/PMPhaseArtifacts.jl`
  - CSV/JSON writing for branch scan and compare tables
- Modify: `src/models/Models.jl`
  - include/export new public entrypoint(s)
- Modify: `src/models/entrypoints.jl`
  - re-export new entrypoint if this file remains the unified public surface
- Create: `scripts/pnjl/diagnose_pm_phase.jl`
  - thin CLI wrapper for compare-only diagnostics
- Create: `tests/unit/models/test_pm_phase_diagnostic.jl`
  - unit coverage for seed handling, status normalization, transition interpolation, endpoint extraction, comparison summaries
- Create: `tests/integration/models/test_pm_phase_diagnostic_smoke.jl`
  - smoke coverage for artifact generation and end-to-end bounded runs
- Modify: `docs/api/models/phase/Overview.md`
  - add the compare-only diagnostic entrypoint
- Create: `docs/api/models/phase/PMPhaseDiagnostic.md`
  - focused usage/API doc

## Existing Sources To Reuse Explicitly

- Solver/thermo source: `src/models/constraint_solver.jl`
  - existing fixed-state results already provide `omega`, `pressure`, `rho_norm`, `residual_norm`
- Maxwell reference source: `src/models/phase/PhaseCore.jl` and current phase-path helpers in `src/models/phase/CEPDetector.jl` / `src/models/phase/ProductionPhasePipeline.jl`
- Seed derivation source candidates:
  - low-`mu` / high-`mu` equilibria from the same solve stack exposed through `src/models/constraint_solver.jl`
  - if a direct helper is needed, add a thin local wrapper in `src/models/phase/PMPhaseSeeds.jl` instead of spreading seed logic into the diagnostic loop
- Artifact conventions: current phase artifact style in `src/models/phase/PhasePipeline.jl`

## Frozen v1 Defaults

- `mu_grid_default = 289.5:0.1:292.5 MeV`
- `mu_refine_step = 0.01 MeV`
- `comparison_mu_tol = 0.05 MeV`
- `residual_accept_tol = 1e-6`
- `continuity_x_tol = 0.25`
- `continuity_rho_tol = 0.15`
- `continuity_mode = :branch_local`
- `fallback_mode = :none`

## Frozen v1 Artifact Contract

### Branch status enum

- `accepted`
- `nonconverged`
- `branch_jump`
- `invalid_thermo`

### Seed source enum

- `seed0`
- `previous_same_branch`
- `manual_override`

### Endpoint cause enum

- `physical_loss_candidate`
- `nonconvergence`
- `branch_jump`
- `out_of_grid`

### Comparison status enum

- `both`
- `pm_only`
- `maxwell_only`
- `neither`

### Required `pm_branch_scan.csv` fields

- `T_MeV`
- `mu_MeV`
- `branch`
- `branch_status`
- `seed_source`
- `continuity_ok`
- `converged`
- `residual_norm`
- `rho_norm`
- `pressure`
- `omega`
- `endpoint_cause`

### Required `pm_phase_summary.json` fields per temperature

- `T_MeV`
- `mu_transition_pm_MeV`
- `hadron_endpoint_mu_MeV`
- `quark_endpoint_mu_MeV`
- `bistable_window_width_MeV`
- `comparison_status`
- `comparison_mu_tol_MeV`
- `residual_accept_tol`
- `continuity_x_tol`
- `continuity_rho_tol`

### Required `pm_vs_maxwell.csv` fields

- `T_MeV`
- `mu_transition_pm_MeV`
- `mu_transition_maxwell_MeV`
- `delta_mu_pm_minus_maxwell_MeV`
- `comparison_status`
- `branch_disappears_first`
- `hadron_endpoint_mu_MeV`
- `quark_endpoint_mu_MeV`
- `bistable_window_width_MeV`

## Chunk 1: Seed Contract And Public Type Surface

### Task 1: Add failing unit tests for seed and enum contracts

**Files:**
- Create: `tests/unit/models/test_pm_phase_diagnostic.jl`
- Create: `src/models/phase/PMPhaseTypes.jl`
- Create: `src/models/phase/PMPhaseSeeds.jl`

- [ ] **Step 1: Write the failing test for enum contracts**

Add assertions that the v1 enum/status values exist exactly as frozen above.

- [ ] **Step 2: Write the failing test for `seed_pair` schema**

Add a test that expects a normalized `seed_pair` object with:
- `hadron_seed0`
- `quark_seed0`
- `continuity_mode`
- `fallback_mode`

- [ ] **Step 3: Write the failing test for branch seed propagation**

Add a test that accepted branch results reuse `previous_same_branch`, while failed points do not overwrite the cached branch seed.

- [ ] **Step 4: Run tests to verify they fail**

Run: `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_pm_phase_diagnostic.jl"; include("tests/unit/runtests.jl")'`
Expected: FAIL because the new type/seed files do not exist.

- [ ] **Step 5: Create `PMPhaseTypes.jl` with minimal type surface**

Implement:
- schema constants
- tiny structs or normalized named-tuple constructors for:
  - branch rows
  - temperature summaries
  - comparison summaries

- [ ] **Step 6: Create `PMPhaseSeeds.jl` with minimal `seed_pair` normalization**

Implement helpers for:
- seed-pair normalization
- branch-local propagation
- seed-source tagging

- [ ] **Step 7: Wire exports in `src/models/Models.jl` and `src/models/entrypoints.jl`**

Add `include(...)` and `export ...` for the new public entrypoint/types that already exist at this stage.

- [ ] **Step 8: Run tests to verify they pass**

Run the same unit command.
Expected: PASS for enum and seed-contract tests.

- [ ] **Step 9: Commit**

```bash
git add tests/unit/models/test_pm_phase_diagnostic.jl src/models/phase/PMPhaseTypes.jl src/models/phase/PMPhaseSeeds.jl src/models/Models.jl src/models/entrypoints.jl
git commit -m "feat: add P-mu diagnostic seed and type contracts"
```

## Chunk 2: Pure Summary And Comparison Helpers

### Task 2: Add failing unit tests for transition and endpoint extraction

**Files:**
- Modify: `tests/unit/models/test_pm_phase_diagnostic.jl`
- Create: `src/models/phase/PMPhaseDiagnostic.jl`

- [ ] **Step 1: Write the failing test for `Delta Omega` zero crossing**

Use synthetic branch tables to assert linear interpolation of `mu_transition_pm_MeV`.

- [ ] **Step 2: Write the failing test for endpoint extraction**

Assert correct `hadron_endpoint_mu_MeV`, `quark_endpoint_mu_MeV`, and `bistable_window_width_MeV`.

- [ ] **Step 3: Write the failing test for `no_bistable_window`**

Assert that disjoint accepted branches report the correct no-window condition.

- [ ] **Step 4: Write the failing test for Maxwell comparison summary**

Assert generation of `delta_mu_pm_minus_maxwell_MeV` and `comparison_status`.

- [ ] **Step 5: Write the failing test for `Delta P` cross-check and branch-loss ordering**

Assert that:
- `Delta P` sign change is reported consistently with `Delta Omega`
- summary records which branch disappears first

- [ ] **Step 6: Run tests to verify they fail**

Run: `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_pm_phase_diagnostic.jl"; include("tests/unit/runtests.jl")'`
Expected: FAIL on missing helpers.

- [ ] **Step 7: Implement transition interpolation helper**

Add a small pure helper in `PMPhaseDiagnostic.jl` that computes `mu_transition_pm` from adjacent `Delta Omega` points.

- [ ] **Step 8: Implement endpoint and window extraction helpers**

Add pure helpers for branch overlap, endpoint detection, and no-window summaries.

- [ ] **Step 9: Implement `Delta P` cross-check and branch-loss ordering helper**

Add a helper that:
- compares `Delta Omega` and `Delta P` zero-cross behavior
- records `branch_disappears_first`

- [ ] **Step 10: Implement Maxwell comparison helper**

Add a pure helper that compares PM and Maxwell values using `comparison_mu_tol`.

- [ ] **Step 11: Run tests to verify they pass**

Run the same unit command.
Expected: PASS.

- [ ] **Step 12: Commit**

```bash
git add tests/unit/models/test_pm_phase_diagnostic.jl src/models/phase/PMPhaseDiagnostic.jl
git commit -m "feat: add P-mu transition and comparison summaries"
```

## Chunk 3: Branch Acceptance And Continuity Rules

### Task 3: Add failing unit tests for acceptance and continuity thresholds

**Files:**
- Modify: `tests/unit/models/test_pm_phase_diagnostic.jl`
- Modify: `src/models/phase/PMPhaseDiagnostic.jl`

- [ ] **Step 1: Write the failing test for acceptance rules**

Assert that a branch point is accepted only when:
- `converged=true`
- thermo values are finite
- `residual_norm <= residual_accept_tol`

- [ ] **Step 2: Write the failing test for continuity success**

Assert that a small `x_state` and `rho_norm` change gives `continuity_ok=true`.

- [ ] **Step 3: Write the failing test for branch jump**

Assert that large jumps produce `branch_status = branch_jump` and `endpoint_cause = branch_jump`.

- [ ] **Step 4: Run tests to verify they fail**

Run the unit command again.
Expected: FAIL on missing acceptance/continuity helpers.

- [ ] **Step 5: Implement acceptance helper with explicit residual tolerance**

Add helper(s) that normalize solver outputs into v1 branch statuses.

- [ ] **Step 6: Implement continuity helper with fixed defaults**

Add helper(s) using `continuity_x_tol` and `continuity_rho_tol`.

- [ ] **Step 7: Run tests to verify they pass**

Run the same unit command.
Expected: PASS.

- [ ] **Step 8: Commit**

```bash
git add tests/unit/models/test_pm_phase_diagnostic.jl src/models/phase/PMPhaseDiagnostic.jl
git commit -m "feat: add P-mu branch acceptance thresholds"
```

## Chunk 4: Real Dual-Branch Solve Loop

### Task 4: Add failing integration smoke coverage for real branch scans

**Files:**
- Create: `tests/integration/models/test_pm_phase_diagnostic_smoke.jl`
- Modify: `src/models/phase/PMPhaseSeeds.jl`
- Modify: `src/models/phase/PMPhaseDiagnostic.jl`

- [ ] **Step 1: Write the failing smoke test for one temperature**

Use:
- `T_values=[130.9]`
- a short `mu_grid`
- `xi=0.0`
- `solver_backend=:legacy`
- `p_num=24`, `t_num=8`
- temp output directory

Assert that the call returns structured rows for both branches.

- [ ] **Step 2: Write the failing smoke test for artifact file existence**

Assert existence of:
- `pm_branch_scan.csv`
- `pm_phase_summary.json`
- `pm_vs_maxwell.csv`

- [ ] **Step 3: Run the smoke test to verify it fails**

Run: `julia --project=. tests/integration/models/test_pm_phase_diagnostic_smoke.jl`
Expected: FAIL because the real solve loop is not implemented.

- [ ] **Step 4: Implement default hadron/quark seed derivation**

In `PMPhaseSeeds.jl`, add concrete v1 helpers that derive `hadron_seed0` and `quark_seed0` from nearby low-`mu` and high-`mu` equilibria at the same `T`.

- [ ] **Step 5: Implement single-branch solve wrapper**

In `PMPhaseDiagnostic.jl`, add a focused helper that runs one branch solve at one `(T, mu)` and returns normalized branch-row data.

- [ ] **Step 6: Implement dual-branch loop for one temperature**

In `PMPhaseDiagnostic.jl`, add the fixed-`mu` loop that solves both branches, applies acceptance/continuity rules, and accumulates rows.

- [ ] **Step 7: Implement artifact writing shim call sites**

Return enough data so later artifact serialization can be added without changing the solve loop shape.

- [ ] **Step 8: Run the smoke test to verify it passes**

Run the same integration command.
Expected: PASS on the bounded scenario.

- [ ] **Step 9: Commit**

```bash
git add tests/integration/models/test_pm_phase_diagnostic_smoke.jl src/models/phase/PMPhaseSeeds.jl src/models/phase/PMPhaseDiagnostic.jl
git commit -m "feat: add P-mu dual-branch solve loop"
```

## Chunk 5: Artifact Writers And Schema Verification

### Task 5: Add failing tests for required artifact fields

**Files:**
- Modify: `tests/unit/models/test_pm_phase_diagnostic.jl`
- Modify: `tests/integration/models/test_pm_phase_diagnostic_smoke.jl`
- Create: `src/models/phase/PMPhaseArtifacts.jl`

- [ ] **Step 1: Write the failing unit test for branch-row schema**

Assert that serialized branch rows contain every required field in the v1 contract.

- [ ] **Step 2: Write the failing integration test for summary schema**

Assert that `pm_phase_summary.json` includes threshold metadata and endpoint values.

- [ ] **Step 3: Write the failing integration test for compare schema**

Assert that `pm_vs_maxwell.csv` contains the required side-by-side comparison columns.

- [ ] **Step 4: Run tests to verify they fail**

Run:
- `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_pm_phase_diagnostic.jl"; include("tests/unit/runtests.jl")'`
- `julia --project=. tests/integration/models/test_pm_phase_diagnostic_smoke.jl`

Expected: FAIL because schema-complete artifact writers do not exist.

- [ ] **Step 5: Implement branch-scan CSV writer**

In `PMPhaseArtifacts.jl`, add serialization for `pm_branch_scan.csv` with the exact required columns.

- [ ] **Step 6: Implement summary JSON writer**

Serialize one summary object per temperature, including threshold metadata.

- [ ] **Step 7: Implement side-by-side comparison CSV writer**

Serialize `pm_vs_maxwell.csv` with exact required comparison columns.

- [ ] **Step 8: Wire artifact writing into the public entrypoint**

Make `analyze_pm_branch_competition(...)` call the new artifact helpers.

- [ ] **Step 9: Run tests to verify they pass**

Run the same unit and integration commands.
Expected: PASS.

- [ ] **Step 10: Commit**

```bash
git add tests/unit/models/test_pm_phase_diagnostic.jl tests/integration/models/test_pm_phase_diagnostic_smoke.jl src/models/phase/PMPhaseArtifacts.jl src/models/phase/PMPhaseDiagnostic.jl
git commit -m "feat: add P-mu diagnostic artifacts"
```

## Chunk 6: Maxwell Reference Wiring

### Task 6: Add failing tests for same-config Maxwell comparison

**Files:**
- Modify: `tests/unit/models/test_pm_phase_diagnostic.jl`
- Modify: `src/models/phase/PMPhaseDiagnostic.jl`

- [ ] **Step 1: Write the failing unit test for Maxwell reference lookup contract**

Assert that comparison requires the same `T`, `xi`, backend, `p_num`, and `t_num`.

- [ ] **Step 2: Write the failing unit test for comparison tolerance handling**

Assert classification into `both`, `pm_only`, `maxwell_only`, and `neither`.

- [ ] **Step 3: Run tests to verify they fail**

Run the unit command.
Expected: FAIL on missing same-config Maxwell wiring.

- [ ] **Step 4: Implement Maxwell reference getter using current phase logic**

Add a focused helper in `PMPhaseDiagnostic.jl` that derives Maxwell reference values from the same physical configuration rather than stale historical artifacts.

Use explicitly:
- `detect_s_shape(...)` from `src/models/phase/PhaseCore.jl`
- `maxwell_construction(...)` from `src/models/phase/PhaseCore.jl`

on branch-comparable fixed-`T` data generated under the same `xi`, backend, `p_num`, and `t_num`.

- [ ] **Step 5: Implement side-by-side comparison classification**

Use `comparison_mu_tol` and the frozen comparison-status enum.

- [ ] **Step 6: Run tests to verify they pass**

Run the unit command.
Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add tests/unit/models/test_pm_phase_diagnostic.jl src/models/phase/PMPhaseDiagnostic.jl
git commit -m "feat: align P-mu diagnostics with Maxwell references"
```

## Chunk 7: Script Entry And API Docs

### Task 7: Add failing script smoke coverage

**Files:**
- Create: `scripts/pnjl/diagnose_pm_phase.jl`
- Modify: `tests/integration/models/test_pm_phase_diagnostic_smoke.jl`
- Modify: `docs/api/models/phase/Overview.md`
- Create: `docs/api/models/phase/PMPhaseDiagnostic.md`

- [ ] **Step 1: Write the failing integration test for the script entrypoint**

Shell out to the thin script with one tiny CEP-neighborhood input.

- [ ] **Step 2: Run the integration test to verify it fails**

Run: `julia --project=. tests/integration/models/test_pm_phase_diagnostic_smoke.jl`
Expected: FAIL because the script does not exist yet.

- [ ] **Step 3: Implement the thin script wrapper**

Add argument parsing only for the bounded diagnostic use case and delegate to the public entrypoint.

- [ ] **Step 4: Document the compare-only diagnostic in API docs**

Add:
- concise overview mention in `docs/api/models/phase/Overview.md`
- focused doc in `docs/api/models/phase/PMPhaseDiagnostic.md`

- [ ] **Step 5: Run the integration test to verify it passes**

Run the same integration command.
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add scripts/pnjl/diagnose_pm_phase.jl tests/integration/models/test_pm_phase_diagnostic_smoke.jl docs/api/models/phase/Overview.md docs/api/models/phase/PMPhaseDiagnostic.md
git commit -m "docs: add P-mu phase diagnostic usage docs"
```

## Chunk 7.5: Local Mu Refinement Around Sign Change

### Task 7.5: Implement the v1 local refinement contract

**Files:**
- Modify: `tests/unit/models/test_pm_phase_diagnostic.jl`
- Modify: `tests/integration/models/test_pm_phase_diagnostic_smoke.jl`
- Modify: `src/models/phase/PMPhaseDiagnostic.jl`

- [ ] **Step 1: Write the failing unit test for refinement trigger**

Assert that when adjacent accepted points bracket a `Delta Omega` sign change, the code requests local refinement toward `mu_refine_step = 0.01 MeV`.

- [ ] **Step 2: Write the failing integration test for refined transition output**

Assert that the diagnostic output reports a refined transition value rather than only the coarse-grid bracket.

- [ ] **Step 3: Run tests to verify they fail**

Run the relevant unit and integration tests.
Expected: FAIL because local refinement is not implemented yet.

- [ ] **Step 4: Implement local sign-change refinement**

Add a small helper that refines only around accepted `Delta Omega` or `Delta P` sign-change brackets down to the frozen `mu_refine_step` target.

- [ ] **Step 5: Run tests to verify they pass**

Run the same unit and integration tests.
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add tests/unit/models/test_pm_phase_diagnostic.jl tests/integration/models/test_pm_phase_diagnostic_smoke.jl src/models/phase/PMPhaseDiagnostic.jl
git commit -m "feat: refine P-mu transition brackets locally"
```

## Chunk 8: Full Verification On The Agreed Scope

### Task 8: Verify the full slice end-to-end

**Files:**
- No new files required unless fixes are needed

- [ ] **Step 1: Run focused unit tests**

Run: `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_pm_phase_diagnostic.jl,models/test_phase_pipeline.jl,models/test_cep_detector.jl,models/test_production_phase_pipeline.jl"; include("tests/unit/runtests.jl")'`
Expected: all selected tests pass.

- [ ] **Step 2: Run integration smoke tests**

Run:
- `julia --project=. tests/integration/models/test_pm_phase_diagnostic_smoke.jl`
- `julia --project=. tests/integration/models/test_phase_core_algorithms_smoke.jl`

Expected: both pass.

- [ ] **Step 3: Run docs consistency**

Run: `julia --project=. scripts/dev/check_docs_consistency.jl`
Expected: `OK`

- [ ] **Step 4: Run one real compare-only diagnostic command**

Run a bounded real command for:
- `xi=0.0`
- `legacy+p24/t8`
- `T≈130.8:0.1:131.1`
- default `mu_grid_default`

Expected: `pm_branch_scan.csv`, `pm_phase_summary.json`, and `pm_vs_maxwell.csv` are generated and readable.

- [ ] **Step 5: Inspect artifact semantics manually**

Check that:
- threshold metadata is present
- endpoint causes are populated
- `comparison_status` values are sensible

- [ ] **Step 6: Commit final verification or fixes**

```bash
git add <modified files>
git commit -m "test: verify P-mu diagnostic compare-only slice"
```

## Notes For The Implementer

- Stay compare-only in v1. Do not alter the current production CEP selection path.
- Do not compare against historical artifacts with mismatched grids or configs.
- If branch identity is ambiguous, record it explicitly in artifacts rather than silently forcing continuity.
- Keep each new file narrow in responsibility so later promotion to a primary phase criterion does not require disentangling one monolithic module.

## Plan Review Status

This plan is based on the approved spec in `docs/superpowers/specs/2026-03-16-pmu-phase-diagnostic-design.md`. It freezes the v1 defaults, seed contract, artifact schema, and Maxwell reference rule that the reviewer asked to make explicit.

Plan complete and saved to `docs/superpowers/plans/2026-03-16-pmu-phase-diagnostic-implementation-plan.md`. Ready to execute?
