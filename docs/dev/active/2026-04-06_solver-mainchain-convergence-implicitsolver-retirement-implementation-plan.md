# Solver Mainchain Convergence and ImplicitSolver Retirement Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove non-mainchain solver branches and fully retire `ImplicitSolver`, so runtime solver behavior only follows the PR56 mainchain (`Models.solve/solve_constraint` -> `ProblemSpec` -> governed selection).

**Architecture:** Keep all solve paths inside `src/models/solver/Solver.jl` + `src/models/solver/ProblemSpec.jl` + `src/models/constraint_solver.jl`. Remove legacy fallback calls to `Main.Models.ImplicitSolver.solve(...)`, then remove module inclusion/export surface and delete `src/models/solver/ImplicitSolver.jl`. Replace legacy-comparison tests with mainchain contract/regression tests at fixed points.

**Tech Stack:** Julia, NLsolve, ForwardDiff, existing `Models` solver stack, unit + integration + regression suites, docs governance scripts.

---

## File Structure Map (planned changes)

- Modify: `src/models/constraint_solver.jl`
  - Remove 4 legacy fallback branches that call `Main.Models.ImplicitSolver.solve(...)`.
  - Keep failure diagnostics in mainchain shape (`converged`, `failed_constraints`, `legacy_fallback_used=false`).
- Modify: `src/models/solver/Solver.jl`
  - Remove compatibility bridge branches that exist only for non-mainchain routing.
  - Keep one canonical path for non-`FixedMu` solve/solve_multi through `solve_constraint` mainchain.
- Modify: `src/models/Models.jl`
  - Remove `include(joinpath(@__DIR__, "solver", "ImplicitSolver.jl"))`.
  - Clean exports/compat aliases if any still imply `ImplicitSolver` surface.
- Delete: `src/models/solver/ImplicitSolver.jl`
- Modify tests:
  - `tests/unit/models/test_solver_dimension_agnostic.jl`
  - `tests/unit/pnjl/test_solver_implicit.jl`
  - `tests/regression/models/test_problem_spec_fixedrho_parity_regression.jl`
  - (as needed) `tests/unit/models/test_solver.jl`, `tests/unit/models/test_problem_spec_contract.jl`
- Modify docs:
  - `docs/dev/active/` add/update one execution task doc for this retirement wave
  - `docs/api/models/solver/*` if exported API wording changed

---

## Chunk 1: Guardrails First (TDD + baseline evidence)

### Task 1: Add failing tests for “no ImplicitSolver runtime dependency”

**Files:**
- Modify: `tests/unit/models/test_solver_dimension_agnostic.jl`
- Modify: `tests/regression/models/test_problem_spec_fixedrho_parity_regression.jl`
- Modify: `tests/unit/pnjl/test_solver_implicit.jl`

- [x] **Step 1: Replace legacy-API assumptions with retirement contract tests**
  - Add assertions that `Models` does not expose `ImplicitSolver` runtime path.
  - Replace parity tests using `Models.ImplicitSolver.solve(...)` with fixed-point mainchain baseline checks.

- [x] **Step 2: Run targeted tests and confirm failure before implementation**
  - Run: `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_solver_dimension_agnostic.jl,pnjl/test_solver_implicit.jl"; include("tests/unit/runtests.jl")'`
  - Run: `julia --project=. -e 'ENV["REGRESSION_FILES"]="models/test_problem_spec_fixedrho_parity_regression.jl"; include("tests/regression/runtests.jl")'`
  - Result: unit guardrails failed as expected before implementation; regression guardrails passed.

- [ ] **Step 3: Commit test-only guardrails**
  - `git commit -m "test(models/solver): add retirement guardrails for ImplicitSolver removal"`

---

## Chunk 2: Remove Legacy Fallback Runtime Branches

### Task 2: Drop `constraint_solver` fallback calls to `Main.Models.ImplicitSolver`

**Files:**
- Modify: `src/models/constraint_solver.jl`
- Test: `tests/unit/models/test_problem_spec_contract.jl`
- Test: `tests/unit/models/test_solver.jl`

- [x] **Step 1: Remove fallback block in `_solve_constraint_fixedmu`**
  - Remove `if allow_legacy_fallback && !s.converged ... Main.Models.ImplicitSolver.solve(...)` branch.
  - Keep unified return payload with explicit mainchain diagnostics.

- [x] **Step 2: Remove fallback blocks in non-`FixedMu` solvers**
  - `_solve_constraint_fixedentropy`
  - `_solve_constraint_fixedsigma`
  - `_solve_constraint_fixedasymrho`
  - Ensure `legacy_fallback_used` remains deterministic (`false`) or field removed consistently with callers.

- [x] **Step 3: Run focused unit tests**
  - Run: `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_problem_spec_contract.jl,models/test_solver.jl"; include("tests/unit/runtests.jl")'`
  - Result: PASS after updating solver parity assertions to mainchain-consistent contract checks.

- [ ] **Step 4: Commit runtime fallback cleanup**
  - `git commit -m "refactor(models/solver): remove constraint_solver legacy fallback branches"`

---

## Chunk 3: Converge Solver Routing and Delete Module

### Task 3: Remove bridge-only paths and retire module include

**Files:**
- Modify: `src/models/solver/Solver.jl`
- Modify: `src/models/Models.jl`
- Delete: `src/models/solver/ImplicitSolver.jl`
- Test: `tests/unit/models/test_solver_dimension_agnostic.jl`

- [x] **Step 1: Simplify non-mainchain bridge logic in `Solver.jl`**
  - Remove compatibility branches that gate between old/non-old route.
  - Keep single route for solve/solve_multi into `solve_constraint` governed path.

- [x] **Step 2: Remove module include and references in `Models.jl`**
  - Delete `include(... "ImplicitSolver.jl")`.
  - Clean any direct exports depending on `Models.ImplicitSolver` symbol.

- [x] **Step 3: Delete retired source file**
  - Delete: `src/models/solver/ImplicitSolver.jl`

- [x] **Step 4: Run targeted tests for symbol surface and routing**
  - Run: `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_solver_dimension_agnostic.jl,models/test_solver.jl,pnjl/test_solver_implicit.jl"; include("tests/unit/runtests.jl")'`
  - Run: `julia --project=. -e 'ENV["REGRESSION_FILES"]="models/test_problem_spec_fixedrho_parity_regression.jl"; include("tests/regression/runtests.jl")'`
  - Result: PASS and no runtime `Models.ImplicitSolver` references under `src/`.

- [ ] **Step 5: Commit solver routing convergence + module retirement**
  - `git commit -m "refactor(models/solver): retire ImplicitSolver and converge to mainchain routing"`

---

## Chunk 4: Full Verification + Docs/Regression Closure

### Task 4: Validate full gate and update docs status

**Files:**
- Modify: `docs/dev/active/` retirement task doc (create/update)
- Modify: `docs/api/models/solver/*` (if export wording changed)

- [x] **Step 1: Run smoke suites**
  - `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'` -> PASS (781/781)
  - `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'` -> PASS (455/455)
  - `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'` -> PASS (468 pass, 1 expected broken)

- [x] **Step 2: Run governance checks**
  - `julia --project=. scripts/dev/check_docs_consistency.jl` -> PASS
  - `julia --project=. scripts/dev/check_active_docs_governance.jl` -> PASS after filename normalization to `YYYY-MM-DD_*.md`
  - `julia --project=. scripts/dev/check_legacy_solver_switch_leakage.jl` -> PASS

- [x] **Step 3: Verify no code references remain**
  - Search `src`: no `ImplicitSolver` references.
  - Search `tests`: only retirement guardrail assertions and historical test title text remain.

- [x] **Step 4: Update task doc and risk notes**
  - Record verification outputs and known trade-offs.
  - Mark completed checklist items and next follow-up if any.

- [ ] **Step 5: Final commit**
  - `git commit -m "docs(dev): close ImplicitSolver retirement with full verification evidence"`

---

## Rollback Strategy

- If non-`FixedMu` convergence drops on smoke/regression, immediately revert the specific mode cleanup commit and retain only passing chunks.
- Do not reintroduce `ImplicitSolver` as broad fallback; use mode-scoped mitigation with explicit diagnostics.

## Done Criteria (DoD)

- [x] `src/models/solver/ImplicitSolver.jl` removed.
- [x] `src/models/Models.jl` no longer includes `ImplicitSolver`.
- [x] `constraint_solver` has no `Main.Models.ImplicitSolver.solve(...)` fallback calls.
- [x] Solver runtime paths are uniformly mainchain (`solve_constraint`/`ProblemSpec` governed flow).
- [x] Unit + integration + regression smoke all pass.
- [x] Governance scripts pass and docs status is updated.
