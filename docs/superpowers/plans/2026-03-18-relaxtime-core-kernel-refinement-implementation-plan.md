# Relaxtime Core Kernel Refinement Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Deliver a first behavior-preserving `relaxtime` core-kernel PR that stabilizes `AverageScatteringRate` and `RelaxationTime` with explicit contract tests, without absorbing script/workflow/output scope.

**Architecture:** Use test-first refinement in small, reviewable slices. First lock observable kernel contracts in unit tests, then perform file-internal layering cleanup in `src/relaxtime/AverageScatteringRate.jl` and `src/relaxtime/RelaxationTime.jl`, and finally run a regression gate only when fixture-facing semantics are touched. Keep public entrypoints stable and keep CSV/script/workflow paths out of active refactoring.

**Tech Stack:** Julia 1.10+, existing `src/relaxtime` modules, `Test`, repo unit/regression harness (`tests/unit/runtests.jl`, `tests/regression/runtests.jl`).

---

## File Map

- Modify: `src/relaxtime/AverageScatteringRate.jl`
  - tighten file-internal helper boundaries for cache/interpolation/grid/refinement logic
- Modify: `src/relaxtime/RelaxationTime.jl`
  - tighten alias/rate aggregation boundaries while keeping outward behavior stable
- Modify: `tests/unit/relaxtime/test_average_scattering_rate.jl`
  - lock interpolation/cache/layering/refinement contracts
- Modify: `tests/unit/relaxtime/test_relaxation_time.jl`
  - lock tau algebra, alias handling, and warn-and-clamp contracts
- Conditional modify: `tests/regression/relaxtime/test_tau_xi_probe_regression.jl`
  - only if fixture-facing semantics change

## Scope Guardrails (Must Keep)

- Do not actively refactor `scripts/relaxtime/*` in this plan.
- Do not actively refactor `scripts/analysis/*` in this plan.
- Do not regenerate full `data/outputs/*` in this plan.
- Do not actively refactor CSV loader helpers in `src/relaxtime/RelaxationTime.jl` unless a minimal mechanical edit is strictly needed to keep the file valid.
- Preserve public callable surface:
  - `average_scattering_rate`
  - `build_w0cdf_pchip_cache`
  - `precompute_cross_section!`
  - `number_density`
  - `compute_average_rates`
  - `relaxation_rates`
  - `relaxation_times`
  - `relaxation_rate_contribution_rows`

## Chunk 1: Contract Locking Tests First

### Task 1: Strengthen `AverageScatteringRate` contract tests

**Files:**
- Modify: `tests/unit/relaxtime/test_average_scattering_rate.jl`
- Test target: `src/relaxtime/AverageScatteringRate.jl`

- [ ] **Step 1: Add/adjust failing tests for explicit cache-range contract**

Ensure tests explicitly verify out-of-range `sigma(s)` lookup returns `0.0` and does not clamp to edge values.

- [ ] **Step 2: Add/adjust failing tests for threshold and resonance injection invariants**

Assert that injected grids remain sorted/unique and add points only under trigger conditions.

- [ ] **Step 3: Add/adjust failing tests for adaptive-refinement contract**

Assert refinement reduces interpolation error in a controlled under-resolved case while preserving finite/ordered sample tables.

- [ ] **Step 4: Add/adjust failing tests for struct vs `NamedTuple` equivalence**

Cover preserved entrypoints touched by refinement (`average_scattering_rate`, `build_w0cdf_pchip_cache` paths where practical).

- [ ] **Step 5: Add/adjust failing tests for default-grid vs explicit-grid equivalence**

Add a dedicated assertion for scenarios where the code intends equivalence between internally built default grids and explicit grid inputs, with explicit tolerance.

- [ ] **Step 6: Run targeted unit test file to confirm initial failures are meaningful**

Run: `julia --project=. -e 'include("tests/unit/relaxtime/test_average_scattering_rate.jl")'`
Expected: FAIL only on newly added/refined contract assertions.

- [ ] **Step 7: Commit test-only contract locks**

```bash
git add tests/unit/relaxtime/test_average_scattering_rate.jl
git commit -m "test(relaxtime): lock average scattering cache and refinement contracts"
```

### Task 2: Strengthen `RelaxationTime` contract tests

**Files:**
- Modify: `tests/unit/relaxtime/test_relaxation_time.jl`
- Test target: `src/relaxtime/RelaxationTime.jl`

- [ ] **Step 1: Add/adjust failing tests for alias-resolution contract**

Verify canonical and compatibility keys produce consistent rate lookup behavior.

- [ ] **Step 2: Add/adjust failing tests for tau aggregation algebra**

Verify contribution totals map exactly onto `tau_inv` species totals.

- [ ] **Step 3: Add/adjust failing tests for negative-rate warn-and-clamp behavior**

Add focused test inputs that force a negative intermediate and assert final clamped non-negative outputs (and warning behavior if test harness captures logs).

- [ ] **Step 4: Run targeted unit test file to confirm initial failures are meaningful**

Run: `julia --project=. -e 'include("tests/unit/relaxtime/test_relaxation_time.jl")'`
Expected: FAIL only on newly added/refined contract assertions.

- [ ] **Step 5: Commit test-only contract locks**

```bash
git add tests/unit/relaxtime/test_relaxation_time.jl
git commit -m "test(relaxtime): lock tau aggregation and alias contracts"
```

## Chunk 2: `AverageScatteringRate` Internal Layering Refinement

### Task 3: Refine helper boundaries without changing public behavior

**Files:**
- Modify: `src/relaxtime/AverageScatteringRate.jl`
- Test: `tests/unit/relaxtime/test_average_scattering_rate.jl`

- [ ] **Step 1: Isolate parameter/distribution helper section boundaries**

Reorganize local helper grouping for flavor lookup and isotropic/anisotropic distribution dispatch.

- [ ] **Step 2: Isolate cache/interpolation helper section boundaries**

Reorganize cache storage, slope maintenance, and interpolation helper grouping without changing behavior.

- [ ] **Step 3: Isolate grid-design and sampling helper section boundaries**

Reorganize w0cdf design, threshold layering, resonance layering, and adaptive refinement helper grouping.

- [ ] **Step 4: Isolate cached quadrature helper section boundaries**

Reorganize interval/semi-infinite grid caching helper grouping and default-vs-explicit grid plumbing.

- [ ] **Step 5: Preserve existing keyword defaults and call-chain semantics**

Do not change public signatures unless strictly required for correctness; prefer internal helper extraction and naming clarity.

- [ ] **Step 6: Keep cache and interpolation invariants stable**

Ensure:
- PCHIP behavior stays intact
- out-of-range `sigma(s)` remains `0.0`
- local-window averaging behavior remains consistent

- [ ] **Step 7: Keep sampling/refinement contracts stable**

Ensure threshold layering, resonance layering, and adaptive refinement remain behaviorally equivalent to test-defined expectations.

- [ ] **Step 8: Run focused unit tests for `AverageScatteringRate`**

Run: `julia --project=. -e 'include("tests/unit/relaxtime/test_average_scattering_rate.jl")'`
Expected: PASS.

- [ ] **Step 9: Run relaxtime unit smoke selector to detect local coupling breakage**

Run: `julia --project=. -e 'ENV["UNIT_FILES"]="relaxtime/test_average_scattering_rate.jl,relaxtime/test_relaxation_time.jl"; include("tests/unit/runtests.jl")'`
Expected: PASS for selected files.

- [ ] **Step 10: Commit `AverageScatteringRate` refinement**

```bash
git add src/relaxtime/AverageScatteringRate.jl tests/unit/relaxtime/test_average_scattering_rate.jl
git commit -m "refactor(relaxtime): clarify average scattering kernel boundaries"
```

## Chunk 3: `RelaxationTime` Internal Layering Refinement

### Task 4: Refine rate and tau aggregation boundaries without widening scope

**Files:**
- Modify: `src/relaxtime/RelaxationTime.jl`
- Modify: `tests/unit/relaxtime/test_relaxation_time.jl`
- Optional mechanical-only touch: CSV helper area in same file (no functional broadening)

- [ ] **Step 1: Clarify alias-resolution and normalization helper boundaries**

Group and simplify internal helper flow from input normalization to rate lookup.

- [ ] **Step 2: Clarify average-rate reuse orchestration boundaries**

Keep `compute_average_rates` call behavior stable while making missing-rate and cache handoff behavior easier to audit.

- [ ] **Step 3: Clarify tau aggregation core and contribution-row derivation**

Keep formula behavior stable and keep warn-and-clamp contract explicit.

- [ ] **Step 4: Avoid active CSV-helper refactor**

Leave `_read_sigma_table_csv` and `load_cross_section_caches_from_dir` unchanged except for minimal mechanical consistency edits if required.

- [ ] **Step 5: Run focused unit tests for `RelaxationTime`**

Run: `julia --project=. -e 'include("tests/unit/relaxtime/test_relaxation_time.jl")'`
Expected: PASS.

- [ ] **Step 6: Re-run paired relaxtime unit selector**

Run: `julia --project=. -e 'ENV["UNIT_FILES"]="relaxtime/test_average_scattering_rate.jl,relaxtime/test_relaxation_time.jl"; include("tests/unit/runtests.jl")'`
Expected: PASS.

- [ ] **Step 7: Commit `RelaxationTime` refinement**

```bash
git add src/relaxtime/RelaxationTime.jl tests/unit/relaxtime/test_relaxation_time.jl
git commit -m "refactor(relaxtime): isolate tau aggregation and rate contracts"
```

## Chunk 4: Regression Gate, Governance, and PR Packaging

### Task 5: Apply conditional regression gate and finalize PR contract statement

**Files:**
- Conditional modify: `tests/regression/relaxtime/test_tau_xi_probe_regression.jl`
- Optional notes/doc touch if needed for PR clarity

- [ ] **Step 1: Determine whether fixture-facing semantics changed**

Use this gate:
- interpolation semantics changed?
- cache support-range semantics changed?
- rate-key resolution semantics changed?
- tau aggregation semantics changed?

- [ ] **Step 2: Run or skip probe regression based on gate**

If any gate item is yes:
- Run: `julia --project=. -e 'ENV["REGRESSION_FILES"]="relaxtime/test_tau_xi_probe_regression.jl"; include("tests/regression/runtests.jl")'`
- Expected: PASS (or intentional fixture update with explicit rationale)

If all gate items are no:
- Record explicit PR note: "probe regression skipped because no fixture-facing contract changed"

- [ ] **Step 3: Run docs consistency guard if any user-facing contract doc changed**

Run (only when docs touched): `julia --project=. scripts/dev/check_docs_consistency.jl`
Expected: PASS.

- [ ] **Step 4: Verify clean scope boundary before opening PR**

Run:
- `git status --short`
- `git diff --name-only --cached`
- `git diff --name-only`

Expected: both staged and unstaged sets remain within agreed kernel/test files for this first PR.

- [ ] **Step 5: Commit conditional regression updates if any**

```bash
git add tests/regression/relaxtime/test_tau_xi_probe_regression.jl
git commit -m "test(relaxtime): align probe guard with kernel contract refinements"
```

If no regression file changed, skip this commit step.

### Task 6: Prepare PR message aligned with first-round scope

**Files:**
- No mandatory file edits; PR body content only

- [ ] **Step 1: Draft PR summary with strict scope statement**

Must state:
- this round stabilizes kernel behavior only
- scripts/workflows/outputs are intentionally deferred

- [ ] **Step 2: List behavior contracts now explicitly guarded by tests**

Include cache-range, interpolation/refinement, alias/tau algebra, and clamp policy.

- [ ] **Step 3: Include regression-gate decision in PR body**

State whether probe regression was run or intentionally skipped and why.

- [ ] **Step 4: Final verification sweep before PR**

Run:
- `julia --project=. -e 'include("tests/unit/relaxtime/test_average_scattering_rate.jl")'`
- `julia --project=. -e 'include("tests/unit/relaxtime/test_relaxation_time.jl")'`
- plus conditional regression command from Task 5 if gate says required

Expected: PASS.

- [ ] **Step 5: Build DoD evidence checklist in PR notes**

Map each spec DoD item to concrete evidence:
- scope boundary evidence: file-list checks from Task 5 Step 4
- behavior-preserving refinement evidence: targeted unit test passes
- default-grid vs explicit-grid equivalence evidence: targeted unit assertion pass with stated tolerance
- conditional regression decision evidence: gate decision + command result or explicit skip rationale
- no script/workflow/artifact expansion evidence: PR file list and scope statement

## Notes For Implementers

- Keep changes incremental and commit frequently.
- Avoid broad formatting-only edits in large kernel files.
- If a step reveals hidden coupling to scripts/artifacts, stop and split that work into a follow-up plan rather than widening this PR.
- If behavior-preserving confidence drops, add/strengthen tests before further refactor.
