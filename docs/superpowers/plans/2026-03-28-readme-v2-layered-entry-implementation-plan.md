# README v2 Layered Entry Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Deliver a concise README v2 that is new-user-first, script-first, and governance-aligned with stable entrypoint policy.

**Architecture:** Rewrite `README.md` into 8 fixed sections from the approved design spec, keep deep content out of README, and route details to existing docs. Validate through docs/script governance checks and unit smoke, then confirm README quickstart artifact contract and git-hygiene behavior.

**Tech Stack:** Julia 1.10+, markdown docs, governance scripts under `scripts/dev/`, unit smoke entry under `tests/unit/runtests.jl`.

---

## File Structure (Planned)

- Modify: `README.md`
- Modify: `.gitignore` (only if choosing ignore-path strategy for quickstart outputs)
- Modify: `docs/guides/STATUS.md` (only if README wording change requires sync)
- Create: `docs/dev/active/readme_v2_content_migration_map.md` (temporary implementation artifact, can be archived/removed after merge)
- Test/Check: `scripts/dev/check_docs_consistency.jl`
- Test/Check: `scripts/dev/check_script_entrypoints.jl`
- Test/Check: `tests/unit/runtests.jl` (smoke profile)

## Execution Safety Baseline

- Record `BASE_SHA` before implementation starts:
  - `git rev-parse HEAD`
- Use this `BASE_SHA` for any file-scoped rollback in this plan.

## Chunk 0: Mandatory Deliverables From Spec

### Task 0: Produce migration map, text budget, and verification matrix before README edits

**Files:**
- Create: `docs/dev/active/readme_v2_content_migration_map.md`
- Reference: `README.md`
- Reference: `docs/superpowers/specs/2026-03-28-readme-v2-layered-entry-design.md`

- [ ] **Step 1: Write old README to docs migration map**
  - For each major old README block, record: source section, target path, action (`keep-short-link`, `move`, `drop`).

- [ ] **Step 2: Write 8-section line budget table**
  - Include per-section max lines and total max <=250 lines.

- [ ] **Step 3: Write verification matrix**
  - Include checks for AC #1-#8, command/URL exact match expectations, and pass/fail evidence fields.

- [ ] **Step 4: Run self-audit against spec mandatory deliverables**
  - Expected: PASS for mandatory deliverables #1/#2/#3 coverage.

- [ ] **Step 5: Commit**
  - `git add docs/dev/active/readme_v2_content_migration_map.md`
  - `git commit -m "docs: add README v2 migration map and validation matrix"`

## Chunk 1: README v2 Content Skeleton + Quickstart Contract

### Task 1: Build failing contract checklist for README v2 shape

**Files:**
- Modify: `README.md`
- Test/Check: manual checklist based on `docs/superpowers/specs/2026-03-28-readme-v2-layered-entry-design.md`

- [ ] **Step 1: Define explicit acceptance checklist before editing README**
  - Include: 8 fixed sections, <=250 lines target, `Latest Release` within first 20 lines, exact URL `https://github.com/w5851/Julia_RelaxTime/releases/latest`, one primary script-chain quickstart block, one optional server command, and explicit "5-10 minutes" quickstart expectation text.

- [ ] **Step 2: Run pre-edit baseline check against current README**
  - Run: manual checklist pass on current `README.md`.
  - Expected: FAIL on section shape and concision constraints.

- [ ] **Step 3: Draft README v2 skeleton only (headings + placeholders)**
  - Keep all 8 required sections in target order.

- [ ] **Step 4: Re-check skeleton contract**
  - Expected: PASS for structure, pending for detailed content.

- [ ] **Step 5: Commit**
  - `git add README.md`
  - `git commit -m "docs: scaffold README v2 layered-entry structure"`

### Task 2: Implement quickstart and stable-entry matrix with script-first priority

**Files:**
- Modify: `README.md`
- Reference: `docs/guides/scripts/README.md`

- [ ] **Step 1: Add quickstart failing checks (self-audit list)**
  - Validate command executability, artifact verification command, and secondary server path wording.

- [ ] **Step 2: Add minimal executable quickstart content**
  - Step 1 env prep command.
  - Step 2 `calculate_phase_structure.jl --preset=smoke` command with output directory.
  - Step 3 artifact verification command.
  - Add expected completion window text: `5-10 minutes` after environment preparation.

- [ ] **Step 3: Add stable script-entry matrix from whitelist authority**
  - State `docs/guides/scripts/README.md` as canonical stable-entry source.
  - Include PNJL and Server/API stable entrypoints only.

- [ ] **Step 4: Re-run quickstart self-audit checks**
  - Expected: PASS on executable wording and boundary consistency.

- [ ] **Step 5: Commit**
  - `git add README.md`
  - `git commit -m "docs: add script-first quickstart and stable entry matrix"`

## Chunk 2: Boundary/Governance Sections + Verification Closure

### Task 3: Add capability boundaries, docs navigation, and concise repo map

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Add failing boundary checklist**
  - Must mention non-user-entry directories (`scripts/dev`, `scripts/analysis`, `scripts/debug`, `scripts/perf`) and compatibility-only historical path policy.

- [ ] **Step 2: Create section placeholders for boundary/docs-navigation/repo-map**
  - Keep architecture narrative placeholder constrained to navigation-only intent.

- [ ] **Step 3: Implement capability boundary section only and self-check**
  - PASS condition: all non-user-entry directories are explicitly named.

- [ ] **Step 4: Implement docs navigation section only and self-check**
  - PASS condition: user/dev/API routes each have at least one link target.

- [ ] **Step 5: Implement minimal repo map section only and self-check**
  - PASS condition: top-level map remains concise and no deep architecture explanation is introduced.
  - PASS condition: include at least one deep-doc link (prefer `docs/architecture/`) for architecture details.

- [ ] **Step 6: Enforce concision budget**
  - Ensure README remains <=250 lines.

- [ ] **Step 7: Enforce architecture-navigation hard cap**
  - PASS condition: architecture/navigation block in README is <=12 lines.

- [ ] **Step 8: Re-check boundary and concision checklist**
  - Expected: PASS.

- [ ] **Step 9: Commit**
  - `git add README.md`
  - `git commit -m "docs: tighten README boundaries and docs navigation"`

### Task 3B: Fill README verification section and contribution/governance section explicitly

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Implement README Section 5 (Verification Commands)**
  - Include artifact-level verification command and maintainer checks (unit smoke + docs consistency).

- [ ] **Step 2: Implement README Section 8 (Contribution + Governance links)**
  - Add links to `.github/CONTRIBUTING.md`, `.github/SECURITY.md`, and code-of-conduct entry.

- [ ] **Step 3: Verify section coverage**
  - PASS condition: Section 5 and Section 8 are both non-empty and actionable.

- [ ] **Step 4: Commit**
  - `git add README.md`
  - `git commit -m "docs: complete README verification and governance sections"`

### Task 4: Resolve quickstart output git-hygiene policy

**Files:**
- Modify: `.gitignore` (optional path A)
- Modify: `README.md` (required if choosing path B cleanup)

- [ ] **Step 1: Validate existing ignore coverage for `data/outputs/results/phase_smoke`**
  - Run: inspect `.gitignore` and check whether quickstart output path is ignored.
  - Expected: likely not ignored in current repo state.

- [ ] **Step 2: Choose one strategy and implement minimally**
  - Path A: add ignore rule for quickstart output directory.
  - Path B: keep repository ignore policy unchanged, add explicit cleanup command in README quickstart section.

- [ ] **Step 3: Verify git-hygiene outcome**
  - Run quickstart output simulation or path check to ensure demo instructions do not leave confusing dirty status.
  - PASS condition: either output path is ignored or cleanup command reliably returns clean status for demo artifacts.
  - Expected: PASS for either ignore or cleanup contract.

- [ ] **Step 4: Update README verification section to reflect chosen strategy**
  - Ensure wording is concrete and actionable.

- [ ] **Step 5: Commit**
  - `git add README.md .gitignore`
  - `git commit -m "docs: enforce README quickstart git-hygiene contract"`

### Task 5: Run governance checks, smoke checks, and finalize docs sync

**Files:**
- Modify: `README.md`
- Modify: `docs/guides/STATUS.md` (only if consistency checker requires synchronized wording)

- [ ] **Step 1: Run docs consistency gate**
  - Run: `julia --project=. scripts/dev/check_docs_consistency.jl`
  - Expected: PASS.

- [ ] **Step 2: Run script entrypoint governance gate**
  - Run: `julia --project=. scripts/dev/check_script_entrypoints.jl`
  - Expected: PASS.

- [ ] **Step 3: Run unit smoke profile for README verify command validity**
  - Run: `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
  - Expected: PASS.

- [ ] **Step 4: Run link validation gate (auto or fallback)**
  - Preferred auto command:
    - `julia --project=. scripts/dev/check_docs_consistency.jl`
    - (treat this as link/consistency gate in this repository)
  - Fallback: manually verify all internal README links resolve and all referenced files/paths exist.
  - Expected: PASS.

- [ ] **Step 5: Apply minimal sync fixes if any check reports drift**
  - Update README and `docs/guides/STATUS.md` only as needed.
  - Re-run only failed checks until PASS.

- [ ] **Step 6: Final commit**
  - `git add README.md docs/guides/STATUS.md`
  - `git commit -m "docs: finalize README v2 with governance-verified entry contract"`

## Rollback Procedure

- If rewrite introduces drift/noise, revert only files touched by this plan (never full-tree rollback), then re-run:
  - `git restore --source=<BASE_SHA> -- README.md docs/guides/STATUS.md .gitignore docs/dev/active/readme_v2_content_migration_map.md`
  - `julia --project=. scripts/dev/check_docs_consistency.jl`
  - `julia --project=. scripts/dev/check_script_entrypoints.jl`

## Definition of Done

- README satisfies approved 8-section design and line budget.
- Script-first quickstart is executable and artifact-verified.
- Stable-entry boundary and non-user-entry policy are explicit.
- Release anchor policy is correctly reflected.
- Governance checks and unit smoke pass after rewrite.
