# README v2 Layered Entry Design

## Context

Current `README.md` mixes unit convention, feature details, historical updates, architecture description, migration notes, and quickstart guidance in one long page. This makes first-run onboarding slower and weakens the role of README as an execution-oriented entrypoint.

The user has explicitly confirmed these priorities for README v2:

1. Audience priority: new users first.
2. Workflow priority: script-based path first (not frontend-first).
3. Scope policy: keep README concise; move long background/history into docs.
4. Structure preference: layered entry layout (option A).

## Goals

1. Make README operational for first-time users, with quickstart phase run targeted at 5-10 minutes on a typical developer machine after environment preparation.
2. Freeze stable script-entry messaging to avoid "too many entrypoints" ambiguity.
3. Clarify capability boundaries (what is stable vs what is internal/experimental).
4. Keep release anchor and verification commands visible at top-level.

## Non-Goals

1. Do not turn README into a complete architecture manual.
2. Do not duplicate long theory/history content already maintained in `docs/`.
3. Do not present frontend as the primary production capability path.
4. Do not copy large architecture/principle blocks from old README into README v2.
5. Architecture/deep-principle content should be migrated to `docs/architecture/` (or existing guide/reference docs), with README keeping only links.

## Selected Approach

Adopt a layered-entry README with 8 fixed sections:

1. One-line project positioning + Latest Release anchor.
2. 3-step quickstart (script workflow first).
3. Stable script-entry matrix (PNJL + server as secondary).
4. Capability boundaries and non-target scope.
5. Recommended verification commands.
6. Documentation navigation (user/dev/API).
7. Minimal repository structure map.
8. Contribution and governance links.

This approach keeps README as an operational index and pushes deep content to guided documents.

## Information Architecture

### Section 1: Positioning + Release

- One sentence for repository identity (PNJL/NJL + transport + reproducible workflows).
- Keep `Latest Release` link within first 20 lines for capability contract alignment.
- Canonical target is the non-404-safe URL:
  - `https://github.com/w5851/Julia_RelaxTime/releases/latest`
- If using version-pinned links elsewhere, update them in the same change batch as release creation.

### Section 2: Quickstart (Script First)

- Quickstart is fixed to exactly 3 steps, all runnable from repo root.
- Step 1 (environment):
  - `julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'`
- Step 2 (primary script workflow):
  - `julia --project=. scripts/pnjl/calculate_phase_structure.jl --preset=smoke --output_dir=data/outputs/results/phase_smoke`
- Step 3 (artifact validation):
  - `julia --project=. -e 'println(isfile("data/outputs/results/phase_smoke/phase_summary.json") && isfile("data/outputs/results/phase_smoke/phase_report.md"))'`
- Output directory hygiene requirement:
  - ensure quickstart output path is covered by `.gitignore` policy, or provide cleanup command in README so sample runs do not pollute `git status`.
- Secondary path is listed as optional:
  - `julia --project=. scripts/server/server_full.jl` as API/web framework entry, not core workflow claim.

### Section 3: Stable Entry Matrix

- State that stable script entrypoints are governed by:
  - `docs/guides/scripts/README.md`
- Surface only whitelist-level items, grouped by domain:
  - PNJL
  - Server/API

### Section 4: Capability Boundaries

- Stable user entrypoints are whitelist-based.
- `scripts/dev/`, `scripts/analysis/`, `scripts/debug/`, `scripts/perf/` are not default user entrypoints.
- Historical paths are compatibility-only and not recommended for new workflows.

### Section 5: Verification Commands

- Provide concise "trust but verify" command set:
  - quickstart output contract check (artifact-level)
  - optional maintainer checks: unit smoke + docs consistency
- Keep list short to reduce friction.

### Section 6: Docs Navigation

- Route users by intent:
  - user guides / quickstart
  - developer process docs
  - API docs

### Section 7: Minimal Repo Map

- Show only top-level directories with role labels.
- Avoid long per-submodule explanations.

### Section 8: Contribution + Governance

- Keep pointers to `.github/CONTRIBUTING.md`, `.github/SECURITY.md`, code of conduct.
- Add one line on docs-governance principle (README as index, docs as source of detail).

## Data and Flow Contract

1. First-time user opens README.
2. User chooses one of stable quick paths (default: script-based PNJL workflow).
3. User verifies output artifacts (`phase_summary.json`, `phase_report.md`).
4. User deep-dives using linked docs by intent.

This flow intentionally minimizes branching at entry and postpones advanced detail.

## Error-Handling and Drift Control

1. If new script entrypoints are introduced, README only references them after they are added to `docs/guides/scripts/README.md` stable list.
2. If release links change, README release anchor must update in same docs/governance change batch.
3. Stable-entrypoint drift gate is required before README entrypoint update is merged:
   - run script-entry governance checker (`scripts/dev/check_script_entrypoints.jl`).
   - if gate fails, merge must be blocked; in CI this check should be configured as required.
4. When docs and README wording diverge, docs consistency checker is required gate.

## Testing and Verification Plan

1. Run docs consistency check after README rewrite:
   - `julia --project=. scripts/dev/check_docs_consistency.jl`
2. Run unit smoke to confirm recommended command remains valid:
   - `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
3. Run script entrypoint governance check:
   - `julia --project=. scripts/dev/check_script_entrypoints.jl`
4. Quickstart artifact verification:
   - confirm `data/outputs/results/phase_smoke/phase_summary.json` exists
   - confirm `data/outputs/results/phase_smoke/phase_report.md` exists
5. Git hygiene verification:
   - confirm quickstart output path is ignored by `.gitignore` (directory or pattern-level ignore).
6. Automated link validation in CI (preferred):
   - run markdown link check tool on `README.md` (or equivalent repository-approved link checker).
7. Manual markdown validation fallback (only if automated checker is unavailable):
   - confirm all internal links resolve
   - confirm quickstart commands match actual scripts and current paths

## Risks and Mitigations

1. Risk: oversimplified README hides important constraints.
   - Mitigation: explicit capability-boundary section + direct links to status and scripts guide.
2. Risk: stale entrypoint references.
   - Mitigation: whitelist authority in `docs/guides/scripts/README.md` and consistency checks.
3. Risk: quickstart path drifts from real workflows.
   - Mitigation: keep one canonical script-first minimal path and verify in routine docs updates.

## Acceptance Criteria

1. README main body is structured into the 8 approved sections and remains within 250 lines.
2. Quickstart section contains exactly one primary script-chain command block and one optional server entry command.
3. Stable entrypoint policy and non-user-entry directories are explicitly documented.
4. `Latest Release` anchor appears within first 20 lines and points to `https://github.com/w5851/Julia_RelaxTime/releases/latest`.
5. Verification section includes artifact-level check plus maintainer checks, and commands are executable in current repository layout.
6. README architecture content is constrained to a navigation-only block of at most 12 lines; deep explanations are linked to `docs/architecture/` or other docs pages.
7. Quickstart output directory path is covered by `.gitignore` policy, or README includes an explicit cleanup step to keep `git status` clean after demo runs.
8. Quickstart primary command uses smoke preset and documents expected completion window (5-10 minutes on a typical developer machine with warm Julia environment).

## Implementation Deliverables (Mandatory)

The implementation plan must include the following deliverables:

1. Exact content cut/move map from old README to destination docs.
2. Section-by-section target text budget for README v2.
3. Verification procedure (artifact checks + governance checks + link checks).
4. Rollback procedure (restore previous README from git and re-run docs checks).
5. `.gitignore` verification note for quickstart output directories.
