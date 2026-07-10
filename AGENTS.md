# AGENTS.md

Repository guidance for coding agents working in `Julia_RelaxTime`.

## Working Context

- Primary language: Julia.
- Main environment: `Project.toml` in the repo root.
- Benchmark environment: `benchmark/Project.toml`.
- CI workflows pin Julia `1.12.5`.
- Root `Project.toml` declares Julia `1.10` compatibility.
- This repo is include-driven, not a standard packaged `src/PackageName.jl` layout.
- Prefer unified entrypoints exposed through `Models` and `src/models/entrypoints.jl`.
- `src/pnjl/PNJL.jl` 已从主线移除；若历史文档仍提及该路径，应按“兼容层历史说明”理解，不作为当前实现入口。
- 统一入口请使用 `Models` 与 `src/models/entrypoints.jl`。
- Post-PR76 contract alignment: keep mixed-meson governance and non-fixedmu unified joint-solve semantics unchanged unless explicitly required by task scope.
- Default user-facing communication language is Chinese.

## Repo Rules From Copilot / Cursor

- `.github/copilot-instructions.md` is present and must be followed.
- No `.cursorrules` file exists.
- No `.cursor/rules/` directory exists.
- New core Julia modules must include unit tests.
- Stable public entrypoints should also update `docs/api/`.
- Model config must be read from `config/models/<model>/`.
- Shared physics constants belong in `config/physics/`.
- Test layering must remain `tests/unit/`, `tests/integration/`, `tests/regression/`, `tests/validation/`, and `benchmark/`.
- Non-test scripts do not belong under `tests/`.
- Analysis scripts belong in `scripts/analysis/`.
- Performance probes belong in `scripts/perf/`.
- Smoke tests must stay deterministic, fast, and free of external dependencies.
- Performance work should be justified by profiling or benchmark evidence.

## Codex Collaboration Rules

### Task Execution

- For changes touching `src/models/`, `src/relaxtime/`, or `src/simulation/`, read the relevant entrypoints and the corresponding test layer before editing.
- For changes on unified solver or workflow paths, keep mixed-meson governance and non-fixedmu unified joint-solve semantics unchanged unless the task explicitly requires otherwise.
- For cross-module refactors, prefer the smallest convergent change first; only introduce new shared abstractions when they are clearly required by the task.
- For stable CLI workflows documented in `docs/guides/scripts/README.md`, prefer the `scripts/dev/run_with_sysimage.ps1` or `scripts/dev/run_with_sysimage.sh` wrappers so local sysimage reuse follows repo convention.
- Unless the user explicitly asks for analysis only, prefer implementing, validating, and reporting instead of stopping at a proposal.

### Scope, Verification, And Knowledge Governance

- User approval or task consent applies only to the current stated target, scope, files, approach, risk level, and external impact; if any of these materially changes, stop and re-state the plan before proceeding.
- For changes affecting stable entrypoints, unified solver behavior, numerical results, script contracts, server/API behavior, or user-visible workflows, distinguish "written/configured" from "effective/usable" in the final report.
- When a task intentionally changes semantics, verify both that the target change works and that the explicitly non-changing constraints still hold.
- Keep long-lived facts, stable entrypoints, directory responsibilities, and governance rules in the existing `docs/` structure; keep task process notes in `docs/dev/active/` or archived task records rather than mixing process history into API or architecture references.
- Parallel read-only investigation is acceptable, but edits to the same file, module, or logical block must have a single active writer and be reconciled by the main agent before validation or commit.

### Search Scope And Historical Archives

- Default content search should prioritize current fact sources: `README.md`, `AGENTS.md`, `.github/copilot-instructions.md`, `docs/api/`, `docs/guides/`, `docs/dev/active/`, `docs/dev/backlog/`, `src/`, `scripts/`, `tests/`, and `config/`.
- `docs/dev/archived/**` is a historical audit and completed-task record. Do not include it in broad default content searches because it may preserve obsolete paths, old status, or superseded implementation plans.
- Only search `docs/dev/archived/**` when the task explicitly needs historical traceability, provenance, archived task evidence, old PR context, or an explanation of why a prior decision was made. When using archived content, label it as historical evidence rather than current repository state.

### Numerical And Regression Governance

- For changes affecting numerical results, scan logic, solver paths, cross sections, or relaxation-time workflows, explicitly decide which test layer must cover the change: `tests/unit/`, `tests/integration/`, `tests/regression/`, or `tests/validation/`.
- If a change may introduce numerical drift, run or update the corresponding regression coverage unless there is a concrete reason not to; when regression coverage is skipped, explain the reason and risk.
- Do not loosen tolerances just to make tests pass. Any tolerance adjustment must explain the error source, physical meaning, and validation evidence.
- Performance-oriented changes should be backed by profiling, benchmarks, or existing performance evidence before broad hot-path rewrites are made.

### Documentation And Paper Collaboration

- When stable entrypoints, script contracts, or front/back-end data contracts change, explicitly evaluate whether `docs/api/`, `docs/dev/`, and `docs/guides/` must be updated.
- Paper-writing assistance must not invent experimental conclusions, numerical trends, figure provenance, or literature conclusions. Only write from code, result files, figures, or clearly provided references.
- For tasks involving method description, experimental setup, captions, or result interpretation, check terminology, variable names, units, and implementation consistency before drafting prose.

### Frontend Contract First

- For frontend-facing work where UI implementation is not explicitly requested, prioritize backend contract deliverables: field definitions, units, example payloads, error semantics, and boundary-state behavior.
- Any numerical field exposed to frontend consumers should document units, value domain, missing-value semantics, error semantics, and ordering assumptions where applicable.
- If an interface is not yet stable, prefer producing a mock or intermediate contract document before hard-coding frontend behavior.

### Reusable Workflow Capture

- If the same collaboration pattern or high-frequency task appears three or more times, prefer capturing it as a skill, script entrypoint, or reusable template instead of relying on repeated ad hoc prompting.
- After substantial tasks, when a reusable process or governance rule becomes clear, prefer updating `docs/dev/` or the relevant skill instructions rather than leaving the workflow only in chat history.
- Repository-specific skills under `.agents/skills/` are the canonical project copies. Do not install another skill with the same name in a user-level skill root; keep generic cross-project skills in one user-level root or vendor them here, but never both.
- Keep positive trigger conditions in SKILL.md frontmatter `description`; keep the body focused on execution, boundaries, resources, and validation.
- Claude Code adapters under `.claude/skills/` are generated wrappers; regenerate them with `scripts/dev/sync_claude_skills.jl` and do not edit their bodies manually.
- Keep `.claude/commands/` for true parameterized action macros only. Do not duplicate SKILL workflows in slash-command files.

## Setup Commands

Instantiate and precompile the main Julia environment:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

Run the standard package test entrypoint:

```sh
julia --project=. test/runtests.jl
```

Run layered wrapper entrypoints (preferred for focused local loops):

```sh
julia --project=. test/unit.jl
julia --project=. test/integration.jl
julia --project=. test/regression.jl
julia --project=. test/validation.jl
```

Alternative package test command:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

Start the API plus web server:

```sh
julia --project=. scripts/server/server_full.jl
```

Prefer the stable CLI wrapper on Windows / PowerShell when available:

```sh
powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 scripts/pnjl/calculate_phase_structure.jl --preset=smoke
```

Prefer the stable CLI wrapper on Linux / macOS when available:

```sh
sh scripts/dev/run_with_sysimage.sh scripts/pnjl/calculate_phase_structure.jl --preset=smoke
```

Bootstrap a matching prebuilt sysimage on Windows / PowerShell:

```sh
powershell -ExecutionPolicy Bypass -File scripts/dev/bootstrap_sysimage.ps1
```

Bootstrap a matching prebuilt sysimage on Linux / macOS:

```sh
sh scripts/dev/bootstrap_sysimage.sh
```

Run phase CLI with optional precompile warmup (enabled by default):

```sh
julia --project=. scripts/pnjl/calculate_phase_structure.jl --preset=smoke
```

Disable phase CLI warmup if needed:

```sh
julia --project=. -e 'ENV["PHASE_PRECOMPILE_WARMUP"]="0"; include("scripts/pnjl/calculate_phase_structure.jl")'
```

Select phase CLI precompile profile (default `scan`):

```sh
julia --project=. -e 'ENV["PHASE_PRECOMPILE_PROFILE"]="core"; include("scripts/pnjl/calculate_phase_structure.jl")'
```

Install JS tooling used for dependency graph rendering:

```sh
npm install
```

Render the dependency graph SVG:

```sh
npm run deps:render
```

Regenerate the dependency graph markdown:

```sh
julia --project=. scripts/dev/gen_deps.jl
```

## Test Commands

Unit smoke profile (ultra-fast local edit-run loop):

```sh
julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'
```

Note: local non-smoke profiles emit a warning to encourage smoke-first edit-run loops.

Unit core profile (broader pre-merge gate):

```sh
julia --project=. -e 'ENV["UNIT_PROFILE"]="core"; include("tests/unit/runtests.jl")'
```

Unit full profile:

```sh
julia --project=. -e 'ENV["UNIT_PROFILE"]="full"; include("tests/unit/runtests.jl")'
```

Run a single unit test file through the supported selector:

```sh
julia --project=. -e 'ENV["UNIT_FILES"]="relaxtime/test_transport_coefficients.jl"; include("tests/unit/runtests.jl")'
```

Direct single-file unit iteration:

```sh
julia --project=. -e 'include("tests/unit/relaxtime/test_transport_coefficients.jl")'
```

Integration smoke profile (ultra-fast local edit-run loop):

```sh
julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'
```

Integration core profile (broader pre-merge gate):

```sh
julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="core"; include("tests/integration/runtests.jl")'
```

Integration full profile:

```sh
julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="full"; include("tests/integration/runtests.jl")'
```

Run a single integration test file through the supported selector:

```sh
julia --project=. -e 'ENV["INTEGRATION_FILES"]="relaxtime/test_transport_workflow_smoke.jl"; include("tests/integration/runtests.jl")'
```

Direct single integration test file:

```sh
julia --project=. -e 'include("tests/integration/relaxtime/test_transport_workflow_smoke.jl")'
```

Regression smoke profile (ultra-fast local edit-run loop):

```sh
julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'
```

Regression core profile (broader pre-merge gate):

```sh
julia --project=. -e 'ENV["REGRESSION_PROFILE"]="core"; include("tests/regression/runtests.jl")'
```

Regression full profile:

```sh
julia --project=. -e 'ENV["REGRESSION_PROFILE"]="full"; include("tests/regression/runtests.jl")'
```

Run a single regression file through the supported selector:

```sh
julia --project=. -e 'ENV["REGRESSION_FILES"]="relaxtime/test_transport_fixedpoint_regression.jl"; include("tests/regression/runtests.jl")'
```

Validation suite:

```sh
julia --project=. tests/validation/runtests.jl
```

Validation smoke profile:

```sh
julia --project=. -e 'ENV["VALIDATION_PROFILE"]="smoke"; include("tests/validation/runtests.jl")'
```

Validation core profile:

```sh
julia --project=. -e 'ENV["VALIDATION_PROFILE"]="core"; include("tests/validation/runtests.jl")'
```

Enable optional precompile warmup for test runners (helps AD-heavy cold start):

```sh
julia --project=. -e 'ENV["TEST_PRECOMPILE_WARMUP"]="1"; ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'
```

Choose precompile capability profile for test runners (default `test`):

```sh
julia --project=. -e 'ENV["TEST_PRECOMPILE_WARMUP"]="1"; ENV["TEST_PRECOMPILE_PROFILE"]="core"; ENV["INTEGRATION_PROFILE"]="core"; include("tests/integration/runtests.jl")'
```

Run a single validation file through the supported selector:

```sh
julia --project=. -e 'ENV["VALIDATION_FILES"]="relaxtime/test_mott_reference_mapping.jl"; include("tests/validation/runtests.jl")'
```

VS Code / TestItems wrapper entrypoint:

```sh
julia --project=. test/vscode_testitems.jl
```

## Governance And Audit Commands

There is no repo-wide `JuliaFormatter` or dedicated linter config in the repository.
Preserve local style and use the governance scripts as the nearest equivalent to lint / policy checks.

```sh
julia --project=. scripts/dev/check_unit_skip_policy.jl
julia --project=. scripts/dev/check_docs_consistency.jl
julia --project=. scripts/dev/check_active_docs_governance.jl
julia --project=. scripts/dev/check_script_entrypoints.jl
julia --project=. scripts/dev/check_models_entry_contract.jl
julia --project=. scripts/dev/check_solver_contract_leakage.jl
julia --project=. scripts/dev/check_relaxtime_script_governance.jl
julia --project=. scripts/dev/check_pnjl_migration_guard.jl
julia --project=. scripts/dev/analyze_deps.jl
julia --project=. scripts/dev/check_model_profile_matrix.jl
julia --project=. scripts/dev/check_data_output_path_guard.jl
julia --project=. scripts/dev/check_precompile_profile_coverage.jl
julia --project=. scripts/dev/check_precompile_trace_budget.jl
julia --project=. scripts/dev/check_skill_governance.jl
julia --project=. scripts/dev/sync_claude_skills.jl --check
```

Build persistent sysimage for AD-heavy test/script workflows:

```sh
julia --project=. scripts/dev/build_sysimage.jl
```

Use the sysimage for faster cold-start test runs:

```sh
julia --sysimage=build/JuliaRelaxTime.dll --project=. -e 'ENV["INTEGRATION_PROFILE"]="core"; include("tests/integration/runtests.jl")'
```

## Benchmark Commands

```sh
julia --project=benchmark -e 'using PkgBenchmark; benchmarkpkg(".")'
julia --project=benchmark benchmark/relaxtime/bench_total_cross_section.jl
```

## Repository Layout

- `src/`: reusable Julia modules.
- `src/models/`: unified model architecture and workflow entrypoints.
- `src/relaxtime/`: transport, scattering, propagators, and relaxation-time logic.
- `src/simulation/`: HTTP server and simulation helpers.
- `tests/unit/`: deterministic fast unit tests.
- `tests/integration/`: cross-module workflow tests.
- `tests/regression/`: internal numerical baseline regression tests.
- `tests/validation/`: external reference and acceptance tests.
- `benchmark/`: performance benchmarks only.
- `scripts/analysis/`: analysis and diagnostics, not automated tests.
- `scripts/perf/`: profiling and performance probes, not correctness tests.
- `docs/api/`: public API docs for stable entrypoints.
- `docs/dev/active/`: active task docs that should be archived after completion.

## Code Style

- Follow the style already present in the edited file; do not apply broad reformatting.
- Keep modules and types in `PascalCase`.
- Keep functions in `snake_case`.
- Keep test files named `test_*.jl`.
- Keep benchmark files named `bench_*.jl` or `benchmark_*.jl` under `benchmark/`.
- Prefer small focused files with one clear responsibility.
- Do not move exploratory code into `src/`.

## Imports, Types, And Naming

- Place `using` and `import` near the top of each module, after any module docstring.
- Prefer explicit imports such as `using Main.ParameterTypes: QuarkParams, ThermoParams`.
- Use `import` when extending methods and `using` when consuming names.
- Relative imports like `using ..GaussLegendre` are normal in this codebase.
- Preserve `if !isdefined(Main, :ModuleName)` guards when shared modules may already be loaded into `Main`.
- Prefer concrete structs for stable API inputs, but preserve existing `NamedTuple` compatibility where supported.
- Use multiple dispatch instead of large flag-heavy functions when extending behavior.
- Avoid `Any`-typed containers in hot paths and use `const` for global bindings.
- Internal physics quantities use natural units; for new fm^-1 values prefer the `_inv_fm` suffix.
- Keep external MeV-facing names explicit, such as `T_MeV` or `mu_MeV`.
- Higher-order units should follow local conventions like `sigma_fm4` or `coupling_inv_fm4`.

## Error Handling, Performance, And Tests

- Validate public inputs early.
- Prefer clear `ArgumentError` messages that name the offending parameter.
- Reuse shared validation helpers from `src/utils/ValidationUtils.jl` when applicable.
- Use `@warn` for recoverable fallback behavior such as config fallback paths.
- Do not silently ignore invalid physical or numerical inputs.
- Keep numerical safeguards explicit around finiteness, positivity, and grid compatibility.
- Avoid unnecessary allocations in hot loops and integration kernels.
- Avoid repeated recalculation inside solver or quadrature loops.
- Only add `@inline` or `@inbounds` when justified by hot-path use.
- Put performance comparisons in `benchmark/`, not in smoke tests.
- Every new core module should include unit tests, and stable public entrypoints should also update `docs/api/`.

## Commit Message Governance (Mandatory)

When creating any git commit in this repository, the agent MUST follow historical commit style.

1) Style sampling (required)
- Run: git log -10 --oneline
- Extract prefix patterns from recent history (e.g. docs:, fix:, refactor:, ci:, feat:, or type(scope):).

2) Style matching (required)
- The new commit message MUST use an existing prefix pattern already present in sampled history.
- Prefer matching the nearest 3 commits in the same change category.
- If multiple styles exist, choose the dominant recent style; do not invent a new format.

3) Message quality (required)
- Keep message concise and intent-driven (focus on why, not only what).
- Subject line should be one line, clear scope, no filler words.
- Avoid vague messages like “update files”, “fix issues”, “misc changes”.

4) Pre-commit self-check (required)
- Before commit, the agent must verify:
  - Prefix matches historical style
  - Scope/category is accurate
  - Message is concise and meaningful

5) Fallback rule
- If style cannot be determined confidently from recent history, use:
  - docs: ... for documentation-only changes
  - fix: ... for bug fixes
  - refactor: ... for non-functional code restructuring
  - ci: ... for workflow/pipeline changes
  and note that fallback was used.

6) Non-compliance
- If the message does not match repository style, the agent must revise it before committing.

## Literature-To-Implementation Subagents

This repository may use project subagents for literature-backed computational decisions. The goal is not to maintain a full bibliography here; the goal is to turn literature evidence into formulas, implementation choices, validation plans, regression targets, and paper handoff notes.

Use these project agents when the user explicitly asks for subagents or parallel agents on literature/research-engineering tasks:

- `relax-literature-search-strategist`: targeted searches for formulas, algorithms, parameterizations, validation points, and reproducibility signals.
- `relax-method-reviewer`: equation assumptions, notation, units, reproducibility, implementation risk, and validation requirements.
- `relax-evidence-synthesizer`: method comparison, project-fit synthesis, recommended path, test layer implications, and documentation impacts.
- `relax-gap-analyst`: missing evidence, formula gaps, validation gaps, regression needs, docs tasks, and paper-workspace handoffs.

Citation and BibTeX cleanup belongs in `D:\Desktop\paper`, especially under `D:\Desktop\paper\bib`. Do not edit that master bibliography from this repository. If citation metadata or duplicate cleanup is needed, return a handoff request for the `paper-citation-curator`.

For protocol details, see `docs/analysis/literature_to_implementation_protocol.md`.
