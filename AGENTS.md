# AGENTS.md

Repository guidance for coding agents working in `Julia_RelaxTime`.

## Working Context

- Primary language: Julia.
- Main environment: `Project.toml` in the repo root.
- Benchmark environment: `benchmark/Project.toml`.
- CI workflows pin Julia `1.12.5`; root compatibility remains Julia `1.10`.
- This repo is include-driven, not a standard packaged `src/PackageName.jl` layout.
- Prefer unified entrypoints exposed through `Models` and `src/models/entrypoints.jl`.
- `src/pnjl/PNJL.jl` is historical compatibility context, not a current implementation entrypoint.
- Keep mixed-meson governance and non-fixedmu unified joint-solve semantics unchanged unless task scope explicitly requires otherwise.
- Default user-facing communication language is Chinese.
- QuadGK is not a root runtime/test dependency. It is allowed only through the explicitly activated/stacked `benchmark/` environment as an external comparison oracle; see `docs/decisions/0006-isolate-optional-numerical-oracles.md`.

## Repo Rules From Copilot / Cursor

- `.github/copilot-instructions.md` is present and must be followed.
- No `.cursorrules` file or `.cursor/rules/` directory exists.
- New core Julia modules require unit tests; stable public entrypoints also require `docs/api/` updates.
- Model config belongs in `config/models/<model>/`; shared physics constants belong in `config/physics/`.
- Preserve test layering: `tests/unit/`, `tests/integration/`, `tests/regression/`, `tests/validation/`, and `benchmark/`.
- Non-test scripts do not belong under `tests/`; analysis scripts belong in `scripts/analysis/`; performance probes belong in `scripts/perf/`.
- Smoke tests must be deterministic, fast, and free of external dependencies.
- Performance work must be justified by profiling or benchmark evidence.

## Codex Collaboration Rules

### Task Execution

- Before editing `src/models/`, `src/relaxtime/`, or `src/simulation/`, read the relevant entrypoints and corresponding test layer.
- For unified solver or workflow changes, preserve mixed-meson governance and non-fixedmu unified joint-solve semantics unless explicitly in scope.
- For cross-module refactors, prefer the smallest convergent change and introduce shared abstractions only when required.
- For stable CLI workflows, prefer `scripts/dev/run_with_sysimage.ps1` or `scripts/dev/run_with_sysimage.sh` as documented in `docs/guides/scripts/README.md`.
- Unless the user asks for analysis only, implement, validate, and report rather than stopping at a proposal.

### Scope, Verification, And Knowledge Governance

- User approval applies only to the stated target, scope, files, approach, risk level, and external impact. If any materially changes, stop and re-state the plan.
- For changes affecting stable entrypoints, solver behavior, numerical results, script contracts, server/API behavior, or user-visible workflows, distinguish written/configured behavior from effective/usable behavior.
- When semantics change, verify both the target change and the explicitly non-changing constraints.
- Keep long-lived facts and governance in the existing `docs/` structure; keep task process notes in `docs/dev/active/` or archived task records.
- Parallel read-only investigation is acceptable, but edits to one file/module/logical block require a single writer and main-agent reconciliation before validation or commit.

### Search Scope And Historical Archives

- Default searches prioritize `README.md`, `AGENTS.md`, `.github/copilot-instructions.md`, `docs/api/`, `docs/guides/`, `docs/dev/active/`, `docs/dev/backlog/`, `src/`, `scripts/`, `tests/`, and `config/`.
- Do not include `docs/dev/archived/**` in broad default searches; it may preserve obsolete paths and superseded plans.
- Search archives only for explicit historical traceability, provenance, old PR context, or prior-decision evidence, and label those results as historical.

### Numerical And Regression Governance

- For numerical, scan, solver, cross-section, or relaxation-time changes, explicitly select the required unit/integration/regression/validation layer.
- If numerical drift is possible, run or update corresponding regression coverage; if skipped, state the reason and risk.
- Never loosen tolerances merely to pass tests. Explain the error source, physical meaning, and evidence for any tolerance change.
- Back performance changes with profiling, benchmarks, or existing performance evidence before broad hot-path rewrites.

### Documentation And Paper Collaboration

- When stable entrypoints, scripts, or data contracts change, evaluate `docs/api/`, `docs/dev/`, and `docs/guides/` impacts.
- Do not invent experimental conclusions, numerical trends, figure provenance, or literature conclusions.
- For methods, setup, captions, or interpretation, verify terminology, variables, units, and implementation consistency first.

### Frontend Contract First

- If UI implementation is not explicitly requested, prioritize backend field definitions, units, examples, error semantics, and boundary behavior.
- Numerical fields exposed to frontend consumers must document units, domain, missing/error semantics, and ordering assumptions.
- Prefer a mock/intermediate contract before hard-coding an unstable interface.

### Reusable Workflow Capture

- After the same collaboration pattern appears three or more times, capture it as a skill, script entrypoint, or template.
- Put durable process rules in `docs/dev/` or the relevant skill rather than leaving them only in chat history.
- Repository skills under `.agents/skills/` are canonical; do not duplicate the same skill name in a user-level root.
- Keep positive triggers in SKILL.md frontmatter; keep bodies focused on execution, boundaries, resources, and validation.
- Regenerate `.claude/skills/` adapters with `scripts/dev/sync_claude_skills.jl`; never edit generated bodies manually.
- Reserve `.claude/commands/` for true parameterized action macros, not duplicated skill workflows.

## Command Reference

- Environment setup, test profiles, governance checks, sysimage commands, and benchmarks: `docs/dev/agent_command_reference.md`.
- Stable user-facing script entrypoints and wrappers: `docs/guides/scripts/README.md`.
- Test organization and smoke/core/full semantics: `docs/dev/testing_governance.md`.
- Use focused selectors first; expand validation in proportion to numerical and workflow risk.

## Repository Layout

- `src/`: reusable Julia modules.
- `src/models/`: unified model architecture and workflow entrypoints.
- `src/relaxtime/`: transport, scattering, propagators, and relaxation-time logic.
- `src/simulation/`: HTTP server and simulation helpers.
- `tests/unit/`: deterministic fast unit tests.
- `tests/integration/`: cross-module workflow tests.
- `tests/regression/`: internal numerical baseline regression tests.
- `tests/validation/`: external reference and acceptance tests.
- `benchmark/`: performance benchmarks and isolated comparison-oracle dependencies only.
- `scripts/analysis/`: analysis and diagnostics, not automated tests.
- `scripts/perf/`: profiling and performance probes, not correctness tests.
- `docs/api/`: public API docs for stable entrypoints.
- `docs/dev/active/`: active task docs; archive them only when the full task is complete.

## Code Style

- Follow local file style; do not apply broad reformatting.
- Modules and types use `PascalCase`; functions use `snake_case`.
- Test files use `test_*.jl`; benchmark files use `bench_*.jl` or `benchmark_*.jl` under `benchmark/`.
- Prefer small focused files with one clear responsibility; do not move exploratory code into `src/`.

## Imports, Types, And Naming

- Place `using` and `import` near the module top after any module docstring.
- Prefer explicit imports such as `using Main.ParameterTypes: QuarkParams, ThermoParams`.
- Use `import` when extending methods and `using` when consuming names.
- Relative imports such as `using ..GaussLegendre` are normal.
- Preserve `if !isdefined(Main, :ModuleName)` guards where shared modules may already be loaded into `Main`.
- Prefer concrete structs for stable API inputs while preserving existing `NamedTuple` compatibility.
- Prefer multiple dispatch over large flag-heavy functions; avoid `Any` containers in hot paths and use `const` for global bindings.
- Internal quantities use natural units. Use `_inv_fm` for new fm^-1 names and explicit MeV-facing names such as `T_MeV` or `mu_MeV`.
- Higher-order units follow local forms such as `sigma_fm4` or `coupling_inv_fm4`.

## Error Handling, Performance, And Tests

- Validate public inputs early and use clear `ArgumentError` messages naming the invalid parameter.
- Reuse `src/utils/ValidationUtils.jl` where applicable.
- Use `@warn` for recoverable fallback behavior; never silently ignore invalid physical or numerical inputs.
- Keep finiteness, positivity, and grid-compatibility safeguards explicit.
- Avoid allocations and repeated calculations in hot loops; add `@inline` or `@inbounds` only with hot-path justification.
- Put performance comparisons in `benchmark/`, not smoke tests.
- Every new core module needs unit tests; stable public entrypoints also need `docs/api/` updates.

## Commit Message Governance (Mandatory)

Before every commit:

1. Run `git log -10 --oneline` and sample recent prefix patterns.
2. Use an existing recent prefix (`docs:`, `fix:`, `refactor:`, `ci:`, `feat:`, or an observed scoped form).
3. Prefer the nearest three commits in the same category and keep the subject concise, intent-driven, and single-line.
4. Verify that prefix, scope/category, and wording match the actual staged change.
5. If history is inconclusive, use `docs:` for docs-only, `fix:` for bugs, `refactor:` for non-functional restructuring, or `ci:` for workflow changes, and note the fallback.
6. Revise any non-compliant message before committing.

## Literature-To-Implementation Subagents

Use these project agents only when the user explicitly requests subagents or parallel literature/research-engineering work:

- `relax-literature-search-strategist`: formulas, algorithms, parameterizations, validation points, and reproducibility signals.
- `relax-method-reviewer`: assumptions, notation, units, reproducibility, implementation risk, and validation requirements.
- `relax-evidence-synthesizer`: method comparison, project fit, recommended path, test implications, and documentation impacts.
- `relax-gap-analyst`: evidence/formula/validation gaps, regression needs, docs tasks, and paper handoffs.

Citation and BibTeX cleanup belongs in `D:\Desktop\paper`, especially `D:\Desktop\paper\bib`; do not edit that master bibliography here. Return citation cleanup to the `paper-citation-curator`.

Protocol: `docs/analysis/literature_to_implementation_protocol.md`.
