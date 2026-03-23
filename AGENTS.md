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

## Setup Commands

Instantiate and precompile the main Julia environment:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

Run the standard package test entrypoint:

```sh
julia --project=. test/runtests.jl
```

Alternative package test command:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

Start the API plus web server:

```sh
julia --project=. scripts/server/server_full.jl
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

Unit smoke profile used by CI:

```sh
julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'
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

Integration smoke profile:

```sh
julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'
```

Integration full profile:

```sh
julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="full"; include("tests/integration/runtests.jl")'
```

Direct single integration test file:

```sh
julia --project=. -e 'include("tests/integration/relaxtime/test_transport_workflow_smoke.jl")'
```

Regression smoke profile:

```sh
julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'
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
julia --project=. scripts/dev/check_pnjl_migration_guard.jl
julia --project=. scripts/dev/analyze_deps.jl
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
