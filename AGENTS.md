# AGENTS.md

Repository guidance for `Julia_RelaxTime`. Default user-facing language: Chinese.

## Working Context

- Julia project, include-driven rather than a packaged `src/PackageName.jl` layout.
- Root environment: `Project.toml`; CI uses Julia `1.12.5`, root compatibility is Julia `1.10`.
- Use unified entrypoints in `Models` and `src/models/entrypoints.jl`. `src/pnjl/PNJL.jl` is historical compatibility context.
- Preserve mixed-meson governance and non-fixedmu unified joint-solve semantics unless explicitly in scope.
- QuadGK is an optional comparison oracle in the explicitly activated/stacked `benchmark/` environment. See [ADR 0006](docs/decisions/0006-isolate-optional-numerical-oracles.md).
- `AGENTS.md` is the shared repository instruction source; `.github/copilot-instructions.md` points here. No Cursor-specific rules are maintained.

## Task Scope and Records

- Read and validate in proportion to the change. Before editing `src/models/`, `src/relaxtime/` or `src/simulation/`, inspect the affected entrypoint and test layer.
- Inspect Git state and preserve unrelated work. Use one writer per file or logical block when work is delegated.
- Continue explicit tasks directly. For ambiguous continuations or track changes, run `scripts/dev/check_task_ledger.jl --preflight` and read the relevant items in `config/governance/task_tracks.toml`.
- Use the conversation or an existing issue/PR for ordinary planning. Maintain a task document only when it serves a shared design, user request or existing ledger contract. See [task tracking](docs/dev/task_tracking_governance.md).
- Update current summaries in place. Retain the inputs, code and key results needed to support scientific conclusions; use Git history for superseded tracked material. Repeated execution logs, source copies and cleanup receipts need not become permanent project assets.
- Search current source, scripts, tests, config, API/guides and active/backlog docs first. Read `docs/dev/archived/` for a specific historical question.
- Complete implementation, proportionate validation and a concise report. Existing authorization covers the corresponding steps; ask again only if scope or external impact materially changes.

## Correctness and Verification

- New core behavior needs meaningful coverage; reuse existing tests when sufficient. Prefer behavior, formulas, limits, units, error paths and numerical drift over source wording or file-count assertions.
- For solver, scan, cross-section and transport changes, select the affected unit/integration/regression/validation tests. Explain skipped numerical checks and their remaining risk.
- Change tolerances only with an explanation of the numerical error and evidence for the new criterion. An algorithm upgrade alone does not justify replacing a baseline.
- Validate public inputs with clear `ArgumentError` messages. Reuse `src/utils/ValidationUtils.jl`; keep finiteness, positivity, grid compatibility and recoverable warnings explicit.
- Justify performance changes with profiling or benchmarks. Keep smoke deterministic, fast and independent of external services.
- Verify effective behavior when changing entrypoints or configuration. Public API, CLI or data-contract changes need the corresponding `docs/api/` or `docs/guides/` update.

## Repository Layout

| Location | Responsibility |
| --- | --- |
| `src/models/`, `src/relaxtime/`, `src/simulation/` | Unified models, transport, HTTP service and simulation helpers |
| `config/models/<model>/`, `config/physics/` | Model configuration and shared physical constants |
| `tests/unit/`, `tests/integration/` | Focused behavior and cross-module workflows |
| `tests/regression/`, `tests/validation/` | Internal numerical baselines and external references |
| `benchmark/`, `scripts/perf/` | Performance comparisons/optional oracles and profiling |
| `scripts/analysis/` | Research and diagnostic scripts |
| `data/outputs/results/`, `data/outputs/figures/` | Numerical results/metadata and rendered figures |
| `docs/api/`, `docs/guides/`, `docs/analysis/` | Public contracts, usage and scientific interpretation |
| `docs/dev/active/`, `backlog/`, `archived/` | Shared ongoing work, deferred work and completed/superseded records |

## Julia Conventions

- Follow local style; use `PascalCase` for modules/types and `snake_case` for functions. Tests use `test_*.jl`; benchmarks use `bench_*.jl` or `benchmark_*.jl`.
- Put imports near the module top. Prefer explicit imports; use `import` to extend a method. Relative imports and `if !isdefined(Main, :ModuleName)` guards are normal in this include-driven project.
- Use concrete structs for stable inputs while preserving existing `NamedTuple` support. Prefer multiple dispatch; avoid `Any` and repeated allocations in hot paths. Add `@inline`/`@inbounds` with a demonstrated hot-path benefit.
- Internal values use natural units. New fm^-1 names use `_inv_fm`; MeV-facing values use names such as `T_MeV`. Follow existing higher-order forms such as `sigma_fm4` and `coupling_inv_fm4`. Domain Unicode names are supported.

## Scientific and Paper Work

- Check formulas, variables, units, implementation and actual evidence before writing a result or caption. Separate observations, conditional derivations and mechanism hypotheses.
- Use [research evidence guidance](docs/analysis/governance/research_evidence.md) for retained studies, and the relevant production/plotting SOP for formal deliveries.
- Document numerical API units, domains, missing/error values and ordering; map HTTP requests to implemented handlers before describing availability.
- Literature-to-implementation agents are available when explicitly requested; see [the protocol](docs/analysis/governance/literature_to_implementation_protocol.md). Citation/BibTeX curation belongs in `D:/Desktop/paper`, especially `paper/bib`.

## Reusable Guidance

- Keep a rule or helper where it has an ongoing consumer and reduces maintenance. Consolidate duplicated instructions; retire obsolete scripts with their dedicated tests.
- Repository skills in `.agents/skills/` are canonical. Keep them focused on project decisions and reusable operations; link shared rules instead of copying them.
- Generate `.claude/skills/` adapters with `scripts/dev/sync_claude_skills.jl`. `.claude/commands/` is for parameterized action macros.
- Commit subjects are concise, single-line descriptions of staged changes; familiar prefixes include `docs:`, `fix:`, `refactor:`, `ci:` and `feat:`.

## Command Reference

- [Environment, focused tests, governance and benchmarks](docs/dev/agent_command_reference.md)
- [Stable script entrypoints and sysimage wrappers](docs/guides/scripts/README.md)
- [Test layers and profile selection](docs/dev/testing_governance.md)

Start with focused selectors; expand checks to resolve a concrete remaining risk.
