# Agent command reference

This page is the detailed command catalog referenced by the root `AGENTS.md`. It keeps frequently changing command syntax out of the always-loaded instruction file while preserving one reviewable repository source.

## Environment and setup

Instantiate and precompile the main environment:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

Run the standard package test entrypoint or the layered wrappers:

```sh
julia --project=. test/runtests.jl
julia --project=. test/unit.jl
julia --project=. test/integration.jl
julia --project=. test/regression.jl
julia --project=. test/validation.jl
julia --project=. -e 'using Pkg; Pkg.test()'
```

Start the API and web server:

```sh
julia --project=. scripts/server/server_full.jl
```

Prefer the stable sysimage wrappers for documented CLI workflows:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 scripts/pnjl/calculate_phase_structure.jl --preset=smoke
powershell -ExecutionPolicy Bypass -File scripts/dev/bootstrap_sysimage.ps1
```

```sh
sh scripts/dev/run_with_sysimage.sh scripts/pnjl/calculate_phase_structure.jl --preset=smoke
sh scripts/dev/bootstrap_sysimage.sh
```

Run phase CLI warmup controls:

```sh
julia --project=. scripts/pnjl/calculate_phase_structure.jl --preset=smoke
julia --project=. -e 'ENV["PHASE_PRECOMPILE_WARMUP"]="0"; include("scripts/pnjl/calculate_phase_structure.jl")'
julia --project=. -e 'ENV["PHASE_PRECOMPILE_PROFILE"]="core"; include("scripts/pnjl/calculate_phase_structure.jl")'
```

Install and use dependency-graph rendering tools:

```sh
npm install
npm run deps:render
julia --project=. scripts/dev/gen_deps.jl
```

## Focused test profiles

Unit profiles:

```sh
julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'
julia --project=. -e 'ENV["UNIT_PROFILE"]="core"; include("tests/unit/runtests.jl")'
julia --project=. -e 'ENV["UNIT_PROFILE"]="full"; include("tests/unit/runtests.jl")'
julia --project=. -e 'ENV["UNIT_FILES"]="relaxtime/test_transport_coefficients.jl"; include("tests/unit/runtests.jl")'
julia --project=. -e 'include("tests/unit/relaxtime/test_transport_coefficients.jl")'
```

Integration profiles:

```sh
julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'
julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="core"; include("tests/integration/runtests.jl")'
julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="full"; include("tests/integration/runtests.jl")'
julia --project=. -e 'ENV["INTEGRATION_FILES"]="relaxtime/test_transport_workflow_smoke.jl"; include("tests/integration/runtests.jl")'
julia --project=. -e 'include("tests/integration/relaxtime/test_transport_workflow_smoke.jl")'
```

Regression profiles:

```sh
julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'
julia --project=. -e 'ENV["REGRESSION_PROFILE"]="core"; include("tests/regression/runtests.jl")'
julia --project=. -e 'ENV["REGRESSION_PROFILE"]="full"; include("tests/regression/runtests.jl")'
julia --project=. -e 'ENV["REGRESSION_FILES"]="relaxtime/test_transport_fixedpoint_regression.jl"; include("tests/regression/runtests.jl")'
```

Validation profiles:

```sh
julia --project=. tests/validation/runtests.jl
julia --project=. -e 'ENV["VALIDATION_PROFILE"]="smoke"; include("tests/validation/runtests.jl")'
julia --project=. -e 'ENV["VALIDATION_PROFILE"]="core"; include("tests/validation/runtests.jl")'
julia --project=. -e 'ENV["VALIDATION_FILES"]="relaxtime/test_mott_reference_mapping.jl"; include("tests/validation/runtests.jl")'
```

Optional test precompile warmup and VS Code/TestItems entrypoint:

```sh
julia --project=. -e 'ENV["TEST_PRECOMPILE_WARMUP"]="1"; ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'
julia --project=. -e 'ENV["TEST_PRECOMPILE_WARMUP"]="1"; ENV["TEST_PRECOMPILE_PROFILE"]="core"; ENV["INTEGRATION_PROFILE"]="core"; include("tests/integration/runtests.jl")'
julia --project=. test/vscode_testitems.jl
```

Local non-smoke profiles intentionally emit a smoke-first warning.

## Governance and audit commands

There is no repository-wide JuliaFormatter or dedicated linter configuration. Preserve local style and use these checks as the nearest policy gate:

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
julia --project=. scripts/dev/check_dependency_policy.jl
julia --project=. scripts/dev/check_agent_instructions.jl
julia --project=. scripts/dev/check_model_profile_matrix.jl
julia --project=. scripts/dev/check_data_output_path_guard.jl
julia --project=. scripts/dev/check_precompile_profile_coverage.jl
julia --project=. scripts/dev/check_precompile_trace_budget.jl
julia --project=. scripts/dev/check_skill_governance.jl
julia --project=. scripts/dev/sync_claude_skills.jl --check
julia --project=. scripts/dev/check_task_ledger.jl
julia --project=. scripts/dev/check_task_ledger.jl --preflight
julia --project=. scripts/dev/check_task_ledger.jl --preflight --track rs-transport
```

Task ledger preflight is read-only. It reports the selected primary track, current task,
branch, HEAD and dirty paths; it never changes the ledger or worktree. `--track ID` selects
an existing track and is also checked during the default validation.

Build and use a persistent sysimage:

```sh
julia --project=. scripts/dev/build_sysimage.jl
julia --sysimage=build/JuliaRelaxTime.dll --project=. -e 'ENV["INTEGRATION_PROFILE"]="core"; include("tests/integration/runtests.jl")'
```

## Benchmark commands

The benchmark environment is isolated from the normal root runtime/test environment. Optional comparison oracles such as QuadGK are declared only in `benchmark/Project.toml`. Instantiate that environment first. Because this repository is include-driven, benchmark scripts that also need root dependencies use an explicit per-process environment stack.

```sh
julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'
```

Windows / PowerShell:

```powershell
$env:JULIA_LOAD_PATH = "@;benchmark;@stdlib"
julia --project=. benchmark/relaxtime/benchmark_quadgk_oracle_smoke.jl
```

Linux / macOS:

```sh
JULIA_LOAD_PATH='@:benchmark:@stdlib' julia --project=. benchmark/relaxtime/benchmark_quadgk_oracle_smoke.jl
```

Use `benchmark/` for comparative performance results and `scripts/perf/` for focused profiling probes. Smoke tests must not contain performance thresholds.
