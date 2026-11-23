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

### Frontend tests

The default Node suite runs without a server; the live integration test is skipped:

```sh
node --test web/js/*.test.mjs
```

For live integration, start a dedicated disposable API server, then explicitly select it:

```powershell
$env:JRT_RUN_LIVE_API_TESTS = "1"
$env:JRT_TEST_API_BASE_URL = "http://127.0.0.1:8080"
node --test web/js/api.scan.integration.test.mjs
```

The live test submits a small scan and cancels the jobs it created in `finally`.
An already running numerical scan may continue until it returns, and output files can remain;
stop the disposable server after testing. Queue saturation is a separate opt-in
(`JRT_TEST_QUEUE_FULL=1`): it can submit up to 60 scans of 2,304 points each.
Use a small `PNJL_SCAN_MAX_PENDING` on that dedicated server and never target a shared research queue.

## Unified scan entrypoint

`scripts/models/run_unified_scan.jl` is the unified entrypoint for `Models` scans (`scan tmu`, `scan trho`, `scan magnetic`) and `workflow phase`; transport-coefficient scans keep using `scripts/relaxtime/run_gap_transport_scan.jl`. `--help` must be the first argument and prints the subcommand list:

```sh
julia --project=. scripts/models/run_unified_scan.jl --help
julia --project=. scripts/models/run_unified_scan.jl scan tmu --model_kind=PNJL --T_values=150 --mu_values=0,100 --xi_values=0.0 --output_path=data/outputs/results/tmu_smoke.csv --overwrite=true
julia --project=. scripts/models/run_unified_scan.jl scan magnetic --model_kind=PNJLMagnetic --solver_mode=fixed_mu --T_values=150 --mu_values=0 --eB_values=20000 --xi_values=0 --output_path=data/outputs/results/pnjl/scan/magnetic/selected.csv
```

Full parameter examples and the sysimage wrappers for these entrypoints are in the [stable script entrypoints guide](../guides/scripts/README.md).

## Governance and audit commands

There is no repository-wide JuliaFormatter or dedicated linter configuration. Preserve local style and use these checks as the nearest policy gate:

The script-entrypoint check parses Julia syntax without executing the scripts and checks their static relative references.

```sh
julia --project=. scripts/dev/check_unit_skip_policy.jl
julia --project=. scripts/dev/check_docs_consistency.jl
julia --project=. scripts/dev/check_active_docs_governance.jl
julia --project=. scripts/dev/check_sop_governance.jl
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
julia --project=. scripts/dev/check_task_ledger.jl --preflight --full-paths
```

The scheduled precompile check runs on Ubuntu with Julia 1.12.5. It caps the trace at
1100 lines without a sysimage and 450 with a fresh sysimage. Statements found only
in the trace without the sysimage, including the focused AD subset, are reported as
coverage diagnostics.

Task ledger preflight is read-only. It reports branch, HEAD, the dirty-path count and up
to 10 paths; `--full-paths` prints every path. A valid ledger also supplies the selected
track/current task, task file, blockers, next actions and up to 5 evidence entries.
Invalid ledgers still produce Git context and a failing exit status, without a trusted
task selection. `--track ID` selects an existing track and is also checked during default
validation. Recover explicit chat/file/issue targets first; use preflight for unresolved
track selection or coordination. It never changes the ledger or worktree.

Build and use a persistent sysimage:

```sh
julia --project=. scripts/dev/build_sysimage.jl
julia --sysimage=build/JuliaRelaxTime.dll --project=. -e 'ENV["INTEGRATION_PROFILE"]="core"; include("tests/integration/runtests.jl")'
```

## Benchmark commands

PNJL single-point/scan runs, scheduled trend comparisons, and result lookup: [benchmark trend guide](../guides/benchmark_trends.md).

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
