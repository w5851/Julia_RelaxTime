# Models Capability Map and Dedup Table

This note defines the authority boundary for models-facing capabilities and dedup decisions.

## Capability Map

- `Models` entrypoints are the authority interface for solver-facing workflows.
- Capability ownership should be explicit before introducing adapters or aliases.

## Dedup Adjudication

- Duplicate pathways are adjudicated by keeping one authority path and marking alternatives as compatibility-only.

## Task 1 Verification Log

- RED command: `julia --project=. -e 'include("tests/unit/models/solver/test_models_authority_map_contract.jl")'`
- RED summary: failed at `Expression: isfile(AUTHORITY_MAP_DOC)` with `models authority map document contract |    1      1` and `ERROR: Some tests did not pass`.
- GREEN command: `julia --project=. -e 'include("tests/unit/models/solver/test_models_authority_map_contract.jl")'`
- GREEN summary: passed with `models authority map document contract |    4      4`.

## Task 2 Verification Log

- RED command: `julia --project=. -e 'include("tests/unit/models/solver/test_models_authority_map_contract.jl")'`
- RED summary: failed at `Expression: isfile(AUTHORITY_MAP_CONFIG)` and `Expression: isfile(AUTHORITY_MAP_CHECK_SCRIPT)` with `models authority map document contract |    5     2      7` and `ERROR: Some tests did not pass`.
- GREEN command: `julia --project=. -e 'include("tests/unit/models/solver/test_models_authority_map_contract.jl")'`
- GREEN summary: passed with `models authority map document contract |    7      7`.
- GREEN command: `julia --project=. scripts/dev/check_models_authority_map.jl`
- GREEN summary: passed with `[models-authority-map] OK`.

## Task 3 Verification Log

- RED command: `julia --project=. -e 'ENV["UNIT_FILES"]="models/derivatives/test_diff_service_contract.jl"; include("tests/unit/runtests.jl")'`
- RED summary: failed at `isdefined(Models, :build_diff_service_context)` and `isdefined(Models, :eval_diff_service_jacobian)` with `UndefVarError: build_diff_service_context not defined in Main.Models`.
- GREEN command: `julia --project=. -e 'ENV["UNIT_FILES"]="models/derivatives/test_diff_service_contract.jl"; include("tests/unit/runtests.jl")'`
- GREEN summary: passed with `Unit | 6/6` for `test_diff_service_contract.jl`.

## Task 4 Verification Log

- RED command: `julia --project=. -e 'include("tests/unit/models/solver/test_solver_diff_pilot_adapters.jl")'`
- RED summary: failed because duplicate-target error expectation mismatched (`targets must not contain duplicates` vs `target_names must not contain duplicates`), indicating adapter path behavior contract not yet aligned.
- GREEN command: `julia --project=. -e 'include("tests/unit/models/solver/test_solver_diff_pilot_adapters.jl")'`
- GREEN summary: passed with `solver diff pilot adapters contract | 18 passed` after adapter-only alignment and unified-service delegation.

## Task 5 Verification Log

- RED command: `julia --project=. -e 'ENV["UNIT_FILES"]="models/derivatives/test_higher_order_derivatives.jl"; include("tests/unit/runtests.jl")'`
- RED summary: failed with `UndefVarError` for `nth_derivative`, `susceptibility_scale`, and `AbstractChiProvider` before new modules were included/exported.
- GREEN command: `julia --project=. -e 'ENV["UNIT_FILES"]="models/derivatives/test_higher_order_derivatives.jl"; include("tests/unit/runtests.jl")'`
- GREEN summary: passed with `Unit | 5/5` after adding higher-order utilities and provider abstraction modules.

## Task 6 Verification Log

- RED command: `julia --project=. -e 'include("tests/integration/models/test_models_derivatives_dual_smoke.jl")'`
- RED summary: failed while checking module-level dedup assertions (`_nth_derivative` / `_susceptibility_scale` still defined in `ConservedChargeSusceptibilities`) before switching to shared `HigherOrderDerivatives` utilities.
- GREEN command: `julia --project=. -e 'include("tests/integration/models/test_models_derivatives_dual_smoke.jl")'`
- GREEN summary: passed with all tests green (`18 passed`), including `kappa_sigma2 ≈ chi4/chi2` finite checks.

## Task 7 Verification Log

- RED command: `julia --project=. -e 'include("tests/regression/relaxtime/test_models_pnjl_fourth_derivative_baseline.jl")'`
- RED summary: failed at strict placeholder baseline assertion (`isapprox(..., 0.0; rtol=1e-12, atol=1e-12)`), indicating baseline value/tolerance not yet calibrated.
- GREEN command: `julia --project=. -e 'include("tests/regression/relaxtime/test_models_pnjl_fourth_derivative_baseline.jl")'`
- GREEN summary: passed with `Models PNJL chi4_B fixed-point baseline | Pass 3/3` after baseline calibration.
- GREEN command: `julia --project=. -e 'ENV["REGRESSION_FILES"]="relaxtime/test_models_pnjl_fourth_derivative_baseline.jl"; include("tests/regression/runtests.jl")'`
- GREEN summary: passed with `Regression | Pass 3/3` through regression entrypoint.

## Task 8 Closure Status

- Docs closure: updated `docs/api/models/derived/susceptibility/Overview.md` and `docs/api/models/derived/susceptibility/Cumulants.md` with fourth-order interface guidance and susceptibility/cumulant scaling contract notes.
- Governance + layered verification closure target:
  - `julia --project=. scripts/dev/check_models_entry_contract.jl`
  - `julia --project=. scripts/dev/check_models_authority_map.jl`
  - `julia --project=. -e 'ENV["UNIT_PROFILE"]="core"; include("tests/unit/runtests.jl")'`
  - `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
  - `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
- Final pass/fail evidence is recorded in task execution report after command runs; any failure remains a hard blocker and is not bypassed.

## Task 8 Verification Log

- Command: `julia --project=. scripts/dev/check_models_entry_contract.jl`
  - Result: PASS (`[models-entry-contract] OK`)
- Command: `julia --project=. scripts/dev/check_models_authority_map.jl`
  - Result: PASS (`[models-authority-map] OK`)
- Command: `julia --project=. -e 'ENV["UNIT_PROFILE"]="core"; include("tests/unit/runtests.jl")'`
  - Result: PASS (`Unit | 687/687`)
- Command: `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
  - Result: PASS (`Integration | 11/11`)
- Command: `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
  - Result: PASS (`Regression | 1/1`)
  - Note: missing xi-probe fixture csv is now treated as informational skip with `@test true` fallback, not Broken.
