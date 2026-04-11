# Solver Domain Split and Declarative Pipeline Design

Date: 2026-04-10
Status: Proposed (approved in-session, pending implementation planning)
Scope: `src/models/solver`, `src/models/workflow`, `src/models/entrypoints.jl`, `scripts/pnjl/calculate_phase_structure.jl`

## 1. Context and Goal

Current `src/models/solver` has grown into a high-density directory with oversized files (e.g., `Solver.jl`, `ProblemSpecOrchestrator.jl`) and mixed responsibilities. Existing architecture already converges toward `ProblemSpec -> Orchestrator -> RootEngine`, but workflow orchestration is still scattered across entrypoints and scripts.

This design targets two outcomes:

1. Split solver by domain boundaries for high cohesion and low coupling.
2. Introduce a minimal declarative pipeline layer to orchestrate full research workflows (model build -> solve -> diagnostics -> analysis -> export -> reproducibility manifest).

Compatibility policy for this work is **breaking upgrade allowed**.

## 2. Design Decisions

### 2.1 Solver Domain Restructure (Wave-B compatible)

Restructure `src/models/solver` into subdomains, using move + re-export first (no algorithm changes in the first wave):

- `src/models/solver/api/`
  - `SolverAPI.jl`
  - Public entry facade only: `solve`, `solve_multi`, `solve_constraint`, `solve_vec`, `solve_named`
- `src/models/solver/spec/`
  - `ProblemSpec.jl`, `ConstraintComponents.jl`, `ConstraintModes.jl`, `Conditions.jl`
  - Constraint contracts and mode semantics
- `src/models/solver/orchestrator/`
  - `ProblemSpecOrchestrator.jl`, `PrimaryStrategy.jl`, `SeedStrategies.jl`
  - Attempt planning and fallback orchestration
- `src/models/solver/governance/`
  - `CandidateGovernance.jl`, `WeightedFallback.jl`
  - Candidate evaluation, hard rules, selectors
- `src/models/solver/runtime/`
  - `GenericRootEngine.jl`, `GapSolver.jl`, `ConstraintSolver*.jl`
  - Numerical solve kernels only
- `src/models/solver/diagnostics/`
  - `SolverDiagnostics.jl`, `SolverDiagnosticsTypes.jl`, `ThermoPostprocess.jl`
  - Diagnostics and postprocessing
- `src/models/solver/compat/`
  - `ImplicitAdapters.jl`, `ImplicitGapLegacy.jl`, `SchemaAdapter.jl`
  - Migration compatibility adapters
- `src/models/solver/config/`
  - `SolverRuntimeConfig.jl`, `StateSchema.jl`
  - Runtime and schema configuration

Boundary rules:

1. `api -> orchestrator/spec` allowed; reverse dependency forbidden.
2. `runtime` cannot depend on `api`.
3. `governance` cannot call `ConstraintSolver*` directly.
4. `compat` can only be referenced explicitly by `api/orchestrator`, not leaked into `runtime`.

Governance-runtime interaction contract:

- `orchestrator` is the only caller that can invoke both `runtime` and `governance`.
- `runtime` returns raw typed candidate payloads (solution, thermo, residual, convergence metadata).
- `governance` accepts candidate payloads plus hard-rule context, then returns normalized candidates and selection result.
- `runtime` must not import selector policy or hard-rule logic.
- `governance` must not call `ConstraintSolver*` or invoke numerical solves.

### 2.2 Declarative Workflow Minimal Contract

Add a workflow orchestration layer under `src/models/workflow/`:

- `PipelineTypes.jl`
  - `PipelineSpec`, `PipelineStage`, `PipelineContext`, `PipelineArtifact`
- `PipelineRunner.jl`
  - Stage dependency validation and ordered execution
- `StageCatalog.jl`
  - Standard stage registration

Minimal contract (linear stage chain for first release):

- `PipelineSpec`
  - `name`, `version`, `model_kind`, `stages::Vector{Symbol}`, `params::NamedTuple`, `io_contract`
- `PipelineStage`
  - `id::Symbol`, `requires::Vector{Symbol}`, `provides::Vector{Symbol}`, `run!::Function`
- `PipelineContext`
  - `state::Dict{Symbol,Any}` (orchestration-only dynamic container)
  - `provenance` (`git_commit`, `config_hash`, `run_id`, `timestamp`)
- `PipelineRunner.run(spec)`
  - Validates dependency chain
  - Executes stage sequence
  - Records per-stage artifact metadata (`path`, `hash`, `schema_version`)

Dependency and output validation rules:

- Cyclic dependencies are invalid and must raise `ArgumentError` before execution.
- If two stages provide the same symbol and no explicit merge policy is configured, runner must raise `ArgumentError` before execution.
- `provides` symbols must be unique across stage set in first release.
- Stage `requires` must be satisfied by initial context or predecessors in resolved order.

Stage execution interface (mandatory):

- Signature: `run!(ctx::PipelineContext, spec::PipelineSpec, stage::PipelineStage) -> StageResult`
- Stage input rule:
  - Read only from `ctx.state` keys listed in `stage.requires`.
  - Read-only access to `spec.params`.
- Stage output rule:
  - Must return `StageResult` with:
    - `produced::Dict{Symbol,Any}` (keys must be subset of `stage.provides`)
    - `artifacts::Vector{PipelineArtifact}`
    - `metrics::Dict{Symbol,Any}` (optional)
  - Runner merges `produced` into `ctx.state` after stage success.
- Side-effect rule:
  - External side effects allowed only in `export_artifacts` and `emit_repro_manifest` stages in first release.
  - Other stages should be pure with respect to filesystem/network.

First-stage catalog for full research flow:

1. `build_model`
2. `prepare_grid`
3. `solve_points`
4. `collect_diagnostics`
5. `analyze_phase`
6. `export_artifacts`
7. `emit_repro_manifest`

`io_contract` minimal definition (required in first implementation):

- Fields:
  - `contract_version::Symbol` (start with `:v1`)
  - `required_inputs::Vector{Symbol}`
  - `required_outputs::Vector{Symbol}`
  - `artifact_schema_version::Symbol` (start with `:artifact_v1`)
  - `manifest_schema_version::Symbol` (start with `:manifest_v1`)
- Validation:
  - `required_inputs` must be provided by initial context or previous stages.
  - `required_outputs` must exist in final context after last stage.
  - Duplicate stage `id` is invalid and must raise `ArgumentError` before execution.

Runner error semantics (required):

- Execution policy: fail-fast (default and only mode in first release).
- On stage failure:
  - Mark failing stage status as `:failed`.
  - Keep artifacts emitted by completed stages unchanged.
  - Do not execute downstream stages.
  - Return structured failure result with `failed_stage`, `error_kind`, `error_msg`, `completed_stages`.
  - Persist failure details into manifest stage record (`error_kind`, `error_msg`).
  - For downstream non-executed stages, runner must emit stage records with `status="skipped"` and null timestamps.
- CLI behavior on failure:
  - non-zero exit code.
  - print concise stage-local error summary.

Failure-manifest ownership rule:

- Runner owns manifest writing on both success and failure paths.
- `emit_repro_manifest` stage is treated as logical enrichment only; its downstream skip must not prevent manifest persistence.
- Minimum failure manifest must always include pipeline metadata, completed/failed/skipped stage records, and error summary.

### 2.3 Entrypoint and CLI Integration

- `Models.run_phase_pipeline` becomes a thin wrapper around `PipelineRunner`.
- `scripts/pnjl/calculate_phase_structure.jl` becomes: parse args -> build `PipelineSpec` -> run pipeline.
- Keep `Models` as the unified user-facing entrypoint.

CLI mapping contract (first release):

- Existing `PhaseCliConfig` fields map 1:1 into `PipelineSpec.params` under stable keys.
- `--config`, `--preset`, and explicit CLI overrides precedence is fixed as:
  1. base defaults
  2. model default config (`config/models/<model>/phase_pipeline_default.toml`, if present)
  3. user `--config`
  4. `--preset`
  5. explicit CLI flags
- Unknown CLI options remain hard errors.
- Backward behavior target: existing smoke preset command lines continue to run without extra required flags.

## 3. Migration Plan and Quality Gates

### Wave-1: Structural split only

- Move files into new subdirectories.
- Rebuild include order and re-export compatibility shell.
- No behavior change.

Gate:

- `tests/unit` smoke + core green.
- Old-to-new file path migration table committed and reviewed.

### Wave-2: API facade hardening

- Centralize public entry APIs in `solver/api/SolverAPI.jl`.
- Remove cross-layer direct dependencies violating boundary rules.

Gate:

- `tests/unit` + `tests/integration` core green.
- For stable public entrypoints changed at `Models` boundary, update `docs/api/` in same wave.

### Wave-3: Pipeline minimal contract

- Add `PipelineTypes`, `PipelineRunner`, `StageCatalog`.
- Route `Models.run_phase_pipeline` through `PipelineRunner` thin wrapper.

Gate:

- `tests/integration` + `tests/regression` smoke/core green.

### Wave-4: CLI migration and reproducibility outputs

- Migrate `scripts/pnjl/calculate_phase_structure.jl` to `PipelineSpec` construction.
- Emit `run_manifest` + artifact schema version metadata.

Gate:

- End-to-end script smoke.
- Fixed-point regression checks on key reference points.
- Manifest schema validation (`manifest_v1`) green.

## 4. Required New Tests

1. Stage contract tests (`requires/provides` completeness and type checks).
2. Runner behavior tests (missing dependency, stage failure, artifact logging).
3. Consistency tests (legacy `run_phase_pipeline` vs new pipeline, smoke-level key metrics within tolerance policy).
4. Reproducibility tests (same `PipelineSpec + config_hash` yields comparable manifest fields).

Reproducibility comparison matrix:

- Exact-match fields:
  - `pipeline_name`, `pipeline_version`, `model_kind`, `config_hash`, `artifact_schema_version`, `manifest_schema_version`
  - stage `id` sequence and stage `status` sequence
  - artifact relative filenames and per-artifact SHA-256
- Normalized-compare fields:
  - `artifact_paths` compared after path-root normalization to relative paths
  - `argv` compared as token list preserving order
- Ignored fields for reproducibility assertion:
  - `generated_at`, stage timestamps, `run_id`
  - `git_commit` only ignored when repo is dirty; exact-match otherwise

Tolerance policy for consistency tests:

- Scalar metrics (`pressure`, `omega`, `rho_norm`, `entropy`, `energy`): `rtol=1e-6`, `atol=1e-8`.
- Vector metrics (`x_state`, `mu_vec`, `masses`): elementwise `rtol=1e-6`, `atol=1e-8`.
- Categorical fields (`converged`, stage status, selection reason): exact match.

## 5. Risks and Controls

- Risk: refactor introduces subtle numerical path drift.
  - Control: wave split + regression gates + fixed-point baselines.
- Risk: orchestration layer leaks dynamic typing into hot kernels.
  - Control: keep `Dict{Symbol,Any}` limited to orchestration context; runtime kernels remain typed.
- Risk: compatibility code leaks back into core runtime.
  - Control: explicit dependency boundary checks in review and tests.

## 6. Acceptance Criteria

1. Solver directory matches target subdomain layout.
2. Public solve APIs are exposed through API facade only.
3. `run_phase_pipeline` executes via declarative stage runner.
4. CLI computes via `PipelineSpec` and emits reproducibility metadata.
5. Required tests pass at each migration wave gate.
6. Stage interface and manifest failure semantics are implemented exactly as specified.

## 7. Planning Boundary

To keep implementation planning estimable and executable, split into two linked plans:

1. Plan-A: solver domain split and boundary hardening (Wave-1 and Wave-2).
2. Plan-B: declarative pipeline runner + CLI/manifest migration (Wave-3 and Wave-4).

Integration checkpoint:

- After Plan-A completion, run integration tests and freeze API facade contract before starting Plan-B.

Plan handoff contract (Plan-A -> Plan-B):

- Frozen artifacts required before Plan-B starts:
  - `solver/api/SolverAPI.jl` public signatures and exported symbols.
  - include topology for `src/models/solver/*` subdomains.
  - boundary rule tests for forbidden dependencies.
  - migration mapping table finalized and committed.
- Plan-B may not change Plan-A frozen API signatures without explicit re-baseline and review.

## 8. Appendix A: Old-to-New Solver File Mapping (Initial)

- `src/models/solver/Solver.jl` -> `src/models/solver/api/SolverAPI.jl`
- `src/models/solver/ProblemSpec.jl` -> `src/models/solver/spec/ProblemSpec.jl`
- `src/models/solver/ConstraintComponents.jl` -> `src/models/solver/spec/ConstraintComponents.jl`
- `src/models/solver/ConstraintModes.jl` -> `src/models/solver/spec/ConstraintModes.jl`
- `src/models/solver/Conditions.jl` -> `src/models/solver/spec/Conditions.jl`
- `src/models/solver/ProblemSpecOrchestrator.jl` -> `src/models/solver/orchestrator/ProblemSpecOrchestrator.jl`
- `src/models/solver/PrimaryStrategy.jl` -> `src/models/solver/orchestrator/PrimaryStrategy.jl`
- `src/models/solver/SeedStrategies.jl` -> `src/models/solver/orchestrator/SeedStrategies.jl`
- `src/models/solver/CandidateGovernance.jl` -> `src/models/solver/governance/CandidateGovernance.jl`
- `src/models/solver/WeightedFallback.jl` -> `src/models/solver/governance/WeightedFallback.jl`
- `src/models/solver/GenericRootEngine.jl` -> `src/models/solver/runtime/GenericRootEngine.jl`
- `src/models/solver/GapSolver.jl` -> `src/models/solver/runtime/GapSolver.jl`
- `src/models/solver/ConstraintSolver.jl` -> `src/models/solver/runtime/ConstraintSolver.jl`
- `src/models/solver/ConstraintSolverCommon.jl` -> `src/models/solver/runtime/ConstraintSolverCommon.jl`
- `src/models/solver/ConstraintSolverFixedMu.jl` -> `src/models/solver/runtime/ConstraintSolverFixedMu.jl`
- `src/models/solver/ConstraintSolverFixedRho.jl` -> `src/models/solver/runtime/ConstraintSolverFixedRho.jl`
- `src/models/solver/ConstraintSolverFixedEntropy.jl` -> `src/models/solver/runtime/ConstraintSolverFixedEntropy.jl`
- `src/models/solver/ConstraintSolverFixedSigma.jl` -> `src/models/solver/runtime/ConstraintSolverFixedSigma.jl`
- `src/models/solver/ConstraintSolverFixedAsymmetricRho.jl` -> `src/models/solver/runtime/ConstraintSolverFixedAsymmetricRho.jl`
- `src/models/solver/SolverDiagnostics.jl` -> `src/models/solver/diagnostics/SolverDiagnostics.jl`
- `src/models/solver/SolverDiagnosticsTypes.jl` -> `src/models/solver/diagnostics/SolverDiagnosticsTypes.jl`
- `src/models/solver/ThermoPostprocess.jl` -> `src/models/solver/diagnostics/ThermoPostprocess.jl`
- `src/models/solver/ImplicitAdapters.jl` -> `src/models/solver/compat/ImplicitAdapters.jl`
- `src/models/solver/ImplicitGapLegacy.jl` -> `src/models/solver/compat/ImplicitGapLegacy.jl`
- `src/models/solver/SchemaAdapter.jl` -> `src/models/solver/compat/SchemaAdapter.jl`
- `src/models/solver/SolverRuntimeConfig.jl` -> `src/models/solver/config/SolverRuntimeConfig.jl`
- `src/models/solver/StateSchema.jl` -> `src/models/solver/config/StateSchema.jl`

## 9. Appendix B: `run_manifest` Schema (`manifest_v1`)

Required top-level keys:

- `generated_at::String`
- `git_commit::Union{String,Nothing}`
- `run_id::String`
- `pipeline_name::String`
- `pipeline_version::String`
- `model_kind::String`
- `config_hash::String`
- `argv::Vector{String}`
- `effective_config::Dict{String,Any}`
- `stages::Vector{Dict{String,Any}}` (each with `id`, `status`, `started_at`, `ended_at`)
- `artifact_paths::Dict{String,String}`
- `artifact_schema_version::String`
- `manifest_schema_version::String` (must equal `"manifest_v1"`)

Stage item required keys:

- `id::String`
- `status::String` (`"completed" | "failed" | "skipped"`)
- `started_at::Union{String,Nothing}` (UTC ISO-8601; null when `skipped`)
- `ended_at::Union{String,Nothing}` (UTC ISO-8601; null when `skipped`)
- `error_kind::Union{String,Nothing}`
- `error_msg::Union{String,Nothing}`

Hash and timestamp normalization:

- `config_hash` uses SHA-256 over canonical serialized effective config bytes.
- Artifact `hash` uses SHA-256 over file bytes.
- All timestamps in manifest use UTC ISO-8601 format.

## 10. Appendix C: Type Normalization Rules

Canonical in-memory types:

- `PipelineSpec.version::String`
- `io_contract.contract_version::Symbol`
- `io_contract.artifact_schema_version::Symbol`
- `io_contract.manifest_schema_version::Symbol`

Canonical persisted manifest types:

- `pipeline_version::String`
- `artifact_schema_version::String`
- `manifest_schema_version::String`

Normalization rule table:

- Persist: `PipelineSpec.version` is already canonical `String`
- Load: `String` -> `Symbol(string)` only for known enum-like fields
- Comparison tests:
  - in-memory compare on `Symbol`
  - manifest compare on normalized `String`
