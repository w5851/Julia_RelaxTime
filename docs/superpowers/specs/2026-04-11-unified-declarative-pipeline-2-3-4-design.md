# Unified Declarative Pipeline (2+3+4) Design

Date: 2026-04-11

Scope: `src/models/workflow/*`, `src/models/entrypoints.jl`, `src/models/workflows/*`, `scripts/models/run_unified_scan.jl`, `scripts/relaxtime/run_relaxtime_orchestrator.jl`, `config/workflows/relaxtime/*`

## 1. Background and Goal

Current repository has already landed a declarative phase pipeline around `Models.run_phase_pipeline` (`src/models/workflow/PipelineTypes.jl`, `src/models/workflow/PipelineRunner.jl`, `src/models/workflow/StageCatalog.jl`).

However, three high-value capability chains are still partially fragmented in orchestration style:

1. Models workflow entrypoints (`solve_gap_and_transport`, `solve_gap_and_meson_point`, rotation/gas-liquid workflow entrys)
2. Scan chains (`run_tmu_scan`, `run_trho_scan`, `scripts/models/run_unified_scan.jl`)
3. Relaxtime config-driven orchestrator (`config/workflows/relaxtime/*`, `scripts/relaxtime/run_relaxtime_orchestrator.jl`)

The design target is to unify these three chains under one declarative orchestration contract (without breaking current public APIs), so that numerical debugging/fixing and regression governance can rely on common manifests, stage records, and reproducibility metadata.

## 2. Non-Goals

- No large-scale rewrite of physics kernels in `src/relaxtime/*` or solver internals in `src/models/solver/*`
- No forced migration of all analysis/perf scripts into runtime pipeline in this wave
- No immediate removal of legacy script flags/aliases in this wave

## 3. Chosen Architecture (Option B)

Adopt a dual-layer kernel architecture:

- Keep `src/models/workflow/*` as the execution kernel (types, runner, stage validation, fail-fast semantics, manifest/hash)
- Add a domain adapter layer to map chains (2+3+4) to one shared stage semantics skeleton
- Keep script layer thin: parse CLI/profile only, then dispatch into stable `Models` entrypoints

Why this option:

- Better risk/cost balance than forcing hard unification immediately
- Preserves existing user-facing API contracts while unifying internal orchestration observability
- Provides direct foundation for relaxtime numerical issue triage via common manifests and regression tags

## 4. Logical Components and Responsibilities

### 4.1 Execution Kernel (existing, extended but stable)

- `PipelineTypes`: canonical pipeline contracts (`PipelineSpec`, `PipelineStage`, `PipelineContext`, `PipelineIOContract`)
- `PipelineRunner`: dependency validation, runtime execution, fail-fast behavior, skipped stage marking, manifest persistence
- Shared hash and reproducibility primitives (`config_hash`, `artifact_hash`, `run_id`, `git_commit`)

### 4.2 New Domain Adapter Layer

Recommended path: `src/models/workflow/adapters/`

- `WorkflowAdapter`: maps workflow entrypoints (transport/meson/rotation/gas-liquid) into canonical pipeline inputs/outputs
- `ScanAdapter`: maps scan entrypoints (`run_tmu_scan`, `run_trho_scan`, unified scan script) into canonical pipeline stages
- `RelaxtimeOrchestratorAdapter`: maps relaxtime toml-driven orchestrator into canonical pipeline stages

Adapters must only perform normalization and wiring, not embed heavy physics logic.

### 4.3 Stage Catalog by Domain

Recommended path: `src/models/workflow/catalog/`

Each chain provides a stage list with explicit `requires/provides`, reusing one semantic skeleton:

1. `prepare_inputs`
2. `solve_core`
3. `postprocess`
4. `export_artifacts`
5. `emit_repro_manifest`

### 4.4 Pipeline IO and Artifact Policy

Recommended path: `src/models/workflow/io/`

- Path strategy for artifacts and run directories
- Manifest schema consolidation and versioning
- Promotion/archive helper contracts reused by workflow/scan/relaxtime-orchestrator families

## 5. Unified Data Flow for Chains 2+3+4

For all three chains, use identical high-level flow:

1. `prepare_inputs`
   - Parse CLI/kwargs/profile/aliases
   - Produce canonical config payload and `config_hash`
2. `solve_core`
   - Invoke existing stable business kernels (`TransportWorkflow.solve_*`, `MesonMassWorkflow.solve_*`, `TmuScan.run_*`, orchestrated scan unit)
3. `postprocess`
   - Extract diagnostics and numerical health indicators
   - Attach baseline regression keys (for later fixed-point checks)
4. `export_artifacts`
   - Persist outputs and produce metadata index for generated files
5. `emit_repro_manifest`
   - Persist unified `run_manifest.json` with stage records and provenance fields

## 6. Compatibility Strategy

- Keep `src/models/entrypoints.jl` exported APIs stable; route internally through adapters + runner
- Keep existing script flags and parameter aliases in this wave; normalize aliases inside adapters
- Reuse existing alias governance from `config/workflows/relaxtime/schema/aliases_v1.toml`
- Migrate scripts to thin wrappers first; avoid behavior drift during wave transition

## 7. Error Model and Observability

### 7.1 Runtime Failure Semantics

- Maintain fail-fast semantics from current runner
- First failed stage recorded as `failed`; later pending stages marked `skipped`

### 7.2 Standardized Error Kind Taxonomy

Adapters and catalogs should normalize errors into:

- `input_validation_error`
- `numerical_convergence_error`
- `artifact_io_error`
- `unexpected_error`

### 7.3 Manifest Extension Fields

Extend pipeline-level manifest metadata with:

- `pipeline_family` (`workflow`, `scan`, `relaxtime_orchestrator`)
- `baseline_suite`
- `physics_profile`
- `adapter_version`

These fields are required for future numerical regression triage and CI gating.

## 8. Test Strategy

### 8.1 Unit

- Adapter normalization and alias mapping contracts
- Stage `requires/provides` completeness and uniqueness checks
- Manifest extension fields present and typed correctly

### 8.2 Integration (minimum closure)

- One workflow smoke route (transport or meson)
- One scan smoke route (tmu or trho)
- One relaxtime orchestrator smoke route (cross-section orchestrated path)

### 8.3 Regression

- Legacy route vs declarative route consistency tests on key observable fields
- Reuse and expand existing baseline-style consistency tests (including relaxtime workflow-vs-direct paths)
- Explicit tolerances must be codified in tests (rtol/atol)

## 9. Migration Waves and Gates

### Wave-1: Adapter + Catalog Extraction

- Add adapters and catalog declarations
- Keep old runtime path as default
- Validate no behavior drift via smoke + targeted consistency checks

### Wave-2: Models Entrypoint Switch

- Route relevant `Models` entrypoints to declarative orchestration path
- Keep fallback switch to old route for controlled rollback

### Wave-3: Script Thin-Wrapper Migration

- Convert `scripts/models/run_unified_scan.jl` and `scripts/relaxtime/run_relaxtime_orchestrator.jl` to thin dispatchers
- Keep CLI interface behavior compatible

### Wave-4: Regression Baseline and Docs Finalization

- Refresh baselines where needed under governance
- Update API docs for stable entrypoint behavior and declarative workflow contracts
- Close migration with governance checks and manifest schema confirmation

## 10. Definition of Done

The initiative is considered complete when all are true:

1. Workflow/Scan/Relaxtime-orchestrator chains all produce unified-structure manifests
2. At least one smoke integration path for each family passes
3. Legacy-vs-new key metrics stay within codified tolerances
4. Public entrypoint compatibility remains intact at `Models` boundary
5. Documentation reflects new declarative orchestration contracts

## 11. Risks and Mitigations

### Risk A: Hidden behavior drift during script migration

Mitigation: dual-path period + consistency regression tests before switching defaults.

### Risk B: Adapter layer accumulates business logic

Mitigation: enforce adapter-only normalization role; keep physics logic in existing domain modules.

### Risk C: Manifest schema fragmentation

Mitigation: single manifest schema versioning policy in workflow IO layer.

### Risk D: Numerical issue tracing still fragmented

Mitigation: standardize `pipeline_family`, baseline tags, and error taxonomy in all chain manifests.

## 12. Immediate Next Step

After this design is approved, generate an implementation plan that decomposes Wave-1..4 into test-first, bite-sized tasks and assigns exact files/tests/commands for each task.
