# Models Full Homomorphism Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make all physical model folders under `src/models` structurally and API-homomorphic in one round, and remove old model-path compatibility includes.

**Architecture:** Introduce canonical per-model anchors (`api.jl`, `capabilities.jl`, `core/`, `adapters/`, `workflows/`) for all 7 model kinds, then route `Models` entrypoints through these anchors. Capability checks become explicit and missing capabilities fail with `UnsupportedCapabilityError` instead of implicit method gaps.

**Tech Stack:** Julia 1.10+ (repo baseline), include-driven module organization, existing governance scripts and layered test runners.

---

### Task 1: Add Failing Homomorphism Governance Test

**Files:**
- Create: `tests/unit/models/test_model_structure_homomorphism.jl`
- Modify: `tests/unit/runtests.jl`
- Test: `tests/unit/models/test_model_structure_homomorphism.jl`

- [ ] **Step 1: Write the failing test**

```julia
using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

const MODEL_LAYOUT = Dict(
    :njl => joinpath(PROJECT_ROOT, "src", "models", "njl"),
    :njl2 => joinpath(PROJECT_ROOT, "src", "models", "njl2"),
    :pnjl => joinpath(PROJECT_ROOT, "src", "models", "pnjl"),
    :rpnjl => joinpath(PROJECT_ROOT, "src", "models", "rpnjl"),
    :pnjl_magnetic => joinpath(PROJECT_ROOT, "src", "models", "pnjl_magnetic"),
    :rotation => joinpath(PROJECT_ROOT, "src", "models", "rotation"),
    :gas_liquid => joinpath(PROJECT_ROOT, "src", "models", "gas_liquid"),
)

@testset "model structure homomorphism" begin
    for (kind, base) in MODEL_LAYOUT
        @test isdir(base)
        @test isfile(joinpath(base, "api.jl"))
        @test isfile(joinpath(base, "capabilities.jl"))
        @test isdir(joinpath(base, "core"))
        @test isdir(joinpath(base, "adapters"))
        @test isdir(joinpath(base, "workflows"))
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'include("tests/unit/models/test_model_structure_homomorphism.jl")'`
Expected: FAIL with missing directories/files for multiple model kinds.

- [ ] **Step 3: Register test in unit runner**

```julia
# add to unit model test include list
include(joinpath(@__DIR__, "models", "test_model_structure_homomorphism.jl"))
```

- [ ] **Step 4: Run unit selector to ensure runner wiring works**

Run: `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_model_structure_homomorphism.jl"; ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
Expected: FAIL for structure assertions (not include/wiring errors).

- [ ] **Step 5: Commit**

```bash
git add tests/unit/models/test_model_structure_homomorphism.jl tests/unit/runtests.jl
git commit -m "test(models): add failing homomorphism structure guard"
```

### Task 2: Add Failing API Homomorphism/Capability Test

**Files:**
- Create: `tests/unit/models/test_model_api_homomorphism.jl`
- Modify: `tests/unit/runtests.jl`
- Test: `tests/unit/models/test_model_api_homomorphism.jl`

- [ ] **Step 1: Write the failing test**

```julia
using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "model API homomorphism" begin
    for kind in (:NJL, :NJL2, :PNJL, :RPNJL, :PNJLMagnetic, :Rotation, :GasLiquid)
        model = Models.create_model(kind)
        caps = Models.model_capabilities(model)
        @test hasproperty(caps, :supports_solve_gap)
        @test hasproperty(caps, :supports_model_thermo)
        @test hasproperty(caps, :supports_number_densities)
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'include("tests/unit/models/test_model_api_homomorphism.jl")'`
Expected: FAIL because `model_capabilities` and/or capability schema does not exist yet.

- [ ] **Step 3: Register in unit runner**

```julia
include(joinpath(@__DIR__, "models", "test_model_api_homomorphism.jl"))
```

- [ ] **Step 4: Run selected unit test**

Run: `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_model_api_homomorphism.jl"; ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
Expected: FAIL at capability API assertions.

- [ ] **Step 5: Commit**

```bash
git add tests/unit/models/test_model_api_homomorphism.jl tests/unit/runtests.jl
git commit -m "test(models): add failing capability homomorphism contract"
```

### Task 3: Create Canonical Model Directories and Bridge Existing Implementations

**Files:**
- Create: `src/models/njl/api.jl`, `src/models/njl/capabilities.jl`, `src/models/njl/adapters/entrypoint_adapter.jl`, `src/models/njl/workflows/noop.jl`
- Create: `src/models/njl2/api.jl`, `src/models/njl2/capabilities.jl`, `src/models/njl2/adapters/entrypoint_adapter.jl`, `src/models/njl2/workflows/noop.jl`
- Create: `src/models/pnjl/api.jl`, `src/models/pnjl/capabilities.jl`, `src/models/pnjl/adapters/entrypoint_adapter.jl`, `src/models/pnjl/workflows/noop.jl`
- Create: `src/models/rpnjl/api.jl`, `src/models/rpnjl/capabilities.jl`, `src/models/rpnjl/adapters/entrypoint_adapter.jl`, `src/models/rpnjl/workflows/noop.jl`
- Create: `src/models/pnjl_magnetic/api.jl`, `src/models/pnjl_magnetic/capabilities.jl`, `src/models/pnjl_magnetic/adapters/entrypoint_adapter.jl`, `src/models/pnjl_magnetic/workflows/noop.jl`
- Create: `src/models/rotation/api.jl`, `src/models/rotation/capabilities.jl`, `src/models/rotation/adapters/entrypoint_adapter.jl`, `src/models/rotation/workflows/noop.jl`
- Create: `src/models/gas_liquid/api.jl`, `src/models/gas_liquid/capabilities.jl`, `src/models/gas_liquid/adapters/entrypoint_adapter.jl`, `src/models/gas_liquid/workflows/noop.jl`
- Modify: `src/models/Models.jl`
- Test: `tests/unit/models/test_model_structure_homomorphism.jl`

- [ ] **Step 1: Add generic capability error type and structs in `Models` domain**

```julia
struct UnsupportedCapabilityError <: Exception
    model_kind::Symbol
    capability::Symbol
end

Base.showerror(io::IO, e::UnsupportedCapabilityError) =
    print(io, "unsupported capability ", e.capability, " for model kind ", e.model_kind)

Base.@kwdef struct ModelCapabilities
    supports_solve_gap::Bool = true
    supports_model_thermo::Bool = true
    supports_number_densities::Bool = true
end
```

- [ ] **Step 2: Create per-model `capabilities.jl` files**

```julia
# example: src/models/njl/capabilities.jl
@inline njl_capabilities() = ModelCapabilities(
    supports_solve_gap=true,
    supports_model_thermo=true,
    supports_number_densities=true,
)
```

- [ ] **Step 3: Create per-model `api.jl` wrappers that delegate to existing implementation**

```julia
# example: src/models/njl/api.jl
@inline function njl_solve_gap(model, T, mu_vec; kwargs...)
    return solve_gap(model, T, mu_vec; kwargs...)
end

@inline function njl_number_densities(model, x_state, T, mu_vec; kwargs...)
    return number_densities(model, x_state, T, mu_vec; kwargs...)
end
```

- [ ] **Step 4: Add includes in `Models.jl` for all canonical anchors**

```julia
include(joinpath(@__DIR__, "njl", "capabilities.jl"))
include(joinpath(@__DIR__, "njl", "api.jl"))
# repeat for njl2/pnjl/rpnjl/pnjl_magnetic/rotation/gas_liquid
```

- [ ] **Step 5: Run structure homomorphism test and ensure green**

Run: `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_model_structure_homomorphism.jl"; ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/models/Models.jl src/models/njl src/models/njl2 src/models/pnjl src/models/rpnjl src/models/pnjl_magnetic src/models/rotation src/models/gas_liquid
git commit -m "refactor(models): add canonical homomorphic model anchors"
```

### Task 4: Implement Unified Capability Routing API

**Files:**
- Modify: `src/models/abstract_model.jl`
- Modify: `src/models/Models.jl`
- Modify: `src/models/entrypoints.jl`
- Test: `tests/unit/models/test_model_api_homomorphism.jl`

- [ ] **Step 1: Add `model_capabilities(model)` dispatches**

```julia
@inline model_capabilities(::NJLModel) = njl_capabilities()
@inline model_capabilities(::NJL2Model) = njl2_capabilities()
@inline model_capabilities(::PNJLModel) = pnjl_capabilities()
@inline model_capabilities(::RPNJLModel) = rpnjl_capabilities()
@inline model_capabilities(::PNJLMagneticModel) = pnjl_magnetic_capabilities()
@inline model_capabilities(::RotationModel) = rotation_capabilities()
@inline model_capabilities(::GasLiquidModel) = gas_liquid_capabilities()
```

- [ ] **Step 2: Enforce capability checks before unsupported contract use**

```julia
function require_capability(model, capability::Symbol)
    caps = model_capabilities(model)
    key = Symbol("supports_", capability)
    hasproperty(caps, key) || throw(UnsupportedCapabilityError(typeof(model).name.name, capability))
    getproperty(caps, key) || throw(UnsupportedCapabilityError(typeof(model).name.name, capability))
    return nothing
end
```

- [ ] **Step 3: Apply checks inside unified wrappers (`solve_gap`, `model_thermo`, `number_densities`)**

```julia
function solve_gap(model::AbstractQCDModel, T, mu_vec; kwargs...)
    require_capability(model, :solve_gap)
    return _solve_gap_impl(model, T, mu_vec; kwargs...)
end
```

- [ ] **Step 4: Export routing API for tests/governance**

```julia
export ModelCapabilities, UnsupportedCapabilityError, model_capabilities, require_capability
```

- [ ] **Step 5: Run API homomorphism test and ensure green**

Run: `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_model_api_homomorphism.jl"; ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/models/abstract_model.jl src/models/Models.jl src/models/entrypoints.jl tests/unit/models/test_model_api_homomorphism.jl
git commit -m "feat(models): add explicit capability routing for homomorphic model API"
```

### Task 5: Remove Old Model-Path Compatibility Includes Repo-Wide

**Files:**
- Modify: all scripts/tests with `include(... "src", "models", "<old path>")` references
- Modify: `scripts/dev/check_models_entry_contract.jl` (boundary warnings if needed)
- Test: targeted script help/smoke + integration selectors

- [ ] **Step 1: Replace old paths with canonical homomorphic paths**

```julia
# before
include(joinpath(PROJECT_ROOT, "src", "models", "workflows", "TransportWorkflow.jl"))

# after (example)
include(joinpath(PROJECT_ROOT, "src", "models", "workflow_apps", "TransportWorkflow.jl"))
```

- [ ] **Step 2: Remove temporary compatibility include stubs (if any were introduced during migration)**

```julia
# delete legacy include bridges that only forward old path -> new path
```

- [ ] **Step 3: Verify no old path references remain**

Run: `julia --project=. -e 'using Glob; println("manual path grep done in CI scripts")'`
Expected: old compatibility path patterns absent in repo code references.

- [ ] **Step 4: Commit**

```bash
git add src scripts tests
git commit -m "refactor(models): remove old model-path compatibility includes"
```

### Task 6: Governance and Layered Verification Sweep

**Files:**
- Modify: none (verification only unless failures)
- Test: governance + layered tests

- [ ] **Step 1: Run governance scripts**

Run:

```bash
julia --project=. scripts/dev/check_models_entry_contract.jl
julia --project=. scripts/dev/check_solver_contract_leakage.jl
julia --project=. scripts/dev/check_pnjl_migration_guard.jl
```

Expected: all `OK`.

- [ ] **Step 2: Run focused layered suites**

Run:

```bash
julia --project=. -e 'ENV["UNIT_PROFILE"]="core"; include("tests/unit/runtests.jl")'
julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="core"; include("tests/integration/runtests.jl")'
julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'
julia --project=. -e 'ENV["VALIDATION_PROFILE"]="smoke"; include("tests/validation/runtests.jl")'
```

Expected: PASS (or documented pre-existing skips only).

- [ ] **Step 3: Commit verification evidence update (if docs/issue notes changed)**

```bash
git add docs/dev tests
git commit -m "docs(models): record homomorphism verification evidence"
```

### Task 7: Final Closeout for #96

**Files:**
- Modify: issue #96 comments (via `gh issue comment`)
- Optional Modify: `docs/dev/active/*` if tracked execution log exists

- [ ] **Step 1: Post before/after model structure matrix to #96**

Run: `gh issue comment 96 --body "<structure matrix + capability matrix + verification summary>"`

- [ ] **Step 2: Post CI run evidence and unresolved debt (if any)**

Run: `gh run list --limit 20`
Expected: relevant workflows green.

- [ ] **Step 3: Decide closure state with user**

If all DoD met, propose close; otherwise list exact remaining blockers.

- [ ] **Step 4: Commit (only if local files changed)**

```bash
git add .
git commit -m "docs(models): finalize issue 96 homomorphism status"
```
