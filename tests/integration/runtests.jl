# Integration test entrypoint
#
# These tests verify cross-module workflows, smoke tests, and
# multi-component interactions. They complement unit tests by
# testing at a higher level of integration.
#
# Run:
#   julia --project=. --eval 'include("tests/integration/runtests.jl")'
#
# Optional ENV knobs:
#   INTEGRATION_PROFILE=smoke|core|full
#     - smoke (default): ultra-fast core subset for edit-run loop
#     - core: previous broader smoke gate used by CI / pre-merge
#     - full: all test_*.jl across all integration subdirectories
#   INTEGRATION_FILES="models/test_a.jl;relaxtime/test_b.jl"
#     - runs only selected files (absolute path or relative to tests/integration)

using Test

const INTEGRATION_DIR = @__DIR__
const PROJECT_ROOT = normpath(joinpath(INTEGRATION_DIR, "..", ".."))

function _maybe_precompile_warmup()
    enabled = get(ENV, "TEST_PRECOMPILE_WARMUP", "0") in ("1", "true", "TRUE", "yes", "YES")
    enabled || return
    try
        if !isdefined(Main, :Models)
            Base.include(Main, joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
        end
        profile = Symbol(lowercase(get(ENV, "TEST_PRECOMPILE_PROFILE", "test")))
        Models.run_precompile_profile(profile)
    catch err
        @warn "Integration precompile warmup failed; continuing without warmup" exception=(err, catch_backtrace())
    end
end

function _warn_local_non_smoke(profile::String)
    is_ci = get(ENV, "CI", "") in ("1", "true", "TRUE", "yes", "YES")
    is_ci && return
    profile == "smoke" && return
    @warn "Local integration test run uses non-smoke profile; prefer smoke for edit-run loop" profile=profile recommended="INTEGRATION_PROFILE=smoke"
end

function _selected_integration_files()
    raw = strip(get(ENV, "INTEGRATION_FILES", ""))
    isempty(raw) && return nothing

    parts = split(replace(raw, ',' => ';'), ';'; keepempty=false)
    files = String[]
    for part in parts
        entry = strip(part)
        isempty(entry) && continue
        path = isabspath(entry) ? entry : joinpath(INTEGRATION_DIR, entry)
        isfile(path) || error("INTEGRATION_FILES entry does not exist: $(path)")
        push!(files, path)
    end

    isempty(files) && return nothing
    return files
end

function _should_include_integration_file(path::String)
    file = lowercase(basename(path))
    endswith(file, ".jl") || return false
    startswith(file, "test_") || return false

    # Exclude utility files (test_utils.jl is a helper, not a standalone test)
    file == "test_utils.jl" && return false

    return true
end

function _include_integration_dir(dir::String)
    isdir(dir) || return
    files = sort(readdir(dir; join=true))
    for f in files
        _should_include_integration_file(f) || continue
        rel = relpath(f, INTEGRATION_DIR)
        println("[integration-full] including: " * rel)
        try
            include(f)
        catch err
            println(stderr, "[integration] include failed: " * rel)
            rethrow(err)
        end
    end
end

# Ultra-fast subset for local edit-run loop (<~1 min target on warm cache).
# Full-only budget exclusions:
# - config/test_config_profile_smoke.jl: local subprocess/config load cost exceeded smoke budget.
# - relaxtime/test_plot_contract_smoke.jl: launches orchestrator Julia subprocesses.
# - relaxtime/test_cross_section_orchestrated_smoke.jl: orchestrated cross-section path is too slow for smoke.
# - relaxtime/test_xi_smoothness_batch_runner_smoke.jl: batch runner subprocess path is too slow for smoke.
# - relaxtime/test_phase_guided_transport_scan_smoke.jl: dry-run CLI subprocess path is too slow for smoke.
const INTEGRATION_SMOKE_FILES = [
    joinpath(INTEGRATION_DIR, "relaxtime", "test_xi_smoothness_plot_smoke.jl"),
    joinpath(INTEGRATION_DIR, "relaxtime", "test_xi_smoothness_evaluation_smoke.jl"),
]

# Broader gate retained for CI / pre-merge confidence.
const INTEGRATION_CORE_FILES = [
    # Models: 1 e2e + 1 contract guard
    joinpath(INTEGRATION_DIR, "models", "test_models_phase0_smoke.jl"),
    joinpath(INTEGRATION_DIR, "models", "test_models_dispatch_interface_smoke.jl"),
    joinpath(INTEGRATION_DIR, "models", "test_phase_pipeline_runner_smoke.jl"),

    # PNJL: 1 e2e + 1 contract guard
    joinpath(INTEGRATION_DIR, "pnjl", "test_tmu_scan_smoke.jl"),
    joinpath(INTEGRATION_DIR, "pnjl", "test_trho_scan_semantic_modes_smoke.jl"),

    # Variants: representative e2e
    joinpath(INTEGRATION_DIR, "models", "test_rotation_workflow_smoke.jl"),
    joinpath(INTEGRATION_DIR, "models", "test_gas_liquid_workflow_smoke.jl"),

    # Config: contract guard
    joinpath(INTEGRATION_DIR, "config", "test_config_profile_smoke.jl"),

    # RelaxTime: 1 e2e + 1 contract guard
    joinpath(INTEGRATION_DIR, "relaxtime", "test_transport_workflow_smoke.jl"),
    joinpath(INTEGRATION_DIR, "relaxtime", "test_orchestrator_smoke.jl"),
]

@testset "Integration" begin
    _maybe_precompile_warmup()

    selected = _selected_integration_files()

    if selected !== nothing
        @testset "Selected" begin
            for file in selected
                rel = relpath(file, INTEGRATION_DIR)
                println("[integration-selected] including: " * rel)
                include(file)
            end
        end
        return
    end

    profile = lowercase(get(ENV, "INTEGRATION_PROFILE", "smoke"))
    _warn_local_non_smoke(profile)

    if profile == "smoke"
        @testset "Smoke" begin
            for f in INTEGRATION_SMOKE_FILES
                if isfile(f)
                    include(f)
                else
                    @warn "Integration smoke file not found" path=f
                end
            end
        end
    elseif profile == "core"
        @testset "Core" begin
            for f in INTEGRATION_CORE_FILES
                if isfile(f)
                    include(f)
                else
                    @warn "Integration core file not found" path=f
                end
            end
        end
    elseif profile == "full"
        @testset "Models" begin
            _include_integration_dir(joinpath(INTEGRATION_DIR, "models"))
        end

        @testset "PNJL" begin
            _include_integration_dir(joinpath(INTEGRATION_DIR, "pnjl"))
        end

        @testset "RelaxTime" begin
            # Include test_utils.jl first (helper, not a test)
            utils_path = joinpath(INTEGRATION_DIR, "relaxtime", "test_utils.jl")
            if isfile(utils_path)
                include(utils_path)
            end
            _include_integration_dir(joinpath(INTEGRATION_DIR, "relaxtime"))
        end

        @testset "Config" begin
            _include_integration_dir(joinpath(INTEGRATION_DIR, "config"))
        end
    else
        error("Unknown INTEGRATION_PROFILE=$(profile). Use smoke, core, or full")
    end
end
