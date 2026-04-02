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
#   INTEGRATION_PROFILE=smoke|full
#     - smoke (default): curated fast subset
#     - full: all test_*.jl across all integration subdirectories

using Test

const INTEGRATION_DIR = @__DIR__

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
        try
            include(f)
        catch err
            println(stderr, "[integration] include failed: " * rel)
            rethrow(err)
        end
    end
end

# Curated smoke subset — fast-running integration tests for CI.
const INTEGRATION_SMOKE_FILES = [
    # Models contracts
    joinpath(INTEGRATION_DIR, "models", "test_models_phase0_smoke.jl"),
    joinpath(INTEGRATION_DIR, "models", "test_models_dispatch_interface_smoke.jl"),
    joinpath(INTEGRATION_DIR, "models", "test_models_implicitdiff_flavor_mu_smoke.jl"),
    joinpath(INTEGRATION_DIR, "models", "test_rpnjl_model_factory_smoke.jl"),
    joinpath(INTEGRATION_DIR, "models", "test_phase_cli_smoke.jl"),
    joinpath(INTEGRATION_DIR, "models", "test_solver_backend_semantic_parity_guard.jl"),

    # PNJL scans
    joinpath(INTEGRATION_DIR, "pnjl", "test_tmu_scan_smoke.jl"),
    joinpath(INTEGRATION_DIR, "pnjl", "test_conserved_charge_susceptibilities_smoke.jl"),
    joinpath(INTEGRATION_DIR, "pnjl", "test_solver_random_physical_smoke.jl"),

    # Rotation workflow
    joinpath(INTEGRATION_DIR, "models", "test_rotation_workflow_smoke.jl"),

    # Config
    joinpath(INTEGRATION_DIR, "config", "test_config_profile_smoke.jl"),

    # RelaxTime workflows
    joinpath(INTEGRATION_DIR, "relaxtime", "test_transport_workflow_smoke.jl"),
    joinpath(INTEGRATION_DIR, "relaxtime", "test_meson_mass_workflow_smoke.jl"),
    joinpath(INTEGRATION_DIR, "relaxtime", "test_orchestrator_smoke.jl"),
    joinpath(INTEGRATION_DIR, "relaxtime", "test_cross_section_orchestrated_smoke.jl"),
    joinpath(INTEGRATION_DIR, "relaxtime", "test_plot_contract_smoke.jl"),
    joinpath(INTEGRATION_DIR, "relaxtime", "test_task_center_point_mode_smoke.jl"),

    # Gas-Liquid workflow
    joinpath(INTEGRATION_DIR, "models", "test_gas_liquid_workflow_smoke.jl"),
]

@testset "Integration" begin
    profile = lowercase(get(ENV, "INTEGRATION_PROFILE", "smoke"))

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
        error("Unknown INTEGRATION_PROFILE=$(profile). Use smoke or full")
    end
end
