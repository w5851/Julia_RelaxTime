# Unit test entrypoint
#
# Design goals:
# - Default: fast, deterministic, CI-friendly ("smoke" profile).
# - Each test_*.jl maps 1:1 to a standalone src module.
# - Integration/smoke/workflow tests live under tests/integration/.
# - Baseline/validation tests live under tests/validation/.
#
# Run:
#   julia --project=. --eval 'include("tests/unit/runtests.jl")'
#
# Optional ENV knobs:
#   UNIT_INCLUDE_PERF=1   # also include files whose name contains "performance"
#   UNIT_INCLUDE_WIP=1    # include entries in DEFAULT_SKIP (for migration/debug only)
#   UNIT_PROFILE=smoke|full
#     - smoke (default): curated fast subset, should be green
#     - full: all test_*.jl across all unit subdirectories

using Test

const UNIT_DIR = @__DIR__

# Blacklist for tests temporarily excluded.
const DEFAULT_SKIP = Set([
    # (currently empty — all tests should pass)
])

# Curated fast subset for CI / daily smoke runs.
# Only references files that exist under tests/unit/.
const SMOKE_FILES = [
    # [Core Numerics] 基础数值模块
    joinpath(UNIT_DIR, "numerics", "test_gausslegendre.jl"),
    joinpath(UNIT_DIR, "numerics", "test_cauchypv.jl"),
    joinpath(UNIT_DIR, "numerics", "test_phase_space_sampling.jl"),

    # [Types & Config] 基础类型与配置
    joinpath(UNIT_DIR, "types", "test_parameter_types.jl"),
    joinpath(UNIT_DIR, "config", "test_config_loader.jl"),

    # [NJL Core] NJL 模型核心
    joinpath(UNIT_DIR, "njl", "test_njl_core.jl"),
    joinpath(UNIT_DIR, "njl", "test_njl2_core.jl"),

    # [Models] 模型子系统（工厂/配置/工具）
    joinpath(UNIT_DIR, "models", "test_njl_model_factory.jl"),
    joinpath(UNIT_DIR, "models", "test_scan_config.jl"),
    joinpath(UNIT_DIR, "models", "test_adaptive_rho_refinement.jl"),
    joinpath(UNIT_DIR, "models", "test_workflow_param_adapters.jl"),

    # [PNJL Solver] 求解器核心
    joinpath(UNIT_DIR, "pnjl", "test_solver_conditions.jl"),
    joinpath(UNIT_DIR, "pnjl", "test_solver_seed_strategies.jl"),
    joinpath(UNIT_DIR, "pnjl", "test_solver_constraint_modes.jl"),
    joinpath(UNIT_DIR, "pnjl", "test_pnjl_core.jl"),
    joinpath(UNIT_DIR, "pnjl", "test_core_integrals.jl"),

    # [RelaxTime Core] 输运核心模块
    joinpath(UNIT_DIR, "relaxtime", "test_oneloopintegrals.jl"),
    joinpath(UNIT_DIR, "relaxtime", "test_effective_couplings.jl"),
    joinpath(UNIT_DIR, "relaxtime", "test_meson_mass.jl"),
    joinpath(UNIT_DIR, "relaxtime", "test_mott_transition.jl"),
    joinpath(UNIT_DIR, "relaxtime", "test_particle_symbols.jl"),
    joinpath(UNIT_DIR, "relaxtime", "test_afieldbuilder.jl"),

    # [Simulation] 模拟子系统
    joinpath(UNIT_DIR, "simulation", "test_frame_transformations.jl"),
    joinpath(UNIT_DIR, "simulation", "test_momentum_mapping.jl"),
]

function _selected_unit_files()
    raw = strip(get(ENV, "UNIT_FILES", ""))
    isempty(raw) && return nothing

    parts = split(replace(raw, ',' => ';'), ';'; keepempty=false)
    files = String[]
    for p in parts
        f = strip(p)
        isempty(f) && continue
        path = isabspath(f) ? f : joinpath(UNIT_DIR, f)
        isfile(path) || error("UNIT_FILES entry does not exist: $(path)")
        push!(files, path)
    end
    isempty(files) && return nothing
    return files
end

function _should_include_unit_file(path::String)
    file = lowercase(basename(path))

    endswith(file, ".jl") || return false
    startswith(file, "test_") || return false

    include_all = get(ENV, "UNIT_INCLUDE_ALL", "0") in ("1", "true", "TRUE", "yes", "YES")
    include_wip = get(ENV, "UNIT_INCLUDE_WIP", "0") in ("1", "true", "TRUE", "yes", "YES")

    if !include_all
        if (file in DEFAULT_SKIP) && !include_wip
            return false
        end
    end

    # Exclude performance-style tests by default.
    include_perf = get(ENV, "UNIT_INCLUDE_PERF", "0") in ("1", "true", "TRUE", "yes", "YES")
    if occursin("performance", file) && !include_perf
        return false
    end

    return true
end

function _include_dir(dir::String)
    files = sort(readdir(dir; join=true))
    verbose_include = get(ENV, "UNIT_VERBOSE_INCLUDE", "0") in ("1", "true", "TRUE", "yes", "YES")
    for f in files
        _should_include_unit_file(f) || continue
        rel = relpath(f, UNIT_DIR)
        if verbose_include
            println("[unit-full] including: " * rel)
        end
        try
            include(f)
        catch err
            println(stderr, "[unit-full] include failed: " * rel)
            println(stderr, "[unit-full] reproduce: UNIT_PROFILE=full UNIT_FILES='" * replace(rel, '\\' => '/') * "' julia --project=. tests/unit/runtests.jl")
            rethrow(err)
        end
    end
end

@testset "Unit" begin
    selected = _selected_unit_files()

    if selected !== nothing
        @testset "Selected" begin
            for f in selected
                include(f)
            end
        end
        return
    end

    profile = lowercase(get(ENV, "UNIT_PROFILE", "smoke"))

    if profile == "smoke"
        @testset "Smoke" begin
            for f in SMOKE_FILES
                include(f)
            end
        end
    elseif profile == "full"
        @testset "Types" begin
            _include_dir(joinpath(UNIT_DIR, "types"))
        end

        @testset "Config" begin
            _include_dir(joinpath(UNIT_DIR, "config"))
        end

        @testset "Numerics" begin
            _include_dir(joinpath(UNIT_DIR, "numerics"))
        end

        @testset "NJL" begin
            _include_dir(joinpath(UNIT_DIR, "njl"))
        end

        @testset "Models" begin
            _include_dir(joinpath(UNIT_DIR, "models"))
        end

        @testset "PNJL" begin
            _include_dir(joinpath(UNIT_DIR, "pnjl"))
        end

        @testset "RelaxTime" begin
            _include_dir(joinpath(UNIT_DIR, "relaxtime"))
        end

        @testset "Simulation" begin
            _include_dir(joinpath(UNIT_DIR, "simulation"))
        end
    else
        error("Unknown UNIT_PROFILE=$(profile). Use UNIT_PROFILE=smoke or UNIT_PROFILE=full")
    end
end
