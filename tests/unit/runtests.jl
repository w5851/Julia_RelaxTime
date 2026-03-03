# Unit test entrypoint
#
# Design goals:
# - Default: fast, deterministic, CI-friendly ("smoke" profile).
# - The repo currently contains a mix of true unit tests, scripts, perf probes, and legacy/WIP tests.
#   The entrypoint therefore defaults to a curated set, and offers opt-in broader runs.
#
# Run:
#   julia --project=. --eval 'include("tests/unit/runtests.jl")'
#
# Optional ENV knobs:
#   UNIT_INCLUDE_PERF=1   # also include files whose name contains "performance"
#   UNIT_INCLUDE_WIP=1    # include entries in DEFAULT_SKIP (for migration/debug only)
#   UNIT_INCLUDE_SLOW=1   # include slower/IO-touching smoke tests
#   UNIT_PROFILE=smoke|full
#     - smoke (default): curated, should be green
#     - full: include most test_*.jl for migration work; script/perf style tests should be moved out of unit

using Test

const UNIT_DIR = @__DIR__

# A small blacklist for tests that are currently WIP / outdated / intentionally non-unit.
# Rule of thumb:
# - keep true unit tests here only;
# - analysis/perf/script-style checks should move to tests/analysis, tests/perf, or scripts.
# Default `include("tests/unit/runtests.jl")` should be green and reasonably fast.
const DEFAULT_SKIP = Set([
    # phase-2: DEFAULT_SKIP 已清零；若需临时豁免，请同步更新 config/ci/unit_skip_policy.toml
])

const SMOKE_FILES = [
    # [Solver Robustness] 保留理由：高价值主链回归，覆盖随机采样但可确定复现。
    joinpath(UNIT_DIR, "pnjl", "test_solver_random_physical_smoke.jl"),

    # [Models Contracts] 保留理由：覆盖 models 子系统的最小可运行契约与后端切换。
    # Models phase-0 contract (solve_gap -> omega)
    joinpath(UNIT_DIR, "models", "test_models_phase0_smoke.jl"),

    # Models implicit differentiation wiring (NJL)
    joinpath(UNIT_DIR, "models", "test_models_implicitdiff_smoke.jl"),

    # Models legacy adapter (PNJL) wiring
    joinpath(UNIT_DIR, "models", "test_models_legacy_adapter_smoke.jl"),

    # Models legacy adapter (NJL) wiring
    joinpath(UNIT_DIR, "models", "test_models_legacy_njl_smoke.jl"),

    # Models PNJL minimal params injection (stage-2)
    joinpath(UNIT_DIR, "models", "test_pnjl_params_injection_smoke.jl"),
    joinpath(UNIT_DIR, "models", "test_pnjl_lambda_injection_smoke.jl"),
    joinpath(UNIT_DIR, "models", "test_pnjl_gk_polyakov_injection_smoke.jl"),
    joinpath(UNIT_DIR, "models", "test_pnjl_models_integrals_smoke.jl"),
    joinpath(UNIT_DIR, "models", "test_pnjl_integrals_forwarddiff_smoke.jl"),
    joinpath(UNIT_DIR, "models", "test_gap_residual_generic_smoke.jl"),
    joinpath(UNIT_DIR, "models", "test_pnjl_solve_gap_generic_smoke.jl"),
    joinpath(UNIT_DIR, "models", "test_models_dispatch_interface_smoke.jl"),
    joinpath(UNIT_DIR, "models", "test_models_unified_entrypoints_smoke.jl"),
    joinpath(UNIT_DIR, "models", "test_models_native_solver_phase1_smoke.jl"),
    joinpath(UNIT_DIR, "models", "test_models_derivatives_dual_smoke.jl"),
    joinpath(UNIT_DIR, "models", "test_phase_core_algorithms_smoke.jl"),
    joinpath(UNIT_DIR, "models", "test_phase_pipeline_smoke.jl"),
    joinpath(UNIT_DIR, "models", "test_phase_artifacts_promotion_smoke.jl"),
    joinpath(UNIT_DIR, "models", "test_phase_legacy_path_detach_smoke.jl"),
    joinpath(UNIT_DIR, "models", "test_pnjl_solve_gap_backend_switch_smoke.jl"),
    joinpath(UNIT_DIR, "models", "test_pnjl_thermo_bridge_multipoint_smoke.jl"),
    joinpath(UNIT_DIR, "models", "test_rpnjl_model_factory_smoke.jl"),

    # [Core Numerics] 保留理由：基础数值模块，变动少但影响范围大。
    # Core numerics / integration utils
    joinpath(UNIT_DIR, "integration", "test_gausslegendre.jl"),
    joinpath(UNIT_DIR, "integration", "test_cauchypv.jl"),

    # [RelaxTime Numerics] 保留理由：输运底层数值稳定性哨兵。
    # RelaxTime core numerics
    joinpath(UNIT_DIR, "relaxtime", "test_b0_correction.jl"),
    joinpath(UNIT_DIR, "relaxtime", "test_aniso_A_switch_smoke.jl"),

    # [Transport Workflow & Backend Bridge] 保留理由：阶段4/5主交付回归面。
    # Transport workflow (gap -> densities -> transport) wiring
    joinpath(UNIT_DIR, "relaxtime", "test_transport_workflow_smoke.jl"),
    joinpath(UNIT_DIR, "relaxtime", "test_transport_workflow_solver_backend_switch_smoke.jl"),
    joinpath(UNIT_DIR, "relaxtime", "test_transport_legacy_models_bridge_smoke.jl"),

    # [Workflow Cross-Checks] 保留理由：跨工作流交叉影响监控。
    # Meson mass workflow (gap -> meson mass/width -> Mott) wiring
    joinpath(UNIT_DIR, "relaxtime", "test_meson_mass_workflow_smoke.jl"),
    joinpath(UNIT_DIR, "pnjl", "test_tmu_scan_smoke.jl"),
    joinpath(UNIT_DIR, "pnjl", "test_tmu_scan_solver_backend_models_smoke.jl"),
    joinpath(UNIT_DIR, "pnjl", "test_trho_scan_smoke.jl"),
    joinpath(UNIT_DIR, "pnjl", "test_trho_scan_solver_backend_models_smoke.jl"),
    joinpath(UNIT_DIR, "pnjl", "test_constraint_modes_parity_smoke.jl"),
    joinpath(UNIT_DIR, "pnjl", "test_constraint_fixedpoint_baseline_smoke.jl"),
    joinpath(UNIT_DIR, "pnjl", "test_scan_fixedpoint_baseline_smoke.jl"),
    joinpath(UNIT_DIR, "pnjl", "test_magnetic_fixedpoint_baseline_smoke.jl"),
    joinpath(UNIT_DIR, "pnjl", "test_magnetic_nmax_convergence.jl"),
    joinpath(UNIT_DIR, "pnjl", "test_solver_constraints_models_backend_smoke.jl"),

    # [Config Injection] 保留理由：配置体系稳定性与可复现性保障。
    # Config profile selection/override rules
    joinpath(UNIT_DIR, "config", "test_config_profile_smoke.jl"),
    joinpath(UNIT_DIR, "config", "test_pnjl_profile_dynamic_constants_smoke.jl"),
    joinpath(UNIT_DIR, "config", "test_pnjl_rpnjl_critical_params_smoke.jl"),
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

            include_slow = get(ENV, "UNIT_INCLUDE_SLOW", "0") in ("1", "true", "TRUE", "yes", "YES")
            if include_slow
                include(joinpath(UNIT_DIR, "models", "test_rpnjl_bridge_smoke.jl"))
                include(joinpath(UNIT_DIR, "relaxtime", "test_transport_workflow_toml_prefer_energy_aniso_smoke.jl"))
            end
        end
    elseif profile == "full"
        @testset "Integration" begin
            _include_dir(joinpath(UNIT_DIR, "integration"))
        end

        @testset "PNJL" begin
            _include_dir(joinpath(UNIT_DIR, "pnjl"))
        end

        @testset "RelaxTime" begin
            _include_dir(joinpath(UNIT_DIR, "relaxtime"))
        end
    else
        error("Unknown UNIT_PROFILE=$(profile). Use UNIT_PROFILE=smoke or UNIT_PROFILE=full")
    end
end
