using Test

const REGRESSION_DIR = @__DIR__
const PROJECT_ROOT = normpath(joinpath(REGRESSION_DIR, "..", ".."))

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
        @warn "Regression precompile warmup failed; continuing without warmup" exception=(err, catch_backtrace())
    end
end

function _warn_local_non_smoke(profile::String)
    is_ci = get(ENV, "CI", "") in ("1", "true", "TRUE", "yes", "YES")
    is_ci && return
    profile == "smoke" && return
    @warn "Local regression test run uses non-smoke profile; prefer smoke for edit-run loop" profile=profile recommended="REGRESSION_PROFILE=smoke"
end

const EXPECTED_OPTIONAL_FIXTURE_SKIPS = [
    (
        name = "tau xi probe regression fixtures",
        file = "tests/regression/relaxtime/test_tau_xi_probe_regression.jl",
        reason = "缺少 data/outputs/results/relaxtime/scan/_xi_probe_T190_summary.csv 或 _xi_probe_T200_summary.csv 时会触发 @test_skip",
    ),
]

const SMOKE_FILES = [
    joinpath(REGRESSION_DIR, "relaxtime", "test_tau_xi_probe_regression.jl"),
]

const CORE_SMOKE_FILES = [
    # NJL: fixed-point guard
    joinpath(REGRESSION_DIR, "njl", "test_njl_gap_fixedpoint_regression.jl"),

    # RPNJL: fixed-point guard
    joinpath(REGRESSION_DIR, "rpnjl", "test_rpnjl_gap_fixedpoint_regression.jl"),

    # PNJL: fixed-point + semantic selection guard
    joinpath(REGRESSION_DIR, "pnjl", "test_scan_fixedpoint_regression.jl"),
    joinpath(REGRESSION_DIR, "pnjl", "test_constraint_selection_regression.jl"),

    # Models: one fixed-point + one contract semantic guard
    joinpath(REGRESSION_DIR, "models", "test_dimension_agnostic_solver_regression.jl"),
    joinpath(REGRESSION_DIR, "models", "test_solver_contract_regression.jl"),

    # RelaxTime: transport fixed-point + workflow consistency
    joinpath(REGRESSION_DIR, "relaxtime", "test_transport_fixedpoint_regression.jl"),
    joinpath(REGRESSION_DIR, "relaxtime", "test_t190_mixed_p_chain_regression.jl"),
    joinpath(REGRESSION_DIR, "relaxtime", "test_meson_density_plot_review_case_regression.jl"),
    joinpath(REGRESSION_DIR, "relaxtime", "test_tau_xi_probe_regression.jl"),
    joinpath(REGRESSION_DIR, "relaxtime", "test_workflow_vs_direct_consistency.jl"),
]

function _selected_regression_files()
    raw = strip(get(ENV, "REGRESSION_FILES", ""))
    isempty(raw) && return nothing

    parts = split(replace(raw, ',' => ';'), ';'; keepempty=false)
    files = String[]
    for part in parts
        entry = strip(part)
        isempty(entry) && continue
        path = isabspath(entry) ? entry : joinpath(REGRESSION_DIR, entry)
        isfile(path) || error("REGRESSION_FILES entry does not exist: $(path)")
        push!(files, path)
    end

    isempty(files) && return nothing
    return files
end

function _include_regression_dir(dir::String)
    isdir(dir) || return
    files = sort(readdir(dir; join=true))
    for file in files
        name = lowercase(basename(file))
        endswith(name, ".jl") || continue
        startswith(name, "test_") || continue
        include(file)
    end
end

function _print_expected_regression_skips()
    isempty(EXPECTED_OPTIONAL_FIXTURE_SKIPS) && return
    println("[regression] 预期可选跳过项（用于解释 summary 中的 Broken 计数）:")
    for item in EXPECTED_OPTIONAL_FIXTURE_SKIPS
        println("  - $(item.name)")
        println("    file: $(item.file)")
        println("    reason: $(item.reason)")
    end
end

@testset "Regression" begin
    _maybe_precompile_warmup()

    _print_expected_regression_skips()

    selected = _selected_regression_files()

    if selected !== nothing
        @testset "Selected" begin
            for file in selected
                include(file)
            end
        end
        return
    end

    profile = lowercase(get(ENV, "REGRESSION_PROFILE", "smoke"))
    _warn_local_non_smoke(profile)

    if profile == "smoke"
        @testset "Smoke" begin
            for file in SMOKE_FILES
                include(file)
            end
        end
    elseif profile == "core"
        @testset "Core" begin
            for file in CORE_SMOKE_FILES
                include(file)
            end
        end
    elseif profile == "full"
        @testset "NJL" begin
            _include_regression_dir(joinpath(REGRESSION_DIR, "njl"))
        end

        @testset "PNJL" begin
            _include_regression_dir(joinpath(REGRESSION_DIR, "pnjl"))
        end

        @testset "Models" begin
            _include_regression_dir(joinpath(REGRESSION_DIR, "models"))
        end

        @testset "RPNJL" begin
            _include_regression_dir(joinpath(REGRESSION_DIR, "rpnjl"))
        end

        @testset "RelaxTime" begin
            _include_regression_dir(joinpath(REGRESSION_DIR, "relaxtime"))
        end

        @testset "Phase" begin
            _include_regression_dir(joinpath(REGRESSION_DIR, "phase"))
        end
    else
        error("Unknown REGRESSION_PROFILE=$(profile). Use smoke, core, or full")
    end
end
