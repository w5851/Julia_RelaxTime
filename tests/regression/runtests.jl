using Test

const REGRESSION_DIR = @__DIR__

const EXPECTED_OPTIONAL_FIXTURE_SKIPS = [
    (
        name = "tau xi probe regression fixtures",
        file = "tests/regression/relaxtime/test_tau_xi_probe_regression.jl",
        reason = "缺少 data/outputs/results/relaxtime/scan/_xi_probe_T190_summary.csv 或 _xi_probe_T200_summary.csv 时会触发 @test_skip",
    ),
]

const SMOKE_FILES = [
    joinpath(REGRESSION_DIR, "njl", "test_njl_gap_fixedpoint_regression.jl"),
    joinpath(REGRESSION_DIR, "rpnjl", "test_rpnjl_gap_fixedpoint_regression.jl"),
    joinpath(REGRESSION_DIR, "pnjl", "test_scan_fixedpoint_regression.jl"),
    joinpath(REGRESSION_DIR, "pnjl", "test_constraint_fixedpoint_regression.jl"),
    joinpath(REGRESSION_DIR, "pnjl", "test_constraint_selection_regression.jl"),
    joinpath(REGRESSION_DIR, "pnjl", "test_magnetic_fixedpoint_regression.jl"),
    joinpath(REGRESSION_DIR, "models", "test_dimension_agnostic_solver_regression.jl"),
    joinpath(REGRESSION_DIR, "models", "test_fixedrho_precision_guard_regression.jl"),
    joinpath(REGRESSION_DIR, "models", "test_problem_spec_fixedrho_parity_regression.jl"),
    joinpath(REGRESSION_DIR, "models", "test_solver_attempt_engine_convergence_regression.jl"),
    joinpath(REGRESSION_DIR, "models", "test_solver_diagnostic_exception_regression.jl"),
    joinpath(REGRESSION_DIR, "models", "test_solver_contract_regression.jl"),
    joinpath(REGRESSION_DIR, "models", "test_solver_phase3_fixedpoint_regression.jl"),
    joinpath(REGRESSION_DIR, "models", "test_fixedrho_semantic_equivalence_regression.jl"),
    joinpath(REGRESSION_DIR, "models", "test_firstorder_manifold_branch_stability.jl"),
    joinpath(REGRESSION_DIR, "relaxtime", "test_transport_fixedpoint_regression.jl"),
    joinpath(REGRESSION_DIR, "relaxtime", "test_tau_xi_probe_regression.jl"),
    joinpath(REGRESSION_DIR, "relaxtime", "test_total_cross_section_fixedpoint_regression.jl"),
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

    if profile == "smoke"
        @testset "Smoke" begin
            for file in SMOKE_FILES
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
        error("Unknown REGRESSION_PROFILE=$(profile). Use smoke or full")
    end
end
