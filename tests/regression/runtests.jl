using Test

const REGRESSION_DIR = @__DIR__

const SMOKE_FILES = [
    joinpath(REGRESSION_DIR, "njl", "test_njl_gap_fixedpoint_regression.jl"),
    joinpath(REGRESSION_DIR, "rpnjl", "test_rpnjl_gap_fixedpoint_regression.jl"),
    joinpath(REGRESSION_DIR, "pnjl", "test_scan_fixedpoint_regression.jl"),
    joinpath(REGRESSION_DIR, "pnjl", "test_constraint_fixedpoint_regression.jl"),
    joinpath(REGRESSION_DIR, "pnjl", "test_magnetic_fixedpoint_regression.jl"),
    joinpath(REGRESSION_DIR, "relaxtime", "test_transport_fixedpoint_regression.jl"),
    joinpath(REGRESSION_DIR, "relaxtime", "test_tau_xi_probe_regression.jl"),
    joinpath(REGRESSION_DIR, "relaxtime", "test_total_cross_section_fixedpoint_regression.jl"),
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

@testset "Regression" begin
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
