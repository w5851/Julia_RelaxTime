# Validation / baseline regression test entrypoint
#
# These tests compare numerical results against stored baselines
# to detect unintended numerical drift across versions.
#
# Run:
#   julia --project=. --eval 'include("tests/validation/runtests.jl")'
#
# Optional ENV knobs:
#   VALIDATION_PROFILE=smoke|core|full
#     - smoke (default): ultra-fast subset for edit-run loop
#     - core: representative validation subset for pre-merge confidence
#     - full: all test_*.jl under tests/validation/
#   VALIDATION_FILES="pnjl/reference/test_crossover_legacy_source_consistency.jl;relaxtime/test_mott_reference_mapping.jl"
#     - runs only selected files (absolute path or relative to tests/validation)

using Test

const VALIDATION_DIR = @__DIR__

include(joinpath(VALIDATION_DIR, "common", "data_paths.jl"))

const VALIDATION_SMOKE_FILES = [
    joinpath(VALIDATION_DIR, "pnjl", "reference", "test_crossover_legacy_source_consistency.jl"),
    joinpath(VALIDATION_DIR, "relaxtime", "test_mott_reference_mapping.jl"),
    joinpath(VALIDATION_DIR, "relaxtime", "test_legacy_meson_solver_status_layered.jl"),
]

const VALIDATION_CORE_FILES = [
    joinpath(VALIDATION_DIR, "pnjl", "reference", "test_crossover_legacy_source_consistency.jl"),
    joinpath(VALIDATION_DIR, "relaxtime", "test_mott_reference_mapping.jl"),
    joinpath(VALIDATION_DIR, "relaxtime", "test_legacy_fortran_meson_trend_and_xi_label.jl"),
    joinpath(VALIDATION_DIR, "relaxtime", "test_meson_density_matches_mott_reference_continuation.jl"),
]

function _selected_validation_files()
    raw = strip(get(ENV, "VALIDATION_FILES", ""))
    isempty(raw) && return nothing

    parts = split(replace(raw, ',' => ';'), ';'; keepempty=false)
    files = String[]
    for part in parts
        entry = strip(part)
        isempty(entry) && continue
        path = isabspath(entry) ? entry : joinpath(VALIDATION_DIR, entry)
        isfile(path) || error("VALIDATION_FILES entry does not exist: $(path)")
        push!(files, path)
    end

    isempty(files) && return nothing
    return files
end

function _warn_local_non_smoke(profile::String)
    is_ci = get(ENV, "CI", "") in ("1", "true", "TRUE", "yes", "YES")
    is_ci && return
    profile == "smoke" && return
    @warn "Local validation run uses non-smoke profile; prefer smoke for edit-run loop" profile=profile recommended="VALIDATION_PROFILE=smoke"
end

function _include_validation_dir(dir::String)
    isdir(dir) || return
    entries = sort(readdir(dir; join=true))
    for entry in entries
        if isdir(entry)
            _include_validation_dir(entry)
            continue
        end

        file = lowercase(basename(entry))
        endswith(file, ".jl") || continue
        startswith(file, "test_") || continue
        rel = relpath(entry, VALIDATION_DIR)
        try
            include(entry)
        catch err
            println(stderr, "[validation] include failed: " * rel)
            rethrow(err)
        end
    end
end

@testset "Validation" begin
    selected = _selected_validation_files()
    if selected !== nothing
        @testset "Selected" begin
            for file in selected
                include(file)
            end
        end
        return
    end

    profile = lowercase(get(ENV, "VALIDATION_PROFILE", "smoke"))
    _warn_local_non_smoke(profile)

    if profile == "smoke"
        @testset "Smoke" begin
            for file in VALIDATION_SMOKE_FILES
                include(file)
            end
        end
    elseif profile == "core"
        @testset "Core" begin
            for file in VALIDATION_CORE_FILES
                include(file)
            end
        end
    elseif profile == "full"
        @testset "PNJL Baselines" begin
            _include_validation_dir(joinpath(VALIDATION_DIR, "pnjl"))
        end

        @testset "RelaxTime Baselines" begin
            _include_validation_dir(joinpath(VALIDATION_DIR, "relaxtime"))
        end
    else
        error("Unknown VALIDATION_PROFILE=$(profile). Use smoke, core, or full")
    end
end
