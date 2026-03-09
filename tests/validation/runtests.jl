# Validation / baseline regression test entrypoint
#
# These tests compare numerical results against stored baselines
# to detect unintended numerical drift across versions.
#
# Run:
#   julia --project=. --eval 'include("tests/validation/runtests.jl")'

using Test

const VALIDATION_DIR = @__DIR__

include(joinpath(VALIDATION_DIR, "common", "data_paths.jl"))

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
    @testset "PNJL Baselines" begin
        _include_validation_dir(joinpath(VALIDATION_DIR, "pnjl"))
    end

    @testset "RelaxTime Baselines" begin
        _include_validation_dir(joinpath(VALIDATION_DIR, "relaxtime"))
    end
end
