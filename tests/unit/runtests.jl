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
#   UNIT_PROFILE=smoke|full
#     - smoke (default): curated, should be green
#     - full: include most test_*.jl (still may fail until legacy tests are updated)

using Test

const UNIT_DIR = @__DIR__

# A small blacklist for tests that are currently WIP / outdated / intentionally non-unit.
# Default `include("tests/unit/runtests.jl")` should be green and reasonably fast.
const DEFAULT_SKIP = Set([
    # PNJL: missing symbol / AD singularities / references removed APIs
    "test_core_integrals.jl",
    "test_implicit_jacobian.jl",
    "test_scans.jl",

    # RelaxTime: uses deprecated API signature (should be updated before re-enabling)
    "test_differential_cross_section.jl",
])

const SMOKE_FILES = [
    # PNJL solver physical robustness (includes deterministic random sampling)
    joinpath(UNIT_DIR, "pnjl", "test_solver_random_physical_smoke.jl"),

    # Core numerics / integration utils
    joinpath(UNIT_DIR, "integration", "test_gausslegendre.jl"),
    joinpath(UNIT_DIR, "integration", "test_cauchypv.jl"),

    # RelaxTime core numerics
    joinpath(UNIT_DIR, "relaxtime", "test_b0_correction.jl"),
]

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
    for f in files
        _should_include_unit_file(f) || continue
        include(f)
    end
end

@testset "Unit" begin
    profile = lowercase(get(ENV, "UNIT_PROFILE", "smoke"))

    if profile == "smoke"
        @testset "Smoke" begin
            for f in SMOKE_FILES
                include(f)
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
