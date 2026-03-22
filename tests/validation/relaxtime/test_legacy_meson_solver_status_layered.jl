using Test

if !isdefined(Main, :validation_targets_path)
    include(joinpath(@__DIR__, "..", "common", "data_paths.jl"))
end

const _MOTT_MAP_SCRIPT = joinpath(
    @__DIR__,
    "..", "..", "..",
    "scripts", "analysis", "mott_reference_mapping.jl",
)

isfile(_MOTT_MAP_SCRIPT) || error("mott reference mapping script missing: $(_MOTT_MAP_SCRIPT)")

if !isdefined(Main, :MottReferenceMapping)
    Base.include(Main, _MOTT_MAP_SCRIPT)
end

using .MottReferenceMapping: load_reference_table, validate_reference_schema

const _FORTRAN_FILE = validation_targets_path(
    "relaxtime", "legacy", "meson", "legacy_meson_scan_fortran_muB0_v1.csv",
)

function _status_counts(rows)
    out = Dict{String,Int}()
    for r in rows
        s = String(r[:solver_status])
        out[s] = get(out, s, 0) + 1
    end
    return out
end

function _subset_count(rows, impl::String, mesons::Tuple{Vararg{String}}, statuses::Tuple{Vararg{String}})
    c = 0
    for r in rows
        String(r[:source_impl]) == impl || continue
        String(r[:meson]) in mesons || continue
        String(r[:solver_status]) in statuses || continue
        c += 1
    end
    return c
end

@testset "Legacy meson solver status layered assertions" begin
    rows_fortran = validate_reference_schema(load_reference_table(_FORTRAN_FILE))

    @test length(rows_fortran) == 608

    counts_fortran = _status_counts(rows_fortran)

    @test (
        get(counts_fortran, "proxy_mass_from_threshold", 0) == 456 ||
        get(counts_fortran, "exact_pole_legacy", 0) == 456
    )

    # Light channels (pi/K) should no longer be fully proxy-only.
    light_statuses = ("exact_light_pole", "exact_light_no_bracket", "exact_light_minabs")
    @test _subset_count(rows_fortran, "fortran", ("pi", "K"), light_statuses) > 0

    # Current staged baseline expects no hard failures in solver_status.
    @test get(counts_fortran, "fail", 0) == 0
end
