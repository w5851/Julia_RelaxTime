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

function _rows_for(rows, meson::String, T_MeV::Float64, xi::Float64)
    out = Dict{Symbol,Any}[]
    for r in rows
        String(r[:source_impl]) == "fortran" || continue
        String(r[:meson]) == meson || continue
        isapprox(r[:T_MeV], T_MeV; atol=1e-12) || continue
        isapprox(r[:xi], xi; atol=1e-12) || continue
        push!(out, r)
    end
    return out
end

function _isotropic_rows(rows, meson::String)
    out = Dict{Float64,Dict{Symbol,Any}}()
    for r in rows
        String(r[:source_impl]) == "fortran" || continue
        String(r[:meson]) == meson || continue
        isapprox(r[:xi], 0.0; atol=1e-12) || continue
        T_MeV = Float64(r[:T_MeV])
        if !haskey(out, T_MeV)
            out[T_MeV] = r
        end
    end
    Ts = sort(collect(keys(out)))
    return Dict{Symbol,Any}[out[T] for T in Ts]
end

@testset "Fortran pi/K isotropic trend" begin
    rows = validate_reference_schema(load_reference_table(_FORTRAN_FILE))

    pi_rows = _isotropic_rows(rows, "pi")
    K_rows = _isotropic_rows(rows, "K")

    expected_Ts = [120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0]
    @test [r[:T_MeV] for r in pi_rows] == expected_Ts
    @test [r[:T_MeV] for r in K_rows] == expected_Ts

    pi_masses = Float64[r[:mass_MeV] for r in pi_rows]
    K_masses = Float64[r[:mass_MeV] for r in K_rows]

    @test all(isfinite, pi_masses)
    @test all(isfinite, K_masses)
    @test all(>(0.0), pi_masses)
    @test all(>(0.0), K_masses)

    pi_high = pi_masses[5:8]
    K_high = K_masses[5:8]
    @test all(diff(pi_high) .> 0.0)
    @test all(diff(K_high) .> 0.0)
end

@testset "Legacy Fortran xi is label-only" begin
    rows = validate_reference_schema(load_reference_table(_FORTRAN_FILE))

    xis = collect(-0.4:0.05:0.4)
    Ts = (160.0, 180.0, 200.0, 220.0)
    mesons = ("pi", "K")

    for meson in mesons
        for T_MeV in Ts
            baseline_row = _rows_for(rows, meson, T_MeV, 0.0)
            @test !isempty(baseline_row)
            baseline_masses = unique(Float64[r[:mass_MeV] for r in baseline_row])
            @test length(baseline_masses) == 1
            baseline = baseline_masses[1]

            for xi in xis
                matched = _rows_for(rows, meson, T_MeV, xi)
                @test !isempty(matched)
                masses = unique(Float64[r[:mass_MeV] for r in matched])
                @test length(masses) == 1
                @test isapprox(masses[1], baseline; rtol=0.0, atol=1e-9)
            end
        end
    end
end
