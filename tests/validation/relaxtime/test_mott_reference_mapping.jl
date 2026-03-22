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

using .MottReferenceMapping: load_reference_table, validate_reference_schema, estimate_mott_temperature
using .MottReferenceMapping: supported_theoretical_mesons, collect_mesons_from_rows

const _MOCK_FILE = validation_targets_path(
    "relaxtime", "legacy", "meson", "legacy_meson_scan_mock_muB0_v1.csv",
)

@testset "Mott reference mapping schema" begin
    rows = load_reference_table(_MOCK_FILE)
    @test !isempty(rows)

    validated = validate_reference_schema(rows)
    @test length(validated) == length(rows)
    @test all(r -> haskey(r, :meson), validated)
    @test all(r -> haskey(r, :gap_MeV), validated)
end

@testset "Mott temperature extraction" begin
    rows = load_reference_table(_MOCK_FILE)
    filtered = filter(r -> r[:source_impl] == "fortran" && r[:meson] == "pi", rows)
    @test length(filtered) == 2

    Ts = [r[:T_MeV] for r in filtered]
    gaps = [r[:gap_MeV] for r in filtered]

    result = estimate_mott_temperature(Ts, gaps)
    @test isapprox(result.T_mott_MeV, 156.6666666667; atol=1e-6)
    @test result.method == :linear
end

@testset "Mott fallback minabs" begin
    Ts = [150.0, 160.0]
    gaps = [-5.0, -1.0]
    result = estimate_mott_temperature(Ts, gaps)
    @test result.method == :minabs_approx
    @test result.approx
    @test result.T_mott_MeV == 160.0
end

@testset "Meson channel policy helpers" begin
    theoretical = supported_theoretical_mesons()
    @test :pi in theoretical
    @test :eta_prime in theoretical
    @test :sigma_prime in theoretical
    @test length(theoretical) == 8

    rows = load_reference_table(_MOCK_FILE)
    found = collect_mesons_from_rows(rows)
    @test found == [:K, :pi]
end
