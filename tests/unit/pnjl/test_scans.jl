using Test
using Printf
using DelimitedFiles

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

using .PNJL
using .PNJL.TmuScan
using .PNJL.TrhoScan
using .PNJL.AdaptiveRhoRefinement
using .PNJL.SeedCache: DEFAULT_SEED_PATH

@testset "TrhoScan defaults" begin
    @test first(DEFAULT_RHO_VALUES) == 0.0
    @test last(DEFAULT_RHO_VALUES) == 3.0
    @test issorted(DEFAULT_RHO_VALUES)
    low_count = count(x -> x <= 0.3, DEFAULT_RHO_VALUES)
    @test low_count > 40
end

@testset "Adaptive rho refinement" begin
    rho = [0.0, 0.1, 0.2, 0.4]
    mu = [320.0, 318.0, 317.9, 340.0]
    config = AdaptiveRhoConfig(; slope_tol=2.0, min_gap=0.05, digits=4, max_points=4)
    extra = suggest_refinement_points(rho, mu; config=config)
    @test length(extra) == 1
    @test abs(extra[1] - 0.15) < 1e-4
    merged = merge_rho_values(rho, extra; digits=4)
    @test issorted(merged)
    @test any(x -> abs(x - 0.15) < 1e-4, merged)
    strict_config = AdaptiveRhoConfig(; slope_tol=0.5, min_gap=0.05)
    @test isempty(suggest_refinement_points(rho, mu; config=strict_config))
end

function _read_data(path)
    lines = readlines(path)
    @test length(lines) >= 2
    header = split(lines[1], ',')
    data = split(lines[2], ',')
    return header, data
end

@testset "TmuScan single point" begin
    tmp_dir = mktempdir()
    output = joinpath(tmp_dir, "tmu_scan.csv")
    stats = run_tmu_scan(
        T_values = [90.0],
        mu_values = [10.0],
        xi_values = [0.0],
        output_path = output,
        seed_path = DEFAULT_SEED_PATH,
        overwrite = true,
        resume = false,
        p_num = 12,
        t_num = 6,
    )
    @test isfile(output)
    @test stats.total == 1
    @test stats.success >= 0
    header, data = _read_data(output)
    @test header[1:3] == ["T_MeV", "mu_MeV", "xi"]
    @test length(header) == length(data)
    @test parse(Float64, data[1]) ≈ 90.0
    @test parse(Float64, data[2]) ≈ 10.0
end

@testset "TrhoScan single point" begin
    tmp_dir = mktempdir()
    output = joinpath(tmp_dir, "trho_scan.csv")
    stats = run_trho_scan(
        T_values = [90.0],
        rho_values = [0.2],
        xi_values = [0.0],
        output_path = output,
        seed_path = DEFAULT_SEED_PATH,
        overwrite = true,
        resume = false,
        p_num = 12,
        t_num = 6,
    )
    @test isfile(output)
    @test stats.total == 1
    header, data = _read_data(output)
    @test header[1:3] == ["T_MeV", "rho", "xi"]
    @test length(header) == length(data)
    @test parse(Float64, data[1]) ≈ 90.0
    @test parse(Float64, data[2]) ≈ 0.2
    @test data[18] == "true"
    mu_u = parse(Float64, data[4])
    @test isfinite(mu_u)
end
