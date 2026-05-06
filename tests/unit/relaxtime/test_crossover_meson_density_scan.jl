using Test

const _MODELS_PATH_CMDS = joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    Base.include(Main, _MODELS_PATH_CMDS)
end

using Main.Models: run_crossover_meson_density_scan

@testset "crossover meson density scan argument validation" begin
    @test_throws ArgumentError run_crossover_meson_density_scan(
        mu_max_MeV=0.0,
        T_min_MeV=100.0,
        T_max_MeV=200.0,
    )
    @test_throws ArgumentError run_crossover_meson_density_scan(
        mu_max_MeV=200.0,
        mu_min_MeV=220.0,
        T_min_MeV=100.0,
        T_max_MeV=200.0,
    )
    @test_throws ArgumentError run_crossover_meson_density_scan(
        mu_max_MeV=200.0,
        T_min_MeV=200.0,
        T_max_MeV=100.0,
    )
    @test_throws ArgumentError run_crossover_meson_density_scan(
        mu_max_MeV=200.0,
        T_min_MeV=100.0,
        T_max_MeV=200.0,
        n_mu=2,
    )
end
