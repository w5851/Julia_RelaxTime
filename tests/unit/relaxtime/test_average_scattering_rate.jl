using Test
using BenchmarkTools

include("../../../src/relaxtime/AverageScatteringRate.jl")

using .AverageScatteringRate

# 构造简化参数，降低节点以加快单测
const QUARK_PARAMS = (
    m = (u = 1.52, d = 1.52, s = 2.50),
    μ = (u = 0.3, d = 0.3, s = 0.0)
)
const THERMO_ISO = (T = 0.15, Φ = 0.5, Φbar = 0.5, ξ = 0.0)
const THERMO_ANISO = (T = 0.15, Φ = 0.5, Φbar = 0.5, ξ = 0.2)
const K_COEFFS = (K_σπ=2.0, K_σK=2.0, K_σ=3.0, K_δπ=1.5, K_δK=1.5)

function constant_sigma_cache(process::Symbol; sigma::Float64=1.0)
    cache = CrossSectionCache(process)
    AverageScatteringRate.insert_sigma!(cache, 0.0, sigma)
    AverageScatteringRate.insert_sigma!(cache, 500.0, sigma)
    return cache
end

@testset "average_scattering_rate (isotropic)" begin
    cache = constant_sigma_cache(:uu_to_uu; sigma=1.0)
    ω = average_scattering_rate(
        :uu_to_uu,
        QUARK_PARAMS,
        THERMO_ISO,
        K_COEFFS;
        p_nodes=4,
        angle_nodes=2,
        phi_nodes=2,
        cs_cache=cache,
        n_sigma_points=4,
    )
    @test isfinite(ω)
    @test ω > 0
end

@testset "average_scattering_rate (anisotropic)" begin
    cache = constant_sigma_cache(:uu_to_uu; sigma=1.0)
    ω = average_scattering_rate(
        :uu_to_uu,
        QUARK_PARAMS,
        THERMO_ANISO,
        K_COEFFS;
        p_nodes=4,
        angle_nodes=2,
        phi_nodes=2,
        cs_cache=cache,
        n_sigma_points=4,
    )
    @test isfinite(ω)
    @test ω > 0
end

@testset "CrossSectionCache interpolation" begin
    cache = CrossSectionCache(:uu_to_uu)
    AverageScatteringRate.insert_sigma!(cache, 10.0, 1.0)
    AverageScatteringRate.insert_sigma!(cache, 20.0, 3.0)
    @test AverageScatteringRate.interpolate_sigma(cache, 10.0) == 1.0
    @test AverageScatteringRate.interpolate_sigma(cache, 20.0) == 3.0
    @test AverageScatteringRate.interpolate_sigma(cache, 15.0) ≈ 2.0 atol=1e-12
end

