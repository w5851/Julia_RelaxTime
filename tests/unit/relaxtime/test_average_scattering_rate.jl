using Test
using BenchmarkTools

# Load all relaxtime modules via RelaxTime.jl entry point
const _RELAX_TIME_MODULE_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "RelaxTime.jl"))
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RELAX_TIME_MODULE_PATH)
end

using Main.AverageScatteringRate
using Main.ParameterTypes: QuarkParams, ThermoParams
using Main.GaussLegendre: gauleg

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

@testset "average_scattering_rate struct/NamedTuple equivalence with cache" begin
    cache = constant_sigma_cache(:uu_to_uu; sigma=1.0)
    q_struct = QuarkParams(QUARK_PARAMS)
    t_struct = ThermoParams(THERMO_ISO)

    ω_nt = average_scattering_rate(
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

    ω_struct = average_scattering_rate(
        :uu_to_uu,
        q_struct,
        t_struct,
        K_COEFFS;
        p_nodes=4,
        angle_nodes=2,
        phi_nodes=2,
        cs_cache=cache,
        n_sigma_points=4,
    )

    @test isfinite(ω_nt)
    @test isfinite(ω_struct)
    @test isapprox(ω_struct, ω_nt; rtol=1e-12, atol=0.0)
end

@testset "CrossSectionCache interpolation" begin
    cache = CrossSectionCache(:uu_to_uu)
    AverageScatteringRate.insert_sigma!(cache, 10.0, 1.0)
    AverageScatteringRate.insert_sigma!(cache, 20.0, 3.0)
    @test AverageScatteringRate.interpolate_sigma(cache, 10.0) == 1.0
    @test AverageScatteringRate.interpolate_sigma(cache, 20.0) == 3.0
    @test AverageScatteringRate.interpolate_sigma(cache, 15.0) ≈ 2.0 atol=1e-12
end

@testset "number_density default grid matches explicit grids" begin
    p_nodes = 6
    angle_nodes = 4
    scale = 10.0

    t_grid, t_w = gauleg(0.0, 1.0, p_nodes)
    p_grid = Float64[]
    p_w = Float64[]
    for (t, wt) in zip(t_grid, t_w)
        if t >= 0.9999
            continue
        end
        inv_gap = 1.0 / (1.0 - t)
        push!(p_grid, scale * t * inv_gap)
        push!(p_w, wt * scale * inv_gap^2)
    end
    cos_grid, cos_w = gauleg(0.0, 1.0, angle_nodes)

    ρ_default = AverageScatteringRate.number_density(
        :u,
        QUARK_PARAMS.m.u,
        QUARK_PARAMS.μ.u,
        THERMO_ISO.T,
        THERMO_ISO.Φ,
        THERMO_ISO.Φbar,
        THERMO_ISO.ξ;
        p_nodes=p_nodes,
        angle_nodes=angle_nodes,
        scale=scale,
    )

    ρ_explicit = AverageScatteringRate.number_density(
        :u,
        QUARK_PARAMS.m.u,
        QUARK_PARAMS.μ.u,
        THERMO_ISO.T,
        THERMO_ISO.Φ,
        THERMO_ISO.Φbar,
        THERMO_ISO.ξ;
        p_grid=p_grid,
        p_w=p_w,
        cos_grid=cos_grid,
        cos_w=cos_w,
    )

    @test isfinite(ρ_default)
    @test isapprox(ρ_default, ρ_explicit; rtol=1e-12, atol=0.0)
end

@testset "average_scattering_rate default grids match explicit grids" begin
    cache = constant_sigma_cache(:uu_to_uu; sigma=1.0)
    p_nodes = 4
    angle_nodes = 2
    phi_nodes = 2
    scale = 10.0

    t_grid, t_w = gauleg(0.0, 1.0, p_nodes)
    p_grid = Float64[]
    p_w = Float64[]
    for (t, wt) in zip(t_grid, t_w)
        if t >= 0.9999
            continue
        end
        inv_gap = 1.0 / (1.0 - t)
        push!(p_grid, scale * t * inv_gap)
        push!(p_w, wt * scale * inv_gap^2)
    end
    cos_grid, cos_w = gauleg(-1.0, 1.0, angle_nodes)
    phi_grid, phi_w = gauleg(0.0, 2π, phi_nodes)

    ω_default = average_scattering_rate(
        :uu_to_uu,
        QUARK_PARAMS,
        THERMO_ISO,
        K_COEFFS;
        p_nodes=p_nodes,
        angle_nodes=angle_nodes,
        phi_nodes=phi_nodes,
        scale=scale,
        cs_cache=cache,
        n_sigma_points=4,
    )

    ω_explicit = average_scattering_rate(
        :uu_to_uu,
        QUARK_PARAMS,
        THERMO_ISO,
        K_COEFFS;
        p_grid=p_grid,
        p_w=p_w,
        cos_grid=cos_grid,
        cos_w=cos_w,
        phi_grid=phi_grid,
        phi_w=phi_w,
        cs_cache=cache,
        n_sigma_points=4,
    )

    @test isfinite(ω_default)
    @test isapprox(ω_default, ω_explicit; rtol=1e-12, atol=0.0)
end

