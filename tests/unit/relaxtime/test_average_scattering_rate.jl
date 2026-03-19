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
using Main.TotalCrossSection: total_cross_section
using Main.OneLoopIntegrals: A
using Main.EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using Main.Constants_PNJL: G_fm2, K_fm5

# 构造简化参数，降低节点以加快单测
const QUARK_PARAMS = (
    m = (u = 1.52, d = 1.52, s = 2.50),
    μ = (u = 0.3, d = 0.3, s = 0.0),
    A = (u = 0.2, d = 0.2, s = 0.2)
)
const THERMO_ISO = (T = 0.15, Φ = 0.5, Φbar = 0.5, ξ = 0.0)
const THERMO_ANISO = (T = 0.15, Φ = 0.5, Φbar = 0.5, ξ = 0.2)
const K_COEFFS = begin
    nodes_p, weights_p = gauleg(0.0, 20.0, 16)
    A_u = A(QUARK_PARAMS.m.u, QUARK_PARAMS.μ.u, THERMO_ISO.T, THERMO_ISO.Φ, THERMO_ISO.Φbar, nodes_p, weights_p)
    A_s = A(QUARK_PARAMS.m.s, QUARK_PARAMS.μ.s, THERMO_ISO.T, THERMO_ISO.Φ, THERMO_ISO.Φbar, nodes_p, weights_p)
    G_u = calculate_G_from_A(A_u, QUARK_PARAMS.m.u)
    G_s = calculate_G_from_A(A_s, QUARK_PARAMS.m.s)
    calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
end

function reference_sigma(process::Symbol, s::Float64)
    return total_cross_section(process, s, QUARK_PARAMS, THERMO_ISO, K_COEFFS; n_points=4)
end

function cache_sigma(cache::CrossSectionCache, s::Float64)
    return AverageScatteringRate.get_sigma(cache, s, QUARK_PARAMS, THERMO_ISO, K_COEFFS; n_points=4)
end

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

@testset "CrossSectionCache out-of-support contract" begin
    cache = CrossSectionCache(:uu_to_uu)
    AverageScatteringRate.insert_sigma!(cache, 10.0, 1.0)
    AverageScatteringRate.insert_sigma!(cache, 20.0, 3.0)

    @test AverageScatteringRate.interpolate_sigma(cache, 9.999) === nothing
    @test AverageScatteringRate.interpolate_sigma(cache, 20.001) === nothing

    @test cache_sigma(cache, 9.999) == 0.0
    @test cache_sigma(cache, 20.001) == 0.0
end

@testset "threshold injection invariants: trigger, sorted, unique" begin
    sparse_grid = [10.8, 12.5]
    cache_triggered = CrossSectionCache(:uu_to_uu)
    precompute_cross_section!(
        cache_triggered,
        sparse_grid,
        QUARK_PARAMS,
        THERMO_ISO,
        K_COEFFS;
        n_points=4,
        threshold_subtraction=true,
        asym_window=0.5,
        asym_fit_min_points=32,
        asym_extra_points=12,
    )

    @test length(cache_triggered.s_vals) > length(sparse_grid)
    @test cache_triggered.s_vals == sort(cache_triggered.s_vals)
    @test all(isfinite, cache_triggered.s_vals)
    @test all(isfinite, cache_triggered.sigma_vals)
    @test length(unique(cache_triggered.s_vals)) == length(cache_triggered.s_vals)

    dense_grid = collect(range(9.45, stop=9.85, length=8))
    cache_not_triggered = CrossSectionCache(:uu_to_uu)
    precompute_cross_section!(
        cache_not_triggered,
        dense_grid,
        QUARK_PARAMS,
        THERMO_ISO,
        K_COEFFS;
        n_points=4,
        threshold_subtraction=true,
        asym_window=0.5,
        asym_fit_min_points=1,
        asym_extra_points=12,
    )

    @test length(cache_not_triggered.s_vals) == length(dense_grid)
end

@testset "adaptive refinement improves under-resolved interpolation" begin
    process = :uubar_to_uubar
    s_th = (QUARK_PARAMS.m.u + QUARK_PARAMS.m.u)^2
    coarse_grid = [s_th + 0.01, s_th + 0.04, s_th + 0.15]

    coarse_cache = CrossSectionCache(process)
    precompute_cross_section!(
        coarse_cache,
        coarse_grid,
        QUARK_PARAMS,
        THERMO_ISO,
        K_COEFFS;
        n_points=4,
        threshold_subtraction=false,
    )

    refined_cache = CrossSectionCache(process)
    precompute_cross_section!(
        refined_cache,
        coarse_grid,
        QUARK_PARAMS,
        THERMO_ISO,
        K_COEFFS;
        n_points=4,
        threshold_subtraction=true,
        asym_window=0.2,
        asym_fit_min_points=16,
        asym_extra_points=20,
    )

    validation_s = [s_th + δ for δ in (1e-4, 5e-4, 1e-3, 2e-3, 5e-3, 8e-3, 0.012, 0.02, 0.03, 0.06, 0.10)]
    err_coarse = [abs(cache_sigma(coarse_cache, s) - reference_sigma(process, s)) for s in validation_s]
    err_refined = [abs(cache_sigma(refined_cache, s) - reference_sigma(process, s)) for s in validation_s]

    @test maximum(err_refined) < 0.9 * maximum(err_coarse)
    @test refined_cache.s_vals == sort(refined_cache.s_vals)
    @test all(isfinite, refined_cache.s_vals)
    @test all(isfinite, refined_cache.sigma_vals)
end

@testset "average_scattering_rate struct/NamedTuple equivalence keeps compatibility path" begin
    cache = constant_sigma_cache(:uu_to_uu; sigma=1.0)
    q_struct = QuarkParams(QUARK_PARAMS)
    t_iso_struct = ThermoParams(THERMO_ISO)
    t_aniso_struct = ThermoParams(THERMO_ANISO)

    ω_nt_iso = average_scattering_rate(
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
    ω_struct_iso = average_scattering_rate(
        :uu_to_uu,
        q_struct,
        t_iso_struct,
        K_COEFFS;
        p_nodes=4,
        angle_nodes=2,
        phi_nodes=2,
        cs_cache=cache,
        n_sigma_points=4,
    )

    ω_nt_aniso = average_scattering_rate(
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
    ω_struct_aniso = average_scattering_rate(
        :uu_to_uu,
        q_struct,
        t_aniso_struct,
        K_COEFFS;
        p_nodes=4,
        angle_nodes=2,
        phi_nodes=2,
        cs_cache=cache,
        n_sigma_points=4,
    )

    @test isapprox(ω_struct_iso, ω_nt_iso; rtol=1e-12, atol=0.0)
    @test isapprox(ω_struct_aniso, ω_nt_aniso; rtol=1e-12, atol=0.0)
end

@testset "design_w0cdf_s_grid struct/NamedTuple equivalence" begin
    q_struct = QuarkParams(QUARK_PARAMS)
    t_struct = ThermoParams(THERMO_ISO)

    s_grid_nt = AverageScatteringRate.design_w0cdf_s_grid(
        :uu_to_uu,
        QUARK_PARAMS,
        THERMO_ISO;
        N=12,
        p_nodes=4,
        angle_nodes=2,
        phi_nodes=2,
        p_cutoff=6.0,
    )

    s_grid_struct = AverageScatteringRate.design_w0cdf_s_grid(
        :uu_to_uu,
        q_struct,
        t_struct;
        N=12,
        p_nodes=4,
        angle_nodes=2,
        phi_nodes=2,
        p_cutoff=6.0,
    )

    @test length(s_grid_nt) == length(s_grid_struct)
    @test isapprox(first(s_grid_nt), first(s_grid_struct); rtol=1e-12, atol=0.0)
    @test isapprox(last(s_grid_nt), last(s_grid_struct); rtol=1e-12, atol=0.0)
    @test all(isfinite, s_grid_nt)
    @test all(isfinite, s_grid_struct)
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

@testset "average_scattering_rate default density grids match explicit density grids" begin
    cache = constant_sigma_cache(:uu_to_uu; sigma=1.0)
    p_nodes = 4
    angle_nodes = 2
    phi_nodes = 2
    density_p_nodes = 6
    scale = 10.0

    t_grid, t_w = gauleg(0.0, 1.0, density_p_nodes)
    density_p_grid = Float64[]
    density_p_w = Float64[]
    for (t, wt) in zip(t_grid, t_w)
        if t >= 0.9999
            continue
        end
        inv_gap = 1.0 / (1.0 - t)
        push!(density_p_grid, scale * t * inv_gap)
        push!(density_p_w, wt * scale * inv_gap^2)
    end

    ω_default_density = average_scattering_rate(
        :uu_to_uu,
        QUARK_PARAMS,
        THERMO_ISO,
        K_COEFFS;
        p_nodes=p_nodes,
        angle_nodes=angle_nodes,
        phi_nodes=phi_nodes,
        density_p_nodes=density_p_nodes,
        scale=scale,
        density_scale=scale,
        cs_cache=cache,
        n_sigma_points=4,
    )

    ω_explicit_density = average_scattering_rate(
        :uu_to_uu,
        QUARK_PARAMS,
        THERMO_ISO,
        K_COEFFS;
        p_nodes=p_nodes,
        angle_nodes=angle_nodes,
        phi_nodes=phi_nodes,
        density_p_grid=density_p_grid,
        density_p_w=density_p_w,
        scale=scale,
        density_scale=scale,
        cs_cache=cache,
        n_sigma_points=4,
    )

    @test isfinite(ω_default_density)
    @test isapprox(ω_default_density, ω_explicit_density; rtol=1e-12, atol=0.0)
end

@testset "average_scattering_rate direct mode bypasses cache interpolation" begin
    zero_cache = constant_sigma_cache(:uu_to_uu; sigma=0.0)

    ω_cached = average_scattering_rate(
        :uu_to_uu,
        QUARK_PARAMS,
        THERMO_ISO,
        K_COEFFS;
        p_nodes=4,
        angle_nodes=2,
        phi_nodes=2,
        cs_cache=zero_cache,
        n_sigma_points=4,
        interpolation_mode=:pchip,
    )

    ω_direct = average_scattering_rate(
        :uu_to_uu,
        QUARK_PARAMS,
        THERMO_ISO,
        K_COEFFS;
        p_nodes=4,
        angle_nodes=2,
        phi_nodes=2,
        cs_cache=zero_cache,
        n_sigma_points=4,
        interpolation_mode=:direct,
    )

    @test isfinite(ω_cached)
    @test isfinite(ω_direct)
    @test ω_cached == 0.0
    @test ω_direct > 0.0
end

@testset "average_scattering_rate hybrid_threshold mode stays finite" begin
    zero_cache = CrossSectionCache(:uu_to_uu)
    AverageScatteringRate.insert_sigma!(zero_cache, 0.0, 0.0)
    AverageScatteringRate.insert_sigma!(zero_cache, 1.0e6, 0.0)

    ω_hybrid = average_scattering_rate(
        :uu_to_uu,
        QUARK_PARAMS,
        THERMO_ISO,
        K_COEFFS;
        p_nodes=4,
        angle_nodes=2,
        phi_nodes=2,
        cs_cache=zero_cache,
        n_sigma_points=4,
        interpolation_mode=:hybrid_threshold,
        apply_s_domain_cut=false,
        sigma_cutoff=nothing,
    )

    @test isfinite(ω_hybrid)
    @test ω_hybrid >= 0.0
end

@testset "hybrid_threshold gate triggers only for configured channels" begin
    zero_cache = CrossSectionCache(:udbar_to_udbar)
    AverageScatteringRate.insert_sigma!(zero_cache, 0.0, 0.0)
    AverageScatteringRate.insert_sigma!(zero_cache, 1.0e6, 0.0)

    ω_hybrid = average_scattering_rate(
        :udbar_to_udbar,
        QUARK_PARAMS,
        THERMO_ISO,
        K_COEFFS;
        p_nodes=4,
        angle_nodes=2,
        phi_nodes=2,
        cs_cache=zero_cache,
        n_sigma_points=4,
        interpolation_mode=:hybrid_threshold,
        apply_s_domain_cut=false,
        sigma_cutoff=nothing,
    )

    @test isfinite(ω_hybrid)
    @test ω_hybrid == 0.0
end

@testset "hybrid_threshold mode remains deterministic for repeated evaluation" begin
    process = :uubar_to_uubar
    sparse_grid = [10.8, 12.5, 14.0, 18.0]

    cache = CrossSectionCache(process)
    precompute_cross_section!(
        cache,
        sparse_grid,
        QUARK_PARAMS,
        THERMO_ISO,
        K_COEFFS;
        n_points=4,
        threshold_subtraction=true,
        asym_window=0.8,
        asym_fit_min_points=8,
        asym_extra_points=16,
    )

    cache.peak_ratio = 100.0
    cache.peak_s = cache.asym_s0 + 0.1
    cache.peak_dirty = false

    row_s = cache.asym_s0 + 0.2
    v1 = AverageScatteringRate._get_sigma_core(
        cache,
        row_s,
        QUARK_PARAMS,
        THERMO_ISO,
        K_COEFFS;
        n_points=4,
        interpolation_mode=:hybrid_threshold,
    )

    v2 = AverageScatteringRate._get_sigma_core(
        cache,
        row_s,
        QUARK_PARAMS,
        THERMO_ISO,
        K_COEFFS;
        n_points=4,
        interpolation_mode=:hybrid_threshold,
    )

    @test isfinite(v1)
    @test isfinite(v2)
    @test isapprox(v2, v1; rtol=1e-12, atol=0.0)
end
