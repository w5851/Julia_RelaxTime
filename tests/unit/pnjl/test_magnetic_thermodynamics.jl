# MagneticThermodynamics 模块单元测试
#
# 合并原 test_magnetic_density.jl, test_magnetic_coupling_GB.jl,
# test_magnetic_omega_components.jl

using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const PNJL = Models.pnjl_module()
const MagneticConfig = getproperty(PNJL, :MagneticConfig)
const coupling_GB = getproperty(PNJL, :coupling_GB)
const MAGNETIC_EB_MIN_FM2 = getproperty(PNJL, :MAGNETIC_EB_MIN_FM2)
const energy_landau = getproperty(PNJL, :energy_landau)
const calculate_magnetic_rho = getproperty(PNJL, :calculate_magnetic_rho)
const calculate_magnetic_number_densities = getproperty(PNJL, :calculate_magnetic_number_densities)
const calculate_magnetic_omega_components = getproperty(PNJL, :calculate_magnetic_omega_components)
const default_imc_params = getproperty(PNJL, :default_imc_params)
const default_magnetic_config = getproperty(PNJL, :default_magnetic_config)
const resolve_magnetic_nmax = getproperty(PNJL, :resolve_magnetic_nmax)
const omega_magnetic_mfir = getproperty(PNJL, :omega_magnetic_mfir)

@testset "MagneticThermodynamics" begin
    @testset "magnetic coupling G(B)" begin
        @test_throws ArgumentError coupling_GB(0.0)
        @test_throws ArgumentError coupling_GB(-MAGNETIC_EB_MIN_FM2)

        g1 = coupling_GB(MAGNETIC_EB_MIN_FM2)
        g2 = coupling_GB(0.10)
        @test isfinite(g1)
        @test isfinite(g2)
        @test g1 > 0
        @test g2 > 0
        @test default_imc_params().a ≈ 0.0108805
    end

    @testset "MFIR default and explicit legacy route" begin
        conf = default_magnetic_config(eB_fm2=MAGNETIC_EB_MIN_FM2)
        @test conf.route == :mfir
        @test conf.zeta_num >= 8
        @test isfinite(omega_magnetic_mfir(1.0, 2 / 3, MAGNETIC_EB_MIN_FM2; zeta_num=16))
        legacy = MagneticConfig(eB_fm2=MAGNETIC_EB_MIN_FM2, route=:landau_legacy)
        @test legacy.route == :landau_legacy
        @test_throws ArgumentError MagneticConfig(eB_fm2=MAGNETIC_EB_MIN_FM2, route=:smooth_landau)
    end

    @testset "positive magnetic-field contract" begin
        @test_throws ArgumentError MagneticConfig(eB_fm2=0.0)
        @test_throws ArgumentError MagneticConfig(eB_fm2=0.5 * MAGNETIC_EB_MIN_FM2)
        @test_throws ArgumentError energy_landau(1.0, 0.0, 0, 1 / 3, 0.0)
        conf = MagneticConfig(eB_fm2=MAGNETIC_EB_MIN_FM2)
        @test conf.eB_fm2 == MAGNETIC_EB_MIN_FM2
    end

    @testset "magnetic densities" begin
        x_state = SVector{5, Float64}(-0.03, -0.03, -0.04, 0.2, 0.2)
        mu_vec = SVector{3, Float64}(0.5, 0.5, 0.5)
        T_fm = 0.7
        conf = MagneticConfig(eB_fm2=0.08, p_num=24, pz_max=10.0)

        rho = calculate_magnetic_rho(x_state, mu_vec, T_fm, conf)
        @test all(isfinite.(rho))

        nd = calculate_magnetic_number_densities(x_state, mu_vec, T_fm, conf)
        @test all(isfinite.(nd.quark))
        @test nd.net === nd.quark
        @test nd.antiquark === nothing
        @test isfinite(nd.baryon)
        @test nd.quark ≈ rho rtol=5e-8 atol=1e-10

        rho_minus = calculate_magnetic_rho(x_state, -mu_vec, T_fm, conf)
        @test rho_minus ≈ -rho rtol=5e-8 atol=1e-10
    end

    @testset "magnetic control contract" begin
        x_state = SVector{5, Float64}(-0.03, -0.03, -0.04, 0.2, 0.2)
        mu_vec = SVector{3, Float64}(0.4, 0.4, 0.4)
        conf = MagneticConfig(eB_fm2=0.08, p_num=24, pz_max=10.0)
        @test_throws ArgumentError calculate_magnetic_omega_components(x_state, mu_vec, 0.0, conf)
        @test_throws ArgumentError calculate_magnetic_omega_components(x_state, mu_vec, 0.7, conf; xi=0.2)
        comp_low = calculate_magnetic_omega_components(x_state, mu_vec, 0.7, conf;
            p_num=4, t_num=4, pz_max=5.0, n_max=1)
        @test isfinite(comp_low.omega)
    end

    @testset "magnetic omega components" begin
        x_state = SVector{5, Float64}(-0.03, -0.03, -0.04, 0.2, 0.2)
        mu_vec = SVector{3, Float64}(0.4, 0.4, 0.4)
        T_fm = 0.8

        conf = MagneticConfig(eB_fm2=0.12, p_num=24, pz_max=10.0)
        comp = calculate_magnetic_omega_components(x_state, mu_vec, T_fm, conf)

        @test isfinite(comp.omega)
        @test isfinite(comp.vac)
        @test isfinite(comp.therm)
        @test comp.n_max >= 0
        @test comp.G_B > 0
        @test isapprox(comp.omega, comp.chi + comp.poly + comp.vac + comp.therm; rtol=1e-10, atol=1e-12)
    end

    @testset "shared thermal Landau cutoff profile" begin
        conf = default_magnetic_config(eB_fm2=0.3)
        μ = SVector{3, Float64}(0.4, 0.4, 0.4)
        n = resolve_magnetic_nmax(240.0 / 197.3269804, μ, conf)
        @test conf.n_max_policy == :thermal_tail
        @test conf.thermal_tail_factor == 30.0
        @test conf.n_max_floor == 3
        @test conf.n_max_cap == 10000
        @test n > 3
        @test n <= 10000
        @test resolve_magnetic_nmax(240.0 / 197.3269804, collect(μ), conf) == n

        source = default_magnetic_config(
            eB_fm2=0.2 * (1000.0 / 197.3269804)^2,
            profile="magnetic_source_parity",
        )
        @test source.n_max == 79
        @test source.p_num == 128
        @test source.pz_max == 40.0
        @test source.zeta_num == 256

        fixed = MagneticConfig(eB_fm2=0.1, n_max=79, n_max_policy=:thermal_tail)
        @test resolve_magnetic_nmax(1.0, μ, fixed) == 79
        over_budget = MagneticConfig(eB_fm2=MAGNETIC_EB_MIN_FM2, n_max_cap=10000)
        @test_throws ArgumentError resolve_magnetic_nmax(240.0 / 197.3269804, μ, over_budget)
        @test_throws ArgumentError MagneticConfig(eB_fm2=0.1, n_max_policy=:unknown)
        @test_throws ArgumentError MagneticConfig(eB_fm2=0.1, thermal_tail_factor=0.0)
        @test_throws ArgumentError MagneticConfig(eB_fm2=0.1, n_max_cap=2, n_max_floor=3)
    end
end
