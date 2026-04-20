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

const PNJL_MODEL = Models.create_model(:PNJL)
const PNJL = Models.pnjl_module()
const MagneticConfig = getproperty(PNJL, :MagneticConfig)
const coupling_GB = getproperty(PNJL, :coupling_GB)
const PNJLConstants = getproperty(PNJL, :Constants_PNJL)
const G_fm2 = getproperty(PNJLConstants, :G_fm2)
const calculate_magnetic_rho = getproperty(PNJL, :calculate_magnetic_rho)
const calculate_magnetic_number_densities = getproperty(PNJL, :calculate_magnetic_number_densities)
const cached_nodes = getproperty(PNJL, :cached_nodes)
const calculate_magnetic_omega_components = getproperty(PNJL, :calculate_magnetic_omega_components)
const calculate_magnetic_omega = getproperty(PNJL, :calculate_magnetic_omega)

@testset "MagneticThermodynamics" begin
    @testset "magnetic coupling G(B)" begin
        g0 = coupling_GB(0.0)
        @test isapprox(g0, G_fm2; rtol=1e-12)

        g1 = coupling_GB(0.05)
        g2 = coupling_GB(0.10)
        @test isfinite(g1)
        @test isfinite(g2)
        @test g1 > 0
        @test g2 > 0
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
        @test isfinite(nd.baryon)
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

        thermal_nodes = cached_nodes(24, 8)
        omega_legacy = Models.omega(PNJL_MODEL, x_state, T_fm, mu_vec; thermal_nodes=thermal_nodes, xi=0.0)
        omega_b0 = calculate_magnetic_omega(x_state, mu_vec, T_fm, MagneticConfig(eB_fm2=0.0, p_num=24, pz_max=10.0))
        @test isapprox(omega_b0, omega_legacy; rtol=1e-7, atol=1e-8)

        comp_b0 = calculate_magnetic_omega_components(x_state, mu_vec, T_fm, MagneticConfig(eB_fm2=0.0, p_num=24, pz_max=10.0))
        masses_ref = Models.calculate_mass_vec(PNJL_MODEL, SVector{3, Float64}(x_state[1], x_state[2], x_state[3]))
        vac_ref = Models.vacuum_contribution(PNJL_MODEL, masses_ref)
        therm_ref = Models.thermal_contribution(PNJL_MODEL, masses_ref, x_state[4], x_state[5], mu_vec, T_fm; p_num=24, t_num=8, xi=0.0)

        @test isfinite(comp_b0.vac)
        @test isfinite(comp_b0.therm)
        @test isapprox(comp_b0.vac, vac_ref; rtol=1e-8, atol=1e-10)
        @test isapprox(comp_b0.therm, therm_ref; rtol=1e-8, atol=1e-10)
        @test isapprox(comp_b0.omega, comp_b0.chi + comp_b0.poly + comp_b0.vac + comp_b0.therm; rtol=1e-10, atol=1e-12)
    end
end
