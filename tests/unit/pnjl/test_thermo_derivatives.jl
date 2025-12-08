# PNJL ThermoDerivatives 单元测试

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))

module PNJL
    using ..Constants_PNJL
    include(joinpath(Main.PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))
end

@testset "thermo_derivatives basic finite" begin
    T_mev = 130.0
    mu_mev = 320.0
    derivs = PNJL.ThermoDerivatives.thermo_derivatives(T_mev, mu_mev; xi=0.0, p_num=48, t_num=12)
    @test derivs.converged
    @test isfinite(derivs.pressure)
    @test isfinite(derivs.energy)
    @test isfinite(derivs.dP_dT)
    @test isfinite(derivs.dP_dmu)
    @test isfinite(derivs.dEpsilon_dT)
    @test isfinite(derivs.dEpsilon_dmu)
    @test isfinite(derivs.dn_dT)
    @test isfinite(derivs.dn_dmu)
    @test !isnan(derivs.dP_depsilon_n)
    @test !isnan(derivs.dP_dn_epsilon)
end

@testset "dP/dT matches central diff approximately" begin
    T_mev = 130.0
    mu_mev = 320.0
    seed = PNJL.ThermoDerivatives.solve_equilibrium_mu(T_mev, mu_mev; xi=0.0, p_num=48, t_num=12)
    δT = 0.25
    f_plus = PNJL.ThermoDerivatives.solve_equilibrium_mu(T_mev + δT, mu_mev; xi=0.0, seed_state=seed.x_state, p_num=48, t_num=12).pressure
    f_minus = PNJL.ThermoDerivatives.solve_equilibrium_mu(T_mev - δT, mu_mev; xi=0.0, seed_state=seed.x_state, p_num=48, t_num=12).pressure
    fd_est = (f_plus - f_minus) / (2δT)
    autodiff_val = PNJL.ThermoDerivatives.dP_dT(T_mev, mu_mev; xi=0.0, seed_state=seed.x_state, p_num=48, t_num=12)
    @test isapprox(autodiff_val, fd_est; atol=5e-3, rtol=1e-2)
end

@testset "quasiparticle energy derivatives finite" begin
    T_mev = 130.0
    mu_mev = 320.0
    p_fm = 0.5
    seed = PNJL.ThermoDerivatives.solve_equilibrium_mu(T_mev, mu_mev; xi=0.0, p_num=48, t_num=12)
    E = PNJL.ThermoDerivatives.quasiparticle_energy(T_mev, mu_mev, p_fm; xi=0.0, seed_state=seed.x_state, p_num=48, t_num=12)
    @test isfinite(E)
    @test isfinite(PNJL.ThermoDerivatives.dE_dT(T_mev, mu_mev, p_fm; xi=0.0, seed_state=seed.x_state, p_num=48, t_num=12))
    @test isfinite(PNJL.ThermoDerivatives.dE_dmu(T_mev, mu_mev, p_fm; xi=0.0, seed_state=seed.x_state, p_num=48, t_num=12))
end
