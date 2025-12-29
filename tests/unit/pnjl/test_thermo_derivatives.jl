# PNJL ThermoDerivatives 单元测试
#
# 测试内容：
# 1. 旧接口（MeV 单位）基本功能
# 2. 新接口（fm⁻¹ 单位）基本功能
# 3. 返回类型检查（无 Dual 泄漏）

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))

module PNJL
    using ..Constants_PNJL
    include(joinpath(Main.PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))
end

# ============================================================================
# 旧接口测试（MeV 单位）
# ============================================================================

@testset "thermo_derivatives basic finite (legacy)" begin
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

@testset "dP/dT matches central diff approximately (legacy)" begin
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

@testset "quasiparticle energy derivatives finite (legacy)" begin
    T_mev = 130.0
    mu_mev = 320.0
    p_fm = 0.5
    seed = PNJL.ThermoDerivatives.solve_equilibrium_mu(T_mev, mu_mev; xi=0.0, p_num=48, t_num=12)
    E = PNJL.ThermoDerivatives.quasiparticle_energy(T_mev, mu_mev, p_fm; xi=0.0, seed_state=seed.x_state, p_num=48, t_num=12)
    @test isfinite(E)
    @test isfinite(PNJL.ThermoDerivatives.dE_dT(T_mev, mu_mev, p_fm; xi=0.0, seed_state=seed.x_state, p_num=48, t_num=12))
    @test isfinite(PNJL.ThermoDerivatives.dE_dmu(T_mev, mu_mev, p_fm; xi=0.0, seed_state=seed.x_state, p_num=48, t_num=12))
end

# ============================================================================
# 新接口测试（fm⁻¹ 单位）
# ============================================================================

@testset "mass_derivatives (new interface)" begin
    T_fm = 0.5  # ~100 MeV
    μ_fm = 1.5  # ~300 MeV

    md = PNJL.ThermoDerivatives.mass_derivatives(T_fm, μ_fm)
    
    @test all(isa.(md.masses, Float64))
    @test all(isa.(md.dM_dT, Float64))
    @test all(isa.(md.dM_dmu, Float64))
    @test all(isfinite.(md.masses))
    @test all(isfinite.(md.dM_dT))
    @test all(isfinite.(md.dM_dmu))
end

@testset "thermo_derivatives (new interface)" begin
    T_fm = 0.5
    μ_fm = 1.5

    td = PNJL.ThermoDerivatives.thermo_derivatives(T_fm, μ_fm)
    
    @test isa(td.pressure, Float64)
    @test isa(td.dP_dT, Float64)
    @test isa(td.dP_dmu, Float64)
    @test td.converged
    @test isfinite(td.pressure)
    @test isfinite(td.dP_dT)
    @test isfinite(td.dP_dmu)
end

@testset "bulk_viscosity_coefficients (new interface)" begin
    T_fm = 0.5
    μ_fm = 1.5

    bv = PNJL.ThermoDerivatives.bulk_viscosity_coefficients(T_fm, μ_fm)
    
    # 类型检查
    @test typeof(bv.v_n_sq) == Float64
    @test typeof(bv.dμB_dT_sigma) == Float64
    @test eltype(bv.masses) == Float64
    @test typeof(bv.s) == Float64
    @test typeof(bv.n_B) == Float64
    
    # 数值有限
    @test isfinite(bv.v_n_sq)
    @test isfinite(bv.dμB_dT_sigma)
    @test all(isfinite.(bv.masses))
end

# ============================================================================
# 类型检查测试（确保无 Dual 泄漏）
# ============================================================================

using ForwardDiff
using StaticArrays

@testset "all_quantities type check" begin
    using .PNJL.ThermoDerivatives: IMPLICIT_SOLVER, get_thermal_nodes, set_config, 
        CURRENT_XI, CURRENT_P_NUM, CURRENT_T_NUM, calculate_thermo, calculate_rho,
        compute_masses_from_state

    set_config(xi=0.0)
    
    function all_quantities(θ::AbstractVector)
        (x_out, _) = IMPLICIT_SOLVER(θ)
        x_sv = SVector{5}(Tuple(x_out))
        mu_vec = SVector{3}(θ[2], θ[2], θ[2])
        
        nodes = get_thermal_nodes(CURRENT_P_NUM[], CURRENT_T_NUM[])
        _, _, s, _ = calculate_thermo(x_sv, mu_vec, θ[1], nodes, CURRENT_XI[])
        rho_vec = calculate_rho(x_sv, mu_vec, θ[1], nodes, CURRENT_XI[])
        n = sum(rho_vec) / 3
        masses = compute_masses_from_state(x_out)
        
        return [s, n, masses[1], masses[2], masses[3]]
    end

    θ = [0.5, 1.5]

    # 直接调用
    result = all_quantities(θ)
    @test eltype(result) == Float64

    # 计算 Jacobian
    J = ForwardDiff.jacobian(all_quantities, θ)
    @test size(J) == (5, 2)
    @test eltype(J) == Float64
    @test all(isfinite.(J))
end
