# PNJL ThermoDerivatives 单元测试
#
# 测试内容：
# 1. thermo_derivatives/mass_derivatives 基本功能（fm⁻¹）
# 2. dP/dT 与中心差分的一致性（近似）
# 3. 返回类型检查（无 Dual 泄漏）

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

const _CONSTANTS_PNJL_PATH = normpath(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
const _GAUSS_LEGENDRE_PATH = normpath(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
const _MODELS_PATH = normpath(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, _CONSTANTS_PNJL_PATH)
end
if !isdefined(Main, :GaussLegendre)
    Base.include(Main, _GAUSS_LEGENDRE_PATH)
end
if !isdefined(Main, :Models)
    Base.include(Main, _MODELS_PATH)
end

using .Models
using Main.Constants_PNJL: ħc_MeV_fm
const PNJL = Models.pnjl_module()

# ============================================================================
# 基本功能（fm⁻¹）
# ============================================================================

@testset "thermo_derivatives basic finite" begin
    T_mev = 130.0
    mu_mev = 320.0
    T_fm = T_mev / ħc_MeV_fm
    mu_fm = mu_mev / ħc_MeV_fm

    derivs = PNJL.thermo_derivatives(T_fm, mu_fm; xi=0.0, p_num=16, t_num=6)
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
    T_fm = T_mev / ħc_MeV_fm
    mu_fm = mu_mev / ħc_MeV_fm

    δT_mev = 0.25
    δT_fm = δT_mev / ħc_MeV_fm

    P_plus = PNJL.thermo_derivatives(T_fm + δT_fm, mu_fm; xi=0.0, p_num=16, t_num=6).pressure
    P_minus = PNJL.thermo_derivatives(T_fm - δT_fm, mu_fm; xi=0.0, p_num=16, t_num=6).pressure
    fd_est = (P_plus - P_minus) / (2δT_fm)

    autodiff_val = PNJL.dP_dT(T_fm, mu_fm; xi=0.0, p_num=16, t_num=6)
    @test isapprox(autodiff_val, fd_est; atol=2e-2, rtol=3e-2)
end

# ============================================================================
# 新接口测试（fm⁻¹ 单位）
# ============================================================================

@testset "mass_derivatives (new interface)" begin
    T_fm = 0.5  # ~100 MeV
    μ_fm = 1.5  # ~300 MeV

    md = PNJL.mass_derivatives(T_fm, μ_fm; xi=0.0, p_num=16, t_num=6)
    
    @test all(isa.(md.masses, Float64))
    @test all(isa.(md.dM_dT, Float64))
    @test all(isa.(md.dM_dmu, Float64))
    @test all(isfinite.(md.masses))
    @test all(isfinite.(md.dM_dT))
    @test all(isfinite.(md.dM_dmu))
end

@testset "mass and bulk backend selector" begin
    T_fm = 0.5
    μ_fm = 1.5
    kwargs = (; xi=0.0, p_num=8, t_num=4)

    md_td = PNJL.mass_derivatives(T_fm, μ_fm; order=2, derivative_backend=:taylordiff, kwargs...)
    md_fd = PNJL.mass_derivatives(T_fm, μ_fm; order=2, derivative_backend=:forwarddiff, kwargs...)
    @test all(isapprox.(md_td.masses, md_fd.masses; rtol=1e-10, atol=1e-12))
    @test all(isapprox.(md_td.dM_dT, md_fd.dM_dT; rtol=1e-10, atol=1e-12))
    @test all(isapprox.(md_td.dM_dmu, md_fd.dM_dmu; rtol=1e-10, atol=1e-12))
    @test all(isapprox.(md_td.d2M_dTdmu, md_fd.d2M_dTdmu; rtol=1e-9, atol=1e-11))

    bv_td = PNJL.bulk_viscosity_coefficients(T_fm, μ_fm; derivative_backend=:taylordiff, kwargs...)
    bv_fd = PNJL.bulk_viscosity_coefficients(T_fm, μ_fm; derivative_backend=:forwarddiff, kwargs...)
    @test isapprox(bv_td.v_n_sq, bv_fd.v_n_sq; rtol=1e-10, atol=1e-12)
    @test isapprox(bv_td.dμB_dT_sigma, bv_fd.dμB_dT_sigma; rtol=1e-10, atol=1e-12)
    @test all(isapprox.(bv_td.masses, bv_fd.masses; rtol=1e-10, atol=1e-12))
end

@testset "thermo_derivatives (new interface)" begin
    T_fm = 0.5
    μ_fm = 1.5

    td = PNJL.thermo_derivatives(T_fm, μ_fm; xi=0.0, p_num=16, t_num=6)
    
    @test isa(td.pressure, Float64)
    @test isa(td.dP_dT, Float64)
    @test isa(td.dP_dmu, Float64)
    @test td.converged
    @test isfinite(td.pressure)
    @test isfinite(td.dP_dT)
    @test isfinite(td.dP_dmu)
end

@testset "thermo derivative backend selector" begin
    T_fm = 0.5
    μ_fm = 1.5
    kwargs = (; xi=0.0, p_num=8, t_num=4)

    auto = PNJL.thermo_derivatives(T_fm, μ_fm; derivative_backend=:auto, kwargs...)
    td = PNJL.thermo_derivatives(T_fm, μ_fm; derivative_backend=:taylordiff, kwargs...)
    fd = PNJL.thermo_derivatives(T_fm, μ_fm; derivative_backend=:forwarddiff, kwargs...)

    @test isapprox(auto.dP_dT, td.dP_dT; rtol=1e-12, atol=1e-12)
    @test isapprox(auto.dP_dmu, td.dP_dmu; rtol=1e-12, atol=1e-12)
    @test isapprox(td.dP_dT, fd.dP_dT; rtol=1e-10, atol=1e-12)
    @test isapprox(td.dP_dmu, fd.dP_dmu; rtol=1e-10, atol=1e-12)

    @test isapprox(PNJL.dP_dT(T_fm, μ_fm; derivative_backend=:auto, kwargs...),
        PNJL.dP_dT(T_fm, μ_fm; derivative_backend=:taylordiff, kwargs...);
        rtol=1e-12, atol=1e-12)
    @test isapprox(PNJL.dP_dmu(T_fm, μ_fm; derivative_backend=:taylordiff, kwargs...),
        PNJL.dP_dmu(T_fm, μ_fm; derivative_backend=:forwarddiff, kwargs...);
        rtol=1e-10, atol=1e-12)
end

@testset "bulk_viscosity_coefficients (new interface)" begin
    T_fm = 0.5
    μ_fm = 1.5

    bv = PNJL.bulk_viscosity_coefficients(T_fm, μ_fm; xi=0.0, p_num=16, t_num=6)
    
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

@testset "bulk_viscosity_coefficients supports mu=0" begin
    T_fm = 0.5
    μ_fm = 0.0
    bv = PNJL.bulk_viscosity_coefficients(T_fm, μ_fm; xi=0.0, p_num=12, t_num=4)
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
    local ThermoDerivatives = getproperty(PNJL, :ThermoDerivatives)
    local IMPLICIT_SOLVER = getproperty(ThermoDerivatives, :IMPLICIT_SOLVER)
    local get_thermal_nodes = getproperty(ThermoDerivatives, :get_thermal_nodes)
    local set_config = getproperty(ThermoDerivatives, :set_config)
    local CURRENT_XI = getproperty(ThermoDerivatives, :CURRENT_XI)
    local CURRENT_P_NUM = getproperty(ThermoDerivatives, :CURRENT_P_NUM)
    local CURRENT_T_NUM = getproperty(ThermoDerivatives, :CURRENT_T_NUM)
    local calculate_thermo = getproperty(ThermoDerivatives, :calculate_thermo)
    local calculate_rho = getproperty(ThermoDerivatives, :calculate_rho)
    local compute_masses_from_state = getproperty(ThermoDerivatives, :compute_masses_from_state)

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
