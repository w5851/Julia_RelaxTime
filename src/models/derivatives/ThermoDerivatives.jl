"""
    ThermoDerivatives

热力学导数计算模块，使用 ImplicitDifferentiation.jl 实现隐函数定理。

## 核心思想

对于非线性方程 F(x; θ) = 0，其中 x 是状态变量，θ = (T, μ) 是参数，
ImplicitDifferentiation.jl 自动使用隐函数定理计算 dx/dθ，支持任意阶导数。

策略约定：本模块导数统一采用 AD（ForwardDiff + ImplicitDifferentiation.jl）。
有限差分仅允许作为测试/诊断对照，不作为默认或主线路径。

## 主要函数

- `mass_derivatives`: 计算质量及其对 T/μ 的导数
- `thermo_derivatives`: 计算热力学量及其导数
- `bulk_viscosity_coefficients`: 计算体粘滞系数所需的导数

## 使用示例

```julia
# 一阶导数
md = mass_derivatives(T_fm, mu_fm)
println(md.dM_dT)  # ∂M/∂T

# 二阶导数
md2 = mass_derivatives(T_fm, mu_fm; order=2)
println(md2.d2M_dT2)  # ∂²M/∂T²
```
"""
module ThermoDerivatives

using ForwardDiff
using StaticArrays
using ImplicitDifferentiation
using LinearAlgebra: dot

# 导入默认积分规模常量（保持旧 API 习惯：p_num/t_num 默认值与 legacy 一致）。
using ..PNJLCore: DEFAULT_MOMENTUM_COUNT, DEFAULT_THETA_COUNT, cached_nodes

# Unified equilibrium facade (model-kind mapping, solve_gap helper)
const _EQUILIBRIUM_FACADE_PATH = normpath(joinpath(@__DIR__, "..", "pnjl_physics", "core", "EquilibriumFacade.jl"))
if !isdefined(Main, :EquilibriumFacade)
    Base.include(Main, _EQUILIBRIUM_FACADE_PATH)
end
const EquilibriumFacade = Main.EquilibriumFacade

# 导入常量
const _CONSTANTS_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "constants", "Constants_PNJL.jl"))
if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, _CONSTANTS_PATH)
end
using Main.Constants_PNJL: G_fm2, K_fm5

export mass_derivatives, thermo_derivatives, bulk_derivative_coeffs
export bulk_viscosity_coefficients, compute_B_bracket
export dP_dT, dP_dmu

# -----------------------------------------------------------------------------
# Compatibility helpers (used by tests / debugging)
# -----------------------------------------------------------------------------

"""get_thermal_nodes(p_num, t_num)

返回 legacy 热积分节点缓存（用于测试/诊断）。
"""
@inline get_thermal_nodes(p_num::Int, t_num::Int) = cached_nodes(p_num, t_num)

const _MODEL_CACHE = Dict{Symbol, Any}()

@inline function _get_model(model_kind::Symbol)
    return get!(_MODEL_CACHE, model_kind) do
        Main.Models.create_model(model_kind)
    end
end

@inline calculate_thermo(x_state, mu_vec, T_fm, thermal_nodes, xi) =
    Main.Models.model_thermo(
        _get_model(:PNJL),
        x_state,
        mu_vec,
        T_fm;
        p_num=CURRENT_P_NUM[],
        t_num=CURRENT_T_NUM[],
        xi=xi,
    )

@inline calculate_rho(x_state, mu_vec, T_fm, thermal_nodes, xi) =
    Main.Models.model_rho(
        _get_model(:PNJL),
        x_state,
        mu_vec,
        T_fm;
        p_num=CURRENT_P_NUM[],
        t_num=CURRENT_T_NUM[],
        xi=xi,
    )

"""IMPLICIT_SOLVER(θ)

返回当前 `set_config` 指定的后端/积分规模下的隐函数求解器输出 `(x_state, meta)`。
"""
@inline IMPLICIT_SOLVER(θ) = _implicit_gap_solver()(θ)

# ============================================================================
# Thermo backend selector
# ============================================================================

const _IMPLICIT_GAP_SOLVER_CACHE = Dict{Tuple{Float64, Int, Int, UInt}, Any}()

@inline function _implicit_gap_solver(; model=nothing)
    isdefined(Main, :Models) || error("Models not loaded; expected Main.Models")
    isdefined(Main.Models, :create_implicit_gap_solver) || error("Models.create_implicit_gap_solver is not defined")

    m = model === nothing ? _get_model(:PNJL) : model
    key = (Float64(CURRENT_XI[]), Int(CURRENT_P_NUM[]), Int(CURRENT_T_NUM[]), Base.objectid(m))
    if haskey(_IMPLICIT_GAP_SOLVER_CACHE, key)
        return _IMPLICIT_GAP_SOLVER_CACHE[key]
    end

    f = Main.Models.create_implicit_gap_solver(m;
        xi=CURRENT_XI[],
        p_num=CURRENT_P_NUM[],
        t_num=CURRENT_T_NUM[],
    )
    _IMPLICIT_GAP_SOLVER_CACHE[key] = f
    return f
end

@inline function _thermo_backend(x_state::SVector{5},
    mu_vec,
    T_fm,
    thermal_nodes,
    xi;
    p_num::Int,
    t_num::Int,
    model=nothing)

    m = model === nothing ? _get_model(:PNJL) : model
    return Main.Models.model_thermo(
        m,
        x_state,
        mu_vec,
        T_fm;
        p_num=p_num,
        t_num=t_num,
        xi=xi,
    )
end

@inline function _rho_backend(x_state::SVector{5},
    mu_vec,
    T_fm,
    thermal_nodes,
    xi;
    p_num::Int,
    t_num::Int,
    model=nothing)

    m = model === nothing ? _get_model(:PNJL) : model
    return Main.Models.model_rho(
        m,
        x_state,
        mu_vec,
        T_fm;
        p_num=p_num,
        t_num=t_num,
        xi=xi,
    )
end

# ============================================================================
# 全局配置
# ============================================================================

const CURRENT_XI = Ref{Float64}(0.0)
const CURRENT_P_NUM = Ref{Int}(DEFAULT_MOMENTUM_COUNT)
const CURRENT_T_NUM = Ref{Int}(DEFAULT_THETA_COUNT)

function set_config(; xi::Real=0.0, p_num::Int=DEFAULT_MOMENTUM_COUNT, t_num::Int=DEFAULT_THETA_COUNT)
    CURRENT_XI[] = Float64(xi)
    CURRENT_P_NUM[] = p_num
    CURRENT_T_NUM[] = t_num
end

# ============================================================================
# 质量计算
# ============================================================================

function compute_masses_from_state(x::AbstractVector)
    φ_u, φ_d, φ_s = x[1], x[2], x[3]
    m_u0 = 0.0055 / 0.197327
    m_s0 = 0.140 / 0.197327
    
    m_u = m_u0 - 4G_fm2 * φ_u + 2K_fm5 * φ_d * φ_s
    m_d = m_u0 - 4G_fm2 * φ_d + 2K_fm5 * φ_u * φ_s
    m_s = m_s0 - 4G_fm2 * φ_s + 2K_fm5 * φ_u * φ_d
    
    return SVector{3}(m_u, m_d, m_s)
end

# ============================================================================
# 质量导数
# ============================================================================

"""
    mass_derivatives(T_fm, mu_fm; order=1, xi=0.0, p_num, t_num)

计算夸克有效质量及其对 T/μ 的导数。
"""
function mass_derivatives(T_fm::Real, mu_fm::Real;
                          order::Int=1,
                          xi::Real=0.0,
                          p_num::Int=DEFAULT_MOMENTUM_COUNT,
                          t_num::Int=DEFAULT_THETA_COUNT,
                          model=nothing)
    set_config(xi=xi, p_num=p_num, t_num=t_num)
    θ = [T_fm, mu_fm]
    implicit = _implicit_gap_solver(model=model)
    
    function mass_vec_func(θ_in)
        (x_out, _) = implicit(θ_in)
        return collect(compute_masses_from_state(x_out))
    end
    
    masses_vec = mass_vec_func(θ)
    masses = SVector{3}(masses_vec...)
    
    if order == 1
        dM_dT = SVector{3}(ForwardDiff.derivative(T -> mass_vec_func([T, θ[2]])[i], θ[1]) for i in 1:3)
        dM_dmu = SVector{3}(ForwardDiff.derivative(μ -> mass_vec_func([θ[1], μ])[i], θ[2]) for i in 1:3)
        return (masses=masses, dM_dT=dM_dT, dM_dmu=dM_dmu)
        
    elseif order == 2
        dM_dT = SVector{3}(ForwardDiff.derivative(T -> mass_vec_func([T, θ[2]])[i], θ[1]) for i in 1:3)
        dM_dmu = SVector{3}(ForwardDiff.derivative(μ -> mass_vec_func([θ[1], μ])[i], θ[2]) for i in 1:3)
        
        d2M_dT2 = SVector{3}(
            ForwardDiff.derivative(T -> ForwardDiff.derivative(t -> mass_vec_func([t, θ[2]])[i], T), θ[1])
            for i in 1:3
        )
        d2M_dmu2 = SVector{3}(
            ForwardDiff.derivative(μ -> ForwardDiff.derivative(m -> mass_vec_func([θ[1], m])[i], μ), θ[2])
            for i in 1:3
        )
        d2M_dTdmu = SVector{3}(
            ForwardDiff.derivative(T -> ForwardDiff.derivative(μ -> mass_vec_func([T, μ])[i], θ[2]), θ[1])
            for i in 1:3
        )
        
        return (masses=masses, dM_dT=dM_dT, dM_dmu=dM_dmu,
                d2M_dT2=d2M_dT2, d2M_dTdmu=d2M_dTdmu, d2M_dmu2=d2M_dmu2)
    else
        error("order must be 1 or 2, got $order")
    end
end

# ============================================================================
# 热力学量导数
# ============================================================================

"""
    thermo_derivatives(T_fm, mu_fm; xi=0.0, p_num, t_num)

计算热力学量及其一阶导数。
"""
function thermo_derivatives(T_fm::Real, mu_fm::Real;
                            xi::Real=0.0,
                            p_num::Int=DEFAULT_MOMENTUM_COUNT,
                            t_num::Int=DEFAULT_THETA_COUNT,
                            model=nothing)
    set_config(xi=xi, p_num=p_num, t_num=t_num)
    θ = [T_fm, mu_fm]
    implicit = _implicit_gap_solver(model=model)
    
    function thermo_func(θ_in)
        (x_out, _) = implicit(θ_in)
        x_sv = SVector{5}(Tuple(x_out))
        mu_vec = SVector{3}(θ_in[2], θ_in[2], θ_in[2])

        P, rho_norm, s, ε = _thermo_backend(x_sv, mu_vec, θ_in[1], nothing, CURRENT_XI[];
            p_num=p_num, t_num=t_num, model=model)
        rho_vec = _rho_backend(x_sv, mu_vec, θ_in[1], nothing, CURRENT_XI[];
            p_num=p_num, t_num=t_num, model=model)
        n = sum(rho_vec) / 3
        
        return [P, ε, n, s]
    end
    
    base_vals = thermo_func(θ)
    pressure, energy, rho, entropy = base_vals[1], base_vals[2], base_vals[3], base_vals[4]

    # Derivative routines above treat x_state as implicit function of θ.
    # Here we also compute a consistent rho_norm for reporting.
    (x_base, _) = implicit(θ)
    x_sv = SVector{5}(Tuple(x_base))
    mu_vec = SVector{3}(θ[2], θ[2], θ[2])
    rho_vec0 = _rho_backend(x_sv, mu_vec, θ[1], nothing, CURRENT_XI[];
        p_num=p_num, t_num=t_num, model=model)
    rho_norm = sum(rho_vec0) / (3.0 * Main.Models.ρ0)
    
    P_T = ForwardDiff.derivative(T -> thermo_func([T, θ[2]])[1], θ[1])
    P_mu = ForwardDiff.derivative(μ -> thermo_func([θ[1], μ])[1], θ[2])
    E_T = ForwardDiff.derivative(T -> thermo_func([T, θ[2]])[2], θ[1])
    E_mu = ForwardDiff.derivative(μ -> thermo_func([θ[1], μ])[2], θ[2])
    n_T = ForwardDiff.derivative(T -> thermo_func([T, θ[2]])[3], θ[1])
    n_mu = ForwardDiff.derivative(μ -> thermo_func([θ[1], μ])[3], θ[2])
    
    denom_eps = E_T * n_mu - E_mu * n_T
    denom_n = n_T * E_mu - n_mu * E_T
    dP_depsilon_n = denom_eps == 0 ? NaN : (P_T * n_mu - P_mu * n_T) / denom_eps
    dP_dn_epsilon = denom_n == 0 ? NaN : (P_T * E_mu - P_mu * E_T) / denom_n
    
    md = mass_derivatives(T_fm, mu_fm; order=1, xi=xi, p_num=p_num, t_num=t_num, model=model)
    
    return (
        pressure=pressure,
        energy=energy,
        rho=rho,
        rho_norm=rho_norm,
        entropy=entropy,
        dP_dT=P_T,
        dP_dmu=P_mu,
        dEpsilon_dT=E_T,
        dEpsilon_dmu=E_mu,
        dn_dT=n_T,
        dn_dmu=n_mu,
        dP_depsilon_n=dP_depsilon_n,
        dP_dn_epsilon=dP_dn_epsilon,
        masses=md.masses,
        dM_dT=md.dM_dT,
        dM_dmu=md.dM_dmu,
        converged=true,
        iterations=missing,
        residual_norm=missing,
    )
end

"""
    bulk_derivative_coeffs(T_fm, mu_fm; kwargs...)

返回体粘滞系数公式中常用的导数组合。
"""
function bulk_derivative_coeffs(T_fm::Real, mu_fm::Real; kwargs...)
    derivs = thermo_derivatives(T_fm, mu_fm; kwargs...)
    return (
        dP_depsilon_n=derivs.dP_depsilon_n,
        dP_dn_epsilon=derivs.dP_dn_epsilon,
        dM_dT=derivs.dM_dT,
        dM_dmu=derivs.dM_dmu,
    )
end

# ============================================================================
# 高阶导数
# ============================================================================

function _nth_derivative(f, x, n::Int)
    n >= 1 || error("order must be ≥ 1")
    n == 1 && return ForwardDiff.derivative(f, x)
    inner = y -> _nth_derivative(f, y, n - 1)
    return ForwardDiff.derivative(inner, x)
end

"""
    dP_dT(T_fm, mu_fm; order=1, kwargs...)

总压强对温度的 n 阶导数。
"""
function dP_dT(T_fm::Real, mu_fm::Real; order::Int=1, xi::Real=0.0,
               p_num::Int=DEFAULT_MOMENTUM_COUNT, t_num::Int=DEFAULT_THETA_COUNT,
               model=nothing)
    set_config(xi=xi, p_num=p_num, t_num=t_num)
    implicit = _implicit_gap_solver(model=model)
    
    function P_func(T)
        θ = [T, mu_fm]
        (x_out, _) = implicit(θ)
        x_sv = SVector{5}(Tuple(x_out))
        mu_vec = SVector{3}(mu_fm, mu_fm, mu_fm)
        P, _, _, _ = _thermo_backend(x_sv, mu_vec, T, nothing, CURRENT_XI[];
            p_num=p_num, t_num=t_num, model=model)
        return P
    end
    
    return _nth_derivative(P_func, T_fm, order)
end

"""
    dP_dmu(T_fm, mu_fm; order=1, kwargs...)

总压强对化学势的 n 阶导数。
"""
function dP_dmu(T_fm::Real, mu_fm::Real; order::Int=1, xi::Real=0.0,
                p_num::Int=DEFAULT_MOMENTUM_COUNT, t_num::Int=DEFAULT_THETA_COUNT,
                model=nothing)
    set_config(xi=xi, p_num=p_num, t_num=t_num)
    implicit = _implicit_gap_solver(model=model)
    
    function P_func(μ)
        θ = [T_fm, μ]
        (x_out, _) = implicit(θ)
        x_sv = SVector{5}(Tuple(x_out))
        mu_vec = SVector{3}(μ, μ, μ)
        P, _, _, _ = _thermo_backend(x_sv, mu_vec, T_fm, nothing, CURRENT_XI[];
            p_num=p_num, t_num=t_num, model=model)
        return P
    end
    
    return _nth_derivative(P_func, mu_fm, order)
end

# ============================================================================
# 体粘滞系数
# ============================================================================

"""
    bulk_viscosity_coefficients(T_fm, mu_fm; kwargs...)

计算体粘滞系数公式所需的所有热力学导数系数。

使用 ForwardDiff + ImplicitDifferentiation.jl 实现自动微分。

## 技术说明

通过 ImplicitDifferentiation.jl 计算 dx/dθ，然后使用链式法则计算热力学量的导数。
热力学量统一通过 ModelThermodynamics（内部调用 Models.omega 并使用 ForwardDiff）计算，
因此会出现嵌套 AD（类似 legacy 的结构）。经测试在当前环境下可以正常工作。

链式法则：
- ds/dT = ∂s/∂T|_x + ∂s/∂x · dx/dT
- ds/dμ = ∂s/∂μ|_x + ∂s/∂x · dx/dμ
- dn/dT = ∂n/∂T|_x + ∂n/∂x · dx/dT
- dn/dμ = ∂n/∂μ|_x + ∂n/∂x · dx/dμ

其中 dx/dθ 通过 ImplicitDifferentiation.jl 计算。
"""
function bulk_viscosity_coefficients(T_fm::Real, mu_fm::Real;
                                     xi::Real=0.0,
                                     p_num::Int=DEFAULT_MOMENTUM_COUNT,
                                     t_num::Int=DEFAULT_THETA_COUNT,
                                     model=nothing)
    set_config(xi=xi, p_num=p_num, t_num=t_num)
    implicit = _implicit_gap_solver(model=model)
    
    T_val = T_fm
    μ_val = mu_fm
    θ = [T_val, μ_val]
    
    # 定义只计算状态变量的函数（通过 ImplicitFunction）
    function solve_state(θ_in)
        (x_out, _) = implicit(θ_in)
        return collect(x_out)
    end
    
    # 计算基础状态
    x_base = solve_state(θ)
    x_sv = SVector{5}(Tuple(x_base))
    mu_vec = SVector{3}(μ_val, μ_val, μ_val)
    
    # 计算基础热力学量（内部使用 ForwardDiff；legacy 与 models 都支持嵌套 AD）
    _, _, s, _ = _thermo_backend(x_sv, mu_vec, T_val, nothing, CURRENT_XI[];
        p_num=p_num, t_num=t_num, model=model)
    rho_vec = _rho_backend(x_sv, mu_vec, T_val, nothing, CURRENT_XI[];
        p_num=p_num, t_num=t_num, model=model)
    n_B = sum(rho_vec) / 3
    masses = compute_masses_from_state(x_base)
    
    # 使用 ForwardDiff.jacobian 计算 dx/dθ
    dx_dθ = ForwardDiff.jacobian(solve_state, θ)
    """
    # 中心差分替代方案（仅用于测试/诊断，默认不使用）
    _fd_step = sqrt(eps(Float64))
    ΔT = max(1e-5, _fd_step * abs(T_val))       # fm⁻¹，约 1e-5
    Δμ = max(1e-5, _fd_step * (abs(μ_val) + 1e-6))
    x_T_p = solve_state([T_val + ΔT, μ_val]);  x_T_m = solve_state([T_val - ΔT, μ_val])
    x_μ_p = solve_state([T_val, μ_val + Δμ]);  x_μ_m = solve_state([T_val, μ_val - Δμ])
    dx_dθ = hcat((x_T_p .- x_T_m) ./ (2ΔT), (x_μ_p .- x_μ_m) ./ (2Δμ))
    """
    # 计算 ∂s/∂x, ∂n/∂x（固定 T, μ）
    # 通过 ModelThermodynamics（嵌套 ForwardDiff 可以工作）
    function s_of_x(x_vec)
        x_s = SVector{5}(Tuple(x_vec))
        _, _, s_val, _ = _thermo_backend(x_s, mu_vec, T_val, nothing, CURRENT_XI[];
            p_num=p_num, t_num=t_num, model=model)
        return s_val
    end
    ds_dx = ForwardDiff.gradient(s_of_x, x_base)
    
    function n_of_x(x_vec)
        x_s = SVector{5}(Tuple(x_vec))
        rho = _rho_backend(x_s, mu_vec, T_val, nothing, CURRENT_XI[];
            p_num=p_num, t_num=t_num, model=model)
        return sum(rho) / 3
    end
    dn_dx = ForwardDiff.gradient(n_of_x, x_base)
    
    # 计算 ∂s/∂T, ∂s/∂μ, ∂n/∂T, ∂n/∂μ（固定 x）
    function s_of_T(T)
        _, _, s_val, _ = _thermo_backend(x_sv, mu_vec, T, nothing, CURRENT_XI[];
            p_num=p_num, t_num=t_num, model=model)
        return s_val
    end
    s_T_partial = ForwardDiff.derivative(s_of_T, T_val)
    
    function s_of_mu(μ)
        mu_v = SVector{3}(μ, μ, μ)
        _, _, s_val, _ = _thermo_backend(x_sv, mu_v, T_val, nothing, CURRENT_XI[];
            p_num=p_num, t_num=t_num, model=model)
        return s_val
    end
    s_μ_partial = ForwardDiff.derivative(s_of_mu, μ_val)
    
    function n_of_T(T)
        rho = _rho_backend(x_sv, mu_vec, T, nothing, CURRENT_XI[];
            p_num=p_num, t_num=t_num, model=model)
        return sum(rho) / 3
    end
    n_T_partial = ForwardDiff.derivative(n_of_T, T_val)
    
    function n_of_mu(μ)
        mu_v = SVector{3}(μ, μ, μ)
        rho = _rho_backend(x_sv, mu_v, T_val, nothing, CURRENT_XI[];
            p_num=p_num, t_num=t_num, model=model)
        return sum(rho) / 3
    end
    n_μ_partial = ForwardDiff.derivative(n_of_mu, μ_val)
    
    # 计算 ∂M/∂x（质量只依赖于 x，不直接依赖于 T, μ）
    function M_of_x(x_vec)
        m = compute_masses_from_state(x_vec)
        return collect(m)
    end
    dM_dx = ForwardDiff.jacobian(M_of_x, x_base)
    
    # 应用链式法则
    # ds/dT = ∂s/∂T + ∂s/∂x · dx/dT
    # ds/dμ = ∂s/∂μ + ∂s/∂x · dx/dμ
    ds_dT = s_T_partial + dot(ds_dx, dx_dθ[:, 1])
    ds_dμq = s_μ_partial + dot(ds_dx, dx_dθ[:, 2])
    
    dn_dT = n_T_partial + dot(dn_dx, dx_dθ[:, 1])
    dn_dμq = n_μ_partial + dot(dn_dx, dx_dθ[:, 2])
    
    # dM/dT = ∂M/∂x · dx/dT
    # dM/dμ = ∂M/∂x · dx/dμ
    dM_dT_vec = dM_dx * dx_dθ[:, 1]
    dM_dμq_vec = dM_dx * dx_dθ[:, 2]
    dM_dT = SVector{3}(dM_dT_vec...)
    dM_dμq = SVector{3}(dM_dμq_vec...)
    
    # 转换到重子化学势：∂/∂μ_B = (1/3)·∂/∂μ_q
    ds_dμB = ds_dμq / 3.0
    dn_dμB = dn_dμq / 3.0
    dM_dμB = dM_dμq ./ 3.0
    
    # v_n² = (s·∂n_B/∂μ_B - n_B·∂n_B/∂T) / [T·(∂s/∂T·∂n_B/∂μ_B - ∂s/∂μ_B·∂n_B/∂T)]
    numerator_vn = s * dn_dμB - n_B * dn_dT
    denominator_vn = T_val * (ds_dT * dn_dμB - ds_dμB * dn_dT)
    v_n_sq = numerator_vn / denominator_vn
    
    # ∂μ_B/∂T|_σ = - (∂σ/∂T) / (∂σ/∂μ_B), 其中 σ = s/n_B。
    # 为避免 n_B=0 时的除零问题，将分子分母同乘 n_B^2，得到等价形式：
    #   ∂μ_B/∂T|_σ = - (n_B·ds/dT - s·dn_B/dT) / (n_B·ds/dμ_B - s·dn_B/dμ_B)
    # 该形式在 n_B→0 的极限下通常是良定义的（例如 μ=0 线）。
    num_muB_T = n_B * ds_dT - s * dn_dT
    den_muB_T = n_B * ds_dμB - s * dn_dμB
    dμB_dT_sig = -num_muB_T / den_muB_T
    
    return (
        v_n_sq = v_n_sq,
        dμB_dT_sigma = dμB_dT_sig,
        masses = SVector{3}(masses...),
        dM_dT = dM_dT,
        dM_dμB = dM_dμB,
        s = s,
        n_B = n_B,
    )
end

"""
    compute_B_bracket(p, M, μq, T, v_n_sq, dμB_dT_sigma, dM_dT, dM_dμB; is_antiquark=false)

计算体粘滞公式中的 B 项。
"""
function compute_B_bracket(p::Real, M::Real, μq::Real, T::Real,
                           v_n_sq::Real, dμB_dT_sigma::Real,
                           dM_dT::Real, dM_dμB::Real;
                           is_antiquark::Bool=false)
    E = sqrt(p^2 + M^2)
    dE_dT_val = (M / E) * dM_dT
    dE_dμB = (M / E) * dM_dμB
    b_q = 1.0 / 3.0
    
    if is_antiquark
        dx_dT_sigma = dE_dT_val + (dE_dμB + b_q) * dμB_dT_sigma
        x = E + μq
    else
        dx_dT_sigma = dE_dT_val + (dE_dμB - b_q) * dμB_dT_sigma
        x = E - μq
    end
    
    dxt_dT_sigma = dx_dT_sigma / T - x / T^2
    B = p^2 + 3 * v_n_sq * T^2 * E * dxt_dT_sigma
    
    return B
end

end # module ThermoDerivatives

