"""
    ThermoDerivatives

热力学导数计算模块，默认使用 TaylorDiff 显式 Taylor-series gap Newton 路线。

## 核心思想

对于非线性方程 F(x; θ) = 0，其中 x 是状态变量，θ = (T, μ) 是参数，
路线先求 primal gap 解，再在单变量 Taylor 代数中显式求解
`F(x(δ), T(δ), μ(δ)) = 0`。`derivative_backend=:auto` 与
`:taylordiff` 都使用这条 TD-only 生产路径；旧
`derivative_backend=:forwarddiff` implicit fallback 已下线。

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

using StaticArrays
import ..Models
using ..PNJLChiBTaylorDiff: gap_series_parameter_direction, pressure_series_parameter_direction
using ..PNJLChiBTaylorDiff: nth_derivative_from_series

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

const _DERIVATIVE_BACKENDS = (:auto, :taylordiff)

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
        Models.create_model(model_kind)
    end
end

@inline calculate_thermo(x_state, mu_vec, T_fm, thermal_nodes, xi) =
    Models.model_thermo(
        _get_model(:PNJL),
        x_state,
        mu_vec,
        T_fm;
        p_num=CURRENT_P_NUM[],
        t_num=CURRENT_T_NUM[],
        xi=xi,
    )

@inline calculate_rho(x_state, mu_vec, T_fm, thermal_nodes, xi) =
    Models.model_rho(
        _get_model(:PNJL),
        x_state,
        mu_vec,
        T_fm;
        p_num=CURRENT_P_NUM[],
        t_num=CURRENT_T_NUM[],
        xi=xi,
    )

# ============================================================================
# Thermo backend selector
# ============================================================================

@inline function _thermo_backend(x_state::SVector{5},
    mu_vec,
    T_fm,
    thermal_nodes,
    xi;
    p_num::Int,
    t_num::Int,
    model=nothing)

    m = model === nothing ? _get_model(:PNJL) : model
    return Models.model_thermo(
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
    return Models.model_rho(
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
# TaylorDiff backend helpers
# ============================================================================

@inline function _validate_derivative_backend(derivative_backend::Symbol)
    if derivative_backend === :forwarddiff
        throw(ArgumentError("derivative_backend=:forwarddiff has been retired from ThermoDerivatives; use derivative_backend=:auto or :taylordiff. For residual adapter audits, use build_*_problem(...).forward_solve/conditions instead of retired implicit-solver factories."))
    end
    derivative_backend in _DERIVATIVE_BACKENDS && return derivative_backend
    throw(ArgumentError("derivative_backend must be one of $(_DERIVATIVE_BACKENDS), got $(derivative_backend)"))
end

@inline function _resolve_thermo_backend(model, derivative_backend::Symbol)
    m = model === nothing ? _get_model(:PNJL) : model
    backend = _validate_derivative_backend(derivative_backend)
    backend = backend === :auto ? :taylordiff : backend
    if backend === :taylordiff && !(m isa Models.AbstractPNJLModel)
        throw(ArgumentError("derivative_backend=:taylordiff requires an AbstractPNJLModel, got $(typeof(m))"))
    end
    return m, backend
end

@inline _symmetric_mu_vec(mu_fm::Real) =
    SVector{3, Float64}(Float64(mu_fm), Float64(mu_fm), Float64(mu_fm))

@inline _zero_mu_direction() = SVector{3, Float64}(0.0, 0.0, 0.0)
@inline _symmetric_mu_direction() = SVector{3, Float64}(1.0, 1.0, 1.0)

@inline function _series_kwargs(; xi, p_num, t_num, model, series_iterations, linear_solve, series_residual_tol)
    return (;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        model=model,
        series_iterations=series_iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )
end

@inline function _extract_series_vector(v::SVector{N}, order::Int) where {N}
    return SVector{N, Float64}(ntuple(i -> nth_derivative_from_series(v[i], order), Val(N)))
end

function _gap_state_series_taylordiff(
    model::Models.AbstractPNJLModel,
    T_fm::Real,
    mu_fm::Real,
    T_direction::Real,
    mu_direction;
    order::Int,
    xi::Real,
    p_num::Int,
    t_num::Int,
    series_iterations::Union{Nothing, Int},
    linear_solve::Symbol,
    series_residual_tol::Real,
)
    return gap_series_parameter_direction(
        T_fm,
        _symmetric_mu_vec(mu_fm),
        T_direction,
        mu_direction;
        order=order,
        _series_kwargs(
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            model=model,
            series_iterations=series_iterations,
            linear_solve=linear_solve,
            series_residual_tol=series_residual_tol,
        )...,
    )
end

function _pressure_series_taylordiff(
    model::Models.AbstractPNJLModel,
    T_fm::Real,
    mu_fm::Real,
    T_direction::Real,
    mu_direction;
    order::Int,
    xi::Real,
    p_num::Int,
    t_num::Int,
    series_iterations::Union{Nothing, Int},
    linear_solve::Symbol,
    series_residual_tol::Real,
)
    return pressure_series_parameter_direction(
        T_fm,
        _symmetric_mu_vec(mu_fm),
        T_direction,
        mu_direction;
        order=order,
        _series_kwargs(
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            model=model,
            series_iterations=series_iterations,
            linear_solve=linear_solve,
            series_residual_tol=series_residual_tol,
        )...,
    )
end

function _mass_series_taylordiff(
    model::Models.AbstractPNJLModel,
    T_fm::Real,
    mu_fm::Real,
    T_direction::Real,
    mu_direction;
    order::Int,
    xi::Real,
    p_num::Int,
    t_num::Int,
    series_iterations::Union{Nothing, Int},
    linear_solve::Symbol,
    series_residual_tol::Real,
)
    result = _gap_state_series_taylordiff(
        model,
        T_fm,
        mu_fm,
        T_direction,
        mu_direction;
        order=order,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        series_iterations=series_iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )
    return compute_masses_from_state(result.x_state)
end

function _pressure_derivatives_order2_taylordiff(
    model::Models.AbstractPNJLModel,
    T_fm::Real,
    mu_fm::Real;
    xi::Real,
    p_num::Int,
    t_num::Int,
    series_iterations::Union{Nothing, Int},
    linear_solve::Symbol,
    series_residual_tol::Real,
)
    zero_mu = _zero_mu_direction()
    sym_mu = _symmetric_mu_direction()

    p_T = _pressure_series_taylordiff(
        model, T_fm, mu_fm, 1.0, zero_mu;
        order=2, xi=xi, p_num=p_num, t_num=t_num,
        series_iterations=series_iterations, linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )
    p_mu = _pressure_series_taylordiff(
        model, T_fm, mu_fm, 0.0, sym_mu;
        order=2, xi=xi, p_num=p_num, t_num=t_num,
        series_iterations=series_iterations, linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )
    p_Tmu = _pressure_series_taylordiff(
        model, T_fm, mu_fm, 1.0, sym_mu;
        order=2, xi=xi, p_num=p_num, t_num=t_num,
        series_iterations=series_iterations, linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )

    pressure = nth_derivative_from_series(p_T, 0)
    P_T = nth_derivative_from_series(p_T, 1)
    P_mu = nth_derivative_from_series(p_mu, 1)
    P_TT = nth_derivative_from_series(p_T, 2)
    P_mumu = nth_derivative_from_series(p_mu, 2)
    P_Tmu = (nth_derivative_from_series(p_Tmu, 2) - P_TT - P_mumu) / 2

    return (
        pressure=pressure,
        P_T=P_T,
        P_mu=P_mu,
        P_TT=P_TT,
        P_Tmu=P_Tmu,
        P_mumu=P_mumu,
    )
end

function _mass_derivatives_taylordiff(
    T_fm::Real,
    mu_fm::Real;
    order::Int,
    xi::Real,
    p_num::Int,
    t_num::Int,
    model::Models.AbstractPNJLModel,
    series_iterations::Union{Nothing, Int},
    linear_solve::Symbol,
    series_residual_tol::Real,
)
    order == 1 || order == 2 || throw(ArgumentError("order must be 1 or 2, got $order"))
    zero_mu = _zero_mu_direction()
    sym_mu = _symmetric_mu_direction()

    m_T = _mass_series_taylordiff(
        model, T_fm, mu_fm, 1.0, zero_mu;
        order=order, xi=xi, p_num=p_num, t_num=t_num,
        series_iterations=series_iterations, linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )
    m_mu = _mass_series_taylordiff(
        model, T_fm, mu_fm, 0.0, sym_mu;
        order=order, xi=xi, p_num=p_num, t_num=t_num,
        series_iterations=series_iterations, linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )

    masses = _extract_series_vector(m_T, 0)
    dM_dT = _extract_series_vector(m_T, 1)
    dM_dmu = _extract_series_vector(m_mu, 1)

    if order == 1
        return (masses=masses, dM_dT=dM_dT, dM_dmu=dM_dmu)
    end

    m_Tmu = _mass_series_taylordiff(
        model, T_fm, mu_fm, 1.0, sym_mu;
        order=2, xi=xi, p_num=p_num, t_num=t_num,
        series_iterations=series_iterations, linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )

    d2M_dT2 = _extract_series_vector(m_T, 2)
    d2M_dmu2 = _extract_series_vector(m_mu, 2)
    d2M_dTdmu = (_extract_series_vector(m_Tmu, 2) - d2M_dT2 - d2M_dmu2) / 2

    return (
        masses=masses,
        dM_dT=dM_dT,
        dM_dmu=dM_dmu,
        d2M_dT2=d2M_dT2,
        d2M_dTdmu=d2M_dTdmu,
        d2M_dmu2=d2M_dmu2,
    )
end

function _thermo_derivatives_taylordiff(
    T_fm::Real,
    mu_fm::Real;
    xi::Real,
    p_num::Int,
    t_num::Int,
    model::Models.AbstractPNJLModel,
    series_iterations::Union{Nothing, Int},
    linear_solve::Symbol,
    series_residual_tol::Real,
)
    pd = _pressure_derivatives_order2_taylordiff(
        model,
        T_fm,
        mu_fm;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        series_iterations=series_iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )

    T_val = Float64(T_fm)
    mu_val = Float64(mu_fm)
    rho = pd.P_mu / 3
    entropy = pd.P_T
    energy = -pd.pressure + mu_val * pd.P_mu + T_val * pd.P_T
    E_T = mu_val * pd.P_Tmu + T_val * pd.P_TT
    E_mu = mu_val * pd.P_mumu + T_val * pd.P_Tmu
    n_T = pd.P_Tmu / 3
    n_mu = pd.P_mumu / 3

    denom_eps = E_T * n_mu - E_mu * n_T
    denom_n = n_T * E_mu - n_mu * E_T
    dP_depsilon_n = denom_eps == 0 ? NaN : (pd.P_T * n_mu - pd.P_mu * n_T) / denom_eps
    dP_dn_epsilon = denom_n == 0 ? NaN : (pd.P_T * E_mu - pd.P_mu * E_T) / denom_n

    md = _mass_derivatives_taylordiff(
        T_fm,
        mu_fm;
        order=1,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        model=model,
        series_iterations=series_iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )

    return (
        pressure=pd.pressure,
        energy=energy,
        rho=rho,
        rho_norm=rho / Models.ρ0,
        entropy=entropy,
        dP_dT=pd.P_T,
        dP_dmu=pd.P_mu,
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

function _bulk_viscosity_coefficients_taylordiff(
    T_fm::Real,
    mu_fm::Real;
    xi::Real,
    p_num::Int,
    t_num::Int,
    model::Models.AbstractPNJLModel,
    series_iterations::Union{Nothing, Int},
    linear_solve::Symbol,
    series_residual_tol::Real,
)
    pd = _pressure_derivatives_order2_taylordiff(
        model,
        T_fm,
        mu_fm;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        series_iterations=series_iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )
    md = _mass_derivatives_taylordiff(
        T_fm,
        mu_fm;
        order=1,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        model=model,
        series_iterations=series_iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )

    T_val = Float64(T_fm)
    s = pd.P_T
    n_B = pd.P_mu / 3
    ds_dT = pd.P_TT
    ds_dμq = pd.P_Tmu
    dn_dT = pd.P_Tmu / 3
    dn_dμq = pd.P_mumu / 3

    dM_dμB = md.dM_dmu ./ 3.0
    ds_dμB = ds_dμq / 3.0
    dn_dμB = dn_dμq / 3.0

    numerator_vn = s * dn_dμB - n_B * dn_dT
    denominator_vn = T_val * (ds_dT * dn_dμB - ds_dμB * dn_dT)
    v_n_sq = numerator_vn / denominator_vn

    num_muB_T = n_B * ds_dT - s * dn_dT
    den_muB_T = n_B * ds_dμB - s * dn_dμB
    dμB_dT_sig = -num_muB_T / den_muB_T

    c_p = abs(n_B) <= sqrt(eps(Float64)) ? NaN : T_val * (ds_dT - ds_dμB * s / n_B)

    return (
        v_n_sq=v_n_sq,
        dμB_dT_sigma=dμB_dT_sig,
        masses=md.masses,
        dM_dT=md.dM_dT,
        dM_dμB=dM_dμB,
        ds_dT=ds_dT,
        ds_dμB=ds_dμB,
        dn_dT=dn_dT,
        dn_dμB=dn_dμB,
        c_p=c_p,
        s=s,
        n_B=n_B,
    )
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
                          model=nothing,
                          derivative_backend::Symbol=:auto,
                          series_iterations::Union{Nothing, Int}=nothing,
                          linear_solve::Symbol=:auto,
                          series_residual_tol::Real=1e-7)
    set_config(xi=xi, p_num=p_num, t_num=t_num)
    resolved_model, backend = _resolve_thermo_backend(model, derivative_backend)
    backend === :taylordiff || error("unreachable derivative backend: $backend")
    return _mass_derivatives_taylordiff(
        T_fm,
        mu_fm;
        order=order,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        model=resolved_model,
        series_iterations=series_iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )
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
                            model=nothing,
                            derivative_backend::Symbol=:auto,
                            series_iterations::Union{Nothing, Int}=nothing,
                            linear_solve::Symbol=:auto,
                            series_residual_tol::Real=1e-7)
    set_config(xi=xi, p_num=p_num, t_num=t_num)
    resolved_model, backend = _resolve_thermo_backend(model, derivative_backend)
    backend === :taylordiff || error("unreachable derivative backend: $backend")
    return _thermo_derivatives_taylordiff(
        T_fm,
        mu_fm;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        model=resolved_model,
        series_iterations=series_iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
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

"""
    dP_dT(T_fm, mu_fm; order=1, kwargs...)

总压强对温度的 n 阶导数。
"""
function dP_dT(T_fm::Real, mu_fm::Real; order::Int=1, xi::Real=0.0,
               p_num::Int=DEFAULT_MOMENTUM_COUNT, t_num::Int=DEFAULT_THETA_COUNT,
               model=nothing,
               derivative_backend::Symbol=:auto,
               series_iterations::Union{Nothing, Int}=nothing,
               linear_solve::Symbol=:auto,
               series_residual_tol::Real=1e-7)
    set_config(xi=xi, p_num=p_num, t_num=t_num)
    resolved_model, backend = _resolve_thermo_backend(model, derivative_backend)
    backend === :taylordiff || error("unreachable derivative backend: $backend")
    order >= 1 || throw(ArgumentError("order must be >= 1, got $order"))
    series = _pressure_series_taylordiff(
        resolved_model,
        T_fm,
        mu_fm,
        1.0,
        _zero_mu_direction();
        order=order,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        series_iterations=series_iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )
    return nth_derivative_from_series(series, order)
end

"""
    dP_dmu(T_fm, mu_fm; order=1, kwargs...)

总压强对化学势的 n 阶导数。
"""
function dP_dmu(T_fm::Real, mu_fm::Real; order::Int=1, xi::Real=0.0,
                p_num::Int=DEFAULT_MOMENTUM_COUNT, t_num::Int=DEFAULT_THETA_COUNT,
                model=nothing,
                derivative_backend::Symbol=:auto,
                series_iterations::Union{Nothing, Int}=nothing,
                linear_solve::Symbol=:auto,
                series_residual_tol::Real=1e-7)
    set_config(xi=xi, p_num=p_num, t_num=t_num)
    resolved_model, backend = _resolve_thermo_backend(model, derivative_backend)
    backend === :taylordiff || error("unreachable derivative backend: $backend")
    order >= 1 || throw(ArgumentError("order must be >= 1, got $order"))
    series = _pressure_series_taylordiff(
        resolved_model,
        T_fm,
        mu_fm,
        0.0,
        _symmetric_mu_direction();
        order=order,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        series_iterations=series_iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )
    return nth_derivative_from_series(series, order)
end

# ============================================================================
# 体粘滞系数
# ============================================================================

"""
    bulk_viscosity_coefficients(T_fm, mu_fm; kwargs...)

计算体粘滞系数公式所需的所有热力学导数系数。

使用 TaylorDiff explicit Taylor-series gap Newton 实现导数传播。

## 技术说明

默认后端从压强与质量的 Taylor series 中解析装配 `ds/dT`、
`ds/dμB`、`dn/dT`、`dn/dμB` 等组合，不再调用旧 `ImplicitFunction`
fallback。
"""
function bulk_viscosity_coefficients(T_fm::Real, mu_fm::Real;
                                     xi::Real=0.0,
                                     p_num::Int=DEFAULT_MOMENTUM_COUNT,
                                     t_num::Int=DEFAULT_THETA_COUNT,
                                     model=nothing,
                                     derivative_backend::Symbol=:auto,
                                     series_iterations::Union{Nothing, Int}=nothing,
                                     linear_solve::Symbol=:auto,
                                     series_residual_tol::Real=1e-7)
    set_config(xi=xi, p_num=p_num, t_num=t_num)
    resolved_model, backend = _resolve_thermo_backend(model, derivative_backend)
    backend === :taylordiff || error("unreachable derivative backend: $backend")
    return _bulk_viscosity_coefficients_taylordiff(
        T_fm,
        mu_fm;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        model=resolved_model,
        series_iterations=series_iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
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
