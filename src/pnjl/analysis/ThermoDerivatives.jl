"""
    ThermoDerivatives

热力学导数计算模块，使用 ImplicitDifferentiation.jl 实现隐函数定理。

## 核心思想

对于非线性方程 F(x; θ) = 0，其中 x 是状态变量，θ = (T, μ) 是参数，
ImplicitDifferentiation.jl 自动使用隐函数定理计算 dx/dθ，支持任意阶导数。

## 主要函数

- `solve_equilibrium_mu`: 求解能隙方程
- `mass_derivatives`: 计算质量及其对 T/μ 的导数（支持任意阶）
- `thermo_derivatives`: 计算热力学量及其导数
- `bulk_derivative_coeffs`: 计算体粘滞系数所需的导数组合

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
using NLsolve
using ImplicitDifferentiation
using LinearAlgebra: dot

using ..AnisoGapSolver:
    cached_nodes,
    calculate_thermo,
    calculate_rho,
    calculate_mass_vec,
    calculate_energy_isotropic,
    residual_mu!,
    DEFAULT_MU_GUESS,
    DEFAULT_MOMENTUM_COUNT,
    DEFAULT_THETA_COUNT

using ..Constants_PNJL: G_fm2, K_fm5

export solve_equilibrium_mu, thermo_derivatives, mass_derivatives, bulk_derivative_coeffs
export quasiparticle_energy, dE_dT, dE_dmu
export bulk_viscosity_coefficients, compute_B_bracket, v_n_squared, dmuB_dT_sigma

# ============================================================================
# 全局缓存的积分节点
# ============================================================================

const THERMAL_NODES = Ref{Any}(nothing)
const THERMAL_NODES_CONFIG = Ref{Tuple{Int,Int}}((0, 0))

"""获取或创建缓存的积分节点"""
function get_thermal_nodes(p_num::Int, t_num::Int)
    if THERMAL_NODES[] === nothing || THERMAL_NODES_CONFIG[] != (p_num, t_num)
        THERMAL_NODES[] = cached_nodes(p_num, t_num)
        THERMAL_NODES_CONFIG[] = (p_num, t_num)
    end
    return THERMAL_NODES[]
end

# ============================================================================
# ImplicitDifferentiation.jl 设置
# ============================================================================

# 存储当前求解配置的全局变量（用于 ImplicitFunction）
const CURRENT_XI = Ref{Float64}(0.0)
const CURRENT_P_NUM = Ref{Int}(DEFAULT_MOMENTUM_COUNT)
const CURRENT_T_NUM = Ref{Int}(DEFAULT_THETA_COUNT)

"""
前向求解函数：给定参数 θ = [T, μ]，求解能隙方程并返回 (x, nothing)
其中 x = [φ_u, φ_d, φ_s, Φ, Φbar]
"""
function forward_solve(θ::AbstractVector)
    T_val = Float64(θ[1])
    μ_val = Float64(θ[2])
    mu_vec = SVector{3}(μ_val, μ_val, μ_val)
    thermal_nodes = get_thermal_nodes(CURRENT_P_NUM[], CURRENT_T_NUM[])
    xi = CURRENT_XI[]
    
    seed = [0.2, 0.2, 0.3, 0.1, 0.1]
    f! = (F, x) -> residual_mu!(F, x, mu_vec, T_val, thermal_nodes, xi)
    res = nlsolve(f!, seed; autodiff=:forward, method=:newton, xtol=1e-9, ftol=1e-9)
    
    return (res.zero, nothing)
end

"""
条件函数：返回残差 F(θ, x)，解满足 F = 0
签名必须是 (θ, y, z) -> F
"""
function conditions(θ::AbstractVector, y::AbstractVector, z)
    T_val = θ[1]
    μ_val = θ[2]
    mu_vec = SVector{3}(μ_val, μ_val, μ_val)
    thermal_nodes = get_thermal_nodes(CURRENT_P_NUM[], CURRENT_T_NUM[])
    xi = CURRENT_XI[]
    
    F = similar(y, promote_type(eltype(y), eltype(θ)))
    residual_mu!(F, y, mu_vec, T_val, thermal_nodes, xi)
    return F
end

# 创建全局 ImplicitFunction 实例
# 使用 DirectLinearSolver + MatrixRepresentation 以支持高阶导数（嵌套 Dual 数）
const IMPLICIT_SOLVER = ImplicitFunction(
    forward_solve, conditions; 
    linear_solver=DirectLinearSolver(),
    representation=MatrixRepresentation()
)

"""
设置求解配置（在调用 IMPLICIT_SOLVER 之前调用）
"""
function set_solver_config(; xi::Real=0.0, p_num::Int=DEFAULT_MOMENTUM_COUNT, t_num::Int=DEFAULT_THETA_COUNT)
    CURRENT_XI[] = Float64(xi)
    CURRENT_P_NUM[] = p_num
    CURRENT_T_NUM[] = t_num
    # 预热积分节点缓存
    get_thermal_nodes(p_num, t_num)
end

# ============================================================================
# 基础求解函数
# ============================================================================

"""
    solve_equilibrium_mu(T_fm, mu_fm; xi=0.0, seed_state, p_num, t_num, kwargs...)

在指定 (T, μ, ξ) 下求解能隙方程，并返回热力学量。

# 参数
- `T_fm`: 温度 (fm⁻¹)
- `mu_fm`: 夸克化学势 (fm⁻¹)
- `xi`: 各向异性参数
- `seed_state`: 初始猜测（默认 DEFAULT_MU_GUESS）
- `p_num`, `t_num`: 积分节点数

# 返回
NamedTuple 包含：pressure, energy, rho, entropy, x_state, converged 等
"""
function solve_equilibrium_mu(T_fm::Real, mu_fm::Real; 
                              xi::Real=0.0, 
                              seed_state=DEFAULT_MU_GUESS, 
                              p_num::Int=DEFAULT_MOMENTUM_COUNT, 
                              t_num::Int=DEFAULT_THETA_COUNT, 
                              nlsolve_kwargs...)
    thermal_nodes = get_thermal_nodes(p_num, t_num)
    seed_type = promote_type(typeof(T_fm), typeof(mu_fm), Float64)
    seed_vec = collect(seed_type.(seed_state))
    mu_vec = SVector{3}(mu_fm, mu_fm, mu_fm)

    f! = (F, x) -> residual_mu!(F, x, mu_vec, T_fm, thermal_nodes, xi)
    extra = (; nlsolve_kwargs...)
    res = nlsolve(f!, seed_vec; autodiff=:forward, method=:newton, xtol=1e-9, ftol=1e-9, extra...)

    x_state = SVector{5}(Tuple(res.zero))
    pressure, rho_norm, entropy, energy = calculate_thermo(x_state, mu_vec, T_fm, thermal_nodes, xi)
    rho_vec = calculate_rho(x_state, mu_vec, T_fm, thermal_nodes, xi)
    rho_baryon = sum(rho_vec) / 3

    return (
        pressure=pressure,
        energy=energy,
        rho=rho_baryon,
        rho_norm=rho_norm,
        entropy=entropy,
        x_state=x_state,
        mu_vec=mu_vec,
        converged=res.f_converged,
        iterations=res.iterations,
        residual_norm=res.residual_norm,
    )
end

# ============================================================================
# 质量计算函数
# ============================================================================

"""从状态变量计算夸克有效质量"""
function compute_masses_from_state(x::AbstractVector)
    φ_u, φ_d, φ_s = x[1], x[2], x[3]
    # m_i = m_i0 - 4G*φ_i + 2K*φ_j*φ_k
    m_u0 = 0.0055 / 0.197327  # MeV -> fm⁻¹
    m_s0 = 0.140 / 0.197327
    
    m_u = m_u0 - 4G_fm2 * φ_u + 2K_fm5 * φ_d * φ_s
    m_d = m_u0 - 4G_fm2 * φ_d + 2K_fm5 * φ_u * φ_s
    m_s = m_s0 - 4G_fm2 * φ_s + 2K_fm5 * φ_u * φ_d
    
    return SVector{3}(m_u, m_d, m_s)
end

"""从状态变量计算单个味道的质量（用于标量函数求导）"""
function compute_mass_flavor(x::AbstractVector, flavor::Int)
    masses = compute_masses_from_state(x)
    return masses[flavor]
end

# ============================================================================
# 质量导数（使用 ImplicitDifferentiation.jl）
# ============================================================================

"""
    mass_derivatives(T_fm, mu_fm; order=1, xi=0.0, p_num, t_num)

计算夸克有效质量及其对 T/μ 的导数。

使用 ImplicitDifferentiation.jl 自动计算任意阶导数。

# 参数
- `T_fm`: 温度 (fm⁻¹)
- `mu_fm`: 夸克化学势 (fm⁻¹)
- `order`: 导数阶数（1 或 2）
- `xi`: 各向异性参数

# 返回
- order=1: (masses, dM_dT, dM_dmu)
- order=2: (masses, dM_dT, dM_dmu, d2M_dT2, d2M_dTdmu, d2M_dmu2)
"""
function mass_derivatives(T_fm::Real, mu_fm::Real; 
                          order::Int=1,
                          xi::Real=0.0, 
                          p_num::Int=DEFAULT_MOMENTUM_COUNT, 
                          t_num::Int=DEFAULT_THETA_COUNT)
    # 设置求解配置
    set_solver_config(; xi=xi, p_num=p_num, t_num=t_num)
    
    θ = [Float64(T_fm), Float64(mu_fm)]
    
    # 定义质量函数（返回三味质量的向量）
    function mass_vec_func(θ_in)
        (x_out, _) = IMPLICIT_SOLVER(θ_in)
        return collect(compute_masses_from_state(x_out))
    end
    
    # 计算质量
    masses_vec = mass_vec_func(θ)
    masses = SVector{3}(masses_vec...)
    
    if order == 1
        # 一阶导数
        dM_dT = SVector{3}(ForwardDiff.derivative(T -> mass_vec_func([T, θ[2]])[i], θ[1]) for i in 1:3)
        dM_dmu = SVector{3}(ForwardDiff.derivative(μ -> mass_vec_func([θ[1], μ])[i], θ[2]) for i in 1:3)
        
        return (masses=masses, dM_dT=dM_dT, dM_dmu=dM_dmu)
        
    elseif order == 2
        # 一阶导数
        dM_dT = SVector{3}(ForwardDiff.derivative(T -> mass_vec_func([T, θ[2]])[i], θ[1]) for i in 1:3)
        dM_dmu = SVector{3}(ForwardDiff.derivative(μ -> mass_vec_func([θ[1], μ])[i], θ[2]) for i in 1:3)
        
        # 二阶导数
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

# 返回
NamedTuple 包含：
- 基础量：pressure, energy, rho, entropy
- 一阶导数：dP_dT, dP_dmu, dEpsilon_dT, dEpsilon_dmu, dn_dT, dn_dmu
- 组合导数：dP_depsilon_n, dP_dn_epsilon
- 质量导数：masses, dM_dT, dM_dmu
"""
function thermo_derivatives(T_fm::Real, mu_fm::Real; 
                            xi::Real=0.0, 
                            p_num::Int=DEFAULT_MOMENTUM_COUNT, 
                            t_num::Int=DEFAULT_THETA_COUNT)
    # 设置求解配置
    set_solver_config(; xi=xi, p_num=p_num, t_num=t_num)
    thermal_nodes = get_thermal_nodes(p_num, t_num)
    
    θ = [Float64(T_fm), Float64(mu_fm)]
    
    # 定义热力学量函数
    function thermo_func(θ_in)
        (x_out, _) = IMPLICIT_SOLVER(θ_in)
        x_sv = SVector{5}(Tuple(x_out))
        mu_vec = SVector{3}(θ_in[2], θ_in[2], θ_in[2])
        
        P, rho_norm, s, ε = calculate_thermo(x_sv, mu_vec, θ_in[1], thermal_nodes, CURRENT_XI[])
        rho_vec = calculate_rho(x_sv, mu_vec, θ_in[1], thermal_nodes, CURRENT_XI[])
        n = sum(rho_vec) / 3
        
        return [P, ε, n, s]
    end
    
    # 计算基础量
    base_vals = thermo_func(θ)
    pressure, energy, rho, entropy = base_vals[1], base_vals[2], base_vals[3], base_vals[4]
    
    # 计算一阶导数
    P_T = ForwardDiff.derivative(T -> thermo_func([T, θ[2]])[1], θ[1])
    P_mu = ForwardDiff.derivative(μ -> thermo_func([θ[1], μ])[1], θ[2])
    E_T = ForwardDiff.derivative(T -> thermo_func([T, θ[2]])[2], θ[1])
    E_mu = ForwardDiff.derivative(μ -> thermo_func([θ[1], μ])[2], θ[2])
    n_T = ForwardDiff.derivative(T -> thermo_func([T, θ[2]])[3], θ[1])
    n_mu = ForwardDiff.derivative(μ -> thermo_func([θ[1], μ])[3], θ[2])
    
    # 计算组合导数
    denom_eps = E_T * n_mu - E_mu * n_T
    denom_n = n_T * E_mu - n_mu * E_T
    dP_depsilon_n = denom_eps == 0 ? NaN : (P_T * n_mu - P_mu * n_T) / denom_eps
    dP_dn_epsilon = denom_n == 0 ? NaN : (P_T * E_mu - P_mu * E_T) / denom_n
    
    # 计算质量导数
    md = mass_derivatives(T_fm, mu_fm; order=1, xi=xi, p_num=p_num, t_num=t_num)
    
    # 获取收敛信息
    base = solve_equilibrium_mu(T_fm, mu_fm; xi=xi, p_num=p_num, t_num=t_num)
    
    return (
        pressure=pressure,
        energy=energy,
        rho=rho,
        rho_norm=base.rho_norm,
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
        converged=base.converged,
        iterations=base.iterations,
        residual_norm=base.residual_norm,
    )
end

"""
    bulk_derivative_coeffs(T_fm, mu_fm; kwargs...)

返回体粘滞系数公式中常用的导数组合。
"""
function bulk_derivative_coeffs(T_fm::Real, mu_fm::Real; 
                                xi::Real=0.0, 
                                p_num::Int=DEFAULT_MOMENTUM_COUNT, 
                                t_num::Int=DEFAULT_THETA_COUNT)
    derivs = thermo_derivatives(T_fm, mu_fm; xi=xi, p_num=p_num, t_num=t_num)
    return (
        dP_depsilon_n=derivs.dP_depsilon_n,
        dP_dn_epsilon=derivs.dP_dn_epsilon,
        dM_dT=derivs.dM_dT,
        dM_dmu=derivs.dM_dmu,
    )
end

# ============================================================================
# 单粒子能量相关函数
# ============================================================================

"""
    quasiparticle_energy(T_fm, mu_fm, p_fm; flavor=1, xi=0.0, kwargs...)

给定动量 p_fm 和味道，返回各向同性色散关系 E = sqrt(p² + m²)。
"""
function quasiparticle_energy(T_fm::Real, mu_fm::Real, p_fm::Real; 
                              flavor::Int=1, 
                              xi::Real=0.0, 
                              p_num::Int=DEFAULT_MOMENTUM_COUNT, 
                              t_num::Int=DEFAULT_THETA_COUNT)
    1 <= flavor <= 3 || error("flavor must be 1, 2, or 3")
    base = solve_equilibrium_mu(T_fm, mu_fm; xi=xi, p_num=p_num, t_num=t_num)
    masses = calculate_mass_vec(base.x_state)
    return calculate_energy_isotropic(masses[flavor], p_fm)
end

"""
    dE_dT(T_fm, mu_fm, p_fm; m, dM_dT)

单粒子能量对温度的导数（链式法则）。
"""
function dE_dT(T_fm::Real, mu_fm::Real, p_fm::Real; m::Real, dM_dT::Real)
    E = calculate_energy_isotropic(m, p_fm)
    return (m / E) * dM_dT
end

"""
    dE_dmu(T_fm, mu_fm, p_fm; m, dM_dmu)

单粒子能量对化学势的导数（链式法则）。
"""
function dE_dmu(T_fm::Real, mu_fm::Real, p_fm::Real; m::Real, dM_dmu::Real)
    E = calculate_energy_isotropic(m, p_fm)
    return (m / E) * dM_dmu
end

# ============================================================================
# 高阶导数辅助函数
# ============================================================================

"""
    dP_dT(T_fm, mu_fm; order=1, xi=0.0, kwargs...)

总压强对温度的 n 阶导数。
"""
function dP_dT(T_fm::Real, mu_fm::Real; 
               order::Int=1, 
               xi::Real=0.0, 
               p_num::Int=DEFAULT_MOMENTUM_COUNT, 
               t_num::Int=DEFAULT_THETA_COUNT)
    set_solver_config(; xi=xi, p_num=p_num, t_num=t_num)
    thermal_nodes = get_thermal_nodes(p_num, t_num)
    
    function P_func(T)
        θ = [T, Float64(mu_fm)]
        (x_out, _) = IMPLICIT_SOLVER(θ)
        x_sv = SVector{5}(Tuple(x_out))
        mu_vec = SVector{3}(mu_fm, mu_fm, mu_fm)
        P, _, _, _ = calculate_thermo(x_sv, mu_vec, T, thermal_nodes, CURRENT_XI[])
        return P
    end
    
    return _nth_derivative(P_func, Float64(T_fm), order)
end

"""
    dP_dmu(T_fm, mu_fm; order=1, xi=0.0, kwargs...)

总压强对化学势的 n 阶导数。
"""
function dP_dmu(T_fm::Real, mu_fm::Real; 
                order::Int=1, 
                xi::Real=0.0, 
                p_num::Int=DEFAULT_MOMENTUM_COUNT, 
                t_num::Int=DEFAULT_THETA_COUNT)
    set_solver_config(; xi=xi, p_num=p_num, t_num=t_num)
    thermal_nodes = get_thermal_nodes(p_num, t_num)
    
    function P_func(μ)
        θ = [Float64(T_fm), μ]
        (x_out, _) = IMPLICIT_SOLVER(θ)
        x_sv = SVector{5}(Tuple(x_out))
        mu_vec = SVector{3}(μ, μ, μ)
        P, _, _, _ = calculate_thermo(x_sv, mu_vec, T_fm, thermal_nodes, CURRENT_XI[])
        return P
    end
    
    return _nth_derivative(P_func, Float64(mu_fm), order)
end

"""递归计算 n 阶导数"""
function _nth_derivative(f, x, n::Int)
    n >= 1 || error("order must be ≥ 1")
    n == 1 && return ForwardDiff.derivative(f, x)
    inner = y -> _nth_derivative(f, y, n - 1)
    return ForwardDiff.derivative(inner, x)
end

# 导出高阶导数函数
export dP_dT, dP_dmu

# ============================================================================
# 体粘滞系数热力学导数（等熵声速形式）
# ============================================================================

"""
    v_n_squared(T_fm, mu_fm; xi=0.0, p_num, t_num)

计算热力学速度 v_n²，用于体粘滞系数公式。

公式：
    v_n² = (s·∂n_B/∂μ_B - n_B·∂n_B/∂T) / [T·(∂s/∂T·∂n_B/∂μ_B - ∂s/∂μ_B·∂n_B/∂T)]

其中所有导数都是对重子化学势 μ_B 的导数。

# 参数
- `T_fm`: 温度 (fm⁻¹)
- `mu_fm`: 夸克化学势 (fm⁻¹)

# 返回
- v_n²: 热力学速度平方（无量纲）
"""
function v_n_squared(T_fm::Real, mu_fm::Real;
                     xi::Real=0.0,
                     p_num::Int=DEFAULT_MOMENTUM_COUNT,
                     t_num::Int=DEFAULT_THETA_COUNT)
    set_solver_config(; xi=xi, p_num=p_num, t_num=t_num)
    thermal_nodes = get_thermal_nodes(p_num, t_num)
    
    # 定义热力学量函数
    function thermo_funcs(T, μ)
        θ = [T, μ]
        (x_out, _) = IMPLICIT_SOLVER(θ)
        x_sv = SVector{5}(Tuple(x_out))
        mu_vec = SVector{3}(μ, μ, μ)
        _, _, s, _ = calculate_thermo(x_sv, mu_vec, T, thermal_nodes, CURRENT_XI[])
        rho_vec = calculate_rho(x_sv, mu_vec, T, thermal_nodes, CURRENT_XI[])
        n = sum(rho_vec) / 3
        return (s=s, n=n)
    end
    
    # 基础量
    base = thermo_funcs(T_fm, mu_fm)
    s, n_B = base.s, base.n
    
    # 一阶导数（对夸克化学势）
    ds_dT = ForwardDiff.derivative(T -> thermo_funcs(T, mu_fm).s, Float64(T_fm))
    ds_dμq = ForwardDiff.derivative(μ -> thermo_funcs(T_fm, μ).s, Float64(mu_fm))
    dn_dT = ForwardDiff.derivative(T -> thermo_funcs(T, mu_fm).n, Float64(T_fm))
    dn_dμq = ForwardDiff.derivative(μ -> thermo_funcs(T_fm, μ).n, Float64(mu_fm))
    
    # 转换到重子化学势：∂/∂μ_B = (1/3)·∂/∂μ_q
    ds_dμB = ds_dμq / 3.0
    dn_dμB = dn_dμq / 3.0
    
    # v_n² = (s·∂n_B/∂μ_B - n_B·∂n_B/∂T) / [T·(∂s/∂T·∂n_B/∂μ_B - ∂s/∂μ_B·∂n_B/∂T)]
    numerator = s * dn_dμB - n_B * dn_dT
    denominator = T_fm * (ds_dT * dn_dμB - ds_dμB * dn_dT)
    
    return numerator / denominator
end

"""
    dmuB_dT_sigma(T_fm, mu_fm; xi=0.0, p_num, t_num)

计算固定 σ=s/n_B 时重子化学势对温度的导数 ∂μ_B/∂T|_σ。

公式：
    ∂μ_B/∂T|_σ = -∂σ/∂T / ∂σ/∂μ_B

其中：
    ∂σ/∂T = ∂s/∂T/n_B - s·∂n_B/∂T/n_B²
    ∂σ/∂μ_B = ∂s/∂μ_B/n_B - s·∂n_B/∂μ_B/n_B²

# 参数
- `T_fm`: 温度 (fm⁻¹)
- `mu_fm`: 夸克化学势 (fm⁻¹)

# 返回
- ∂μ_B/∂T|_σ (fm⁻¹/fm⁻¹ = 无量纲)
"""
function dmuB_dT_sigma(T_fm::Real, mu_fm::Real;
                       xi::Real=0.0,
                       p_num::Int=DEFAULT_MOMENTUM_COUNT,
                       t_num::Int=DEFAULT_THETA_COUNT)
    set_solver_config(; xi=xi, p_num=p_num, t_num=t_num)
    thermal_nodes = get_thermal_nodes(p_num, t_num)
    
    # 定义热力学量函数
    function thermo_funcs(T, μ)
        θ = [T, μ]
        (x_out, _) = IMPLICIT_SOLVER(θ)
        x_sv = SVector{5}(Tuple(x_out))
        mu_vec = SVector{3}(μ, μ, μ)
        _, _, s, _ = calculate_thermo(x_sv, mu_vec, T, thermal_nodes, CURRENT_XI[])
        rho_vec = calculate_rho(x_sv, mu_vec, T, thermal_nodes, CURRENT_XI[])
        n = sum(rho_vec) / 3
        return (s=s, n=n)
    end
    
    # 基础量
    base = thermo_funcs(T_fm, mu_fm)
    s, n_B = base.s, base.n
    
    # 一阶导数（对夸克化学势）
    ds_dT = ForwardDiff.derivative(T -> thermo_funcs(T, mu_fm).s, Float64(T_fm))
    ds_dμq = ForwardDiff.derivative(μ -> thermo_funcs(T_fm, μ).s, Float64(mu_fm))
    dn_dT = ForwardDiff.derivative(T -> thermo_funcs(T, mu_fm).n, Float64(T_fm))
    dn_dμq = ForwardDiff.derivative(μ -> thermo_funcs(T_fm, μ).n, Float64(mu_fm))
    
    # 转换到重子化学势
    ds_dμB = ds_dμq / 3.0
    dn_dμB = dn_dμq / 3.0
    
    # ∂σ/∂T 和 ∂σ/∂μ_B
    dσ_dT = ds_dT / n_B - s * dn_dT / n_B^2
    dσ_dμB = ds_dμB / n_B - s * dn_dμB / n_B^2
    
    return -dσ_dT / dσ_dμB
end

"""
    bulk_viscosity_coefficients(T_fm, mu_fm; xi=0.0, p_num, t_num)

计算体粘滞系数公式（等熵声速形式）所需的所有热力学导数系数。

使用重子化学势约定，与 Fortran 代码一致。

## 优化说明
此函数经过优化，只调用一次 IMPLICIT_SOLVER，然后通过 ForwardDiff.jacobian
一次性计算所有需要的导数（s, n, M 对 T 和 μ 的导数）。

# 参数
- `T_fm`: 温度 (fm⁻¹)
- `mu_fm`: 夸克化学势 (fm⁻¹)

# 返回
NamedTuple 包含：
- `v_n_sq`: 热力学速度 v_n²
- `dμB_dT_sigma`: ∂μ_B/∂T|_σ
- `masses`: 夸克有效质量 (u, d, s)
- `dM_dT`: 质量对温度的导数
- `dM_dμB`: 质量对重子化学势的导数
- `s`: 熵密度
- `n_B`: 重子数密度
"""
function bulk_viscosity_coefficients(T_fm::Real, mu_fm::Real;
                                     xi::Real=0.0,
                                     p_num::Int=DEFAULT_MOMENTUM_COUNT,
                                     t_num::Int=DEFAULT_THETA_COUNT)
    set_solver_config(; xi=xi, p_num=p_num, t_num=t_num)
    thermal_nodes = get_thermal_nodes(p_num, t_num)
    
    # 定义统一的计算函数：返回 [s, n, M_u, M_d, M_s]
    # 这样只需要一次 jacobian 调用就能得到所有导数
    function all_quantities(θ::AbstractVector)
        (x_out, _) = IMPLICIT_SOLVER(θ)
        x_sv = SVector{5}(Tuple(x_out))
        mu_vec = SVector{3}(θ[2], θ[2], θ[2])
        
        # 热力学量
        _, _, s, _ = calculate_thermo(x_sv, mu_vec, θ[1], thermal_nodes, CURRENT_XI[])
        rho_vec = calculate_rho(x_sv, mu_vec, θ[1], thermal_nodes, CURRENT_XI[])
        n = sum(rho_vec) / 3
        
        # 质量
        masses = compute_masses_from_state(x_out)
        
        return [s, n, masses[1], masses[2], masses[3]]
    end
    
    θ = [Float64(T_fm), Float64(mu_fm)]
    
    # 计算基础量
    base_vals = all_quantities(θ)
    s = base_vals[1]
    n_B = base_vals[2]
    masses = SVector{3}(base_vals[3], base_vals[4], base_vals[5])
    
    # 一次性计算 Jacobian: ∂[s, n, M_u, M_d, M_s]/∂[T, μ]
    # J[i,j] = ∂(quantity_i)/∂(θ_j)
    J = ForwardDiff.jacobian(all_quantities, θ)
    
    # 提取导数（对夸克化学势）
    ds_dT = J[1, 1]
    ds_dμq = J[1, 2]
    dn_dT = J[2, 1]
    dn_dμq = J[2, 2]
    dM_dT = SVector{3}(J[3, 1], J[4, 1], J[5, 1])
    dM_dμq = SVector{3}(J[3, 2], J[4, 2], J[5, 2])
    
    # 转换到重子化学势：∂/∂μ_B = (1/3)·∂/∂μ_q
    ds_dμB = ds_dμq / 3.0
    dn_dμB = dn_dμq / 3.0
    dM_dμB = dM_dμq ./ 3.0
    
    # v_n² = (s·∂n_B/∂μ_B - n_B·∂n_B/∂T) / [T·(∂s/∂T·∂n_B/∂μ_B - ∂s/∂μ_B·∂n_B/∂T)]
    numerator_vn = s * dn_dμB - n_B * dn_dT
    denominator_vn = T_fm * (ds_dT * dn_dμB - ds_dμB * dn_dT)
    v_n_sq = numerator_vn / denominator_vn
    
    # ∂μ_B/∂T|_σ = -∂σ/∂T / ∂σ/∂μ_B
    dσ_dT = ds_dT / n_B - s * dn_dT / n_B^2
    dσ_dμB = ds_dμB / n_B - s * dn_dμB / n_B^2
    dμB_dT_sig = -dσ_dT / dσ_dμB
    
    return (
        v_n_sq = v_n_sq,
        dμB_dT_sigma = dμB_dT_sig,
        masses = masses,
        dM_dT = dM_dT,
        dM_dμB = dM_dμB,
        s = s,
        n_B = n_B,
    )
end

"""
    compute_B_bracket(p, M, μq, T, v_n_sq, dμB_dT_sigma, dM_dT, dM_dμB; is_antiquark=false)

计算体粘滞公式中的 B 项（括号项）。

公式：
    B = p² + 3·v_n²·T²·E·∂[(E∓μ_q)/T]/∂T|_σ

对于夸克：使用 E - μ_q
对于反夸克：使用 E + μ_q

# 参数
- `p`: 动量 (fm⁻¹)
- `M`: 夸克有效质量 (fm⁻¹)
- `μq`: 夸克化学势 (fm⁻¹)
- `T`: 温度 (fm⁻¹)
- `v_n_sq`: 热力学速度 v_n²
- `dμB_dT_sigma`: ∂μ_B/∂T|_σ
- `dM_dT`: ∂M/∂T
- `dM_dμB`: ∂M/∂μ_B（对重子化学势的导数）
- `is_antiquark`: 是否为反夸克

# 返回
- B: 括号项 (fm⁻²)
"""
function compute_B_bracket(p::Real, M::Real, μq::Real, T::Real,
                           v_n_sq::Real, dμB_dT_sigma::Real,
                           dM_dT::Real, dM_dμB::Real;
                           is_antiquark::Bool=false)
    E = sqrt(p^2 + M^2)
    
    # 能量导数（对重子化学势）
    dE_dT_val = (M / E) * dM_dT
    dE_dμB = (M / E) * dM_dμB
    
    # 夸克重子数 b_q = 1/3
    b_q = 1.0 / 3.0
    
    if is_antiquark
        # 反夸克：x = E + μ_q
        # ∂x/∂T|_σ = ∂E/∂T + (∂E/∂μ_B + 1/3)·∂μ_B/∂T|_σ
        dx_dT_sigma = dE_dT_val + (dE_dμB + b_q) * dμB_dT_sigma
        x = E + μq
    else
        # 夸克：x = E - μ_q
        # ∂x/∂T|_σ = ∂E/∂T + (∂E/∂μ_B - 1/3)·∂μ_B/∂T|_σ
        dx_dT_sigma = dE_dT_val + (dE_dμB - b_q) * dμB_dT_sigma
        x = E - μq
    end
    
    # ∂[x/T]/∂T|_σ = (1/T)·∂x/∂T|_σ - x/T²
    dxt_dT_sigma = dx_dT_sigma / T - x / T^2
    
    # B = p² + 3·v_n²·T²·E·∂[x/T]/∂T|_σ
    B = p^2 + 3 * v_n_sq * T^2 * E * dxt_dT_sigma
    
    return B
end

end # module
