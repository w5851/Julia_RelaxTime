"""
    ImplicitSolver

PNJL 隐函数求解器模块，整合 ImplicitDifferentiation.jl。

## 主要功能
- 统一的求解接口（支持多种约束模式）
- 基于 ImplicitDifferentiation.jl 的自动微分支持
- 灵活的初值策略

## 使用示例
```julia
# 固定化学势求解
result = solve(FixedMu(), T_fm, μ_fm)

# 固定密度求解
result = solve(FixedRho(1.0), T_fm)

# 使用自定义初值策略
result = solve(FixedMu(), T_fm, μ_fm; seed_strategy=MultiSeed())
```
"""
module ImplicitSolver

using StaticArrays
using NLsolve
using ForwardDiff
using ImplicitDifferentiation

# 使用相对路径导入，避免重复定义
using ..ConstraintModes: ConstraintMode, FixedMu, FixedRho, FixedEntropy, FixedSigma, state_dim, param_dim
using ..SeedStrategies: SeedStrategy, DefaultSeed, MultiSeed, ContinuitySeed, get_seed, get_all_seeds, default_omega_selector
using ..Conditions: GapParams, gap_conditions, build_residual!

# 导入 core 模块的 Thermodynamics
const _THERMO_PATH = normpath(joinpath(@__DIR__, "..", "core", "Thermodynamics.jl"))
if !isdefined(@__MODULE__, :Thermodynamics)
    include(_THERMO_PATH)
end
using .Thermodynamics: calculate_pressure, calculate_omega, calculate_rho, calculate_thermo, calculate_mass_vec, ρ0
using .Thermodynamics.Integrals: cached_nodes, DEFAULT_MOMENTUM_COUNT, DEFAULT_THETA_COUNT

export solve, SolverResult
export create_implicit_solver, solve_with_derivatives

# ============================================================================
# 物理性判据与兜底求解（Newton → Trust-Region）
# ============================================================================

@inline function _default_is_physical_solution(x_state::SVector{5, Float64}, masses::SVector{3, Float64}; phi_tol::Float64=1e-8)
    Φ = x_state[4]
    Φbar = x_state[5]
    if !(isfinite(Φ) && isfinite(Φbar) && (-phi_tol <= Φ <= 1 + phi_tol) && (-phi_tol <= Φbar <= 1 + phi_tol))
        return false
    end
    if any(!isfinite, masses) || any(m -> m <= 0.0, masses)
        return false
    end
    return true
end

@inline function _all_finite_thermo(omega::Float64, pressure::Float64, rho_norm::Float64, entropy::Float64, energy::Float64)
    return isfinite(omega) && isfinite(pressure) && isfinite(rho_norm) && isfinite(entropy) && isfinite(energy)
end

function _postprocess_candidate(postprocess_fn::Function, physicality_check::Function, x_sol)
    pp = postprocess_fn(x_sol)
    phys = physicality_check(pp.x_state, pp.masses) && _all_finite_thermo(pp.omega, pp.pressure, pp.rho_norm, pp.entropy, pp.energy)
    return (phys=phys, x_sol=Vector{Float64}(x_sol), pp...)
end

function _choose_candidate(primary_res, primary_cand, fallback_res, fallback_cand; residual_norm_max::Float64)
    primary_good = primary_res.f_converged && isfinite(primary_res.residual_norm) && primary_res.residual_norm <= residual_norm_max && primary_cand.phys
    fallback_good = fallback_res.f_converged && isfinite(fallback_res.residual_norm) && fallback_res.residual_norm <= residual_norm_max && fallback_cand.phys

    if fallback_good && !primary_good
        return fallback_res, fallback_cand
    elseif primary_good && !fallback_good
        return primary_res, primary_cand
    elseif fallback_good && primary_good
        # 同样“好”的情况下：优先 omega 更小（P 更大）；再比 residual_norm
        if fallback_cand.omega < primary_cand.omega
            return fallback_res, fallback_cand
        elseif fallback_cand.omega > primary_cand.omega
            return primary_res, primary_cand
        else
            return (fallback_res.residual_norm < primary_res.residual_norm) ? (fallback_res, fallback_cand) : (primary_res, primary_cand)
        end
    end

    # 都不够好：优先收敛；否则 residual 更小
    if fallback_res.f_converged && !primary_res.f_converged
        return fallback_res, fallback_cand
    elseif primary_res.f_converged && !fallback_res.f_converged
        return primary_res, primary_cand
    end
    if isfinite(fallback_res.residual_norm) && isfinite(primary_res.residual_norm)
        return (fallback_res.residual_norm < primary_res.residual_norm) ? (fallback_res, fallback_cand) : (primary_res, primary_cand)
    end
    return primary_res, primary_cand
end

function _nlsolve_with_tr_fallback(residual_fn!, x0;
    primary_method::Symbol,
    fallback_method::Symbol=:trust_region,
    use_fallback::Bool=true,
    physicality_check::Function=_default_is_physical_solution,
    residual_norm_max::Float64=1e-6,
    postprocess_fn::Function,
    nlsolve_kwargs...)

    primary_res = nlsolve(residual_fn!, x0; autodiff=:forward, method=primary_method, xtol=1e-9, ftol=1e-9, nlsolve_kwargs...)

    local primary_cand
    try
        primary_cand = _postprocess_candidate(postprocess_fn, physicality_check, primary_res.zero)
    catch
        primary_cand = (phys=false,
                        x_sol=Vector{Float64}(primary_res.zero),
                        x_state=SVector{5, Float64}(fill(NaN, 5)),
                        mu_vec=SVector{3, Float64}(fill(NaN, 3)),
                        omega=NaN,
                        pressure=NaN,
                        rho_norm=NaN,
                        entropy=NaN,
                        energy=NaN,
                        masses=SVector{3, Float64}(fill(NaN, 3)))
    end

    need_fallback = use_fallback && (
        !primary_res.f_converged ||
        !isfinite(primary_res.residual_norm) ||
        primary_res.residual_norm > residual_norm_max ||
        !primary_cand.phys
    )

    if !need_fallback
        return primary_res, primary_cand
    end

    local fallback_res
    local fallback_cand
    try
        fallback_res = nlsolve(residual_fn!, x0; autodiff=:forward, method=fallback_method, xtol=1e-9, ftol=1e-9, nlsolve_kwargs...)
        fallback_cand = _postprocess_candidate(postprocess_fn, physicality_check, fallback_res.zero)
    catch
        return primary_res, primary_cand
    end

    return _choose_candidate(primary_res, primary_cand, fallback_res, fallback_cand; residual_norm_max=residual_norm_max)
end

# ============================================================================
# 求解结果结构
# ============================================================================

"""
    SolverResult

求解结果结构体。

# 字段
- `mode::ConstraintMode`: 求解模式
- `converged::Bool`: 是否收敛
- `solution::Vector{Float64}`: 解向量
- `x_state::SVector{5, Float64}`: 状态变量 [φ_u, φ_d, φ_s, Φ, Φ̄]
- `mu_vec::SVector{3, Float64}`: 化学势 [μ_u, μ_d, μ_s]
- `omega::Float64`: 巨热力学势 Ω
- `pressure::Float64`: 压强 P
- `rho_norm::Float64`: 归一化密度 ρ/ρ₀
- `entropy::Float64`: 熵密度 s
- `energy::Float64`: 能量密度 ε
- `masses::SVector{3, Float64}`: 有效质量 [M_u, M_d, M_s]
- `iterations::Int`: 迭代次数
- `residual_norm::Float64`: 残差范数
- `xi::Float64`: 各向异性参数
"""
struct SolverResult
    mode::ConstraintMode
    converged::Bool
    solution::Vector{Float64}
    x_state::SVector{5, Float64}
    mu_vec::SVector{3, Float64}
    omega::Float64
    pressure::Float64
    rho_norm::Float64
    entropy::Float64
    energy::Float64
    masses::SVector{3, Float64}
    iterations::Int
    residual_norm::Float64
    xi::Float64
end

# ============================================================================
# 核心求解函数
# ============================================================================

"""
    solve(mode::FixedMu, T_fm, μ_fm; kwargs...) -> SolverResult

固定化学势模式求解。

# 参数
- `T_fm`: 温度 (fm⁻¹)
- `μ_fm`: 夸克化学势 (fm⁻¹)
- `xi`: 各向异性参数（默认 0.0）
- `seed_strategy`: 初值策略（默认 DefaultSeed()）
- `p_num`, `t_num`: 积分节点数
"""
function solve(::FixedMu, T_fm::Real, μ_fm::Real;
               xi::Real=0.0,
               seed_strategy::SeedStrategy=DefaultSeed(),
               p_num::Int=DEFAULT_MOMENTUM_COUNT,
               t_num::Int=DEFAULT_THETA_COUNT,
               nlsolve_method::Symbol=:newton,
               trust_region_fallback::Bool=true,
               auto_multiseed_fallback::Bool=true,
               fallback_method::Symbol=:trust_region,
               physicality_check::Function=_default_is_physical_solution,
               residual_norm_max::Real=1e-6,
               nlsolve_kwargs...)
    
    mode = FixedMu()

    # 如果用户显式给了 MultiSeed，则直接走 solve_multi（符合文档示例），并避免在内部递归 fallback。
    if seed_strategy isa MultiSeed
        return solve_multi(mode, T_fm, μ_fm;
            seed_strategy=seed_strategy,
            nlsolve_method=nlsolve_method,
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            trust_region_fallback=trust_region_fallback,
            auto_multiseed_fallback=false,
            fallback_method=fallback_method,
            physicality_check=physicality_check,
            residual_norm_max=residual_norm_max,
            nlsolve_kwargs...)
    end
    thermal_nodes = cached_nodes(p_num, t_num)
    params = GapParams(Float64(T_fm), thermal_nodes, Float64(xi))
    mu_vec = SVector{3}(μ_fm, μ_fm, μ_fm)
    
    # 获取初值
    θ = [T_fm, μ_fm]
    seed = get_seed(seed_strategy, θ, mode)
    x0 = Float64.(seed)
    
    # 构建残差函数并求解
    residual_fn! = build_residual!(mode, mu_vec, params)
    postprocess_fn = x_sol -> begin
        x_state = SVector{5}(Tuple(x_sol))
        pressure, rho_norm, entropy, energy = calculate_thermo(x_state, mu_vec, T_fm, thermal_nodes, xi)
        omega = -pressure
        masses = calculate_mass_vec(x_state)
        return (x_state=x_state, mu_vec=mu_vec, omega=omega, pressure=pressure, rho_norm=rho_norm, entropy=entropy, energy=energy, masses=masses)
    end
    res, cand = _nlsolve_with_tr_fallback(residual_fn!, x0;
        primary_method=nlsolve_method,
        fallback_method=fallback_method,
        use_fallback=trust_region_fallback,
        physicality_check=physicality_check,
        residual_norm_max=Float64(residual_norm_max),
        postprocess_fn=postprocess_fn,
        nlsolve_kwargs...)

    converged = res.f_converged && cand.phys && isfinite(res.residual_norm) && (res.residual_norm <= Float64(residual_norm_max))

    single = SolverResult(
        mode,
        converged,
        cand.x_sol,
        cand.x_state,
        mu_vec,
        cand.omega,
        cand.pressure,
        cand.rho_norm,
        cand.entropy,
        cand.energy,
        cand.masses,
        res.iterations,
        res.residual_norm,
        Float64(xi),
    )

    if converged || !auto_multiseed_fallback
        return single
    end

    # 单初值没能得到“可用物理解”时，自动回退到多初值策略。
    try
        multi = solve_multi(mode, T_fm, μ_fm;
            seed_strategy=MultiSeed(),
            nlsolve_method=nlsolve_method,
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            trust_region_fallback=trust_region_fallback,
            auto_multiseed_fallback=false,
            fallback_method=fallback_method,
            physicality_check=physicality_check,
            residual_norm_max=residual_norm_max,
            nlsolve_kwargs...)
        return multi
    catch
        return single
    end
end

"""
    solve(mode::FixedRho, T_fm; kwargs...) -> SolverResult

固定密度模式求解。

# 参数
- `T_fm`: 温度 (fm⁻¹)
- `xi`: 各向异性参数（默认 0.0）
- `seed_strategy`: 初值策略（默认 DefaultSeed()）
"""
function solve(mode::FixedRho, T_fm::Real;
               xi::Real=0.0,
               seed_strategy::SeedStrategy=DefaultSeed(),
               p_num::Int=DEFAULT_MOMENTUM_COUNT,
               t_num::Int=DEFAULT_THETA_COUNT,
               nlsolve_method::Symbol=:newton,
               trust_region_fallback::Bool=true,
               fallback_method::Symbol=:trust_region,
               physicality_check::Function=_default_is_physical_solution,
               residual_norm_max::Real=1e-6,
               nlsolve_kwargs...)
    
    thermal_nodes = cached_nodes(p_num, t_num)
    params = GapParams(Float64(T_fm), thermal_nodes, Float64(xi))
    
    # 获取初值
    θ = [T_fm]
    seed = get_seed(seed_strategy, θ, mode)
    x0 = Float64.(seed)
    
    # 构建残差函数并求解
    residual_fn! = build_residual!(mode, params)
    postprocess_fn = x_sol -> begin
        x_state = SVector{5}(Tuple(x_sol[1:5]))
        mu_vec = SVector{3}(x_sol[6], x_sol[7], x_sol[8])
        pressure, rho_norm, entropy, energy = calculate_thermo(x_state, mu_vec, T_fm, thermal_nodes, xi)
        omega = -pressure
        masses = calculate_mass_vec(x_state)
        return (x_state=x_state, mu_vec=mu_vec, omega=omega, pressure=pressure, rho_norm=rho_norm, entropy=entropy, energy=energy, masses=masses)
    end
    res, cand = _nlsolve_with_tr_fallback(residual_fn!, x0;
        primary_method=nlsolve_method,
        fallback_method=fallback_method,
        use_fallback=trust_region_fallback,
        physicality_check=physicality_check,
        residual_norm_max=Float64(residual_norm_max),
        postprocess_fn=postprocess_fn,
        nlsolve_kwargs...)

    converged = res.f_converged && cand.phys && isfinite(res.residual_norm) && (res.residual_norm <= Float64(residual_norm_max))
    
    return SolverResult(
        mode,
        converged,
        cand.x_sol,
        cand.x_state,
        cand.mu_vec,
        cand.omega,
        cand.pressure,
        cand.rho_norm,
        cand.entropy,
        cand.energy,
        cand.masses,
        res.iterations,
        res.residual_norm,
        Float64(xi),
    )
end

"""
    solve(mode::FixedEntropy, T_fm; kwargs...) -> SolverResult

固定熵密度模式求解。
"""
function solve(mode::FixedEntropy, T_fm::Real;
               xi::Real=0.0,
               seed_strategy::SeedStrategy=DefaultSeed(),
               p_num::Int=DEFAULT_MOMENTUM_COUNT,
               t_num::Int=DEFAULT_THETA_COUNT,
               nlsolve_method::Symbol=:newton,
               trust_region_fallback::Bool=true,
               fallback_method::Symbol=:trust_region,
               physicality_check::Function=_default_is_physical_solution,
               residual_norm_max::Real=1e-6,
               nlsolve_kwargs...)
    
    thermal_nodes = cached_nodes(p_num, t_num)
    params = GapParams(Float64(T_fm), thermal_nodes, Float64(xi))
    
    θ = [T_fm]
    seed = get_seed(seed_strategy, θ, mode)
    x0 = Float64.(seed)
    
    residual_fn! = build_residual!(mode, params)
    postprocess_fn = x_sol -> begin
        x_state = SVector{5}(Tuple(x_sol[1:5]))
        mu_vec = SVector{3}(x_sol[6], x_sol[7], x_sol[8])
        pressure, rho_norm, entropy, energy = calculate_thermo(x_state, mu_vec, T_fm, thermal_nodes, xi)
        omega = -pressure
        masses = calculate_mass_vec(x_state)
        return (x_state=x_state, mu_vec=mu_vec, omega=omega, pressure=pressure, rho_norm=rho_norm, entropy=entropy, energy=energy, masses=masses)
    end
    res, cand = _nlsolve_with_tr_fallback(residual_fn!, x0;
        primary_method=nlsolve_method,
        fallback_method=fallback_method,
        use_fallback=trust_region_fallback,
        physicality_check=physicality_check,
        residual_norm_max=Float64(residual_norm_max),
        postprocess_fn=postprocess_fn,
        nlsolve_kwargs...)

    converged = res.f_converged && cand.phys && isfinite(res.residual_norm) && (res.residual_norm <= Float64(residual_norm_max))
    
    return SolverResult(
        mode,
        converged,
        cand.x_sol,
        cand.x_state,
        cand.mu_vec,
        cand.omega,
        cand.pressure,
        cand.rho_norm,
        cand.entropy,
        cand.energy,
        cand.masses,
        res.iterations,
        res.residual_norm,
        Float64(xi),
    )
end

"""
    solve(mode::FixedSigma, T_fm; kwargs...) -> SolverResult

固定比熵模式求解。
"""
function solve(mode::FixedSigma, T_fm::Real;
               xi::Real=0.0,
               seed_strategy::SeedStrategy=DefaultSeed(),
               p_num::Int=DEFAULT_MOMENTUM_COUNT,
               t_num::Int=DEFAULT_THETA_COUNT,
               nlsolve_method::Symbol=:newton,
               trust_region_fallback::Bool=true,
               fallback_method::Symbol=:trust_region,
               physicality_check::Function=_default_is_physical_solution,
               residual_norm_max::Real=1e-6,
               nlsolve_kwargs...)
    
    thermal_nodes = cached_nodes(p_num, t_num)
    params = GapParams(Float64(T_fm), thermal_nodes, Float64(xi))
    
    θ = [T_fm]
    seed = get_seed(seed_strategy, θ, mode)
    x0 = Float64.(seed)
    
    residual_fn! = build_residual!(mode, params)
    postprocess_fn = x_sol -> begin
        x_state = SVector{5}(Tuple(x_sol[1:5]))
        mu_vec = SVector{3}(x_sol[6], x_sol[7], x_sol[8])
        pressure, rho_norm, entropy, energy = calculate_thermo(x_state, mu_vec, T_fm, thermal_nodes, xi)
        omega = -pressure
        masses = calculate_mass_vec(x_state)
        return (x_state=x_state, mu_vec=mu_vec, omega=omega, pressure=pressure, rho_norm=rho_norm, entropy=entropy, energy=energy, masses=masses)
    end
    res, cand = _nlsolve_with_tr_fallback(residual_fn!, x0;
        primary_method=nlsolve_method,
        fallback_method=fallback_method,
        use_fallback=trust_region_fallback,
        physicality_check=physicality_check,
        residual_norm_max=Float64(residual_norm_max),
        postprocess_fn=postprocess_fn,
        nlsolve_kwargs...)

    converged = res.f_converged && cand.phys && isfinite(res.residual_norm) && (res.residual_norm <= Float64(residual_norm_max))
    
    return SolverResult(
        mode,
        converged,
        cand.x_sol,
        cand.x_state,
        cand.mu_vec,
        cand.omega,
        cand.pressure,
        cand.rho_norm,
        cand.entropy,
        cand.energy,
        cand.masses,
        res.iterations,
        res.residual_norm,
        Float64(xi),
    )
end

# ============================================================================
# 多初值求解
# ============================================================================

"""
    solve_multi(mode, args...; seed_strategy::MultiSeed, kwargs...) -> SolverResult

使用多初值策略求解，返回最优解。
"""
function solve_multi(mode::FixedMu, T_fm::Real, μ_fm::Real;
                     seed_strategy::MultiSeed=MultiSeed(),
                     nlsolve_method::Symbol=:newton,
                     kwargs...)
    θ = [T_fm, μ_fm]
    seeds = get_all_seeds(seed_strategy, θ, mode)
    
    results = SolverResult[]
    for seed in seeds
        try
            result = solve(mode, T_fm, μ_fm; 
                          seed_strategy=DefaultSeed(seed, seed, :hadron),
                          nlsolve_method=nlsolve_method,
                          auto_multiseed_fallback=false,
                          kwargs...)
            push!(results, result)
        catch e
            @warn "Solve failed with seed" seed exception=e
        end
    end

    isempty(results) && error("All seeds failed (exceptions) in solve_multi")

    converged = filter(r -> r.converged, results)
    isempty(converged) && error("All seeds failed to converge to a physical solution")

    return seed_strategy.selector(converged)
end

export solve_multi

# ============================================================================
# ImplicitDifferentiation.jl 集成
# ============================================================================

# 全局配置存储
const IMPLICIT_CONFIG = Ref{NamedTuple}((
    xi = 0.0,
    p_num = DEFAULT_MOMENTUM_COUNT,
    t_num = DEFAULT_THETA_COUNT,
    thermal_nodes = nothing,
))

"""设置隐函数求解器配置"""
function set_implicit_config(; xi::Real=0.0, p_num::Int=DEFAULT_MOMENTUM_COUNT, t_num::Int=DEFAULT_THETA_COUNT)
    thermal_nodes = cached_nodes(p_num, t_num)
    IMPLICIT_CONFIG[] = (xi=Float64(xi), p_num=p_num, t_num=t_num, thermal_nodes=thermal_nodes)
end

"""
前向求解函数（ImplicitDifferentiation.jl 接口）

θ = [T, μ] -> x = [φ_u, φ_d, φ_s, Φ, Φ̄]
"""
function forward_solve_mu(θ::AbstractVector)
    T_fm = Float64(θ[1])
    μ_fm = Float64(θ[2])
    config = IMPLICIT_CONFIG[]
    
    mu_vec = SVector{3}(μ_fm, μ_fm, μ_fm)
    params = GapParams(T_fm, config.thermal_nodes, config.xi)
    
    seed = get_seed(DefaultSeed(), θ, FixedMu())
    residual_fn! = build_residual!(FixedMu(), mu_vec, params)
    res = nlsolve(residual_fn!, seed; autodiff=:forward, method=:newton, xtol=1e-9, ftol=1e-9)
    
    return (res.zero, nothing)
end

"""
条件函数（ImplicitDifferentiation.jl 接口）
"""
function conditions_mu(θ::AbstractVector, x::AbstractVector, z)
    T_fm = θ[1]
    μ_fm = θ[2]
    config = IMPLICIT_CONFIG[]
    
    mu_vec = SVector{3}(μ_fm, μ_fm, μ_fm)
    params = GapParams(T_fm, config.thermal_nodes, config.xi)
    x_state = SVector{5}(Tuple(x))
    
    return Vector(gap_conditions(x_state, mu_vec, params))
end

"""
    create_implicit_solver(; kwargs...) -> ImplicitFunction

创建支持自动微分的隐函数求解器。

# 返回
ImplicitFunction 对象，可用于计算解及其对参数的导数。

# 示例
```julia
solver = create_implicit_solver(xi=0.0)
θ = [T_fm, μ_fm]
x, _ = solver(θ)  # 求解
# 使用 ForwardDiff 计算导数
dx_dθ = ForwardDiff.jacobian(θ -> solver(θ)[1], θ)
```
"""
function create_implicit_solver(; xi::Real=0.0, p_num::Int=DEFAULT_MOMENTUM_COUNT, t_num::Int=DEFAULT_THETA_COUNT)
    set_implicit_config(xi=xi, p_num=p_num, t_num=t_num)
    return ImplicitFunction(
        forward_solve_mu,
        conditions_mu;
        linear_solver=DirectLinearSolver(),
        representation=MatrixRepresentation(),
    )
end

"""
    solve_with_derivatives(T_fm, μ_fm; order=1, kwargs...) -> NamedTuple

求解并计算解对参数的导数。

# 参数
- `T_fm`: 温度 (fm⁻¹)
- `μ_fm`: 化学势 (fm⁻¹)
- `order`: 导数阶数（1 或 2）

# 返回
NamedTuple 包含：
- `x`: 解向量
- `dx_dT`: ∂x/∂T
- `dx_dμ`: ∂x/∂μ
- `d2x_dT2`, `d2x_dμ2`, `d2x_dTdμ`（order=2 时）
"""
function solve_with_derivatives(T_fm::Real, μ_fm::Real;
                                order::Int=1,
                                xi::Real=0.0,
                                p_num::Int=DEFAULT_MOMENTUM_COUNT,
                                t_num::Int=DEFAULT_THETA_COUNT)
    solver = create_implicit_solver(xi=xi, p_num=p_num, t_num=t_num)
    θ = [Float64(T_fm), Float64(μ_fm)]
    
    # 基础解
    x, _ = solver(θ)
    
    if order == 1
        # 一阶导数
        dx_dT = ForwardDiff.derivative(T -> solver([T, θ[2]])[1], θ[1])
        dx_dμ = ForwardDiff.derivative(μ -> solver([θ[1], μ])[1], θ[2])
        return (x=x, dx_dT=dx_dT, dx_dμ=dx_dμ)
        
    elseif order == 2
        # 一阶导数
        dx_dT = ForwardDiff.derivative(T -> solver([T, θ[2]])[1], θ[1])
        dx_dμ = ForwardDiff.derivative(μ -> solver([θ[1], μ])[1], θ[2])
        
        # 二阶导数
        d2x_dT2 = ForwardDiff.derivative(
            T -> ForwardDiff.derivative(t -> solver([t, θ[2]])[1], T),
            θ[1]
        )
        d2x_dμ2 = ForwardDiff.derivative(
            μ -> ForwardDiff.derivative(m -> solver([θ[1], m])[1], μ),
            θ[2]
        )
        d2x_dTdμ = ForwardDiff.derivative(
            T -> ForwardDiff.derivative(μ -> solver([T, μ])[1], θ[2]),
            θ[1]
        )
        
        return (x=x, dx_dT=dx_dT, dx_dμ=dx_dμ, 
                d2x_dT2=d2x_dT2, d2x_dμ2=d2x_dμ2, d2x_dTdμ=d2x_dTdμ)
    else
        error("order must be 1 or 2, got $order")
    end
end

end # module ImplicitSolver
