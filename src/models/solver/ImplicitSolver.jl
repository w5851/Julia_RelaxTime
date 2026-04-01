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

# 从 Models 域导入，避免重复定义
import Main.Models: ConstraintMode, FixedMu, FixedRho, FixedAsymmetricRho, FixedEntropy, FixedSigma, state_dim, param_dim
import Main.Models: AbstractQCDModel, model_thermo, calculate_mass_vec
import Main.Models: RootPolicy, solve_root_with_policy
using ..SeedStrategies: SeedStrategy, DefaultSeed, MultiSeed, ContinuitySeed, HybridContinuitySeed, PhaseAwareContinuitySeed, get_seed, get_all_seeds, default_omega_selector, update!
using ..Conditions: GapParams, gap_conditions, build_residual!

const cached_nodes = Main.Models.cached_nodes
const DEFAULT_MOMENTUM_COUNT = Main.Models.default_momentum_count()
const DEFAULT_THETA_COUNT = Main.Models.default_theta_count()
using Main.Constants_PNJL: ρ0_inv_fm3
const ρ0 = ρ0_inv_fm3

export solve, SolverResult
export create_implicit_solver, solve_with_derivatives
export solve_weighted_block_fallback
export solve_with_root_diagnostics

@inline function _get_model(model_kind::Symbol)
    if model_kind === :PNJL || model_kind === :RPNJL
        return Main.Models.get_cached_model(model_kind)
    end
    error("Unsupported model kind in ImplicitSolver: $(model_kind)")
end

@inline function _postprocess_payload(model_kind::Symbol, x_state, mu_vec, T_fm;
    p_num::Int,
    t_num::Int,
    xi,
)
    model = _get_model(model_kind)
    pressure, rho_norm, entropy, energy = model_thermo(
        model,
        x_state,
        mu_vec,
        T_fm;
        p_num=p_num,
        t_num=t_num,
        xi=xi,
    )
    omega = -pressure
    masses = calculate_mass_vec(model, SVector{3}(x_state[1], x_state[2], x_state[3]))
    return (x_state=x_state, mu_vec=mu_vec, omega=omega, pressure=pressure, rho_norm=rho_norm, entropy=entropy, energy=energy, masses=masses)
end

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

    cache = Dict{Symbol,Tuple{Any,Any}}()

    solve_once = function (method::Symbol, seed::Vector{Float64})
        res = nlsolve(residual_fn!, seed; autodiff=:forward, method=method, xtol=1e-9, ftol=1e-9, nlsolve_kwargs...)

        local cand
        try
            cand = _postprocess_candidate(postprocess_fn, physicality_check, res.zero)
        catch
            cand = (phys=false,
                    x_sol=Vector{Float64}(res.zero),
                    x_state=SVector{5, Float64}(fill(NaN, 5)),
                    mu_vec=SVector{3, Float64}(fill(NaN, 3)),
                    omega=NaN,
                    pressure=NaN,
                    rho_norm=NaN,
                    entropy=NaN,
                    energy=NaN,
                    masses=SVector{3, Float64}(fill(NaN, 3)))
        end

        cache[method] = (res, cand)
        return (
            x=Vector{Float64}(res.zero),
            converged=Bool(res.f_converged) && Bool(cand.phys),
            residual_norm=Float64(res.residual_norm),
            score=isfinite(cand.omega) ? Float64(cand.omega) : NaN,
        )
    end

    policy = RootPolicy(
        primary_method=primary_method,
        fallback_method=fallback_method,
        use_fallback=use_fallback,
        use_multiseed=false,
        residual_norm_max=residual_norm_max,
        require_converged=true,
        diagnostics_level=:basic,
    )

    solved = solve_root_with_policy(solve_once, Vector{Float64}(x0); policy=policy)
    selected_method = solved.diagnostics.attempts[solved.diagnostics.selected_attempt].method

    if haskey(cache, selected_method)
        return cache[selected_method]
    end

    picked = solve_once(selected_method, Vector{Float64}(x0))
    res, cand = cache[selected_method]
    _ = picked
    return res, cand
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

    # PhaseAwareContinuitySeed 可选：第一个点用 MultiSeed 自举（选 Ω 最小的物理解），随后由调用方 update! 进入连续跟踪。
    if seed_strategy isa PhaseAwareContinuitySeed
        s = seed_strategy::PhaseAwareContinuitySeed
        if s.bootstrap_multiseed && s.previous_solution === nothing
            return solve_multi(mode, T_fm, μ_fm;
                seed_strategy=s.bootstrap_strategy,
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
    end

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
    params = GapParams(Float64(T_fm), thermal_nodes, Float64(xi); p_num=p_num, t_num=t_num, model_kind=:PNJL)
    mu_vec = SVector{3}(μ_fm, μ_fm, μ_fm)
    
    # 获取初值
    θ = [T_fm, μ_fm]
    seed = get_seed(seed_strategy, θ, mode)
    x0 = Float64.(seed)
    
    # 构建残差函数并求解
    residual_fn! = build_residual!(mode, mu_vec, params)
    postprocess_fn = x_sol -> begin
        x_state = SVector{5}(Tuple(x_sol))
        return _postprocess_payload(:PNJL, x_state, mu_vec, T_fm;
            p_num=p_num,
            t_num=t_num,
            xi=xi,
        )
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

    if seed_strategy isa MultiSeed
        return solve_multi(mode, T_fm;
            seed_strategy=seed_strategy,
            nlsolve_method=nlsolve_method,
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            trust_region_fallback=trust_region_fallback,
            fallback_method=fallback_method,
            physicality_check=physicality_check,
            residual_norm_max=residual_norm_max,
            nlsolve_kwargs...)
    end
    
    thermal_nodes = cached_nodes(p_num, t_num)
            params = GapParams(Float64(T_fm), thermal_nodes, Float64(xi); p_num=p_num, t_num=t_num, model_kind=:PNJL)
    
    # 获取初值
    θ = [T_fm]
    seed = get_seed(seed_strategy, θ, mode)
    x0 = Float64.(seed)
    
    # 构建残差函数并求解
    residual_fn! = build_residual!(mode, params)
    postprocess_fn = x_sol -> begin
        x_state = SVector{5}(Tuple(x_sol[1:5]))
        mu_vec = SVector{3}(x_sol[6], x_sol[7], x_sol[8])
        return _postprocess_payload(:PNJL, x_state, mu_vec, T_fm;
            p_num=p_num,
            t_num=t_num,
            xi=xi,
        )
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
    solve(mode::FixedAsymmetricRho, T_fm; kwargs...) -> SolverResult

固定非对称约束模式求解。

# 约束
- `sum(rho)/(3ρ0) = mode.rho_target`
- `rho_u/rho_d = mode.ud_ratio_target`
- `rho_s = mode.s_target`
"""
function solve(mode::FixedAsymmetricRho, T_fm::Real;
               xi::Real=0.0,
               seed_strategy::SeedStrategy=DefaultSeed(),
               p_num::Int=DEFAULT_MOMENTUM_COUNT,
               t_num::Int=DEFAULT_THETA_COUNT,
               model_kind::Symbol=:PNJL,
               nlsolve_method::Symbol=:newton,
               trust_region_fallback::Bool=true,
               fallback_method::Symbol=:trust_region,
               enforce_physicality::Bool=false,
               physicality_check::Function=_default_is_physical_solution,
               residual_norm_max::Real=1e-6,
               nlsolve_kwargs...)

    if seed_strategy isa HybridContinuitySeed
        s = seed_strategy::HybridContinuitySeed

        continuity_result = try
            solve(mode, T_fm;
                xi=xi,
                seed_strategy=s.continuity,
                p_num=p_num,
                t_num=t_num,
                model_kind=model_kind,
                nlsolve_method=nlsolve_method,
                trust_region_fallback=trust_region_fallback,
                fallback_method=fallback_method,
                enforce_physicality=enforce_physicality,
                physicality_check=physicality_check,
                residual_norm_max=residual_norm_max,
                nlsolve_kwargs...)
        catch
            nothing
        end

        if continuity_result !== nothing && continuity_result.converged
            update!(s, continuity_result.solution)
            return continuity_result
        end

        fallback_result = solve(mode, T_fm;
            xi=xi,
            seed_strategy=s.fallback,
            p_num=p_num,
            t_num=t_num,
            model_kind=model_kind,
            nlsolve_method=nlsolve_method,
            trust_region_fallback=trust_region_fallback,
            fallback_method=fallback_method,
            enforce_physicality=enforce_physicality,
            physicality_check=physicality_check,
            residual_norm_max=residual_norm_max,
            nlsolve_kwargs...)

        if fallback_result.converged
            update!(s, fallback_result.solution)
        end
        return fallback_result
    end

    if seed_strategy isa MultiSeed
        return solve_multi(mode, T_fm;
            seed_strategy=seed_strategy,
            nlsolve_method=nlsolve_method,
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            model_kind=model_kind,
            trust_region_fallback=trust_region_fallback,
            fallback_method=fallback_method,
            physicality_check=physicality_check,
            residual_norm_max=residual_norm_max,
            nlsolve_kwargs...)
    end

    thermal_nodes = cached_nodes(p_num, t_num)
    params = GapParams(Float64(T_fm), thermal_nodes, Float64(xi);
        p_num=p_num, t_num=t_num, model_kind=model_kind)

    θ = [T_fm]
    seed = get_seed(seed_strategy, θ, mode)
    x0 = Float64.(seed)

    residual_fn! = build_residual!(mode, params)
    effective_physicality_check = enforce_physicality ? physicality_check : ((_, _) -> true)
    postprocess_fn = x_sol -> begin
        x_state = SVector{5}(Tuple(x_sol[1:5]))
        mu_vec = SVector{3}(x_sol[6], x_sol[7], x_sol[8])
        return _postprocess_payload(model_kind, x_state, mu_vec, T_fm;
            p_num=p_num,
            t_num=t_num,
            xi=xi,
        )
    end

    res, cand = _nlsolve_with_tr_fallback(residual_fn!, x0;
        primary_method=nlsolve_method,
        fallback_method=fallback_method,
        use_fallback=trust_region_fallback,
        physicality_check=effective_physicality_check,
        residual_norm_max=Float64(residual_norm_max),
        postprocess_fn=postprocess_fn,
        nlsolve_kwargs...)

    converged = if enforce_physicality
        res.f_converged && cand.phys && isfinite(res.residual_norm) && (res.residual_norm <= Float64(residual_norm_max))
    else
        isfinite(res.residual_norm) && (res.residual_norm <= Float64(residual_norm_max))
    end

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
    params = GapParams(Float64(T_fm), thermal_nodes, Float64(xi); p_num=p_num, t_num=t_num, model_kind=:PNJL)
    
    θ = [T_fm]
    seed = get_seed(seed_strategy, θ, mode)
    x0 = Float64.(seed)
    
    residual_fn! = build_residual!(mode, params)
    postprocess_fn = x_sol -> begin
        x_state = SVector{5}(Tuple(x_sol[1:5]))
        mu_vec = SVector{3}(x_sol[6], x_sol[7], x_sol[8])
        return _postprocess_payload(:PNJL, x_state, mu_vec, T_fm;
            p_num=p_num,
            t_num=t_num,
            xi=xi,
        )
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
    params = GapParams(Float64(T_fm), thermal_nodes, Float64(xi); p_num=p_num, t_num=t_num, model_kind=:PNJL)
    
    θ = [T_fm]
    seed = get_seed(seed_strategy, θ, mode)
    x0 = Float64.(seed)
    
    residual_fn! = build_residual!(mode, params)
    postprocess_fn = x_sol -> begin
        x_state = SVector{5}(Tuple(x_sol[1:5]))
        mu_vec = SVector{3}(x_sol[6], x_sol[7], x_sol[8])
        return _postprocess_payload(:PNJL, x_state, mu_vec, T_fm;
            p_num=p_num,
            t_num=t_num,
            xi=xi,
        )
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
@inline function _solve_multi_collect(
    mode::ConstraintMode,
    θ::AbstractVector,
    seed_strategy::MultiSeed,
    run_with_seed::Function,
)
    seeds = get_all_seeds(seed_strategy, θ, mode)

    results = SolverResult[]
    for seed in seeds
        try
            push!(results, run_with_seed(seed))
        catch e
            @warn "Solve failed with seed" seed exception=e
        end
    end

    isempty(results) && error("All seeds failed (exceptions) in solve_multi")

    converged = filter(r -> r.converged, results)
    isempty(converged) && error("All seeds failed to converge to a physical solution")

    return seed_strategy.selector(converged)
end

@inline function _is_physically_preferred_result(r::SolverResult)
    return _default_is_physical_solution(r.x_state, r.masses) &&
           _all_finite_thermo(r.omega, r.pressure, r.rho_norm, r.entropy, r.energy)
end

function _physical_first_omega_selector(results::AbstractVector{SolverResult})
    converged = filter(r -> r.converged, results)
    isempty(converged) && return first(results)

    physical = filter(_is_physically_preferred_result, converged)
    pool = isempty(physical) ? converged : physical

    return argmin(r -> (r.omega, r.residual_norm), pool)
end

function solve_multi(mode::FixedMu, T_fm::Real, μ_fm::Real;
                     seed_strategy::MultiSeed=MultiSeed(),
                     nlsolve_method::Symbol=:newton,
                     kwargs...)
    θ = [T_fm, μ_fm]
    return _solve_multi_collect(mode, θ, seed_strategy, seed -> begin
        solve(mode, T_fm, μ_fm;
            seed_strategy=DefaultSeed(seed, seed, :hadron),
            nlsolve_method=nlsolve_method,
            auto_multiseed_fallback=false,
            kwargs...)
    end)
end

function solve_multi(mode::FixedRho, T_fm::Real;
                     seed_strategy::MultiSeed=MultiSeed(),
                     nlsolve_method::Symbol=:newton,
                     kwargs...)
    θ = [T_fm]
    return _solve_multi_collect(mode, θ, seed_strategy, seed -> begin
        solve(mode, T_fm;
            seed_strategy=DefaultSeed(seed, seed, :hadron),
            nlsolve_method=nlsolve_method,
            kwargs...)
    end)
end

function solve_multi(mode::FixedAsymmetricRho, T_fm::Real;
                     seed_strategy::MultiSeed=MultiSeed(),
                     nlsolve_method::Symbol=:newton,
                     kwargs...)
    effective_seed_strategy = seed_strategy
    if seed_strategy.selector === default_omega_selector
        effective_seed_strategy = MultiSeed(seed_strategy.candidates, _physical_first_omega_selector)
    end

    θ = [T_fm]
    return _solve_multi_collect(mode, θ, effective_seed_strategy, seed -> begin
        solve(mode, T_fm;
            seed_strategy=DefaultSeed(seed, seed, :hadron),
            nlsolve_method=nlsolve_method,
            kwargs...)
    end)
end

export solve_multi

# ============================================================================
# Weighted-Block fallback（用于 FixedAsymmetricRho 失败点救援）
# ============================================================================

@inline function _weighted_stage_schedule()
    return (
        (w6=0.10, w7=0.10, w8=1.0, method=:trust_region, iterations=600),
        (w6=0.50, w7=0.50, w8=1.0, method=:trust_region, iterations=1000),
        (w6=1.00, w7=1.00, w8=1.0, method=:trust_region, iterations=1400),
        (w6=1.00, w7=1.00, w8=1.0, method=:newton,       iterations=1000),
    )
end

function _run_weighted_stage(raw_residual!::Function, x0::Vector{Float64};
    w6::Float64,
    w7::Float64,
    w8::Float64,
    method::Symbol,
    iterations::Int,
)
    weights = Float64[1.0, 1.0, 1.0, 1.0, 1.0, w6, w7, w8]

    weighted_residual! = (F, x) -> begin
        rawF = similar(F)
        raw_residual!(rawF, x)
        @inbounds for i in 1:8
            F[i] = weights[i] * rawF[i]
        end
        return nothing
    end

    res = nlsolve(weighted_residual!, x0;
        autodiff=:forward,
        method=method,
        xtol=1e-9,
        ftol=1e-9,
        iterations=iterations,
    )

    x = Vector{Float64}(res.zero)
    rawF = zeros(8)
    raw_residual!(rawF, x)
    raw_norm = sqrt(sum(abs2, rawF))

    return (x=x, raw_norm=raw_norm)
end

"""
    solve_weighted_block_fallback(mode::FixedAsymmetricRho, T_fm; initial_seed, kwargs...) -> Union{SolverResult, Nothing}

仅用于 Hybrid 失败点的兜底：
1) 对约束分量做分阶段加权求解（w6/w7 由小到大）；
2) 用阶段最优点作为 seed，再走标准 `solve` 做严格判定。
"""
function solve_weighted_block_fallback(mode::FixedAsymmetricRho, T_fm::Real;
    initial_seed::AbstractVector{<:Real},
    max_seed_candidates::Int=3,
    xi::Real=0.0,
    p_num::Int=DEFAULT_MOMENTUM_COUNT,
    t_num::Int=DEFAULT_THETA_COUNT,
    model_kind::Symbol=:PNJL,
    residual_norm_max::Real=1e-6,
    nlsolve_kwargs...)

    base_seed = if length(initial_seed) >= 8
        Float64.(initial_seed[1:8])
    else
        extend_seed(Float64.(initial_seed), mode)
    end

    seed_candidates = Vector{Vector{Float64}}()
    push!(seed_candidates, copy(base_seed))
    extra = max(max_seed_candidates - 1, 0)
    for seed in Iterators.take(get_all_seeds(MultiSeed(), [T_fm], mode), extra)
        s8 = length(seed) >= 8 ? Float64.(seed[1:8]) : extend_seed(Float64.(seed), mode)
        push!(seed_candidates, s8)
    end

    uniq = Dict{String, Vector{Float64}}()
    for s in seed_candidates
        key = join(round.(s; digits=6), ",")
        if !haskey(uniq, key)
            uniq[key] = s
        end
    end
    seed_candidates = collect(values(uniq))

    thermal_nodes = cached_nodes(p_num, t_num)
    params = GapParams(Float64(T_fm), thermal_nodes, Float64(xi);
        p_num=p_num, t_num=t_num, model_kind=model_kind)
    raw_residual! = build_residual!(mode, params)

    best_x = copy(base_seed)
    best_raw = Inf

    for seed in seed_candidates
        x = copy(seed)
        for cfg in _weighted_stage_schedule()
            stage = try
                _run_weighted_stage(raw_residual!, x;
                    w6=cfg.w6,
                    w7=cfg.w7,
                    w8=cfg.w8,
                    method=cfg.method,
                    iterations=cfg.iterations,
                )
            catch
                nothing
            end
            stage === nothing && continue

            x = stage.x
            if isfinite(stage.raw_norm) && stage.raw_norm < best_raw
                best_raw = stage.raw_norm
                best_x = copy(stage.x)
            end
            if isfinite(best_raw) && best_raw <= 1e-8
                early_seed = DefaultSeed(best_x, best_x, :hadron)
                early_result = try
                    solve(mode, T_fm;
                        xi=xi,
                        seed_strategy=early_seed,
                        p_num=p_num,
                        t_num=t_num,
                        model_kind=model_kind,
                        nlsolve_method=:trust_region,
                        trust_region_fallback=true,
                        residual_norm_max=residual_norm_max,
                        nlsolve_kwargs...)
                catch
                    nothing
                end
                if early_result !== nothing && early_result.converged
                    return early_result
                end
            end
            if best_raw <= 1e-8
                break
            end
        end
        if best_raw <= 1e-8
            break
        end
    end

    final_seed = DefaultSeed(best_x, best_x, :hadron)
    result = try
        solve(mode, T_fm;
            xi=xi,
            seed_strategy=final_seed,
            p_num=p_num,
            t_num=t_num,
            model_kind=model_kind,
            nlsolve_method=:trust_region,
            trust_region_fallback=true,
            residual_norm_max=residual_norm_max,
            nlsolve_kwargs...)
    catch
        nothing
    end

    return result
end

# ============================================================================
# ImplicitDifferentiation.jl 集成
# ============================================================================

# 全局配置存储
const IMPLICIT_CONFIG = Ref{NamedTuple}((
    xi = 0.0,
    p_num = DEFAULT_MOMENTUM_COUNT,
    t_num = DEFAULT_THETA_COUNT,
    model_kind=:PNJL,
    thermal_nodes = cached_nodes(DEFAULT_MOMENTUM_COUNT, DEFAULT_THETA_COUNT),
))

function set_implicit_config(; xi::Real=0.0,
                             p_num::Int=DEFAULT_MOMENTUM_COUNT,
                             t_num::Int=DEFAULT_THETA_COUNT,
                             model_kind::Symbol=:PNJL)
    thermal_nodes = cached_nodes(p_num, t_num)
    IMPLICIT_CONFIG[] = (
        xi=Float64(xi),
        p_num=p_num,
        t_num=t_num,
        model_kind=model_kind,
        thermal_nodes=thermal_nodes,
    )
    return nothing
end

@inline function _forward_solve_mu_with_config(θ::AbstractVector, config::NamedTuple)
    T_fm = Float64(θ[1])
    μ_fm = Float64(θ[2])

    mu_vec = SVector{3, Float64}(μ_fm, μ_fm, μ_fm)
    params = GapParams(T_fm, config.thermal_nodes, config.xi;
        p_num=config.p_num,
        t_num=config.t_num,
        model_kind=config.model_kind,
    )

    seed = get_seed(DefaultSeed(), θ, FixedMu())
    residual_fn! = build_residual!(FixedMu(), mu_vec, params)
    res = nlsolve(residual_fn!, seed; autodiff=:forward, method=:newton, xtol=1e-9, ftol=1e-9)

    return (res.zero, nothing)
end

@inline function _conditions_mu_with_config(θ::AbstractVector, x::AbstractVector, z, config::NamedTuple)
    _ = z
    T_fm = θ[1]
    μ_fm = θ[2]

    mu_vec = SVector{3}(μ_fm, μ_fm, μ_fm)
    params = GapParams(T_fm, config.thermal_nodes, config.xi;
        p_num=config.p_num,
        t_num=config.t_num,
        model_kind=config.model_kind,
    )
    x_state = SVector{5}(x[1], x[2], x[3], x[4], x[5])

    return Vector(gap_conditions(x_state, mu_vec, params))
end

@inline function _build_fixedmu_implicit_adapters(config::NamedTuple)
    return (
        forward_solve=(θ -> _forward_solve_mu_with_config(θ, config)),
        conditions=((θ, x, z) -> _conditions_mu_with_config(θ, x, z, config)),
    )
end

"""
前向求解函数（ImplicitDifferentiation.jl 接口）

θ = [T, μ] -> x = [φ_u, φ_d, φ_s, Φ, Φ̄]
"""
function forward_solve_mu(θ::AbstractVector)
    config = IMPLICIT_CONFIG[]
    return _forward_solve_mu_with_config(θ, config)
end

"""
条件函数（ImplicitDifferentiation.jl 接口）
"""
function conditions_mu(θ::AbstractVector, x::AbstractVector, z)
    config = IMPLICIT_CONFIG[]
    return _conditions_mu_with_config(θ, x, z, config)
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
function create_implicit_solver(; xi::Real=0.0,
                                p_num::Int=DEFAULT_MOMENTUM_COUNT,
                                t_num::Int=DEFAULT_THETA_COUNT,
                                model_kind::Symbol=:PNJL)
    thermal_nodes = cached_nodes(p_num, t_num)
    local_config = (
        xi=Float64(xi),
        p_num=p_num,
        t_num=t_num,
        model_kind=model_kind,
        thermal_nodes=thermal_nodes,
    )
    adapters = _build_fixedmu_implicit_adapters(local_config)

    # Backward-compatible side effect for callers relying on set_implicit_config/global path.
    IMPLICIT_CONFIG[] = local_config

    return ImplicitFunction(
        adapters.forward_solve,
        adapters.conditions;
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
                                model_kind::Symbol=:PNJL,
                                p_num::Int=DEFAULT_MOMENTUM_COUNT,
                                t_num::Int=DEFAULT_THETA_COUNT)
    solver = create_implicit_solver(
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        model_kind=model_kind,
    )
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

"""
    solve_with_root_diagnostics(mode::FixedMu, T_fm, μ_fm; kwargs...) -> NamedTuple

与 `solve(mode, ...)` 行为一致，但额外返回 root diagnostics 结构，
用于迁移期诊断兼容。
"""
function solve_with_root_diagnostics(::FixedMu, T_fm::Real, μ_fm::Real;
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
    thermal_nodes = cached_nodes(p_num, t_num)
    params = GapParams(Float64(T_fm), thermal_nodes, Float64(xi); p_num=p_num, t_num=t_num, model_kind=:PNJL)
    mu_vec = SVector{3}(μ_fm, μ_fm, μ_fm)

    θ = [T_fm, μ_fm]
    seed = get_seed(seed_strategy, θ, mode)
    x0 = Float64.(seed)

    residual_fn! = build_residual!(mode, mu_vec, params)
    postprocess_fn = x_sol -> begin
        x_state = SVector{5}(Tuple(x_sol))
        return _postprocess_payload(:PNJL, x_state, mu_vec, T_fm;
            p_num=p_num,
            t_num=t_num,
            xi=xi,
        )
    end

    cache = Dict{Symbol,Tuple{Any,Any}}()
    solve_once = function (method::Symbol, seedv::Vector{Float64})
        res = nlsolve(residual_fn!, seedv; autodiff=:forward, method=method, xtol=1e-9, ftol=1e-9, nlsolve_kwargs...)
        local cand
        try
            cand = _postprocess_candidate(postprocess_fn, physicality_check, res.zero)
        catch
            cand = (phys=false,
                    x_sol=Vector{Float64}(res.zero),
                    x_state=SVector{5, Float64}(fill(NaN, 5)),
                    mu_vec=SVector{3, Float64}(fill(NaN, 3)),
                    omega=NaN,
                    pressure=NaN,
                    rho_norm=NaN,
                    entropy=NaN,
                    energy=NaN,
                    masses=SVector{3, Float64}(fill(NaN, 3)))
        end

        cache[method] = (res, cand)
        return (
            x=Vector{Float64}(res.zero),
            converged=Bool(res.f_converged) && Bool(cand.phys),
            residual_norm=Float64(res.residual_norm),
            score=isfinite(cand.omega) ? Float64(cand.omega) : NaN,
        )
    end

    policy = RootPolicy(
        primary_method=nlsolve_method,
        fallback_method=fallback_method,
        use_fallback=trust_region_fallback,
        use_multiseed=false,
        residual_norm_max=Float64(residual_norm_max),
        require_converged=true,
        diagnostics_level=:basic,
    )

    solved = solve_root_with_policy(solve_once, x0; policy=policy)
    selected_method = solved.diagnostics.attempts[solved.diagnostics.selected_attempt].method
    if !haskey(cache, selected_method)
        _ = solve_once(selected_method, x0)
    end
    res, cand = cache[selected_method]

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
        attempts = map(solved.diagnostics.attempts) do a
            (
                method=a.method,
                seed_source=a.seed_source,
                converged=a.converged,
                residual_norm=a.residual_norm,
                quality_tag=a.quality_tag,
                score=a.score,
            )
        end
        return (
            result=single,
            root_diagnostics=(
                selected_method=selected_method,
                selected_attempt=solved.diagnostics.selected_attempt,
                attempts=attempts,
            ),
        )
    end

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
    attempts = map(solved.diagnostics.attempts) do a
        (
            method=a.method,
            seed_source=a.seed_source,
            converged=a.converged,
            residual_norm=a.residual_norm,
            quality_tag=a.quality_tag,
            score=a.score,
        )
    end
    return (
        result=multi,
        root_diagnostics=(
            selected_method=:multiseed,
            selected_attempt=solved.diagnostics.selected_attempt,
            attempts=attempts,
        ),
    )
end

@inline function _diagnostics_payload_from_result(result::SolverResult; selected_method::Symbol=:unknown, seed_source::Symbol=:seed)
    attempt = (
        method=selected_method,
        seed_source=seed_source,
        converged=Bool(result.converged),
        residual_norm=Float64(result.residual_norm),
        quality_tag=Bool(result.converged) ? :good : :degraded,
        score=isfinite(result.omega) ? Float64(result.omega) : NaN,
    )
    return (
        result=result,
        root_diagnostics=(
            selected_method=selected_method,
            selected_attempt=1,
            attempts=[attempt],
        ),
    )
end

"""
    solve_with_root_diagnostics(mode::FixedRho, T_fm; kwargs...) -> NamedTuple

固定密度模式的兼容诊断入口。
"""
function solve_with_root_diagnostics(mode::FixedRho, T_fm::Real; kwargs...)
    result = solve(mode, T_fm; kwargs...)
    return _diagnostics_payload_from_result(result; selected_method=:legacy_fallback, seed_source=:seed)
end

"""
    solve_with_root_diagnostics(mode::FixedAsymmetricRho, T_fm; kwargs...) -> NamedTuple

固定非对称密度模式的兼容诊断入口。
"""
function solve_with_root_diagnostics(mode::FixedAsymmetricRho, T_fm::Real; kwargs...)
    result = solve(mode, T_fm; kwargs...)
    return _diagnostics_payload_from_result(result; selected_method=:legacy_fallback, seed_source=:seed)
end

"""
    solve_with_root_diagnostics(mode::FixedEntropy, T_fm; kwargs...) -> NamedTuple

固定熵密度模式的兼容诊断入口。
"""
function solve_with_root_diagnostics(mode::FixedEntropy, T_fm::Real; kwargs...)
    result = solve(mode, T_fm; kwargs...)
    return _diagnostics_payload_from_result(result; selected_method=:legacy_fallback, seed_source=:seed)
end

"""
    solve_with_root_diagnostics(mode::FixedSigma, T_fm; kwargs...) -> NamedTuple

固定比熵模式的兼容诊断入口。
"""
function solve_with_root_diagnostics(mode::FixedSigma, T_fm::Real; kwargs...)
    result = solve(mode, T_fm; kwargs...)
    return _diagnostics_payload_from_result(result; selected_method=:legacy_fallback, seed_source=:seed)
end

end # module ImplicitSolver
