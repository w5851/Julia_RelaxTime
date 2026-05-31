if !isdefined(Main, :EquilibriumFacade)
    @eval module EquilibriumFacade

"""
    EquilibriumFacade

薄封装：统一“平衡态求解 (solve_gap) + 状态向量展开 + 有效质量计算”的入口。

目的：
- workflow/scan 等高层不再各写一套 `thermo_backend/solver_backend` 的分支
- 保持 legacy/models 两后端语义一致

注意：这里的 `thermo_backend` 指派生量/残差评估所用热力学后端（:legacy|:models）；
`solver_backend` 仅在 models 模型支持时控制求解器实现（:legacy|:models）。
"""

using StaticArrays
using LinearAlgebra: norm
import Main: Models

const calculate_mass_vec = Models.calculate_mass_vec
const create_model = Models.create_model
const solve_gap = Models.solve_gap
const state_vector = Models.state_vector
const normalize_mu_vec = Models.normalize_mu_vec
const gap_residual = Models.gap_residual
const PNJLModel = Models.PNJLModel
const FixedMu = Models.FixedMu
const solve = Models.solve
const MeanFieldState = Models.MeanFieldState
const NLsolveGapSolver = Models.NLsolveGapSolver

export pnjl_model_kind
export solve_equilibrium_backend

const _MODEL_CACHE = Dict{Symbol, Any}()

@inline function _build_fixedmu_physicality_check(
    model;
    enforce_mass_floor::Bool,
    phi_max_positive::Float64,
    phi_min_negative::Float64,
)
    if !enforce_mass_floor && !isfinite(phi_max_positive) && !isfinite(phi_min_negative)
        return ((_, _) -> true)
    end

    m0_vec = calculate_mass_vec(model, SVector{3}(0.0, 0.0, 0.0))
    mass_floor = (
        u=Float64(m0_vec[1]),
        d=Float64(m0_vec[2]),
        s=Float64(m0_vec[3]),
    )

    return function (x_state, masses)
        if length(x_state) >= 3
            phi_u = Float64(x_state[1])
            phi_d = Float64(x_state[2])
            phi_s = Float64(x_state[3])
            if isfinite(phi_max_positive)
                (phi_u <= phi_max_positive && phi_d <= phi_max_positive && phi_s <= phi_max_positive) || return false
            end
            if isfinite(phi_min_negative)
                (phi_u >= phi_min_negative && phi_d >= phi_min_negative && phi_s >= phi_min_negative) || return false
            end
        end

        if enforce_mass_floor && length(masses) >= 3
            mu = Float64(masses[1])
            md = Float64(masses[2])
            ms = Float64(masses[3])
            (mu >= mass_floor.u && md >= mass_floor.d && ms >= mass_floor.s) || return false
        end

        return true
    end
end

@inline function _get_model(model_kind::Symbol)
    return get!(_MODEL_CACHE, model_kind) do
        create_model(model_kind)
    end
end

@inline function pnjl_model_kind(thermo_backend::Symbol)::Symbol
    return :PNJL
end

"""solve_equilibrium_backend(T_fm, mu_fm; ...) -> NamedTuple

返回值字段（稳定给 workflow 使用）：
- `converged`: Bool（对 models/legacy 都要求 solve_gap 成功，否则直接抛错）
- `x_state`: SVector{5}
- `mu_vec`: SVector{3}
- `masses`: SVector{3}
- `iterations`, `residual_norm`: 预留字段（目前为 missing）

关键词参数：
- `solver_backend`: :legacy|:models
- `seed_state`: solver_backend=:models 时的 5D 初值
- `solver_kwargs`: 透传给 legacy 求解器（例如 iterations 等）
- `models_solver`, `models_residual_norm_max`: models solver 配置
- `model`: 可选注入 models model（默认取 PNJL 缓存模型）
"""
function solve_equilibrium_backend(
    T_fm::Real,
    mu_fm::Real;
    xi::Real=0.0,
    solver_backend::Symbol=:auto,
    p_num::Int=24,
    t_num::Int=8,
    seed_state=nothing,
    solver_kwargs::NamedTuple=(;),
    models_solver=nothing,
    models_residual_norm_max::Real=1e-4,
    fixedmu_seed_strategy=nothing,
    fixedmu_evaluate_all_attempts::Bool=true,
    enforce_fixedmu_mass_floor::Bool=true,
    fixedmu_phi_max_positive::Real=0.1,
    fixedmu_phi_min_negative::Real=-Inf,
    model=nothing,
)
    kind = :PNJL

    m = model === nothing ? _get_model(kind) : model
    solved_mu_vec = nothing
    solved_iterations = missing
    solved_residual_norm = missing
    solved_converged = true

    effective_solver_backend = if solver_backend === :auto
        :models
    else
        solver_backend
    end

    st = if effective_solver_backend === :legacy
        # Legacy solver backend still accepts legacy-style solver kwargs.
        solve_gap(m, T_fm, mu_fm;
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            solver_kwargs...,
        )
    elseif effective_solver_backend === :models
        # In long-lived include-driven sessions (integration full profile), Models
        # may be re-included by different test files. A cached model created before a
        # re-include can fail strict `isa PNJLModel` even when semantically
        # equivalent. Rebuild from current Models module in that case.
        if !(m isa PNJLModel)
            m = create_model(kind)
        end
        if models_solver === nothing
            fixedmu_physicality_check = _build_fixedmu_physicality_check(
                m;
                enforce_mass_floor=Bool(enforce_fixedmu_mass_floor),
                phi_max_positive=Float64(fixedmu_phi_max_positive),
                phi_min_negative=Float64(fixedmu_phi_min_negative),
            )
            fixed_mode = FixedMu()
            if fixedmu_seed_strategy !== nothing
                solved = solve(m, fixed_mode, T_fm, mu_fm;
                    seed_strategy=fixedmu_seed_strategy,
                    evaluate_all_attempts=Bool(fixedmu_evaluate_all_attempts),
                    xi=xi,
                    p_num=p_num,
                    t_num=t_num,
                    residual_norm_max=models_residual_norm_max,
                    physicality_check=fixedmu_physicality_check,
                )
            elseif seed_state === nothing
                solved = solve(m, fixed_mode, T_fm, mu_fm;
                    xi=xi,
                    p_num=p_num,
                    t_num=t_num,
                    residual_norm_max=models_residual_norm_max,
                    physicality_check=fixedmu_physicality_check,
                )
            else
                solved = solve(m, fixed_mode, T_fm, mu_fm;
                    seed_guess=seed_state,
                    continuity_seed=true,
                    xi=xi,
                    p_num=p_num,
                    t_num=t_num,
                    residual_norm_max=models_residual_norm_max,
                    physicality_check=fixedmu_physicality_check,
                )
            end
            solved_converged = Bool(solved.converged)
            fixedmu_solution_ok = solved_converged && fixedmu_physicality_check(solved.x_state, solved.masses)
            if fixedmu_solution_ok
                st = MeanFieldState(solved.x_state)
                solved_mu_vec = solved.mu_vec
                solved_iterations = solved.iterations
                solved_residual_norm = solved.residual_norm
                st
            else
                fallback_solver = NLsolveGapSolver(method=:trust_region, jacobian=:forward)
                solve_gap(m, T_fm, mu_fm;
                    solver_backend=:models,
                    solver=fallback_solver,
                    initial_guess=seed_state,
                    residual_norm_max=models_residual_norm_max,
                    xi=xi,
                    p_num=p_num,
                    t_num=t_num,
                )
            end
        else
            solve_gap(m, T_fm, mu_fm;
                solver_backend=:models,
                solver=models_solver,
                initial_guess=seed_state,
                residual_norm_max=models_residual_norm_max,
                xi=xi,
                p_num=p_num,
                t_num=t_num,
            )
        end
    else
        throw(ArgumentError("unknown solver_backend=$solver_backend (expected :auto, :legacy or :models)"))
    end

    x_state = state_vector(st)
    mu_vec = solved_mu_vec === nothing ? normalize_mu_vec(mu_fm) : normalize_mu_vec(solved_mu_vec)

    masses = calculate_mass_vec(m, SVector{3}(Tuple(x_state[1:3])))

    return (
        converged=solved_converged,
        x_state=SVector{5}(Tuple(x_state)),
        mu_vec=SVector{3}(Tuple(mu_vec)),
        masses=masses,
        iterations=solved_iterations,
        residual_norm=solved_residual_norm,
    )
end

"""solve_equilibrium_backend(T_fm, mu_vec; ...) -> NamedTuple

Flavor-chemical-potential overload used by workflow/path layers when
`μ_u, μ_d, μ_s` must differ explicitly.
"""
function solve_equilibrium_backend(
    T_fm::Real,
    mu_vec::AbstractVector{<:Real};
    xi::Real=0.0,
    solver_backend::Symbol=:auto,
    p_num::Int=24,
    t_num::Int=8,
    seed_state=nothing,
    solver_kwargs::NamedTuple=(;),
    models_solver=nothing,
    models_residual_norm_max::Real=1e-4,
    model=nothing,
)
    kind = :PNJL
    m = model === nothing ? _get_model(kind) : model
    μ = normalize_mu_vec(mu_vec)

    effective_solver_backend = if solver_backend === :auto
        :models
    else
        solver_backend
    end

    st = if effective_solver_backend === :legacy
        solve_gap(m, T_fm, μ;
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            solver_kwargs...,
        )
    elseif effective_solver_backend === :models
        if models_solver === nothing
            solve_gap(m, T_fm, μ;
                solver_backend=:models,
                initial_guess=seed_state,
                residual_norm_max=models_residual_norm_max,
                xi=xi,
                p_num=p_num,
                t_num=t_num,
            )
        else
            solve_gap(m, T_fm, μ;
                solver_backend=:models,
                solver=models_solver,
                initial_guess=seed_state,
                residual_norm_max=models_residual_norm_max,
                xi=xi,
                p_num=p_num,
                t_num=t_num,
            )
        end
    else
        throw(ArgumentError("unknown solver_backend=$solver_backend (expected :auto, :legacy or :models)"))
    end

    x_state = state_vector(st)
    masses = calculate_mass_vec(m, SVector{3}(Tuple(x_state[1:3])))
    residual = gap_residual(m, x_state, T_fm, μ; xi=xi, p_num=p_num, t_num=t_num)
    residual_norm = norm(residual)

    return (
        converged=isfinite(residual_norm) && residual_norm <= Float64(models_residual_norm_max),
        x_state=SVector{5}(Tuple(x_state)),
        mu_vec=SVector{3}(Tuple(μ)),
        masses=masses,
        iterations=-1,
        residual_norm=Float64(residual_norm),
    )
end

end # module EquilibriumFacade
end
