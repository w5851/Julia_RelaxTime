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

export pnjl_model_kind
export solve_equilibrium_backend

const _MODEL_CACHE = Dict{Symbol, Any}()

@inline function _get_model(model_kind::Symbol)
    return get!(_MODEL_CACHE, model_kind) do
        Main.Models.create_model(model_kind)
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
    model=nothing,
)
    kind = :PNJL

    isdefined(Main, :Models) || error("Models not loaded; expected Main.Models")
    isdefined(Main.Models, :solve_gap) || error("Models.solve_gap is not defined")
    isdefined(Main.Models, :state_vector) || error("Models.state_vector is not defined")
    isdefined(Main.Models, :normalize_mu_vec) || error("Models.normalize_mu_vec is not defined")

    m = model === nothing ? _get_model(kind) : model

    effective_solver_backend = if solver_backend === :auto
        :models
    else
        solver_backend
    end

    st = if effective_solver_backend === :legacy
        # Legacy solver backend still accepts legacy-style solver kwargs.
        Main.Models.solve_gap(m, T_fm, mu_fm;
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            solver_kwargs...,
        )
    elseif effective_solver_backend === :models
        # In long-lived include-driven sessions (integration full profile), Main.Models
        # may be re-included by different test files. A cached model created before a
        # re-include can fail strict `isa Main.Models.PNJLModel` even when semantically
        # equivalent. Rebuild from current Main.Models in that case.
        if !(m isa Main.Models.PNJLModel)
            m = Main.Models.create_model(kind)
        end

        solver = models_solver === nothing ? Main.Models.NLsolveGapSolver(method=:trust_region, jacobian=:forward) : models_solver
        Main.Models.solve_gap(m, T_fm, mu_fm;
            solver_backend=:models,
            solver=solver,
            initial_guess=seed_state,
            residual_norm_max=models_residual_norm_max,
            xi=xi,
            p_num=p_num,
            t_num=t_num,
        )
    else
        throw(ArgumentError("unknown solver_backend=$solver_backend (expected :auto, :legacy or :models)"))
    end

    x_state = Main.Models.state_vector(st)
    mu_vec = Main.Models.normalize_mu_vec(mu_fm)

    masses = Main.Models.calculate_mass_vec(m, SVector{3}(Tuple(x_state[1:3])))

    return (
        converged=true,
        x_state=SVector{5}(Tuple(x_state)),
        mu_vec=SVector{3}(Tuple(mu_vec)),
        masses=masses,
        iterations=missing,
        residual_norm=missing,
    )
end

end # module EquilibriumFacade
end
