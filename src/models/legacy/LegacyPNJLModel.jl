"""LegacyPNJLModel

阶段 1：legacy 适配器（adapter）。

目的：把旧 PNJL core（src/pnjl/core/Thermodynamics.jl + solver）包装成实现 Models 接口的模型类型，
让高层工作流可以只依赖 `Models.omega/solve_gap/number_densities` 等统一入口。

注意：该适配器刻意复用 legacy 的 `calculate_omega` / `calculate_number_densities`，
以便与 legacy backend 做逐点对比。
"""

using StaticArrays

export LegacyPNJLModel

# Include-once helper
const _INCLUDE_ONCE_PATH = normpath(joinpath(@__DIR__, "..", "..", "utils", "IncludeOnce.jl"))
if !isdefined(Main, :IncludeOnce)
    Base.include(Main, _INCLUDE_ONCE_PATH)
end
const IncludeOnce = Main.IncludeOnce

const _THERMO_FACADE_PATH = normpath(joinpath(@__DIR__, "..", "core", "ThermoFacade.jl"))
const ThermoFacade = IncludeOnce.include_once!(Main, :ThermoFacade, _THERMO_FACADE_PATH)

const DEFAULT_MOMENTUM_COUNT = ThermoFacade.default_momentum_count()
const DEFAULT_THETA_COUNT = ThermoFacade.default_theta_count()

const _PNJL_ENTRY_PATH = normpath(joinpath(@__DIR__, "..", "pnjl", "PNJL.jl"))

@inline function _pnjl_solver_module()
    pnjl = IncludeOnce.include_once!(Main, :PNJL, _PNJL_ENTRY_PATH)
    if !(isdefined(pnjl, :solve) && isdefined(pnjl, :FixedMu))
        error("PNJL module loaded but required API (solve/FixedMu) is missing")
    end
    return pnjl
end

"""legacy PNJL 模型适配器（参数仍来自 legacy 常量/配置）。"""
struct LegacyPNJLModel <: AbstractPNJLModel
end

"""calculate_mass_vec(::LegacyPNJLModel, φ) -> SVector{3}

Legacy 适配器的有效质量向量；直接复用 legacy thermo 的质量公式。
需要支持 ForwardDiff.Dual（用于 gap_residual=∇Ω 的 AD）。
"""
@inline function calculate_mass_vec(::LegacyPNJLModel, φ::SVector{3, T}; kwargs...) where {T}
    _ = kwargs
    return ThermoFacade.calculate_mass_vec_backend(φ; model_kind=:PNJL)
end

# -----------------------------------------------------------------------------
# Unified entrypoints
# -----------------------------------------------------------------------------

"""返回 (chi, poly, vac, therm, masses, omega) 的 NamedTuple（legacy 公式）。"""
function omega_components(
    ::LegacyPNJLModel,
    x_state,
    T_fm,
    mu_vec;
    p_num::Int=DEFAULT_MOMENTUM_COUNT,
    t_num::Int=DEFAULT_THETA_COUNT,
    thermal_nodes=nothing,
    xi=0.0,
    kwargs...
)
    _ = kwargs

    μ = normalize_mu_vec(mu_vec)
    st = meanfield_state(x_state)
    x5 = state_vector(st)

    return ThermoFacade.calculate_omega_components_backend(
        x5,
        μ,
        T_fm;
        model_kind=:PNJL,
        p_num=p_num,
        t_num=t_num,
        thermal_nodes=thermal_nodes,
        xi=xi,
    )
end

"""Ω 标量值（legacy 公式）。"""
function omega(
    model::LegacyPNJLModel,
    x_state,
    T_fm,
    mu_vec;
    p_num::Int=DEFAULT_MOMENTUM_COUNT,
    t_num::Int=DEFAULT_THETA_COUNT,
    thermal_nodes=nothing,
    xi=0.0,
    kwargs...
)
    return omega_components(model, x_state, T_fm, mu_vec;
        p_num=p_num,
        t_num=t_num,
        thermal_nodes=thermal_nodes,
        xi=xi,
        kwargs...).omega
end

"""返回 (quark, antiquark) 数密度（legacy 路径）。"""
function number_densities(
    ::LegacyPNJLModel,
    x_state,
    T_fm,
    mu_vec;
    thermal_nodes=nothing,
    p_num::Int=DEFAULT_MOMENTUM_COUNT,
    t_num::Int=DEFAULT_THETA_COUNT,
    xi=0.0,
    kwargs...
)
    _ = kwargs

    μ = normalize_mu_vec(mu_vec)
    st = meanfield_state(x_state)
    x5 = state_vector(st)

    return ThermoFacade.calculate_number_densities_backend(
        x5,
        μ,
        T_fm;
        model_kind=:PNJL,
        p_num=p_num,
        t_num=t_num,
        thermal_nodes=thermal_nodes,
        xi=xi,
    )
end

"""solve_gap(::LegacyPNJLModel, T, mu_vec) -> MeanFieldState

说明：复用 legacy solver 的 FixedMu 模式；默认 thermo_backend=:legacy。
"""
function solve_gap(
    ::LegacyPNJLModel,
    T_fm,
    mu_vec;
    xi::Real=0.0,
    p_num::Int=DEFAULT_MOMENTUM_COUNT,
    t_num::Int=DEFAULT_THETA_COUNT,
    kwargs...
)
    pnjl = _pnjl_solver_module()

    μ = normalize_mu_vec(mu_vec)
    if !(μ[1] == μ[2] == μ[3])
        throw(ArgumentError("LegacyPNJLModel.solve_gap currently requires mu_u == mu_d == mu_s (FixedMu mode). Got mu_vec=$μ"))
    end

    legacy_kwargs = Dict{Symbol,Any}(pairs(kwargs))
    pop!(legacy_kwargs, :thermo_backend, nothing)

    mode = Base.invokelatest(pnjl.FixedMu)
    res = Base.invokelatest(pnjl.solve, mode, T_fm, μ[1];
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        legacy_kwargs...)

    res.converged || error("LegacyPNJLModel.solve_gap did not converge (residual_norm=$(res.residual_norm))")
    return MeanFieldState(res.x_state)
end
