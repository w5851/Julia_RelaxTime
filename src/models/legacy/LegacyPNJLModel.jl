"""LegacyPNJLModel

阶段 1：legacy 适配器（adapter）。

目的：把旧 PNJL core（src/pnjl/core/Thermodynamics.jl + solver）包装成实现 Models 接口的模型类型，
让高层工作流可以只依赖 `Models.omega/solve_gap/number_densities` 等统一入口。

注意：该适配器刻意复用 legacy 的 `calculate_omega` / `calculate_number_densities`，
以便与 legacy backend 做逐点对比。
"""

using StaticArrays

export LegacyPNJLModel

# Load legacy thermodynamics once at module load (avoid runtime include world-age issues)
const _LEGACY_THERMO_PATH = normpath(joinpath(@__DIR__, "..", "..", "pnjl", "core", "Thermodynamics.jl"))
if !isdefined(@__MODULE__, :LegacyThermodynamics)
    include(_LEGACY_THERMO_PATH)
    const LegacyThermodynamics = Thermodynamics
end

const _pnjl_cached_nodes = LegacyThermodynamics.Integrals.cached_nodes
const DEFAULT_MOMENTUM_COUNT = LegacyThermodynamics.Integrals.DEFAULT_MOMENTUM_COUNT
const DEFAULT_THETA_COUNT = LegacyThermodynamics.Integrals.DEFAULT_THETA_COUNT

const _PNJL_SOLVER_PATH = normpath(joinpath(@__DIR__, "..", "..", "pnjl", "solver", "Solver.jl"))

@inline function _ensure_pnjl_solver_loaded()
    if !(isdefined(Main, :Solver) && isdefined(Main.Solver, :solve) && isdefined(Main.Solver, :FixedMu))
        Base.include(Main, _PNJL_SOLVER_PATH)
    end
    return nothing
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
    return LegacyThermodynamics.calculate_mass_vec(φ)
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

    nodes = isnothing(thermal_nodes) ? _pnjl_cached_nodes(p_num, t_num) : thermal_nodes
    thermal_p_mesh, cosθ_mesh, thermal_coefficients = nodes

    φ = SVector{3, eltype(x5)}(x5[1], x5[2], x5[3])
    Φ, Φbar = x5[4], x5[5]

    chi = LegacyThermodynamics.calculate_chiral(φ)
    poly = LegacyThermodynamics.calculate_U(T_fm, Φ, Φbar)

    masses = LegacyThermodynamics.calculate_mass_vec(φ)
    vac = LegacyThermodynamics.Integrals.calculate_energy_sum(masses)
    therm = LegacyThermodynamics.Integrals.calculate_log_sum(
        masses,
        thermal_p_mesh,
        cosθ_mesh,
        thermal_coefficients,
        Φ,
        Φbar,
        μ,
        T_fm,
        xi,
    )

    ω = chi + poly + vac + therm
    return (chi=chi, poly=poly, vac=vac, therm=therm, masses=masses, omega=ω)
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

    nodes = isnothing(thermal_nodes) ? _pnjl_cached_nodes(p_num, t_num) : thermal_nodes
    return LegacyThermodynamics.calculate_number_densities(x5, μ, T_fm, nodes, xi)
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
    thermo_backend::Symbol=:legacy,
    kwargs...
)
    _ensure_pnjl_solver_loaded()

    μ = normalize_mu_vec(mu_vec)
    if !(μ[1] == μ[2] == μ[3])
        throw(ArgumentError("LegacyPNJLModel.solve_gap currently requires mu_u == mu_d == mu_s (FixedMu mode). Got mu_vec=$μ"))
    end

    mode = Base.invokelatest(Main.Solver.FixedMu)
    res = Base.invokelatest(Main.Solver.solve, mode, T_fm, μ[1];
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        thermo_backend=thermo_backend,
        kwargs...)

    res.converged || error("LegacyPNJLModel.solve_gap did not converge (residual_norm=$(res.residual_norm))")
    return MeanFieldState(res.x_state)
end
