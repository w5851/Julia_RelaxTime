module MesonMassWorkflow

"""
串联：PNJL 平衡求解 → 介子质量/宽度（MesonMass）→ Mott 阈值与 gap（MottTransition）。

定位：当前作为 models 侧 workflow 实体承载模块。
若调用方已自行提供 (quark_params, thermo_params)，则可直接调用 src/relaxtime/MesonMass.jl 与
src/relaxtime/MottTransition.jl，无需本工作流。

单位约定：
- 输入与输出均使用 fm⁻¹（T_fm, mu_fm, 质量/动量）
"""

# Include-once helper
const _INCLUDE_ONCE_PATH = normpath(joinpath(@__DIR__, "..", "..", "utils", "IncludeOnce.jl"))
if !isdefined(Main, :IncludeOnce)
    Base.include(Main, _INCLUDE_ONCE_PATH)
end
const IncludeOnce = Main.IncludeOnce

# Load all relaxtime submodules via RelaxTime.jl (ensures Constants_PNJL, MesonMass, MottTransition etc.)
const _RELAX_TIME_MODULE_PATH = normpath(joinpath(@__DIR__, "..", "..", "relaxtime", "RelaxTime.jl"))
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RELAX_TIME_MODULE_PATH)
end

# Unified equilibrium facade (solve_gap + state_vector + masses)
const _EQUILIBRIUM_FACADE_PATH = normpath(joinpath(@__DIR__, "..", "pnjl_physics", "core", "EquilibriumFacade.jl"))
IncludeOnce.include_once!(Main, :EquilibriumFacade, _EQUILIBRIUM_FACADE_PATH)

# Shared parameter structs (QuarkParams/ThermoParams)
const _PARAMETER_TYPES_PATH = normpath(joinpath(@__DIR__, "..", "..", "types", "ParameterTypes.jl"))
IncludeOnce.include_once!(Main, :ParameterTypes, _PARAMETER_TYPES_PATH)
using Main.ParameterTypes: QuarkParams, ThermoParams

const _WORKFLOW_PARAM_ADAPTERS_PATH = normpath(joinpath(@__DIR__, "WorkflowParamAdapters.jl"))
const WorkflowParamAdapters = IncludeOnce.include_once!(Main, :WorkflowParamAdapters, _WORKFLOW_PARAM_ADAPTERS_PATH)
using .WorkflowParamAdapters: normalize_quark_params, normalize_thermo_params

using Main.Constants_PNJL: ħc_MeV_fm
const HADRON_SEED_5 = Main.Models.HADRON_SEED_5
const DEFAULT_MOMENTUM_COUNT = Main.Models.PNJLCore.DEFAULT_MOMENTUM_COUNT
const DEFAULT_THETA_COUNT = Main.Models.PNJLCore.DEFAULT_THETA_COUNT
using Main.MesonMass: solve_meson_mass, default_meson_mass_guess
using Main.MottTransition: mott_threshold_mass, mott_gap, mott_threshold_masses, mott_gaps

export DEFAULT_MESONS
export solve_gap_and_meson_point
export build_equilibrium_params

const DEFAULT_MESONS = (
    :pi,
    :K,
    :eta,
    :eta_prime,
    :sigma_pi,
    :sigma_K,
    :sigma,
    :sigma_prime,
)

@inline function _unique_positive_candidates(vals::Vector{Float64})
    out = Float64[]
    for v in vals
        isfinite(v) || continue
        v > 0 || continue
        any(x -> abs(x - v) ≤ 1e-12, out) && continue
        push!(out, v)
    end
    return out
end

function _solve_meson_mass_with_retries(
    meson::Symbol,
    quark_params,
    thermo_params;
    k_norm::Float64,
    mass_kwargs::NamedTuple,
)
    qp = normalize_quark_params(quark_params)
    tp = normalize_thermo_params(thermo_params)

    # 基础阈值：用于构造更稳健的初值候选。
    thr = if _is_mixed_meson(meson)
        mott_threshold_masses(meson, qp).min
    else
        mott_threshold_mass(meson, qp)
    end

    guess = default_meson_mass_guess(meson, qp)
    mass_candidates = _unique_positive_candidates(Float64[
        guess,
        0.7 * guess,
        0.4 * guess,
        0.2 * guess,
        0.7 * thr,
        0.4 * thr,
        0.2 * thr,
        thr,
        0.7,
        1.0,
    ])

    gamma_candidates = Float64[0.0, 1e-6, 1e-4]

    best = nothing
    best_resid = Inf

    for m0 in mass_candidates
        for g0 in gamma_candidates
            res = try
                solve_meson_mass(meson, qp, tp;
                    k_norm=k_norm,
                    initial_mass=m0,
                    initial_gamma=g0,
                    mass_kwargs...,
                )
            catch
                nothing
            end
            res === nothing && continue

            resid = Float64(res.residual_norm)
            if isfinite(resid) && resid < best_resid
                best = res
                best_resid = resid
            end
            if Bool(res.converged)
                return res
            end
        end
    end

    return best
end

@inline function _is_mixed_meson(meson::Symbol)::Bool
    return meson === :eta || meson === :eta_prime || meson === :sigma || meson === :sigma_prime
end

"""将 PNJL 平衡求解结果转换成 (quark_params, thermo_params)。"""
function build_equilibrium_params(base, T_fm::Real, mu_fm::Real; xi::Real=0.0)
    Φ = Float64(base.x_state[4])
    Φbar = Float64(base.x_state[5])

    masses = base.masses
    quark_params = QuarkParams((
        m=(u=Float64(masses[1]), d=Float64(masses[2]), s=Float64(masses[3])),
        μ=(u=Float64(mu_fm), d=Float64(mu_fm), s=Float64(mu_fm)),
    ))
    thermo_params = ThermoParams((T=Float64(T_fm), Φ=Φ, Φbar=Φbar, ξ=Float64(xi)))
    return (quark_params=quark_params, thermo_params=thermo_params)
end

"""一次性完成：平衡求解 → 介子质量/宽度 → Mott 阈值与 gap。

返回 NamedTuple：
- equilibrium: PNJL.solve 的输出
- quark_params, thermo_params
- meson_results: Dict{Symbol,NamedTuple}，每个介子对应 solve_meson_mass 的输出与阈值/gap

关键词：
- xi: 各向异性参数 ξ
- p_num/t_num: 能隙/密度积分节点（传给 PNJL.solve）
- seed_state: 初值策略（默认 HADRON_SEED_5）
- solver_kwargs: 透传到 PNJL.solve
- mass_kwargs: 透传到 MesonMass.solve_meson_mass（例如 nlsolve 的参数）
"""
function solve_gap_and_meson_point(
    T_fm::Real,
    mu_fm::Real;
    xi::Real=0.0,
    solver_backend::Symbol=:legacy,
    mesons::Tuple{Vararg{Symbol}}=DEFAULT_MESONS,
    k_norm::Real=0.0,
    p_num::Int=DEFAULT_MOMENTUM_COUNT,
    t_num::Int=DEFAULT_THETA_COUNT,
    seed_state=HADRON_SEED_5,
    solver_kwargs::NamedTuple=(;),
    models_solver=nothing,
    models_residual_norm_max::Real=1e-4,
    mass_kwargs::NamedTuple=(;),
)
    seed_guess = if seed_state isa AbstractVector
        length(seed_state) >= 5 || throw(ArgumentError("seed_state must have length >= 5 (got $(length(seed_state)))"))
        Float64.(seed_state[1:5])
    else
        seed_state
    end

    base = Main.EquilibriumFacade.solve_equilibrium_backend(
        T_fm,
        mu_fm;
        xi=xi,
        solver_backend=solver_backend,
        p_num=p_num,
        t_num=t_num,
        seed_state=seed_guess,
        solver_kwargs=solver_kwargs,
        models_solver=models_solver,
        models_residual_norm_max=models_residual_norm_max,
    )

    params = build_equilibrium_params(base, T_fm, mu_fm; xi=xi)
    quark_params = params.quark_params
    thermo_params = params.thermo_params

    meson_results = Dict{Symbol,NamedTuple}()

    for meson in mesons
        res = _solve_meson_mass_with_retries(meson, quark_params, thermo_params;
            k_norm=Float64(k_norm),
            mass_kwargs=mass_kwargs,
        )

        mass = res === nothing ? NaN : Float64(res.mass)
        gamma = res === nothing ? NaN : Float64(res.gamma)
        converged = res === nothing ? false : Bool(res.converged)
        residual = res === nothing ? Inf : Float64(res.residual_norm)

        if _is_mixed_meson(meson)
            qp = normalize_quark_params(quark_params)
            thr = mott_threshold_masses(meson, qp)
            gaps = isfinite(mass) ? mott_gaps(meson, mass, qp) : (uu=NaN, ss=NaN, min=NaN)
            meson_results[meson] = (mass=mass, gamma=gamma, converged=converged, residual=residual,
                                    threshold=thr, gaps=gaps)
        else
            qp = normalize_quark_params(quark_params)
            thr = mott_threshold_mass(meson, qp)
            gapv = isfinite(mass) ? mott_gap(meson, mass, qp) : NaN
            meson_results[meson] = (mass=mass, gamma=gamma, converged=converged, residual=residual,
                                    threshold=thr, gap=gapv)
        end
    end

    return (equilibrium=base, quark_params=quark_params, thermo_params=thermo_params, meson_results=meson_results)
end

end # module
