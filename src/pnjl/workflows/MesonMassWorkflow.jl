module MesonMassWorkflow

"""
串联：PNJL 平衡求解 → 介子质量/宽度（MesonMass）→ Mott 阈值与 gap（MottTransition）。

定位：放在 src/pnjl/workflows 下，是因为“完整流程”需要 PNJL.solve 得到 (m_u,m_s,Φ,Φbar)。
若调用方已自行提供 (quark_params, thermo_params)，则可直接调用 src/relaxtime/MesonMass.jl 与
src/relaxtime/MottTransition.jl，无需本工作流。

单位约定：
- 输入与输出均使用 fm⁻¹（T_fm, mu_fm, 质量/动量）
"""

include("../../Constants_PNJL.jl")
include("../PNJL.jl")
include("../../relaxtime/MesonMass.jl")
include("../../relaxtime/MottTransition.jl")

using .Constants_PNJL: ħc_MeV_fm
using .PNJL: solve, FixedMu
using .PNJL: HADRON_SEED_5, DEFAULT_MOMENTUM_COUNT, DEFAULT_THETA_COUNT
using .MesonMass: solve_meson_mass, default_meson_mass_guess
using .MottTransition: mott_threshold_mass, mott_gap, mott_threshold_masses, mott_gaps

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
    quark_params::NamedTuple,
    thermo_params::NamedTuple;
    k_norm::Float64,
    mass_kwargs::NamedTuple,
)
    # 基础阈值：用于构造更稳健的初值候选。
    thr = if _is_mixed_meson(meson)
        mott_threshold_masses(meson, quark_params).min
    else
        mott_threshold_mass(meson, quark_params)
    end

    guess = default_meson_mass_guess(meson, quark_params)
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
                solve_meson_mass(meson, quark_params, thermo_params;
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
    quark_params = (
        m=(u=Float64(masses[1]), d=Float64(masses[2]), s=Float64(masses[3])),
        μ=(u=Float64(mu_fm), d=Float64(mu_fm), s=Float64(mu_fm)),
    )
    thermo_params = (T=Float64(T_fm), Φ=Φ, Φbar=Φbar, ξ=Float64(xi))
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
    mesons::Tuple{Vararg{Symbol}}=DEFAULT_MESONS,
    k_norm::Real=0.0,
    p_num::Int=DEFAULT_MOMENTUM_COUNT,
    t_num::Int=DEFAULT_THETA_COUNT,
    seed_state=HADRON_SEED_5,
    solver_kwargs::NamedTuple=(;),
    mass_kwargs::NamedTuple=(;),
)
    seed_strategy = if seed_state isa AbstractVector
        s5 = Float64.(seed_state[1:5])
        PNJL.DefaultSeed(s5, s5, :hadron)
    else
        PNJL.DefaultSeed(phase_hint=:auto)
    end

    base = solve(FixedMu(), T_fm, mu_fm;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        seed_strategy=seed_strategy,
        solver_kwargs...,
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
            thr = mott_threshold_masses(meson, quark_params)
            gaps = isfinite(mass) ? mott_gaps(meson, mass, quark_params) : (uu=NaN, ss=NaN, min=NaN)
            meson_results[meson] = (mass=mass, gamma=gamma, converged=converged, residual=residual,
                                    threshold=thr, gaps=gaps)
        else
            thr = mott_threshold_mass(meson, quark_params)
            gapv = isfinite(mass) ? mott_gap(meson, mass, quark_params) : NaN
            meson_results[meson] = (mass=mass, gamma=gamma, converged=converged, residual=residual,
                                    threshold=thr, gap=gapv)
        end
    end

    return (equilibrium=base, quark_params=quark_params, thermo_params=thermo_params, meson_results=meson_results)
end

end # module
