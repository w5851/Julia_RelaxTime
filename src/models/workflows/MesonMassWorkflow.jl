module MesonMassWorkflow

"""
串联：PNJL 平衡求解 → 介子质量/宽度（MesonMass）→ Mott 阈值与 gap（MottTransition）。

定位：当前作为 models 侧 workflow 实体承载模块。
若调用方已自行提供 (quark_params, thermo_params)，则可直接调用 src/relaxtime/MesonMass.jl 与
src/relaxtime/MottTransition.jl，无需本工作流。

单位约定：
- 输入与输出均使用 fm⁻¹（T_fm, mu_fm, 质量/动量）
"""

# Load all relaxtime submodules via RelaxTime.jl (ensures Constants_PNJL, MesonMass, MottTransition etc.)
const _RELAX_TIME_MODULE_PATH = normpath(joinpath(@__DIR__, "..", "..", "relaxtime", "RelaxTime.jl"))
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RELAX_TIME_MODULE_PATH)
end

# Unified equilibrium facade (solve_gap + state_vector + masses)
const _EQUILIBRIUM_FACADE_PATH = normpath(joinpath(@__DIR__, "..", "pnjl_physics", "core", "EquilibriumFacade.jl"))
if !isdefined(Main, :EquilibriumFacade)
    Base.include(Main, _EQUILIBRIUM_FACADE_PATH)
end

# Shared parameter structs (QuarkParams/ThermoParams)
const _PARAMETER_TYPES_PATH = normpath(joinpath(@__DIR__, "..", "..", "types", "ParameterTypes.jl"))
if !isdefined(Main, :ParameterTypes)
    Base.include(Main, _PARAMETER_TYPES_PATH)
end
using Main.ParameterTypes: QuarkParams, ThermoParams

const _WORKFLOW_PARAM_ADAPTERS_PATH = normpath(joinpath(@__DIR__, "WorkflowParamAdapters.jl"))
if !isdefined(Main, :WorkflowParamAdapters)
    Base.include(Main, _WORKFLOW_PARAM_ADAPTERS_PATH)
end
const WorkflowParamAdapters = Main.WorkflowParamAdapters
using .WorkflowParamAdapters: normalize_quark_params, normalize_thermo_params

using Main.Constants_PNJL: ħc_MeV_fm
using Main.Models: RootPolicy, solve_root_with_policy
const HADRON_SEED_5 = Main.Models.HADRON_SEED_5
const DEFAULT_MOMENTUM_COUNT = Main.Models.default_momentum_count()
const DEFAULT_THETA_COUNT = Main.Models.default_theta_count()
using Main.MesonMass: solve_meson_mass, default_meson_mass_guess, ensure_quark_params_has_A
using Main.PolarizationAniso: polarization_with_width
using Main.EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings, mixing_matrix_elements
using Main.Constants_PNJL: G_fm2, K_fm5
using Main.MottTransition: mott_threshold_mass, mott_gap, mott_threshold_masses, mott_gaps

export DEFAULT_MESONS
export solve_gap_and_meson_point
export build_equilibrium_params

const _DEFAULT_MESON_ROOT_POLICY = (
    residual_norm_max=1e-6,
    require_converged=true,
    use_trust_region_fallback=true,
)

@inline function _as_vec2(x)::Union{Nothing,Vector{Float64}}
    x === nothing && return nothing
    if x isa AbstractVector && length(x) >= 2
        return Float64[x[1], x[2]]
    end
    return nothing
end

@inline function _mass_result_good(res, policy)::Bool
    res === nothing && return false
    isfinite(Float64(res.residual_norm)) || return false
    Float64(res.residual_norm) <= Float64(policy.residual_norm_max) || return false
    if Bool(policy.require_converged)
        return Bool(res.converged)
    end
    return true
end

@inline function _mixed_pair(meson::Symbol)
    if meson === :eta || meson === :eta_prime
        return (:eta, :eta_prime, :P)
    elseif meson === :sigma || meson === :sigma_prime
        return (:sigma, :sigma_prime, :S)
    end
    throw(ArgumentError("not a mixed meson: $meson"))
end

@inline function _label_sign_for_mixed(meson::Symbol)::Int
    if meson === :eta || meson === :sigma
        return -1
    elseif meson === :eta_prime || meson === :sigma_prime
        return 1
    end
    throw(ArgumentError("not a mixed meson: $meson"))
end

function _mixed_branch_scores(
    meson::Symbol,
    mass::Float64,
    gamma::Float64,
    quark_params,
    thermo_params,
    k_norm::Float64,
)
    (m_light, m_heavy, channel_sym) = _mixed_pair(meson)
    _ = (m_light, m_heavy)

    qp = normalize_quark_params(quark_params)
    tp = normalize_thermo_params(thermo_params)
    qpA = ensure_quark_params_has_A((
        m=(u=Float64(qp.m.u), d=Float64(qp.m.d), s=Float64(qp.m.s)),
        μ=(u=Float64(qp.μ.u), d=Float64(qp.μ.d), s=Float64(qp.μ.s)),
    ), (
        T=Float64(tp.T),
        Φ=Float64(tp.Φ),
        Φbar=Float64(tp.Φbar),
        ξ=Float64(tp.ξ),
    ))

    G_u = calculate_G_from_A(qpA.A.u, qpA.m.u)
    G_s = calculate_G_from_A(qpA.A.s, qpA.m.s)
    Kc = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

    μq = Float64(qp.μ.u)
    m_u = Float64(qp.m.u)
    m_s = Float64(qp.m.s)

    Πuu_re, Πuu_im = polarization_with_width(
        channel_sym,
        mass,
        gamma,
        k_norm,
        m_u,
        m_u,
        μq,
        μq,
        Float64(tp.T),
        Float64(tp.Φ),
        Float64(tp.Φbar),
        Float64(tp.ξ),
        qpA.A.u,
        qpA.A.u,
        0,
    )
    Πss_re, Πss_im = polarization_with_width(
        channel_sym,
        mass,
        gamma,
        k_norm,
        m_s,
        m_s,
        μq,
        μq,
        Float64(tp.T),
        Float64(tp.Φ),
        Float64(tp.Φbar),
        Float64(tp.ξ),
        qpA.A.s,
        qpA.A.s,
        2,
    )
    Πuu = ComplexF64(Πuu_re, Πuu_im)
    Πss = ComplexF64(Πss_re, Πss_im)
    elems = mixing_matrix_elements(Πuu, Πss, Kc, channel_sym)

    root = sqrt((elems.M00 - elems.M88)^2 + 4.0 * elems.M08^2)
    plus_v = (elems.M00 + elems.M88) + root
    minus_v = (elems.M00 + elems.M88) - root

    score_plus = hypot(real(plus_v), imag(plus_v))
    score_minus = hypot(real(minus_v), imag(minus_v))
    return score_plus, score_minus
end

function _annotate_mixed_branch_scores(
    meson_results::Dict{Symbol,NamedTuple},
    quark_params,
    thermo_params,
    k_norm::Float64,
)
    for meson in keys(meson_results)
        _is_mixed_meson(meson) || continue
        entry = meson_results[meson]
        if !isfinite(entry.mass) || !isfinite(entry.gamma)
            meson_results[meson] = merge(entry, (
                branch_score_plus=Inf,
                branch_score_minus=Inf,
                branch_score_selected=Inf,
                branch_sign=0,
            ))
            continue
        end
        splus, sminus = _mixed_branch_scores(
            meson,
            Float64(entry.mass),
            Float64(entry.gamma),
            quark_params,
            thermo_params,
            k_norm,
        )
        sign_label = _label_sign_for_mixed(meson)
        ssel = sign_label == 1 ? splus : sminus
        meson_results[meson] = merge(entry, (
            branch_score_plus=splus,
            branch_score_minus=sminus,
            branch_score_selected=ssel,
            branch_sign=sign_label,
        ))
    end
end

function _apply_identity_alignment!(meson_results::Dict{Symbol,NamedTuple})
    for pair in ((:eta, :eta_prime), (:sigma, :sigma_prime))
        m1, m2 = pair
        (haskey(meson_results, m1) && haskey(meson_results, m2)) || continue

        e1 = meson_results[m1]
        e2 = meson_results[m2]

        # Option A: label mapping (m1->-, m2->+)
        cost_A = Float64(e1.branch_score_minus) + Float64(e2.branch_score_plus)
        # Option B: swapped identity (m1->+, m2->-)
        cost_B = Float64(e1.branch_score_plus) + Float64(e2.branch_score_minus)

        if cost_B < cost_A
            meson_results[m1], meson_results[m2] = e2, e1
            e1 = meson_results[m1]
            e2 = meson_results[m2]
            meson_results[m1] = merge(e1, (
                branch_score_selected=Float64(e1.branch_score_minus),
                branch_sign=-1,
            ))
            meson_results[m2] = merge(e2, (
                branch_score_selected=Float64(e2.branch_score_plus),
                branch_sign=1,
            ))
        else
            meson_results[m1] = merge(e1, (
                branch_score_selected=Float64(e1.branch_score_minus),
                branch_sign=-1,
            ))
            meson_results[m2] = merge(e2, (
                branch_score_selected=Float64(e2.branch_score_plus),
                branch_sign=1,
            ))
        end
    end
end

@inline function _tracking_keys_for_mixed_pair(meson::Symbol)
    if meson === :eta || meson === :eta_prime
        return (:eta_minus, :eta_plus)
    elseif meson === :sigma || meson === :sigma_prime
        return (:sigma_minus, :sigma_plus)
    end
    throw(ArgumentError("not a mixed meson: $meson"))
end

@inline function _seed_distance(entry, seed::Vector{Float64})
    m = Float64(entry.mass)
    g = Float64(entry.gamma)
    if !(isfinite(m) && isfinite(g) && isfinite(seed[1]) && isfinite(seed[2]))
        return Inf
    end
    return hypot(m - seed[1], g - seed[2])
end

function _normalize_tracking_seed_state(mixed_seed_tracking_state)
    out = Dict{Symbol,Vector{Float64}}()
    mixed_seed_tracking_state isa AbstractDict || return out
    for (k, v) in mixed_seed_tracking_state
        k isa Symbol || continue
        vv = _as_vec2(v)
        vv === nothing && continue
        out[k] = vv
    end
    return out
end

function _apply_identity_tracking_label_output!(
    meson_results::Dict{Symbol,NamedTuple},
    tracking_seed_in::Dict{Symbol,Vector{Float64}},
)
    tracking_out = Dict{Symbol,Vector{Float64}}()

    for pair in ((:eta, :eta_prime), (:sigma, :sigma_prime))
        m1, m2 = pair
        (haskey(meson_results, m1) && haskey(meson_results, m2)) || continue

        e1 = meson_results[m1]
        e2 = meson_results[m2]
        minus_key, plus_key = _tracking_keys_for_mixed_pair(m1)

        prev_minus = get(tracking_seed_in, minus_key, nothing)
        prev_plus = get(tracking_seed_in, plus_key, nothing)

        continuity_A = Inf
        continuity_B = Inf
        if prev_minus !== nothing && prev_plus !== nothing
            continuity_A = _seed_distance(e1, prev_minus) + _seed_distance(e2, prev_plus)
            continuity_B = _seed_distance(e2, prev_minus) + _seed_distance(e1, prev_plus)
        end

        branch_A = Float64(e1.branch_score_minus) + Float64(e2.branch_score_plus)
        branch_B = Float64(e1.branch_score_plus) + Float64(e2.branch_score_minus)

        swap = if isfinite(continuity_A) && isfinite(continuity_B)
            continuity_B < continuity_A || (isapprox(continuity_B, continuity_A; atol=1e-12) && branch_B < branch_A)
        else
            branch_B < branch_A
        end

        minus_entry = swap ? e2 : e1
        plus_entry = swap ? e1 : e2

        meson_results[m1] = merge(minus_entry, (
            branch_score_selected=Float64(minus_entry.branch_score_minus),
            branch_sign=-1,
        ))
        meson_results[m2] = merge(plus_entry, (
            branch_score_selected=Float64(plus_entry.branch_score_plus),
            branch_sign=1,
        ))

        tracking_out[minus_key] = Float64[Float64(minus_entry.mass), Float64(minus_entry.gamma)]
        tracking_out[plus_key] = Float64[Float64(plus_entry.mass), Float64(plus_entry.gamma)]
    end

    return tracking_out
end

function _solve_meson_mass_with_policy(
    meson::Symbol,
    quark_params,
    thermo_params;
    k_norm::Float64,
    mass_kwargs::NamedTuple,
    policy::NamedTuple,
    experiment_mode::Symbol=:default,
    initial_seed::Union{Nothing,Vector{Float64}}=nothing,
)
    qp = normalize_quark_params(quark_params)
    tp = normalize_thermo_params(thermo_params)

    function run_once(method::Union{Nothing,Symbol}; m0::Float64, g0::Float64)
        extra_kwargs = method === nothing ? mass_kwargs : merge(mass_kwargs, (method=method,))
        return try
            solve_meson_mass(meson, qp, tp;
                k_norm=k_norm,
                initial_mass=m0,
                initial_gamma=g0,
                extra_kwargs...,
            )
        catch
            nothing
        end
    end

    if initial_seed !== nothing
        if experiment_mode === :fortran_track
            m0 = initial_seed[1] * 1.001
            g0 = max(initial_seed[2] * 1.001, 0.0)
            seed_res = run_once(nothing; m0=m0, g0=g0)
            if seed_res === nothing
                diag = (selected_method=:fortran_track, selected_attempt=0, attempts=NamedTuple[])
                return nothing, :bad, nothing, diag
            end
            attempt = (
                method=:fortran_track,
                seed_source=:seed,
                converged=Bool(seed_res.converged),
                residual_norm=Float64(seed_res.residual_norm),
                quality_tag=_mass_result_good(seed_res, policy) ? :good : :bad,
                score=NaN,
            )
            diag = (selected_method=:fortran_track, selected_attempt=1, attempts=[attempt])
            if _mass_result_good(seed_res, policy)
                return seed_res, :good, Float64[seed_res.mass, seed_res.gamma], diag
            end
            return seed_res, :bad, Float64[seed_res.mass, seed_res.gamma], diag
        end
    end

    solve_once_callback = function (method::Symbol, seed::Vector{Float64})
        picked = run_once(method; m0=seed[1], g0=max(seed[2], 0.0))
        if picked !== nothing
            qp_local = normalize_quark_params(quark_params)
            tp_local = normalize_thermo_params(thermo_params)
            local_score = try
                splus, sminus = _mixed_branch_scores(
                    meson,
                    Float64(picked.mass),
                    Float64(picked.gamma),
                    qp_local,
                    tp_local,
                    Float64(k_norm),
                )
                min(splus, sminus)
            catch
                NaN
            end
            return merge(picked, (score=local_score,))
        end
        return picked
    end

    guess = default_meson_mass_guess(meson, qp)
    seed0 = initial_seed === nothing ? Float64[guess, 0.0] : Float64[initial_seed[1], initial_seed[2]]

    root_policy = RootPolicy(
        primary_method=:newton,
        fallback_method=:trust_region,
        use_fallback=Bool(policy.use_trust_region_fallback),
        use_multiseed=true,
        residual_norm_max=Float64(policy.residual_norm_max),
        require_converged=Bool(policy.require_converged),
        diagnostics_level=:basic,
    )

    solved = solve_root_with_policy(solve_once_callback, seed0; policy=root_policy)
    root_quality = solved.quality_tag
    diag_attempts = map(solved.diagnostics.attempts) do a
        (
            method=a.method,
            seed_source=a.seed_source,
            converged=a.converged,
            residual_norm=a.residual_norm,
            quality_tag=a.quality_tag,
            score=a.score,
        )
    end
    root_diag = (
        selected_method=solved.diagnostics.attempts[solved.diagnostics.selected_attempt].method,
        selected_attempt=solved.diagnostics.selected_attempt,
        attempts=diag_attempts,
    )

    selected_method = solved.diagnostics.attempts[solved.diagnostics.selected_attempt].method
    selected_seed = Float64[solved.x[1], solved.x[2]]
    selected = run_once(selected_method; m0=selected_seed[1], g0=max(selected_seed[2], 0.0))

    if selected === nothing
        if !isfinite(solved.x[1]) || !isfinite(solved.x[2])
            return nothing, root_quality, nothing, root_diag
        end
        selected = (
            mass=Float64(solved.x[1]),
            gamma=Float64(solved.x[2]),
            converged=Bool(solved.converged),
            residual_norm=Float64(solved.residual_norm),
            solution=nothing,
        )
    end

    out_seed = (isfinite(selected.mass) && isfinite(selected.gamma)) ? Float64[selected.mass, selected.gamma] : nothing
    return selected, root_quality, out_seed, root_diag
end

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
    solver_backend::Symbol=:auto,
    mesons::Tuple{Vararg{Symbol}}=DEFAULT_MESONS,
    k_norm::Real=0.0,
    p_num::Int=DEFAULT_MOMENTUM_COUNT,
    t_num::Int=DEFAULT_THETA_COUNT,
    seed_state=HADRON_SEED_5,
    solver_kwargs::NamedTuple=(;),
    models_solver=nothing,
    models_residual_norm_max::Real=1e-4,
    mass_kwargs::NamedTuple=(;),
    meson_seed_state=nothing,
    meson_root_policy::NamedTuple=(;),
    meson_experiment_mode::Symbol=:default,
    mixed_branch_align::Symbol=:identity_track_label_output,
    mixed_seed_tracking_state=nothing,
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

    effective_policy = merge(_DEFAULT_MESON_ROOT_POLICY, meson_root_policy)

    meson_results = Dict{Symbol,NamedTuple}()
    meson_seed_out = Dict{Symbol,Vector{Float64}}()

    seed_in = Dict{Symbol,Vector{Float64}}()
    tracking_seed_in = _normalize_tracking_seed_state(mixed_seed_tracking_state)
    if meson_seed_state isa AbstractDict
        for (k, v) in meson_seed_state
            k isa Symbol || continue
            vv = _as_vec2(v)
            vv === nothing && continue
            seed_in[k] = vv
        end
    end

    if mixed_branch_align === :identity_track_label_output
        for meson in mesons
            _is_mixed_meson(meson) || continue
            minus_key, plus_key = _tracking_keys_for_mixed_pair(meson)
            key = _label_sign_for_mixed(meson) == -1 ? minus_key : plus_key
            if haskey(tracking_seed_in, key)
                seed_in[meson] = tracking_seed_in[key]
            end
        end
    end

    for meson in mesons
        initial_seed = haskey(seed_in, meson) ? seed_in[meson] : nothing
        res, root_quality, seed_vec, root_diagnostics = _solve_meson_mass_with_policy(
            meson,
            quark_params,
            thermo_params;
            k_norm=Float64(k_norm),
            mass_kwargs=mass_kwargs,
            policy=effective_policy,
            experiment_mode=meson_experiment_mode,
            initial_seed=initial_seed,
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
                                    root_quality=root_quality,
                                    root_diagnostics=root_diagnostics,
                                    threshold=thr, gaps=gaps)
        else
            qp = normalize_quark_params(quark_params)
            thr = mott_threshold_mass(meson, qp)
            gapv = isfinite(mass) ? mott_gap(meson, mass, qp) : NaN
            meson_results[meson] = (mass=mass, gamma=gamma, converged=converged, residual=residual,
                                    root_quality=root_quality,
                                    root_diagnostics=root_diagnostics,
                                    threshold=thr, gap=gapv)
        end

        if seed_vec !== nothing && isfinite(seed_vec[1]) && isfinite(seed_vec[2])
            meson_seed_out[meson] = seed_vec
        end
    end

    _annotate_mixed_branch_scores(meson_results, quark_params, thermo_params, Float64(k_norm))
    mixed_seed_tracking_out = nothing
    if mixed_branch_align === :identity
        _apply_identity_alignment!(meson_results)
    elseif mixed_branch_align === :identity_track_label_output
        mixed_seed_tracking_out = _apply_identity_tracking_label_output!(meson_results, tracking_seed_in)
    elseif mixed_branch_align !== :label
        throw(ArgumentError("mixed_branch_align must be :label, :identity or :identity_track_label_output, got $mixed_branch_align"))
    end

    return (
        equilibrium=base,
        quark_params=quark_params,
        thermo_params=thermo_params,
        meson_results=meson_results,
        meson_seed_state=meson_seed_out,
        mixed_seed_tracking=mixed_seed_tracking_out,
    )
end

end # module
