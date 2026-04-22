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

const _BBO_AVAILABLE = let ok = false
    try
        @eval import BlackBoxOptim
        ok = true
    catch
        ok = false
    end
    ok
end

const _WORKFLOW_PARAM_ADAPTERS_PATH = normpath(joinpath(@__DIR__, "WorkflowParamAdapters.jl"))
if !isdefined(Main, :WorkflowParamAdapters)
    Base.include(Main, _WORKFLOW_PARAM_ADAPTERS_PATH)
end
const WorkflowParamAdapters = Main.WorkflowParamAdapters
using .WorkflowParamAdapters: normalize_quark_params, normalize_thermo_params

using Main.Constants_PNJL: ħc_MeV_fm
using ..Models: RootPolicy, solve_root_with_policy, execute_governance_selector, select_residual_min_candidate
using ..Models: HADRON_SEED_5, default_momentum_count, default_theta_count
const DEFAULT_MOMENTUM_COUNT = default_momentum_count()
const DEFAULT_THETA_COUNT = default_theta_count()
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

@inline function _selected_branch_score_for_meson(meson::Symbol, score_plus::Float64, score_minus::Float64)::Float64
    return _label_sign_for_mixed(meson) == 1 ? score_plus : score_minus
end

@inline function _distance_to_box(x::Float64, xmin::Float64, xmax::Float64)::Float64
    if x < xmin
        return xmin - x
    elseif x > xmax
        return x - xmax
    end
    return 0.0
end

function _blackbox_global_fallback_seed(
    objective::Function,
    lower::Vector{Float64},
    upper::Vector{Float64};
    max_steps::Int=80,
)
    fallback = 0.5 .* (lower .+ upper)

    if _BBO_AVAILABLE
        dim = length(lower)
        search_range = collect(zip(lower, upper))
        wrapped_objective = x -> objective(collect(x))

        res = try
            BlackBoxOptim.bboptimize(
                wrapped_objective;
                SearchRange=search_range,
                NumDimensions=dim,
                MaxSteps=max_steps,
                Method=:adaptive_de_rand_1_bin_radiuslimited,
                TraceMode=:silent,
            )
        catch
            try
                BlackBoxOptim.bboptimize(
                    wrapped_objective;
                    SearchRange=search_range,
                    NumDimensions=dim,
                    MaxSteps=max_steps,
                    TraceMode=:silent,
                )
            catch
                nothing
            end
        end

        if res !== nothing
            best = try
                BlackBoxOptim.best_candidate(res)
            catch
                nothing
            end
            if best !== nothing
                return Float64[best...]
            end
        end
    end

    # Fallback when BlackBoxOptim is unavailable or optimizer fails.
    return fallback
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

function _build_tracking_state_strict_sign_binding(
    meson_results::Dict{Symbol,NamedTuple},
)
    tracking_out = Dict{Symbol,Vector{Float64}}()
    for pair in ((:eta, :eta_prime), (:sigma, :sigma_prime))
        m1, m2 = pair
        (haskey(meson_results, m1) && haskey(meson_results, m2)) || continue
        minus_key, plus_key = _tracking_keys_for_mixed_pair(m1)
        e1 = meson_results[m1]
        e2 = meson_results[m2]
        tracking_out[minus_key] = Float64[Float64(e1.mass), Float64(e1.gamma)]
        tracking_out[plus_key] = Float64[Float64(e2.mass), Float64(e2.gamma)]
    end
    return tracking_out
end

@inline function _quality_from_residual_and_converged(residual::Float64, converged::Bool, policy)
    if converged && isfinite(residual) && residual <= Float64(policy.residual_norm_max)
        return :good
    elseif isfinite(residual)
        return :degraded
    end
    return :bad
end

function _collect_joint_mixed_candidates(
    meson::Symbol,
    current_entry,
    prev_seed::Union{Nothing,Vector{Float64}},
    qp,
    tp,
    quark_params,
    thermo_params,
    k_norm::Float64,
    mass_kwargs::NamedTuple,
    policy::NamedTuple,
)
    seeds = Vector{Vector{Float64}}()
    if isfinite(Float64(current_entry.mass)) && isfinite(Float64(current_entry.gamma))
        push!(seeds, Float64[Float64(current_entry.mass), max(Float64(current_entry.gamma), 0.0)])
    end
    if prev_seed !== nothing
        push!(seeds, Float64[Float64(prev_seed[1]), max(Float64(prev_seed[2]), 0.0)])
    end

    thr = mott_threshold_masses(meson, qp).min
    for s in (0.9, 1.0, 1.1)
        push!(seeds, Float64[max(thr * s, 1e-6), 0.0])
    end

    if meson === :eta_prime
        tp = normalize_thermo_params(thermo_params)
        xi = Float64(tp.ξ)
        T_MeV = Float64(tp.T) * ħc_MeV_fm
        resid0 = Float64(current_entry.residual)
        if xi <= -0.25 && 240.0 <= T_MeV <= 260.0 && isfinite(resid0) && resid0 >= 0.03
            for m0 in range(2.82, 2.92; length=5), g0 in range(2.30, 2.50; length=5)
                push!(seeds, Float64[m0, g0])
            end
        end
    end

    uniq = Vector{Vector{Float64}}()
    seen = Set{String}()
    for s in seeds
        key = join(round.(s; digits=8), ",")
        key in seen && continue
        push!(uniq, s)
        push!(seen, key)
    end

    candidates = NamedTuple[]

    push!(candidates, (
        mass=Float64(current_entry.mass),
        gamma=Float64(current_entry.gamma),
        residual=Float64(current_entry.residual),
        converged=Bool(current_entry.converged),
        quality=Symbol(current_entry.root_quality),
        method=:current,
        seed=nothing,
        score_plus=Float64(current_entry.branch_score_plus),
        score_minus=Float64(current_entry.branch_score_minus),
        source=:current,
    ))

    for s in uniq
        for method in (:newton, :trust_region)
            rr = try
                solve_meson_mass(
                    meson,
                    qp,
                    tp;
                    k_norm=k_norm,
                    initial_mass=s[1],
                    initial_gamma=max(s[2], 0.0),
                    merge(mass_kwargs, (method=method,))...,
                )
            catch
                nothing
            end
            rr === nothing && continue
            resid = Float64(rr.residual_norm)
            mass = Float64(rr.mass)
            gamma = Float64(rr.gamma)
            splus, sminus = _mixed_branch_scores(meson, mass, gamma, quark_params, thermo_params, k_norm)
            push!(candidates, (
                mass=mass,
                gamma=gamma,
                residual=resid,
                converged=Bool(rr.converged),
                quality=_quality_from_residual_and_converged(resid, Bool(rr.converged), policy),
                method=method,
                seed=copy(s),
                score_plus=splus,
                score_minus=sminus,
                source=:joint_probe,
            ))
        end
    end

    return candidates
end

function _joint_refine_mixed_pairs!(
    meson_results::Dict{Symbol,NamedTuple},
    quark_params,
    thermo_params,
    k_norm::Float64,
    tracking_seed_in::Dict{Symbol,Vector{Float64}},
    mass_kwargs::NamedTuple,
    policy::NamedTuple,
)
    qp = normalize_quark_params(quark_params)
    tp = normalize_thermo_params(thermo_params)

    for pair in ((:eta, :eta_prime), (:sigma, :sigma_prime))
        m_minus, m_plus = pair
        (haskey(meson_results, m_minus) && haskey(meson_results, m_plus)) || continue

        e_minus = meson_results[m_minus]
        e_plus = meson_results[m_plus]
        both_good = (e_minus.root_quality === :good) && (e_plus.root_quality === :good)
        both_good && continue

        minus_key, plus_key = _tracking_keys_for_mixed_pair(m_minus)
        prev_minus = get(tracking_seed_in, minus_key, nothing)
        prev_plus = get(tracking_seed_in, plus_key, nothing)

        cand_minus = _collect_joint_mixed_candidates(
            m_minus,
            e_minus,
            prev_minus,
            qp,
            tp,
            quark_params,
            thermo_params,
            k_norm,
            mass_kwargs,
            policy,
        )
        cand_plus = _collect_joint_mixed_candidates(
            m_plus,
            e_plus,
            prev_plus,
            qp,
            tp,
            quark_params,
            thermo_params,
            k_norm,
            mass_kwargs,
            policy,
        )

        if isempty(cand_minus) || isempty(cand_plus)
            continue
        end

        w_branch = 1.0
        w_cont = 0.05
        w_phys = 0.5

        tp_local = normalize_thermo_params(thermo_params)
        xi_local = Float64(tp_local.ξ)
        T_local_MeV = Float64(tp_local.T) * ħc_MeV_fm
        roi_pair_mode = (m_minus === :eta && m_plus === :eta_prime && xi_local <= -0.25 && 240.0 <= T_local_MeV <= 260.0 && isfinite(Float64(e_plus.residual)) && Float64(e_plus.residual) >= 0.03)
        w_roi = 20.0

        combo_cost(cminus, cplus) = begin
            branch = w_branch * (Float64(cminus.score_minus) + Float64(cplus.score_plus))
            cont = 0.0
            if prev_minus !== nothing
                cont += hypot(Float64(cminus.mass) - prev_minus[1], Float64(cminus.gamma) - prev_minus[2])
            end
            if prev_plus !== nothing
                cont += hypot(Float64(cplus.mass) - prev_plus[1], Float64(cplus.gamma) - prev_plus[2])
            end
            cont *= w_cont
            pen = 0.0
            for c in (cminus, cplus)
                if c.quality !== :good
                    pen += 5.0 + min(max(Float64(c.residual), 0.0), 1.0) * 10.0
                end
                pen += 0.02 * abs(Float64(c.gamma))
            end
            if prev_plus !== nothing
                gprev = abs(prev_plus[2])
                gnow = abs(Float64(cplus.gamma))
                pen += 0.8 * abs(gnow - gprev)
            end
            if roi_pair_mode
                m = Float64(cplus.mass)
                g = Float64(cplus.gamma)
                roi_dist = _distance_to_box(m, 2.82, 2.92) + _distance_to_box(g, 2.30, 2.50)
                pen += w_roi * roi_dist
            end
            pen *= w_phys
            return branch + cont + pen
        end

        best = (J=Inf, cminus=nothing, cplus=nothing)
        total_combo = 0
        for cminus in cand_minus
            for cplus in cand_plus
                total_combo += 1
                J = combo_cost(cminus, cplus)
                if isfinite(J) && J < best.J
                    best = (J=J, cminus=cminus, cplus=cplus)
                end
            end
        end

        b2_guard_applied = false
        if m_minus === :eta && m_plus === :eta_prime && prev_plus !== nothing && !(best.cplus === nothing)
            prev_g = abs(Float64(prev_plus[2]))
            best_g = abs(Float64(best.cplus.gamma))
            if prev_g >= 2.3 && best_g < 2.0
                alt = (J=Inf, cminus=nothing, cplus=nothing)
                for cminus in cand_minus
                    for cplus in cand_plus
                        gnow = abs(Float64(cplus.gamma))
                        gnow >= 2.3 || continue
                        J = combo_cost(cminus, cplus)
                        if isfinite(J) && J < alt.J
                            alt = (J=J, cminus=cminus, cplus=cplus)
                        end
                    end
                end
                if !(alt.cplus === nothing) && isfinite(alt.J) && alt.J <= best.J + 5.0
                    best = alt
                    b2_guard_applied = true
                end
            end
        end

        current_J = combo_cost(cand_minus[1], cand_plus[1])
        if !(best.cminus === nothing || best.cplus === nothing) && isfinite(best.J) && best.J + 1e-12 < current_J
            cminus = best.cminus
            cplus = best.cplus

            thr_minus = mott_threshold_masses(m_minus, qp)
            gaps_minus = mott_gaps(m_minus, Float64(cminus.mass), qp)
            diag_minus = merge(e_minus.root_diagnostics, (
                selected_method=Symbol(cminus.method),
                governance_selection_reason=:joint_pair_objective_min,
                joint_pair_refine_applied=true,
                joint_pair_objective=Float64(best.J),
                joint_pair_candidate_count=total_combo,
                joint_pair_selection_reason=(roi_pair_mode ? :pair_min_roi : :pair_min),
                b2_guard_applied=b2_guard_applied,
            ))
            meson_results[m_minus] = merge(e_minus, (
                mass=Float64(cminus.mass),
                gamma=Float64(cminus.gamma),
                converged=Bool(cminus.converged),
                residual=Float64(cminus.residual),
                root_quality=Symbol(cminus.quality),
                root_diagnostics=diag_minus,
                threshold=thr_minus,
                gaps=gaps_minus,
            ))

            thr_plus = mott_threshold_masses(m_plus, qp)
            gaps_plus = mott_gaps(m_plus, Float64(cplus.mass), qp)
            diag_plus = merge(e_plus.root_diagnostics, (
                selected_method=Symbol(cplus.method),
                governance_selection_reason=:joint_pair_objective_min,
                joint_pair_refine_applied=true,
                joint_pair_objective=Float64(best.J),
                joint_pair_candidate_count=total_combo,
                joint_pair_selection_reason=(roi_pair_mode ? :pair_min_roi : :pair_min),
                b2_guard_applied=b2_guard_applied,
            ))
            meson_results[m_plus] = merge(e_plus, (
                mass=Float64(cplus.mass),
                gamma=Float64(cplus.gamma),
                converged=Bool(cplus.converged),
                residual=Float64(cplus.residual),
                root_quality=Symbol(cplus.quality),
                root_diagnostics=diag_plus,
                threshold=thr_plus,
                gaps=gaps_plus,
            ))
        end
    end
end

function _apply_eta_prime_b2_post_guard!(
    meson_results::Dict{Symbol,NamedTuple},
    tracking_seed_in::Dict{Symbol,Vector{Float64}},
    quark_params,
    thermo_params,
    k_norm::Float64,
    mass_kwargs::NamedTuple,
    policy::NamedTuple,
)
    haskey(meson_results, :eta_prime) || return
    prev_plus = get(tracking_seed_in, :eta_plus, nothing)
    prev_plus === nothing && return

    e = meson_results[:eta_prime]
    cur_g = abs(Float64(e.gamma))
    prev_g = abs(Float64(prev_plus[2]))
    (isfinite(cur_g) && isfinite(prev_g)) || return

    # B2 guard condition: previous step stayed in high-gamma basin, but current
    # selected point drops into low-gamma basin.
    if !(prev_g >= 2.3 && cur_g < 2.0)
        return
    end

    qp = normalize_quark_params(quark_params)
    tp = normalize_thermo_params(thermo_params)
    seed_m = max(Float64(prev_plus[1]), 1e-6)
    seed_g = max(Float64(prev_plus[2]), 0.0)

    best = nothing
    for method in (:trust_region, :newton)
        rr = try
            solve_meson_mass(
                :eta_prime,
                qp,
                tp;
                k_norm=k_norm,
                initial_mass=seed_m,
                initial_gamma=seed_g,
                merge(mass_kwargs, (method=method,))...,
            )
        catch
            nothing
        end
        rr === nothing && continue

        mass = Float64(rr.mass)
        gamma = Float64(rr.gamma)
        resid = Float64(rr.residual_norm)
        if !(isfinite(mass) && isfinite(gamma) && isfinite(resid))
            continue
        end

        splus, sminus = _mixed_branch_scores(:eta_prime, mass, gamma, quark_params, thermo_params, k_norm)
        score = _selected_branch_score_for_meson(:eta_prime, splus, sminus)
        cand = (
            mass=mass,
            gamma=gamma,
            residual=resid,
            converged=Bool(rr.converged),
            score=score,
            method=method,
        )

        if best === nothing || cand.residual < best.residual
            best = cand
        end
    end

    best === nothing && return
    best_g = abs(Float64(best.gamma))
    best_g >= 2.3 || return

    cur_resid = Float64(e.residual)
    if !(isfinite(best.residual) && (!isfinite(cur_resid) || best.residual <= cur_resid + 1e-6))
        return
    end

    thr = mott_threshold_masses(:eta_prime, qp)
    gaps = mott_gaps(:eta_prime, Float64(best.mass), qp)
    root_quality = _quality_from_residual_and_converged(Float64(best.residual), Bool(best.converged), policy)
    diag = merge(e.root_diagnostics, (
        selected_method=best.method,
        governance_selection_reason=:b2_post_guard,
        b2_guard_applied=true,
    ))

    meson_results[:eta_prime] = merge(e, (
        mass=Float64(best.mass),
        gamma=Float64(best.gamma),
        residual=Float64(best.residual),
        converged=Bool(best.converged),
        root_quality=root_quality,
        root_diagnostics=diag,
        threshold=thr,
        gaps=gaps,
    ))
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
    force_global_fallback::Bool=false,
)
    qp = normalize_quark_params(quark_params)
    tp = normalize_thermo_params(thermo_params)

    function candidate_from_result(method::Symbol, seed_source::Symbol, seed::Vector{Float64}, result)
        if result === nothing
            return (
                method=method,
                seed_source=seed_source,
                seed=copy(seed),
                mass=NaN,
                gamma=NaN,
                converged=false,
                residual_norm=Inf,
                hard_constraint_ok=false,
                failed_constraints=Symbol[:solver_failed],
                pressure=-Inf,
                selection_score=Inf,
                quality_tag=:bad,
            )
        end
        resid = Float64(result.residual_norm)
        conv = Bool(result.converged)
        ok = conv && isfinite(resid) && resid <= Float64(policy.residual_norm_max)
        sel_score = if _is_mixed_meson(meson)
            try
                splus, sminus = _mixed_branch_scores(
                    meson,
                    Float64(result.mass),
                    Float64(result.gamma),
                    qp,
                    tp,
                    Float64(k_norm),
                )
                _selected_branch_score_for_meson(meson, splus, sminus)
            catch
                isfinite(resid) ? resid : Inf
            end
        else
            isfinite(resid) ? resid : Inf
        end
        return (
            method=method,
            seed_source=seed_source,
            seed=copy(seed),
            mass=Float64(result.mass),
            gamma=Float64(result.gamma),
            converged=conv,
            residual_norm=resid,
            hard_constraint_ok=ok,
            failed_constraints=(ok ? Symbol[] : Symbol[:meson_root_unqualified]),
            pressure=-Inf,
            selection_score=sel_score,
            quality_tag=(ok ? :good : (isfinite(resid) ? :degraded : :bad)),
        )
    end

    function governed_probe_candidates(seed0::Vector{Float64}, selected_seed::Vector{Float64}, selected_quality::Symbol)
        selected_quality == :good && return nothing

        continuity_seed = initial_seed
        continuity_weight_base = 0.02
        continuity_weight_active = continuity_weight_base
        second_pass_eta_prime_fine_gamma = false
        continuity_guard_applied = false

        @inline function continuity_distance_from_seed(mass::Float64, gamma::Float64)
            continuity_seed === nothing && return Inf
            (isfinite(mass) && isfinite(gamma)) || return Inf
            return hypot(mass - continuity_seed[1], gamma - continuity_seed[2])
        end

        @inline function continuity_distance_from_candidate(c)
            continuity_seed === nothing && return Inf
            m = Float64(get(c, :mass, NaN))
            g = Float64(get(c, :gamma, NaN))
            if !(isfinite(m) && isfinite(g))
                return Inf
            end
            return continuity_distance_from_seed(m, g)
        end

        m0 = seed0[1]
        extra_seeds = Vector{Vector{Float64}}()
        push!(extra_seeds, copy(seed0))
        push!(extra_seeds, copy(selected_seed))
        push!(extra_seeds, Float64[max(m0 * 0.7, 1e-6), 0.0])
        push!(extra_seeds, Float64[max(m0 * 1.3, 1e-6), 0.0])
        push!(extra_seeds, Float64[max(m0, 1e-6), 0.05])

        if _is_mixed_meson(meson)
            thr = mott_threshold_masses(meson, qp).min
            base = abs(selected_seed[1])
            push!(extra_seeds, Float64[max(base, 1e-6), 0.0])
            push!(extra_seeds, Float64[max(base * 0.9, 1e-6), 0.0])
            push!(extra_seeds, Float64[max(base * 1.1, 1e-6), 0.0])
            push!(extra_seeds, Float64[max(thr, 1e-6), 0.0])
            push!(extra_seeds, Float64[max(thr * 0.9, 1e-6), 0.0])
            push!(extra_seeds, Float64[max(thr * 1.1, 1e-6), 0.0])
        end

        uniq = Vector{Vector{Float64}}()
        seen = Set{String}()
        for s in extra_seeds
            key = join(round.(s; digits=8), ",")
            key in seen && continue
            push!(uniq, s)
            push!(seen, key)
        end

        raw_candidates = NamedTuple[]
        for s in uniq
            for method in (:newton, :trust_region)
                r = run_once(method; m0=s[1], g0=max(s[2], 0.0))
                push!(raw_candidates, candidate_from_result(method, :governed_probe, s, r))
            end
        end

        selector_input = map(enumerate(raw_candidates)) do (i, c)
            (
                converged=c.converged,
                pressure=c.pressure,
                residual_norm=c.residual_norm,
                hard_constraint_ok=c.hard_constraint_ok,
                failed_constraints=c.failed_constraints,
                seed_index=i,
                selection_score=c.selection_score,
                quality_tag=c.quality_tag,
                method=c.method,
                seed_source=c.seed_source,
                mass=c.mass,
                gamma=c.gamma,
            )
        end

        selector_fn = if _is_mixed_meson(meson)
            candidates -> begin
                eligible = [i for i in eachindex(candidates) if Bool(candidates[i].hard_constraint_ok)]
                search = isempty(eligible) ? collect(eachindex(candidates)) : eligible
                seed_mass_ref = abs(seed0[1])
                w_mass = 0.02
                w_gamma = 0.005
                w_cont = continuity_weight_active

                candidate_primary_metric(c) = begin
                    s = Float64(c.selection_score)
                    r = isfinite(s) ? s : Inf
                    if isempty(eligible)
                        m = abs(Float64(get(c, :mass, NaN)))
                        g = abs(Float64(get(c, :gamma, NaN)))
                        if isfinite(m)
                            r += w_mass * abs(m - seed_mass_ref)
                        end
                        if isfinite(g)
                            r += w_gamma * g
                        end
                        d = continuity_distance_from_candidate(c)
                        if isfinite(d)
                            r += w_cont * d
                        end
                    end
                    return r
                end

                candidate_cont_metric(c) = continuity_distance_from_candidate(c)

                best_i = search[1]
                for i in search[2:end]
                    ci = candidates[i]
                    cb = candidates[best_i]
                    si = candidate_primary_metric(ci)
                    sb = candidate_primary_metric(cb)
                    di = candidate_cont_metric(ci)
                    db = candidate_cont_metric(cb)
                    ri = Float64(ci.residual_norm)
                    rb = Float64(cb.residual_norm)
                    if (isfinite(si) && isfinite(sb) && si < sb) ||
                       (isfinite(si) && !isfinite(sb)) ||
                       (si == sb && isfinite(di) && isfinite(db) && di < db) ||
                       (si == sb && di == db && ri < rb)
                        best_i = i
                    end
                end
                (
                    selected_index=best_i,
                    selected_candidate=candidates[best_i],
                    selection_reason=(isempty(eligible) ? :no_candidate_passed_constraints : :mixed_branch_score_min),
                )
            end
        else
            select_residual_min_candidate
        end

        governance_selected = execute_governance_selector(
            selector_input;
            selector=selector_fn,
            residual_norm_max=Float64(policy.residual_norm_max),
            require_converged=Bool(policy.require_converged),
        )

        second_pass_triggered = false
        second_pass_candidate_count = 0
        nudged_restart_applied = false
        roi_rescue_attempted = false
        roi_rescue_applied = false
        roi_rescue_candidate_count = 0
        global_fallback_attempted = false
        global_fallback_applied = false
        global_fallback_candidate_count = 0
        b2_guard_applied = false

        if _is_mixed_meson(meson) && governance_selected.selection_reason === :no_candidate_passed_constraints
            second_pass_triggered = true
            extra2 = Vector{Vector{Float64}}()

            thr2 = mott_threshold_masses(meson, qp).min
            mbase = abs(selected_seed[1])
            gbase = abs(selected_seed[2])
            centers = Float64[]
            append!(centers, (mbase, mbase * 0.95, mbase * 1.05, thr2, thr2 * 0.95, thr2 * 1.05))
            if continuity_seed !== nothing
                mcont = abs(continuity_seed[1])
                append!(centers, (mcont, mcont * 0.95, mcont * 1.05))
            end
            gammas = Float64[0.0, max(gbase, 0.0), max(gbase * 0.6, 0.0), max(gbase * 1.2, 0.0), 0.05, 0.15]

            if meson === :eta_prime
                second_pass_eta_prime_fine_gamma = true
                append!(gammas, [
                    0.8,
                    1.0,
                    1.2,
                    1.4,
                    1.6,
                    1.8,
                    2.0,
                    2.2,
                    2.4,
                    2.6,
                    2.8,
                ])
            end

            for mc in centers, gc in gammas
                push!(extra2, Float64[max(mc, 1e-6), max(gc, 0.0)])
            end

            uniq2 = Vector{Vector{Float64}}()
            seen2 = Set{String}()
            for s in extra2
                key = join(round.(s; digits=8), ",")
                key in seen2 && continue
                push!(uniq2, s)
                push!(seen2, key)
            end

            raw_candidates_2 = NamedTuple[]
            for s in uniq2
                for method in (:newton, :trust_region)
                    r = run_once(method; m0=s[1], g0=max(s[2], 0.0))
                    push!(raw_candidates_2, candidate_from_result(method, :governed_second_pass, s, r))
                end
            end

            second_pass_candidate_count = length(raw_candidates_2)
            if second_pass_candidate_count > 0
                append!(raw_candidates, raw_candidates_2)

                sel_resid = selected === nothing ? Inf : Float64(selected.residual_norm)
                if isfinite(sel_resid)
                    continuity_weight_active = continuity_weight_base * (1.0 + clamp(sel_resid / Float64(policy.residual_norm_max), 0.0, 5.0))
                else
                    continuity_weight_active = continuity_weight_base * 6.0
                end

                selector_input = map(enumerate(raw_candidates)) do (i, c)
                    (
                        converged=c.converged,
                        pressure=c.pressure,
                        residual_norm=c.residual_norm,
                        hard_constraint_ok=c.hard_constraint_ok,
                        failed_constraints=c.failed_constraints,
                        seed_index=i,
                        selection_score=c.selection_score,
                        quality_tag=c.quality_tag,
                        method=c.method,
                        seed_source=c.seed_source,
                        mass=c.mass,
                        gamma=c.gamma,
                    )
                end

                governance_selected = execute_governance_selector(
                    selector_input;
                    selector=selector_fn,
                    residual_norm_max=Float64(policy.residual_norm_max),
                    require_converged=Bool(policy.require_converged),
                )
            end

            # Minimal solver-side nudge restart to escape local attraction basin.
            # Keep this constrained to mixed/no-candidate scenario only.
            chosen_idx0 = Int(governance_selected.selected_index)
            chosen0 = raw_candidates[chosen_idx0]
            base_m = Float64(chosen0.mass)
            base_g = Float64(chosen0.gamma)
            if isfinite(base_m) && isfinite(base_g)
                nudges = (
                    (0.15, 0.8),
                    (0.20, 0.8),
                    (0.25, 0.9),
                    (0.25, 1.0),
                )
                nudge_candidates = NamedTuple[]
                for (dm, dg) in nudges
                    m0n = max(base_m + dm, 1e-6)
                    g0n = max(base_g + dg, 0.0)
                    for method in (:trust_region, :newton)
                        rn = run_once(method; m0=m0n, g0=g0n)
                        push!(nudge_candidates, candidate_from_result(method, :nudge_restart, Float64[m0n, g0n], rn))
                    end
                end
                if !isempty(nudge_candidates)
                    append!(raw_candidates, nudge_candidates)
                    selector_input = map(enumerate(raw_candidates)) do (i, c)
                        (
                            converged=c.converged,
                            pressure=c.pressure,
                            residual_norm=c.residual_norm,
                            hard_constraint_ok=c.hard_constraint_ok,
                            failed_constraints=c.failed_constraints,
                            seed_index=i,
                            selection_score=c.selection_score,
                            quality_tag=c.quality_tag,
                            method=c.method,
                            seed_source=c.seed_source,
                            mass=c.mass,
                            gamma=c.gamma,
                        )
                    end
                    governance_selected_2 = execute_governance_selector(
                        selector_input;
                        selector=selector_fn,
                        residual_norm_max=Float64(policy.residual_norm_max),
                        require_converged=Bool(policy.require_converged),
                    )
                    if governance_selected_2.selected_index != governance_selected.selected_index
                        nudged_restart_applied = true
                    end
                    governance_selected = governance_selected_2
                end
            end

            chosen_idx_roi_base = Int(governance_selected.selected_index)
            chosen_roi_base = raw_candidates[chosen_idx_roi_base]
            chosen_roi_resid = Float64(chosen_roi_base.residual_norm)
            if meson === :eta_prime && isfinite(chosen_roi_resid) && chosen_roi_resid >= 0.015
                roi_rescue_attempted = true
                roi_m_grid = range(2.82, 2.92; length=5)
                roi_g_grid = range(2.30, 2.50; length=5)
                roi_candidates = NamedTuple[]
                for m0n in roi_m_grid, g0n in roi_g_grid
                    for method in (:trust_region, :newton)
                        rn = run_once(method; m0=Float64(m0n), g0=Float64(g0n))
                        push!(roi_candidates, candidate_from_result(method, :roi_rescue, Float64[m0n, g0n], rn))
                    end
                end
                roi_rescue_candidate_count = length(roi_candidates)
                if roi_rescue_candidate_count > 0
                    append!(raw_candidates, roi_candidates)
                    old_c = raw_candidates[Int(governance_selected.selected_index)]
                    old_score = Float64(old_c.selection_score)
                    old_resid = Float64(old_c.residual_norm)
                    old_cont = continuity_distance_from_candidate(old_c)
                    best_roi_idx = 0
                    best_roi_score = Inf
                    for (i, c) in enumerate(raw_candidates)
                        get(c, :seed_source, :unknown) === :roi_rescue || continue
                        s = Float64(c.selection_score)
                        isfinite(s) || continue
                        if s < best_roi_score
                            best_roi_score = s
                            best_roi_idx = i
                        end
                    end

                    if best_roi_idx > 0
                        new_c = raw_candidates[best_roi_idx]
                        new_score = Float64(new_c.selection_score)
                        new_resid = Float64(new_c.residual_norm)
                        new_cont = continuity_distance_from_candidate(new_c)

                        better_score = isfinite(new_score) && (!isfinite(old_score) || new_score + 1e-12 < old_score)
                        continuity_ok = (!isfinite(old_cont) || !isfinite(new_cont) || new_cont <= 3.0 * old_cont || (isfinite(old_resid) && isfinite(new_resid) && new_resid <= 0.9 * old_resid))
                        residual_ok = !isfinite(new_resid) || !isfinite(old_resid) || new_resid <= old_resid + 0.01

                        if better_score && continuity_ok && residual_ok
                            if best_roi_idx != governance_selected.selected_index
                                roi_rescue_applied = true
                            end
                            governance_selected = (
                                selected_index=best_roi_idx,
                                selected_candidate=raw_candidates[best_roi_idx],
                                selection_reason=:roi_rescue_objective_min,
                            )
                        end
                    end
                end
            end

            chosen_idx_global_base = Int(governance_selected.selected_index)
            chosen_global_base = raw_candidates[chosen_idx_global_base]
            chosen_global_resid = Float64(chosen_global_base.residual_norm)
            if meson === :eta_prime && ((isfinite(chosen_global_resid) && chosen_global_resid >= 0.03) || force_global_fallback)
                global_fallback_attempted = true
                global_candidates = NamedTuple[]
                lower = Float64[2.2, 0.0]
                upper = Float64[3.3, 3.0]
                objective = function (x::Vector{Float64})
                    m = clamp(Float64(x[1]), lower[1], upper[1])
                    g = clamp(Float64(x[2]), lower[2], upper[2])
                    r = run_once(:trust_region; m0=m, g0=g)
                    c = candidate_from_result(:trust_region, :global_fallback_bbo_obj, Float64[m, g], r)
                    s = Float64(c.selection_score)
                    cont = continuity_distance_from_candidate(c)
                    return s + (isfinite(cont) ? 0.01 * cont : 0.0)
                end
                best_seed = _blackbox_global_fallback_seed(objective, lower, upper; max_steps=80)
                for method in (:trust_region, :newton)
                    rn = run_once(method; m0=best_seed[1], g0=best_seed[2])
                    push!(global_candidates, candidate_from_result(method, :global_fallback_bbo, Float64[best_seed[1], best_seed[2]], rn))
                end
                global_fallback_candidate_count = length(global_candidates)
                if global_fallback_candidate_count > 0
                    old_c = raw_candidates[Int(governance_selected.selected_index)]
                    old_resid = Float64(old_c.residual_norm)
                    old_cont = continuity_distance_from_candidate(old_c)
                    best_idx = 0
                    best_score = Inf
                    for (i, c) in enumerate(global_candidates)
                        resid = Float64(c.selection_score)
                        isfinite(resid) || continue
                        cont = continuity_distance_from_candidate(c)
                        score = resid + (isfinite(cont) ? 0.01 * cont : 0.0)
                        if score < best_score
                            best_score = score
                            best_idx = i
                        end
                    end
                    if best_idx > 0
                        best_c = global_candidates[best_idx]
                        new_resid = Float64(best_c.residual_norm)
                        new_cont = continuity_distance_from_candidate(best_c)
                        better = isfinite(new_resid) && (!isfinite(old_resid) || new_resid + 1e-12 < old_resid)
                        continuity_ok = (!isfinite(old_cont) || !isfinite(new_cont) || new_cont <= 3.0 * old_cont || new_resid <= 0.8 * old_resid)
                        if better && continuity_ok
                            push!(raw_candidates, best_c)
                            governance_selected = (
                                selected_index=length(raw_candidates),
                                selected_candidate=best_c,
                                selection_reason=:global_fallback_min,
                            )
                            global_fallback_applied = true
                        end
                    end
                end
            end
        end

        chosen_idx = Int(governance_selected.selected_index)
        chosen = raw_candidates[chosen_idx]

        if _is_mixed_meson(meson) && continuity_seed !== nothing
            chosen_dist = continuity_distance_from_candidate(chosen)
            if isfinite(chosen_dist)
                filtered = [(i, c) for (i, c) in enumerate(raw_candidates) if isfinite(continuity_distance_from_candidate(c))]
                if !isempty(filtered)
                    dists = [(i=i, d=continuity_distance_from_candidate(c)) for (i, c) in filtered]
                    best = dists[1]
                    for t in dists[2:end]
                        if t.d < best.d
                            best = t
                        end
                    end
                    min_dist = best.d
                    min_dist_idx = best.i
                    dist_ratio = chosen_dist / max(min_dist, 1e-12)
                    if dist_ratio > 3.0
                        alt = raw_candidates[min_dist_idx]
                        alt_resid = Float64(alt.residual_norm)
                        chosen_resid = Float64(chosen.residual_norm)
                        if isfinite(alt_resid) && isfinite(chosen_resid) && alt_resid <= 5.0 * chosen_resid
                            chosen = alt
                            chosen_idx = min_dist_idx
                            continuity_guard_applied = true
                        end
                    end
                end
            end
        end

        if force_global_fallback && meson === :eta_prime && continuity_seed !== nothing
            prev_g = abs(Float64(continuity_seed[2]))
            chosen_g = abs(Float64(chosen.gamma))
            if isfinite(prev_g) && isfinite(chosen_g) && prev_g >= 2.3 && chosen_g < 2.0
                b2_candidates = NamedTuple[]
                seed_b2 = Float64[max(continuity_seed[1], 1e-6), max(continuity_seed[2], 0.0)]
                for method in (:trust_region, :newton)
                    rb2 = run_once(method; m0=seed_b2[1], g0=seed_b2[2])
                    push!(b2_candidates, candidate_from_result(method, :b2_continuity_guard, seed_b2, rb2))
                end

                best_i = 0
                best_resid = Inf
                for (i, c) in enumerate(b2_candidates)
                    g = abs(Float64(c.gamma))
                    r = Float64(c.residual_norm)
                    (isfinite(g) && isfinite(r) && g >= 2.3) || continue
                    if r < best_resid
                        best_resid = r
                        best_i = i
                    end
                end

                if best_i > 0
                    chosen_resid = Float64(chosen.residual_norm)
                    if isfinite(best_resid) && (!isfinite(chosen_resid) || best_resid <= chosen_resid + 1e-6)
                        append!(raw_candidates, b2_candidates)
                        chosen_idx = length(raw_candidates) - length(b2_candidates) + best_i
                        chosen = raw_candidates[chosen_idx]
                        governance_selected = (
                            selected_index=chosen_idx,
                            selected_candidate=chosen,
                            selection_reason=:b2_continuity_guard,
                        )
                        b2_guard_applied = true
                    end
                end
            end
        end

        chosen_result = run_once(chosen.method; m0=chosen.seed[1], g0=max(chosen.seed[2], 0.0))
        return (
            probe_selected_result=chosen_result,
            selected_method=chosen.method,
            selected_seed=chosen.seed,
            selection_reason=governance_selected.selection_reason,
            candidate_count=length(raw_candidates),
            second_pass_triggered=second_pass_triggered,
            second_pass_candidate_count=second_pass_candidate_count,
            nudged_restart_applied=nudged_restart_applied,
            roi_rescue_attempted=roi_rescue_attempted,
            roi_rescue_applied=roi_rescue_applied,
            roi_rescue_candidate_count=roi_rescue_candidate_count,
            global_fallback_attempted=global_fallback_attempted,
            global_fallback_applied=global_fallback_applied,
            global_fallback_candidate_count=global_fallback_candidate_count,
            b2_guard_applied=b2_guard_applied,
            continuity_penalty_weight=continuity_weight_active,
            continuity_penalty_distance=continuity_distance_from_seed(Float64(chosen.mass), Float64(chosen.gamma)),
            second_pass_eta_prime_fine_gamma=second_pass_eta_prime_fine_gamma,
            continuity_guard_applied=continuity_guard_applied,
        )
    end

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
                _selected_branch_score_for_meson(meson, splus, sminus)
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

    selected_method = solved.diagnostics.attempts[solved.diagnostics.selected_attempt].method
    selected_seed = Float64[solved.x[1], solved.x[2]]
    selected = run_once(selected_method; m0=selected_seed[1], g0=max(selected_seed[2], 0.0))

    if selected === nothing
        if isfinite(solved.x[1]) && isfinite(solved.x[2])
            selected = (
                mass=Float64(solved.x[1]),
                gamma=Float64(solved.x[2]),
                converged=Bool(solved.converged),
                residual_norm=Float64(solved.residual_norm),
                solution=nothing,
            )
        end
    end

    governance_probe = governed_probe_candidates(seed0, selected_seed, root_quality)

    if governance_probe !== nothing && governance_probe.probe_selected_result !== nothing
        probe_selected = governance_probe.probe_selected_result
        probe_good = _mass_result_good(probe_selected, policy)
        current_good = selected === nothing ? false : _mass_result_good(selected, policy)
        probe_resid = Float64(probe_selected.residual_norm)
        current_resid = selected === nothing ? Inf : Float64(selected.residual_norm)
        if (probe_good && !current_good) || (probe_good == current_good && isfinite(probe_resid) && (!isfinite(current_resid) || probe_resid < current_resid))
            selected = probe_selected
            selected_method = governance_probe.selected_method
            selected_seed = governance_probe.selected_seed
            root_quality = _mass_result_good(selected, policy) ? :good : :degraded
        end
    end

    root_diag = (
        selected_method=selected_method,
        selected_attempt=solved.diagnostics.selected_attempt,
        attempts=diag_attempts,
        governance_candidate_count=(governance_probe === nothing ? 0 : governance_probe.candidate_count),
        governance_selection_reason=(governance_probe === nothing ? :none : governance_probe.selection_reason),
        second_pass_triggered=(governance_probe === nothing ? false : governance_probe.second_pass_triggered),
        second_pass_candidate_count=(governance_probe === nothing ? 0 : governance_probe.second_pass_candidate_count),
        nudged_restart_applied=(governance_probe === nothing ? false : governance_probe.nudged_restart_applied),
        roi_rescue_attempted=(governance_probe === nothing ? false : governance_probe.roi_rescue_attempted),
        roi_rescue_applied=(governance_probe === nothing ? false : governance_probe.roi_rescue_applied),
        roi_rescue_candidate_count=(governance_probe === nothing ? 0 : governance_probe.roi_rescue_candidate_count),
        global_fallback_attempted=(governance_probe === nothing ? false : governance_probe.global_fallback_attempted),
        global_fallback_applied=(governance_probe === nothing ? false : governance_probe.global_fallback_applied),
        global_fallback_candidate_count=(governance_probe === nothing ? 0 : governance_probe.global_fallback_candidate_count),
        continuity_penalty_weight=(governance_probe === nothing ? 0.0 : governance_probe.continuity_penalty_weight),
        continuity_penalty_distance=(governance_probe === nothing ? Inf : governance_probe.continuity_penalty_distance),
        second_pass_eta_prime_fine_gamma=(governance_probe === nothing ? false : governance_probe.second_pass_eta_prime_fine_gamma),
        continuity_guard_applied=(governance_probe === nothing ? false : governance_probe.continuity_guard_applied),
        b2_guard_applied=(governance_probe === nothing ? false : governance_probe.b2_guard_applied),
        pre_joint_mass=(selected === nothing ? NaN : Float64(selected.mass)),
        pre_joint_gamma=(selected === nothing ? NaN : Float64(selected.gamma)),
        pre_joint_residual=(selected === nothing ? Inf : Float64(selected.residual_norm)),
        post_joint_mass=(selected === nothing ? NaN : Float64(selected.mass)),
        post_joint_gamma=(selected === nothing ? NaN : Float64(selected.gamma)),
        post_joint_residual=(selected === nothing ? Inf : Float64(selected.residual_norm)),
        joint_pair_refine_applied=false,
        joint_pair_objective=Inf,
        joint_pair_candidate_count=0,
        joint_pair_selection_reason=:none,
    )

    if selected === nothing
        return nothing, root_quality, nothing, root_diag
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
    force_global_fallback::Bool=false,
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

    if mixed_branch_align === :identity_track_label_output || mixed_branch_align === :strict_sign_binding
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
            force_global_fallback=force_global_fallback,
        )

        mass = res === nothing ? NaN : Float64(res.mass)
        gamma = res === nothing ? NaN : Float64(res.gamma)
        if _is_mixed_meson(meson) && isfinite(mass) && mass < 0.0
            mass = -mass
            if isfinite(gamma)
                gamma = abs(gamma)
            end
        end
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

    _joint_refine_mixed_pairs!(
        meson_results,
        quark_params,
        thermo_params,
        Float64(k_norm),
        tracking_seed_in,
        mass_kwargs,
        effective_policy,
    )

    if force_global_fallback
        _apply_eta_prime_b2_post_guard!(
            meson_results,
            tracking_seed_in,
            quark_params,
            thermo_params,
            Float64(k_norm),
            mass_kwargs,
            effective_policy,
        )
    end

    _annotate_mixed_branch_scores(meson_results, quark_params, thermo_params, Float64(k_norm))
    for meson in mesons
        _is_mixed_meson(meson) || continue
        haskey(meson_results, meson) || continue
        ent = meson_results[meson]
        diag = ent.root_diagnostics
        meson_results[meson] = merge(ent, (
            root_diagnostics=merge(diag, (
                post_joint_mass=Float64(ent.mass),
                post_joint_gamma=Float64(ent.gamma),
                post_joint_residual=Float64(ent.residual),
            )),
        ))
    end
    mixed_seed_tracking_out = nothing
    if mixed_branch_align === :identity
        _apply_identity_alignment!(meson_results)
    elseif mixed_branch_align === :identity_track_label_output
        mixed_seed_tracking_out = _apply_identity_tracking_label_output!(meson_results, tracking_seed_in)
    elseif mixed_branch_align === :strict_sign_binding
        mixed_seed_tracking_out = _build_tracking_state_strict_sign_binding(meson_results)
    elseif mixed_branch_align !== :label
        throw(ArgumentError("mixed_branch_align must be :label, :identity, :identity_track_label_output or :strict_sign_binding, got $mixed_branch_align"))
    end

    # B2 post-alignment label guard:
    # if previous plus branch was in high-gamma basin but eta_prime label output
    # drops to low-gamma while eta keeps a high-gamma candidate, swap labels back.
    if force_global_fallback && mixed_branch_align === :identity_track_label_output &&
       haskey(meson_results, :eta) && haskey(meson_results, :eta_prime)
        prev_plus = get(tracking_seed_in, :eta_plus, nothing)
        if prev_plus !== nothing && isfinite(prev_plus[2]) && abs(Float64(prev_plus[2])) >= 2.3
            eta_e = meson_results[:eta]
            etap_e = meson_results[:eta_prime]
            etap_g = abs(Float64(etap_e.gamma))
            eta_g = abs(Float64(eta_e.gamma))
            etap_r = Float64(etap_e.residual)
            eta_r = Float64(eta_e.residual)
            if isfinite(etap_g) && isfinite(eta_g) && etap_g < 2.0 && eta_g >= 2.3 &&
               isfinite(eta_r) && isfinite(etap_r) && eta_r <= etap_r + 1e-6
                diag_eta = merge(eta_e.root_diagnostics, (
                    b2_guard_applied=true,
                    governance_selection_reason=:b2_label_guard,
                ))
                diag_etap = merge(etap_e.root_diagnostics, (
                    b2_guard_applied=true,
                    governance_selection_reason=:b2_label_guard,
                ))
                meson_results[:eta] = merge(etap_e, (root_diagnostics=diag_etap,))
                meson_results[:eta_prime] = merge(eta_e, (root_diagnostics=diag_eta,))

                if mixed_seed_tracking_out !== nothing
                    mixed_seed_tracking_out[:eta_minus] = Float64[Float64(meson_results[:eta].mass), Float64(meson_results[:eta].gamma)]
                    mixed_seed_tracking_out[:eta_plus] = Float64[Float64(meson_results[:eta_prime].mass), Float64(meson_results[:eta_prime].gamma)]
                end
            end
        end
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
