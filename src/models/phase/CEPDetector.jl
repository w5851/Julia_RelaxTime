function _is_valid_s_curve(mu_vals, rho_vals; maxwell_options::NamedTuple=(;))
    sres = detect_s_shape(mu_vals, rho_vals)
    sres.has_s_shape || return (false, nothing, sres)
    mres = maxwell_construction(mu_vals, rho_vals; spinodal_hint=sres, maxwell_options...)
    if mres.converged && mres.mu_transition !== nothing
        return (true, Float64(mres.mu_transition), sres)
    end
    return (false, nothing, sres)
end

function _classify_interpolate_curve(mu_vals, rho_vals; maxwell_options::NamedTuple=(;))
    cres = _classify_s_curve(mu_vals, rho_vals; maxwell_options=maxwell_options)
    if cres.status == :valid
        return (true, :valid, cres.mu_transition, cres.sres)
    elseif cres.status == :weak_s_shape
        return (true, :weak_s_shape, cres.mu_transition, cres.sres)
    end
    return (false, cres.status, nothing, cres.sres)
end

function _classify_s_curve(mu_vals, rho_vals;
        maxwell_options::NamedTuple=(;),
        area_tol_good::Float64=1e-4,
        area_tol_bad::Float64=5e-4)
    sres = detect_s_shape(mu_vals, rho_vals)
    sres.has_s_shape || return (
        status=:invalid,
        mu_transition=nothing,
        sres=sres,
        maxwell=MaxwellResult(),
        area_residual=nothing,
        reason="no_s_shape",
    )

    mres = maxwell_construction(mu_vals, rho_vals; spinodal_hint=sres, maxwell_options...)
    if !(mres.converged && mres.mu_transition !== nothing)
        reason = get(mres.details, :reason, "maxwell_failed")
        if reason == "no_sign_change"
            weak = _weak_s_shape_metrics(mu_vals, rho_vals, sres)
            if _is_weak_s_shape(weak)
                weak_mu = if sres.mu_spinodal_hadron !== nothing && sres.mu_spinodal_quark !== nothing
                    0.5 * (Float64(sres.mu_spinodal_hadron) + Float64(sres.mu_spinodal_quark))
                else
                    nothing
                end
                return (
                    status=:weak_s_shape,
                    mu_transition=weak_mu,
                    sres=sres,
                    maxwell=mres,
                    area_residual=nothing,
                    reason="weak_s_shape_no_sign_change",
                )
            end
        end
        return (
            status=:invalid,
            mu_transition=nothing,
            sres=sres,
            maxwell=mres,
            area_residual=nothing,
            reason=String(reason),
        )
    end

    area = something(mres.area_residual, Inf)
    if area <= area_tol_good
        return (
            status=:valid,
            mu_transition=Float64(mres.mu_transition),
            sres=sres,
            maxwell=mres,
            area_residual=Float64(area),
            reason="ok",
        )
    elseif area >= area_tol_bad
        return (
            status=:invalid,
            mu_transition=Float64(mres.mu_transition),
            sres=sres,
            maxwell=mres,
            area_residual=Float64(area),
            reason="area_residual_too_large",
        )
    end

    return (
        status=:unknown,
        mu_transition=Float64(mres.mu_transition),
        sres=sres,
        maxwell=mres,
        area_residual=Float64(area),
        reason="area_residual_gray_zone",
    )
end

function _weak_s_shape_metrics(mu_vals, rho_vals, sres)
    mu_sorted, rho_sorted = _sort_curve_by_rho(mu_vals, rho_vals)
    min_negative_slope = Inf
    negative_segment_width = 0.0
    negative_segment_count = 0

    for i in 1:(length(rho_sorted) - 1)
        drho = rho_sorted[i + 1] - rho_sorted[i]
        drho == 0 && continue
        slope = (mu_sorted[i + 1] - mu_sorted[i]) / drho
        if slope < 0
            min_negative_slope = min(min_negative_slope, slope)
            negative_segment_width += abs(drho)
            negative_segment_count += 1
        end
    end

    spinodal_mu_gap = if sres.mu_spinodal_hadron !== nothing && sres.mu_spinodal_quark !== nothing
        abs(Float64(sres.mu_spinodal_hadron) - Float64(sres.mu_spinodal_quark))
    else
        Inf
    end
    spinodal_rho_gap = if sres.rho_spinodal_hadron !== nothing && sres.rho_spinodal_quark !== nothing
        abs(Float64(sres.rho_spinodal_hadron) - Float64(sres.rho_spinodal_quark))
    else
        Inf
    end

    return (
        min_negative_slope=min_negative_slope,
        negative_segment_width=negative_segment_width,
        negative_segment_count=negative_segment_count,
        spinodal_mu_gap=spinodal_mu_gap,
        spinodal_rho_gap=spinodal_rho_gap,
    )
end

function _is_weak_s_shape(metrics;
        min_negative_slope_tol::Float64=0.01,
        negative_segment_width_tol::Float64=0.051,
        negative_segment_count_max::Int=1,
        spinodal_mu_gap_tol::Float64=0.01,
        spinodal_rho_gap_tol::Float64=0.1)
    return metrics.negative_segment_count <= negative_segment_count_max ||
           metrics.negative_segment_width <= negative_segment_width_tol ||
           abs(metrics.min_negative_slope) <= min_negative_slope_tol ||
           metrics.spinodal_mu_gap <= spinodal_mu_gap_tol ||
           metrics.spinodal_rho_gap <= spinodal_rho_gap_tol
end

@inline _is_cep_low_side_status(status::Symbol) = status in (:valid, :weak_s_shape)

function _interpolate_curve(curves::Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}, T_target::Float64)
    temps = sort(collect(keys(curves)))
    T_below = nothing
    T_above = nothing

    for T in temps
        if T < T_target
            T_below = T
        elseif T > T_target && T_above === nothing
            T_above = T
        end
    end

    (T_below === nothing || T_above === nothing) && return nothing, nothing
    mu_below, rho_below = curves[T_below]
    mu_above, rho_above = curves[T_above]
    length(rho_below) == length(rho_above) || return nothing, nothing

    α = (T_target - T_below) / (T_above - T_below)
    mu_interp = mu_below .+ α .* (mu_above .- mu_below)
    return mu_interp, rho_below
end

function _interpolate_or_reevaluate_curve(
        curves::Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}},
        T_target::Float64;
        evaluate_at_T::Union{Nothing, Function}=nothing)
    if haskey(curves, T_target)
        return curves[T_target]
    end
    if evaluate_at_T !== nothing
        curve = evaluate_at_T(T_target, 0)
        curve !== nothing && return curve
    end
    interpolated = _interpolate_curve(curves, T_target)
    # `_interpolate_curve` historically returns `(nothing, nothing)` when
    # the target is outside the supplied temperature span. Normalize that
    # sentinel here so the three-state classifier records a genuine unknown
    # rather than attempting `detect_s_shape(nothing, nothing)`.
    interpolated[1] === nothing && return nothing
    return interpolated
end

function _uniform_rho_grid(curves::Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}; atol::Float64=1e-10)
    isempty(curves) && return true
    temps = sort(collect(keys(curves)))
    _, rho_ref = curves[temps[1]]
    for T in temps[2:end]
        _, rho = curves[T]
        length(rho) == length(rho_ref) || return false
        for i in eachindex(rho)
            abs(Float64(rho[i]) - Float64(rho_ref[i])) <= atol || return false
        end
    end
    return true
end

function _lookup_curve(curves::Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}, T::Float64; atol::Float64=1e-10)
    if haskey(curves, T)
        return curves[T]
    end
    for (Tk, curve) in curves
        abs(Tk - T) <= atol && return curve
    end
    return nothing
end

function _evaluate_direct_at_T(curves::Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}, T::Float64;
        evaluate_at_T::Union{Nothing, Function}=nothing,
        max_refine_level::Int=2,
        maxwell_options::NamedTuple=(;),
        area_tol_good::Float64=1e-4,
        area_tol_bad::Float64=5e-4)
    last_unknown = (
        status=:unknown,
        mu_transition=nothing,
        area_residual=nothing,
        level=0,
        reason="no_curve",
    )

    for level in 0:max_refine_level
        curve = nothing
        if level == 0
            curve = _lookup_curve(curves, T)
        end
        if curve === nothing && evaluate_at_T !== nothing
            curve = evaluate_at_T(T, level)
        end
        curve === nothing && continue

        mu_vals, rho_vals = curve
        cres = _classify_s_curve(mu_vals, rho_vals;
            maxwell_options=maxwell_options,
            area_tol_good=area_tol_good,
            area_tol_bad=area_tol_bad)
        if cres.status == :unknown
            last_unknown = (
                status=:unknown,
                mu_transition=cres.mu_transition,
                area_residual=cres.area_residual,
                level=level,
                reason=String(cres.reason),
            )
            continue
        end
        return (
            status=cres.status,
            mu_transition=cres.mu_transition,
            area_residual=cres.area_residual,
            level=level,
            reason=String(cres.reason),
        )
    end

    return last_unknown
end

function _auto_initial_step(temperatures::Vector{Float64}, tol::Float64)
    if length(temperatures) < 2
        return max(0.5, tol * 4)
    end
    diffs = Float64[]
    sizehint!(diffs, max(0, length(temperatures) - 1))
    for i in 1:(length(temperatures) - 1)
        Δ = temperatures[i + 1] - temperatures[i]
        Δ > 0 && push!(diffs, Δ)
    end
    isempty(diffs) && return max(0.5, tol * 4)
    return max(minimum(diffs), tol)
end

function _scan_bracket(temperatures::Vector{Float64}, eval_fn::Function)
    last_valid = nothing
    first_invalid = nothing
    for T in temperatures
        er = eval_fn(T)
        if _is_cep_low_side_status(er.status)
            last_valid = (Float64(T), er.mu_transition)
            first_invalid = nothing
        elseif er.status == :invalid
            if last_valid !== nothing && first_invalid === nothing
                first_invalid = Float64(T)
            end
        end
    end
    return last_valid, first_invalid
end

function _directional_bracket(temperatures::Vector{Float64}, eval_fn::Function;
        direct_start::Symbol=:low,
        direct_initial_step::Float64=NaN,
        direct_expand_factor::Float64=2.0,
    direct_max_expand_steps::Int=8,
    tol::Float64=0.01)
    isempty(temperatures) && return nothing

    T_min = first(temperatures)
    T_max = last(temperatures)
    T_start = if direct_start == :high
        T_max
    elseif direct_start == :mid
        0.5 * (T_min + T_max)
    else
        T_min
    end

    start_res = eval_fn(T_start)
    start_res.status == :unknown && return nothing

    direction = _is_cep_low_side_status(start_res.status) ? 1.0 : -1.0
    step0 = isfinite(direct_initial_step) && direct_initial_step > 0 ? direct_initial_step : _auto_initial_step(temperatures, tol)
    expand_factor = direct_expand_factor > 1.0 ? direct_expand_factor : 2.0

    if direction > 0
        current_valid_T = Float64(T_start)
        current_valid_mu = start_res.mu_transition
        for i in 1:max(1, direct_max_expand_steps)
            step = step0 * expand_factor^(i - 1)
            T_candidate = min(T_max, T_start + step)
            T_candidate <= current_valid_T && break
            er = eval_fn(T_candidate)
            if _is_cep_low_side_status(er.status)
                current_valid_T = Float64(T_candidate)
                current_valid_mu = er.mu_transition
                continue
            elseif er.status == :invalid
                return (
                    T_low=current_valid_T,
                    mu_low=current_valid_mu,
                    T_high=Float64(T_candidate),
                )
            end
        end
    else
        current_invalid_T = Float64(T_start)
        for i in 1:max(1, direct_max_expand_steps)
            step = step0 * expand_factor^(i - 1)
            T_candidate = max(T_min, T_start - step)
            T_candidate >= current_invalid_T && break
            er = eval_fn(T_candidate)
            if er.status == :invalid
                current_invalid_T = Float64(T_candidate)
                continue
            elseif _is_cep_low_side_status(er.status)
                return (
                    T_low=Float64(T_candidate),
                    mu_low=er.mu_transition,
                    T_high=current_invalid_T,
                )
            end
        end
    end

    return nothing
end

function _find_cep_direct_legacy(curves::Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}};
        evaluate_at_T::Union{Nothing, Function}=nothing,
        maxwell_options::NamedTuple=(;),
        tol::Float64=0.01,
        max_bisect_iter::Int=20,
        area_tol_good::Float64=1e-4,
        area_tol_bad::Float64=5e-4,
        max_refine_level::Int=2,
        direct_bracket_mode::Symbol=:directional,
        direct_start::Symbol=:low,
        direct_initial_step::Float64=NaN,
        direct_expand_factor::Float64=2.0,
        direct_max_expand_steps::Int=8,
        direct_fallback_scan::Bool=true)
    isempty(curves) && return CEPResult(method=:no_valid_s_shape)

    if area_tol_good > area_tol_bad
        area_tol_good, area_tol_bad = area_tol_bad, area_tol_good
    end

    temperatures = sort(collect(keys(curves)))
    eval_count = 0
    unknown_count = 0
    cache = Dict{Float64, NamedTuple}()

    eval_fn = function (T::Float64)
        if haskey(cache, T)
            return cache[T]
        end
        er = _evaluate_direct_at_T(curves, T;
            evaluate_at_T=evaluate_at_T,
            max_refine_level=max_refine_level,
            maxwell_options=maxwell_options,
            area_tol_good=area_tol_good,
            area_tol_bad=area_tol_bad)
        eval_count += er.level + 1
        er.status == :unknown && (unknown_count += 1)
        cache[T] = er
        return er
    end

    bracket = nothing
    if direct_bracket_mode == :directional
        bracket = _directional_bracket(temperatures, eval_fn;
            direct_start=direct_start,
            direct_initial_step=direct_initial_step,
            direct_expand_factor=direct_expand_factor,
            direct_max_expand_steps=direct_max_expand_steps,
            tol=tol)
    end

    if bracket === nothing && direct_fallback_scan
        last_valid, first_invalid = _scan_bracket(temperatures, eval_fn)
        if last_valid !== nothing
            T_low0, mu_low0 = last_valid
            if first_invalid === nothing
                return CEPResult(
                    found=true,
                    T_cep_MeV=T_low0,
                    mu_cep_MeV=Float64(something(mu_low0, NaN)),
                    uncertainty_T_MeV=0.0,
                    eval_count=eval_count,
                    unknown_count=unknown_count,
                    reason="directional_bracket_failed_missing_upper",
                    method=:fallback_last_valid_missing_upper,
                )
            end
            bracket = (
                T_low=T_low0,
                mu_low=mu_low0,
                T_high=first_invalid,
            )
        else
            return CEPResult(
                method=:no_valid_s_shape,
                reason="directional_and_scan_no_valid",
                eval_count=eval_count,
                unknown_count=unknown_count,
            )
        end
    end

    bracket === nothing && return CEPResult(
        method=:no_valid_s_shape,
        reason="directional_bracket_not_found",
        eval_count=eval_count,
        unknown_count=unknown_count,
    )

    T_low = Float64(bracket.T_low)
    T_high = Float64(bracket.T_high)
    last_mu = Float64(something(bracket.mu_low, NaN))
    bisect_count = 0

    while (T_high - T_low) > tol && bisect_count < max_bisect_iter
        T_mid = 0.5 * (T_low + T_high)
        eval_mid = eval_fn(T_mid)

        if _is_cep_low_side_status(eval_mid.status)
            T_low = T_mid
            eval_mid.mu_transition !== nothing && (last_mu = Float64(eval_mid.mu_transition))
        elseif eval_mid.status == :invalid
            T_high = T_mid
        else
            unknown_count += 1
            return CEPResult(
                found=true,
                T_cep_MeV=T_low,
                mu_cep_MeV=last_mu,
                uncertainty_T_MeV=0.5 * (T_high - T_low),
                T_bracket_low_MeV=T_low,
                T_bracket_high_MeV=T_high,
                bracket_width_T_MeV=T_high - T_low,
                eval_count=eval_count,
                unknown_count=unknown_count,
                reason=eval_mid.reason,
                method=:fallback_last_valid_unknown_mid,
            )
        end

        bisect_count += 1
    end

    T_cep = 0.5 * (T_low + T_high)
    return CEPResult(
        found=true,
        T_cep_MeV=T_cep,
        mu_cep_MeV=last_mu,
        uncertainty_T_MeV=0.5 * (T_high - T_low),
        T_bracket_low_MeV=T_low,
        T_bracket_high_MeV=T_high,
        bracket_width_T_MeV=T_high - T_low,
        eval_count=eval_count,
        unknown_count=unknown_count,
        method=:direct_bisect_last_valid_maxwell,
    )
end

# Keep the historical private helper name available to downstream diagnostic
# scripts while the public `find_cep` entrypoint uses the three-state path.
function _find_cep_direct(curves::Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}; kwargs...)
    return _find_cep_direct_legacy(curves; kwargs...)
end

function _find_cep_legacy(curves::Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}};
        maxwell_options::NamedTuple=(;),
        tol::Float64=0.01,
        max_bisect_iter::Int=20,
        strategy::Symbol=:interpolate,
        evaluate_at_T::Union{Nothing, Function}=nothing,
        area_tol_good::Float64=1e-4,
        area_tol_bad::Float64=5e-4,
    max_refine_level::Int=2,
    direct_bracket_mode::Symbol=:directional,
    direct_start::Symbol=:low,
    direct_initial_step::Float64=NaN,
    direct_expand_factor::Float64=2.0,
    direct_max_expand_steps::Int=8,
    direct_fallback_scan::Bool=true)
    isempty(curves) && return CEPResult()

    if strategy == :direct
        return _find_cep_direct_legacy(curves;
            evaluate_at_T=evaluate_at_T,
            maxwell_options=maxwell_options,
            tol=tol,
            max_bisect_iter=max_bisect_iter,
            area_tol_good=area_tol_good,
            area_tol_bad=area_tol_bad,
                max_refine_level=max_refine_level,
                direct_bracket_mode=direct_bracket_mode,
                direct_start=direct_start,
                direct_initial_step=direct_initial_step,
                direct_expand_factor=direct_expand_factor,
                direct_max_expand_steps=direct_max_expand_steps,
                direct_fallback_scan=direct_fallback_scan)
    end

    temperatures = sort(collect(keys(curves)))
    last_with_s = nothing
    first_without_s = nothing
    last_status = :invalid

    for T in temperatures
        mu_vals, rho_vals = curves[T]
        valid, status, mu_transition, _ = _classify_interpolate_curve(mu_vals, rho_vals; maxwell_options=maxwell_options)
        if valid
            last_with_s = (Float64(T), status, Float64(something(mu_transition, NaN)))
            last_status = status
            first_without_s = nothing
        elseif last_with_s !== nothing && first_without_s === nothing
            first_without_s = Float64(T)
        end
    end

    last_with_s === nothing && return CEPResult(method=:no_valid_s_shape)

    T_low, low_status, mu_low = last_with_s
    T_high = first_without_s === nothing ? T_low + max(0.1, tol) : first_without_s
    last_mu = mu_low
    bisect_count = 0

    if !_uniform_rho_grid(curves)
        return CEPResult(found=true, T_cep_MeV=T_low, mu_cep_MeV=last_mu, method=:fallback_last_valid_nonuniform_grid)
    end

    while (T_high - T_low) > tol && bisect_count < max_bisect_iter
        T_mid = 0.5 * (T_low + T_high)
        mu_vals, rho_vals = _interpolate_or_reevaluate_curve(curves, T_mid; evaluate_at_T=evaluate_at_T)

        if mu_vals === nothing || rho_vals === nothing
            return CEPResult(found=true, T_cep_MeV=T_low, mu_cep_MeV=last_mu, method=:fallback_last_valid_interpolation_failed)
        end

        valid, status, mu_transition, _ = _classify_interpolate_curve(mu_vals, rho_vals; maxwell_options=maxwell_options)
        curves[T_mid] = (mu_vals, rho_vals)

        if valid
            T_low = T_mid
            low_status = status
            mu_transition !== nothing && (last_mu = Float64(mu_transition))
        else
            T_high = T_mid
        end
        bisect_count += 1
    end

    T_cep = 0.5 * (T_low + T_high)
    method = low_status == :weak_s_shape ? :bisect_weak_s_shape_disappearance : :bisect_last_valid_maxwell
    return CEPResult(
        found=true,
        T_cep_MeV=T_cep,
        mu_cep_MeV=last_mu,
        uncertainty_T_MeV=0.5 * (T_high - T_low),
        T_bracket_low_MeV=T_low,
        T_bracket_high_MeV=T_high,
        bracket_width_T_MeV=T_high - T_low,
        method=method,
    )
end

# ---------------------------------------------------------------------------
# Three-state CEP contract
# ---------------------------------------------------------------------------

@inline function _research_slice_status(cres, T::Float64, monotone_certificate)
    if cres.status == :valid
        return :confirmed_first_order
    elseif cres.status == :invalid && cres.reason == "no_s_shape" && monotone_certificate !== nothing
        return monotone_certificate(T, cres) ? :confirmed_monotone : :ambiguous_near_critical
    end
    return :ambiguous_near_critical
end

function _three_state_cep_result(low, high;
        eval_count::Int=0,
        unknown_count::Int=0,
        resolution_target::Float64=NaN,
        reason::Union{Nothing, String}=nothing,
        method::Symbol=:three_state_frontier)
    if low === nothing
        return CEPResult(
            result_status=:not_found,
            eval_count=eval_count,
            unknown_count=unknown_count,
            reason=isnothing(reason) ? "no_confirmed_first_order" : reason,
            method=method,
            temperature_resolution_target_MeV=resolution_target,
        )
    end

    T_low = Float64(low.T)
    mu_low = Float64(something(low.mu, NaN))
    T_high = high === nothing ? NaN : Float64(high.T)
    width = isfinite(T_high) ? T_high - T_low : NaN
    return CEPResult(
        result_status=:ambiguous,
        T_bracket_low_MeV=T_low,
        T_bracket_high_MeV=T_high,
        bracket_width_T_MeV=width,
        T_last_first_order_MeV=T_low,
        mu_last_first_order_MeV=mu_low,
        T_first_monotone_MeV=T_high,
        ambiguity_width_T_MeV=width,
        temperature_resolution_target_MeV=resolution_target,
        eval_count=eval_count,
        unknown_count=unknown_count,
        reason=isnothing(reason) ? "ambiguous_interval_between_phase_evidence" : reason,
        method=method,
    )
end

function find_cep(curves::Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}};
        maxwell_options::NamedTuple=(;),
        tol::Float64=0.01,
        max_bisect_iter::Int=20,
        strategy::Symbol=:interpolate,
        evaluate_at_T::Union{Nothing, Function}=nothing,
        area_tol_good::Float64=1e-4,
        area_tol_bad::Float64=5e-4,
        max_refine_level::Int=2,
        direct_bracket_mode::Symbol=:directional,
        direct_start::Symbol=:low,
        direct_initial_step::Float64=NaN,
        direct_expand_factor::Float64=2.0,
        direct_max_expand_steps::Int=8,
        direct_fallback_scan::Bool=true,
        monotone_certificate::Union{Nothing, Function}=nothing,
        unknown_budget::Int=typemax(Int))
    isfinite(tol) && tol > 0 || throw(ArgumentError(
        "temperature resolution target must be finite and positive, got $(tol)",
    ))
    unknown_budget >= 0 || throw(ArgumentError(
        "unknown_budget must be nonnegative, got $(unknown_budget)",
    ))
    max_bisect_iter >= 0 || throw(ArgumentError(
        "max_bisect_iter must be nonnegative, got $(max_bisect_iter)",
    ))
    isempty(curves) && return _three_state_cep_result(nothing, nothing;
        resolution_target=tol,
        method=Symbol("three_state_", strategy),
        reason="empty_curve_set")

    area_tol_good > area_tol_bad && ((area_tol_good, area_tol_bad) = (area_tol_bad, area_tol_good))
    temperatures = sort(collect(keys(curves)))
    eval_count = 0
    unknown_count = 0
    evaluations = Dict{Float64, NamedTuple}()
    statuses = Dict{Float64, Symbol}()

    evaluate_at = function (T::Float64)
        key = Float64(T)
        haskey(evaluations, key) && return evaluations[key]
        er = if strategy == :direct
            _evaluate_direct_at_T(curves, key;
                evaluate_at_T=evaluate_at_T,
                max_refine_level=max_refine_level,
                maxwell_options=maxwell_options,
                area_tol_good=area_tol_good,
                area_tol_bad=area_tol_bad)
        else
            curve = _interpolate_or_reevaluate_curve(curves, key; evaluate_at_T=evaluate_at_T)
            if curve === nothing
                (status=:unknown, mu_transition=nothing, area_residual=nothing, level=0, reason="curve_unavailable")
            else
                mu_vals, rho_vals = curve
                cres = _classify_s_curve(mu_vals, rho_vals;
                    maxwell_options=maxwell_options,
                    area_tol_good=area_tol_good,
                    area_tol_bad=area_tol_bad)
                curves[key] = curve
                merge(cres, (level=0,))
            end
        end
        eval_count += er.level + 1
        semantic = _research_slice_status(er, key, monotone_certificate)
        # `unknown_budget` applies only to solver/classification unknowns. A
        # weak-S or Maxwell-gray slice remains ambiguous and is never counted
        # as a budget event or relabelled as monotone.
        er.status == :unknown && (unknown_count += 1)
        statuses[key] = semantic
        evaluations[key] = er
        return er
    end

    for T in temperatures
        evaluate_at(Float64(T))
    end

    low = nothing
    high = nothing
    for T in sort(collect(keys(evaluations)))
        semantic = statuses[T]
        if semantic == :confirmed_first_order
            low = (T=T, mu=evaluations[T].mu_transition)
        elseif semantic == :confirmed_monotone && low !== nothing && high === nothing
            high = (T=T,)
        end
    end

    if low === nothing
        return _three_state_cep_result(nothing, nothing;
            eval_count=eval_count,
            unknown_count=unknown_count,
            resolution_target=tol,
            method=Symbol("three_state_", strategy),
            reason=unknown_count > unknown_budget ?
                "unknown_budget_exceeded:no_confirmed_first_order" : "no_confirmed_first_order")
    end

    if high !== nothing
        low_anchor = low
        high_anchor = high
        low_search_hi = high
        high_search_lo = low
        for _ in 1:max_bisect_iter
            # The initial supplied curves are always classified, but once the
            # budget is exceeded no further expensive re-evaluation is
            # requested. The current evidence interval is retained verbatim.
            unknown_count > unknown_budget && break
            changed = false
            if low_search_hi.T - low_anchor.T > tol
                T_mid = 0.5 * (low_anchor.T + low_search_hi.T)
                er = evaluate_at(T_mid)
                semantic = statuses[T_mid]
                if semantic == :confirmed_first_order
                    low_anchor = (T=T_mid, mu=er.mu_transition)
                else
                    low_search_hi = (T=T_mid,)
                    semantic == :confirmed_monotone && (high_anchor = (T=T_mid,))
                end
                changed = true
            end
            unknown_count > unknown_budget && break
            if high_anchor.T - high_search_lo.T > tol
                T_mid = 0.5 * (high_search_lo.T + high_anchor.T)
                er = evaluate_at(T_mid)
                semantic = statuses[T_mid]
                if semantic == :confirmed_monotone
                    high_anchor = (T=T_mid,)
                else
                    high_search_lo = (T=T_mid,)
                    semantic == :confirmed_first_order && (low_anchor = (T=T_mid, mu=er.mu_transition))
                end
                changed = true
            end
            changed || break
            low = low_anchor
            high = high_anchor
            low_search_hi.T - low_anchor.T <= tol && high_anchor.T - high_search_lo.T <= tol && break
        end
        low, high = low_anchor, high_anchor
    end

    budget_exhausted = unknown_count > unknown_budget
    result_reason = if budget_exhausted
        "unknown_budget_exceeded"
    elseif high === nothing
        "confirmed_first_order_without_monotone_anchor"
    else
        nothing
    end
    return _three_state_cep_result(low, high;
        eval_count=eval_count,
        unknown_count=unknown_count,
        resolution_target=tol,
        method=Symbol("three_state_", strategy),
        reason=result_reason)
end
