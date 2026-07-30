"""Density-support cascade used by the opt-in PNJL phase production path.

The module deliberately contains only curve evidence and point-selection
logic.  Equilibrium solving, Maxwell construction, and solver governance stay
in the Models/scan layers.  This keeps the router usable by both the
production shadow runner and the historical pilot without introducing a
second physical solver.
"""
module RhoSupportRefinement

using LinearAlgebra
using Statistics

using ..Models: detect_s_shape

export RhoSupportConfig,
       RhoSupportPrior,
       RhoSupportAssessment,
       analyze_rho_support_cascade

Base.@kwdef struct RhoSupportConfig
    support_slope_tol::Float64 = 4e-2
    positive_slope_margin::Float64 = 1e-2
    negative_slope_margin::Float64 = 1e-2
    minimum_negative_secant_run::Int = 2
    target_point_count::Int = 9
    max_extra_points::Int = 12
    support_expansion_gaps::Float64 = 0.5
    local_fit_rmse_tol::Float64 = 2e-2
    near_critical_slope_tol::Float64 = 1e-2
end

struct RhoSupportPrior
    T::Float64
    center::Float64
    width::Float64
end

struct RhoSupportAssessment
    status::Symbol
    stage::Symbol
    reason::String
    baseline_has_s_shape::Bool
    support_origin::Symbol
    support_low::Union{Nothing, Float64}
    support_high::Union{Nothing, Float64}
    extra_point_count::Int
    polynomial_fit_count::Int
    targeted_min_secant_slope::Union{Nothing, Float64}
    fit_min_derivative::Union{Nothing, Float64}
    fit_rmse::Union{Nothing, Float64}
    fit_has_s_topology::Bool
    spinodal_rho_center::Union{Nothing, Float64}
    spinodal_rho_gap::Union{Nothing, Float64}
    coarse_point_count::Int
    total_point_count::Int
end

function _validate(config::RhoSupportConfig)
    config.support_slope_tol > 0 || throw(ArgumentError("support_slope_tol must be positive"))
    config.positive_slope_margin > 0 || throw(ArgumentError("positive_slope_margin must be positive"))
    config.negative_slope_margin > 0 || throw(ArgumentError("negative_slope_margin must be positive"))
    config.minimum_negative_secant_run >= 2 || throw(ArgumentError("minimum_negative_secant_run must be at least 2"))
    config.target_point_count >= 5 && isodd(config.target_point_count) ||
        throw(ArgumentError("target_point_count must be an odd integer of at least 5"))
    config.max_extra_points >= config.target_point_count ||
        throw(ArgumentError("max_extra_points must cover the targeted grid"))
    config.support_expansion_gaps >= 0 || throw(ArgumentError("support_expansion_gaps must be nonnegative"))
    config.local_fit_rmse_tol > 0 || throw(ArgumentError("local_fit_rmse_tol must be positive"))
    config.near_critical_slope_tol > 0 || throw(ArgumentError("near_critical_slope_tol must be positive"))
    return config
end

function _sorted_curve(rho_values, mu_values)
    length(rho_values) == length(mu_values) || throw(ArgumentError("rho and mu lengths must match"))
    length(rho_values) >= 5 || throw(ArgumentError("at least 5 curve points are required"))
    rho = Float64.(rho_values)
    mu = Float64.(mu_values)
    all(isfinite, rho) && all(isfinite, mu) || throw(ArgumentError("rho and mu values must be finite"))
    order = sortperm(rho)
    rho = rho[order]
    mu = mu[order]
    all(diff(rho) .> 0) || throw(ArgumentError("rho values must be unique and strictly increasing"))
    return rho, mu
end

function _normalize_curve(rho::Vector{Float64}, mu::Vector{Float64})
    rho_center = (first(rho) + last(rho)) / 2
    rho_scale = (last(rho) - first(rho)) / 2
    rho_scale > 0 || throw(ArgumentError("rho span must be positive"))
    mu_center = median(mu)
    mu_scale = max(maximum(abs.(mu .- mu_center)), 100eps(Float64))
    return (; x=(rho .- rho_center) ./ rho_scale,
            y=(mu .- mu_center) ./ mu_scale,
            rho_center, rho_scale, mu_center, mu_scale)
end

@inline function _polyval(coefficients, x::Float64)
    value = 0.0
    for coefficient in Iterators.reverse(coefficients)
        value = muladd(value, x, coefficient)
    end
    return value
end

@inline function _polyderivative(coefficients, x::Float64)
    value = 0.0
    for power in reverse(1:(length(coefficients) - 1))
        value = muladd(value, x, power * coefficients[power + 1])
    end
    return value
end

function _fit_cubic(x::Vector{Float64}, y::Vector{Float64})
    length(x) >= 4 || return nothing
    matrix = Matrix{Float64}(undef, length(x), 4)
    for row in eachindex(x)
        matrix[row, 1] = 1.0
        matrix[row, 2] = x[row]
        matrix[row, 3] = x[row]^2
        matrix[row, 4] = x[row]^3
    end
    coefficients = qr(matrix) \ y
    all(isfinite, coefficients) || return nothing
    return coefficients
end

function _has_secant_topology(x::Vector{Float64}, y::Vector{Float64}, margin::Float64, minimum_run::Int)
    slopes = diff(y) ./ diff(x)
    for start_index in eachindex(slopes)
        slopes[start_index] < -margin || continue
        stop_index = start_index
        while stop_index < lastindex(slopes) && slopes[stop_index + 1] < -margin
            stop_index += 1
        end
        stop_index - start_index + 1 >= minimum_run || continue
        left_stop = start_index - 1
        right_start = stop_index + 1
        left_stop >= firstindex(slopes) && any(>(margin), @view slopes[firstindex(slopes):left_stop]) || continue
        right_start <= lastindex(slopes) && any(>(margin), @view slopes[right_start:lastindex(slopes)]) || continue
        return true
    end
    return false
end

function _low_slope_support(rho::Vector{Float64}, mu::Vector{Float64}, config::RhoSupportConfig)
    normalized = _normalize_curve(rho, mu)
    slopes = diff(normalized.y) ./ diff(normalized.x)
    best_start = 0
    best_stop = -1
    best_length = 0
    best_minimum = Inf
    index = firstindex(slopes)
    while index <= lastindex(slopes)
        if abs(slopes[index]) > config.support_slope_tol
            index += 1
            continue
        end
        run_start = index
        while index < lastindex(slopes) && abs(slopes[index + 1]) <= config.support_slope_tol
            index += 1
        end
        run_stop = index
        run_length = run_stop - run_start + 1
        run_minimum = minimum(abs, @view slopes[run_start:run_stop])
        is_interior = run_start > firstindex(slopes) && run_stop < lastindex(slopes)
        if is_interior && (run_length > best_length || (run_length == best_length && run_minimum < best_minimum))
            best_start = run_start
            best_stop = run_stop
            best_length = run_length
            best_minimum = run_minimum
        end
        index += 1
    end
    best_length > 0 || return (available=false, low=nothing, high=nothing, width=0.0, min_abs_slope=minimum(abs, slopes))
    low = rho[best_start]
    high = rho[best_stop + 1]
    return (available=true, low=low, high=high, width=high - low, min_abs_slope=best_minimum)
end

function _support_window(rho, support, prior, config::RhoSupportConfig)
    gap = maximum(diff(rho))
    expansion = config.support_expansion_gaps * gap
    current = support.available ?
        (low=max(first(rho), Float64(support.low) - expansion), high=min(last(rho), Float64(support.high) + expansion)) : nothing
    predicted = prior === nothing ? nothing :
        (low=max(first(rho), prior.center - 0.5prior.width - expansion), high=min(last(rho), prior.center + 0.5prior.width + expansion))
    if current !== nothing && predicted !== nothing
        low = max(current.low, predicted.low)
        high = min(current.high, predicted.high)
        high > low && return (available=true, low=low, high=high, origin=:coarse_and_prior)
        return (available=true, low=current.low, high=current.high, origin=:coarse_low_slope)
    elseif current !== nothing
        return (available=true, low=current.low, high=current.high, origin=:coarse_low_slope)
    elseif predicted !== nothing && predicted.high > predicted.low
        return (available=true, low=predicted.low, high=predicted.high, origin=:temperature_prior)
    end
    return (available=false, low=nothing, high=nothing, origin=:none)
end

function _targeted_points(rho::Vector{Float64}, low::Float64, high::Float64, config::RhoSupportConfig)
    candidates = collect(range(low, high; length=config.target_point_count))
    scale = max(last(rho) - first(rho), 1.0)
    atol = 100eps(Float64) * scale
    additions = Float64[]
    for value in candidates
        any(existing -> abs(existing - value) <= atol, rho) && continue
        push!(additions, value)
        length(additions) >= config.max_extra_points && break
    end
    return additions
end

function _local_cubic_evidence(rho::Vector{Float64}, mu::Vector{Float64}, config::RhoSupportConfig)
    normalized = _normalize_curve(rho, mu)
    coefficients = _fit_cubic(normalized.x, normalized.y)
    coefficients === nothing && return (available=false, min_derivative=nothing, rmse=nothing, has_s_topology=false)
    fitted = [_polyval(coefficients, value) for value in normalized.x]
    residuals = fitted .- normalized.y
    derivative_grid = range(-1.0, 1.0; length=201)
    derivatives = [_polyderivative(coefficients, value) for value in derivative_grid]
    all(isfinite, residuals) && all(isfinite, derivatives) ||
        return (available=false, min_derivative=nothing, rmse=nothing, has_s_topology=false)
    minimum_derivative = minimum(derivatives)
    return (
        available=true,
        min_derivative=minimum_derivative,
        rmse=sqrt(mean(abs2, residuals)),
        has_s_topology=first(derivatives) > config.positive_slope_margin &&
            last(derivatives) > config.positive_slope_margin &&
            minimum_derivative < -config.negative_slope_margin,
    )
end

@inline function _spinodal_gap(sres)
    sres.rho_spinodal_hadron === nothing || sres.rho_spinodal_quark === nothing ? nothing :
        abs(Float64(sres.rho_spinodal_hadron) - Float64(sres.rho_spinodal_quark))
end

@inline function _spinodal_center(sres)
    sres.rho_spinodal_hadron === nothing || sres.rho_spinodal_quark === nothing ? nothing :
        0.5(Float64(sres.rho_spinodal_hadron) + Float64(sres.rho_spinodal_quark))
end

function _assessment(status, stage, reason, baseline, origin, low, high, extra, fit_count,
        targeted_slope, fit_slope, fit_rmse, fit_s, center, gap, coarse_count, total_count)
    return RhoSupportAssessment(status, stage, reason, baseline, origin, low, high, extra,
        fit_count, targeted_slope, fit_slope, fit_rmse, fit_s, center, gap, coarse_count, total_count)
end

"""Analyze a sampled curve and optionally request targeted rho points.

The returned assessment is a router/evidence record.  Callers must re-run
their physical S-shape, Maxwell, and geometry gates on the augmented curve;
this function never upgrades a weak or unresolved curve by itself.
"""
function analyze_rho_support_cascade(rho_values, mu_values;
        sample_mu::Union{Nothing, Function}=nothing,
        prior::Union{Nothing, RhoSupportPrior}=nothing,
        config::RhoSupportConfig=RhoSupportConfig())
    _validate(config)
    rho, mu = _sorted_curve(rho_values, mu_values)
    if prior !== nothing
        all(isfinite, (prior.T, prior.center, prior.width)) || throw(ArgumentError("rho-support prior fields must be finite"))
        prior.width > 0 || throw(ArgumentError("rho-support prior width must be positive"))
    end
    normalized = _normalize_curve(rho, mu)
    baseline = detect_s_shape(mu, rho)
    coarse_topology = _has_secant_topology(normalized.x, normalized.y,
        config.negative_slope_margin, config.minimum_negative_secant_run)
    if baseline.has_s_shape && coarse_topology
        return _assessment(:resolved_s_shape, :fast_topology,
            "coarse sampled curve resolves persistent + -> - -> + topology", true,
            :coarse_topology, nothing, nothing, 0, 0, nothing, nothing, nothing, false,
            _spinodal_center(baseline), _spinodal_gap(baseline), length(rho), length(rho))
    end

    support = _low_slope_support(rho, mu, config)
    window = _support_window(rho, support, prior, config)
    if !window.available
        return _assessment(:unresolved, :fast_unresolved,
            baseline.has_s_shape ? "discrete topology lacks persistent density support" :
            "no interior low-slope density support or temperature prior", baseline.has_s_shape,
            :none, nothing, nothing, 0, 0, nothing, nothing, nothing, false, nothing, nothing,
            length(rho), length(rho))
    end
    if sample_mu === nothing
        return _assessment(:near_critical, :rho_support_warning,
            "rho support exists but no targeted sampler is available", baseline.has_s_shape,
            window.origin, window.low, window.high, 0, 0, nothing, nothing, nothing, false,
            nothing, nothing, length(rho), length(rho))
    end

    targeted_rho = _targeted_points(rho, window.low, window.high, config)
    isempty(targeted_rho) && return _assessment(:near_critical, :rho_support_warning,
        "rho support contains no unevaluated targeted points", baseline.has_s_shape,
        window.origin, window.low, window.high, 0, 0, nothing, nothing, nothing, false,
        nothing, nothing, length(rho), length(rho))
    targeted_mu = try
        Float64[Float64(sample_mu(value)) for value in targeted_rho]
    catch
        return _assessment(:unresolved, :rho_support_targeted,
            "targeted sampler failed while probing rho support", baseline.has_s_shape,
            window.origin, window.low, window.high, 0, 0, nothing, nothing, nothing, false,
            nothing, nothing, length(rho), length(rho))
    end
    if !all(isfinite, targeted_mu)
        return _assessment(:unresolved, :rho_support_targeted,
            "targeted sampler returned non-finite mu values", baseline.has_s_shape,
            window.origin, window.low, window.high, length(targeted_rho), 0, nothing, nothing, nothing, false,
            nothing, nothing, length(rho), length(rho) + length(targeted_rho))
    end
    local_indices = findall(value -> window.low <= value <= window.high, rho)
    local_rho = vcat(rho[local_indices], targeted_rho)
    local_mu = vcat(mu[local_indices], targeted_mu)
    order = sortperm(local_rho)
    local_rho = local_rho[order]
    local_mu = local_mu[order]
    local_normalized = _normalize_curve(local_rho, local_mu)
    targeted_topology = _has_secant_topology(local_normalized.x, local_normalized.y,
        config.negative_slope_margin, config.minimum_negative_secant_run)
    targeted_slope = minimum(diff(local_normalized.y) ./ diff(local_normalized.x))
    fit = _local_cubic_evidence(local_rho, local_mu, config)
    targeted_sres = detect_s_shape(local_mu, local_rho)
    fit_usable = fit.available && fit.rmse <= config.local_fit_rmse_tol
    status, stage, reason = if targeted_topology && fit_usable && fit.has_s_topology
        (:resolved_s_shape, :rho_support_targeted,
         "targeted rho samples and local cubic agree on persistent + -> - -> + topology")
    elseif prior !== nothing && fit_usable && targeted_slope > config.positive_slope_margin &&
           fit.min_derivative > config.positive_slope_margin && !baseline.has_s_shape
        (:resolved_monotone, :rho_support_targeted,
         "previous-temperature density window is positively sloped after targeted validation")
    elseif fit_usable && abs(fit.min_derivative) <= config.near_critical_slope_tol
        (:near_critical, :rho_support_targeted, "targeted local cubic has a near-zero minimum slope")
    else
        (:unresolved, :rho_support_targeted,
         "targeted topology and local cubic do not form a sufficient certificate")
    end
    return _assessment(status, stage, reason, baseline.has_s_shape, window.origin,
        window.low, window.high, length(targeted_rho), fit.available ? 1 : 0,
        targeted_slope, fit.min_derivative, fit.rmse, fit.has_s_topology,
        status == :resolved_s_shape ? _spinodal_center(targeted_sres) : nothing,
        status == :resolved_s_shape ? _spinodal_gap(targeted_sres) : nothing,
        length(rho), length(rho) + length(targeted_rho))
end

end # module RhoSupportRefinement
