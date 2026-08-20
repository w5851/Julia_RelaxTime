module PNJLCriticalityFeasibility

using CSV
using JSON3
using LinearAlgebra
using SHA
using Statistics

export CriticalityConfig,
       CriticalityAssessment,
       RhoSupportConfig,
       RhoSupportPrior,
       CascadeAssessment,
       CepDensityExtrapolation,
       SyntheticCase,
       analyze_curve,
       analyze_curve_cascade,
       default_cases,
       evaluate_case,
       evaluate_cascade_case,
       run_temperature_sequence,
       extrapolate_cep_from_density_width,
       run_suite,
       write_evidence_package

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const VALID_STATES = (
    :resolved_s_shape,
    :near_critical,
    :resolved_monotone,
    :unresolved,
)

Base.@kwdef struct CriticalityConfig
    near_slope_tol::Float64 = 5e-3
    resolved_slope_margin::Float64 = 1e-2
    holdout_rmse_tol::Float64 = 1e-2
    holdout_max_error_tol::Float64 = 3e-2
    derivative_spread_tol::Float64 = 2e-2
    derivative_grid_points::Int = 801
    minimum_validation_points::Int = 33
    validation_gap_ratio::Float64 = 0.5
    fit_window_trim_fraction::Float64 = 0.1
    minimum_negative_secant_run::Int = 3
end

struct FitEvidence
    degree::Int
    window::Symbol
    split::Symbol
    min_derivative::Float64
    holdout_rmse::Float64
    holdout_max_error::Float64
    holdout_count::Int
    has_s_topology::Bool
end

struct CriticalityAssessment
    status::Symbol
    reason::String
    baseline_has_s_shape::Bool
    validation_available::Bool
    validation_min_secant_slope::Union{Nothing, Float64}
    validation_max_secant_slope::Union{Nothing, Float64}
    validation_negative_secant_run::Int
    validation_has_s_topology::Bool
    coarse_negative_secant_run::Int
    coarse_has_s_topology::Bool
    fit_min_derivative_median::Float64
    fit_min_derivative_spread::Float64
    holdout_rmse_max::Float64
    holdout_max_error::Float64
    fit_stable::Bool
    coarse_point_count::Int
    validation_point_count::Int
    fit_count::Int
    fit_holdout_evaluation_count::Int
    fit_s_topology_all::Bool
    baseline_time_ns::Int
    evidence_time_ns::Int
end

Base.@kwdef struct SyntheticCase
    name::String
    linear_coefficient::Float64
    coarse_points::Int
    validation_points::Int = 0
    noise_amplitude::Float64 = 0.0
    expected_status::Symbol
    expected_cascade_status::Union{Nothing, Symbol} = nothing
    description::String
end

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

struct CascadeAssessment
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

struct CepDensityExtrapolation
    found::Bool
    T_cep::Union{Nothing, Float64}
    slope::Union{Nothing, Float64}
    intercept::Union{Nothing, Float64}
    rmse::Union{Nothing, Float64}
    point_count::Int
    reason::String
end

function _ensure_models_loaded!()
    if !isdefined(Main, :Models)
        Base.include(Main, joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
    end
    return Base.invokelatest(() -> getfield(Main, :Models))
end

function _detect_s_shape(mu::Vector{Float64}, rho::Vector{Float64})
    models = _ensure_models_loaded!()
    detector = Base.invokelatest(() -> getfield(models, :detect_s_shape))
    return Base.invokelatest(detector, mu, rho)
end

function _models_rho_support_config(config::RhoSupportConfig)
    models = _ensure_models_loaded!()
    return getfield(models, :RhoSupportConfig)(;
        support_slope_tol=config.support_slope_tol,
        positive_slope_margin=config.positive_slope_margin,
        negative_slope_margin=config.negative_slope_margin,
        minimum_negative_secant_run=config.minimum_negative_secant_run,
        target_point_count=config.target_point_count,
        max_extra_points=config.max_extra_points,
        support_expansion_gaps=config.support_expansion_gaps,
        local_fit_rmse_tol=config.local_fit_rmse_tol,
        near_critical_slope_tol=config.near_critical_slope_tol,
    )
end

function _delegate_rho_support_assessment(result)
    return CascadeAssessment(
        result.status,
        result.stage,
        result.reason,
        result.baseline_has_s_shape,
        result.support_origin,
        result.support_low,
        result.support_high,
        result.extra_point_count,
        result.polynomial_fit_count,
        result.targeted_min_secant_slope,
        result.fit_min_derivative,
        result.fit_rmse,
        result.fit_has_s_topology,
        result.spinodal_rho_center,
        result.spinodal_rho_gap,
        result.coarse_point_count,
        result.total_point_count,
    )
end

function _validate_config(config::CriticalityConfig)
    config.near_slope_tol > 0 || throw(ArgumentError("near_slope_tol must be positive"))
    config.resolved_slope_margin > config.near_slope_tol || throw(ArgumentError(
        "resolved_slope_margin must exceed near_slope_tol",
    ))
    config.holdout_rmse_tol > 0 || throw(ArgumentError("holdout_rmse_tol must be positive"))
    config.holdout_max_error_tol > 0 || throw(ArgumentError("holdout_max_error_tol must be positive"))
    config.derivative_spread_tol > 0 || throw(ArgumentError("derivative_spread_tol must be positive"))
    config.derivative_grid_points >= 101 || throw(ArgumentError("derivative_grid_points must be at least 101"))
    config.minimum_validation_points >= 3 || throw(ArgumentError("minimum_validation_points must be at least 3"))
    0 < config.validation_gap_ratio < 1 || throw(ArgumentError("validation_gap_ratio must lie in (0, 1)"))
    0 < config.fit_window_trim_fraction < 0.25 || throw(ArgumentError(
        "fit_window_trim_fraction must lie in (0, 0.25)",
    ))
    config.minimum_negative_secant_run >= 2 || throw(ArgumentError(
        "minimum_negative_secant_run must be at least 2",
    ))
    return config
end

function _validate_rho_support_config(config::RhoSupportConfig)
    config.support_slope_tol > 0 || throw(ArgumentError("support_slope_tol must be positive"))
    config.positive_slope_margin > 0 || throw(ArgumentError("positive_slope_margin must be positive"))
    config.negative_slope_margin > 0 || throw(ArgumentError("negative_slope_margin must be positive"))
    config.minimum_negative_secant_run >= 2 || throw(ArgumentError(
        "minimum_negative_secant_run must be at least 2",
    ))
    config.target_point_count >= 5 && isodd(config.target_point_count) || throw(ArgumentError(
        "target_point_count must be an odd integer of at least 5",
    ))
    config.max_extra_points >= config.target_point_count || throw(ArgumentError(
        "max_extra_points must cover the targeted grid",
    ))
    config.support_expansion_gaps >= 0 || throw(ArgumentError(
        "support_expansion_gaps must be nonnegative",
    ))
    config.local_fit_rmse_tol > 0 || throw(ArgumentError("local_fit_rmse_tol must be positive"))
    config.near_critical_slope_tol > 0 || throw(ArgumentError(
        "near_critical_slope_tol must be positive",
    ))
    return config
end

function _sorted_finite_curve(rho_values, mu_values; minimum_points::Int=6)
    length(rho_values) == length(mu_values) || throw(ArgumentError("rho and mu lengths must match"))
    length(rho_values) >= minimum_points || throw(ArgumentError("at least $minimum_points curve points are required"))
    all(isfinite, rho_values) || throw(ArgumentError("rho values must be finite"))
    all(isfinite, mu_values) || throw(ArgumentError("mu values must be finite"))

    order = sortperm(rho_values)
    rho = Float64.(rho_values[order])
    mu = Float64.(mu_values[order])
    all(diff(rho) .> 0) || throw(ArgumentError("rho values must be unique"))
    return rho, mu
end

function _normalize_curve(rho::Vector{Float64}, mu::Vector{Float64})
    rho_center = (first(rho) + last(rho)) / 2
    rho_scale = (last(rho) - first(rho)) / 2
    rho_scale > 0 || throw(ArgumentError("rho span must be positive"))
    mu_center = median(mu)
    mu_scale = maximum(abs.(mu .- mu_center))
    mu_scale = max(mu_scale, 100 * eps(Float64))
    x = (rho .- rho_center) ./ rho_scale
    y = (mu .- mu_center) ./ mu_scale
    return (; x, y, rho_center, rho_scale, mu_center, mu_scale)
end

@inline function _polyval(coefficients::Vector{Float64}, x::Float64)
    value = 0.0
    for coefficient in Iterators.reverse(coefficients)
        value = muladd(value, x, coefficient)
    end
    return value
end

@inline function _polyderivative(coefficients::Vector{Float64}, x::Float64)
    value = 0.0
    for power in reverse(1:(length(coefficients) - 1))
        value = muladd(value, x, power * coefficients[power + 1])
    end
    return value
end

function _fit_polynomial(x::Vector{Float64}, y::Vector{Float64}, degree::Int)
    length(x) > degree || throw(ArgumentError("polynomial degree requires more training points"))
    matrix = Matrix{Float64}(undef, length(x), degree + 1)
    for row in eachindex(x)
        matrix[row, 1] = 1.0
        for column in 2:(degree + 1)
            matrix[row, column] = matrix[row, column - 1] * x[row]
        end
    end
    return qr(matrix) \ y
end

function _fit_evidence(
        x::Vector{Float64},
        y::Vector{Float64},
        training_indices::Vector{Int},
        holdout_indices::Vector{Int},
        degree::Int,
        window::Symbol,
        split::Symbol,
        derivative_grid_points::Int,
        derivative_min::Float64,
        derivative_max::Float64,
        slope_margin::Float64)
    coefficients = _fit_polynomial(x[training_indices], y[training_indices], degree)
    all(isfinite, coefficients) || error("polynomial fit produced non-finite coefficients")
    holdout_errors = [
        _polyval(coefficients, x[index]) - y[index]
        for index in holdout_indices
    ]
    derivative_grid = range(derivative_min, derivative_max; length=derivative_grid_points)
    derivatives = [_polyderivative(coefficients, value) for value in derivative_grid]
    all(isfinite, derivatives) || error("polynomial derivative check produced non-finite values")
    all(isfinite, holdout_errors) || error("polynomial holdout check produced non-finite values")
    min_derivative = minimum(derivatives)
    has_s_topology = first(derivatives) > slope_margin &&
                     last(derivatives) > slope_margin &&
                     min_derivative < -slope_margin
    return FitEvidence(
        degree,
        window,
        split,
        min_derivative,
        sqrt(mean(abs2, holdout_errors)),
        maximum(abs, holdout_errors),
        length(holdout_errors),
        has_s_topology,
    )
end

function _collect_fit_evidence(x::Vector{Float64}, y::Vector{Float64}, config::CriticalityConfig)
    evidence = FitEvidence[]
    windows = [(:full, collect(eachindex(x)))]
    trim = max(1, floor(Int, config.fit_window_trim_fraction * length(x)))
    central = collect((firstindex(x) + trim):(lastindex(x) - trim))
    length(central) >= 8 && push!(windows, (:central, central))

    for (window, indices) in windows
        odd = indices[1:2:end]
        even = indices[2:2:end]
        for (training, holdout, split) in ((odd, even, :odd), (even, odd, :even))
            isempty(holdout) && continue
            for degree in (3, 5)
                length(training) > degree || continue
                push!(evidence, _fit_evidence(
                    x,
                    y,
                    training,
                    holdout,
                    degree,
                    window,
                    Symbol(split, :_degree_, degree),
                    config.derivative_grid_points,
                    minimum(x[indices]),
                    maximum(x[indices]),
                    config.resolved_slope_margin,
                ))
            end
        end
    end
    windows_seen = unique(getfield.(evidence, :window))
    degrees_seen = unique(getfield.(evidence, :degree))
    length(windows_seen) >= 2 || throw(ArgumentError(
        "curve does not support full and central fitting-window checks",
    ))
    all(degree -> degree in degrees_seen, (3, 5)) || throw(ArgumentError(
        "curve does not support both cubic and quintic fit checks",
    ))
    return evidence
end

function _maximum_negative_secant_run(
        x::Vector{Float64},
        y::Vector{Float64},
        margin::Float64)
    maximum_run = 0
    current_run = 0
    for slope in diff(y) ./ diff(x)
        if slope < -margin
            current_run += 1
            maximum_run = max(maximum_run, current_run)
        else
            current_run = 0
        end
    end
    return maximum_run
end

function _has_s_secant_topology(
        x::Vector{Float64},
        y::Vector{Float64},
        margin::Float64,
        minimum_negative_run::Int)
    slopes = diff(y) ./ diff(x)
    for start_index in eachindex(slopes)
        slopes[start_index] < -margin || continue
        stop_index = start_index
        while stop_index < lastindex(slopes) && slopes[stop_index + 1] < -margin
            stop_index += 1
        end
        stop_index - start_index + 1 >= minimum_negative_run || continue
        any(>(margin), @view slopes[firstindex(slopes):(start_index - 1)]) || continue
        any(>(margin), @view slopes[(stop_index + 1):lastindex(slopes)]) || continue
        return true
    end
    return false
end

function _longest_low_slope_support(
        rho::Vector{Float64},
        mu::Vector{Float64},
        config::RhoSupportConfig)
    normalization = _normalize_curve(rho, mu)
    slopes = diff(normalization.y) ./ diff(normalization.x)
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
        while index < lastindex(slopes) &&
              abs(slopes[index + 1]) <= config.support_slope_tol
            index += 1
        end
        run_stop = index
        run_length = run_stop - run_start + 1
        run_minimum = minimum(abs, @view slopes[run_start:run_stop])
        is_interior = run_start > firstindex(slopes) && run_stop < lastindex(slopes)
        if is_interior &&
           (run_length > best_length || (run_length == best_length && run_minimum < best_minimum))
            best_start = run_start
            best_stop = run_stop
            best_length = run_length
            best_minimum = run_minimum
        end
        index += 1
    end
    best_length > 0 || return (
        available=false,
        low=nothing,
        high=nothing,
        width=0.0,
        min_abs_slope=minimum(abs, slopes),
    )
    low = rho[best_start]
    high = rho[best_stop + 1]
    return (
        available=true,
        low=low,
        high=high,
        width=high - low,
        min_abs_slope=best_minimum,
    )
end

function _rho_support_window(
        rho::Vector{Float64},
        support,
        prior::Union{Nothing, RhoSupportPrior},
        config::RhoSupportConfig)
    gap = maximum(diff(rho))
    expansion = config.support_expansion_gaps * gap
    current = support.available ? (
        low=max(first(rho), Float64(support.low) - expansion),
        high=min(last(rho), Float64(support.high) + expansion),
    ) : nothing
    predicted = prior === nothing ? nothing : (
        low=max(first(rho), prior.center - 0.5 * prior.width - expansion),
        high=min(last(rho), prior.center + 0.5 * prior.width + expansion),
    )

    if current !== nothing && predicted !== nothing
        overlap_low = max(current.low, predicted.low)
        overlap_high = min(current.high, predicted.high)
        if overlap_high > overlap_low
            return (available=true, low=overlap_low, high=overlap_high, origin=:coarse_and_prior)
        end
        return (available=true, low=current.low, high=current.high, origin=:coarse_low_slope)
    elseif current !== nothing
        return (available=true, low=current.low, high=current.high, origin=:coarse_low_slope)
    elseif predicted !== nothing && predicted.high > predicted.low
        return (available=true, low=predicted.low, high=predicted.high, origin=:temperature_prior)
    end
    return (available=false, low=nothing, high=nothing, origin=:none)
end

function _targeted_rho_points(
        rho::Vector{Float64},
        low::Float64,
        high::Float64,
        config::RhoSupportConfig)
    candidates = collect(range(low, high; length=config.target_point_count))
    scale = max(last(rho) - first(rho), 1.0)
    atol = 100 * eps(Float64) * scale
    additions = Float64[]
    for value in candidates
        any(existing -> abs(existing - value) <= atol, rho) && continue
        push!(additions, value)
        length(additions) >= config.max_extra_points && break
    end
    return additions
end

function _local_cubic_evidence(
        rho::Vector{Float64},
        mu::Vector{Float64},
        config::RhoSupportConfig)
    length(rho) >= 4 || return (
        available=false,
        min_derivative=nothing,
        rmse=nothing,
        has_s_topology=false,
    )
    normalization = _normalize_curve(rho, mu)
    coefficients = _fit_polynomial(normalization.x, normalization.y, 3)
    fitted = [_polyval(coefficients, value) for value in normalization.x]
    residuals = fitted .- normalization.y
    derivative_grid = range(-1.0, 1.0; length=201)
    derivatives = [_polyderivative(coefficients, value) for value in derivative_grid]
    all(isfinite, coefficients) && all(isfinite, residuals) && all(isfinite, derivatives) ||
        error("local cubic evidence produced non-finite values")
    min_derivative = minimum(derivatives)
    return (
        available=true,
        min_derivative=min_derivative,
        rmse=sqrt(mean(abs2, residuals)),
        has_s_topology=first(derivatives) > config.positive_slope_margin &&
            last(derivatives) > config.positive_slope_margin &&
            min_derivative < -config.negative_slope_margin,
    )
end

function _spinodal_rho_gap(sres)
    if sres.rho_spinodal_hadron === nothing || sres.rho_spinodal_quark === nothing
        return nothing
    end
    return abs(Float64(sres.rho_spinodal_hadron) - Float64(sres.rho_spinodal_quark))
end

function _spinodal_rho_center(sres)
    if sres.rho_spinodal_hadron === nothing || sres.rho_spinodal_quark === nothing
        return nothing
    end
    return 0.5 * (
        Float64(sres.rho_spinodal_hadron) + Float64(sres.rho_spinodal_quark)
    )
end

"""
Run the staged rho-support criticality prototype.

The fast path uses only the existing sampled curve. The targeted path calls
`sample_mu(rho)` only inside a density window selected from the current
low-slope support and/or a previous-temperature prior. Local cubic evidence is
a router; `resolved_s_shape` additionally requires sampled `+ -> - -> +`
secant topology.
"""
function analyze_curve_cascade(
        rho_values,
        mu_values;
        sample_mu::Union{Nothing, Function}=nothing,
        prior::Union{Nothing, RhoSupportPrior}=nothing,
        config::RhoSupportConfig=RhoSupportConfig())
    models = _ensure_models_loaded!()
    if isdefined(models, :analyze_rho_support_cascade)
        modern_prior = prior === nothing ? nothing : getfield(models, :RhoSupportPrior)(prior.T, prior.center, prior.width)
        modern_result = Base.invokelatest(
            getfield(models, :analyze_rho_support_cascade),
            rho_values, mu_values;
            sample_mu=sample_mu,
            prior=modern_prior,
            config=_models_rho_support_config(config),
        )
        return _delegate_rho_support_assessment(modern_result)
    end
    _validate_rho_support_config(config)
    rho, mu = _sorted_finite_curve(rho_values, mu_values; minimum_points=5)
    if prior !== nothing
        all(isfinite, (prior.T, prior.center, prior.width)) || throw(ArgumentError(
            "rho-support prior fields must be finite",
        ))
        prior.width > 0 || throw(ArgumentError("rho-support prior width must be positive"))
    end
    normalization = _normalize_curve(rho, mu)
    baseline = _detect_s_shape(mu, rho)
    coarse_topology = _has_s_secant_topology(
        normalization.x,
        normalization.y,
        config.negative_slope_margin,
        config.minimum_negative_secant_run,
    )
    if baseline.has_s_shape && coarse_topology
        return CascadeAssessment(
            :resolved_s_shape,
            :fast_topology,
            "coarse sampled curve resolves persistent + -> - -> + topology",
            true,
            :coarse_topology,
            nothing,
            nothing,
            0,
            0,
            nothing,
            nothing,
            nothing,
            false,
            _spinodal_rho_center(baseline),
            _spinodal_rho_gap(baseline),
            length(rho),
            length(rho),
        )
    end

    support = _longest_low_slope_support(rho, mu, config)
    window = _rho_support_window(rho, support, prior, config)
    if !window.available
        return CascadeAssessment(
            :unresolved,
            :fast_unresolved,
            baseline.has_s_shape ?
                "discrete topology lacks persistent density support" :
                "no interior low-slope density support or temperature prior",
            baseline.has_s_shape,
            :none,
            nothing,
            nothing,
            0,
            0,
            nothing,
            nothing,
            nothing,
            false,
            nothing,
            nothing,
            length(rho),
            length(rho),
        )
    end

    if sample_mu === nothing
        return CascadeAssessment(
            :near_critical,
            :rho_support_warning,
            "rho support exists but no targeted sampler is available",
            baseline.has_s_shape,
            window.origin,
            window.low,
            window.high,
            0,
            0,
            nothing,
            nothing,
            nothing,
            false,
            nothing,
            nothing,
            length(rho),
            length(rho),
        )
    end

    targeted_rho = _targeted_rho_points(rho, window.low, window.high, config)
    isempty(targeted_rho) && return CascadeAssessment(
        :near_critical,
        :rho_support_warning,
        "rho support contains no unevaluated targeted points",
        baseline.has_s_shape,
        window.origin,
        window.low,
        window.high,
        0,
        0,
        nothing,
        nothing,
        nothing,
        false,
        nothing,
        nothing,
        length(rho),
        length(rho),
    )
    targeted_mu = Float64[Float64(sample_mu(value)) for value in targeted_rho]
    all(isfinite, targeted_mu) || throw(ArgumentError("targeted sampler returned non-finite mu values"))
    coarse_local = findall(value -> window.low <= value <= window.high, rho)
    local_rho = vcat(rho[coarse_local], targeted_rho)
    local_mu = vcat(mu[coarse_local], targeted_mu)
    order = sortperm(local_rho)
    local_rho = local_rho[order]
    local_mu = local_mu[order]
    targeted_normalization = _normalize_curve(local_rho, local_mu)
    targeted_topology = _has_s_secant_topology(
        targeted_normalization.x,
        targeted_normalization.y,
        config.negative_slope_margin,
        config.minimum_negative_secant_run,
    )
    targeted_min_secant_slope = minimum(
        diff(targeted_normalization.y) ./ diff(targeted_normalization.x),
    )
    local_fit = _local_cubic_evidence(local_rho, local_mu, config)
    targeted_sres = _detect_s_shape(local_mu, local_rho)

    fit_usable = local_fit.available && local_fit.rmse <= config.local_fit_rmse_tol
    status, stage, reason = if targeted_topology && fit_usable && local_fit.has_s_topology
        (
            :resolved_s_shape,
            :rho_support_targeted,
            "targeted rho samples and local cubic agree on persistent + -> - -> + topology",
        )
    elseif prior !== nothing && fit_usable &&
           targeted_min_secant_slope > config.positive_slope_margin &&
           local_fit.min_derivative > config.positive_slope_margin &&
           !baseline.has_s_shape
        (
            :resolved_monotone,
            :rho_support_targeted,
            "previous-temperature density window is positively sloped after targeted validation",
        )
    elseif fit_usable && abs(local_fit.min_derivative) <= config.near_critical_slope_tol
        (
            :near_critical,
            :rho_support_targeted,
            "targeted local cubic has a near-zero minimum slope",
        )
    else
        (
            :unresolved,
            :rho_support_targeted,
            "targeted topology and local cubic do not form a sufficient certificate",
        )
    end
    return CascadeAssessment(
        status,
        stage,
        reason,
        baseline.has_s_shape,
        window.origin,
        window.low,
        window.high,
        length(targeted_rho),
        local_fit.available ? 1 : 0,
        targeted_min_secant_slope,
        local_fit.min_derivative,
        local_fit.rmse,
        local_fit.has_s_topology,
        status == :resolved_s_shape ? _spinodal_rho_center(targeted_sres) : nothing,
        status == :resolved_s_shape ? _spinodal_rho_gap(targeted_sres) : nothing,
        length(rho),
        length(rho) + length(targeted_rho),
    )
end

function _validation_evidence(
        normalization,
        validation_rho,
        validation_mu,
        coarse_max_gap::Float64,
        config::CriticalityConfig)
    validation_rho === nothing && return (
        available=false,
        minimum_slope=nothing,
        maximum_slope=nothing,
        negative_run=0,
        has_s_topology=false,
    )
    validation_mu === nothing && throw(ArgumentError("validation_mu is required with validation_rho"))
    rho, mu = _sorted_finite_curve(validation_rho, validation_mu; minimum_points=3)
    unavailable = (
        available=false,
        minimum_slope=nothing,
        maximum_slope=nothing,
        negative_run=0,
        has_s_topology=false,
    )
    length(rho) >= config.minimum_validation_points || return unavailable
    first(rho) <= normalization.rho_center - normalization.rho_scale || return unavailable
    last(rho) >= normalization.rho_center + normalization.rho_scale || return unavailable
    maximum(diff(rho)) <= config.validation_gap_ratio * coarse_max_gap || return unavailable
    x = (rho .- normalization.rho_center) ./ normalization.rho_scale
    y = (mu .- normalization.mu_center) ./ normalization.mu_scale
    slopes = diff(y) ./ diff(x)
    return (
        available=true,
        minimum_slope=minimum(slopes),
        maximum_slope=maximum(slopes),
        negative_run=_maximum_negative_secant_run(
            x,
            y,
            config.resolved_slope_margin,
        ),
        has_s_topology=_has_s_secant_topology(
            x,
            y,
            config.resolved_slope_margin,
            config.minimum_negative_secant_run,
        ),
    )
end

function _classify(
        baseline_has_s_shape::Bool,
        fit_evidence::Vector{FitEvidence},
        validation,
        config::CriticalityConfig)
    minima = getfield.(fit_evidence, :min_derivative)
    median_minimum = median(minima)
    spread = maximum(minima) - minimum(minima)
    rmse_max = maximum(getfield.(fit_evidence, :holdout_rmse))
    error_max = maximum(getfield.(fit_evidence, :holdout_max_error))
    fit_stable = spread <= config.derivative_spread_tol &&
                 rmse_max <= config.holdout_rmse_tol &&
                 error_max <= config.holdout_max_error_tol
    all_negative = maximum(minima) < -config.resolved_slope_margin
    all_s_topology = all(getfield.(fit_evidence, :has_s_topology))
    all_positive = minimum(minima) > config.resolved_slope_margin
    near_zero = abs(median_minimum) <= config.near_slope_tol ||
                (minimum(minima) <= config.near_slope_tol && maximum(minima) >= -config.near_slope_tol)

    validation_s_topology = validation.available && validation.has_s_topology
    validation_positive = validation.available &&
        validation.minimum_slope > config.resolved_slope_margin

    status, reason = if fit_stable && all_s_topology && validation_s_topology
        :resolved_s_shape, baseline_has_s_shape ?
            "discrete folds, cross-window fits, and targeted validation agree" :
            "targeted validation confirms a persistent fitted negative-derivative interval"
    elseif fit_stable && validation_positive && all_positive && !baseline_has_s_shape
        :resolved_monotone, "independent dense validation and stable fits give positive slope margins"
    elseif fit_stable && all_negative &&
           (!validation.available || validation.maximum_slope > config.resolved_slope_margin)
        :near_critical, "stable cross-window fits warn of a negative-slope interval without an independent persistence certificate"
    elseif fit_stable && near_zero && !baseline_has_s_shape
        :near_critical, "stable fitted minimum slope is within the criticality tolerance"
    else
        :unresolved, "fit, discrete topology, and independent validation do not form a sufficient certificate"
    end
    status in VALID_STATES || error("invalid criticality state: $status")
    return (;
        status,
        reason,
        median_minimum,
        spread,
        rmse_max,
        error_max,
        fit_stable,
        all_s_topology,
    )
end

"""
Assess one sampled `mu(rho)` curve.

`resolved_monotone` requires an independent, denser validation curve. A fitted
minimum slope is diagnostic evidence only and never certifies monotonicity by
itself.
"""
function analyze_curve(
        rho_values,
        mu_values;
        validation_rho=nothing,
        validation_mu=nothing,
        config::CriticalityConfig=CriticalityConfig())
    _validate_config(config)
    rho, mu = _sorted_finite_curve(rho_values, mu_values; minimum_points=11)
    normalization = _normalize_curve(rho, mu)
    baseline_start = time_ns()
    baseline = _detect_s_shape(mu, rho)
    baseline_time_ns = time_ns() - baseline_start

    evidence_start = time_ns()
    fits = _collect_fit_evidence(normalization.x, normalization.y, config)
    validation = _validation_evidence(
        normalization,
        validation_rho,
        validation_mu,
        maximum(diff(rho)),
        config,
    )
    classification = _classify(baseline.has_s_shape, fits, validation, config)
    evidence_time_ns = time_ns() - evidence_start

    return CriticalityAssessment(
        classification.status,
        classification.reason,
        baseline.has_s_shape,
        validation.available,
        validation.minimum_slope,
        validation.maximum_slope,
        validation.negative_run,
        validation.has_s_topology,
        _maximum_negative_secant_run(
            normalization.x,
            normalization.y,
            config.resolved_slope_margin,
        ),
        _has_s_secant_topology(
            normalization.x,
            normalization.y,
            config.resolved_slope_margin,
            config.minimum_negative_secant_run,
        ),
        classification.median_minimum,
        classification.spread,
        classification.rmse_max,
        classification.error_max,
        classification.fit_stable,
        length(rho),
        validation.available ? length(validation_rho) : 0,
        length(fits),
        sum(getfield.(fits, :holdout_count)),
        classification.all_s_topology,
        baseline_time_ns,
        evidence_time_ns,
    )
end

@inline _cusp_mu(x::Float64, linear_coefficient::Float64) = x^3 + linear_coefficient * x

function _case_curve(case::SyntheticCase; validation::Bool=false)
    count = validation ? case.validation_points : case.coarse_points
    count > 0 || return Float64[], Float64[]
    x = collect(range(-1.0, 1.0; length=count))
    rho = x .+ 1.0
    mu = [_cusp_mu(value, case.linear_coefficient) for value in x]
    if !validation && case.noise_amplitude > 0
        for index in eachindex(mu)
            mu[index] += case.noise_amplitude * sin(17 * x[index] + 0.37)
        end
    end
    return rho, mu
end

@inline function _case_mu(case::SyntheticCase, rho::Float64; include_noise::Bool=true)
    x = rho - 1.0
    value = _cusp_mu(x, case.linear_coefficient)
    if include_noise && case.noise_amplitude > 0
        value += case.noise_amplitude * sin(17 * x + 0.37)
    end
    return value
end

function default_cases()
    return SyntheticCase[
        SyntheticCase(
            name="resolved_visible_s",
            linear_coefficient=-0.10,
            coarse_points=25,
            validation_points=401,
            expected_status=:resolved_s_shape,
            expected_cascade_status=:resolved_s_shape,
            description="well-resolved cusp with explicit discrete folds",
        ),
        SyntheticCase(
            name="hidden_s_warning",
            linear_coefficient=-0.02,
            coarse_points=13,
            expected_status=:near_critical,
            expected_cascade_status=:resolved_s_shape,
            description="coarse samples miss the narrow negative-slope interval",
        ),
        SyntheticCase(
            name="hidden_s_targeted_validation",
            linear_coefficient=-0.02,
            coarse_points=13,
            validation_points=401,
            expected_status=:resolved_s_shape,
            expected_cascade_status=:resolved_s_shape,
            description="targeted dense validation resolves the hidden interval",
        ),
        SyntheticCase(
            name="critical_cusp",
            linear_coefficient=0.0,
            coarse_points=13,
            validation_points=401,
            expected_status=:near_critical,
            expected_cascade_status=:near_critical,
            description="critical cubic with a zero minimum slope",
        ),
        SyntheticCase(
            name="monotone_without_certificate",
            linear_coefficient=0.08,
            coarse_points=13,
            expected_status=:unresolved,
            expected_cascade_status=:unresolved,
            description="positive coarse secants are insufficient to certify monotonicity",
        ),
        SyntheticCase(
            name="monotone_with_certificate",
            linear_coefficient=0.08,
            coarse_points=13,
            validation_points=401,
            expected_status=:resolved_monotone,
            expected_cascade_status=:unresolved,
            description="independent dense validation certifies positive slope margin",
        ),
        SyntheticCase(
            name="noisy_monotone_minimum_bias",
            linear_coefficient=0.08,
            coarse_points=17,
            validation_points=401,
            noise_amplitude=0.03,
            expected_status=:unresolved,
            expected_cascade_status=:unresolved,
            description="coarse noise can create negative secants but must not certify an S shape",
        ),
    ]
end


function evaluate_cascade_case(
        case::SyntheticCase;
        config::RhoSupportConfig=RhoSupportConfig(),
        prior::Union{Nothing, RhoSupportPrior}=nothing)
    rho, mu = _case_curve(case)
    assessment = analyze_curve_cascade(
        rho,
        mu;
        sample_mu=value -> _case_mu(case, value),
        prior=prior,
        config=config,
    )
    expected = something(case.expected_cascade_status, case.expected_status)
    return (
        case_name=case.name,
        expected_status=String(expected),
        observed_status=String(assessment.status),
        matches_expected=assessment.status == expected,
        stage=String(assessment.stage),
        baseline_has_s_shape=assessment.baseline_has_s_shape,
        support_origin=String(assessment.support_origin),
        support_low=something(assessment.support_low, missing),
        support_high=something(assessment.support_high, missing),
        targeted_function_evaluations=assessment.extra_point_count,
        polynomial_fits=assessment.polynomial_fit_count,
        targeted_min_secant_slope=something(assessment.targeted_min_secant_slope, missing),
        fit_min_derivative=something(assessment.fit_min_derivative, missing),
        fit_rmse=something(assessment.fit_rmse, missing),
        fit_has_s_topology=assessment.fit_has_s_topology,
        spinodal_rho_center=something(assessment.spinodal_rho_center, missing),
        spinodal_rho_gap=something(assessment.spinodal_rho_gap, missing),
        coarse_point_count=assessment.coarse_point_count,
        total_point_count=assessment.total_point_count,
        reason=assessment.reason,
        description=case.description,
    )
end


function extrapolate_cep_from_density_width(T_values, rho_widths)
    length(T_values) == length(rho_widths) || throw(ArgumentError(
        "temperature and density-width lengths must match",
    ))
    temperatures = Float64.(collect(T_values))
    widths = Float64.(collect(rho_widths))
    all(isfinite, temperatures) || throw(ArgumentError("temperatures must be finite"))
    all(isfinite, widths) || throw(ArgumentError("density widths must be finite"))
    all(value -> value > 0, widths) || throw(ArgumentError("density widths must be positive"))
    point_count = length(temperatures)
    point_count >= 3 || return CepDensityExtrapolation(
        false,
        nothing,
        nothing,
        nothing,
        nothing,
        point_count,
        "at least three resolved S-shape temperatures are required",
    )
    temperature_span = maximum(temperatures) - minimum(temperatures)
    temperature_span > 100 * eps(Float64) * max(maximum(abs, temperatures), 1.0) ||
        throw(ArgumentError("temperatures must span more than one numerical value"))
    design = hcat(temperatures, ones(point_count))
    slope, intercept = design \ abs2.(widths)
    fitted = slope .* temperatures .+ intercept
    rmse = sqrt(mean(abs2, fitted .- abs2.(widths)))
    all(isfinite, (slope, intercept, rmse)) || error(
        "density-width CEP extrapolation produced non-finite values",
    )
    if slope >= 0
        return CepDensityExtrapolation(
            false,
            nothing,
            slope,
            intercept,
            rmse,
            point_count,
            "density-width squared does not decrease with temperature",
        )
    end
    T_cep = -intercept / slope
    isfinite(T_cep) || error("density-width CEP extrapolation produced a non-finite CEP")
    return CepDensityExtrapolation(
        true,
        T_cep,
        slope,
        intercept,
        rmse,
        point_count,
        "linear density-width-squared extrapolation has a finite decreasing branch",
    )
end


function run_temperature_sequence(;
        T_values::AbstractVector{<:Real}=[0.0, 0.25, 0.5, 0.7, 0.82, 0.9, 0.96, 1.0, 1.04, 1.15],
        critical_temperature::Real=1.0,
        coefficient_scale::Real=0.12,
        coarse_points::Int=13,
        config::RhoSupportConfig=RhoSupportConfig())
    coarse_points >= 5 || throw(ArgumentError("coarse_points must be at least 5"))
    coefficient_scale > 0 || throw(ArgumentError("coefficient_scale must be positive"))
    temperatures = Float64.(collect(T_values))
    all(isfinite, temperatures) || throw(ArgumentError("T_values must be finite"))
    issorted(temperatures) || throw(ArgumentError("T_values must be sorted"))
    allunique(temperatures) || throw(ArgumentError("T_values must be unique"))
    T_critical = Float64(critical_temperature)
    isfinite(T_critical) || throw(ArgumentError("critical_temperature must be finite"))
    scale = Float64(coefficient_scale)
    rho = collect(range(0.0, 2.0; length=coarse_points))
    rows = NamedTuple[]
    prior = nothing
    extrapolation_temperatures = Float64[]
    extrapolation_widths = Float64[]
    critical_atol = 100 * eps(Float64) * max(abs(T_critical), 1.0)

    for T in temperatures
        coefficient = scale * (T - T_critical)
        mu = [_cusp_mu(value - 1.0, coefficient) for value in rho]
        prior_used = prior !== nothing
        assessment = analyze_curve_cascade(
            rho,
            mu;
            sample_mu=value -> _cusp_mu(value - 1.0, coefficient),
            prior=prior,
            config=config,
        )
        expected = if T < T_critical - critical_atol
            :resolved_s_shape
        elseif abs(T - T_critical) <= critical_atol
            :near_critical
        else
            :resolved_monotone
        end
        analytic_width = coefficient < 0 ? 2 * sqrt(-coefficient / 3) : missing
        push!(rows, (
            T=T,
            linear_coefficient=coefficient,
            expected_status=String(expected),
            observed_status=String(assessment.status),
            matches_expected=assessment.status == expected,
            prior_used=prior_used,
            support_origin=String(assessment.support_origin),
            targeted_function_evaluations=assessment.extra_point_count,
            polynomial_fits=assessment.polynomial_fit_count,
            spinodal_rho_center=something(assessment.spinodal_rho_center, missing),
            spinodal_rho_gap=something(assessment.spinodal_rho_gap, missing),
            analytic_spinodal_rho_gap=analytic_width,
            reason=assessment.reason,
        ))
        if assessment.status == :resolved_s_shape &&
           assessment.spinodal_rho_center !== nothing &&
           assessment.spinodal_rho_gap !== nothing
            prior = RhoSupportPrior(
                T,
                assessment.spinodal_rho_center,
                assessment.spinodal_rho_gap,
            )
            push!(extrapolation_temperatures, T)
            push!(extrapolation_widths, assessment.spinodal_rho_gap)
        end
    end
    extrapolation = extrapolate_cep_from_density_width(
        extrapolation_temperatures,
        extrapolation_widths,
    )
    return (; rows, extrapolation)
end

function evaluate_case(case::SyntheticCase; config::CriticalityConfig=CriticalityConfig())
    rho, mu = _case_curve(case)
    validation_rho, validation_mu = _case_curve(case; validation=true)
    assessment = analyze_curve(
        rho,
        mu;
        validation_rho=isempty(validation_rho) ? nothing : validation_rho,
        validation_mu=isempty(validation_mu) ? nothing : validation_mu,
        config=config,
    )
    return (
        case_name=case.name,
        expected_status=String(case.expected_status),
        observed_status=String(assessment.status),
        matches_expected=assessment.status == case.expected_status,
        baseline_has_s_shape=assessment.baseline_has_s_shape,
        validation_available=assessment.validation_available,
        coarse_point_count=assessment.coarse_point_count,
        validation_point_count=assessment.validation_point_count,
        fit_count=assessment.fit_count,
        fit_min_derivative_median=assessment.fit_min_derivative_median,
        fit_min_derivative_spread=assessment.fit_min_derivative_spread,
        validation_min_secant_slope=something(assessment.validation_min_secant_slope, missing),
        validation_max_secant_slope=something(assessment.validation_max_secant_slope, missing),
        validation_negative_secant_run=assessment.validation_negative_secant_run,
        validation_has_s_topology=assessment.validation_has_s_topology,
        coarse_negative_secant_run=assessment.coarse_negative_secant_run,
        coarse_has_s_topology=assessment.coarse_has_s_topology,
        holdout_rmse_max=assessment.holdout_rmse_max,
        holdout_max_error=assessment.holdout_max_error,
        fit_stable=assessment.fit_stable,
        fit_s_topology_all=assessment.fit_s_topology_all,
        fit_holdout_evaluation_count=assessment.fit_holdout_evaluation_count,
        baseline_time_ns=assessment.baseline_time_ns,
        evidence_time_ns=assessment.evidence_time_ns,
        reason=assessment.reason,
        description=case.description,
    )
end

function run_suite(;
        cases::Vector{SyntheticCase}=default_cases(),
        config::CriticalityConfig=CriticalityConfig(),
        cascade_config::RhoSupportConfig=RhoSupportConfig(),
        repetitions::Int=100)
    repetitions >= 1 || throw(ArgumentError("repetitions must be positive"))
    _ensure_models_loaded!()
    for case in cases
        evaluate_case(case; config=config)
        evaluate_cascade_case(case; config=cascade_config)
    end
    rows = [evaluate_case(case; config=config) for case in cases]
    cascade_rows = [evaluate_cascade_case(case; config=cascade_config) for case in cases]
    sequence = run_temperature_sequence(; config=cascade_config)
    cep = sequence.extrapolation
    cep_expected_temperature = 1.0
    cep_absolute_error = cep.found ? abs(cep.T_cep - cep_expected_temperature) : missing
    cep_row = (
        found=cep.found,
        expected_T_cep=cep_expected_temperature,
        extrapolated_T_cep=something(cep.T_cep, missing),
        absolute_error=cep_absolute_error,
        within_synthetic_gate=cep.found && cep_absolute_error <= 0.05,
        slope=something(cep.slope, missing),
        intercept=something(cep.intercept, missing),
        rmse=something(cep.rmse, missing),
        resolved_width_point_count=cep.point_count,
        reason=cep.reason,
    )

    # Models loading and compilation are reported by the separate startup probe.
    baseline_start = time_ns()
    for _ in 1:repetitions, case in cases
        rho, mu = _case_curve(case)
        _detect_s_shape(mu, rho)
    end
    baseline_elapsed_ns = time_ns() - baseline_start

    fast_indices = findall(row -> row.targeted_function_evaluations == 0, cascade_rows)
    fast_cases = cases[fast_indices]
    fast_start = time_ns()
    for _ in 1:repetitions, case in fast_cases
        evaluate_cascade_case(case; config=cascade_config)
    end
    fast_elapsed_ns = time_ns() - fast_start

    cascade_start = time_ns()
    for _ in 1:repetitions, case in cases
        evaluate_cascade_case(case; config=cascade_config)
    end
    cascade_elapsed_ns = time_ns() - cascade_start

    evidence_start = time_ns()
    for _ in 1:repetitions, case in cases
        evaluate_case(case; config=config)
    end
    evidence_elapsed_ns = time_ns() - evidence_start
    total_cases = repetitions * length(cases)
    total_coarse_evaluations = repetitions * sum(case.coarse_points for case in cases)
    total_validation_evaluations = repetitions * sum(case.validation_points for case in cases)
    total_targeted_evaluations = repetitions * sum(
        row.targeted_function_evaluations for row in cascade_rows
    )
    total_cascade_fits = repetitions * sum(row.polynomial_fits for row in cascade_rows)
    total_polynomial_fits = repetitions * sum(row.fit_count for row in rows)
    total_fit_holdout_evaluations = repetitions * sum(
        row.fit_holdout_evaluation_count for row in rows
    )
    performance = [
        (
            method="current_secant_detector",
            evidence_scope="synthetic_same_process",
            comparable=true,
            repetitions=repetitions,
            case_evaluations=total_cases,
            coarse_function_evaluations=total_coarse_evaluations,
            validation_function_evaluations=0,
            targeted_function_evaluations=0,
            polynomial_fits=0,
            fit_holdout_evaluations=0,
            classification_matches=missing,
            performance_verdict="baseline_reference",
            pnjl_residual_evaluations=missing,
            pnjl_jacobian_evaluations=missing,
            pnjl_newton_iterations=missing,
            pnjl_anchor_solves=missing,
            pnjl_branch_points=missing,
            elapsed_ns=baseline_elapsed_ns,
            mean_time_ns=baseline_elapsed_ns / total_cases,
            note="in-process warm synthetic curve classification only",
        ),
        (
            method="rho_support_fast_path",
            evidence_scope="synthetic_same_process",
            comparable=true,
            repetitions=repetitions,
            case_evaluations=repetitions * length(fast_cases),
            coarse_function_evaluations=repetitions * sum(case.coarse_points for case in fast_cases),
            validation_function_evaluations=0,
            targeted_function_evaluations=0,
            polynomial_fits=0,
            fit_holdout_evaluations=0,
            classification_matches=repetitions * sum(
                cascade_rows[index].matches_expected for index in fast_indices
            ),
            performance_verdict="synthetic_candidate",
            pnjl_residual_evaluations=missing,
            pnjl_jacobian_evaluations=missing,
            pnjl_newton_iterations=missing,
            pnjl_anchor_solves=missing,
            pnjl_branch_points=missing,
            elapsed_ns=fast_elapsed_ns,
            mean_time_ns=isempty(fast_cases) ? missing :
                fast_elapsed_ns / (repetitions * length(fast_cases)),
            note="cases resolved without any targeted function evaluation or polynomial fit",
        ),
        (
            method="rho_support_cascade",
            evidence_scope="synthetic_same_process",
            comparable=true,
            repetitions=repetitions,
            case_evaluations=total_cases,
            coarse_function_evaluations=total_coarse_evaluations,
            validation_function_evaluations=0,
            targeted_function_evaluations=total_targeted_evaluations,
            polynomial_fits=total_cascade_fits,
            fit_holdout_evaluations=0,
            classification_matches=repetitions * sum(row.matches_expected for row in cascade_rows),
            performance_verdict="synthetic_candidate",
            pnjl_residual_evaluations=missing,
            pnjl_jacobian_evaluations=missing,
            pnjl_newton_iterations=missing,
            pnjl_anchor_solves=missing,
            pnjl_branch_points=missing,
            elapsed_ns=cascade_elapsed_ns,
            mean_time_ns=cascade_elapsed_ns / total_cases,
            note="coarse rho support routes at most twelve new fixed-rho samples and one local cubic fit",
        ),
        (
            method="multi_evidence_criticality",
            evidence_scope="synthetic_same_process",
            comparable=true,
            repetitions=repetitions,
            case_evaluations=total_cases,
            coarse_function_evaluations=total_coarse_evaluations,
            validation_function_evaluations=total_validation_evaluations,
            targeted_function_evaluations=0,
            polynomial_fits=total_polynomial_fits,
            fit_holdout_evaluations=total_fit_holdout_evaluations,
            classification_matches=repetitions * sum(row.matches_expected for row in rows),
            performance_verdict="performance_rejected",
            pnjl_residual_evaluations=missing,
            pnjl_jacobian_evaluations=missing,
            pnjl_newton_iterations=missing,
            pnjl_anchor_solves=missing,
            pnjl_branch_points=missing,
            elapsed_ns=evidence_elapsed_ns,
            mean_time_ns=evidence_elapsed_ns / total_cases,
            note="includes polynomial fits and optional independent validation",
        ),
        (
            method="historical_multibranch_palc_cold",
            evidence_scope="historical_noncomparable",
            comparable=false,
            repetitions=1,
            case_evaluations=1,
            coarse_function_evaluations=missing,
            validation_function_evaluations=missing,
            targeted_function_evaluations=missing,
            polynomial_fits=missing,
            fit_holdout_evaluations=missing,
            classification_matches=missing,
            performance_verdict="historical_noncomparable",
            pnjl_residual_evaluations=missing,
            pnjl_jacobian_evaluations=missing,
            pnjl_newton_iterations=missing,
            pnjl_anchor_solves=missing,
            pnjl_branch_points=16,
            elapsed_ns=125_400_000_000,
            mean_time_ns=125_400_000_000.0,
            note="historical PNJL cold start; includes BifurcationKit compilation and is not comparable",
        ),
        (
            method="historical_multibranch_palc_warm",
            evidence_scope="historical_noncomparable",
            comparable=false,
            repetitions=1,
            case_evaluations=1,
            coarse_function_evaluations=missing,
            validation_function_evaluations=missing,
            targeted_function_evaluations=missing,
            polynomial_fits=missing,
            fit_holdout_evaluations=missing,
            classification_matches=missing,
            performance_verdict="historical_noncomparable",
            pnjl_residual_evaluations=missing,
            pnjl_jacobian_evaluations=missing,
            pnjl_newton_iterations=missing,
            pnjl_anchor_solves=missing,
            pnjl_branch_points=16,
            elapsed_ns=1_370_000_000,
            mean_time_ns=1_370_000_000.0,
            note="historical PNJL warm process; different scenario and not comparable",
        ),
    ]
    return (;
        rows,
        cascade_rows,
        temperature_sequence=sequence.rows,
        cep_extrapolation=cep_row,
        performance,
    )
end

function _sha256_file(path::String)
    return bytes2hex(open(sha256, path))
end

function _assert_finite_rows(rows, label::String)
    for (row_index, row) in enumerate(rows), field in propertynames(row)
        value = getproperty(row, field)
        if value isa AbstractFloat && !isfinite(value)
            error("$label row $row_index field $field is non-finite: $value")
        end
    end
    return rows
end

function _git_head()
    try
        return readchomp(`git -C $PROJECT_ROOT rev-parse HEAD`)
    catch
        return "unknown"
    end
end

function _git_worktree_dirty()
    try
        return !isempty(readchomp(`git -C $PROJECT_ROOT status --porcelain --untracked-files=all`))
    catch
        return true
    end
end

function _write_claim_ledger(path::String)
    rows = [
        (
            claim_id="C1",
            status="supported_synthetic",
            claim="A rho-support cascade resolves the hidden analytic S interval with at most twelve targeted evaluations.",
            evidence_file="tables/cascade_case_results.csv",
            evidence_selector="case_name=hidden_s_warning",
            boundary="Synthetic cusp evidence only; the fixed-rho solve cost has not been measured in PNJL.",
        ),
        (
            claim_id="C2",
            status="supported_by_contract",
            claim="Fast-path cases require zero targeted evaluations and zero polynomial fits.",
            evidence_file="tables/accuracy_performance_frontier.csv",
            evidence_selector="method=rho_support_fast_path",
            boundary="The timing is an in-memory analytic comparison; the zero-extra-work count is the portable result.",
        ),
        (
            claim_id="C3",
            status="supported_synthetic",
            claim="Resolved density widths provide a finite CEP estimate through a linear delta-rho-squared extrapolation.",
            evidence_file="tables/cep_extrapolation_summary.csv",
            evidence_selector="within_synthetic_gate=true",
            boundary="The critical normal-form scaling is synthetic evidence, not a PNJL CEP result.",
        ),
        (
            claim_id="C4",
            status="performance_rejected",
            claim="The original full-window multi-evidence method is suitable for the production hot path.",
            evidence_file="tables/accuracy_performance_frontier.csv",
            evidence_selector="method=multi_evidence_criticality",
            boundary="It remains a diagnostic comparator because its analytic classification overhead is excessive.",
        ),
        (
            claim_id="C5",
            status="author_check",
            claim="The rho-support cascade improves PNJL CEP accuracy and total runtime.",
            evidence_file="tables/accuracy_performance_frontier.csv",
            evidence_selector="method=rho_support_cascade",
            boundary="Requires a same-physics Actions pilot with residual/Jacobian/runner cost counters.",
        ),
    ]
    CSV.write(path, rows)
end

function _write_manifest(
        output_dir::String,
        config::CriticalityConfig,
        cascade_config::RhoSupportConfig,
        case_count::Int,
        all_cases_match_expected::Bool,
        cascade_case_count::Int,
        all_cascade_cases_match_expected::Bool,
        temperature_sequence_count::Int,
        all_temperature_sequence_match_expected::Bool,
        cep_extrapolation_within_synthetic_gate::Bool)
    tables_dir = joinpath(output_dir, "tables")
    case_path = joinpath(tables_dir, "case_results.csv")
    cascade_path = joinpath(tables_dir, "cascade_case_results.csv")
    sequence_path = joinpath(tables_dir, "temperature_sequence_results.csv")
    cep_path = joinpath(tables_dir, "cep_extrapolation_summary.csv")
    performance_path = joinpath(tables_dir, "accuracy_performance_frontier.csv")
    claim_path = joinpath(tables_dir, "claim_ledger.csv")
    for path in (case_path, cascade_path, sequence_path, cep_path, performance_path, claim_path)
        isfile(path) || error("manifest input is missing: $path")
    end
    script_path = abspath(@__FILE__)
    perf_script_path = joinpath(
        PROJECT_ROOT,
        "scripts",
        "perf",
        "pnjl_criticality_feasibility_probe.jl",
    )
    isfile(perf_script_path) || error("performance producer script is missing: $perf_script_path")
    readme_path = joinpath(output_dir, "README.md")
    generated_files = [
        (path="tables/case_results.csv", sha256=_sha256_file(case_path)),
        (path="tables/cascade_case_results.csv", sha256=_sha256_file(cascade_path)),
        (path="tables/temperature_sequence_results.csv", sha256=_sha256_file(sequence_path)),
        (path="tables/cep_extrapolation_summary.csv", sha256=_sha256_file(cep_path)),
        (path="tables/accuracy_performance_frontier.csv", sha256=_sha256_file(performance_path)),
        (path="tables/claim_ledger.csv", sha256=_sha256_file(claim_path)),
    ]
    startup_path = joinpath(tables_dir, "process_startup_timings.csv")
    if isfile(startup_path)
        push!(generated_files, (path="tables/process_startup_timings.csv", sha256=_sha256_file(startup_path)))
    end
    documentation = isfile(readme_path) ? [
        (path="README.md", sha256=_sha256_file(readme_path)),
    ] : NamedTuple[]
    phase_core_path = joinpath(PROJECT_ROOT, "src", "models", "phase", "PhaseCore.jl")
    palc_readme_path = joinpath(PROJECT_ROOT, "scripts", "analysis", "pnjl_bifurcation_spike", "README.md")
    palc_design_path = joinpath(PROJECT_ROOT, "docs", "dev", "archived", "2026-06-18_PALC连续路径求解接口设计.md")
    palc_design_rel_path = "docs/dev/archived/2026-06-18_PALC连续路径求解接口设计.md"
    manifest = (
        schema_version=1,
        artifact_kind="diagnostic_feasibility",
        production_authority=false,
        repo_head=_git_head(),
        repo_head_role="generation_base_head",
        generation_worktree_dirty=_git_worktree_dirty(),
        script_path=replace(relpath(script_path, PROJECT_ROOT), '\\' => '/'),
        script_sha256=_sha256_file(script_path),
        producer_scripts=[
            (
                path=replace(relpath(script_path, PROJECT_ROOT), '\\' => '/'),
                sha256=_sha256_file(script_path),
            ),
            (
                path=replace(relpath(perf_script_path, PROJECT_ROOT), '\\' => '/'),
                sha256=_sha256_file(perf_script_path),
            ),
        ],
        source_evidence=[
            (path="src/models/phase/PhaseCore.jl", sha256=_sha256_file(phase_core_path)),
            (path="scripts/analysis/pnjl_bifurcation_spike/README.md", sha256=_sha256_file(palc_readme_path)),
            (path=palc_design_rel_path, sha256=_sha256_file(palc_design_path)),
        ],
        config=(;
            near_slope_tol=config.near_slope_tol,
            resolved_slope_margin=config.resolved_slope_margin,
            holdout_rmse_tol=config.holdout_rmse_tol,
            holdout_max_error_tol=config.holdout_max_error_tol,
            derivative_spread_tol=config.derivative_spread_tol,
            derivative_grid_points=config.derivative_grid_points,
            minimum_validation_points=config.minimum_validation_points,
            validation_gap_ratio=config.validation_gap_ratio,
            fit_window_trim_fraction=config.fit_window_trim_fraction,
            minimum_negative_secant_run=config.minimum_negative_secant_run,
            rho_support=(;
                support_slope_tol=cascade_config.support_slope_tol,
                positive_slope_margin=cascade_config.positive_slope_margin,
                negative_slope_margin=cascade_config.negative_slope_margin,
                minimum_negative_secant_run=cascade_config.minimum_negative_secant_run,
                target_point_count=cascade_config.target_point_count,
                max_extra_points=cascade_config.max_extra_points,
                support_expansion_gaps=cascade_config.support_expansion_gaps,
                local_fit_rmse_tol=cascade_config.local_fit_rmse_tol,
                near_critical_slope_tol=cascade_config.near_critical_slope_tol,
            ),
        ),
        case_count=case_count,
        all_cases_match_expected=all_cases_match_expected,
        cascade_case_count=cascade_case_count,
        all_cascade_cases_match_expected=all_cascade_cases_match_expected,
        temperature_sequence_count=temperature_sequence_count,
        all_temperature_sequence_match_expected=all_temperature_sequence_match_expected,
        cep_extrapolation_within_synthetic_gate=cep_extrapolation_within_synthetic_gate,
        files=generated_files,
        documentation=documentation,
        limitations=[
            "synthetic cusp family only",
            "no PNJL residual, Jacobian, Newton, or runner-cost measurement",
            "historical PALC timings are explicitly noncomparable",
            "no production phase classifier or reference was modified",
        ],
    )
    open(joinpath(output_dir, "manifest.json"), "w") do io
        JSON3.pretty(io, manifest)
        println(io)
    end
    return manifest
end

function _refresh_manifest(
        output_dir::String;
        config::CriticalityConfig=CriticalityConfig(),
        cascade_config::RhoSupportConfig=RhoSupportConfig())
    case_path = joinpath(output_dir, "tables", "case_results.csv")
    rows = collect(CSV.File(case_path))
    isempty(rows) && error("case_results.csv is empty")
    all_match = all(Bool(row.matches_expected) for row in rows)
    all_match || error("synthetic feasibility suite did not match expected states")
    cascade_rows = collect(CSV.File(joinpath(
        output_dir,
        "tables",
        "cascade_case_results.csv",
    )))
    sequence_rows = collect(CSV.File(joinpath(
        output_dir,
        "tables",
        "temperature_sequence_results.csv",
    )))
    cep_rows = collect(CSV.File(joinpath(
        output_dir,
        "tables",
        "cep_extrapolation_summary.csv",
    )))
    isempty(cascade_rows) && error("cascade_case_results.csv is empty")
    isempty(sequence_rows) && error("temperature_sequence_results.csv is empty")
    length(cep_rows) == 1 || error("cep_extrapolation_summary.csv must contain one row")
    all_cascade_match = all(Bool(row.matches_expected) for row in cascade_rows)
    all_sequence_match = all(Bool(row.matches_expected) for row in sequence_rows)
    cep_gate = Bool(only(cep_rows).within_synthetic_gate)
    all_cascade_match || error("rho-support cascade cases did not match expected states")
    all_sequence_match || error("temperature sequence did not match expected states")
    cep_gate || error("density-width CEP extrapolation missed the synthetic gate")
    return _write_manifest(
        output_dir,
        config,
        cascade_config,
        length(rows),
        all_match,
        length(cascade_rows),
        all_cascade_match,
        length(sequence_rows),
        all_sequence_match,
        cep_gate,
    )
end

function write_evidence_package(
        output_dir::String;
        cases::Vector{SyntheticCase}=default_cases(),
        config::CriticalityConfig=CriticalityConfig(),
        cascade_config::RhoSupportConfig=RhoSupportConfig(),
        repetitions::Int=100)
    suite = run_suite(;
        cases=cases,
        config=config,
        cascade_config=cascade_config,
        repetitions=repetitions,
    )
    tables_dir = joinpath(output_dir, "tables")
    mkpath(tables_dir)
    case_path = joinpath(tables_dir, "case_results.csv")
    cascade_path = joinpath(tables_dir, "cascade_case_results.csv")
    sequence_path = joinpath(tables_dir, "temperature_sequence_results.csv")
    cep_path = joinpath(tables_dir, "cep_extrapolation_summary.csv")
    performance_path = joinpath(tables_dir, "accuracy_performance_frontier.csv")
    claim_path = joinpath(tables_dir, "claim_ledger.csv")
    startup_path = joinpath(tables_dir, "process_startup_timings.csv")
    isfile(startup_path) && rm(startup_path)
    _assert_finite_rows(suite.rows, "case results")
    _assert_finite_rows(suite.cascade_rows, "cascade case results")
    _assert_finite_rows(suite.temperature_sequence, "temperature sequence results")
    _assert_finite_rows([suite.cep_extrapolation], "CEP extrapolation summary")
    _assert_finite_rows(suite.performance, "accuracy/performance frontier")
    CSV.write(case_path, suite.rows)
    CSV.write(cascade_path, suite.cascade_rows)
    CSV.write(sequence_path, suite.temperature_sequence)
    CSV.write(cep_path, [suite.cep_extrapolation])
    CSV.write(performance_path, suite.performance)
    _write_claim_ledger(claim_path)

    all_match = all(row.matches_expected for row in suite.rows)
    all_match || error("synthetic feasibility suite did not match expected states")
    all_cascade_match = all(row.matches_expected for row in suite.cascade_rows)
    all_sequence_match = all(row.matches_expected for row in suite.temperature_sequence)
    all_cascade_match || error("rho-support cascade cases did not match expected states")
    all_sequence_match || error("temperature sequence did not match expected states")
    suite.cep_extrapolation.within_synthetic_gate || error(
        "density-width CEP extrapolation missed the synthetic gate",
    )
    _write_manifest(
        output_dir,
        config,
        cascade_config,
        length(suite.rows),
        all_match,
        length(suite.cascade_rows),
        all_cascade_match,
        length(suite.temperature_sequence),
        all_sequence_match,
        suite.cep_extrapolation.within_synthetic_gate,
    )
    return suite
end

function _parse_cli(args)
    output_dir = joinpath(
        PROJECT_ROOT,
        "docs",
        "analysis",
        "pnjl",
        "algorithmic_feasibility",
        "criticality_feasibility_v1",
    )
    repetitions = 100
    manifest_only = false
    for arg in args
        if startswith(arg, "--output-dir=")
            output_dir = abspath(split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--repetitions=")
            repetitions = parse(Int, split(arg, "="; limit=2)[2])
        elseif arg == "--manifest-only"
            manifest_only = true
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. scripts/analysis/pnjl_criticality_feasibility.jl [--output-dir=PATH] [--repetitions=N] [--manifest-only]")
            return nothing
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
    end
    return (; output_dir, repetitions, manifest_only)
end

function main(args=ARGS)
    options = _parse_cli(args)
    options === nothing && return 0
    if options.manifest_only
        manifest = _refresh_manifest(options.output_dir)
        println("criticality feasibility manifest refreshed: ", options.output_dir)
        println("synthetic cases matched: ", manifest.case_count, "/", manifest.case_count)
    else
        suite = write_evidence_package(options.output_dir; repetitions=options.repetitions)
        println("criticality feasibility output: ", options.output_dir)
        println("synthetic cases matched: ", count(row -> row.matches_expected, suite.rows), "/", length(suite.rows))
    end
    return 0
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    exit(PNJLCriticalityFeasibility.main())
end
