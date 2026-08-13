const EPS_SLOPE = 0.0
const DEFAULT_AREA_TOL = 1e-4
const DEFAULT_MAXWELL_SOLVER_TOL_FACTOR = 0.1
const DEFAULT_MAX_ITER = 60
const DEFAULT_CANDIDATE_STEPS = 64
const MAX_CANDIDATE_STEPS = 1024
const BRACKET_SHRINK_REL = 1e-3
const BRACKET_SHRINK_ABS = 1e-3
const MAXWELL_CANDIDATE_POLICIES = (
    :unique_three_crossing_topology_v1,
    :unique_three_crossing_sign_change_v2,
)

"""Derive an internal Maxwell stopping tolerance from active acceptance gates.

The public `maxwell_construction` default remains `DEFAULT_AREA_TOL` for
backward compatibility.  Production callers use this helper so the internal
bisection is strictly tighter than every active area acceptance tolerance;
the outer rho/temperature geometry gate remains a separate certificate.
"""
function _derive_maxwell_solver_tol(active_tolerances::AbstractVector{<:Real};
        factor::Real=DEFAULT_MAXWELL_SOLVER_TOL_FACTOR)
    isempty(active_tolerances) && throw(ArgumentError("at least one active Maxwell tolerance is required"))
    isfinite(factor) && factor > 0 || throw(ArgumentError(
        "Maxwell solver tolerance factor must be finite and positive, got $(factor)",
    ))
    values = Float64[]
    for tolerance in active_tolerances
        value = Float64(tolerance)
        isfinite(value) && value > 0 || throw(ArgumentError(
            "active Maxwell acceptance tolerances must be finite and positive, got $(tolerance)",
        ))
        push!(values, value)
    end
    derived = Float64(factor) * minimum(values)
    isfinite(derived) && derived > 0 || throw(ArgumentError(
        "derived Maxwell solver tolerance is invalid: $(derived)",
    ))
    return derived
end

@inline function _maxwell_result_satisfies_tol(result, tolerance::Real)
    result.converged && result.area_residual !== nothing &&
        isfinite(Float64(result.area_residual)) &&
        Float64(result.area_residual) <= Float64(tolerance)
end

struct SShapeResult
    has_s_shape::Bool
    mu_spinodal_hadron::Union{Nothing, Float64}
    mu_spinodal_quark::Union{Nothing, Float64}
    rho_spinodal_hadron::Union{Nothing, Float64}
    rho_spinodal_quark::Union{Nothing, Float64}
    derivative_sign_changes::Int
end

SShapeResult() = SShapeResult(false, nothing, nothing, nothing, nothing, 0)

struct MaxwellResult
    converged::Bool
    mu_transition::Union{Nothing, Float64}
    rho_hadron::Union{Nothing, Float64}
    rho_quark::Union{Nothing, Float64}
    area_residual::Union{Nothing, Float64}
    iterations::Int
    details::Dict{Symbol, Any}
end

MaxwellResult() = MaxwellResult(false, nothing, nothing, nothing, nothing, 0, Dict{Symbol, Any}())

@inline _slope_sign(value) = abs(value) <= EPS_SLOPE ? 0 : (value > 0 ? 1 : -1)

function _sort_curve_by_rho(mu_vals::AbstractVector, rho_vals::AbstractVector)
    n = min(length(mu_vals), length(rho_vals))
    order = sortperm(rho_vals[1:n])
    return Float64.(mu_vals[order]), Float64.(rho_vals[order])
end

function _refine_extremum(mu_sorted::Vector{Float64}, rho_sorted::Vector{Float64}, idx::Int; is_maximum::Bool=true)
    n = length(rho_sorted)
    if idx < 2 || idx >= n
        return (rho_sorted[idx] + rho_sorted[idx + 1]) / 2,
               (mu_sorted[idx] + mu_sorted[idx + 1]) / 2
    end

    i1, i2, i3 = idx - 1, idx, idx + 1
    if i3 > n
        i1, i2, i3 = idx - 2, idx - 1, idx
    end
    if i1 < 1
        i1, i2, i3 = 1, 2, 3
    end

    r1, r2, r3 = rho_sorted[i1], rho_sorted[i2], rho_sorted[i3]
    m1, m2, m3 = mu_sorted[i1], mu_sorted[i2], mu_sorted[i3]
    denom = (r1 - r2) * (r1 - r3) * (r2 - r3)

    if abs(denom) < 1e-15
        return (rho_sorted[idx] + rho_sorted[idx + 1]) / 2,
               (mu_sorted[idx] + mu_sorted[idx + 1]) / 2
    end

    a = (r3 * (m2 - m1) + r2 * (m1 - m3) + r1 * (m3 - m2)) / denom
    b = (r3^2 * (m1 - m2) + r2^2 * (m3 - m1) + r1^2 * (m2 - m3)) / denom
    c = (r2 * r3 * (r2 - r3) * m1 + r3 * r1 * (r3 - r1) * m2 + r1 * r2 * (r1 - r2) * m3) / denom

    if (is_maximum && a >= 0) || (!is_maximum && a <= 0)
        return (rho_sorted[idx] + rho_sorted[idx + 1]) / 2,
               (mu_sorted[idx] + mu_sorted[idx + 1]) / 2
    end

    rho_ext = -b / (2 * a)
    rho_min = min(r1, r2, r3)
    rho_max = max(r1, r2, r3)
    if rho_ext < rho_min || rho_ext > rho_max
        return (rho_sorted[idx] + rho_sorted[idx + 1]) / 2,
               (mu_sorted[idx] + mu_sorted[idx + 1]) / 2
    end

    mu_ext = a * rho_ext^2 + b * rho_ext + c
    return rho_ext, mu_ext
end

function detect_s_shape(mu_vals::AbstractVector, rho_vals::AbstractVector; eps::Real=EPS_SLOPE, min_points::Int=5)
    n = min(length(mu_vals), length(rho_vals))
    n >= min_points || return SShapeResult()

    mu_sorted, rho_sorted = _sort_curve_by_rho(mu_vals, rho_vals)

    slopes = Float64[]
    slope_indices = Int[]
    for i in 1:(length(rho_sorted) - 1)
        drho = rho_sorted[i + 1] - rho_sorted[i]
        abs(drho) <= eps && continue
        push!(slopes, (mu_sorted[i + 1] - mu_sorted[i]) / drho)
        push!(slope_indices, i)
    end
    isempty(slopes) && return SShapeResult()

    signs = Int[]
    sign_to_slope_idx = Int[]
    last_sign = 0
    for (i, slope) in enumerate(slopes)
        sign = _slope_sign(slope)
        sign == 0 && continue
        if sign != last_sign
            push!(signs, sign)
            push!(sign_to_slope_idx, i)
            last_sign = sign
        end
    end

    if length(signs) < 3
        return SShapeResult(false, nothing, nothing, nothing, nothing, max(0, length(signs) - 1))
    end

    max_idx = nothing
    min_idx = nothing
    for i in 1:(length(signs) - 1)
        if signs[i] == 1 && signs[i + 1] == -1 && max_idx === nothing
            max_idx = i + 1
        elseif signs[i] == -1 && signs[i + 1] == 1 && max_idx !== nothing && min_idx === nothing
            min_idx = i + 1
        end
    end

    if max_idx === nothing || min_idx === nothing
        return SShapeResult(false, nothing, nothing, nothing, nothing, length(signs) - 1)
    end

    slope_idx_max = sign_to_slope_idx[max_idx]
    slope_idx_min = sign_to_slope_idx[min_idx]
    data_idx_max = slope_indices[slope_idx_max]
    data_idx_min = slope_indices[slope_idx_min]

    rho_hadron, mu_hadron = _refine_extremum(mu_sorted, rho_sorted, data_idx_max; is_maximum=true)
    rho_quark, mu_quark = _refine_extremum(mu_sorted, rho_sorted, data_idx_min; is_maximum=false)

    return SShapeResult(true, mu_hadron, mu_quark, rho_hadron, rho_quark, length(signs) - 1)
end

function _prepare_curve(mu_vals::AbstractVector, rho_vals::AbstractVector)
    n = min(length(mu_vals), length(rho_vals))
    pairs = Vector{Tuple{Float64, Float64}}()
    sizehint!(pairs, n)
    for i in 1:n
        mu = Float64(mu_vals[i])
        rho = Float64(rho_vals[i])
        (isfinite(mu) && isfinite(rho)) || continue
        push!(pairs, (rho, mu))
    end
    sort!(pairs; by=first)
    return Float64[first(p) for p in pairs], Float64[last(p) for p in pairs]
end

function _mu_bracket(hint::SShapeResult)
    mu_max = hint.mu_spinodal_hadron
    mu_min = hint.mu_spinodal_quark
    (mu_max === nothing || mu_min === nothing) && return nothing
    return mu_min < mu_max ? (mu_min, mu_max) : (mu_max, mu_min)
end

function _shrink_bracket(mu_lo::Float64, mu_hi::Float64)
    width = mu_hi - mu_lo
    width > 0 || return nothing
    δ = max(width * BRACKET_SHRINK_REL, BRACKET_SHRINK_ABS)
    if mu_lo + δ >= mu_hi - δ
        δ = width / 4
        mu_lo + δ < mu_hi - δ || return nothing
    end
    return (mu_lo + δ, mu_hi - δ)
end

"""Enumerate all crossings of a horizontal chemical-potential line.

The old implementation intentionally returned only the first and last
crossing.  That is unsafe for Maxwell construction: a two-crossing endpoint
topology can otherwise be mistaken for a valid outer pair.  Exact vertices
are included once and linear crossings are interpolated between adjacent
samples.  A plateau is represented by its endpoints and is therefore still
rejected by the strict three-crossing contract below.
"""
function _find_intersections(mu0::Float64, rho_vals::Vector{Float64}, mu_vals::Vector{Float64};
        atol::Real=1e-9)
    n = min(length(rho_vals), length(mu_vals))
    n >= 2 || return Float64[]
    crossings = Float64[]
    for i in 1:(n - 1)
        r1, r2 = rho_vals[i], rho_vals[i + 1]
        f1, f2 = mu_vals[i] - mu0, mu_vals[i + 1] - mu0
        abs(f1) <= atol && push!(crossings, r1)
        if f1 * f2 < 0.0
            α = f1 / (f1 - f2)
            push!(crossings, r1 + α * (r2 - r1))
        end
        abs(f2) <= atol && push!(crossings, r2)
    end
    sort!(crossings)
    unique_crossings = Float64[]
    merge_atol = max(10.0 * Float64(atol), 32eps(Float64))
    for value in crossings
        isempty(unique_crossings) || abs(value - last(unique_crossings)) > merge_atol || continue
        push!(unique_crossings, value)
    end
    return unique_crossings
end

"""Compatibility helper returning the outermost pair of all crossings."""
function _find_outer_intersections(mu0::Float64, rho_vals::Vector{Float64}, mu_vals::Vector{Float64};
        atol::Real=1e-9)
    crossings = _find_intersections(mu0, rho_vals, mu_vals; atol=atol)
    length(crossings) >= 2 || return nothing, nothing
    return first(crossings), last(crossings)
end

@inline function _interp(r1::Float64, r2::Float64, m1::Float64, m2::Float64, target::Float64)
    r2 == r1 && return m1
    return m1 + (target - r1) / (r2 - r1) * (m2 - m1)
end

function _integrate_difference(rho_vals::Vector{Float64}, mu_vals::Vector{Float64}, rho_left::Float64, rho_right::Float64, mu0::Float64)
    total = 0.0
    for i in 1:(length(rho_vals) - 1)
        r1, r2 = rho_vals[i], rho_vals[i + 1]
        (r2 <= rho_left || r1 >= rho_right) && continue
        left, right = max(r1, rho_left), min(r2, rho_right)
        right <= left && continue
        mu_left = _interp(r1, r2, mu_vals[i], mu_vals[i + 1], left)
        mu_right = _interp(r1, r2, mu_vals[i], mu_vals[i + 1], right)
        total += 0.5 * ((mu_left - mu0) + (mu_right - mu0)) * (right - left)
    end
    return total
end

function _area_probe(mu0::Float64, rho_vals::Vector{Float64}, mu_vals::Vector{Float64};
        endpoint_atol::Real=1e-9)
    crossings = _find_intersections(mu0, rho_vals, mu_vals; atol=endpoint_atol)
    if length(crossings) != 3
        return (valid=false, area=nothing, crossings=crossings,
            endpoint_dependent=false, reason="crossing_count_not_three")
    end
    rho_left, rho_right = first(crossings), last(crossings)
    rho_right > rho_left || return (valid=false, area=nothing, crossings=crossings,
        endpoint_dependent=false, reason="degenerate_crossings")
    endpoint_cell = length(rho_vals) >= 2 ? abs(rho_vals[2] - rho_vals[1]) : Inf
    endpoint_dependent = rho_vals[1] == 0.0 && rho_left <= endpoint_cell + 32eps(Float64)
    return (valid=true,
        area=_integrate_difference(rho_vals, mu_vals, rho_left, rho_right, mu0),
        crossings=crossings,
        endpoint_dependent=endpoint_dependent,
        reason="three_crossings")
end

"""Compatibility scalar area helper.

Only a strict three-crossing topology has an area.  Returning `nothing` for a
topology gap is important because sign-change scans must reset across it.
"""
function _area_difference(mu0::Float64, rho_vals::Vector{Float64}, mu_vals::Vector{Float64})
    probe = _area_probe(mu0, rho_vals, mu_vals)
    return probe.valid ? probe.area : nothing
end

function _bisection_solve(rho_vals::Vector{Float64}, mu_vals::Vector{Float64},
        mu_a::Float64, mu_b::Float64, area_a::Real, area_b::Real,
        tol_area::Real, max_iter::Int)
    max_iter > 0 || return (converged=false, reason="max_iter_nonpositive",
        mu=nothing, area=nothing, iterations=0, crossings=Float64[], endpoint_dependent=false)
    if mu_a == mu_b
        probe = _area_probe(mu_a, rho_vals, mu_vals)
        # A zero-width bracket is a diagnostic probe only.  It is not a
        # Maxwell root because no area sign change was established.
        return (converged=false, reason="grid_hit_not_candidate", mu=mu_a,
            area=probe.area, iterations=0, crossings=probe.crossings,
            endpoint_dependent=probe.endpoint_dependent)
    end
    a, b, fa, fb = mu_a, mu_b, Float64(area_a), Float64(area_b)
    fa * fb <= 0.0 || return (converged=false, reason="invalid_area_bracket",
        mu=nothing, area=nothing, iterations=0, crossings=Float64[], endpoint_dependent=false)
    best = nothing
    for iter in 1:max_iter
        mid = 0.5 * (a + b)
        probe = _area_probe(mid, rho_vals, mu_vals)
        if !probe.valid
            return (converged=false, reason="topology_changed_inside_bisection",
                mu=mid, area=nothing, iterations=iter, crossings=probe.crossings,
                endpoint_dependent=false)
        end
        best = probe
        if abs(probe.area) <= tol_area
            return (converged=true, reason="ok", mu=mid, area=probe.area,
                iterations=iter, crossings=probe.crossings,
                endpoint_dependent=probe.endpoint_dependent)
        end
        if fa * probe.area < 0.0
            b, fb = mid, probe.area
        else
            a, fa = mid, probe.area
        end
    end
    # Exhausting max_iter is a numerical failure, even if a finite best probe
    # is available.  Callers must not publish it as a converged certificate.
    return (converged=false, reason="solver_tolerance_not_met",
        mu=0.5 * (a + b), area=best === nothing ? nothing : best.area,
        iterations=max_iter, crossings=best === nothing ? Float64[] : best.crossings,
        endpoint_dependent=best === nothing ? false : best.endpoint_dependent)
end

function _candidate_roots(rho_vals::Vector{Float64}, mu_vals::Vector{Float64},
        mu_lo::Float64, mu_hi::Float64, steps::Int, tol_area::Real, max_iter::Int;
        candidate_policy::Symbol=:unique_three_crossing_sign_change_v2)
    candidate_policy in MAXWELL_CANDIDATE_POLICIES || throw(ArgumentError(
        "unsupported Maxwell candidate policy $(candidate_policy); expected one of $(MAXWELL_CANDIDATE_POLICIES)",
    ))
    include_grid_hit_candidates = candidate_policy === :unique_three_crossing_topology_v1
    attempt = max(steps, 3)
    valid_intervals = NamedTuple[]
    roots = NamedTuple[]
    failed_candidates = NamedTuple[]
    near_zero_grid_hits = NamedTuple[]
    while attempt <= MAX_CANDIDATE_STEPS
        previous = nothing
        interval_start = nothing
        interval_last = nothing
        best_near_zero = nothing
        interval_has_root = false
        for mu0 in range(mu_lo, mu_hi; length=attempt)
            probe = _area_probe(Float64(mu0), rho_vals, mu_vals)
            if !probe.valid
                if interval_start !== nothing
                    push!(valid_intervals, (mu_low=interval_start, mu_high=interval_last,
                        status="three_crossing_interval"))
                    include_grid_hit_candidates && !interval_has_root && best_near_zero !== nothing &&
                        push!(roots, merge(best_near_zero,
                            (bracket=(best_near_zero.mu, best_near_zero.mu),)))
                end
                previous = nothing
                interval_start = nothing
                interval_last = nothing
                best_near_zero = nothing
                interval_has_root = false
                continue
            end
            interval_start === nothing && (interval_start = Float64(mu0))
            interval_last = Float64(mu0)
            if abs(probe.area) <= tol_area
                hit = (mu=Float64(mu0), area=probe.area, crossings=probe.crossings,
                    endpoint_dependent=probe.endpoint_dependent, reason="near_zero_grid_probe",
                    bracket=(Float64(mu0), Float64(mu0)))
                push!(near_zero_grid_hits, hit)
                if include_grid_hit_candidates &&
                   (best_near_zero === nothing || abs(probe.area) < abs(best_near_zero.area))
                    best_near_zero = merge(probe, (converged=true, reason="grid_hit",
                        mu=Float64(mu0), iterations=0))
                end
            end
            if previous !== nothing && previous.area * probe.area < 0.0
                solved = _bisection_solve(rho_vals, mu_vals,
                    previous.mu, Float64(mu0), previous.area, probe.area,
                    tol_area, max_iter)
                if solved.converged
                    push!(roots, merge(solved, (bracket=(previous.mu, Float64(mu0)),)))
                    interval_has_root = true
                else
                    push!(failed_candidates, merge(solved, (bracket=(previous.mu, Float64(mu0)),)))
                end
            end
            previous = merge(probe, (mu=Float64(mu0),))
        end
        if interval_start !== nothing
            push!(valid_intervals,
                (mu_low=interval_start, mu_high=interval_last, status="three_crossing_interval"))
            include_grid_hit_candidates && !interval_has_root && best_near_zero !== nothing &&
                push!(roots, merge(best_near_zero,
                    (bracket=(best_near_zero.mu, best_near_zero.mu),)))
        end
        # Increasing scan density can expose a topology interval that was
        # skipped at a coarser level.  Duplicate roots from adjacent scans
        # are merged below by overlapping brackets.
        attempt *= 2
    end
    sort!(roots; by=root -> root.mu)
    unique_roots = NamedTuple[]
    for root in roots
        duplicate = false
        for previous in unique_roots
            lo_a, hi_a = root.bracket
            lo_b, hi_b = previous.bracket
            if max(lo_a, lo_b) <= min(hi_a, hi_b) + 32eps(Float64) ||
               abs(root.mu - previous.mu) <= 32eps(Float64) * max(1.0, abs(root.mu))
                duplicate = true
                break
            end
        end
        duplicate || push!(unique_roots, root)
    end
    # The same grid probe can be seen at multiple scan resolutions.  Keep one
    # deterministic diagnostic record per zero-width bracket.
    unique_grid_hits = NamedTuple[]
    for hit in sort(near_zero_grid_hits; by=item -> (item.mu, abs(item.area)))
        duplicate = any(existing -> existing.mu == hit.mu, unique_grid_hits)
        duplicate || push!(unique_grid_hits, hit)
    end
    return (roots=unique_roots, failed_candidates=failed_candidates,
        valid_intervals=valid_intervals, near_zero_grid_hits=unique_grid_hits)
end

function _failure(reason::AbstractString; kwargs...)
    details = Dict{Symbol, Any}(:reason => reason)
    for (k, v) in kwargs
        details[k] = v
    end
    return MaxwellResult(false, nothing, nothing, nothing, nothing, 0, details)
end

function maxwell_construction(mu_vals::AbstractVector, rho_vals::AbstractVector;
        min_samples::Int=12, detect_min_points::Int=6, detect_eps::Real=1e-6,
        candidate_steps::Int=DEFAULT_CANDIDATE_STEPS,
        max_iter::Int=DEFAULT_MAX_ITER, tol_area::Real=DEFAULT_AREA_TOL,
        spinodal_hint::Union{Nothing, SShapeResult}=nothing,
        candidate_policy::Symbol=:unique_three_crossing_sign_change_v2)

    candidate_policy in MAXWELL_CANDIDATE_POLICIES || throw(ArgumentError(
        "unsupported Maxwell candidate policy $(candidate_policy); expected one of $(MAXWELL_CANDIDATE_POLICIES)",
    ))

    rho_sorted, mu_sorted = _prepare_curve(mu_vals, rho_vals)
    length(rho_sorted) < min_samples && return _failure("insufficient_points"; count=length(rho_sorted))

    hint = isnothing(spinodal_hint) ?
        detect_s_shape(mu_vals, rho_vals; eps=detect_eps, min_points=detect_min_points) : spinodal_hint
    hint.has_s_shape || return _failure("no_s_shape")

    mu_bracket = _mu_bracket(hint)
    mu_bracket === nothing && return _failure("invalid_mu_bracket")
    mu_lo, mu_hi = mu_bracket

    tightened = _shrink_bracket(mu_lo, mu_hi)
    tightened === nothing && return _failure("degenerate_bracket"; bracket=(mu_lo, mu_hi))
    mu_lo, mu_hi = tightened

    candidates = _candidate_roots(rho_sorted, mu_sorted, mu_lo, mu_hi,
        candidate_steps, tol_area, max_iter; candidate_policy=candidate_policy)
    roots = candidates.roots
    observed_crossing_count = if !isempty(roots)
        length(first(roots).crossings)
    elseif !isempty(candidates.failed_candidates)
        length(first(candidates.failed_candidates).crossings)
    elseif !isempty(candidates.near_zero_grid_hits)
        length(first(candidates.near_zero_grid_hits).crossings)
    else
        0
    end
    base_details = Dict{Symbol, Any}(
        :mu_search_bracket => (mu_lo, mu_hi),
        :spinodal_hint => (hint.mu_spinodal_hadron, hint.mu_spinodal_quark),
        :candidate_count => length(roots),
        :crossing_count => observed_crossing_count,
        :valid_intervals => candidates.valid_intervals,
        :candidate_policy => candidate_policy,
        :near_zero_grid_hits => candidates.near_zero_grid_hits,
        :near_zero_grid_probe_count => length(candidates.near_zero_grid_hits),
    )
    if isempty(roots)
        if !isempty(candidates.failed_candidates)
            return _failure("bisection_failed"; base_details...,
                bisection_failures=candidates.failed_candidates,
                failure_reason=first(candidates.failed_candidates).reason)
        end
        return _failure("no_sign_change"; base_details...,
            failure_reason="no_three_crossing_sign_change")
    end
    if length(roots) > 1
        base_details[:candidate_roots] = [(
            mu=root.mu, area=abs(root.area), crossings=root.crossings,
            converged=root.converged, endpoint_dependent=root.endpoint_dependent,
            bracket=root.bracket) for root in roots]
        return _failure("multiple_maxwell_candidates"; base_details...,
            failure_reason="multiple_maxwell_candidates")
    end

    root = first(roots)
    root.converged || return _failure("bisection_failed"; base_details...,
        bracket=root.bracket, bisection_reason=root.reason)
    length(root.crossings) == 3 || return _failure("crossing_count_not_three";
        base_details..., crossing_count=length(root.crossings))
    rho_left, rho_right = first(root.crossings), last(root.crossings)
    details = merge(base_details, Dict{Symbol, Any}(
        :mu_bracket => root.bracket,
        :rho_interval => (rho_left, rho_right),
        :crossings => root.crossings,
        :endpoint_dependent => root.endpoint_dependent,
        :endpoint_bracket => root.endpoint_dependent ?
            (rho_sorted[1], min(rho_sorted[2], rho_right)) : nothing,
        :bisection_reason => root.reason,
    ))
    return MaxwellResult(true, root.mu, rho_left, rho_right,
        abs(root.area), root.iterations, details)
end
