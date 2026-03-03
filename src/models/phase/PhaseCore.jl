const EPS_SLOPE = 0.0
const DEFAULT_AREA_TOL = 1e-4
const DEFAULT_MAX_ITER = 60
const DEFAULT_CANDIDATE_STEPS = 64
const MAX_CANDIDATE_STEPS = 1024
const BRACKET_SHRINK_REL = 1e-3
const BRACKET_SHRINK_ABS = 1e-3

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

function _find_outer_intersections(mu0::Float64, rho_vals::Vector{Float64}, mu_vals::Vector{Float64}; atol::Real=1e-9)
    n = length(rho_vals)
    n >= 2 || return nothing, nothing
    left, right = nothing, nothing
    for i in 1:(n - 1)
        r1, r2 = rho_vals[i], rho_vals[i + 1]
        f1, f2 = mu_vals[i] - mu0, mu_vals[i + 1] - mu0
        if abs(f1) < atol
            left === nothing && (left = r1)
            right = r1
        end
        f1 == f2 && continue
        if f1 * f2 < 0
            α = f1 / (f1 - f2)
            crossing = r1 + α * (r2 - r1)
            left === nothing && (left = crossing)
            right = crossing
        elseif abs(f2) < atol
            left === nothing && (left = r2)
            right = r2
        end
    end
    return left, right
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

function _area_difference(mu0::Float64, rho_vals::Vector{Float64}, mu_vals::Vector{Float64})
    rho_left, rho_right = _find_outer_intersections(mu0, rho_vals, mu_vals)
    (rho_left === nothing || rho_right === nothing || rho_right - rho_left <= 1e-9) && return nothing
    return _integrate_difference(rho_vals, mu_vals, rho_left, rho_right, mu0)
end

function _find_mu_bracket(rho_vals::Vector{Float64}, mu_vals::Vector{Float64}, mu_lo::Float64, mu_hi::Float64, steps::Int, tol_area::Real)
    attempt = max(steps, 3)
    while attempt <= MAX_CANDIDATE_STEPS
        prev_mu, prev_area = nothing, nothing
        for mu0 in range(mu_lo, mu_hi; length=attempt)
            area = _area_difference(mu0, rho_vals, mu_vals)
            area === nothing && continue
            abs(area) <= tol_area && return (mu0, mu0, area, area)
            if prev_area !== nothing && area * prev_area < 0
                return (prev_mu, mu0, prev_area, area)
            end
            prev_mu, prev_area = mu0, area
        end
        attempt *= 2
    end
    return nothing
end

function _bisection_solve(rho_vals::Vector{Float64}, mu_vals::Vector{Float64}, mu_a::Float64, mu_b::Float64, area_a::Real, area_b::Real, tol_area::Real, max_iter::Int)
    mu_a == mu_b && return mu_a, area_a, 0
    a, b, fa, fb = mu_a, mu_b, area_a, area_b
    fa * fb <= 0 || return nothing, nothing, 0

    for iter in 1:max_iter
        mid = 0.5 * (a + b)
        area_mid = _area_difference(mid, rho_vals, mu_vals)
        area_mid === nothing && return nothing, nothing, iter
        abs(area_mid) <= tol_area && return mid, area_mid, iter
        if fa * area_mid < 0
            b, fb = mid, area_mid
        else
            a, fa = mid, area_mid
        end
    end

    mid = 0.5 * (a + b)
    area_mid = _area_difference(mid, rho_vals, mu_vals)
    area_mid === nothing && return nothing, nothing, max_iter
    return mid, area_mid, max_iter
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
        spinodal_hint::Union{Nothing, SShapeResult}=nothing)

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

    bracket = _find_mu_bracket(rho_sorted, mu_sorted, mu_lo, mu_hi, candidate_steps, tol_area)
    bracket === nothing && return _failure("no_sign_change"; bracket=(mu_lo, mu_hi))
    mu_a, mu_b, area_a, area_b = bracket

    mu_root, area_root, iterations = _bisection_solve(rho_sorted, mu_sorted, mu_a, mu_b, area_a, area_b, tol_area, max_iter)
    mu_root === nothing && return _failure("bisection_failed"; bracket=(mu_a, mu_b))

    rho_left, rho_right = _find_outer_intersections(mu_root, rho_sorted, mu_sorted)
    (rho_left === nothing || rho_right === nothing || !(rho_left < rho_right)) &&
        return _failure("no_crossings"; mu_transition=mu_root, bracket=(mu_a, mu_b))

    details = Dict(
        :mu_bracket => (mu_a, mu_b),
        :rho_interval => (rho_left, rho_right),
        :spinodal_hint => (hint.mu_spinodal_hadron, hint.mu_spinodal_quark),
    )
    return MaxwellResult(true, mu_root, rho_left, rho_right, abs(area_root), iterations, details)
end
