module MaxwellRhoMu

using ..PhaseTransition: SShapeResult, detect_s_shape, group_curves_by_temperature
import ..CEPFinder: build_curves

export MaxwellResult, maxwell_rho_mu, build_phase_boundary, phase_boundary_from_rows

const DEFAULT_AREA_TOL = 1e-4
const DEFAULT_MAX_ITER = 60
const DEFAULT_CANDIDATE_STEPS = 64
const MAX_CANDIDATE_STEPS = 1024
const BRACKET_SHRINK_REL = 1e-3
const BRACKET_SHRINK_ABS = 1e-3

struct MaxwellResult
    converged::Bool
    mu_coex_MeV::Union{Nothing, Float64}
    rho_gas::Union{Nothing, Float64}
    rho_liquid::Union{Nothing, Float64}
    area_residual::Union{Nothing, Float64}
    iterations::Int
    details::Dict{Symbol, Any}
end

MaxwellResult() = MaxwellResult(false, nothing, nothing, nothing, nothing, 0, Dict{Symbol, Any}())

function maxwell_rho_mu(mu_vals::AbstractVector, rho_vals::AbstractVector;
        min_samples::Int=12, detect_min_points::Int=6, detect_eps::Real=1e-6,
        candidate_steps::Int=DEFAULT_CANDIDATE_STEPS,
        max_iter::Int=DEFAULT_MAX_ITER, tol_area::Real=DEFAULT_AREA_TOL,
        spinodal_hint::Union{Nothing, SShapeResult}=nothing)

    rho_sorted, mu_sorted = _prepare_curve(mu_vals, rho_vals)
    if length(rho_sorted) < min_samples
        return _failure("insufficient_points"; count=length(rho_sorted))
    end

    hint = isnothing(spinodal_hint) ?
        detect_s_shape(mu_vals, rho_vals; eps=detect_eps, min_points=detect_min_points) :
        spinodal_hint
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

    mu_root, area_root, iterations = _bisection_solve(rho_sorted, mu_sorted,
        mu_a, mu_b, area_a, area_b, tol_area, max_iter)
    mu_root === nothing && return _failure("bisection_failed"; bracket=(mu_a, mu_b))

    rho_left, rho_right = _find_outer_intersections(mu_root, rho_sorted, mu_sorted)
    if rho_left === nothing || rho_right === nothing || !(rho_left < rho_right)
        return _failure("no_crossings"; mu_coex=mu_root, bracket=(mu_a, mu_b))
    end

    details = Dict(
        :mu_bracket => (mu_a, mu_b),
        :rho_interval => (rho_left, rho_right),
        :spinodal_hint => (hint.mu_spinodal_low, hint.mu_spinodal_high),
    )
    return MaxwellResult(true, mu_root, rho_left, rho_right, abs(area_root), iterations, details)
end

function build_phase_boundary(curves::Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}; kwargs...)
    results = Dict{Float64, MaxwellResult}()
    for T in sort(collect(keys(curves)))
        mu_vals, rho_vals = curves[T]
        results[T] = maxwell_rho_mu(mu_vals, rho_vals; kwargs...)
    end
    return results
end

function phase_boundary_from_rows(rows; xi::Real=0.0, tol::Real=1e-6, kwargs...)
    grouped = group_curves_by_temperature(rows; xi=xi, tol=tol)
    curves = build_curves(grouped)
    return build_phase_boundary(curves; kwargs...)
end

# ---------------- internal helpers ----------------

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
    sort!(pairs; by = first)
    rho_sorted = Float64[first(p) for p in pairs]
    mu_sorted = Float64[last(p) for p in pairs]
    return rho_sorted, mu_sorted
end

function _mu_bracket(hint::SShapeResult)
    low = hint.mu_spinodal_low
    high = hint.mu_spinodal_high
    if low === nothing || high === nothing
        return nothing
    end
    low < high || return nothing
    return (low, high)
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

function _find_mu_bracket(rho_vals::Vector{Float64}, mu_vals::Vector{Float64},
        mu_lo::Float64, mu_hi::Float64, steps::Int, tol_area::Real)
    attempt = max(steps, 3)
    while attempt <= MAX_CANDIDATE_STEPS
        prev_mu = nothing
        prev_area = nothing
        for mu0 in range(mu_lo, mu_hi; length=attempt)
            area = _area_difference(mu0, rho_vals, mu_vals)
            area === nothing && continue
            if abs(area) <= tol_area
                return (mu0, mu0, area, area)
            end
            if prev_area !== nothing && area * prev_area < 0
                return (prev_mu, mu0, prev_area, area)
            end
            prev_mu = mu0
            prev_area = area
        end
        attempt *= 2
    end
    return nothing
end

function _bisection_solve(rho_vals::Vector{Float64}, mu_vals::Vector{Float64},
        mu_a::Float64, mu_b::Float64, area_a::Real, area_b::Real,
        tol_area::Real, max_iter::Int)

    if mu_a == mu_b
        return mu_a, area_a, 0
    end

    a, b = mu_a, mu_b
    fa, fb = area_a, area_b
    fa * fb <= 0 || return nothing, nothing, 0

    for iter in 1:max_iter
        mid = 0.5 * (a + b)
        area_mid = _area_difference(mid, rho_vals, mu_vals)
        area_mid === nothing && return nothing, nothing, iter
        if abs(area_mid) <= tol_area
            return mid, area_mid, iter
        end
        if fa * area_mid < 0
            b = mid
            fb = area_mid
        else
            a = mid
            fa = area_mid
        end
    end
    mid = 0.5 * (a + b)
    area_mid = _area_difference(mid, rho_vals, mu_vals)
    area_mid === nothing && return nothing, nothing, max_iter
    return mid, area_mid, max_iter
end

function _area_difference(mu0::Float64, rho_vals::Vector{Float64}, mu_vals::Vector{Float64})
    rho_left, rho_right = _find_outer_intersections(mu0, rho_vals, mu_vals)
    if rho_left === nothing || rho_right === nothing || rho_right - rho_left <= 1e-9
        return nothing
    end
    return _integrate_difference(rho_vals, mu_vals, rho_left, rho_right, mu0)
end

function _find_outer_intersections(mu0::Float64, rho_vals::Vector{Float64}, mu_vals::Vector{Float64}; atol::Real=1e-9)
    n = length(rho_vals)
    n >= 2 || return nothing, nothing
    left = nothing
    right = nothing
    for i in 1:(n - 1)
        r1 = rho_vals[i]
        r2 = rho_vals[i + 1]
        f1 = mu_vals[i] - mu0
        f2 = mu_vals[i + 1] - mu0
        if abs(f1) < atol
            left === nothing && (left = r1)
            right = r1
        end
        if f1 == f2
            continue
        end
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

function _integrate_difference(rho_vals::Vector{Float64}, mu_vals::Vector{Float64},
        rho_left::Float64, rho_right::Float64, mu0::Float64)
    total = 0.0
    for i in 1:(length(rho_vals) - 1)
        r1 = rho_vals[i]
        r2 = rho_vals[i + 1]
        if r2 <= rho_left
            continue
        end
        if r1 >= rho_right
            continue
        end
        left = max(r1, rho_left)
        right = min(r2, rho_right)
        if right <= left
            continue
        end
        m1 = mu_vals[i]
        m2 = mu_vals[i + 1]
        mu_left = _interp(r1, r2, m1, m2, left)
        mu_right = _interp(r1, r2, m1, m2, right)
        f_left = mu_left - mu0
        f_right = mu_right - mu0
        total += 0.5 * (f_left + f_right) * (right - left)
    end
    return total
end

@inline function _interp(r1::Float64, r2::Float64, m1::Float64, m2::Float64, target::Float64)
    r2 == r1 && return m1
    t = (target - r1) / (r2 - r1)
    return m1 + t * (m2 - m1)
end

function _failure(reason::AbstractString; kwargs...)
    details = Dict{Symbol, Any}(:reason => reason)
    for (k, v) in kwargs
        details[k] = v
    end
    return MaxwellResult(false, nothing, nothing, nothing, nothing, 0, details)
end

end # module
