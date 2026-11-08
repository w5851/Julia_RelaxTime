module MatchedPALCPilotSupport

using LinearAlgebra

"""Seeds from every crossing of a target density, including folded branches.

Interpolation only supplies a Newton seed; the caller must solve at `target`.
An endpoint seed is retained when a continuation stopped just short of a target.
"""
function target_seeds(points, target; dedup_tol=1e-8)
    isempty(points) && return Vector{Vector{Float64}}()
    seeds = Vector{Vector{Float64}}()
    function add!(x)
        any(s -> norm(s - x) <= dedup_tol, seeds) || push!(seeds, Float64.(x))
    end
    for p in points
        abs(p.rho - target) <= 1e-12 && add!(p.solution)
    end
    for i in 2:length(points)
        a, b = points[i - 1], points[i]
        a.rho == b.rho && continue
        if min(a.rho, b.rho) < target < max(a.rho, b.rho)
            weight = (target - a.rho) / (b.rho - a.rho)
            add!((1 - weight) .* a.solution .+ weight .* b.solution)
        end
    end
    if isempty(seeds)
        nearest = points[argmin(abs(p.rho - target) for p in points)]
        add!(nearest.solution)
    end
    return seeds
end

"""Observe every BK Newton attempt, including retries and duplicate final callbacks."""
Base.@kwdef mutable struct NewtonWork
    attempts::Int = 0
    iterations::Int = 0
    nonconverged::Int = 0
    history::Any = nothing
    last_step::Int = 0
    last_residual::Float64 = Inf
    last_tol::Float64 = 0.0
end

function finish_attempt!(work::NewtonWork)
    work.history === nothing && return
    work.nonconverged += !(isfinite(work.last_residual) && work.last_residual < work.last_tol)
    work.history = nothing
end

function observe_newton!(work::NewtonWork, state)
    if work.history !== state.residuals
        finish_attempt!(work)
        work.history = state.residuals
        work.attempts += 1
        work.last_step = 0
    end
    work.iterations += max(0, state.step - work.last_step)
    work.last_step = state.step
    work.last_residual = state.residual
    work.last_tol = hasproperty(state, :contparams) ? state.contparams.newton_options.tol : state.options.tol
    return true
end

end
