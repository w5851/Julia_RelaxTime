module AdaptiveRhoRefinement

export AdaptiveRhoConfig, suggest_refinement_points, merge_rho_values

"""
    AdaptiveRhoConfig(; slope_tol=5.0, min_gap=0.002, max_points=64, digits=6)

Configuration container that controls how new ρ samples are proposed.
- `slope_tol`: treat |Δμ/Δρ| below this threshold (MeV per normalized ρ) as "flat".
- `min_gap`: only consider segments wider than this spacing.
- `max_points`: cap the number of proposed samples per curve.
- `digits`: rounding precision when deduplicating ρ values.
"""
struct AdaptiveRhoConfig
    slope_tol::Float64
    min_gap::Float64
    max_points::Int
    digits::Int
end

function AdaptiveRhoConfig(; slope_tol::Real=5.0, min_gap::Real=0.002, max_points::Integer=64, digits::Integer=6)
    slope_tol > 0 || error("slope_tol must be positive")
    min_gap > 0 || error("min_gap must be positive")
    max_points > 0 || error("max_points must be positive")
    digits >= 1 || error("digits must be >= 1")
    return AdaptiveRhoConfig(Float64(slope_tol), Float64(min_gap), Int(max_points), Int(digits))
end

"""
    suggest_refinement_points(rho_vals, mu_vals; config=AdaptiveRhoConfig()) -> Vector{Float64}

Given existing samples (ρ, μ), return additional ρ values that should be resampled
because their neighboring slope |Δμ/Δρ| is nearly flat.
"""
function suggest_refinement_points(rho_vals::AbstractVector, mu_vals::AbstractVector; config::AdaptiveRhoConfig=AdaptiveRhoConfig())
    n = min(length(rho_vals), length(mu_vals))
    n >= 2 || return Float64[]
    pairs = Vector{Tuple{Float64, Float64}}(undef, n)
    for i in 1:n
        pairs[i] = (Float64(rho_vals[i]), Float64(mu_vals[i]))
    end
    sort!(pairs; by=first)
    suggestions = Float64[]
    prev_mid = nothing
    for i in 1:(n - 1)
        rho_left, mu_left = pairs[i]
        rho_right, mu_right = pairs[i + 1]
        gap = rho_right - rho_left
        gap > config.min_gap || continue
        slope = abs((mu_right - mu_left) / gap)
        slope <= config.slope_tol || continue
        midpoint = round((rho_left + rho_right) / 2; digits=config.digits)
        if prev_mid === nothing || abs(midpoint - prev_mid) > 10.0^(-config.digits)
            push!(suggestions, midpoint)
            prev_mid = midpoint
            length(suggestions) >= config.max_points && break
        end
    end
    return suggestions
end

"""
    merge_rho_values(existing, additions; digits=6) -> Vector{Float64}

Combine two ρ grids while keeping them sorted and deduplicated at the desired precision.
"""
function merge_rho_values(existing::AbstractVector, additions::AbstractVector; digits::Integer=6)
    digits >= 1 || error("digits must be >= 1")
    rounded = Set{Float64}()
    for value in existing
        push!(rounded, round(Float64(value); digits=digits))
    end
    for value in additions
        push!(rounded, round(Float64(value); digits=digits))
    end
    merged = collect(rounded)
    sort!(merged)
    return merged
end

end # module
