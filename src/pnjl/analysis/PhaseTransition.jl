module PhaseTransition

export SShapeResult, detect_s_shape, group_curves_by_temperature

const EPS_SLOPE = 0.0

struct SShapeResult
    has_s_shape::Bool
    mu_spinodal_low::Union{Nothing, Float64}
    mu_spinodal_high::Union{Nothing, Float64}
    rho_spinodal_low::Union{Nothing, Float64}
    rho_spinodal_high::Union{Nothing, Float64}
    derivative_sign_changes::Int
end

SShapeResult() = SShapeResult(false, nothing, nothing, nothing, nothing, 0)

"""Normalize curve samples by sorting in ascending μ order."""
function _sort_curve(mu_vals::AbstractVector, rho_vals::AbstractVector)
    n = min(length(mu_vals), length(rho_vals))
    order = sortperm(mu_vals[1:n])
    return Float64.(mu_vals[order]), Float64.(rho_vals[order])
end

@inline _slope_sign(value) = abs(value) <= EPS_SLOPE ? 0 : (value > 0 ? 1 : -1)

"""
    detect_s_shape(mu_vals, rho_vals; eps=EPS_SLOPE, min_points=5)

Check whether a ρ(μ) curve exhibits an S-shape (i.e. the derivative changes sign from
positive → negative → positive). Returns an `SShapeResult` containing spinodal estimates.
"""
function detect_s_shape(mu_vals::AbstractVector, rho_vals::AbstractVector; eps::Real=EPS_SLOPE, min_points::Int=5)
    n = min(length(mu_vals), length(rho_vals))
    n >= min_points || return SShapeResult()
    mu_sorted, rho_sorted = _sort_curve(mu_vals, rho_vals)
    slopes = Float64[]
    for i in 1:(length(mu_sorted) - 1)
        dmu = mu_sorted[i + 1] - mu_sorted[i]
        abs(dmu) <= eps && continue
        push!(slopes, (rho_sorted[i + 1] - rho_sorted[i]) / dmu)
    end
    isempty(slopes) && return SShapeResult()

    signs = Int[]
    last = 0
    for slope in slopes
        sign = abs(slope) < eps ? 0 : (slope > 0 ? 1 : -1)
        sign == 0 && continue
        sign == last && continue
        push!(signs, sign)
        last = sign
    end
    length(signs) < 3 && return SShapeResult(false, nothing, nothing, nothing, nothing, length(signs) - 1)

    first_neg_idx = findfirst(x -> x == -1, signs)
    last_pos_after = nothing
    if first_neg_idx !== nothing
        for idx in (first_neg_idx + 1):length(signs)
            idx > length(signs) && break
            if signs[idx] == 1
                last_pos_after = idx
            end
        end
    end
    if first_neg_idx === nothing || last_pos_after === nothing
        return SShapeResult(false, nothing, nothing, nothing, nothing, length(signs) - 1)
    end

    # map sign indices back to slope indices
    slope_indices = Int[]
    prev_sign = 0
    for (i, slope) in pairs(slopes)
        s = _slope_sign(slope)
        s == 0 && continue
        if s != prev_sign
            push!(slope_indices, i)
            prev_sign = s
        end
    end

    low_idx = slope_indices[first_neg_idx]
    high_idx = slope_indices[last_pos_after]
    mu_low = (mu_sorted[low_idx] + mu_sorted[low_idx + 1]) / 2
    mu_high = (mu_sorted[high_idx] + mu_sorted[high_idx + 1]) / 2
    rho_low = (rho_sorted[low_idx] + rho_sorted[low_idx + 1]) / 2
    rho_high = (rho_sorted[high_idx] + rho_sorted[high_idx + 1]) / 2

    return SShapeResult(true, mu_low, mu_high, rho_low, rho_high, length(signs) - 1)
end

"""Group raw CSV-like rows by temperature for downstream analysis."""
function group_curves_by_temperature(rows; xi::Real=0.0, tol::Real=1e-6)
    grouped = Dict{Float64, Vector{Tuple{Float64, Float64}}}()
    for row in rows
        T = try
            parse(Float64, row["T_MeV"])
        catch
            continue
        end
        xi_val = haskey(row, "xi") ? try
                parse(Float64, row["xi"])
            catch
                NaN
            end : NaN
        if isnan(xi_val) || abs(xi_val - xi) > tol
            continue
        end
        rho = try
            parse(Float64, row["rho"])
        catch
            continue
        end
        mu = haskey(row, "mu_avg_MeV") ? try
                parse(Float64, row["mu_avg_MeV"])
            catch
                NaN
            end : NaN
        if isnan(mu)
            mu = haskey(row, "mu_MeV") ? try
                    parse(Float64, row["mu_MeV"])
                catch
                    NaN
                end : NaN
        end
        if isnan(mu) || !isfinite(rho)
            continue
        end
        bucket = get!(grouped, T) do
            Vector{Tuple{Float64, Float64}}()
        end
        push!(bucket, (mu, rho))
    end
    return grouped
end

end # module
