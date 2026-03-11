module ValidationUtils

export require_finite
export require_positive_finite
export require_nonnegative_finite
export require_positive_integer
export validate_grid_weight_pair

@inline function require_finite(name::AbstractString, value)
    isfinite(value) || throw(ArgumentError("$name must be finite"))
    return value
end

@inline function require_positive_finite(name::AbstractString, value)
    isfinite(value) && value > 0 || throw(ArgumentError("$name must be finite and > 0"))
    return value
end

@inline function require_nonnegative_finite(name::AbstractString, value)
    isfinite(value) && value >= 0 || throw(ArgumentError("$name must be finite and >= 0"))
    return value
end

@inline function require_positive_integer(name::AbstractString, value::Integer)
    value > 0 || throw(ArgumentError("$name must be > 0"))
    return value
end

function validate_grid_weight_pair(
    context::AbstractString,
    grid_name::AbstractString,
    grid,
    weight_name::AbstractString,
    weights,
)
    if (grid === nothing) != (weights === nothing)
        throw(ArgumentError("$context: $grid_name and $weight_name must be provided together"))
    end

    grid === nothing && return nothing

    length(grid) == length(weights) || throw(ArgumentError("$context: $grid_name and $weight_name length mismatch"))
    isempty(grid) && throw(ArgumentError("$context: $grid_name must not be empty"))
    all(isfinite, grid) || throw(ArgumentError("$context: $grid_name must be finite"))
    all(isfinite, weights) || throw(ArgumentError("$context: $weight_name must be finite"))
    return nothing
end

end # module ValidationUtils