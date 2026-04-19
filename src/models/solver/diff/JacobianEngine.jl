@inline function _coerce_jacobian_scalar(raw)
    raw isa Real || throw(ArgumentError("jacobian backend for single target must return Real, got $(typeof(raw))"))
    value = Float64(raw)
    isfinite(value) || throw(ArgumentError("jacobian backend returned non-finite scalar: $(raw)"))
    return value
end

@inline function _coerce_jacobian_vector(raw, n_targets::Int)
    raw isa AbstractVector || throw(ArgumentError("jacobian backend for multi-target call must return AbstractVector, got $(typeof(raw))"))
    length(raw) == n_targets || throw(ArgumentError("jacobian backend vector length mismatch: expected $(n_targets), got $(length(raw))"))
    vec = Float64.(raw)
    all(isfinite, vec) || throw(ArgumentError("jacobian backend returned non-finite entries for multi-target call"))
    return vec
end

@inline function jacobian(ctx::ThermoDiffContext, target::DiffTarget, params::ParamSpec)
    raw = ctx.jacobian_backend(ctx, target, params)
    value = _coerce_jacobian_scalar(raw)
    return reshape([value], 1, 1)
end

@inline function jacobian(ctx::ThermoDiffContext, targets::AbstractVector{<:DiffTarget}, params::ParamSpec)
    isempty(targets) && throw(ArgumentError("targets must be non-empty"))
    raw = ctx.jacobian_backend(ctx, targets, params)
    vec = _coerce_jacobian_vector(raw, length(targets))
    return reshape(vec, length(targets), 1)
end
