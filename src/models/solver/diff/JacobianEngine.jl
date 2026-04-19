@inline function _coerce_jacobian_scalar(raw)
    raw isa Real || throw(ArgumentError("jacobian backend for single target must return Real, got $(typeof(raw))"))
    value = Float64(raw)
    isfinite(value) || throw(ArgumentError("jacobian backend returned non-finite scalar: $(raw)"))
    return value
end

@inline function _coerce_jacobian_matrix(raw, n_targets::Int, n_params::Int)
    raw isa AbstractMatrix || throw(ArgumentError("jacobian backend for multi-parameter call must return AbstractMatrix, got $(typeof(raw))"))
    size(raw) == (n_targets, n_params) || throw(ArgumentError("jacobian backend matrix shape mismatch: expected ($(n_targets), $(n_params)), got $(size(raw))"))
    mat = Float64.(raw)
    all(isfinite, mat) || throw(ArgumentError("jacobian backend returned non-finite matrix entries"))
    return mat
end

@inline function _coerce_jacobian_vector(raw, n_targets::Int)
    raw isa AbstractVector || throw(ArgumentError("jacobian backend for multi-target call must return AbstractVector, got $(typeof(raw))"))
    length(raw) == n_targets || throw(ArgumentError("jacobian backend vector length mismatch: expected $(n_targets), got $(length(raw))"))
    vec = Float64.(raw)
    all(isfinite, vec) || throw(ArgumentError("jacobian backend returned non-finite entries for multi-target call"))
    return vec
end

@inline function jacobian(ctx::ThermoDiffContext, target::DiffTarget, params::ParamSpec)
    n_params = length(params.names)
    raw = ctx.jacobian_backend(ctx, target, params)
    if n_params == 1
        value = _coerce_jacobian_scalar(raw)
        return reshape([value], 1, 1)
    end
    return _coerce_jacobian_matrix(raw, 1, n_params)
end

@inline function jacobian(ctx::ThermoDiffContext, targets::AbstractVector{<:DiffTarget}, params::ParamSpec)
    isempty(targets) && throw(ArgumentError("targets must be non-empty"))
    n_params = length(params.names)
    raw = ctx.jacobian_backend(ctx, targets, params)
    if n_params == 1
        vec = _coerce_jacobian_vector(raw, length(targets))
        return reshape(vec, length(targets), 1)
    end
    return _coerce_jacobian_matrix(raw, length(targets), n_params)
end
