@inline function _normalize_param_alias(name::Symbol)
    name === Symbol("μ_fm") && return :mu_fm
    return name
end

@inline function _resolve_target_names(target_names::AbstractVector{Symbol}, targets::AbstractVector{Symbol})
    if !isempty(target_names) && !isempty(targets)
        target_names == targets || throw(ArgumentError("target_names and targets must match when both are provided"))
        return collect(target_names)
    end
    if !isempty(target_names)
        return collect(target_names)
    end
    if !isempty(targets)
        return collect(targets)
    end
    throw(ArgumentError("target_names (or targets) must be non-empty"))
end

@inline function _normalize_param_aliases(names::AbstractVector{Symbol})
    return [_normalize_param_alias(name) for name in names]
end

@inline function _resolve_param_names(param_names::AbstractVector{Symbol}, params::AbstractVector{Symbol})
    normalized_param_names = _normalize_param_aliases(param_names)
    normalized_params = _normalize_param_aliases(params)
    if !isempty(param_names) && !isempty(params)
        normalized_param_names == normalized_params || throw(ArgumentError("param_names and params must match when both are provided"))
        return normalized_param_names
    end
    if !isempty(param_names)
        return normalized_param_names
    end
    if !isempty(params)
        return normalized_params
    end
    throw(ArgumentError("param_names (or params) must be non-empty"))
end

@inline function build_pilot_diff_context(
    result;
    mode,
    model,
    theta::NamedTuple,
    spec_override=nothing,
    jacobian_backend=nothing,
)
    return build_diff_service_context(
        result;
        mode=mode,
        model=model,
        theta=theta,
        spec_override=spec_override,
        jacobian_backend=jacobian_backend,
    )
end

@inline function eval_pilot_derivatives(
    ctx::ThermoDiffContext;
    target_names::AbstractVector{Symbol}=Symbol[],
    param_names::AbstractVector{Symbol}=Symbol[],
    targets::AbstractVector{Symbol}=Symbol[],
    params::AbstractVector{Symbol}=Symbol[],
)
    resolved_targets = _resolve_target_names(target_names, targets)
    resolved_params = _resolve_param_names(param_names, params)

    payload = eval_diff_service_jacobian(
        ctx;
        target_names=resolved_targets,
        param_names=resolved_params,
    )

    return (jacobian=payload.jacobian, by_name=payload.by_name)
end
