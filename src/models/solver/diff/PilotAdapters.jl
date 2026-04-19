@inline function _normalize_pilot_theta(theta::NamedTuple)
    has_mu_ascii = hasproperty(theta, :mu_fm)
    has_mu_unicode = hasproperty(theta, :μ_fm)
    if has_mu_ascii && has_mu_unicode
        throw(ArgumentError("theta cannot contain both :mu_fm and :μ_fm"))
    end
    if has_mu_unicode
        names = propertynames(theta)
        normalized_names = Tuple((name === :μ_fm ? :mu_fm : name) for name in names)
        values = Tuple(getproperty(theta, name) for name in names)
        return NamedTuple{normalized_names}(values)
    end
    return theta
end

@inline function _normalize_param_alias(name::Symbol)
    name === Symbol("μ_fm") && return :mu_fm
    return name
end

@inline function _require_nonempty_names(names::AbstractVector{Symbol}, kind::String)
    isempty(names) && throw(ArgumentError("$(kind) must be non-empty"))
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

@inline function _build_jacobian_by_name(jac::AbstractMatrix{<:Real}, targets::Vector{DiffTarget}, params::Vector{Symbol})
    by_name = Dict{Symbol, Float64}()
    for i in eachindex(targets)
        target_name = targets[i].name
        for j in eachindex(params)
            key = Symbol("$(target_name)__d$(params[j])")
            by_name[key] = Float64(jac[i, j])
        end
    end
    return by_name
end

@inline function build_pilot_diff_context(
    result;
    mode,
    model,
    theta::NamedTuple,
    spec_override=nothing,
    jacobian_backend=nothing,
)
    theta_norm = _normalize_pilot_theta(theta)
    if jacobian_backend === nothing
        return build_thermo_diff_context(
            result;
            mode=mode,
            model=model,
            theta=theta_norm,
            spec_override=spec_override,
        )
    end
    return build_thermo_diff_context(
        result;
        mode=mode,
        model=model,
        theta=theta_norm,
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

    _require_nonempty_names(resolved_targets, "targets")
    _require_nonempty_names(resolved_params, "params")
    length(unique(resolved_targets)) == length(resolved_targets) || throw(ArgumentError("targets must not contain duplicates"))

    target_defs = DiffTarget[diff_target(name) for name in resolved_targets]
    param_defs = ParamSpec(resolved_params)
    jac = jacobian(ctx, target_defs, param_defs)
    by_name = _build_jacobian_by_name(jac, target_defs, param_defs.names)
    return (jacobian=jac, by_name=by_name)
end
