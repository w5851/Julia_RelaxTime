@inline function build_diff_service_context(
    result;
    mode,
    model,
    theta,
    spec_override=nothing,
    jacobian_backend=nothing,
)
    theta isa NamedTuple || throw(ArgumentError("theta must be NamedTuple, got $(typeof(theta))"))

    if jacobian_backend === nothing
        return build_thermo_diff_context(
            result;
            mode=mode,
            model=model,
            theta=theta,
            spec_override=spec_override,
        )
    end

    return build_thermo_diff_context(
        result;
        mode=mode,
        model=model,
        theta=theta,
        spec_override=spec_override,
        jacobian_backend=jacobian_backend,
    )
end

@inline function _normalize_diff_service_param_name(name::Symbol)
    name === Symbol("μ_fm") && return :mu_fm
    return name
end

@inline function eval_diff_service_jacobian(
    ctx::ThermoDiffContext;
    target_names::AbstractVector{Symbol},
    param_names::AbstractVector{Symbol},
)
    isempty(target_names) && throw(ArgumentError("target_names must be non-empty"))
    length(unique(target_names)) == length(target_names) || throw(ArgumentError("target list must not contain duplicates"))
    isempty(param_names) && throw(ArgumentError("param_names must be non-empty"))

    normalized_param_names = _normalize_diff_service_param_name.(param_names)
    length(unique(normalized_param_names)) == length(normalized_param_names) || throw(ArgumentError("param_names contains duplicated names after normalization"))

    resolved_target_names = collect(target_names)
    target_defs = DiffTarget[diff_target(name) for name in resolved_target_names]
    param_spec = ParamSpec(normalized_param_names)
    jac = jacobian(ctx, target_defs, param_spec)

    by_name = Dict{Symbol, Float64}()
    for i in eachindex(resolved_target_names)
        target_name = resolved_target_names[i]
        for j in eachindex(param_spec.names)
            key = Symbol("$(target_name)__d$(param_spec.names[j])")
            by_name[key] = Float64(jac[i, j])
        end
    end

    return (
        jacobian=jac,
        target_names=resolved_target_names,
        param_names=param_spec.names,
        by_name=by_name,
    )
end
