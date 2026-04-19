const _THERMO_DIFF_ALLOWED_PARAMS = (:T_fm, :mu_fm, :xi)

struct ParamSpec
    names::Vector{Symbol}
    function ParamSpec(names::AbstractVector{Symbol})
        return new(_normalize_param_names(collect(names)))
    end
end

@inline function _normalize_param_names(names::AbstractVector{Symbol})
    isempty(names) && throw(ArgumentError("ParamSpec requires at least one parameter"))
    allowed = Set(_THERMO_DIFF_ALLOWED_PARAMS)
    for name in names
        name in allowed || throw(ArgumentError("unsupported parameter in ParamSpec: $(name)"))
    end
    dedup = unique(names)
    length(dedup) == length(names) || throw(ArgumentError("ParamSpec contains duplicated parameter names"))
    order = Dict(name => i for (i, name) in enumerate(_THERMO_DIFF_ALLOWED_PARAMS))
    sort!(dedup; by=name -> order[name])
    return dedup
end

@inline ParamSpec(names::Tuple{Vararg{Symbol}}) = ParamSpec(Symbol[names...])

struct ThermoDiffContext{R<:SolverResult,M<:ConstraintMode,MD,TH,SP,JB}
    result::R
    mode::M
    model::MD
    theta::TH
    spec_override::SP
    jacobian_backend::JB
end

@inline function _require_finite_theta_value(theta::NamedTuple, key::Symbol)
    hasproperty(theta, key) || return
    value = getproperty(theta, key)
    value isa Real || throw(ArgumentError("theta.$(key) must be Real, got $(typeof(value))"))
    isfinite(Float64(value)) || throw(ArgumentError("theta.$(key) must be finite, got $(value)"))
end

@inline function _normalize_theta_mu_key(theta::NamedTuple)
    has_mu_ascii = hasproperty(theta, :mu_fm)
    has_mu_unicode = hasproperty(theta, :μ_fm)

    if has_mu_ascii && has_mu_unicode
        throw(ArgumentError("theta cannot contain both :mu_fm and :μ_fm"))
    end

    if has_mu_unicode
        return merge(theta, (mu_fm=getproperty(theta, :μ_fm),))
    end

    return theta
end

@inline function _normalize_diff_theta(theta::NamedTuple)
    theta = _normalize_theta_mu_key(theta)
    _require_finite_theta_value(theta, :T_fm)
    _require_finite_theta_value(theta, :mu_fm)
    _require_finite_theta_value(theta, :xi)
    return theta
end

@inline function _default_jacobian_backend(ctx::ThermoDiffContext, target, params::ParamSpec)
    _ = ctx
    _ = params
    if target isa DiffTarget
        throw(ErrorException("jacobian backend is not implemented for target $(target.name)"))
    elseif target isa AbstractVector{<:DiffTarget}
        names = join(string.(getfield.(target, :name)), ",")
        throw(ErrorException("jacobian backend is not implemented for targets [$(names)]"))
    end
    throw(ErrorException("jacobian backend is not implemented for target input type $(typeof(target))"))
end

@inline function build_thermo_diff_context(
    result::SolverResult;
    mode,
    model,
    theta,
    spec_override=nothing,
    jacobian_backend=_default_jacobian_backend,
)
    mode isa ConstraintMode || throw(ArgumentError("mode must be ConstraintMode, got $(typeof(mode))"))
    theta isa NamedTuple || throw(ArgumentError("theta must be NamedTuple, got $(typeof(theta))"))

    theta_norm = _normalize_diff_theta(theta)
    return ThermoDiffContext(result, mode, model, theta_norm, spec_override, jacobian_backend)
end
