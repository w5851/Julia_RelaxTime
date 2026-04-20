struct DiffTarget{F}
    name::Symbol
    evaluator::F
end

# Diff target layer only maps symbolic targets to existing services/results.
# It must not host standalone derivative kernels or duplicate differentiation logic.

@inline DiffTarget(name::Symbol) = DiffTarget(name, _ctx -> throw(ErrorException("DiffTarget $(name) evaluator is not implemented")))

@inline function _theta_value(theta::NamedTuple, key::Symbol)
    hasproperty(theta, key) || throw(ArgumentError("theta is missing required key $(key)"))
    value = getproperty(theta, key)
    value isa Real || throw(ArgumentError("theta.$(key) must be Real, got $(typeof(value))"))
    return Float64(value)
end

@inline function _finite_target_value(name::Symbol, value)
    value isa Real || throw(ArgumentError("target $(name) must evaluate to Real, got $(typeof(value))"))
    out = Float64(value)
    isfinite(out) || throw(ArgumentError("target $(name) evaluated to non-finite value: $(value)"))
    return out
end

@inline function _resolve_target_xi(ctx::ThermoDiffContext)
    if hasproperty(ctx.theta, :xi)
        return _theta_value(ctx.theta, :xi)
    end
    if ctx.spec_override !== nothing && haskey(ctx.spec_override, :xi)
        xi = Float64(get(ctx.spec_override, :xi, 0.0))
        isfinite(xi) || throw(ArgumentError("spec_override.xi must be finite, got $(xi)"))
        return xi
    end
    xi = Float64(ctx.result.xi)
    isfinite(xi) || throw(ArgumentError("context result xi must be finite, got $(xi)"))
    return xi
end

@inline _eval_pressure(ctx::ThermoDiffContext) = _finite_target_value(:pressure, ctx.result.pressure)
@inline _eval_entropy(ctx::ThermoDiffContext) = _finite_target_value(:entropy, ctx.result.entropy)
@inline _eval_rho_norm(ctx::ThermoDiffContext) = _finite_target_value(:rho_norm, ctx.result.rho_norm)
@inline _eval_energy(ctx::ThermoDiffContext) = _finite_target_value(:energy, ctx.result.energy)

@inline function _eval_dP_dT(ctx::ThermoDiffContext)
    deriv = ThermoDerivatives
    T_fm = _theta_value(ctx.theta, :T_fm)
    mu_fm = _theta_value(ctx.theta, :mu_fm)
    xi = _resolve_target_xi(ctx)
    p_num = ctx.spec_override === nothing ? default_momentum_count() : Int(get(ctx.spec_override, :p_num, default_momentum_count()))
    t_num = ctx.spec_override === nothing ? default_theta_count() : Int(get(ctx.spec_override, :t_num, default_theta_count()))
    value = deriv.dP_dT(T_fm, mu_fm; xi=xi, p_num=p_num, t_num=t_num, model=ctx.model)
    return _finite_target_value(:dP_dT, value)
end

@inline function _eval_dP_dmu(ctx::ThermoDiffContext)
    deriv = ThermoDerivatives
    T_fm = _theta_value(ctx.theta, :T_fm)
    mu_fm = _theta_value(ctx.theta, :mu_fm)
    xi = _resolve_target_xi(ctx)
    p_num = ctx.spec_override === nothing ? default_momentum_count() : Int(get(ctx.spec_override, :p_num, default_momentum_count()))
    t_num = ctx.spec_override === nothing ? default_theta_count() : Int(get(ctx.spec_override, :t_num, default_theta_count()))
    value = deriv.dP_dmu(T_fm, mu_fm; xi=xi, p_num=p_num, t_num=t_num, model=ctx.model)
    return _finite_target_value(:dP_dmu, value)
end

const _DIFF_TARGET_REGISTRY = Dict{Symbol, DiffTarget}(
    :pressure => DiffTarget(:pressure, _eval_pressure),
    :entropy => DiffTarget(:entropy, _eval_entropy),
    :rho_norm => DiffTarget(:rho_norm, _eval_rho_norm),
    :energy => DiffTarget(:energy, _eval_energy),
    :dP_dT => DiffTarget(:dP_dT, _eval_dP_dT),
    :dP_dmu => DiffTarget(:dP_dmu, _eval_dP_dmu),
)

@inline function diff_target(name::Symbol)
    haskey(_DIFF_TARGET_REGISTRY, name) || throw(ArgumentError("unknown diff target: $(name)"))
    return _DIFF_TARGET_REGISTRY[name]
end
