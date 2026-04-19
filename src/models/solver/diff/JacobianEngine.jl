@inline function _fd_step(value::Float64)
    base = max(abs(value), 1.0)
    return max(1e-5, sqrt(eps(Float64)) * base)
end

@inline function _resolve_theta_value(theta::NamedTuple, key::Symbol)
    hasproperty(theta, key) || throw(ArgumentError("theta is missing required key $(key)"))
    value = getproperty(theta, key)
    value isa Real || throw(ArgumentError("theta.$(key) must be Real, got $(typeof(value))"))
    return Float64(value)
end

@inline function _theta_with_perturb(theta::NamedTuple, key::Symbol, delta::Float64)
    value = _resolve_theta_value(theta, key)
    return merge(theta, (key => value + delta,))
end

@inline function _runtime_option(ctx::ThermoDiffContext, key::Symbol, default)
    if ctx.spec_override === nothing
        return default
    end
    return get(ctx.spec_override, key, default)
end

@inline function _solve_result_at_theta(ctx::ThermoDiffContext, theta::NamedTuple)
    if ctx.mode isa FixedMu
        T_fm = _resolve_theta_value(theta, :T_fm)
        mu_fm = _resolve_theta_value(theta, :mu_fm)
        xi = hasproperty(theta, :xi) ? Float64(getproperty(theta, :xi)) : Float64(_runtime_option(ctx, :xi, ctx.result.xi))
        p_num = Int(_runtime_option(ctx, :p_num, Main.Models.default_momentum_count()))
        t_num = Int(_runtime_option(ctx, :t_num, Main.Models.default_theta_count()))
        residual_norm_max = Float64(_runtime_option(ctx, :residual_norm_max, 1e-6))
        return Main.Models.solve(ctx.model, ctx.mode, T_fm, mu_fm;
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            residual_norm_max=residual_norm_max,
        )
    end
    throw(ArgumentError("default jacobian backend currently supports FixedMu mode only, got $(typeof(ctx.mode))"))
end

@inline function _evaluate_target_at_theta(ctx::ThermoDiffContext, target::DiffTarget, theta::NamedTuple)
    solved = _solve_result_at_theta(ctx, theta)
    local_ctx = ThermoDiffContext(solved, ctx.mode, ctx.model, theta, ctx.spec_override, ctx.jacobian_backend)
    value = target.evaluator(local_ctx)
    value isa Real || throw(ArgumentError("target $(target.name) evaluator must return Real, got $(typeof(value))"))
    out = Float64(value)
    isfinite(out) || throw(ArgumentError("target $(target.name) evaluator returned non-finite value"))
    return out
end

@inline function _evaluate_targets_at_theta(ctx::ThermoDiffContext, targets::AbstractVector{<:DiffTarget}, theta::NamedTuple)
    solved = _solve_result_at_theta(ctx, theta)
    local_ctx = ThermoDiffContext(solved, ctx.mode, ctx.model, theta, ctx.spec_override, ctx.jacobian_backend)
    out = Vector{Float64}(undef, length(targets))
    for i in eachindex(targets)
        value = targets[i].evaluator(local_ctx)
        value isa Real || throw(ArgumentError("target $(targets[i].name) evaluator must return Real, got $(typeof(value))"))
        out[i] = Float64(value)
        isfinite(out[i]) || throw(ArgumentError("target $(targets[i].name) evaluator returned non-finite value"))
    end
    return out
end

@inline function _default_numeric_jacobian(ctx::ThermoDiffContext, target::DiffTarget, params::ParamSpec)
    n_params = length(params.names)
    if n_params == 1
        key = params.names[1]
        center = _resolve_theta_value(ctx.theta, key)
        h = _fd_step(center)
        f_p = _evaluate_target_at_theta(ctx, target, _theta_with_perturb(ctx.theta, key, h))
        f_m = _evaluate_target_at_theta(ctx, target, _theta_with_perturb(ctx.theta, key, -h))
        return (f_p - f_m) / (2h)
    end

    out = Matrix{Float64}(undef, 1, n_params)
    for (j, key) in enumerate(params.names)
        center = _resolve_theta_value(ctx.theta, key)
        h = _fd_step(center)
        f_p = _evaluate_target_at_theta(ctx, target, _theta_with_perturb(ctx.theta, key, h))
        f_m = _evaluate_target_at_theta(ctx, target, _theta_with_perturb(ctx.theta, key, -h))
        out[1, j] = (f_p - f_m) / (2h)
    end
    return out
end

@inline function _default_numeric_jacobian(ctx::ThermoDiffContext, targets::AbstractVector{<:DiffTarget}, params::ParamSpec)
    n_targets = length(targets)
    n_targets > 0 || throw(ArgumentError("targets must be non-empty"))
    n_params = length(params.names)

    if n_params == 1
        key = params.names[1]
        center = _resolve_theta_value(ctx.theta, key)
        h = _fd_step(center)
        out = Vector{Float64}(undef, n_targets)
        theta_p = _theta_with_perturb(ctx.theta, key, h)
        theta_m = _theta_with_perturb(ctx.theta, key, -h)
        values_p = _evaluate_targets_at_theta(ctx, targets, theta_p)
        values_m = _evaluate_targets_at_theta(ctx, targets, theta_m)
        @inbounds for i in eachindex(targets)
            out[i] = (values_p[i] - values_m[i]) / (2h)
        end
        return out
    end

    out = Matrix{Float64}(undef, n_targets, n_params)
    for (j, key) in enumerate(params.names)
        center = _resolve_theta_value(ctx.theta, key)
        h = _fd_step(center)
        theta_p = _theta_with_perturb(ctx.theta, key, h)
        theta_m = _theta_with_perturb(ctx.theta, key, -h)
        values_p = _evaluate_targets_at_theta(ctx, targets, theta_p)
        values_m = _evaluate_targets_at_theta(ctx, targets, theta_m)
        @inbounds for i in eachindex(targets)
            out[i, j] = (values_p[i] - values_m[i]) / (2h)
        end
    end
    return out
end

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
    raw = if ctx.jacobian_backend === _default_jacobian_backend
        _default_numeric_jacobian(ctx, target, params)
    else
        ctx.jacobian_backend(ctx, target, params)
    end
    if n_params == 1
        value = _coerce_jacobian_scalar(raw)
        return reshape([value], 1, 1)
    end
    return _coerce_jacobian_matrix(raw, 1, n_params)
end

@inline function jacobian(ctx::ThermoDiffContext, targets::AbstractVector{<:DiffTarget}, params::ParamSpec)
    isempty(targets) && throw(ArgumentError("targets must be non-empty"))
    n_params = length(params.names)
    raw = if ctx.jacobian_backend === _default_jacobian_backend
        _default_numeric_jacobian(ctx, targets, params)
    else
        ctx.jacobian_backend(ctx, targets, params)
    end
    if n_params == 1
        vec = _coerce_jacobian_vector(raw, length(targets))
        return reshape(vec, length(targets), 1)
    end
    return _coerce_jacobian_matrix(raw, length(targets), n_params)
end
