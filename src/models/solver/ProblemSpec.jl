"""
    ProblemSpec

约束问题契约：统一描述 mode、维度以及求解/后处理回调。
"""

# Contract Core: constraints extension and ProblemSpec data contract.
struct ExtraConstraints{R,F,S}
    residual!::R
    feasible::F
    seed_extend::S
end

@inline _default_extra_residual!(F, x, theta, cfg, mode) = nothing
@inline _default_extra_feasible(candidate, params, mode) = true
@inline _default_seed_extend(seed, mode) = Float64.(seed)

@inline function default_extra_constraints()
    return ExtraConstraints(
        _default_extra_residual!,
        _default_extra_feasible,
        _default_seed_extend,
    )
end

@inline function _resolve_extra_constraints(kwargs)
    extra = get(kwargs, :extra_constraints, default_extra_constraints())
    extra isa ExtraConstraints || throw(ArgumentError("extra_constraints must be ExtraConstraints, got $(typeof(extra))"))
    return extra
end

@inline function _extend_seed_with_extra(seed::AbstractVector{<:Real}, mode::ConstraintMode, extra::ExtraConstraints)
    extended = extra.seed_extend(Float64.(seed), mode)
    return Float64.(extended)
end

@inline function _append_extra_feasible_rule(rules, extra::ExtraConstraints, mode::ConstraintMode)
    extra_rule = c -> (Bool(extra.feasible(c, nothing, mode)), :extra_constraint_failed)
    return vcat(collect(rules), [extra_rule])
end



struct ProblemSpec{M,R,F,C,U,P,H,S,E}
    mode::M
    x_dim::Int
    theta_dim::Int
    residual!::R
    forward_solve::F
    conditions::C
    unpack_solution::U
    postprocess::P
    hard_rules::H
    selector::S
    extra_constraints::E
end

@inline function ProblemSpec(
    mode::ConstraintMode;
    x_dim::Int=state_dim(mode),
    theta_dim::Int=param_dim(mode),
    residual! = _unimplemented_residual,
    forward_solve=_unimplemented_forward_solve,
    conditions=_unimplemented_conditions,
    unpack_solution=identity,
    postprocess=identity,
    hard_rules=default_hard_constraint_rules(),
    selector=select_pressure_max_candidate,
    extra_constraints=default_extra_constraints(),
)
    x_dim > 0 || throw(ArgumentError("x_dim must be positive, got $(x_dim)"))
    theta_dim > 0 || throw(ArgumentError("theta_dim must be positive, got $(theta_dim)"))
    wrapped_forward_solve = function (model, T_fm; fwd_kwargs...)
        if haskey(fwd_kwargs, :extra_constraints)
            return forward_solve(model, T_fm; fwd_kwargs...)
        end
        return forward_solve(model, T_fm; extra_constraints=extra_constraints, fwd_kwargs...)
    end

    return ProblemSpec(mode, x_dim, theta_dim, residual!, wrapped_forward_solve, conditions, unpack_solution, postprocess, hard_rules, selector, extra_constraints)
end

@inline _unimplemented_residual(F, x, theta, cfg, mode) = throw(ArgumentError("residual! is not configured for $(typeof(mode))"))
@inline _unimplemented_forward_solve(theta, cfg, mode) = throw(ArgumentError("forward_solve is not configured for $(typeof(mode))"))
@inline _unimplemented_conditions(theta, x, meta, cfg, mode) = throw(ArgumentError("conditions is not configured for $(typeof(mode))"))



# Registration Core: build mode -> ProblemSpec.forward_solve 唯一注册入口。
@inline function build_problem_spec(mode::ConstraintMode; kwargs...)
    if isempty(kwargs)
        components = build_constraint_components(mode)
        extra_constraints = default_extra_constraints()
        forward_solve = if mode isa FixedRho
            (model, T_fm; fwd_kwargs...) -> _fixedrho_problem_spec_forward_solve(model, mode, T_fm; fwd_kwargs...)
        elseif mode isa FixedMu
            (model, T_fm; fwd_kwargs...) -> _fixedmu_problem_spec_forward_solve(model, mode, T_fm; fwd_kwargs...)
        elseif mode isa FixedEntropy
            (model, T_fm; fwd_kwargs...) -> _fixedentropy_problem_spec_forward_solve(model, mode, T_fm; fwd_kwargs...)
        elseif mode isa FixedSigma
            (model, T_fm; fwd_kwargs...) -> _fixedsigma_problem_spec_forward_solve(model, mode, T_fm; fwd_kwargs...)
        elseif mode isa FixedAsymmetricRho
            (model, T_fm; fwd_kwargs...) -> _fixedasymrho_problem_spec_forward_solve(model, mode, T_fm; fwd_kwargs...)
        else
            (model, T_fm; fwd_kwargs...) -> solve_constraint(model, mode, T_fm; fwd_kwargs...)
        end
        return ProblemSpec(
            mode;
            x_dim=constraint_total_dim(components),
            conditions=params -> build_conditions(mode, params),
            forward_solve=forward_solve,
            extra_constraints=extra_constraints,
        )
    end
    return ProblemSpec(mode; kwargs...)
end
