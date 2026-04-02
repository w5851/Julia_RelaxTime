"""
    ProblemSpec

约束问题契约：统一描述 mode、维度以及求解/后处理回调。
"""
struct ProblemSpec{M,R,F,C,U,P,H,S}
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
)
    x_dim > 0 || throw(ArgumentError("x_dim must be positive, got $(x_dim)"))
    theta_dim > 0 || throw(ArgumentError("theta_dim must be positive, got $(theta_dim)"))
    return ProblemSpec(mode, x_dim, theta_dim, residual!, forward_solve, conditions, unpack_solution, postprocess, hard_rules, selector)
end

@inline _unimplemented_residual(F, x, theta, cfg, mode) = throw(ArgumentError("residual! is not configured for $(typeof(mode))"))
@inline _unimplemented_forward_solve(theta, cfg, mode) = throw(ArgumentError("forward_solve is not configured for $(typeof(mode))"))
@inline _unimplemented_conditions(theta, x, meta, cfg, mode) = throw(ArgumentError("conditions is not configured for $(typeof(mode))"))

@inline function build_problem_spec(mode::ConstraintMode; kwargs...)
    if isempty(kwargs)
        components = build_constraint_components(mode)
        return ProblemSpec(
            mode;
            x_dim=constraint_total_dim(components),
            conditions=params -> build_conditions(mode, params),
            forward_solve=(model, T_fm; fwd_kwargs...) -> solve_constraint(model, mode, T_fm; fwd_kwargs...),
        )
    end
    return ProblemSpec(mode; kwargs...)
end
