"""
    build_implicit_solver(problem::ImplicitProblem, cfg::ImplicitSolverConfig=ImplicitSolverConfig())

legacy `ImplicitDifferentiation` builder 已下线。
"""
@inline function build_implicit_solver(
    problem::ImplicitProblem,
    cfg::ImplicitSolverConfig=ImplicitSolverConfig(),
)
    _ = problem
    _ = cfg
    throw(ArgumentError("build_implicit_solver has been retired with the ImplicitDifferentiation backend. Use build_*_problem(...).forward_solve/conditions for residual adapter audits, or solve_pnjl_with_derivatives(...; derivative_backend=:auto/:taylordiff) for PNJL derivatives."))
end
