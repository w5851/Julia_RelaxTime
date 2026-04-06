"""
    build_implicit_solver(problem::ImplicitProblem, cfg::ImplicitSolverConfig=ImplicitSolverConfig())

根据最小隐式问题契约构建 `ImplicitFunction`。
"""
@inline function build_implicit_solver(
    problem::ImplicitProblem,
    cfg::ImplicitSolverConfig=ImplicitSolverConfig(),
)
    return ImplicitFunction(
        problem.forward_solve,
        problem.conditions;
        linear_solver=cfg.linear_solver,
        representation=cfg.representation,
    )
end
