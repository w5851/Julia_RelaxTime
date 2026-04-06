using ImplicitDifferentiation

"""
    ImplicitProblem

最小隐式问题契约：
- `forward_solve(θ) -> (x, meta)`
- `conditions(θ, x, z) -> residual`

其中 `z` 为隐式微分库预留上下文参数。
"""
struct ImplicitProblem{F,C}
    forward_solve::F
    conditions::C
    x_dim::Int
    theta_dim::Int

    function ImplicitProblem(forward_solve::F, conditions::C, x_dim::Int, theta_dim::Int) where {F,C}
        x_dim > 0 || throw(ArgumentError("x_dim must be positive, got $(x_dim)"))
        theta_dim > 0 || throw(ArgumentError("theta_dim must be positive, got $(theta_dim)"))
        return new{F,C}(forward_solve, conditions, x_dim, theta_dim)
    end
end

@inline function ImplicitProblem(; forward_solve, conditions, x_dim::Int, theta_dim::Int)
    return ImplicitProblem(forward_solve, conditions, x_dim, theta_dim)
end

Base.@kwdef struct ImplicitSolverConfig
    linear_solver = DirectLinearSolver()
    representation = MatrixRepresentation()
end
