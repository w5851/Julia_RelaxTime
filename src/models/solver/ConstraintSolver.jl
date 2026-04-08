"""constraint_solver.jl

Models 侧约束求解聚合入口。

调用拓扑：
- `Solver.jl` 通过 `solve_constraint` 进入 ProblemSpec 主链分发。
- 各 mode 具体实现在 `ConstraintSolverFixed*.jl`。
- 共性工具在 `ConstraintSolverCommon.jl` 单点维护。
"""

include(joinpath(@__DIR__, "ConstraintSolverCommon.jl"))
include(joinpath(@__DIR__, "ConstraintSolverFixedMu.jl"))
include(joinpath(@__DIR__, "ConstraintSolverFixedRho.jl"))
include(joinpath(@__DIR__, "ConstraintSolverFixedEntropy.jl"))
include(joinpath(@__DIR__, "ConstraintSolverFixedSigma.jl"))
include(joinpath(@__DIR__, "ConstraintSolverFixedAsymmetricRho.jl"))
