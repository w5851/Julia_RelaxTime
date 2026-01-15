"""
    Solver

PNJL 求解器模块，包含约束模式定义、条件函数和求解器。

## 子模块
- `ConstraintModes`: 求解模式类型定义
- `SeedStrategies`: 初值策略
- `Conditions`: 条件函数构建
- `ImplicitSolver`: 隐函数求解器
"""
module Solver

# 按依赖顺序加载子模块
include("ConstraintModes.jl")
include("SeedStrategies.jl")
include("Conditions.jl")
include("ImplicitSolver.jl")

# 使用子模块
using .ConstraintModes
using .SeedStrategies
using .Conditions
using .ImplicitSolver

# 重新导出 ConstraintModes
export ConstraintMode, FixedMu, FixedRho, FixedEntropy, FixedSigma
export state_dim, param_dim, constraint_description

# 重新导出 SeedStrategies
export SeedStrategy, DefaultSeed, MultiSeed, ContinuitySeed, PhaseAwareSeed, PhaseAwareContinuitySeed
export get_seed, update!, reset!, get_all_seeds, set_phase!
export HADRON_SEED_5, QUARK_SEED_5, HADRON_SEED_8, QUARK_SEED_8
export PhaseBoundaryData, load_phase_boundary, interpolate_mu_c, get_phase_hint

# 重新导出 Conditions
export gap_conditions, build_conditions, build_residual!
export GapParams

# 重新导出 ImplicitSolver
export solve, solve_multi, SolverResult
export create_implicit_solver, solve_with_derivatives

end # module Solver
