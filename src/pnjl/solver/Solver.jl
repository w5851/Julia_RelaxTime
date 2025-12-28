"""
    Solver

PNJL 求解器模块，包含约束模式定义和条件函数。

## 子模块
- `ConstraintModes`: 求解模式类型定义
- `Conditions`: 条件函数构建
"""
module Solver

include("ConstraintModes.jl")
include("Conditions.jl")

using .ConstraintModes
using .Conditions

# 重新导出
export ConstraintMode, FixedMu, FixedRho, FixedEntropy, FixedSigma
export state_dim, param_dim, constraint_description

export gap_conditions, build_conditions, build_residual!
export GapParams

end # module Solver
