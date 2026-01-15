# 求解器聚合模块 `PNJL.Solver`

代码位置：`src/pnjl/solver/Solver.jl`

`PNJL.Solver` 是一个“聚合/重导出”模块：按依赖顺序 include 子模块，并把常用类型与函数统一导出，方便外部通过 `PNJL.*` 调用。

## 子模块

- `ConstraintModes`：约束模式（`FixedMu/FixedRho/FixedEntropy/FixedSigma`）
- `SeedStrategies`：初值策略（`DefaultSeed/MultiSeed/...`）
- `Conditions`：残差与方程构建
- `ImplicitSolver`：主求解器（`solve/solve_multi`）

## 推荐入口

多数场景直接使用顶层 `PNJL.solve`：

```julia
using PNJL
res = PNJL.solve(PNJL.FixedMu(), T_fm, μ_fm; xi=xi)
```

如需更强鲁棒性：

```julia
res = PNJL.solve(PNJL.FixedMu(), T_fm, μ_fm; xi=xi, seed_strategy=PNJL.MultiSeed())
```
