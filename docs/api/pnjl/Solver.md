# 求解器入口 `PNJL.*`

代码位置：`src/pnjl/PNJL.jl` + `src/pnjl/solver/{ConstraintModes,SeedStrategies,Conditions,ImplicitSolver}.jl`

当前 `PNJL` 顶层按依赖顺序加载 solver 子模块，并直接重导出常用类型与函数，外部统一通过 `PNJL.*` 调用。

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
