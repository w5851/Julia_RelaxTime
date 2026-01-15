# 隐式求解器模块 `PNJL.ImplicitSolver`

代码位置：`src/pnjl/solver/ImplicitSolver.jl`

本模块提供 PNJL 平衡方程的统一求解接口（基于 `NLsolve.jl`），并封装物理性判据、Newton/Trust-Region 兜底，以及多初值求解。

## 核心 API

### `solve`

固定化学势：

```julia
res = PNJL.solve(PNJL.FixedMu(), T_fm, μ_fm; xi=0.0, seed_strategy=PNJL.DefaultSeed())
```

固定密度/熵/比熵：

```julia
res = PNJL.solve(PNJL.FixedRho(ρnorm), T_fm; xi=0.0, seed_strategy=PNJL.DefaultSeed())
```

### `solve_multi`

对 `FixedMu()` 使用 `MultiSeed` 的多初值求解，并按 Ω 最小选择：

```julia
res = PNJL.solve_multi(PNJL.FixedMu(), T_fm, μ_fm; xi=xi, seed_strategy=PNJL.MultiSeed())
```

## 与 SeedStrategies 的关系

- 若 `seed_strategy isa MultiSeed`：`solve` 会直接走 `solve_multi`。
- 若 `seed_strategy isa PhaseAwareContinuitySeed` 且 `bootstrap_multiseed=true` 且当前为“首点”（`previous_solution === nothing`）：`solve` 会对首点调用 `solve_multi` 完成自举。
- 其他策略：通过 `get_seed(seed_strategy, θ, mode)` 生成单一初值，再执行 `nlsolve`。

## 连续性跟踪的注意点

`PhaseAwareContinuitySeed` 仅提供“初值”，不会自动更新其内部状态；扫描脚本需要在收敛后调用：

```julia
PNJL.update!(tracker, res.solution, T_MeV, μ_MeV)
```

其中 `update!` 的 `T_MeV/μ_MeV` 使用 MeV 单位（与相变线 CSV 数据一致）。
