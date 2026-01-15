# 初值策略模块 `PNJL.SeedStrategies`

代码位置：`src/pnjl/solver/SeedStrategies.jl`

本模块负责为 PNJL 隐式求解器提供“初值（seed）”生成与管理策略，用于提升收敛、避免落入错误分支，并在参数扫描时实现连续性跟踪。

## 关键概念

- **seed**：求解器未知量的初始猜测向量。
  - `FixedMu()` 模式：5 维 `x = [φ_u, φ_d, φ_s, Φ, Φ̄]`
  - `FixedRho/FixedEntropy/FixedSigma()` 模式：8 维（额外包含 `μ_u, μ_d, μ_s`）
- **多值解**：一阶相变区域存在强子/夸克等多个可收敛分支，需要多初值或分支追踪。

## 导出符号（摘要）

- 策略类型：`SeedStrategy`, `DefaultSeed`, `MultiSeed`, `ContinuitySeed`, `PhaseAwareSeed`, `PhaseAwareContinuitySeed`
- 核心接口：`get_seed`, `get_all_seeds`, `update!`, `reset!`, `set_phase!`
- 相变线工具：`PhaseBoundaryData`, `load_phase_boundary`, `interpolate_mu_c`, `get_phase_hint`

## 推荐用法

### 单点求解（更稳健：MultiSeed）

```julia
using PNJL

res = PNJL.solve(PNJL.FixedMu(), T_fm, μ_fm; xi=xi, seed_strategy=PNJL.MultiSeed())
```

### 扫描（推荐：首点 MultiSeed + 后续相变感知连续跟踪）

`PhaseAwareContinuitySeed` 本身是“有状态”的跟踪器；连续性沿你调用 `solve` 的顺序推进（见下方“连续性路径”）。

```julia
using PNJL

tracker = PNJL.PhaseAwareContinuitySeed(xi; bootstrap_multiseed=true)

for (T_fm, μ_fm, T_MeV, μ_MeV) in scan_points
    res = PNJL.solve(PNJL.FixedMu(), T_fm, μ_fm; xi=xi, seed_strategy=tracker)
    res.converged && PNJL.update!(tracker, res.solution, T_MeV, μ_MeV)
end
```

注意：`bootstrap_multiseed=true` 仅在 `FixedMu()` 下由求解器触发 `solve_multi`（多初值求解并选择 Ω 最小者）。

## 连续性路径（二维/多维扫描）

`ContinuitySeed/PhaseAwareContinuitySeed` 只保存“上一个点”的解，因此：

- 在 **二维网格 (T, μ)** 上，连续性沿你遍历网格的 **路径/顺序** 生效。
- 实践中常见两种方式：
  1) 固定 `T`，沿 `μ` 扫描：每个 `T` 维护一个 tracker（进入新 T 时 reset）
  2) 固定 `μ`，沿 `T` 扫描：每个 `μ` 维护一个 tracker

仓库示例：`PNJL.TmuScan` 采用“每个温度一条 μ 扫描线”，即每个 `T` 开始前重置 tracker。

## 相关模块

- `PNJL.ImplicitSolver`：真正消费这些策略并执行求解
- `PNJL.DualBranchScan`：一阶相变区域的双分支扫描（显式追踪两分支并比较 Ω）
