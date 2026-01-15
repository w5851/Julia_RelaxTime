# 双分支扫描模块 `PNJL.DualBranchScan`

代码位置：`src/pnjl/scans/DualBranchScan.jl`

本模块用于一阶相变区域的“多值解”处理：显式追踪强子/夸克两条分支，并通过比较 Ω 选择物理解，同时定位相变点（Ω 交叉）与自旋odal 边界。

## 入口

```julia
res = PNJL.run_dual_branch_scan(; T_mev=100.0, mu_range=0.0:10.0:400.0, xi=0.0)
info = PNJL.find_phase_transition(res)
```

## 适用场景

- 你需要明确得到两条分支（而不是“跳过相变线”）。
- 你需要找相变点 `μ_c`、共存区与 spinodal 位置。

## 与 PhaseAwareContinuitySeed 的区别

- `PhaseAwareContinuitySeed`：更像工程化扫描的“稳健初值器”，目标是让单分支扫描不中断/不乱跳。
- `DualBranchScan`：更偏物理分析用途，目标是得到两条分支并作自由能比较。
