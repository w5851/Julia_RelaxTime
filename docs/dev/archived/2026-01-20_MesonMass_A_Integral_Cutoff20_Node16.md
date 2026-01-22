---
title: MesonMass A Integral Cutoff=20 & Nodes=16
archived: true
original: docs/dev/active/2026_1_20介子质量计算结果不同原因分析.md
archived_date: 2026-01-20
---

概要：

- 背景：`ak_compare.csv` 中 Julia 的 A 与 Fortran 差异显著（典型点 `dA_u≈-0.164`），但 `tests/analysis/convergence/test_A_convergence.jl` 的既有参数点仍可通过，产生“测试通过但实际差异大”的现象。
- 结论（基于新增的 ak_compare-like 收敛/报告输出）：在高温/轻质量/Φ≈0.84 的敏感参数区，`pmax=10` 截断不足；使用 `pmax=20` 时，`n=16` 已能把节点误差压到远低于 `TOL=5e-3` 的量级。
- 决策：将介子质量链路与 Fortran A/K 对比脚本中的 A 积分默认配置调整为 `p ∈ [0, 20] fm⁻¹` 且 `n=16`，以在保证精度的同时兼顾性能。

实现变更：

- [src/relaxtime/MesonMass.jl](src/relaxtime/MesonMass.jl)：
  - `ensure_quark_params_has_A` 默认参数改为 `p_nodes=16`、新增 `p_max=20.0`，并用 `gauleg(0.0, p_max, p_nodes)` 生成动量节点。
- `compare_ak_fortran.jl`（已从 scripts/relaxtime 清理；代码见归档）：
  - [docs/dev/archived/2026-01-22_RelaxTime_Fortran_Comparison_Scripts.md](docs/dev/archived/2026-01-22_RelaxTime_Fortran_Comparison_Scripts.md)
  - Julia 侧 A/K 对比改用 `gauleg(0.0, 20.0, 16)`。
- [src/relaxtime/RelaxationTime.jl](src/relaxtime/RelaxationTime.jl)：
  - 计算散射链路中 on-demand 的 `quark_params.A` 改用 `pmax=20,n=16`（更稳健覆盖高温参数区）。
- [src/relaxtime/OneLoopIntegrals.jl](src/relaxtime/OneLoopIntegrals.jl)：
  - 更新 `A` 的文档说明：明确高温/轻质量参数区需要更大的热积分上限，给出 `gauleg(0.0, 20.0, 16)` 的实用折中推荐。

验证方式：

- 运行收敛性测试：`julia --project=. tests/analysis/convergence/test_A_convergence.jl`
- 可选打印报告（显示 pmax=15 vs 20 的差异以及达到 TOL 的最小节点数）：
  - `pwsh`: `$env:A_CONV_REPORT=1; julia --project=. tests/analysis/convergence/test_A_convergence.jl`

备注与后续：

- 该归档聚焦于“Plan Step 3：对齐 A 的积分上限/节点数”的落地；B0/PV 的差异仍需后续按计划继续推进。
- 若未来希望进一步降低计算量，推荐在不改变积分方法的前提下引入“自适应热积分上限”策略（仅在尾部贡献显著的参数点增大 pmax）。

以下为原始内容（保留，以便审阅与历史参考）：

***

目前计算结果见data/outputs/results根目录下，可以看到data\outputs\results\ak_compare.csv文件中A/K计算结果不同，需要分析原因
以及全替换脚本（已从 `scripts/relaxtime` 清理；代码见 docs/dev/archived/2026-01-22_RelaxTime_Fortran_Comparison_Scripts.md）的计算结果也不同，需要分析原因.
## Plan: 介子质量差异排查与修复

本计划目标：定位并修复导致 Julia 与 Fortran 介子质量（A/K 比较及全替换脚本）结果不一致的根本原因。方法：先对齐关键数值实现（B0/PV、对称化、积分节点/截断、常量与映射逻辑），再逐步验证缩小差异，最终回归测试并更新脚本与文档。

### Steps
1. 对齐 B0/PV 实现：检查并统一 src/relaxtime/OneLoopIntegrals.jl 与 RelaxTime/EPNJL_shear-iso-sn-rh-M2/Meson-mass-0703/Cauchy principal value.f90 的主值处理与 gap/节点 策略。  
2. 统一对称化条件：比对并对齐 src/relaxtime/PolarizationAniso.jl 与 RelaxTime/.../polarization.f90 的对称化分支。  
3. 对齐 A 积分节点与上限：在全替换实验脚本（代码见 docs/dev/archived/2026-01-22_RelaxTime_Fortran_Comparison_Scripts.md）与 Julia 主实现中使用相同 Gauss-Legendre 节点数与动量上限。  
4. 检查并修正数据映射：改进 experiment_full_fortran_params.jl 的 key 生成/匹配逻辑，避免舍入导致错行映射（统计 key 命中率并改为精确索引或更严格容差）。  
5. 标准化常量与耦合：比对并统一 src/Constants_PNJL.jl 与 Fortran 常量（hc/hbarc、Λ 等）。  
6. 添加中间量输出与单点比对：在脚本中生成单点 B0、Π、A、K 的完整中间值 CSV，逐列比对并定位残差来源。  
7. 回归验证与文档：更新 docs/dev/active/2026_1_20介子质量计算结果不同原因分析.md，记录变更与验证结果。

### Further Considerations
1. 验证方式：先对齐 PV（最可能原因）→ 若无效则对称化→ 积分节点→ 数据映射。优先按步骤执行以缩小调试范围。  
2. 需要样本：如果可以提供 Fortran 端 A 与 B0 的单点输出样本，可在 1–2 小时内完成首轮比对。  
3. 选项：选择统一 PV 实现（A）或在脚本中为 Fortran/Julia 保持各自实现但增加转换层（B）。推荐先试 A（更可维护）。
