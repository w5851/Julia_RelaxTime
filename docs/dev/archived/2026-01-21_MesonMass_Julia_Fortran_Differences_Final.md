---
title: MesonMass Julia vs Fortran Differences (Final Archive)
archived: true
original: docs/dev/active/2026_1_20介子质量计算结果不同原因分析.md
archived_date: 2026-01-21
---

结论：

- A 的主要差异已定位为“热积分动量上限截断不一致”，并已在 Julia 侧关键链路中统一为 `pmax=20, n=16` 的默认折中；对应实现与验证已单独归档在 [docs/dev/archived/2026-01-20_MesonMass_A_Integral_Cutoff20_Node16.md](docs/dev/archived/2026-01-20_MesonMass_A_Integral_Cutoff20_Node16.md)。
- 在采用上述 A 配置后，π 通道误差已很小；剩余少量点上 K 的误差更显著。
- 基于 `data/outputs/results/meson_full_fortran_params_experiment.csv` 的“交叉组合残差”归因与统计：剩余质量差异与 `B0`（尤其 PV/数值策略差异导致的 `B0` 实部偏差）高度相关；当仅筛选 `|dM_K| > |dM_π|` 的点时，K 的 `|dB0_K|` 明显大于 π 的 `|dB0_π|`。
- 因此本任务在当前阶段以“接受 Fortran 与 Julia 的 `B0/PV` 数值策略差异”为结论收束：暂不继续投入对齐 Fortran 的 PV 细节实现；以 Julia 侧当前实现作为后续工作的基线。

产物与复现入口：

- 结果文件：
  - `data/outputs/results/ak_compare.csv`
  - `data/outputs/results/meson_fortran_compare.csv`
  - `data/outputs/results/meson_full_fortran_params_experiment.csv`
- 相关脚本入口：
  - 已集中归档：
    - docs/dev/archived/2026-01-22_RelaxTime_Fortran_Comparison_Scripts.md

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
