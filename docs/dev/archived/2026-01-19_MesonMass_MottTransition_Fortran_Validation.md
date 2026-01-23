---
title: MesonMass/MottTransition Fortran 对比验证
archived: true
original: docs/dev/active/2026_1_19新增功能测试.md
archived_date: 2026-01-19
---

以下为原始内容（保留，以便审阅与历史参考）：

***

1.添加按 docs/reference/formula/relaxtime/propagator/MesonMass_RPA_Pole.md 中记录的完整流程计算介子质量和 Mott 相变位置的扫描脚本。

已实现（2026-01-21）：

- 工作流模块：src/pnjl/workflows/MesonMassWorkflow.jl
	- 负责“PNJL 平衡求解 → (quark_params, thermo_params) → 介子质量/宽度 → Mott 阈值/gap”的串联。
	- 放在 pnjl/workflows 下的理由：该链路需要 PNJL.solve 提供 (m_u,m_s,Φ,Φbar) 才能算到完整流程。
- 扫描脚本：scripts/relaxtime/run_gap_meson_mass_scan.jl
	- 扫描 (T, μ_B, ξ) 网格并输出 CSV（每点给出阈值与 gap；Mott 点可通过 gap 过零在后处理中定位）。

2.正确性检查：与 Fortran 介子质量输出对照（不要求完全一致，但应大致相同）

已完成（2026-01-21）：

- 对照脚本（已从 `scripts/relaxtime` 清理；代码见归档）：
	- docs/dev/archived/2026-01-22_RelaxTime_Fortran_Comparison_Scripts.md
- 参考数据：D:/Desktop/fortran代码/输运系数/RelaxTime/PNJL-mott-mu_T/PNJL-mu-T/quark_phase/meson.dat
- 示例命令（μ_B=0, ξ=0）：

	# 从归档文档复制 compare_meson_masses_with_fortran.jl 到本地临时文件后运行
	julia --project=. D:/tmp/compare_meson_masses_with_fortran.jl --fortran-meson-file "D:/Desktop/fortran代码/输运系数/RelaxTime/PNJL-mott-mu_T/PNJL-mu-T/quark_phase/meson.dat" --t-list 150,200,210

- 结论（抽样点 T=150/200/210 MeV）：
	- π、K、σ(及相关标量通道)、σ′ 的质量与 Fortran 基本一致（通常 < 1%）。
	- η/η′ 在部分温度点可能出现明显差异：脚本会自动尝试“本征态交换匹配”。在 T=210 MeV 的抽样中，η′ 可对到 ~2% 量级，但 η 仍有 ~30% 的差异。
	- 该差异更可能来自混合通道选根/排序、数值积分与截断、以及 PV/B0 策略差异，而非实现错误；若需要进一步收敛诊断，可增加对该点的输出（混合角、残差、阈值/gap、初值重试路径）并与 Fortran 同步记录。

***

归档补充（2026-01-22）：

- 上述对照脚本依赖外部 Fortran 工作区（非仓库默认依赖），已从 `scripts/relaxtime` 清理；脚本代码集中归档在：
	- docs/dev/archived/2026-01-22_RelaxTime_Fortran_Comparison_Scripts.md
- API 文档已更新为引用归档入口，而非脚本路径。

复现提示：

- 若需复现对照，请从脚本归档文档中复制 `compare_meson_masses_with_fortran.jl` 到本地临时文件再运行（避免把外部依赖脚本长期留在仓库默认脚本集合中）。

补充结论（早期诊断摘要，保留）：

- A 函数收敛性测试通过，仅能说明默认上限与节点数在代表点足够稳定，但不足以单独证明介子质量差异来源。
- B0 差异为主要来源：通过 Fortran B0 输出与 Julia B0 对比、以及 B0 注入实验可见残差差异显著。
- QuadGK 交叉验证结果显示 Julia B0 与独立数值方法的偏差显著小于 Fortran B0，支持 Julia B0 更接近数值收敛极限。
- 介子质量差异仍可能包含混合 B0（负 λ 缺失时回退）、以及 A/K 系数与求根策略等因素。