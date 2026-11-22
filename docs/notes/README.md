# 历史笔记（已退役，只读）

**状态：`retired_readonly`。不再新增文件。** 原 18 份未归类笔记已按 2026-11 的迁移方案分流：13 份迁入长期域（正文零改动，git 记录为 100% rename），5 份因审计路径锁定或已被现行规范取代而原地只读保留。本目录只作历史追溯，不作为现行权威依据。

迁移规则：只移动路径、不改正文语义；需改正文时另起提交并说明依据。

## 已迁出（13 份）

| 原路径 | 主题 | 现址 |
| --- | --- | --- |
| `forwarddiff/ForwardDiff_nested_dual_issue.md` | AD 嵌套双数 | [reference/domain-knowledge/ForwardDiff_nested_dual_issue.md](../reference/domain-knowledge/ForwardDiff_nested_dual_issue.md) |
| `forwarddiff/ForwardDiff_NLsolve_issue.md` | AD 与 NLsolve 交互 | [reference/domain-knowledge/ForwardDiff_NLsolve_issue.md](../reference/domain-knowledge/ForwardDiff_NLsolve_issue.md) |
| `integrand_method/对数奇点(B0)_奇点减法积分.md` | 奇点减法积分 | [reference/domain-knowledge/对数奇点(B0)_奇点减法积分.md](<../reference/domain-knowledge/对数奇点(B0)_奇点减法积分.md>) |
| `integrand_method/含奇点积分处理方法参考.md` | 奇点积分方法 | [reference/domain-knowledge/含奇点积分处理方法参考.md](../reference/domain-knowledge/含奇点积分处理方法参考.md) |
| `integrand_method/稀疏网格与张量积对比_average_rate.md` | 积分网格对比 | [reference/domain-knowledge/稀疏网格与张量积对比_average_rate.md](../reference/domain-knowledge/稀疏网格与张量积对比_average_rate.md) |
| `integrand_method/B0_k_positive_积分方法.md` | B0 正 k 积分 | [reference/domain-knowledge/B0_k_positive_积分方法.md](../reference/domain-knowledge/B0_k_positive_积分方法.md) |
| `integrand_method/cdf_sampling_for_expensive_functions.md` | 采样策略 | [reference/domain-knowledge/cdf_sampling_for_expensive_functions.md](../reference/domain-knowledge/cdf_sampling_for_expensive_functions.md) |
| `pnjl/2026-01-22_PNJL二阶导符号现象说明.md` | 二阶导符号现象 | [reference/domain-knowledge/2026-01-22_PNJL二阶导符号现象说明.md](../reference/domain-knowledge/2026-01-22_PNJL二阶导符号现象说明.md) |
| `pnjl/一阶相变双分支扫描策略.md` | 双分支扫描 | [reference/domain-knowledge/一阶相变双分支扫描策略.md](../reference/domain-knowledge/一阶相变双分支扫描策略.md) |
| `relaxtime/散射截面在积分下限下的发散以及积分的可积性.md` | 积分可积性 | [reference/domain-knowledge/散射截面在积分下限下的发散以及积分的可积性.md](../reference/domain-knowledge/散射截面在积分下限下的发散以及积分的可积性.md) |
| `relaxtime/2026-03-13_T190_T200_阈值邻域与峰值邻域补点复测记录.md` | 补点复测记录 | [analysis/historical/2026-03-13_T190_T200_阈值邻域与峰值邻域补点复测记录.md](../analysis/historical/2026-03-13_T190_T200_阈值邻域与峰值邻域补点复测记录.md) |
| `relaxtime/meson/2026-04-15-meson-debug-progress.md` | 调试进度 | [analysis/historical/2026-04-15-meson-debug-progress.md](../analysis/historical/2026-04-15-meson-debug-progress.md) |
| `structure/Julia模块配置文件设计指南.md` | 模块与配置设计 | [architecture/Julia模块配置文件设计指南.md](../architecture/Julia模块配置文件设计指南.md) |

## 原地只读（5 份）

| 文件 | 主题 | 保留原因 |
| --- | --- | --- |
| [pnjl/一阶相变方法对比.md](pnjl/一阶相变方法对比.md) | 一阶相变边界方法对比 | issue130 v1–v4 审计表把该路径写成历史契约，移动会改变证据路径与文件字节 |
| [pnjl/PNJL计算工作流.md](pnjl/PNJL计算工作流.md) | 历史导航 | 已被现行架构总览与指南取代 |
| [structure/结构体使用的必要性分析.md](structure/结构体使用的必要性分析.md) | 结构体迁移论证 | 已被参数容器契约与迁移指南取代 |
| [structure/结构体使用指南.md](structure/结构体使用指南.md) | 重构前问题分析 | 已被参数容器契约取代 |
| [md/VScode以(md)为后缀的文档的注意事项.md](<md/VScode以(md)为后缀的文档的注意事项.md>) | 编辑器提示 | 无域归属，且不应写入只读归档区 |

## 维护约定

- 本目录**不接受**新文件；新笔记按 [AGENTS.md](../../AGENTS.md) 的文档落点判据直接写进对应长期域。
- 「原地只读」的文件只作历史追溯，不得作为现行权威依据引用。
