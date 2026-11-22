# 历史笔记（已冻结）

**状态：`frozen`。不再新增文件。** 这里的 18 份笔记是早期未归类的推导、方法与排障记录，逐步按下方「计划去向」迁入长期域；迁空后本目录退役。

迁移规则：只移动路径、不改正文语义；每次迁移前先全仓 grep 引用（含 `tests/`、`scripts/`、`config/`、`.agents/`），并同步更新域索引。

| 文件 | 主题 | 计划去向 |
| --- | --- | --- |
| `forwarddiff/ForwardDiff_nested_dual_issue.md` | AD 嵌套双数问题 | `docs/reference/domain-knowledge/` |
| `forwarddiff/ForwardDiff_NLsolve_issue.md` | AD 与 NLsolve 交互 | `docs/reference/domain-knowledge/` |
| `integrand_method/对数奇点(B0)_奇点减法积分.md` | 奇点减法积分 | `docs/reference/domain-knowledge/` |
| `integrand_method/含奇点积分处理方法参考.md` | 奇点积分方法 | `docs/reference/domain-knowledge/` |
| `integrand_method/稀疏网格与张量积对比_average_rate.md` | 积分网格对比 | `docs/reference/domain-knowledge/` |
| `integrand_method/B0_k_positive_积分方法.md` | B0 正 k 积分 | `docs/reference/domain-knowledge/` |
| `integrand_method/cdf_sampling_for_expensive_functions.md` | 采样策略 | `docs/reference/domain-knowledge/` |
| `pnjl/2026-01-22_PNJL二阶导符号现象说明.md` | 二阶导符号现象 | `docs/reference/domain-knowledge/` |
| `pnjl/一阶相变方法对比.md` | 一阶相变方法 | `docs/reference/domain-knowledge/` |
| `pnjl/一阶相变双分支扫描策略.md` | 双分支扫描 | `docs/reference/domain-knowledge/` |
| `pnjl/PNJL计算工作流.md` | 工作流草图 | `docs/reference/domain-knowledge/` |
| `relaxtime/2026-03-13_T190_T200_阈值邻域与峰值邻域补点复测记录.md` | 补点复测记录 | `docs/analysis/historical/` |
| `relaxtime/散射截面在积分下限下的发散以及积分的可积性.md` | 积分可积性 | `docs/reference/domain-knowledge/` |
| `relaxtime/meson/2026-04-15-meson-debug-progress.md` | 调试进度 | `docs/analysis/historical/` |
| `structure/结构体使用的必要性分析.md` | 结构体与类型 | `docs/architecture/` |
| `structure/结构体使用指南.md` | 结构体使用 | `docs/architecture/` |
| `structure/Julia模块配置文件设计指南.md` | 模块与配置设计 | `docs/architecture/` |
| `md/VScode以(md)为后缀的文档的注意事项.md` | 编辑器工具笔记 | `docs/dev/archived/` |

## 维护约定

- 本目录**不接受**新文件；新笔记按 [落点判据](../../AGENTS.md) 直接写进对应长期域。
- 迁移时保持原文件名，避免破坏历史链接；确需改名时在提交说明里列出旧名。
