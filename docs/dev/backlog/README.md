# backlog 目录说明

本目录用于存放长期路线图、候选需求盘点、能力清单与未进入当前执行批次的规划文档。

## 适用场景

- 长期滚动维护的 roadmap
- 功能候选池与优先级盘点
- 暂未拆成 active 任务单的规划文档

## 不适用场景

- 正在执行、需要持续回填状态的任务单
- 已完成且应归档的任务文档

## 与 active / archived 的关系

- `active`：当前正在推进的任务
- `backlog`：尚未正式拉起执行的规划与候选项
- `archived`：已完成任务的历史归档

## 当前规划文档

- [PNJL 可选功能盘点与优先级任务单](2026-03-01_PNJL可选功能盘点与优先级任务单.md)
- [sysimage 产品化与 package/app/service 演进任务单](2026-05-04_sysimage产品化与package-app-service演进任务单.md)
- [介子热力学 regression 与 validation 治理方案](2026-05-09_介子热力学regression与validation治理方案.md)
- [PNJL 全温区参考生产 PR3 规划](2026-07-17_PNJL全温区相图热积分与参考生产.md)
- [带电 K/π 介子数密度物理约束与验证中期路线](2026-07-25_带电KPi介子数密度物理约束与验证中期路线.md)

## 当前恢复入口（2026-11-10 核对）

介子冻结线数值、更新实验图件、矢量交付及论文接入均已完成；
[Arena Issue #313](https://github.com/w5851/Julia_RelaxTime/issues/313) 已按作者接受的现行定义关闭。
论文后续文字审阅由 `D:/Desktop/paper` 管理。PALC 的候选搜索优化、通用生产接入和预算对照
已由 [PR #338](https://github.com/w5851/Julia_RelaxTime/pull/338)、
[#339](https://github.com/w5851/Julia_RelaxTime/pull/339)、
[#340](https://github.com/w5851/Julia_RelaxTime/pull/340) 完成；
任务合同已[归档](../archived/2026-08-18_PALC同口径性能与SeedContinuation一致性诊断.md)。
较长的标准 PNJL Tρ/Tμ 扫描默认 PALC；Tμ 完整旧策略/warm 对照和 RPNJL 默认资格
留待后续实际需要时再安排，当前没有新增计算任务。

| 待评估研究 | 实际剩余工作 | 当前状态 |
| --- | --- | --- |
| 全温区参考生产 PR3 | 确定新生产温区与输入，完成积分/相线/oracle 收敛及新参考产物验收 | `deferred`；算法 PR #137/#139 已完成，生产尚未开展 |

其余文档是按需求重开的候选研究与能力规划。旧阶段勾选保留当时的语境，
实际进度由各页当前摘要、现有实现及 `config/governance/task_tracks.toml` 确定。
P1/BU2020 的条件性扩展见[候选清单](2026-08-03_P1与BU2020后续治理待办.md)。

## 使用规则

- 命名建议沿用 `YYYY-MM-DD_描述.md`
- 启动时先明确剩余范围与验收条件，通常使用对话计划或已有 issue。
- 需要共享设计、用户要求文档或现有 ledger 合同时，再维护 `active/` 文档；参见[任务跟踪治理](../task_tracking_governance.md)。
