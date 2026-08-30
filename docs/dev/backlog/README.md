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
- [带电 K/π 介子数密度物理约束与验证中期路线](2026-07-25_带电KPi介子数密度物理约束与验证中期路线.md)

## 使用规则

- 命名建议沿用 `YYYY-MM-DD_描述.md`
- 当 backlog 条目进入实际执行时，应拆成新的 `docs/dev/active/` 任务单
- backlog 文档不替代实现验收；真正执行时仍需有独立 active 任务单承接
