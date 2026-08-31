# Diagnostic workflow retirement wave 1 v1

这是 `diagnostic_workflow_inventory_v1` 经作者于 2026-08-31 审核后的首批执行包。
本轮只退役作者明确批准的四个一次性 workflow，并对剩余纯手动 workflow 做第二轮
结构审计；不触发数值计算，不删除 GitHub run/artifact，也不改变 solver、reference 或
production 结果。

## 首批退役语义

“退役”在本包中表示：

1. YAML 从 `.github/workflows/` 移出，因此不再出现在 Actions 手动触发入口；
2. 原字节定义保存在 `definitions/`，并记录 source commit、Git blob SHA 和 SHA-256；
3. workflow-specific tests 改为读取 versioned retired definition，继续验证历史合同；
4. 历史 run、artifact、分析包和失败状态保持不变，不把退役解释为历史成功或失败被抹除。

首批四个文件见 `retired_workflows.csv`。活动 workflow 数量由 47 降为 43，含
`workflow_dispatch` 的入口由 43 降为 39；其中纯手动入口由 29 降为 25。

## 第二轮手动入口审计

`manual_dispatch_inventory_v2.csv` 只统计剩余 25 个纯 `workflow_dispatch` workflow，
不把 CI/schedule/release workflow 的手动重跑入口计为冗余。结论分为：

- 8 个 durable production/recovery 入口继续保留；
- 8 个已闭环 Issue #130/legacy 专用入口进入下一批直接退役审核；
- 7 个入口进入“合并为参数化入口或整体退役”审核；
- 2 个 CEP oracle/pilot 入口保留为历史重放审核，待确认未来是否仍需复用。

结构相似度最高的三组见 `consolidation_groups.csv`。相似度只证明 YAML/job/script
结构重复，不证明物理合同可以无条件共用；任何合并都必须保留 scope、输入 SHA、failed-only、
artifact schema 和 solver_called 语义。

## 后续边界

本 PR 不退役第二轮候选。首批合并后，作者可按组选择：

- 对已完成的一次性诊断直接做 wave 2 退役；或
- 对确有未来复用价值的组，新建一个 versioned 参数化 workflow，再迁移测试和入口。

在选择前，不修改剩余 25 个手动 workflow。
