# Diagnostic workflow retirement wave 2 review v1

这是 wave 1 合并后的第二轮纯手动 workflow 审阅包。源树固定为
`main@15e02b8e8f9b23629908b6ee164d1e858f7ac463`。本包只读取 active YAML、已有
inventory/evidence 和仓库引用，不改变 `.github/workflows/`，不删除 run/artifact，不调用
PNJL solver，也不改变数值结果。

## 审阅范围

wave 1 后共有 25 个纯 `workflow_dispatch` 入口。本包把其中 17 个需要作者判断的入口拆成
三类：

- **直接退役候选（8）**：已闭环的 Issue #130/CEP/legacy 一次性入口，历史证据已有替代物。
- **合并或退役候选（7）**：后端脚本/job 结构重复，但仍需判断是否保留一个参数化入口。
- **历史重放审核（2）**：昂贵 CEP oracle/pilot，除非未来有明确复用需求，否则不应长期保留
  active trigger。

没有放入本轮候选的 8 个入口继续保留：dense-reference 主/replay/resume、Zenodo
archive/restore，以及 relaxtime manual/meson/phase-guided production。

## 如何阅读

`manual_workflow_function_review.csv` 对每个入口记录：

1. 它真正执行的 numerical、aggregate 或 replay 功能；
2. 使用的 Julia/Python 后端脚本和输入范围；
3. 可退役的具体理由；
4. 退役或合并前必须保留的输入 SHA、failed-only、artifact schema 和 provenance。

`consolidation_groups.csv` 来自 wave 1 的结构比较。相似度只表示 YAML/job/script 重复，
不表示物理合同等价。任何参数化合并都必须显式保留 scope、source/calculation SHA、
failed-only、oracle/diagnostic 语义和 `solver_called` 字段。

## 作者审核问题

请逐组选择：

- 直接退役并保留 versioned definition；
- 保留一个新的参数化 workflow，再迁移测试和入口；
- 暂时保留，因为存在明确的未来重放需求。

本包不执行第二波 workflow 迁移。作者批准后，另立精确 allowlist 的 retirement 或
parameterization PR。

## 计数与边界

- active workflow：43；带 `workflow_dispatch`：39；纯手动：25。
- 本包审阅：17；明确保留：8。
- `solver_called=false`；`numerical_results_changed=false`；`wave2_changes_executed=false`。
