---
title: PNJL Crossover Legacy Dual-Source Mini Family
archived: true
original: docs/dev/active/2026-03-09_pnjl_crossover_legacy_dual_source_mini_family.md
archived_date: 2026-03-09
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# PNJL Crossover Legacy Dual-Source Mini Family

更新日期：2026-03-09

## 0. 任务定位

本文档把当前 PNJL crossover legacy 双源工作收敛成一个正式的小 family 说明，目标不是再扩大 legacy 扫描覆盖面，而是固定三件事：

- 当前 dual-source crossover family 的组成
- acceptance / evidence 的准入顺序
- 进入收尾前应保持不再回退的约束

## 1. 当前 mini family 组成

当前 `tests/validation/pnjl/reference/` 下的 crossover family 由三层组成：

### 1.1 Reference Acceptance

- target: `tests/validation/data/targets/pnjl/reference/pnjl_crossover_reference_targets_v1.csv`
- test: `tests/validation/pnjl/reference/test_crossover_reference_targets.jl`
- 语义：以仓库内 `data/reference/pnjl/crossover.csv` 提取出的 Julia reference fixed points 为 acceptance 主锚点。

### 1.2 Legacy Dual-Source Guardrail

- target: `tests/validation/data/targets/pnjl/reference/pnjl_crossover_legacy_dual_source_targets_v1.csv`
- test: `tests/validation/pnjl/reference/test_crossover_legacy_dual_source_guardrail.jl`
- 语义：在 legacy C++ 与 legacy Fortran 都能给出 machine-readable fixed points 的前提下，用共享固定点约束 Julia `detect_crossover(...)` 的 branch-safe 输出窗口。

### 1.3 Legacy Source-Consistency Gate

- target: `tests/validation/data/targets/pnjl/reference/pnjl_crossover_legacy_source_consistency_targets_v1.csv`
- test: `tests/validation/pnjl/reference/test_crossover_legacy_source_consistency.jl`
- 语义：该层不是再测 Julia，而是先回答 legacy Fortran 与 legacy C++ 是否彼此一致。只有通过这道 gate 的共享点，才有资格进入 dual-source acceptance guardrail。

## 2. 当前共享 fixed points

当前已正式承认的 shared legacy fixed points 只有 4 个：

- `muB = 2 MeV`，Fortran `T_c = 200.0 MeV`，C++ `T_c = 200.1 MeV`，`|ΔT| = 0.1 MeV`
- `muB = 105 MeV`，Fortran `T_c = 199.0 MeV`，C++ `T_c = 199.1 MeV`，`|ΔT| = 0.1 MeV`
- `muB = 209 MeV`，Fortran `T_c = 196.0 MeV`，C++ `T_c = 196.4 MeV`，`|ΔT| = 0.4 MeV`
- `muB = 291 MeV`，Fortran `T_c = 193.0 MeV`，C++ `T_c = 192.9 MeV`，`|ΔT| = 0.1 MeV`

当前 source-consistency 约束固定为：

- `|ΔT| <= 0.5 MeV`

这 4 个点目前已经满足“共享来源 + 明确字段语义 + source-consistent”三项最低要求，因此可以保留在 dual-source guardrail 中。

## 3. 准入规则

### 3.1 进入 evidence 的条件

以下对象先进入 `provenance/.../evidence/`，而不是直接变成 acceptance target：

- 新发现的 legacy 导出文件
- 只存在单源、尚无第二来源交叉核对的 fixed points
- 需要额外解释字段语义、横轴换算或单位链路的 legacy 结果
- 涉及明显 branch divergence 或热力学稳定性不清的 comparison point

### 3.2 进入 source-consistency gate 的条件

只有同时满足以下条件，legacy crossover 点才进入 source-consistency target：

- 至少两个 legacy 来源可 machine-readable 读取
- 横轴与物理字段已经统一到相同语义
- shared fixed points 不依赖人工读图或手工脚本推断
- 它们的比较目标是“源间一致性”，而不是“强迫 Julia 贴某一旧分支”

### 3.3 进入 acceptance 的条件

只有同时满足以下条件，legacy crossover 点才进入 dual-source guardrail acceptance：

- 已通过 source-consistency gate
- 对应 comparison point 是 branch-unique 或至少没有现成 branch ambiguity 证据
- Julia 侧 quantity 定义与 legacy 侧字段语义一一对应
- 不需要整条 raw scan 才能解释测试含义

## 4. 明确排除项

以下对象当前不应进入 crossover acceptance：

- 原始 legacy `T-mu` 整条扫描曲线
- 仅靠单个 legacy 源存在的 fixed points
- 高温低密区里存在 branch mismatch 证据的 legacy 点
- 需要先比较不同 seed/stability 才能解释的点

这条规则与本轮对 `lorentz_legacy` 的处理保持一致：

- 数值可导出，不等于可作为 acceptance anchor
- branch ambiguity 一旦进入主 compare，默认测试就会退化成“比较分支选择”，而不是比较稳定物理量

## 5. 当前 family 的完成判断

以 2026-03-09 当前状态看，这个 mini family 已满足 archive-ready 的核心条件：

- reference acceptance 已独立落地
- legacy dual-source acceptance 已只保留 shared fixed points
- source-consistency gate 已独立落地
- acceptance / evidence 边界已不再依赖口头说明

仍然保留为 follow-up 的事项只有两类：

- 若未来出现新的 shared fixed points，可按同一准入顺序扩容
- 若未来出现 branch-safe 的更密 fixed-point 子集，可再新增 target 版本，而不是回退到 raw scan compare

## 6. 收尾建议

当前 crossover legacy 双源工作不再需要继续扩大范围。更合理的收尾方式是：

1. 保持现有 3 层结构不再混并
2. 以后新增点先过 source-consistency，再谈 acceptance
3. 不再把 branch-ambiguous legacy 点当作“还差一点就能收紧”的默认 backlog

在此约束下，本 mini family 已可视为完成 formalization，后续只需按版本化 target 继续增量维护。