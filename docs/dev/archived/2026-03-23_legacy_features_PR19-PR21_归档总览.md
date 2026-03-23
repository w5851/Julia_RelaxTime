---
title: Legacy 功能吸入（三个 PR）归档总览
archived: true
original: docs/dev/active/2026-03-23_legacy_features_PR19-PR21_归档总览.md
archived_date: 2026-03-23
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Legacy 功能吸入（三个 PR）归档总览

更新日期：2026-03-23

## 1. 目标与范围

本总览用于补齐“前三个 PR（#19/#20/#21）”的开发过程记录，统一沉淀：

- 目标与拆分策略（先文档、后模型壳）
- 每个 PR 的核心变更与边界
- 合并前验证与冲突处理记录
- 当前主线状态与后续衔接计划

本总览对应的三个 PR：

- PR #19: https://github.com/w5851/Julia_RelaxTime/pull/19
- PR #20: https://github.com/w5851/Julia_RelaxTime/pull/20
- PR #21: https://github.com/w5851/Julia_RelaxTime/pull/21

## 2. 拆分策略（执行口径）

- PR1（#19）仅做“公式/方法文档内化”，不夹带模型实现。
- PR2（#20）仅做 Gas-Liquid 最小核心与统一入口接入。
- PR3（#21）仅做 Rotation 适配壳与最小 workflow；随后在主线合并时处理与 PR2 的入口冲突。

该拆分目标是降低评审耦合，确保每个 PR 都可独立验证和回滚。

## 3. PR 逐项归档

### 3.1 PR #19（文档内化）

标题：`docs: internalize gas-liquid/rotation formula foundations and MFA Omega workflow`

核心产物：

- `docs/reference/formula/models/gas_liquid/GasLiquid_RMF_CoreEquations.md`
- `docs/reference/formula/models/rotation/Rotation_PNJL_CoreEquations.md`
- `docs/reference/formula/models/shared/Lagrangian_to_GrandPotential_MFA_Workflow.md`
- `docs/reference/formula/models/README.md`

目的：

- 先把 legacy 理论内容转化为主线可实现语义。
- 固化“从拉氏量到巨热力学势”的统一推导路径。

状态：已合并（MERGED）。

### 3.2 PR #20（Gas-Liquid 最小核心）

标题：`feat(models): add gas-liquid minimal core with unified Models workflow entry`

核心产物：

- `src/models/variants/gas_liquid/GasLiquidModel.jl`
- `src/models/variants/gas_liquid/core/EquationSet.jl`
- `src/models/variants/gas_liquid/core/Thermodynamics.jl`
- `src/models/variants/gas_liquid/workflows/GasLiquidWorkflow.jl`
- `tests/unit/models/test_gas_liquid_model.jl`
- `tests/unit/models/test_gas_liquid_workflow.jl`
- `tests/integration/models/test_gas_liquid_workflow_smoke.jl`

目的：

- 在 `Models` 统一入口下接入 Gas-Liquid 的最小可运行壳。
- 先打通模型-工作流-测试闭环，再逐步补物理细节。

状态：已合并（MERGED）。

### 3.3 PR #21（Rotation 适配壳）

标题：`feat(models): add rotation adapter shell with T-mu-omega workflow entry`

核心产物：

- `src/models/variants/rotation/RotationModel.jl`
- `src/models/variants/rotation/core/RotationThermo.jl`
- `src/models/variants/rotation/workflows/RotationWorkflow.jl`
- `tests/unit/models/test_rotation_model.jl`
- `tests/unit/models/test_rotation_workflow.jl`
- `tests/integration/models/test_rotation_workflow_smoke.jl`

冲突与处置：

- 与 PR #20 在 `Models` 汇聚入口处发生合并冲突（`Models.jl` / `entrypoints.jl` / `factory.jl` / `test_factory.jl`）。
- 采用人工合并，保留 `GasLiquid` 与 `Rotation` 双方导出和 workflow 路由。

状态：已合并（MERGED）。

## 4. 关键验证记录

在冲突修复与合并前后，执行并通过以下门禁：

- Unit smoke（`UNIT_PROFILE=smoke`）
- Integration smoke（`INTEGRATION_PROFILE=smoke`）
- 迁移守卫（`scripts/dev/check_pnjl_migration_guard.jl`）

补充：文档治理检查 `scripts/dev/check_docs_consistency.jl` 已在后续文献对齐提交中通过。

## 5. 后续承接（已在主线推进）

- 依据 `2307.14402v1.pdf` 对 Rotation 公式文档做了文献对齐增强。
- 依据 `thesis.pdf` 第 2.2 节与 `tmp/公式文档temp.md` 对 Gas-Liquid 公式文档做了逐式对齐、符号映射与参数口径收敛（`f_sigma/f_omega/f_rho`、`rho_B=rho_p+rho_n`）。

这些后续更新已进入 `main`，用于为下一步参数模板与实现增强提供一致依据。

## 6. 结论

前三个 PR 的主线目标已完成：

- 文档先行内化（PR #19）
- Gas-Liquid 最小核心接入（PR #20）
- Rotation 适配壳接入并完成冲突收口（PR #21）

后续工作应转入“参数模板落地 + 文档与实现一致性校核”。
