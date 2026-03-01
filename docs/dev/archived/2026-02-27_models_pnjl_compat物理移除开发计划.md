---
title: models/pnjl/compat 物理移除开发计划
archived: true
original: docs/dev/active/2026-02-27_models_pnjl_compat物理移除开发计划.md
archived_date: 2026-03-01
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# models/pnjl/compat 物理移除开发计划

更新日期：2026-02-27

> 最终目标：在保证回归门禁全绿的前提下，**物理移除** `src/models/pnjl/compat` 实体目录，并确保主路径仅依赖 `src/models/**`。

---

## 1. 背景与目标

- [x] 目标 A：完成 `src/models/pnjl/compat` 残余能力盘点与迁移落位映射。
- [x] 目标 B：将 compat 残余能力迁移/收敛到 `src/models/**` 主域。
- [x] 目标 C：断开所有源码与配置对 `pnjl/compat` 的路径依赖。
- [x] 目标 D：物理删除 `src/models/pnjl/compat` 后通过四重门禁。

---

## 2. 范围与非范围

### 2.1 范围（In Scope）

- [x] `src/models/pnjl/compat/**` 全目录能力迁移与删除。
- [x] 入口与桥接路径收敛（`src/models/pnjl`, `src/models/core`, `src/models/scans`, `src/models/derivatives`, `src/models/legacy`）。
- [x] `config/ci/models_invokelatest_allowlist.toml` 与迁移后实际路径对齐。
- [x] 相关单测/门禁脚本校验与台账记录。

### 2.2 非范围（Out of Scope）

- [ ] 不新增业务功能，不改变对外数值语义。
- [ ] 不做与本次物理移除无关的大规模重构。
- [ ] 不扩展新 UI/新 API，仅做路径与实现收敛。

---

## 3. 现状盘点（基于已完成分析）

### 3.1 现存实体

- [ ] `src/models/pnjl/compat` 仍有 27 个 `.jl` 文件（待移除）。

### 3.2 已识别硬依赖（需断开）

- [ ] `src/models/pnjl/PNJL.jl` 仍 include `compat/PNJL.jl`。
- [ ] `src/models/pnjl/PNJLMagneticModel.jl` 仍指向 `compat/core/MagneticThermodynamics.jl`。
- [ ] `src/models/core/ThermoFacade.jl` 仍 include `pnjl/compat/core/ThermoFacade.jl`。
- [ ] `src/models/core/EquilibriumFacade.jl` 仍 include `pnjl/compat/core/EquilibriumFacade.jl`。
- [ ] `config/ci/models_invokelatest_allowlist.toml` 仍含 `src/models/pnjl/compat/derivatives/ThermoDerivatives.jl` 键。

### 3.3 无同名替代清单（职责迁移主工作量）

- [ ] analysis: `PhaseTransition.jl`
- [ ] core: `Core.jl`, `Integrals.jl`, `MagneticIntegrals.jl`, `MagneticThermodynamics.jl`, `ModelThermodynamics.jl`, `Thermodynamics.jl`
- [ ] derivatives: `ThermoDerivatives.jl`
- [ ] solver: `ConstraintModes.jl`, `SeedStrategies.jl`, `Conditions.jl`, `ImplicitSolver.jl`, `Solver.jl`
- [ ] entry: `DualBranchScanEntry.jl`

---

## 4. 执行策略

- [x] 规则 R1：先“迁移能力并切主路径”，后“物理删除 compat”。
- [x] 规则 R2：优先最小改动，保持 API/调用签名兼容。
- [x] 规则 R3：按批次推进，每批次落地后立即跑局部验证。
- [x] 规则 R4：删除前必须满足“全仓无 compat 源码/配置依赖”。

---

## 5. 任务分解（可勾选）

### Batch D1：残余能力迁移落位设计与占位

- [x] 为 `compat/core/*` 建立主域落位（建议：`src/models/pnjl/core/*`）。
- [x] 为 `compat/solver/*` 建立主域落位（建议：`src/models/pnjl/solver/*`）。
- [x] 为 `compat/derivatives/*` 建立主域落位（建议：`src/models/pnjl/derivatives/*`）。
- [x] 为 `compat/analysis/*` 建立主域落位（建议：`src/models/pnjl/analysis/*`）。
- [x] 明确 `DualBranchScanEntry.jl` 的去向（并入 `src/models/pnjl` 入口或删除）。

验收：
- [x] 输出“文件级迁移映射表（旧路径 -> 新路径）”。

### Batch D2：核心模块迁移与入口切换

- [x] 迁移 `core/*` 到目标主域并修正 include 相对路径。
- [x] 迁移 `solver/*` 与 `derivatives/*` 到目标主域并修正依赖。
- [x] 迁移 `analysis/*`，保证 `PNJL` 主入口可加载。
- [x] 逐步将 `src/models/core/*Facade.jl`、`PNJLMagneticModel.jl` 指向新主域实现。

验收：
- [x] `src/models/**` 主路径可在不依赖 `compat/*` 的情况下完成加载与求解主链路。

### Batch D3：桥接/配置/测试收口

- [x] 将 `src/models/pnjl/PNJL.jl` 从转发 compat 改为直接主域组装。
- [x] 清理桥接层与 legacy 适配器中残余 compat 路径。
- [x] 更新 `models_invokelatest_allowlist.toml` 键到新路径（或移除不再需要项）。
- [x] 修正受影响测试引用路径。

验收：
- [x] `grep` 全仓无 `pnjl/compat` 路径命中（文档历史除外）。

### Batch D4：物理删除与终验

- [x] 物理删除 `src/models/pnjl/compat` 目录。
- [x] 运行 smoke + migration/data-output/prune-wave 四重门禁。
- [x] 记录产物并更新开发文档与执行台账。

验收：
- [x] 目录删除后主流程与 CI 门禁全绿。

---

## 6. 测试与验收标准

### 6.1 必过

- [x] `UNIT_PROFILE=smoke julia --project=. tests/unit/runtests.jl`
- [x] `julia --project=. scripts/dev/check_pnjl_migration_guard.jl --base HEAD --head HEAD`
- [x] `julia --project=. scripts/dev/check_data_output_path_guard.jl --base HEAD --head HEAD`
- [x] `julia --project=. scripts/dev/run_prune_wave_gate.jl --base HEAD --head HEAD`

### 6.2 必稳

- [x] migration guard 中 models-invokelatest-audit 与 allowlist 对齐。
- [x] 不新增任何 `src/models -> compat` 反向依赖。

### 6.3 可追溯

- [x] 每个批次在执行台账 append-only 追加命令、产物、结果。
- [x] 关键产物落档到 `data/outputs/results/*`。

---

## 7. 里程碑

- [x] M-D1：完成迁移映射与主域落位方案。
- [x] M-D2：完成 core/solver/derivatives/analysis 迁移并切换主入口。
- [x] M-D3：完成桥接与配置收口，源码无 compat 依赖。
- [x] M-D4：完成目录物理删除并门禁全绿。

---

## 8. DoD（完成定义）

- [x] `src/models/pnjl/compat` 已物理删除。
- [x] 全仓源码与配置中不再依赖 `pnjl/compat` 路径（文档历史除外）。
- [x] smoke + 三门禁全部通过，且台账与产物完备。

---

## 10. 执行摘要（2026-02-27）

- 已完成实体迁移：`compat/{core,solver,derivatives,analysis}` 与 `PNJL.jl`、`DualBranchScanEntry.jl` 全部迁入 `src/models/pnjl/**`。
- 已完成依赖收口：`src/models/core/{ThermoFacade,EquilibriumFacade}.jl`、`src/models/pnjl/PNJLMagneticModel.jl`、`config/ci/models_invokelatest_allowlist.toml` 已切到新路径。
- 已完成物理删除：`src/models/pnjl/compat` 目录已删除，源码/测试/配置中 `pnjl/compat` 路径命中为 0。
- 迁移中修复两处路径层级问题：`src/models/pnjl/core/Thermodynamics.jl` 与 `src/models/pnjl/PNJL.jl` 的相对路径已修正。
- 回归结果：smoke 845/845，migration/data-output/prune-wave 全绿。

---

## 9. 风险与回退

- [ ] 风险 R1：跨模块迁移引发 include 层级错误。
  - 回退：按 Batch 粒度回滚；保留迁移映射表逐项恢复。
- [ ] 风险 R2：求解器/导数模块迁移后出现 world-age 或重复加载问题。
  - 回退：优先统一 IncludeOnce 入口，必要时临时保留桥接壳层。
- [ ] 风险 R3：allowlist 与审计口径短期不一致导致门禁失败。
  - 回退：先更新 allowlist 到新路径，再做最终删除。
