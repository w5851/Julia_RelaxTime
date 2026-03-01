---
title: models/pnjl/compat 同名替代与职责收敛开发计划
archived: true
original: docs/dev/active/2026-02-27_models_pnjl_compat同名替代与职责收敛开发计划.md
archived_date: 2026-03-01
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# models/pnjl/compat 同名替代与职责收敛开发计划

更新日期：2026-02-27

> 目标：对 `src/models/pnjl/compat/**` 逐文件检查，优先按“同名文件替代”完成收敛；若无同名文件，再按“模块职责等价”完成替代与路径切换，最终缩减 compat 表面面积并保持回归全绿。

---

## 1. 背景与目标

- [x] 目标 A：完成 `src/models/pnjl/compat/**` 全量文件盘点（文件级）。
- [x] 目标 B：按“同名优先”完成替代决策矩阵（保留/替代/删除）。
- [x] 目标 C：对无同名文件执行“职责等价”匹配并完成替代。
- [x] 目标 D：替代后通过 smoke + migration/data-output/prune-wave 门禁。

---

## 2. 范围与非范围

### 2.1 范围（In Scope）

- [x] `src/models/pnjl/compat/**` 全目录（含 `core/derivatives/scans/solver/workflows/analysis`）。
- [x] `src/models/**` 中作为替代目标的同名/等价实现文件。
- [x] 路径常量、include 入口、allowlist、相关测试引用的同步收敛。

### 2.2 非范围（Out of Scope）

- [ ] 不新增业务功能，仅做实现收敛与替代。
- [ ] 不做与本次替代无关的大规模数值重构。

---

## 3. 执行策略（用户确认口径）

- [x] 规则 R1（主规则）：优先按“文件名同名”判断可替代关系。
- [x] 规则 R2（兜底规则）：若未找到同名文件，再按“模块职责等价”判断。
- [x] 规则 R3（安全规则）：替代以最小改动为原则，先替路径/入口，再删冗余实现。

---

## 4. 任务分解（可勾选）

### Batch C1：同名清单与替代矩阵

- [x] 输出 `compat` 文件清单（含子目录）。
- [x] 枚举 `src/models/**` 同名文件并建立映射表。
- [x] 将每个文件标记为：`同名可替代` / `无同名待职责匹配` / `保留`。

验收：
- [x] 形成文件级映射清单（可追溯到路径）。

### Batch C2：同名替代落地

- [x] 对“同名可替代”项完成入口/路径切换。
- [x] 删除或停用已被替代的 compat 文件（按最小风险分批）。
- [x] 修正受影响测试与 allowlist 路径。

验收：
- [x] `pnjl_compat`（或 `pnjl/compat`）在源码中的主路径引用显著下降。

### Batch C3：职责等价替代（无同名项）

- [x] 为每个“无同名”文件确定职责等价目标（函数/模块级）。
- [x] 完成迁移或转发收敛（保证 API/调用签名兼容优先）。
- [x] 逐项补充最小回归证据（单点/小网格）。

验收：
- [x] 形成“职责等价替代清单 + 证据”。

### Batch C4：收口与清理

- [x] 清理残留兼容路径常量与历史注释口径。
- [x] 同步文档口径（开发计划、执行台账）。
- [x] 评估 `src/models/pnjl/compat` 是否可继续裁剪。

验收：
- [x] 目录结构与 `models` 主分层语义一致（`pnjl` 子域优先）。

---

## 5. 测试与验收标准

### 5.1 必过

- [x] `UNIT_PROFILE=smoke julia --project=. tests/unit/runtests.jl`
- [x] `julia --project=. scripts/dev/check_pnjl_migration_guard.jl --base HEAD --head HEAD`
- [x] `julia --project=. scripts/dev/check_data_output_path_guard.jl --base HEAD --head HEAD`
- [x] `julia --project=. scripts/dev/run_prune_wave_gate.jl --base HEAD --head HEAD`

### 5.2 必稳

- [x] 迁移后 `models-invokelatest-audit` 与 allowlist 对齐。
- [x] 不引入新的 legacy 反向依赖。

### 5.3 可追溯

- [x] 每批次在执行台账追加记录（命令、产物、结果）。
- [x] 关键门禁产物落档 `data/outputs/results/*`。

---

## 6. 里程碑

- [x] M-C1：完成同名映射清单。
- [x] M-C2：完成同名替代落地。
- [x] M-C3：完成职责等价替代。
- [x] M-C4：完成收口清理并门禁全绿。

---

## 7. DoD（完成定义）

- [x] `src/models/pnjl/compat` 中重复实现按规则完成替代或有明确保留理由。
- [x] 主路径优先使用 `src/models/**` 现有实现。
- [x] smoke + 三门禁通过，且台账与产物完备。

---

## 9. 执行摘要（2026-02-27 Batch C1~C4）

- C1 映射结论：`compat` 共 27 个 `.jl` 文件；其中同名可替代 12 个（`scans/*` 7 个、`workflows/*` 3 个、`core/*Facade` 2 个），无同名 15 个。
- C2 落地：
  - `src/models/pnjl/compat/PNJL.jl` 的扫描模块 include 已切到 `src/models/scans/*` 同名实现。
  - `tests/unit/models/test_pnjl_thermo_bridge_multipoint_smoke.jl` 已切到 `src/models/core/ThermoFacade.jl`。
  - `src/models/legacy/LegacyPNJLModel.jl` 中 `ThermoFacade` 路径已切到 `src/models/core/ThermoFacade.jl`。
- C3 落地（职责等价）：新增 `src/models/pnjl/PNJL.jl` 作为主入口转发，`PNJLModel`、`ScanEntrypointsBridge`、`ThermoDerivativesBridge`、`LegacyPNJLModel` 改为优先指向该入口。
- C4 验证：smoke + migration/data-output/prune-wave 全部通过。

---

## 8. 风险与回退

- [ ] 风险 R1：职责误判导致行为漂移。
  - 回退：按 Batch 粒度回滚；保留映射清单与门禁证据。
- [ ] 风险 R2：一次性删除过多导致连锁失败。
  - 回退：分批替代、分批验证，不做大爆炸清理。
