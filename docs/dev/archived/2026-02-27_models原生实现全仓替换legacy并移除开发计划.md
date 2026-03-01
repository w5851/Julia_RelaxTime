---
title: Models 原生实现全仓替换 legacy 并移除开发计划
archived: true
original: docs/dev/active/2026-02-27_models原生实现全仓替换legacy并移除开发计划.md
archived_date: 2026-03-01
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Models 原生实现全仓替换 legacy 并移除开发计划

更新日期：2026-02-27

> 目标：全仓统一使用 `src/models` 原生实现；若 `legacy` 存在 `models` 缺失能力，先补写到 `models`，最终移除 `src/legacy_pnjl/**`。

---

## 1. 背景与目标

- [x] 目标 A：所有对外入口（扫描/导数/workflow）不再转发 `legacy_pnjl`。
- [x] 目标 B：`tests/**` 与 `scripts/**` 不再调用 `Models.legacy_pnjl_module()` 或 `Main.PNJL`。
- [x] 目标 C：`legacy` 能力缺口全部补齐到 `models`，并通过回归。
- [x] 目标 D：删除 `src/legacy_pnjl/**` 后 smoke 与门禁稳定通过。

---

## 2. 范围与非范围

### 2.1 范围（In Scope）

- [x] 按顺序替换：入口层 → 调用方 → 缺口补齐 → 物理删除。
- [x] 将 `legacy` 中 `models` 缺失实现迁移/重写到 `src/models/**`（当前兼容实现已下沉至 `src/models/pnjl/compat/**`）。
- [x] 更新相关测试、脚本、文档与门禁脚本的路径口径。

### 2.2 非范围（Out of Scope）

- [ ] 不新增业务功能（仅做架构替换与等价迁移）。
- [ ] 不做与迁移无关的大规模数值重构。

---

## 3. 执行顺序（强约束）

### Batch N1：入口层去 legacy（先做）

- [x] `src/models/entrypoints.jl` 改为直接调用 `src/models/scans/*` 与 `src/models/workflows/*`，移除 `legacy_pnjl` 转发。
- [x] `src/models/derivative_entrypoints.jl` 改为 `models` 原生导数链路，移除 `PNJL.ThermoDerivatives` 依赖。
- [x] `src/models/workflows/*` 去除 `legacy_pnjl` 路径 include（如 `TransportWorkflow` / `MesonMassWorkflow`）。

验收：
- `UNIT_PROFILE=smoke julia --project=. tests/unit/runtests.jl`
- `julia --project=. scripts/dev/check_pnjl_migration_guard.jl --base HEAD --head HEAD`

### Batch N2：调用方切换（第二步）

- [x] `tests/**` 全量替换 `Models.legacy_pnjl_module()`、`Main.PNJL` 调用为 `Models` 原生 API。
- [x] `scripts/**` 全量替换 legacy 入口为 `Models` 原生 API。
- [x] 修复迁移过程中出现的 include 顺序/world-age/命名冲突。

验收：
- `UNIT_PROFILE=smoke julia --project=. tests/unit/runtests.jl`
- `julia --project=. scripts/dev/check_data_output_path_guard.jl --base HEAD --head HEAD`

### Batch N3：能力缺口补齐（第三步，发现即补）

- [x] 建立缺口清单：`legacy` 有但 `models` 无（功能、参数、返回结构、数值口径）。
- [x] 缺口实现策略：优先在 `src/models/**` 原生实现；禁止新增对 `legacy` 的反向依赖。
- [x] 高风险能力优先：
  - [x] 导数与隐函数求解链路
  - [x] 扫描配置与约束模式（FixedMu/FixedRho/FixedAsymmetricRho/...）
  - [x] workflow（transport/meson）
  - [x] 相变分析（S-shape/Maxwell/crossover）
  - [x] 磁场相关能力（若仍依赖 `legacy`）

验收：
- `UNIT_PROFILE=smoke julia --project=. tests/unit/runtests.jl`
- 对每个补齐点增加最小等价回归（单点或小网格对比）

### Batch N4：移除 legacy（最后一步）

- [x] 删除 `src/legacy_pnjl/**`。
- [x] 清除 `src/**`, `tests/**`, `scripts/**` 中所有 `legacy_pnjl`/`Main.PNJL` 直接依赖。
- [x] 提交级验证：新 HEAD 下 legacy 路径不存在。

验收：
- `UNIT_PROFILE=smoke julia --project=. tests/unit/runtests.jl`
- `julia --project=. scripts/dev/run_prune_wave_gate.jl --base HEAD --head HEAD`
- `git ls-tree -r --name-only HEAD src/legacy_pnjl`

---

## 4. 缺口治理规则（必须遵守）

- [x] 发现缺口时，不回退到 legacy；直接在 `models` 实现。
- [x] 保持 API 兼容优先：先兼容调用签名，再优化内部实现。
- [x] 每个缺口附带验证证据（测试/脚本输出/基线比对）。
- [x] 对可能引起口径漂移的改动，执行固定点回归与门禁校验。

---

## 5. 测试与验收标准

### 5.1 必过

- [x] `UNIT_PROFILE=smoke` 全绿。
- [x] migration/data-output/prune-wave 三门禁通过。

### 5.2 必稳

- [x] 无 `legacy_pnjl` 运行时依赖。
- [x] 默认扫描 include/runtime dependency 审计为 0。
- [x] `models-invokelatest-audit` 与 allowlist 对齐（若有漂移需同步治理）。

### 5.3 可追溯

- [x] 每批次有执行台账追加记录。
- [x] 关键产物归档到 `data/outputs/results/*`。

---

## 6. 里程碑

- [x] M1：完成 N1（入口层去 legacy）。
- [x] M2：完成 N2（调用方切换）。
- [x] M3：完成 N3（缺口补齐并回归通过）。
- [x] M4：完成 N4（删除 legacy + 提交级验证）。

---

## 7. DoD（完成定义）

- [x] 仓库中不再存在 `src/legacy_pnjl/**`。
- [x] 主路径（入口/测试/脚本）全部走 `models` 原生实现。
- [x] smoke + 门禁全通过，且有执行台账与产物证据。

---

## 8. 风险与回退

- [ ] 风险 R1：迁移后数值漂移。
  - 回退：按批次回滚；保留固定点基线与对比产物。
- [ ] 风险 R2：高耦合模块迁移导致连锁失败。
  - 回退：按 N1→N4 批次拆分提交，逐批回退。
- [ ] 风险 R3：遗漏隐式依赖（测试未覆盖）。
  - 回退：补充 grep 门禁与 targeted smoke。

---

## 9. 执行记录（分离维护）

- [x] 本文档仅维护计划与状态。
- [x] 执行细节追加到：
  - `docs/dev/active/2026-02-27_models原生实现全仓替换legacy并移除执行台账.md`
