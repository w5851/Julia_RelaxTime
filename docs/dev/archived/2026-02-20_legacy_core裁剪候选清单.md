---
title: legacy core 裁剪候选清单（仅候选，不立即删除）
archived: true
original: docs/dev/active/2026-02-20_legacy_core裁剪候选清单.md
archived_date: 2026-02-21
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# legacy core 裁剪候选清单（仅候选，不立即删除）

更新日期：2026-02-20

---

## 1. 目的

- 给出可追踪的“候选裁剪项”，避免一次性大删导致主链回归失控。
- 每个候选项必须满足：调用证据 + 回退方案 + 最小验证集。

---

## 2. 裁剪门槛（必须同时满足）

- `models` 路径具备等价能力（接口与数值行为）。
- smoke 与定向回归通过。
- 脚本/工作流无隐式依赖，或已迁移到明确新入口。
- 出现问题时可通过 `solver_backend=:legacy` 或 `thermo_backend=:legacy` 回退。

---

## 3. 候选项清单

## C1（低风险）：`src/pnjl/core/ModelThermodynamics.jl` 内重复桥接逻辑

状态：已推进（转发收敛已完成，未删文件）

原因：
- 与 `src/pnjl/core/ThermoFacade.jl` 在 models backend 分发上存在职责重叠。

前置条件：
- `ThermoFacade` 内完成统一 dispatch，不再需要重复桥接入口。

本轮推进（2026-02-21）：
- 已新增 `ThermoFacade.ensure_models_loaded()` / `ThermoFacade.rho0()` 统一入口；
- `EquilibriumFacade`、`ThermoDerivatives` 已不再直接引用 `ThermoFacade.ModelThermodynamics.*` 内部细节。

回退方案：
- 保留现文件并仅做内部转发；先不删除实体文件。

最小验证：
- `models/test_pnjl_thermo_bridge_multipoint_smoke.jl`
- `pnjl/test_tmu_scan_solver_backend_models_smoke.jl`

## C2（中风险）：`src/models/pnjl/PNJLModel.jl` 中 legacy solver 直连默认路径

状态：已推进（新增 `:auto` 能力，默认仍 `:legacy`）

原因：
- 默认 `solver_backend=:legacy` 仍是依赖点；当前已支持 `fallback_legacy_on_failure`，可作为迁移过渡。

本轮推进（2026-02-21）：
- `PNJLModel.solve_gap` 已支持 `solver_backend=:auto`，内部优先走 models 求解；
- 失败时按既有规则受控回退到 legacy（仅对称 `mu` 可回退）；
- 默认值保持 `:legacy`，避免主链行为突变。

前置条件：
- models solver 在关键点位收敛稳定。
- bridge/workflow baseline 不因默认切换而偏移。

回退方案：
- 保留 `solver_backend=:legacy`。
- 默认行为切换失败时回滚到 legacy 默认。

最小验证：
- `models/test_pnjl_solve_gap_backend_switch_smoke.jl`
- `relaxtime/test_transport_legacy_models_bridge_smoke.jl`
- `UNIT_PROFILE=smoke`

## C3（中风险）：`src/pnjl/core/EquilibriumFacade.jl` 默认后端策略

状态：已推进（scan 链路已接入 `:auto`）

原因：
- 当前默认仍为 `:legacy`，但已支持 `:auto`；后续可按场景逐步切换调用方。

前置条件：
- `:auto` 在 workflow/scan 场景验证充分。

本轮推进（2026-02-21）：
- `TmuScan`/`TrhoScan` 已支持 `solver_backend=:auto`；
- 解析规则：`:auto` 在 `thermo_backend=:models` 时转为 `:models`，否则转为 `:legacy`；
- 默认值保持 `:legacy`，避免行为突变。

回退方案：
- 保持默认 `:legacy` 不变，只允许显式启用 `:auto`。

最小验证：
- `pnjl/test_tmu_scan_solver_backend_models_smoke.jl`
- `pnjl/test_trho_scan_solver_backend_models_smoke.jl`
- `UNIT_PROFILE=smoke`

新增验证：
- `pnjl/test_tmu_scan_solver_backend_models_smoke.jl`（新增 `solver_backend=:auto` 子用例）
- `pnjl/test_trho_scan_solver_backend_models_smoke.jl`（新增 `solver_backend=:auto` 子用例）

---

## 4. 本轮结论

- 本轮不删除 legacy 实体文件。
- 先完成“能力并存 + 默认兼容 + 显式开关验证”，下一轮再进入实际删除任务单。