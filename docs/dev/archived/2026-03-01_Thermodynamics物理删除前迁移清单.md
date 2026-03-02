---
title: Thermodynamics 物理删除前迁移清单
archived: true
original: docs/dev/active/2026-03-01_Thermodynamics物理删除前迁移清单.md
archived_date: 2026-03-02
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Thermodynamics 物理删除前迁移清单

更新日期：2026-03-01

> 目标：在保持外部 API 与回归结果稳定的前提下，完成 src/models/pnjl/core/Thermodynamics.jl 的物理删除准备。

---

## 1. 当前状态

- 已完成：Thermodynamics.jl 已收缩为纯转发壳（不含主流程物理计算细节）。
- 已完成：PNJL.jl、solver、derivatives、ThermoFacade 的节点缓存入口已切换到 PNJLCore 命名空间。
- 未完成：legacy backend 仍是主线兼容路径之一，多个模块与测试仍显式依赖 Thermodynamics 模块符号或 :legacy 语义。

结论：当前可进入“删除前迁移”阶段，但不应直接删除 Thermodynamics.jl 文件。

---

## 2. 阻塞项清单（按类型）

### 2.1 直接模块依赖（include/using Thermodynamics）

- src/models/pnjl/core/ThermoFacade.jl
  - 显式 include legacy Thermodynamics 路径（当前唯一保留的兼容壳入口），并在 `:legacy` 分支转发到 Thermodynamics.*

### 2.2 代码路径阻塞（legacy 语义仍在）

- src/models/legacy/LegacyPNJLModel.jl
  - 作为 legacy 哨兵模型，内部显式固定走 `thermo_backend=:legacy`（mass/omega_components/number_densities）
  - `solve_gap` 默认参数仍为 `thermo_backend=:legacy`
- src/models/pnjl/core/ThermoFacade.jl
  - 兼容壳仍 include `Thermodynamics.jl`，并保留 `:legacy` 转发分支
- src/models/pnjl/derivatives/ThermoDerivatives.jl
  - 兼容辅助函数仍直接调用 `ThermoFacade.Thermodynamics.calculate_thermo/rho`
- src/models/pnjl/PNJLModel.jl
  - legacy 求解器分支与 fallback 路径仍显式传 `thermo_backend=:legacy`（兼容策略保留）

### 2.3 测试阻塞（显式覆盖 :legacy）

- tests/unit/models/test_pnjl_thermo_bridge_multipoint_smoke.jl
- tests/unit/pnjl/test_thermo_derivatives.jl
- tests/unit/pnjl/test_bulk_viscosity.jl
- tests/unit/pnjl/test_scan_fixedpoint_baseline_smoke.jl
- tests/unit/pnjl/test_scan_solver_boundary_error_smoke.jl
- tests/unit/relaxtime/test_transport_workflow.jl
- tests/unit/relaxtime/test_workflow_paramtypes_equivalence_smoke.jl
- tests/unit/relaxtime/test_workflow_paramtypes_mixedmode_smoke.jl

### 2.4 批次 C 前置检查结果（2026-03-01）

- [x] 迁移门禁：`scripts/dev/check_pnjl_migration_guard.jl --base HEAD --head HEAD`（exit=0）
- [x] 输出路径门禁：`scripts/dev/check_data_output_path_guard.jl --base HEAD --head HEAD`（exit=0）
- [x] 代码扫描确认：`src/models/**` 中对 `Thermodynamics` 的直接 include/using 已收敛到 `ThermoFacade` 兼容壳链路
- [ ] 测试侧仍有 12 处显式 `thermo_backend=:legacy` 对照用法（需作为批次 C 删除前的哨兵保留或迁移策略输入）

---

## 3. 建议分批实施（最小回滚单元）

### 批次 A：消除直接 include/using 依赖（不删文件）

目标：让业务模块不再直接 import Thermodynamics 模块符号。

动作：
- 将 PNJL.jl 与 core/Core.jl 导出面切到 ThermoFacade 与 Main.Models/PNJLCore 对应接口。
- 将 solver/Conditions.jl 与 solver/ImplicitSolver.jl 里 using .Thermodynamics:* 替换为 ThermoFacade 对应入口。
- 将 core/MagneticThermodynamics.jl 对 Thermodynamics 的函数依赖改为 PNJLCore 或 Main.Models 对应函数。

进展（2026-03-01）：
- [x] `solver/Conditions.jl` 已移除对 `Thermodynamics` 的直接 include/using，改为仅经 `ThermoFacade`；`ρ0` 改由 `Constants_PNJL` 提供。
- [x] `solver/ImplicitSolver.jl` 已移除对 `Thermodynamics` 的直接 include/using，改为仅经 `ThermoFacade`；`ρ0` 改由 `Constants_PNJL` 提供。
- [x] `PNJL.jl` 已移除对 `Thermodynamics` 的直接 include/using，改为通过 `ThermoFacade` 绑定 legacy API（对外导出不变）。
- [x] `core/Core.jl` 已移除对 `Thermodynamics` 的直接 include/using，改为通过 `ThermoFacade` 间接绑定 legacy API。
- [x] `core/MagneticThermodynamics.jl` 已移除对 `Thermodynamics` 的直接 include/using，改为 `PNJLCore + Integrals` 本地组装（并消除加载环）。

验收：
- 无模块再使用 include("Thermodynamics.jl")（仅允许兼容壳保留）。
- 现有 smoke 通过。

本批次已验证：
- `tests/unit/pnjl/test_solver_constraints_models_backend_smoke.jl`（43/43, exit=0）
- `tests/unit/pnjl/test_thermo_derivatives.jl`（全组通过, exit=0）
- `tests/unit/models/test_pnjl_thermo_bridge_multipoint_smoke.jl`（93/93, exit=0）

### 批次 B：legacy 后端退场策略

目标：将默认路径从 :legacy 切换到 :models，并把 :legacy 保留为可选兼容开关。

动作：
- 调整 solver/derivatives/workflows/scans 默认 thermo_backend 为 :models（分文件逐步切）。
- LegacyPNJLModel 保留但不作为默认工厂路径。
- ThermoFacade 保留 :legacy 分支，但内部尽量转发到 models 等价实现。

进展（2026-03-01，首批）：
- [x] `solver/Conditions.jl`：`GapParams` 默认 `thermo_backend` 已切换为 `:models`。
- [x] `solver/ImplicitSolver.jl`：`solve*` / implicit config / `create_implicit_solver` / `solve_with_derivatives` 默认 `thermo_backend` 已切换为 `:models`。
- [x] `derivatives/ThermoDerivatives.jl`：`CURRENT_THERMO_BACKEND` 与导数主接口默认 `thermo_backend` 已切换为 `:models`。
- [x] `scans/TmuScan.jl`、`scans/TrhoScan.jl`、`scans/DualBranchScan.jl` 默认 `thermo_backend` 已切换为 `:models`。
- [x] `workflows/MesonMassWorkflow.jl`、`workflows/TransportWorkflow.jl` 默认 `thermo_backend` 已切换为 `:models`。
- [x] `core/ThermoFacade.jl` 所有 backend 统一入口默认 `thermo_backend` 已切换为 `:models`（legacy 分支保留）。
- [x] `core/EquilibriumFacade.jl` 默认已切换为 `thermo_backend=:models`、`solver_backend=:auto`。
- [x] `LegacyPNJLModel` 保留为显式兼容模型（`:LegacyPNJL`），非默认工厂路径。

验收：
- models 默认路径下 unit smoke 全绿。
- 关键 legacy 对照测试保留少量哨兵用例。

本批次首批已验证：
- `tests/unit/pnjl/test_solver_constraints_models_backend_smoke.jl`（43/43, exit=0）
- `tests/unit/pnjl/test_thermo_derivatives.jl`（全组通过, exit=0）
- `tests/unit/models/test_pnjl_thermo_bridge_multipoint_smoke.jl`（93/93, exit=0）
- `tests/unit/models/test_pnjl_solve_gap_backend_switch_smoke.jl`（15/15, exit=0）

本批次第二波已验证：
- `tests/unit/pnjl/test_tmu_scan_solver_backend_models_smoke.jl`（9/9, exit=0）
- `tests/unit/pnjl/test_trho_scan_solver_backend_models_smoke.jl`（18/18, exit=0）
- `tests/unit/relaxtime/test_transport_workflow_solver_backend_switch_smoke.jl`（11/11, exit=0）

本批次收口已验证：
- `tests/unit/pnjl/test_tmu_scan_smoke.jl`（12/12, exit=0）
- `tests/unit/pnjl/test_trho_scan_smoke.jl`（18/18, exit=0）
- `tests/unit/relaxtime/test_transport_workflow_smoke.jl`（35/35, exit=0）
- `tests/unit/relaxtime/test_meson_mass_workflow_smoke.jl`（15/15, exit=0）
- `tests/unit/models/test_models_unified_entrypoints_smoke.jl`（7/7, exit=0）

### 批次 C：物理删除 Thermodynamics.jl

目标：删除文件并保留最小兼容层（如需）。

动作：
- 移除所有 include/using Thermodynamics 的最后残余。
- 删除 src/models/pnjl/core/Thermodynamics.jl。
- 若外部脚本仍需旧路径，可用极薄 shim 文件或 deprecate 指向 ThermoFacade。

验收：
- 全量门禁通过。
- 无 runtime world-age/redefinition 问题。

---

## 4. 每批次最小验证命令

- julia --project=. tests/unit/models/test_pnjl_thermo_bridge_multipoint_smoke.jl
- julia --project=. tests/unit/pnjl/test_thermo_derivatives.jl
- julia --project=. tests/unit/pnjl/test_solver_constraints_models_backend_smoke.jl
- julia --project=. tests/unit/models/test_pnjl_solve_gap_backend_switch_smoke.jl
- julia --project=. scripts/dev/check_pnjl_migration_guard.jl --base HEAD --head HEAD
- julia --project=. scripts/dev/check_data_output_path_guard.jl --base HEAD --head HEAD

---

## 5. 建议提交拆分

- 提交 1：去除直接 include/using Thermodynamics（无语义改动）
- 提交 2：默认 backend 切换与 legacy 保留策略
- 提交 3：删除 Thermodynamics 文件与兼容 shim（如需要）

每个提交都应可独立回滚，并附对应 smoke/guard 结果。

---

## 6. 批次 C 实施草案（可直接执行）

> 目标：在不破坏现有外部 API 的前提下，完成 `src/models/pnjl/core/Thermodynamics.jl` 物理删除。

### C-0 执行边界（本草案约束）

- 不改物理公式口径，仅做“调用重定向 + 兼容层替换”。
- `:legacy` 语义保留为显式兼容开关；删除的是文件依赖，不是功能开关。
- 任何一步失败都可回滚到上一提交，不跨步叠加修复。

### C-1 兼容 shim 设计（先建替身，再删旧文件）

建议新增一个极薄兼容模块（文件名可用 `src/models/pnjl/core/ThermodynamicsCompat.jl`）：

- 对外导出旧符号名：
  - `calculate_mass_vec`, `calculate_chiral`, `calculate_U`, `calculate_U_derivative_T`
  - `calculate_pressure`, `calculate_omega`, `calculate_omega_components`
  - `calculate_rho`, `calculate_thermo`, `calculate_number_densities`, `ρ0`
- 内部全部转发到 `ThermoFacade` / `PNJLCore` / `Main.Models`。
- 不再包含 `Integrals` 子模块（若必须保兼容，提供最小别名到 `PNJLCore.cached_nodes` 相关能力）。

验收：
- 旧调用方无需改符号名即可运行。
- shim 不包含主流程物理实现细节。

执行进展（2026-03-01）：
- [x] 已新增 `src/models/pnjl/core/ThermodynamicsCompat.jl`，提供旧 API 形状并转发至 `PNJLCore/Main.Models`。
- [x] `src/models/pnjl/core/ThermoFacade.jl` 的 legacy include 已切换到 `ThermodynamicsCompat.jl`。
- [x] 旧 `src/models/pnjl/core/Thermodynamics.jl` 暂保留（尚未删除，符合 C-1 范围）。

C-1 验证：
- `tests/unit/pnjl/test_solver_constraints_models_backend_smoke.jl`（43/43, exit=0）
- `tests/unit/pnjl/test_thermo_derivatives.jl`（全组通过, exit=0）
- `tests/unit/models/test_pnjl_thermo_bridge_multipoint_smoke.jl`（93/93, exit=0）
- `tests/unit/models/test_models_unified_entrypoints_smoke.jl`（7/7, exit=0）

### C-2 调用方切换（从 Thermodynamics 壳切到 shim/facade）

按优先级修改：

1. `src/models/pnjl/core/ThermoFacade.jl`
  - 将 `_LEGACY_THERMO_PATH` 指向 shim（不再指向 `Thermodynamics.jl`）。
2. `src/models/pnjl/PNJL.jl` 与 `src/models/pnjl/core/Core.jl`
  - `const Thermodynamics = ...` 绑定改为 shim 模块。
3. `src/models/pnjl/derivatives/ThermoDerivatives.jl`
  - `ThermoFacade.Thermodynamics.*` 的兼容辅助调用改为 shim/facade 对应入口。
4. 其余 residual 引用
  - 全局扫描 `Thermodynamics.`，确保不再硬依赖被删除文件路径。

验收：
- `src/models/**` 中不再存在 `include("Thermodynamics.jl")`。
- legacy 对照测试仍可显式调用 `thermo_backend=:legacy`。

执行进展（2026-03-01）：
- [x] `src/models/**` 已无 `include("Thermodynamics.jl")` 或旧文件路径级硬引用。
- [x] `PNJL.jl`、`core/Core.jl` 仍通过 `ThermoFacade.Thermodynamics` 暴露兼容别名（来源已切换为 shim）。
- [x] `derivatives/ThermoDerivatives.jl` 的兼容调用保持通过 `ThermoFacade.Thermodynamics`，不依赖旧文件实体。

### C-3 物理删除与清理

在 C-1/C-2 通过后执行：

- 删除 `src/models/pnjl/core/Thermodynamics.jl`。
- 清理文档中的“仍依赖 Thermodynamics.jl”描述，改为“兼容 shim + facade”。

验收：
- 删除后无 world-age/redefinition 错误。
- 所有批次 B 门禁 + 关键 legacy 哨兵测试通过。

执行进展（2026-03-01）：
- [x] 已删除 `src/models/pnjl/core/Thermodynamics.jl`。
- [x] 关键回归通过（`UNIT_FILES=pnjl/test_solver_constraints_models_backend_smoke.jl,pnjl/test_thermo_derivatives.jl,models/test_pnjl_thermo_bridge_multipoint_smoke.jl,models/test_models_unified_entrypoints_smoke.jl`）：`195/195`，`exit=0`。

### C-4 门禁命令（建议最小集合）

- `julia --project=. tests/unit/pnjl/test_solver_constraints_models_backend_smoke.jl`
- `julia --project=. tests/unit/pnjl/test_thermo_derivatives.jl`
- `julia --project=. tests/unit/models/test_pnjl_thermo_bridge_multipoint_smoke.jl`
- `julia --project=. tests/unit/models/test_pnjl_solve_gap_backend_switch_smoke.jl`
- `julia --project=. tests/unit/pnjl/test_tmu_scan_smoke.jl`
- `julia --project=. tests/unit/pnjl/test_trho_scan_smoke.jl`
- `julia --project=. tests/unit/relaxtime/test_transport_workflow_smoke.jl`
- `julia --project=. tests/unit/relaxtime/test_meson_mass_workflow_smoke.jl`
- `julia --project=. scripts/dev/check_pnjl_migration_guard.jl --base HEAD --head HEAD`
- `julia --project=. scripts/dev/check_data_output_path_guard.jl --base HEAD --head HEAD`

执行进展（2026-03-01）：
- [x] 关键单测门禁通过（`UNIT_FILES=pnjl/test_solver_constraints_models_backend_smoke.jl,pnjl/test_thermo_derivatives.jl,models/test_pnjl_thermo_bridge_multipoint_smoke.jl,models/test_pnjl_solve_gap_backend_switch_smoke.jl,pnjl/test_tmu_scan_smoke.jl,pnjl/test_trho_scan_smoke.jl,relaxtime/test_transport_workflow_smoke.jl,relaxtime/test_meson_mass_workflow_smoke.jl`）：`283/283`，`exit=0`。
- [x] `scripts/dev/check_pnjl_migration_guard.jl --base HEAD --head HEAD`：`OK`。
- [x] `scripts/dev/check_data_output_path_guard.jl --base HEAD --head HEAD`：`OK`。
- [x] `ThermoFacade` 已收敛为单执行路径：`thermo_backend=:legacy` 仅作为兼容别名，统一映射到 `:models` 的 `:PNJL` 实现（不再存在 legacy 专用计算分支）。

### C-5 提交建议（收口版）

- 提交 A：新增 `ThermodynamicsCompat` + `ThermoFacade` 改指向（不删旧文件）
- 提交 B：调用方切换（PNJL/Core/Derivatives 等）
- 提交 C：删除 `Thermodynamics.jl` + 文档更新

### C-6 回滚点

- 回滚点 1：提交 A 前失败，直接回退 shim 引入。
- 回滚点 2：提交 B 后失败，保留 shim，回退调用方切换。
- 回滚点 3：提交 C 后失败，恢复 `Thermodynamics.jl` 文件并暂时保留双路径。

### C-7 环境备注（执行期）

- `nightly-full-regression.yml` 在 VS Code 中出现 “Unable to resolve action ...” 时，优先检查本机网络/代理解析；该文件与仓库其他 workflow 的 action 写法一致。
- 若 `git ls-remote https://github.com/actions/checkout.git refs/tags/v4` 可达，则通常属于编辑器扩展缓存或本地诊断状态；先 `Reload Window` 再复查。
