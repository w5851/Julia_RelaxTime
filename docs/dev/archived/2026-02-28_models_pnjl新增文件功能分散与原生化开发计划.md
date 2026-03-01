---
title: models/pnjl 新增文件功能分散与原生化开发计划
archived: true
original: docs/dev/active/2026-02-28_models_pnjl新增文件功能分散与原生化开发计划.md
archived_date: 2026-03-01
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# models/pnjl 新增文件功能分散与原生化开发计划

更新日期：2026-02-28

> 目标：以 Git 新增文件为边界，逐个审查 `src/models/pnjl` 新增文件的职责归属；将可通用能力下沉到 `src/models` 主域原生模块，减少/移除桥接转发（如 `ThermoDerivativesBridge`）并保持数值与门禁稳定。

---

## 1. 背景与目标

- [x] 目标 A：基于 Git 准确识别 `src/models/pnjl` 新增文件清单（新增而非修改）。
- [x] 目标 B：对新增文件逐个形成“保留 pnjl / 分散到主域 / 合并删除”的归属判定。
- [x] 目标 C：优先将可通用能力原生化到 `src/models`（示例：导数能力收口到 `src/models/derivatives/ThermoDerivatives.jl`）。
- [x] 目标 D：桥接层最小化，仅保留迁移期必要壳层并给出退役路径。
- [x] 目标 E：全流程通过 smoke + migration/data-output/prune-wave 门禁。

---

## 2. 范围与非范围

### 2.1 范围（In Scope）

- [x] `src/models/pnjl/**` 中“Git 新增”文件的逐文件审查。
- [x] 关联主域模块：`src/models/derivatives/**`, `src/models/core/**`, `src/models/scans/**`, `src/models/entrypoints.jl`。
- [x] 桥接层审查：`src/models/derivatives/ThermoDerivativesBridge.jl` 等转发层。
- [x] 迁移落地后的引用收口、门禁验证与台账记录。

### 2.2 非范围（Out of Scope）

- [ ] 不新增物理模型功能，不改变对外数值语义。
- [ ] 不做与本次分散目标无关的大规模目录重构。
- [ ] 不在同一批次内引入多项实验性优化。

---

## 3. 判定标准（逐文件必填）

对每个新增文件输出以下结论字段：

- [x] 文件定位：入口 / 领域实现 / 通用算法 / 桥接壳层 / 配置适配。
- [x] 归属判定：
  - 保留 `src/models/pnjl/**`（模型特异）
  - 分散到 `src/models/**`（跨模型可复用）
  - 合并删除（重复或纯转发）
- [x] 依赖评估：是否引入 PNJL 专属常量、类型、边界条件。
- [x] 可迁移性等级：高 / 中 / 低（给出理由）。
- [x] 变更方案：目标路径、调用方改造点、兼容策略。

判定原则：

- [x] P1：凡跨模型复用能力，优先主域原生实现，不长期停留在桥接层。
- [x] P2：凡 PNJL 物理语义强耦合实现，保留在 `src/models/pnjl/**`。
- [x] P3：桥接文件默认视为迁移期资产，需标注退役条件与时间点。

---

## 4. 执行批次（可勾选）

### Batch E1：Git 新增清单冻结

- [x] 使用 Git 命令生成 `src/models/pnjl/**` 新增文件清单（相对基线分支）。
- [x] 过滤重命名/移动，确保“新增”语义准确。
- [x] 产出审查输入表（文件路径、类型、调用入口、依赖关系）。

验收：
- [x] 得到可复核的“新增文件唯一清单”，并可追溯到 Git 结果。

### Batch E2：逐文件职责审查

- [x] 对 E1 全量文件逐个给出归属判定与迁移优先级。
- [x] 标注“可立即分散”“需分步迁移”“暂不迁移”三类。
- [x] 识别重复能力与桥接层可删除项。

验收：
- [x] 形成“逐文件审查矩阵（文件 -> 判定 -> 理由 -> 下一步）”。

### Batch E3：导数链路原生化优先改造（首批）

- [x] 评估并改造 `src/models/pnjl/derivatives/ThermoDerivatives.jl` 能力下沉路径。
- [x] 将通用导数能力收口到 `src/models/derivatives/ThermoDerivatives.jl`。
- [x] 缩减 `src/models/derivatives/ThermoDerivativesBridge.jl` 到最小必要壳层或移除。
- [x] 修正调用方到主域原生入口。

验收：
- [x] 导数主链路不再依赖桥接转发完成核心计算。

### Batch E4：其余可分散能力迁移与收口

- [x] 按优先级迁移 core/analysis/solver 中可复用能力到主域。
- [x] 保留 PNJL 特异模块在 `src/models/pnjl/**`，并补边界注释/文档。
- [x] 移除已完成替代的桥接/重复文件与引用。

验收：
- [x] 迁移后目录职责清晰，桥接数量明显减少且可解释。

### Batch E5：终验与留痕

- [x] 运行 smoke：`UNIT_PROFILE=smoke julia --project=. tests/unit/runtests.jl`
- [x] 运行 migration guard：`julia --project=. scripts/dev/check_pnjl_migration_guard.jl --base HEAD --head HEAD`
- [x] 运行 data-output guard：`julia --project=. scripts/dev/check_data_output_path_guard.jl --base HEAD --head HEAD`
- [x] 运行 prune-wave gate：`julia --project=. scripts/dev/run_prune_wave_gate.jl --base HEAD --head HEAD`
- [x] append-only 追加执行台账（命令、产物、结果、主线映射）。

验收：
- [x] 四重门禁全绿，且迁移证据与台账完备。

---

## 5. 交付物

- [x] 交付 1：Git 新增文件冻结清单。
- [x] 交付 2：逐文件职责审查矩阵（含可迁移性等级）。
- [x] 交付 3：首批原生化改造 PR 级变更（导数链路优先）。
- [x] 交付 4：门禁结果与产物路径汇总。
- [x] 交付 5：执行台账追加记录。

---

## 6. 风险与回退

- [ ] 风险 R1：误将 PNJL 特异逻辑下沉为“通用逻辑”导致语义污染。
  - 回退：按文件粒度回滚，恢复至 pnjl 子域并补判定依据。
- [ ] 风险 R2：桥接移除过快导致调用链断裂。
  - 回退：短期保留兼容壳层，逐调用方迁移后再删除。
- [ ] 风险 R3：导数链路原生化引入数值漂移。
  - 回退：基于现有 baseline 与 smoke 回归，分批启用。

---

## 7. 里程碑

- [x] M-E1：完成 Git 新增清单冻结。
- [x] M-E2：完成逐文件审查矩阵。
- [x] M-E3：完成导数链路首批原生化。
- [x] M-E4：完成其余可分散能力迁移收口。
- [x] M-E5：完成门禁与台账闭环。

---

## 8. DoD（完成定义）

- [x] `src/models/pnjl` 新增文件均有明确归属判定与处理结论。
- [x] 至少一条主链路（导数）完成“桥接转发 -> 主域原生”收口。
- [x] 不引入新的 `src/models/** -> 临时桥接` 长期依赖。
- [x] 四重门禁通过，执行台账可审计。

---

## 9. E1 新增文件冻结清单（Git）

`src/models/pnjl` 新增文件（17）：

- `src/models/pnjl/DualBranchScanEntry.jl`
- `src/models/pnjl/PNJL.jl`
- `src/models/pnjl/analysis/PhaseTransition.jl`
- `src/models/pnjl/core/Core.jl`
- `src/models/pnjl/core/EquilibriumFacade.jl`
- `src/models/pnjl/core/Integrals.jl`
- `src/models/pnjl/core/MagneticIntegrals.jl`
- `src/models/pnjl/core/MagneticThermodynamics.jl`
- `src/models/pnjl/core/ModelThermodynamics.jl`
- `src/models/pnjl/core/ThermoFacade.jl`
- `src/models/pnjl/core/Thermodynamics.jl`
- `src/models/pnjl/derivatives/ThermoDerivatives.jl`
- `src/models/pnjl/solver/Conditions.jl`
- `src/models/pnjl/solver/ConstraintModes.jl`
- `src/models/pnjl/solver/ImplicitSolver.jl`
- `src/models/pnjl/solver/SeedStrategies.jl`
- `src/models/pnjl/solver/Solver.jl`

---

## 10. E2 逐文件审查矩阵（结论版）

### 10.1 可立即分散（已执行）

- `src/models/pnjl/derivatives/ThermoDerivatives.jl`
  - 判定：分散到 `src/models/**`（高）
  - 理由：对外调用面是 `Models.mass_derivatives/thermo_derivatives`，不应继续经 bridge 跳转。
  - 落地：新增 `src/models/derivatives/ThermoDerivatives.jl`，`src/models/derivative_entrypoints.jl` 改为加载该模块；删除 `src/models/derivatives/ThermoDerivativesBridge.jl`。

### 10.2 需分步迁移（本批完成判定，暂不搬迁）

- `src/models/pnjl/core/{Integrals,Thermodynamics,MagneticIntegrals,MagneticThermodynamics,ModelThermodynamics}.jl`
  - 判定：需分步迁移（中）
  - 理由：当前与 PNJL 常量、求解状态和磁场语义强耦合；直接平移会带来接口震荡。
  - 下一步：先抽公共数值核，再做 `src/models/core/**` 原生化。

- `src/models/pnjl/solver/{ConstraintModes,Conditions,ImplicitSolver,SeedStrategies,Solver}.jl`
  - 判定：需分步迁移（中）
  - 理由：约束模式和 seed 策略直接绑定 PNJL 五维状态与物理约束。
  - 下一步：先抽象通用 solver 接口，再下沉共性模块。

### 10.3 暂不迁移（保留 pnjl 子域）

- `src/models/pnjl/PNJL.jl`
- `src/models/pnjl/DualBranchScanEntry.jl`
- `src/models/pnjl/analysis/PhaseTransition.jl`
- `src/models/pnjl/core/{Core,ThermoFacade,EquilibriumFacade}.jl`

判定：保留 `src/models/pnjl/**`（低）
理由：作为 PNJL 子域入口/聚合/相图分析承载层，短期内跨模型复用价值不足。

---

## 11. E3/E4 实际改动摘要

- 新增：`src/models/derivatives/ThermoDerivatives.jl`
- 修改：`src/models/derivative_entrypoints.jl`（Bridge -> ThermoDerivatives）
- 删除：`src/models/derivatives/ThermoDerivativesBridge.jl`
- 修改：`config/ci/models_invokelatest_allowlist.toml`（`bridge.derivatives_module_ref` -> `td.derivatives_module_ref`）

补充（F2/F3/F4 第一批最小可迁移子集）：

- 新增：`src/models/scans/ScanEntrypoints.jl`（主域 scan 入口壳层）
- 新增：`src/models/solver/ConstraintModes.jl`（主域 solver 接口壳层）
- 新增：`src/models/solver/Solver.jl`（主域 solver 模块壳层）
- 修改：`src/models/entrypoints.jl`（`ScanEntrypointsBridge` -> `ScanEntrypoints`）
- 修改：`src/models/workflows/TransportWorkflow.jl`（scan 入口切换到 `ModelsScanEntrypoints`）
- 修改：`src/models/workflows/MesonMassWorkflow.jl`（scan 入口切换到 `ModelsScanEntrypoints`）
- 删除：`src/models/scans/ScanEntrypointsBridge.jl`
- 修改：`config/ci/models_invokelatest_allowlist.toml`
  - `entrypoints`：`bridge.pnjl_module_ref` -> `scan.pnjl_module_ref`
  - `workflow`：`ModelsScanEntrypointsBridge.pnjl_module_ref` -> `ModelsScanEntrypoints.pnjl_module_ref`

---

## 12. E5 门禁结果与产物

- smoke：`845/845` 通过
- migration guard：通过（`observed=15, allowlisted=15`）
- data-output guard：通过
- prune-wave gate：通过
  - `data/outputs/results/pnjl_prune_wave_snapshot_20260228_190255.txt`
  - `data/outputs/results/models_invokelatest_allowlist_delta_20260228_190255.txt`
  - `data/outputs/results/pnjl_scan_default_include_audit_20260228_190255.txt`

补充（F2/F3/F4 后复验）：

- smoke：`845/845` 通过
- migration guard：通过（`observed=15, allowlisted=15`）
- data-output guard：通过
- prune-wave gate：通过
  - `data/outputs/results/pnjl_prune_wave_snapshot_20260228_203222.txt`
  - `data/outputs/results/models_invokelatest_allowlist_delta_20260228_203222.txt`
  - `data/outputs/results/pnjl_scan_default_include_audit_20260228_203222.txt`

---

## 13. F1-F5 执行结论（本轮）

- F1（接口矩阵）：完成，识别 `core` 主域空壳已具备，`solver` 主域空壳缺失。
- F2（主域空壳接口）：完成，补齐 `src/models/solver/*` 最小壳层。
- F3（调用方收口）：完成，`entrypoints + workflows` 的 scan 入口均已从旧 bridge 收口到主域 `ScanEntrypoints`。
- F4（首批内部替换）：完成低风险替换，未触及高耦合物理算法实现。
- F5（门禁与留痕）：完成，smoke 与三门禁均通过，产物与台账可追溯。

---

## 14. G1 执行结论（solver 真下沉：类型 + 通用接口）

- G1.1（类型下沉）：完成，将 `src/models/solver/ConstraintModes.jl` 从 `module_ref` 壳层替换为主域真实约束类型定义（`ConstraintMode`、`Fixed*`、`state_dim/param_dim/constraint_description`）。
- G1.2（通用接口下沉）：完成，将 `src/models/solver/Solver.jl` 从 `module_ref` 壳层替换为主域 `solve_constraint(model, mode, T_fm; kwargs...)` 多重派发接口，并映射到 `solve_fixed*` 内核。
- G1.3（主域接入）：完成，`src/models/Models.jl` 新增 include 接入 `solver/ConstraintModes.jl` 与 `solver/Solver.jl`，并导出 `solve_constraint`。
- G1.4（兼容性修正）：完成，移除 `Models.jl` 对 `Fixed*` 的直接导出，避免 `using .Models` 与 Main 现有标识符冲突告警；保持 `Models.FixedMu` 等命名空间访问。
- G1.5（验证）：通过。
  - smoke：`845/845` 通过
  - migration guard（range HEAD...HEAD）：通过
  - data-output guard（range HEAD...HEAD）：通过
  - prune-wave gate（range HEAD...HEAD）：通过
  - 产物：
    - `data/outputs/results/pnjl_prune_wave_snapshot_20260228_224507.txt`
    - `data/outputs/results/models_invokelatest_allowlist_delta_20260228_224507.txt`
    - `data/outputs/results/pnjl_scan_default_include_audit_20260228_224507.txt`

---

## 15. G2 执行结论（TrhoScan 调用链闭环到主域 solve_constraint）

- G2.1（调用路径切换）：完成，在 `src/models/scans/TrhoScan.jl` 中新增 models 路径适配：当 `solver_backend=:models` 时，不再走 `PNJL.ImplicitSolver.solve`，改为调用 `Main.Models.solve_constraint`。
- G2.2（兼容层设计）：完成，保持 seed 策略与 scan glue 仍使用 PNJL 侧约束类型分派（避免破坏 `get_seed/extend_seed`），仅在实际求解前映射到主域 `Models.Fixed*`，并将返回值适配回 `SolverResult`，因此 `finalize/refine/promote` 链路无侵入。
- G2.3（覆盖范围）：完成，覆盖 `TrhoScan` 的两条入口：`_attempt_with_strategy` 与 `_solve_point`。
- G2.4（验证）：通过。
  - `tests/unit/pnjl/test_trho_scan_solver_backend_models_smoke.jl`：通过（含 fixed_rho / fixed_asymmetric_rho / solver_backend=:auto+models）
  - smoke：`845/845` 通过
  - migration guard（range HEAD...HEAD）：通过
  - data-output guard（range HEAD...HEAD）：通过

---

## 16. G3 执行结论（TmuScan + DualBranchScan 调用链收口）

- G3.1（TmuScan 调用路径切换）：完成，在 `src/models/scans/TmuScan.jl` 中新增 models 路径适配：当 `solver_backend=:models` 时，固定化学势点不再走 `PNJL.ImplicitSolver.solve`，改为调用 `Main.Models.solve_constraint`。
- G3.2（TmuScan 兼容层设计）：完成，保留 scan 侧 seed 组合与 `SolverResult` 后处理链，新增 `FixedMu -> Main.Models.FixedMu` 映射与返回值适配，保证 `finalize/refine/promote` 链路无侵入。
- G3.3（DualBranchScan 收口）：完成，为 `run_dual_branch_scan` 增加 `solver_backend`（支持 `:legacy|:models|:auto`）；在 `_solve_point` 中引入同样的 models 求解适配（默认行为保持 legacy 不变）。
- G3.4（兼容性修正）：完成，针对 fixed-μ models 内核不接受 `iterations` 等 legacy `nlsolve` 参数的问题，增加关键词白名单透传（仅保留 `solver/residual_norm_max/physicality_check`）。
- G3.5（验证）：通过。
  - `tests/unit/pnjl/test_tmu_scan_solver_backend_models_smoke.jl`：通过（含 `solver_backend=:models` 与 `:auto+models`）
  - smoke：`845/845` 通过
  - migration guard（range HEAD...HEAD）：通过
  - data-output guard（range HEAD...HEAD）：通过
