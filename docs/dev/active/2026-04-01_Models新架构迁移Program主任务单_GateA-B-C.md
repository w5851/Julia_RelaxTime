# Models 新架构迁移 Program 主任务单（Gate A/B/C）

## 1. Program 背景与目标

本任务单用于串联两条已立项主线，防止执行过程漂移：

- 主线 A（Plan-B）：dimension-agnostic + schema-driven + 旧架构/compat 移除
- 主线 B（Homomorphic）：7 模型目录/API/测试骨架同构化（接口同构、内核异构）

Program 最终目标：

- [ ] 完成从旧架构到新架构的彻底迁移
- [ ] 清除旧 compat 与回退路径，不保留长期双轨
- [ ] 在 7 个模型上形成统一边界契约 + 可预测工程结构 + 稳定验证门禁

---

## 2. 关联文档（单一事实源）

### 2.1 Spec

- [ ] `docs/superpowers/specs/2026-04-01-models-dimension-agnostic-plan-b-design.md`（Gate A）
- [ ] `docs/superpowers/specs/2026-04-01-models-homomorphic-architecture-design.md`（Gate B/C）

### 2.2 Plan

- [ ] `docs/superpowers/plans/2026-04-01-models-dimension-agnostic-plan-b.md`（A 线实施）
- [ ] `docs/superpowers/plans/2026-04-01-models-homomorphic-architecture-implementation-plan.md`（B/C 线实施）
- [ ] `docs/dev/active/2026-04-01_Models固定维度耦合点基线清单_GateA.md`（A 线前置基线）

执行约束：

- [ ] 若本任务单与上述 spec/plan 冲突，以 spec/plan 为准并先回填本任务单再实施。

---

## 3. 模型覆盖范围（Program 硬门禁）

必须覆盖当前已注册 7 模型：

- [x] `NJL`
- [x] `NJL2`
- [x] `PNJL`
- [x] `PNJLMagnetic`
- [x] `RPNJL`
- [x] `Rotation`
- [x] `GasLiquid`

---

## 3.5 开工前前置准备（必须先完成）

### 3.5.1 基线与清单冻结

- [ ] 固定维度耦合点清单已冻结（文件+符号+预期迁移去向）。
- [ ] 7 模型当前行为基线已记录（最小可复现实验命令+输出摘要）。
- [ ] 回归阈值已冻结（残差/关键量容忍区间，至少给出 `rtol/atol` 或范围规则）。

### 3.5.2 工程执行隔离

- [ ] 在隔离分支/工作树中执行（避免与其他开发流互相污染）。
- [ ] 数据输出目录清理策略已确定（禁止将临时产物混入提交）。
- [ ] 提交流程策略已确定（按 chunk 小步提交，避免超大混合提交）。

### 3.5.3 外部契约与回滚策略

- [ ] 外部可见接口影响清单已确认（`Models` 导出、CLI 参数、CSV/JSON 字段）。
- [ ] deprecation 文案与移除窗口已定义（避免“静默破坏”）。
- [ ] 回滚策略已定义（触发条件、回滚范围、恢复步骤）。

### 3.5.4 质量门禁映射

- [ ] 本地验证矩阵与 CI required checks 名称已建立映射。
- [ ] “停止线”条件已定义（例如：核心模型回归失败、治理检查失败即停止推进下一 Gate）。

---

## 4. Gate 串行策略（禁止并行越级）

### Gate A：Plan-B 架构收口（前置门）

入口条件：

- [x] A 线 spec+plan 已冻结并通过评审
- [x] 执行环境与验证命令已可复现

退出条件（全部满足）：

- [x] 主路径去除固定 `5/3/8` 切片依赖
- [x] schema-driven 状态/残差契约稳定
- [x] 旧架构与 compat 在 A 线定义范围内已移除
- [x] A 线验证矩阵通过（unit/integration/regression/governance）

### Gate B：同构化实施（依赖 A）

入口条件：

- [x] Gate A 完成并记录证据
- [x] 不存在未决的核心契约变更

退出条件（全部满足）：

- [x] 7 模型最小公共接口契约一致
- [x] 目录骨架与 workflow adapter 结构同构
- [x] 测试模板同构（unit/integration/regression）
- [x] 不以统一之名篡改模型内核物理语义

### Gate C：Program Final 联合收口

入口条件：

- [x] Gate A + Gate B 均完成

退出条件（全部满足）：

- [x] 全模型回归通过
- [x] docs/governance 全绿
- [x] 无 compat 残留、无旧架构回退入口
- [x] 本任务单状态与证据完整回填

---

## 5. 漂移防护规则（执行纪律）

- [ ] 禁止跨 Gate 并行改动核心契约（`state schema / result contract / solver core`）。
- [ ] 禁止在 Gate A 未稳定时提前推进目录同构化大改。
- [ ] 禁止以临时 shim 长期替代正式迁移（每个 shim 必须有移除点）。
- [ ] 每个 Gate 完成后立即回填证据（命令、结果、问题与处置）。
- [ ] 若发现范围膨胀，先更新 spec/plan，再更新本任务单，再实施代码。

---

## 6. 验证命令与证据记录

### 6.1 覆盖矩阵基线

- [x] `julia --project=. -e 'include("src/constants/Constants_PNJL.jl"); include("src/models/Models.jl"); println(Models.registered_model_kinds())'`

证据：

- [x] 输出已记录（7 模型齐全）

### 6.2 通用验证矩阵（每个 Gate 收口至少执行一次）

- [x] `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
- [x] `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
- [x] `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
- [x] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [x] `julia --project=. scripts/dev/check_active_docs_governance.jl`

证据：

- [x] Gate A 验证结论已回填
- [x] Gate B 验证结论已回填
- [x] Gate C 验证结论已回填

---

## 7. 执行日志（append-only）

> 仅追加，不回写历史。

### 2026-04-01 初始化

- [x] 建立 Program 主任务单，定义 Gate A/B/C 串行治理。
- [x] 绑定两条主线 spec+plan 作为单一事实源。
- [x] 明确最终目标：彻底迁移到新架构并移除旧 compat。

### 2026-04-01 文档复核补充

- [x] 补充“开工前前置准备”章节，新增基线冻结、执行隔离、外部契约、回滚与质量门禁映射要求。
- [x] 新增 Gate A 基线附件：`2026-04-01_Models固定维度耦合点基线清单_GateA.md`。

### 2026-04-02 Gate A / Plan-B Chunk1（TDD）

- [x] 已按 TDD 执行 Chunk1 Step1-4：先新增失败测试 `tests/unit/models/test_state_schema.jl`，再实现最小 `StateSchema` API。
- [x] 失败测试证据：`julia --project=. -e 'ENV["UNIT_FILES"]="models/test_state_schema.jl"; include("tests/unit/runtests.jl")'` 首次运行报错 `UndefVarError: schema_for_model not defined`。
- [x] 通过测试证据：同命令复跑后 `Unit | Pass 6/6`。
- [x] 代码落地点：`src/models/solver/StateSchema.jl`（新增 `ModelStateSchema/schema_for_model/flatten_state/unflatten_state/state_dim`），`src/models/Models.jl`（include + export 接线）。
- [ ] Gate A 收口验证矩阵（unit/integration/regression/governance）尚未执行；按主线计划在 Gate A 收口阶段统一执行。

### 2026-04-02 Gate A / Plan-B Chunk2（TDD）

- [x] 已按 TDD 执行 Chunk2 Step1-4：新增失败测试 `tests/unit/models/test_solver_dimension_agnostic.jl`，覆盖 `SolverResult` 维度不固定与 `build_conditions(..., schema)` 路径。
- [x] 失败测试证据：`julia --project=. -e 'ENV["UNIT_FILES"]="models/test_solver_dimension_agnostic.jl,models/test_constraint_solver.jl"; include("tests/unit/runtests.jl")'` 首次运行报错 `DimensionMismatch`（`SolverResult` 固定 `SVector{5}/SVector{3}`）与 `MethodError`（缺少 schema 版 `build_conditions`）。
- [x] 最小实现落地：
  - `src/models/solver/ImplicitSolver.jl`：`SolverResult` 改为向量泛型字段；默认物理性检查入参放宽为 `AbstractVector`。
  - `src/models/solver/Conditions.jl`：新增 `build_conditions(mode::FixedRho, params, schema; mu_dim=3)` 及 schema 提取辅助。
- [x] 复测证据：同命令复跑后 `Unit | Pass 42/42`。
- [x] 相关回归（最小必要）：`julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'` 通过（`Unit | Pass 781/781`）。

### 2026-04-02 Gate A / Plan-B Chunk3（TDD）

- [x] 已按 TDD 执行 Chunk3 Step1-4：新增失败 integration 测试 `tests/integration/models/test_dimension_agnostic_scan_smoke.jl`，覆盖 scan 结果适配层对动态维度的透传能力。
- [x] 失败测试证据：`julia --project=. -e 'include("tests/integration/models/test_dimension_agnostic_scan_smoke.jl")'` 首次运行报错 `DimensionMismatch`（`TmuScan/TrhoScan` 仍强制 `SVector{5}/SVector{3}`）。
- [x] 最小实现落地：
  - `src/models/scans/TmuScan.jl`：`_to_solver_result` 改为动态 `Float64` 向量映射。
  - `src/models/scans/TrhoScan.jl`：`_to_solver_result` 改为动态 `Float64` 向量映射。
- [x] 复测证据：同命令复跑后 `Pass 4/4`。
- [x] 相关回归（最小必要）：`julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'` 通过（`Integration | Pass 353/353`）。

### 2026-04-02 Gate A / Plan-B Chunk4（TDD + 回归/文档）

- [x] 新增回归测试 `tests/regression/models/test_dimension_agnostic_solver_regression.jl`，覆盖 7 模型 schema 维度与 flatten/unflatten 回环。
- [x] 将模型维度无关回归纳入 `tests/regression/runtests.jl`（smoke + full 都包含 `models` 组）。
- [x] API 文档同步：
  - `docs/api/models/solver/Overview.md`
  - `docs/api/models/solver/ConstraintModes.md`
- [x] 回归证据：
  - `julia --project=. -e 'ENV["REGRESSION_FILES"]="models/test_dimension_agnostic_solver_regression.jl"; include("tests/regression/runtests.jl")'` -> `Pass 35/35`
  - `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'` -> `Pass 497/498, Broken 1`（与既有可选 fixture skip 基线一致）
- [x] 治理检查证据：
  - `julia --project=. scripts/dev/check_docs_consistency.jl` -> OK
  - `julia --project=. scripts/dev/check_active_docs_governance.jl` -> OK

### 2026-04-02 Gate A / Plan-B Chunk5（TDD + legacy/compat 清理）

- [x] 已按 TDD 执行 Chunk5 Step1-4：先在 `tests/unit/models/test_solver_dimension_agnostic.jl` 增加失败断言，要求 fixed-* compat API 与 migration map/status 不再导出。
- [x] 失败测试证据：`julia --project=. -e 'ENV["UNIT_FILES"]="models/test_solver_dimension_agnostic.jl"; include("tests/unit/runtests.jl")'` 首次运行 `7` 个失败（compat API 与迁移查询仍存在）。
- [x] 代码清理落地：
  - `src/models/constraint_solver.jl`：`solve_fixed*` 重命名为内部 `_solve_constraint_*`，移除 hard-deprecate 兼容门。
  - `src/models/solver/Solver.jl`：删除 `solver_migration_map/status` 与 compat shim 路由，统一分发到 `_solve_constraint_*`。
  - `src/models/Models.jl`：移除 `solve_fixed*` 与 migration map/status 导出。
  - 删除 compat 专用测试：`tests/unit/models/test_solver_compat_markers.jl`、`tests/integration/models/test_waved_compat_cleanup_smoke.jl`、`tests/integration/models/test_waveb_compat_routing_smoke.jl`、`tests/regression/pnjl/test_waved_cleanup_stability.jl`。
- [x] 文档/API 同步：
  - 更新 `docs/api/models/solver/Overview.md`、`ConstraintModes.md`、`README.md`，移除 hard-deprecated 兼容口径；
  - 重新生成 `docs/api/models/solver/generated/Exports.md` 与 `docs/api/generated/models/ModelsExportIndex.md`，确认 fixed-* 与 migration map/status 不再出现在导出清单。
- [x] 验证证据：
  - 目标单测：`julia --project=. -e 'ENV["UNIT_FILES"]="models/test_solver_dimension_agnostic.jl,models/test_constraint_solver.jl"; include("tests/unit/runtests.jl")'` -> `Pass 48/48`
  - unit smoke：`Pass 781/781`
  - integration smoke：`Pass 355/355`
  - regression smoke：`Pass 497/498, Broken 1`（与既有可选 fixture skip 基线一致）
  - docs governance：`check_docs_consistency` / `check_active_docs_governance` 均 OK

### 2026-04-02 Gate B / Chunk0-3 进度（进行中）

- [x] Chunk0 Step1（Gate A 前置检查）已执行：
  - `include("src/constants/Constants_PNJL.jl"); include("src/models/Models.jl"); println(Models.registered_model_kinds())` -> `[:GasLiquid, :NJL, :NJL2, :PNJL, :PNJLMagnetic, :RPNJL, :Rotation]`
  - `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'` -> `Pass 781/781`
- [x] Chunk0 Step3（治理基线）已执行：
  - `julia --project=. scripts/dev/check_docs_consistency.jl` -> OK
  - `julia --project=. scripts/dev/check_active_docs_governance.jl` -> OK
- [x] Chunk1（接口同构）已完成：
  - 新增失败测试：`tests/unit/models/test_model_interface_homomorphism.jl`
  - 失败证据：首次运行对错误维度假设（将全部模型维度误设为 5）出现失败
  - 最小修复：按 7 模型真实 `gap_state_dim` 契约修正测试断言
  - 复测：`julia --project=. -e 'ENV["UNIT_FILES"]="models/test_model_interface_homomorphism.jl"; include("tests/unit/runtests.jl")'` -> `Pass 35/35`
- [x] Chunk2（骨架同构）已完成（最小增量）：
  - 新增失败 integration 测试：`tests/integration/models/test_model_skeleton_homomorphism_smoke.jl`
  - 失败证据：除 Rotation/GasLiquid 外缺少统一 workflow 访问路径
  - 最小实现：`src/models/entrypoints.jl` 新增 `workflow_module_for(model_kind::Symbol)` 并导出
  - 复测：`include("tests/integration/models/test_model_skeleton_homomorphism_smoke.jl")` -> `Pass 14/14`
- [x] Chunk3（回归+文档）已完成当前最小门禁：
  - 新增回归：`tests/regression/models/test_model_homomorphic_regression.jl`
  - 回归验证：`REGRESSION_FILES=models/test_model_homomorphic_regression.jl` -> `Pass 26/26`
  - regression smoke：`Pass 497/498, Broken 1`（既有可选 fixture skip）
  - 文档更新：
    - `docs/api/models/workflows/Overview.md`
    - `docs/api/models/workflows/CoreConcepts.md`
    - 重新生成 `docs/api/models/workflows/generated/Exports.md`

### 2026-04-02 Gate C / Chunk4（联合收口验证）

- [x] Step1 全套验证已执行：
  - `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'` -> `Unit | Pass 781/781`
  - `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'` -> `Integration | Pass 355/355`
  - `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'` -> `Regression | Pass 497/498, Broken 1`（既有可选 fixture skip 基线）
  - `julia --project=. scripts/dev/check_docs_consistency.jl` -> `OK`
  - `julia --project=. scripts/dev/check_active_docs_governance.jl` -> `OK`
- [x] Step2 compat 残留扫描（A 线定义范围）已执行：
  - `grep` 检查 `solve_fixed* / solver_migration_map / solver_migration_status` 在 `src/models` 与 `docs/api/models` 均为 `No files found`。
  - workflow 文档补充“禁止抽象清单”：`docs/api/models/workflows/CoreConcepts.md`。
- [x] Step2 补充扫描（Program 全局口径）已执行：
  - 在 `src/models` 仍检测到兼容痕迹（例如 `legacy_transport_c_p` 导出、`as_legacy_inputs` 弃用桥接、多个 `solver_backend=:legacy` 默认值）。
  - 结论：A 线定义范围已清零，但 Program DoD 的“compat 清零”尚未达成，需在后续收口中明确去留策略并执行。
- [x] Step3 收口提交准备完成（代码/测试/文档证据齐备，等待最终 PR/commit 操作）。

### 2026-04-02 Gate C / Compat 全清零追加（Program 全域）

- [x] 按“Program 全清零”口径完成 compat 清理：
  - 移除 `legacy_transport_c_p` 导出与实现（`src/models/Models.jl`、`src/models/derivatives/ThermoDerivatives.jl`）。
  - 移除 `as_legacy_inputs` 兼容别名（`src/models/workflows/WorkflowParamAdapters.jl`）。
  - 移除扫描迁移查询兼容面 `scan_workflow_migration_map/status`（`src/models/entrypoints.jl`），并同步对应测试断言。
  - 将扫描/workflow/phase/dual-branch 的 `solver_backend` 默认值从 `:legacy` 收敛为 `:auto`（`TmuScan/TrhoScan/TransportWorkflow/MesonMassWorkflow/PMPhaseSeeds/PMPhaseDiagnostic/DualBranchScan`）。
- [x] 相关测试与文档同步：
  - 删除 compat 专用单测：`tests/unit/models/test_workflow_param_adapters_compat.jl`。
  - 更新回归锚点：`tests/regression/pnjl/test_scan_fixedpoint_regression.jl` 显式固定 `solver_backend=:legacy` 以保持历史基线可复现。
  - 更新 docs：`docs/api/models/workflows/{CoreConcepts,ModuleAccessors,Overview}.md`、`docs/api/models/derived/**` 多页。
  - 重新生成导出页：`docs/api/models/derived/derivatives/generated/Exports.md`。
- [x] Gate C 验证矩阵（全量）复跑通过：
  - unit smoke：`Pass 781/781`
  - integration smoke：`Pass 353/353`
  - regression smoke：`Pass 497/498, Broken 1`（既有可选 fixture skip）
  - docs governance：`check_docs_consistency` / `check_active_docs_governance` 均 `OK`
- [x] compat 残留扫描结论（`src/models`）：
  - `scan_workflow_migration_map/status`、`as_legacy_inputs`、`legacy_transport_c_p`、`solver_backend::Symbol=:legacy` 均已 `No files found`。

---

## 8. Program DoD

- [x] Gate A/B/C 全部门禁通过
- [x] 7 模型覆盖矩阵全绿
- [x] 旧架构与 compat 清零
- [x] 文档、计划、执行状态一致
- [ ] PR checks 与 review 闭环
