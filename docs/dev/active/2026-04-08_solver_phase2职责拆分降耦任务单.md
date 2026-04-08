# Solver Phase 2 职责拆分与降耦任务单（编排/后处理/诊断解耦）

## 1. 目标

- 在 Phase 1 单主链已建立的基础上，拆分 solver 巨石职责，落实“高内聚、低耦合”。
- 将 `ProblemSpecOrchestrator` 从“调度 + 物理后处理 + 诊断拼装 + kwargs 管道”收敛为“纯调度编排层”。
- 建立稳定的数据契约，降低 `phase/scans` 对 solver 内部字段的脆弱依赖。

### 1.1 前置与后置

- 前置：`docs/dev/active/2026-04-08_solver_phase1主链收敛任务单.md` 全部 DoD 达成。
- 后置：为 `docs/dev/active/2026-04-08_solver_phase3能力边界固化与治理门禁任务单.md` 提供稳定契约基础。

## 2. 范围与非范围

### 2.1 本期范围（Phase 2）

- [ ] 抽取统一热力学后处理模块，消除重复实现。
- [ ] 引入运行时配置类型，替代关键路径 `Dict{Symbol,Any}` 管道。
- [ ] 诊断结构类型化并提供兼容映射层。
- [ ] 固化 attempt/governance 的策略边界，mode executor 仅保留残差与约束求解职责。

### 2.2 非范围（本期不做）

- [ ] 不在本期完成所有历史字段清理（先保兼容，Phase 3 再收口）。
- [ ] 不在本期调整 physics 公式或数值策略阈值（避免与架构改动耦合）。
- [ ] 不在本期删除 `ConstraintSolverFixed*.jl`（只做职责瘦身和调用对齐）。

## 3. 目标架构与依赖方向

- 目标依赖：
  - `mode executor` 只产出“原始求解解 + 必要元信息”。
  - `ThermoPostprocess` 统一负责 `omega/pressure/rho/entropy/energy/masses/residual_norm` 计算。
  - `SolverDiagnosticsTypes` 统一负责诊断对象构造与序列化映射。
  - `ProblemSpecOrchestrator` 只负责调度、候选治理、调用上述模块。
- 期望方向：`solver_core <- strategy/governance <- phase|scans|workflows`。

## 4. 文件级任务分解（可执行）

### Batch-A：统一后处理抽取（先解重复）

- [x] A1 新增 `src/models/solver/ThermoPostprocess.jl`
  - 提供统一入口（建议）：
    - `compute_thermo_from_solution(...)`
    - `compute_residual_norm_from_solution(...)`
    - `build_solver_candidate(...)`
  - 约束：输出字段命名与当前 `solve_constraint` 结果保持一致。

- [ ] A2 改造 `src/models/solver/ConstraintSolverCommon.jl`
  - 保留纯共性工具（约束规则、候选选择、seed 辅助）。
  - 将重复热力学后处理逻辑迁移到 `ThermoPostprocess.jl`。

- [x] A3 改造 `src/models/solver/ProblemSpecOrchestrator.jl`
  - 移除内联 thermodynamic 计算块，统一改为调用 `ThermoPostprocess`。
  - 确保 `FixedRho` joint 路径与 non-rho 路径后处理一致。

### Batch-B：运行时配置类型化（先试点后推广）

- [x] B1 新增 `src/models/solver/SolverRuntimeConfig.jl`
  - 第一批类型：`FixedRhoRuntimeConfig`, `FixedEntropyRuntimeConfig`。
  - 包含字段校验与默认值解析函数，减少 `Dict` 删除/重写。

- [x] B2 改造 `src/models/solver/ProblemSpecOrchestrator.jl`
  - 在 `FixedRho`、`FixedEntropy` 路径先使用类型化 config。
  - 其他 mode 暂保留旧字典管道，但通过适配器统一入口。

- [ ] B3 改造 `src/models/solver/Solver.jl`
  - 入口处把 kwargs 解析成运行时 config，再传给 orchestrator。
  - 保持对旧 kwargs 的兼容（含报错文案和别名映射）。

### Batch-C：诊断结构类型化（稳定对外契约）

- [x] C1 新增 `src/models/solver/SolverDiagnosticsTypes.jl`
  - 建议类型：
    - `SolverDiagnosticSummary`
    - `SolverDiagnosticCandidate`
    - `SolverDiagnosticFull`
  - 提供 `to_namedtuple`（兼容旧调用方）和 `from_candidate` 构造函数。

- [x] C2 改造 `src/models/solver/SolverDiagnostics.jl`
  - 将当前 NamedTuple 拼装改为“类型对象 -> 兼容映射”。
  - `diagnostic_level=:summary|:full` 的输出语义保持不变。

- [x] C3 适配上层读法（只做兼容，不做行为变化）
  - 文件：`src/models/phase/PMPhaseDiagnostic.jl`, `src/models/scans/TmuScan.jl`, `src/models/scans/TrhoScan.jl`
  - 目标：允许读取类型化诊断或兼容 NamedTuple，不绑定 solver 内部临时字段名。

### Batch-D：策略边界收敛（治理与求解解耦）

- [x] D1 改造 `src/models/solver/CandidateGovernance.jl`
  - 固化 selector 输入契约（必须字段集）。
  - 失败候选和成功候选构造统一使用同一工厂函数。

- [x] D2 改造 `src/models/solver/ProblemSpecOrchestrator.jl`
  - attempt 评估只关心“调用 solve + 收到候选 + 交给 governance”，
    不在 orchestrator 内混入业务字段拼装细节。

- [x] D3 校对 `src/models/solver/ConstraintSolverFixed*.jl`
  - 确保 mode executor 只做 residual/求根，不再重复治理与诊断拼装。

## 5. 逐函数改动清单（函数名 + 预期 diff 类型）

### 5.1 `src/models/solver/ProblemSpecOrchestrator.jl`

- [x] `_fixedrho_joint_problem_spec_forward_solve(...)`
  - diff 类型：`后处理下沉`
  - 从函数内部剥离热力学计算，改为调用 `ThermoPostprocess`。

- [x] `_governed_nonrho_problem_spec_forward_solve(...)`
  - diff 类型：`候选构造收敛`
  - 候选字段尽量通过统一工厂构造，减少手工散装字段。

- [x] `_execute_governed_attempt_plan(...)`
  - diff 类型：`边界收敛`
  - 只负责流程控制，候选标准化交由治理/诊断模块。

### 5.2 `src/models/solver/SolverDiagnostics.jl`

- [x] `_solver_diagnostic_from_candidate(...)`
  - diff 类型：`类型构造`
  - 改为先构造诊断类型，再输出兼容视图。

- [x] `_attach_solver_diagnostic(...)`
  - diff 类型：`兼容层`
  - 对外可继续返回旧字段结构，但内部统一类型化。

### 5.3 `src/models/solver/ConstraintSolverCommon.jl`

- [ ] `_compute_mode_thermo_quantities(...)`
  - diff 类型：`迁移/瘦身`
  - 若迁移后可删除，则保留过渡 wrapper 并标注来源模块。

- [ ] `_compose_mode_residual_norm(...)`
  - diff 类型：`复用统一入口`
  - 与 `ThermoPostprocess` 中 residual 逻辑对齐。

## 6. 验收标准（DoD）

- [ ] `ProblemSpecOrchestrator.jl` 不再包含重复 thermodynamic 公式实现。
- [ ] `FixedRho/FixedEntropy` 关键路径已使用类型化 runtime config。
- [x] diagnostics 内部已类型化，外部行为兼容（phase/scans 无需大改）。
- [x] 上层不依赖 solver 私有临时字段名（至少完成兼容隔离层）。

## 7. 验证计划（每个 Batch 完成后执行）

- [ ] `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
- [ ] `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
- [ ] `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
- [ ] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [ ] `julia --project=. scripts/dev/check_active_docs_governance.jl`

建议附加点测（Phase 2 推荐）：

- [ ] `tests/integration/relaxtime/test_transport_workflow_smoke.jl`
- [ ] `tests/integration` 中涉及 `PMPhase` 与 `Tmu/Trho` 的 smoke 用例（按仓库现有入口挑选）。

## 8. 风险与回滚点

- [ ] 风险 A：后处理迁移导致数值漂移（尤其 entropy/rho_norm）。
  - 回滚点：保留旧实现 wrapper，对单 mode 快速切回。

- [ ] 风险 B：诊断类型化引入上层读取断裂。
  - 回滚点：保留 `to_namedtuple` 兼容输出，先不中断旧字段。

- [ ] 风险 C：配置类型化覆盖不全，导致 kwargs 丢失。
  - 回滚点：未覆盖字段通过 passthrough 适配层转发。

- [ ] 风险 D：治理逻辑收敛改变候选排序。
  - 回滚点：冻结 selector 比较器与质量标签规则，逐步替换实现而不改规则。

## 9. 建议提交粒度

- [ ] Commit 1：新增 `ThermoPostprocess.jl` + 接入一个 mode（最小闭环）。
- [ ] Commit 2：`ProblemSpecOrchestrator.jl` 后处理全面改接入。
- [ ] Commit 3：新增 `SolverRuntimeConfig.jl` 并完成 FixedRho/FixedEntropy 试点。
- [ ] Commit 4：新增 `SolverDiagnosticsTypes.jl` + 兼容映射层。
- [ ] Commit 5：phase/scans 兼容适配 + 全量 smoke/regression 验证。

## 10. 执行记录

- [x] 2026-04-08：基于 Phase 1 任务单扩展，产出 Phase 2 职责拆分任务单。
- [x] 2026-04-08：与 Phase 1 / Phase 3 文档完成衔接字段校对（前置与后置关系明确）。
- [x] 2026-04-08：完成 Batch-A 最小闭环（A1 + A3）。
  - 新增 `src/models/solver/ThermoPostprocess.jl`，提供 `compute_thermo_from_solution(...)`、`compute_residual_norm_from_solution(...)`、`build_solver_candidate(...)`。
  - `ProblemSpecOrchestrator._fixedrho_joint_problem_spec_forward_solve(...)` 去除内联热力学计算，统一接入 `ThermoPostprocess`。
  - 测试：`tests/unit/models/test_solver.jl` 新增统一后处理入口可用性测试并通过；`tests/integration/models/test_solver_auto_backend_semantic_parity.jl` 通过。
- [x] 2026-04-08：完成 Batch-B 试点（B1 + B2，FixedRho/FixedEntropy）。
  - 新增 `src/models/solver/SolverRuntimeConfig.jl`：`FixedRhoRuntimeConfig`、`FixedEntropyRuntimeConfig` 及解析函数。
  - `ProblemSpecOrchestrator` 在 `FixedRho` 与 `FixedEntropy` 路径接入类型化 config，并保持 legacy kwargs 兼容。
  - 修复参数泄漏：`evaluate_all_attempts` 由 ProblemSpec 前向 kwargs 剔除，避免落入 `nlsolve`。
  - 测试：`tests/unit/models/test_solver.jl`（56/56）与 `tests/integration/models/test_solver_auto_backend_semantic_parity.jl`（含 fixedmu parity guard）通过。
- [x] 2026-04-08：完成 Batch-C 核心（C1 + C2）。
  - 新增 `src/models/solver/SolverDiagnosticsTypes.jl`，引入 `SolverDiagnosticSummary` / `SolverDiagnosticCandidate` / `SolverDiagnosticFull` 与 `to_namedtuple` 兼容映射。
  - `src/models/solver/SolverDiagnostics.jl` 改为“内部类型化构造 + 对外 NamedTuple 兼容输出”，`diagnostic_level=:summary|:full` 语义保持不变。
  - 测试：`tests/unit/models/test_solver.jl` 新增类型化诊断兼容测试并通过（60/60）；`tests/integration/models/test_solver_auto_backend_semantic_parity.jl` 通过。
- [x] 2026-04-08：完成 Batch-D（D1 + D2 + D3）。
  - `CandidateGovernance.jl` 新增 selector 输入字段契约常量与校验，新增统一候选工厂 `build_governance_candidate(...)`。
  - `ProblemSpecOrchestrator.jl` attempt 执行路径改为统一使用候选工厂，削弱编排层字段拼装职责。
  - `Solver.jl` multi-seed 尝试路径同步接入统一候选工厂；`ConstraintSolverFixed*.jl` 校对确认仍聚焦 residual/求根。
  - 测试：`tests/unit/models/test_candidate_governance_contract.jl`（61/61）、`tests/unit/models/test_solver.jl`（60/60）、`tests/integration/models/test_solver_auto_backend_semantic_parity.jl`（通过）。
- [x] 2026-04-08：完成 C3 上层兼容适配收口。
  - `src/models/phase/PMPhaseDiagnostic.jl` 适配类型化诊断输入：支持 `SolverDiagnostic*` 类型与 NamedTuple 统一视图，不绑定 solver 内部临时字段。
  - 补充回归：`tests/unit/models/test_pm_phase_diagnostic.jl` 新增 typed diagnostic 兼容用例。
  - 验证：`tests/integration/models/test_phase_solver_diagnostic_adapter_smoke.jl` 通过，确保 phase 侧诊断读取保持兼容。
- [ ] 待补：各 Batch 的改动 PR/commit、测试结果、回滚演练记录。
