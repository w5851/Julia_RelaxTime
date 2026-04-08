# Solver Phase 1 主链收敛任务单（ProblemSpec 单主链）

> 当前状态：**已完成（含收口项）**

## 1. 目标

- 在不改变外部能力语义的前提下，将 solver 收敛到单一主链：
  `solve* -> solve_constraint -> ProblemSpec.forward_solve -> mode executor`。
- 消除 `FixedMu` 的双路径分流（ProblemSpec 路径 vs direct `_solve_constraint_fixedmu`）。
- 对遗留旁路开关做显式治理（报错/弃用提示），避免静默进入非主链。

### 1.1 与后续阶段衔接

- 本文档是三阶段计划的入口阶段，后续衔接文档：
  - `docs/dev/active/2026-04-08_solver_phase2职责拆分降耦任务单.md`
  - `docs/dev/active/2026-04-08_solver_phase3能力边界固化与治理门禁任务单.md`
- Phase 1 完成前，不推进 Phase 2 的结构迁移与 Phase 3 的兼容清理。

## 2. 范围与非范围

### 2.1 本期范围（Phase 1）

- [x] 收敛 `src/models/solver/Solver.jl` 的 `solve_constraint(::FixedMu, ...)` 到 ProblemSpec 主链。
- [x] 统一 legacy 旁路开关拦截策略（含错误文案）。
- [x] 同步更新 solver 主流程注释，确保代码内文档与实现一致。
- [x] 验证 `phase/scans` 仅依赖 `solve_constraint` 公共入口，不依赖 direct `_solve_constraint_fixed*` 调用。

### 2.2 非范围（本期不做）

- [ ] 不在 Phase 1 重写 `ProblemSpecOrchestrator.jl` 的职责拆分（留给 Phase 2）。
- [ ] 不在 Phase 1 引入新的诊断类型结构（留给 Phase 2/3）。
- [ ] 不在 Phase 1 删除 `ConstraintSolverFixed*.jl` 文件本体（先断主入口再清理）。

## 3. 逐函数改动清单（函数名 + 预期 diff 类型）

### 3.1 `src/models/solver/Solver.jl`

- [x] `solve_constraint(model::AbstractQCDModel, mode::FixedMu, T_fm::Real; ...)`
  - diff 类型：`逻辑删除 + 主链强制 + 兼容字段保留`
  - 改动要点：
    1) 删除 `fixedmu_use_problem_spec=false` 时直达 `_solve_constraint_fixedmu(...)` 的分支。
    2) 无论是否传 `problem_spec`，均走 `_solve_with_problem_spec_default(...)`。
    3) 保留返回字段 `fixedmu_problem_spec_active`（兼容上游读取），并确保值为 `true` 或按 raw 回传。

- [ ] `_solve_with_problem_spec_default(model, mode, T_fm, kwargs)`
  - diff 类型：`防御性增强 + 文案统一`
  - 改动要点：
    1) 统一 legacy 开关错误文案：明确“已移除 + 请使用 ProblemSpec 主链”。
    2) 确保 `problem_spec` 类型校验与转发规则保持稳定。

- [ ] `_strip_forward_kwargs(kwargs, blocked)`（如需）
  - diff 类型：`小幅兼容`
  - 改动要点：仅在必要时调整被剔除键，避免遗留键泄漏到下游造成行为分叉。

- [x] 文件头或 `solve_constraint` 相关 docstring
  - diff 类型：`注释对齐`
  - 改动要点：明确唯一主流程，不再描述可选 legacy 分流。

### 3.2 `src/models/solver/ConstraintSolver.jl`

- [x] 模块说明注释
  - diff 类型：`注释对齐`
  - 改动要点：强调该文件是 mode executor 聚合入口，但主调度入口是 `Solver.jl -> ProblemSpec`。

### 3.3 `src/models/solver/ProblemSpec.jl`

- [x] `build_problem_spec(mode::ConstraintMode; kwargs...)` 上方注释
  - diff 类型：`注释增强`
  - 改动要点：补充“作为唯一 forward_solve 注册入口”的职责声明。

## 4. 兼容性护栏（上层模块）

### 4.1 `src/models/scans/TmuScan.jl`

- [x] 检查 `_solve_with_models(...)` 调用路径
  - 目标：确认仅通过 `Main.Models.solve_constraint(...)`；不依赖 direct `_solve_constraint_fixed*`。
- [x] 校验对返回字段的依赖
  - 目标：若读取 `fixedmu_problem_spec_active`，Phase 1 保持兼容回传。

### 4.2 `src/models/scans/TrhoScan.jl`

- [x] 检查 `_solve_with_models(...)` 中 `problem_spec`/`semantic_mode` 分支
  - 目标：确保 Phase 1 后行为不变（仍由 `solve_constraint` 接管分发）。

### 4.3 `src/models/phase/PMPhaseDiagnostic.jl`

- [x] 检查 `solve_constraint(FixedMu)` 结果字段消费
  - 目标：确认其依赖字段在 Phase 1 保持不变，尤其是 `diagnostic`、`converged`、`residual_norm`。

## 5. 验证计划（必须执行）

- [x] `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
- [x] `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
- [x] `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
- [x] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [x] `julia --project=. scripts/dev/check_active_docs_governance.jl`

重点观察点：

- [x] `TmuScan/TrhoScan/PMPhase*` 的 `converged` 分布是否异常漂移。
- [x] `residual_norm` 与 `selection_reason` 是否出现系统性偏移。
- [x] legacy 开关触发时是否按预期报错（而非静默降级）。

## 6. 风险、回滚点与应对

- [x] 风险 A：历史调用显式传 `fixedmu_use_problem_spec=false`。
  - 应对：抛出 `ArgumentError`，文案给出替代用法。
  - 回滚点：仅临时恢复参数容忍（warning），不恢复 direct 分流逻辑。

- [ ] 风险 B：自定义 `problem_spec` 依赖旧 kwargs 清洗顺序。
  - 应对：在 `_solve_with_problem_spec_default` 内补兼容映射。
  - 回滚点：只回滚参数映射，不回滚主链统一。

- [x] 风险 C：phase/scans 对过渡字段的隐式依赖。
  - 应对：保留过渡字段回传（如 `fixedmu_problem_spec_active`）。
  - 回滚点：加适配层映射，不改求解内核。

## 7. 建议提交粒度（便于 bisect/回滚）

- [x] Commit 1：`Solver.jl` 主链收敛（只动 FixedMu 入口分流）。
- [ ] Commit 2：legacy 开关治理 + 注释对齐。
- [ ] Commit 3：上层兼容修补（若测试暴露依赖）。
- [ ] Commit 4：文档/治理脚本收口。

## 8. DoD（Phase 1）

- [x] `solve_constraint` 对全部 mode 仅保留 ProblemSpec 主链入口。
- [x] 无调用方可通过 legacy 开关进入 direct `_solve_constraint_fixed*` 路径。
- [x] `TmuScan/TrhoScan/PMPhase*` smoke 与回归在容差内通过。
- [x] 代码内主流程注释与实现一致。
- [x] 为 Phase 2 预留兼容字段（如 `fixedmu_problem_spec_active`）仍可稳定回传。

## 9. 执行记录

- [x] 2026-04-08：根据 solver 全量分析结论，生成 Phase 1 可执行任务单（逐函数粒度）。
- [x] 2026-04-08：完成 Phase 1 首轮实现并通过验证。
  - 代码：`src/models/solver/Solver.jl`, `src/models/solver/ConstraintSolver.jl`, `src/models/solver/ProblemSpec.jl`
  - 测试：`tests/unit/models/test_solver.jl`, `tests/unit/models/test_solver_diagnostic_contract.jl`, `tests/integration/models/test_solver_auto_backend_semantic_parity.jl`
  - 验证摘要：Unit smoke 772/772；Integration smoke 456/456；Regression smoke 511 pass + 1 broken（既有可选跳过项）
- [x] 2026-04-08：完成 Phase 1 收口项。
  - `_solve_with_problem_spec_default(...)` legacy 开关报错文案统一为“已移除 + ProblemSpec 替代路径”。
  - 完成风险 B 专项回归（`test_problem_spec_contract.jl` + solver 相关 unit/integration 目标集）。
  - docs 治理检查通过：`check_docs_consistency.jl`、`check_active_docs_governance.jl`。

## 10. 剩余项（未完成即不宣称 Phase 1 完结）

- [x] `_solve_with_problem_spec_default(...)` 文案治理与兼容映射收口（任务 3.1 第二项已完成）。
- [x] `_strip_forward_kwargs(...)` 评估结论：Phase 1 不需要额外兼容调整（保持现状，记录为按需关闭）。
- [x] 风险 B（自定义 `problem_spec` kwargs 顺序依赖）已完成专项回归验证（见验证与执行记录）。
- [x] 提交粒度收口：本阶段代码/测试/文档已形成可提交批次（后续按仓库节奏提交）。

结论：Phase 1 已达到 DoD，可进入 Phase 2。
