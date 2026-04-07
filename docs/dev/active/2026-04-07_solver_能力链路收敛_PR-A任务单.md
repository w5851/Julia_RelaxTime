# Solver 能力链路收敛 PR-A 任务单（attempt 引擎归一）

## 1. 背景与目标

- 背景：当前 solver 主链在 `Solver.jl`、`ProblemSpec.jl`、`ConstraintSolver.jl`、`ScanCommon.jl` 之间仍存在“同类能力分散实现”现象，尤其是多候选尝试编排（attempt orchestration）存在平行实现。
- 目标：先完成 PR-A（Phase 1），将“attempt 执行能力”收敛为唯一引擎，降低行为漂移与后续演进成本。
- 目标类型：能力归一化（平台化），非新增业务功能。

## 2. 范围与非范围

### 2.1 本期范围（PR-A）

- [ ] 统一 attempt 执行入口：所有多候选尝试统一走 `execute_attempt_pool`（或同职责升级接口）。
- [ ] `Solver.solve_multi` 改为“构建 attempt plan + 调统一引擎”，移除本地平行 for/try/select 编排。
- [ ] `ProblemSpec` 保留 plan 构建与参数注入职责，避免重复执行编排细节。
- [ ] 统一 early-stop / evaluate-all 语义，确保单点定义与跨层一致。

### 2.2 非范围（延后到 PR-B / PR-C）

- [ ] 不在 PR-A 内完成 seed pool 组装入口全面收敛（Phase 3）。
- [ ] 不在 PR-A 内完成 success 判定与 candidate schema 全面归一（Phase 2 全量）。
- [ ] 不在 PR-A 内完成诊断输出契约与 scan message 格式全面统一（Phase 4）。

## 3. 现状盘点（已确认）

- 已具备统一能力雏形：`src/models/solver/CandidateGovernance.jl` 提供 `execute_attempt_pool`。
- 已完成关键异常语义修复：`InterruptException` 重抛、普通异常透传 `err` 到 `on_error`。
- 当前分散点：
  - `Solver.solve_multi` 仍保留本地候选循环与选择逻辑。
  - `ProblemSpec` 与 `solve_multi` 存在并行 attempt 编排路径。
  - early-stop/evaluate-all 语义存在多处定义风险。

## 4. PR-A 最小改动清单（实施级）

### 4.1 T0 基线与护栏（先落测试）

- [x] T0.1 选取 6 个代表点（覆盖 `FixedMu/FixedRho/FixedAsymmetricRho/FixedEntropy/FixedSigma`）。
- [x] T0.2 固化关键字段基线：`converged/residual_norm/selection_reason/seed_index/failed_constraints`。
- [x] T0.3 新增守护测试文件：`tests/regression/models/test_solver_attempt_engine_convergence_regression.jl`。
- [x] T0.4 将该回归纳入 `tests/regression/runtests.jl` smoke 列表。

#### T0 基线快照（2026-04-07）

- `fixedmu_100_250`：`(converged=true, residual_norm=2.3562268645305217e-14, selection_reason=:pressure_max_under_constraints, seed_index=-1, failed_constraints=[])`
- `fixedmu_140_320`：`(converged=true, residual_norm=2.2392747047217338e-15, selection_reason=:pressure_max_under_constraints, seed_index=-1, failed_constraints=[])`
- `fixedrho_110_0p2`：`(converged=true, residual_norm=1.8366093460261402e-15, selection_reason=:pressure_max_under_constraints, seed_index=-1, failed_constraints=[])`
- `fixedasym_110_0p05`：`(converged=false, residual_norm=Inf, selection_reason=:no_candidate_passed_constraints, seed_index=-1, failed_constraints=[:solver_failed])`
- `fixedentropy_120_0p5`：`(converged=false, residual_norm=Inf, selection_reason=:no_candidate_passed_constraints, seed_index=-1, failed_constraints=[:solver_failed])`
- `fixedsigma_120_10`：`(converged=false, residual_norm=Inf, selection_reason=:no_candidate_passed_constraints, seed_index=-1, failed_constraints=[:solver_failed])`

### 4.2 T1 attempt 引擎归一

- [ ] T1.1 在 `CandidateGovernance` 明确统一接口契约：
  - 输入：attempt plan（seed/method/fallback metadata）
  - 回调：`evaluate_attempt`、`on_error`
  - 输出：标准 candidate 列表（含 diagnostics）
- [x] T1.2 重写 `Solver.solve_multi` 路径，改用 `execute_attempt_pool`。
- [x] T1.3 清理 `Solver.solve_multi` 中重复编排逻辑（for/try/catch/本地选择分支）。
- [ ] T1.4 `ProblemSpec` 侧只保留 plan 构建和参数映射，执行逻辑复用统一引擎。
- [x] T1.5 对齐 early-stop/evaluate-all 默认值与透传规则，确保跨层单义。

#### T1 当前进展记录

- 已完成：`src/models/solver/Solver.jl` 的 `solve_multi(FixedMu/FixedRho/FixedAsymmetricRho)` 两条主路径改为统一调用 `execute_attempt_pool`，移除本地 for/try/catch 重复编排。
- 已完成：`solve_multi` 新增 `evaluate_all_attempts` 统一透传（默认 `true`），并在 forwarded kwargs 中单点剔除，避免多层重复解释。
- 验证通过：
  - `tests/unit/models/test_candidate_governance_contract.jl`
  - `tests/unit/models/test_solver.jl`
  - `tests/regression/models/test_solver_attempt_engine_convergence_regression.jl`
  - `REGRESSION_FILES=models/test_solver_attempt_engine_convergence_regression.jl;pnjl/test_constraint_selection_regression.jl`

### 4.3 T1 验证与回归

- [ ] T1.6 更新/新增单测：
  - `tests/unit/models/test_candidate_governance_contract.jl`
  - `tests/unit/models/test_solver.jl`
  - 必要时补 `tests/unit/models/test_problem_spec_contract.jl`
- [ ] T1.7 跑 targeted regression（至少含新建收敛回归 + constraint selection 相关回归）。
- [ ] T1.8 对比基线样本，确认关键字段一致（或记录可解释差异）。

## 5. 影响文件列表（PR-A 预期）

- 核心实现：
  - `src/models/solver/CandidateGovernance.jl`
  - `src/models/solver/Solver.jl`
  - `src/models/solver/ProblemSpec.jl`
- 可能联动：
  - `src/models/solver/ConstraintSolver.jl`
  - `src/models/scans/ScanCommon.jl`
- 测试：
  - `tests/unit/models/test_candidate_governance_contract.jl`
  - `tests/unit/models/test_solver.jl`
  - `tests/unit/models/test_problem_spec_contract.jl`
  - `tests/regression/models/test_solver_attempt_engine_convergence_regression.jl`（新增）
  - `tests/regression/runtests.jl`

## 6. 验证命令（PR-A）

```powershell
julia --project=. -e 'ENV["UNIT_FILES"]="models/test_candidate_governance_contract.jl;models/test_solver.jl;models/test_problem_spec_contract.jl"; include("tests/unit/runtests.jl")'
julia --project=. -e 'ENV["REGRESSION_FILES"]="models/test_solver_attempt_engine_convergence_regression.jl;pnjl/test_constraint_selection_regression.jl"; include("tests/regression/runtests.jl")'
```

如需快速点对点验证，可额外使用：

```powershell
julia --project=. -e 'include("tests/unit/models/test_candidate_governance_contract.jl")'
```

## 7. 验收标准（PR-A DoD）

- [ ] 代码层：`solve_multi` 不再保留平行 attempt 编排实现。
- [ ] 能力层：attempt 执行路径唯一（统一引擎），无重复语义分支。
- [ ] 测试层：新增守护回归 + 相关 unit/regression 全通过。
- [ ] 行为层：基线样本关键字段保持一致，差异有书面解释。
- [ ] 文档层：本任务单状态与关键决策同步更新。

## 8. 里程碑与提交建议

- M1（护栏就绪）：完成 T0.x，提交 `test(regression): add solver attempt-engine convergence guard`
- M2（主改完成）：完成 T1.1~T1.5，提交 `refactor(models/solver): unify solve_multi attempt orchestration`
- M3（验证收口）：完成 T1.6~T1.8，提交 `test(models/regression): align attempt-engine convergence assertions`

## 9. 风险与回退方案

- 风险 R1：attempt 顺序变化引入 selection 差异。
  - 缓解：固定 plan 顺序 + 基线样本对比。
- 风险 R2：early-stop/evaluate-all 语义收敛导致历史隐式行为变化。
  - 缓解：显式参数透传测试，补契约断言。
- 风险 R3：回归覆盖不足导致漂移未被捕获。
  - 缓解：新增 regression 守护并纳入 smoke 或常用 REGRESSION_FILES。
- 回退：保留小步提交，若行为不可解释可按里程碑回退到 M1 基线。

## 10. 后续路线（占位）

- PR-B：Phase 2 + 3（治理口径与 seed pool 归一）
- PR-C：Phase 4（错误语义 + 诊断契约）
