# Solver 能力链路收敛 PR-C 任务单（错误语义 + 诊断契约）

## 1. 背景与目标

- 衔接关系：PR-A 完成 attempt 主线归一，PR-B 负责治理口径与 seed pool 归一。
- 本期目标：完成 Phase 4，统一错误语义与诊断契约，确保“失败可解释、诊断可回归”。
- 目标类型：可靠性与可观测性收敛。

## 2. 范围与非范围

### 2.1 本期范围（PR-C）

- [ ] C4.1 全链统一异常语义：`InterruptException` 重抛；普通异常降级为候选失败。
- [ ] C4.2 统一 diagnostics 字段：`error_kind/error_msg/attempt_origin/selection_reason` 等最小契约。
- [ ] C4.3 贯通 scan 输出消息，避免多套字符串拼接语义。
- [ ] C4.4 新增“异常注入回归 + 诊断契约回归”。

### 2.2 非范围

- [ ] 不新增新的业务模式或求解算法分支。
- [ ] 不在 PR-C 中再调整治理选择策略本身（只做语义与可观测性统一）。

## 3. 实施分解（可执行）

### 3.1 错误语义统一

- [ ] 在统一 attempt 执行链上固化错误分类（至少区分 interrupt / solver_error / constraint_error）。
- [ ] on_error 回调统一返回结构，不允许调用方私有格式。
- [ ] 清理历史 catch 分支中的语义不一致字段。

### 3.2 诊断契约统一

- [ ] 定义最小诊断字段集合与含义（文档化）。
- [ ] `ProblemSpec` / `Solver` / `ScanCommon` 输出对齐该集合。
- [ ] `_attach_solver_diagnostic`（若继续沿用）字段来源与优先级统一。

### 3.3 回归守护

- [ ] 新增异常注入测试（验证 interrupt 不被吞、普通异常被降级）。
- [ ] 新增诊断契约测试（字段存在性 + 语义稳定性）。
- [ ] 关键回归纳入 smoke 或确保可通过 `REGRESSION_FILES` 单跑。

## 4. 影响文件（预期）

- `src/models/solver/CandidateGovernance.jl`
- `src/models/solver/ProblemSpec.jl`
- `src/models/solver/Solver.jl`
- `src/models/scans/ScanCommon.jl`
- `tests/unit/models/test_candidate_governance_contract.jl`
- `tests/unit/models/test_solver_diagnostic_contract.jl`
- `tests/unit/models/test_scan_common.jl`
- `tests/regression/models/` 下新增诊断/异常回归文件

## 5. 验证命令（建议）

```powershell
julia --project=. -e 'ENV["UNIT_FILES"]="models/test_candidate_governance_contract.jl;models/test_solver_diagnostic_contract.jl;models/test_scan_common.jl"; include("tests/unit/runtests.jl")'
julia --project=. -e 'ENV["REGRESSION_FILES"]="models/test_solver_attempt_engine_convergence_regression.jl;pnjl/test_constraint_selection_regression.jl"; include("tests/regression/runtests.jl")'
```

## 6. PR-C DoD

- [ ] 错误语义全链一致，`InterruptException` 行为受测试保护。
- [ ] diagnostics 字段集合稳定、可追踪、可回归。
- [ ] scan 输出消息与 solver 诊断语义一致。
- [ ] 相关 unit/regression 全通过，任务单记录与实现一致。

## 7. 风险与回退

- 风险：诊断字段统一导致现有脚本解析失败。
  - 缓解：保持向后兼容窗口（字段别名或兼容映射），并在任务单注明迁移点。
- 风险：异常分类过细引入实现复杂度。
  - 缓解：先最小可用分类，后续增量细化。
- 回退：优先回退输出层适配，不回退核心 interrupt 语义保障。
