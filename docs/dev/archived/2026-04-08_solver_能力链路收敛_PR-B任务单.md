---
title: Solver 能力链路收敛 PR-B 任务单（治理口径 + seed pool 归一）
archived: true
original: docs/dev/active/2026-04-08_solver_能力链路收敛_PR-B任务单.md
archived_date: 2026-04-08
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Solver 能力链路收敛 PR-B 任务单（治理口径 + seed pool 归一）

## 1. 背景与目标

- 衔接关系：PR-A 已完成 attempt 引擎主线收敛（`execute_attempt_pool` 单点执行）。
- 本期目标：完成 Phase 2 + Phase 3，消除“治理口径分散”和“seed pool 组装分散”。
- 目标类型：能力归一化与契约稳定化，不新增业务功能。

## 2. 范围与非范围

### 2.1 本期范围（PR-B）

- [x] B2.1 统一 candidate 最小 schema：`converged/residual_norm/hard_constraint_ok/failed_constraints/pressure/seed_index`。
- [x] B2.2 统一 success 判定入口，各层不再各写一套判定。
- [ ] B2.3 统一 selector 调用点，避免二次筛选与语义漂移。
- [x] B2.4 清理 `Solver.jl` 与 `ScanCommon.jl` 中重复治理判断。
- [x] B3.1 提取 mode-aware seed pool builder（单一入口）。
- [x] B3.2 `solve_multi`、`ProblemSpec`、scan 侧统一复用 builder。
- [x] B3.3 统一 seed 去重/顺序/continuity 注入规则并文档化。

### 2.2 非范围（延后到 PR-C）

- [ ] 不在 PR-B 内完成 error_kind/error_msg 诊断字段全链路标准化。
- [ ] 不在 PR-B 内完成 scan 输出消息格式统一与异常注入回归全量覆盖。

## 3. 实施分解（可执行）

### 3.1 治理契约收敛

- [x] 提供 `normalize_candidate_for_governance(...)`（或等价入口），把不同来源 candidate 归一到最小 schema。
- [x] 提供 `evaluate_candidate_success(...)`（或等价入口）作为唯一 success 判定。
- [x] 在 `Solver` / `ProblemSpec` / `ScanCommon` 替换本地治理拼装逻辑，统一走上述入口。

#### 3.1 进展记录（2026-04-08）

- 新增统一入口：`CandidateGovernance.evaluate_candidate_success` 与 `CandidateGovernance.normalize_governance_candidate`。
- `Solver.solve_multi(FixedMu/FixedRho/FixedAsymmetricRho)` 已改为统一调用归一化与 success 入口，移除局部口径分叉。
- 已补单测：`tests/unit/models/test_candidate_governance_contract.jl` 新增 helper 契约测试。
- 已回归验证：
  - `tests/unit/models/test_candidate_governance_contract.jl`
  - `tests/unit/models/test_solver.jl`
  - `REGRESSION_FILES=models/test_solver_attempt_engine_convergence_regression.jl;pnjl/test_constraint_selection_regression.jl`
- 追加收敛（本轮）：
  - `ProblemSpec._execute_governed_attempt_plan` 接入统一 `normalize_governance_candidate/evaluate_candidate_success`。
  - `ScanCommon.attempt_with_candidates` 接入统一 governance helper，移除本地 success 口径分叉。

### 3.2 seed pool 收敛

- [ ] 提取 `build_seed_pool(mode, seed_guess, kwargs...)`（命名可按现有风格调整）。
- [ ] 合并当前 `_build_default_seed_candidates` 的重复使用场景与散点 fallback。
- [ ] 固化顺序稳定性（主 seed、显式候选、默认候选、去重策略）。

#### 3.2 进展记录（2026-04-08）

- 新增统一入口：`CandidateGovernance.build_seed_pool(mode; primary_seed, extra_seed_pool, provided_seed_pool, default_seed_pool, seed_extend)`。
- `ProblemSpec._build_governed_attempt_plan` 已改为复用 `build_seed_pool`，固定顺序为主 seed -> extra -> provided -> default，并统一去重键。
- `Solver.solve_multi(FixedMu/FixedRho/FixedAsymmetricRho)` 与 `TmuScan/TrhoScan` 已改为复用 `build_seed_pool`。
- 去重/顺序/continuity 规则：continuity(若有) 作为 primary；随后 explicit/provided；最后默认候选；统一按 12 位 round-key 去重。

### 3.3 文档与契约同步

- [x] 在本任务单中补“能力 -> 唯一实现点”矩阵。
- [ ] 若导出 API 发生变化，更新 `docs/api/` 对应页面（仅在必要时）。

## 4. 影响文件（预期）

- `src/models/solver/CandidateGovernance.jl`
- `src/models/solver/Solver.jl`
- `src/models/solver/ProblemSpec.jl`
- `src/models/scans/ScanCommon.jl`
- `src/models/solver/SeedStrategies.jl`（或新增 seed builder 承载文件）
- `tests/unit/models/test_candidate_governance_contract.jl`
- `tests/unit/models/test_solver.jl`
- `tests/unit/models/test_problem_spec_contract.jl`
- `tests/unit/models/test_scan_common.jl`

## 5. 验证命令（建议）

```powershell
julia --project=. -e 'ENV["UNIT_FILES"]="models/test_candidate_governance_contract.jl;models/test_solver.jl;models/test_problem_spec_contract.jl;models/test_scan_common.jl"; include("tests/unit/runtests.jl")'
julia --project=. -e 'ENV["REGRESSION_FILES"]="models/test_solver_attempt_engine_convergence_regression.jl;pnjl/test_constraint_selection_regression.jl"; include("tests/regression/runtests.jl")'
```

## 6. PR-B DoD

- [ ] 治理口径单点实现，调用方不再自行拼装 success/hard-constraint 判定。
- [ ] seed pool 单点组装，调用方不再重复拼接候选。
- [ ] 关键 unit/regression 全通过，无不可解释行为漂移。
- [x] 关键 unit/regression 全通过，无不可解释行为漂移。
- [ ] 本任务单状态、验证记录与代码实际一致。
- [x] 本任务单状态、验证记录与代码实际一致。

#### 本轮验证（2026-04-08，seed 收敛后）

- unit:
  - `ENV["UNIT_FILES"]="models/test_candidate_governance_contract.jl;models/test_solver.jl;models/test_problem_spec_contract.jl;models/test_scan_common.jl;models/test_tmu_scan.jl;models/test_trho_scan.jl"`
  - 结果：`310 passed`
- regression:
  - `ENV["REGRESSION_FILES"]="models/test_solver_attempt_engine_convergence_regression.jl;pnjl/test_constraint_selection_regression.jl"`
  - 结果：`36 passed`

## 8. 能力 -> 唯一实现点矩阵（PR-B）

| 能力 | 唯一实现点 | 说明 |
|---|---|---|
| candidate success 判定 | `CandidateGovernance.evaluate_candidate_success` | 各层统一调用，避免局部口径分叉。 |
| candidate 最小 schema 归一 | `CandidateGovernance.normalize_governance_candidate` | 统一 `converged/residual_norm/hard_constraint_ok/failed_constraints/pressure/seed_index`。 |
| seed pool 组装与去重 | `CandidateGovernance.build_seed_pool` | 统一 primary -> extra/provided -> default 顺序与去重规则。 |
| ProblemSpec 治理执行 | `ProblemSpec._execute_governed_attempt_plan` | 统一调用 governance helper，不再本地判定 success。 |
| scan 治理候选聚合 | `ScanCommon.attempt_with_candidates` | 统一调用 governance helper，selection 仍由 selector 决定。 |

## 7. 风险与回退

- 风险：候选归一化可能改变边缘路径选择。
  - 缓解：保持 selector 不变，先收敛结构后收敛策略。
- 风险：seed 去重顺序导致历史路径微调。
  - 缓解：先固定 deterministic 顺序并做基线对比。
- 回退：按小步提交拆分（治理收敛 / seed 收敛），必要时按子阶段回退。
