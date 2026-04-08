# Solver 五大痛点治理 PR-E 任务单（Entry Unification + Governance Unification）

## 1. 目标

- [ ] FixedMu 与非 FixedMu 统一走 ProblemSpec 主链，消除“特殊直通”心智分叉。
- [ ] attempt 执行与 selector 选择形成单一路径，调用方不再重复筛选。

## 2. 范围

### 2.1 本期范围

- [ ] E1.1 `solve_constraint(::FixedMu)` 默认改为 ProblemSpec 主链（保留短期开关以便回退）。
- [ ] E1.2 统一 kwargs 解释语义（`seed_guess/seed_candidates/selector/semantic_mode/diagnostic_level`）。
- [ ] E2.1 明确 selector 单点调用位置（建议只在治理尾部调用一次）。
- [ ] E2.2 清理 `Solver.jl` 局部二次筛选逻辑，避免“先治理再本地再筛选”。
- [ ] E2.3 统一 candidate 最小字段契约在各层传递，不再隐式补字段。

### 2.2 非范围

- [ ] 不在本 PR 内拆分大文件（留给 PR-F）。
- [ ] 不新增新的 selector 策略类型。

## 2.3 深层治理判据（防止保留双轨）

- [ ] FixedMu 旧直通路径不能长期与 ProblemSpec 主链并存；最多保留短期回退开关并设置退役期限。
- [ ] selector 不得在 `Solver`、`ProblemSpec`、scan 侧重复调用；必须定义且验证“唯一调用点”。
- [ ] kwargs 语义映射需形成单一映射表，避免调用层各自删 key/补 key。

## 2.4 退役计划（必须落文档）

- [ ] 标注 `fixedmu_use_problem_spec` 迁移阶段：默认开启日期、旧路径移除条件、最终删除目标版本。
- [ ] 标注历史兼容参数的退役窗口与 warning 文案。

## 3. 实施任务（可勾选）

### 3.1 入口统一

- [ ] 给 `FixedMu` 增加 ProblemSpec 主链回归测试（功能等价 + 诊断字段等价）。
- [ ] 将 `fixedmu_use_problem_spec` 从“可选路径”升级为“默认路径”（保留兼容退路开关，默认关闭旧路）。
- [ ] 清理与旧直通路径绑定的特殊 kwargs 过滤逻辑。

### 3.2 治理统一

- [ ] 审计 `solve_multi` / `ProblemSpec` / scan 侧 selector 调用次数，确保单次选择。
- [ ] 删除重复的 `_select_*` 局部 helper 或将其降为 selector 适配器。
- [ ] 在诊断输出中固化 `selection_reason` 来源字段，保证可追溯。

## 4. 验证清单

- [ ] 单测：`tests/unit/models/test_solver.jl` 增加 FixedMu ProblemSpec 默认链路断言。
- [ ] 单测：`tests/unit/models/test_problem_spec_contract.jl` 增加 selector 单次调用断言（可通过 mock selector 计数）。
- [ ] 单测：增加“FixedMu 主链默认生效”断言，防止回滚为双轨默认。
- [ ] 单测：增加“selector 调用次数 = 1”断言（覆盖 solve_multi 与 ProblemSpec 路径）。
- [ ] 回归：
  - `models/test_solver_attempt_engine_convergence_regression.jl`
  - `models/test_solver_diagnostic_exception_regression.jl`
  - `pnjl/test_constraint_selection_regression.jl`

## 5. 风险与缓解

- 风险：FixedMu 切主链后，历史 kwargs 某些边缘行为变化。
  - 缓解：建立“兼容参数映射表”并加迁移提示 warning（短窗口）。
- 风险：selector 单点化后改变个别分支 tie-break 行为。
  - 缓解：显式固化 tie-break（residual -> pressure -> seed_index）并回归守护。

## 6. PR-E DoD

- [ ] FixedMu 与非 FixedMu 在编排层走统一 ProblemSpec 主路径。
- [ ] selector 只在治理尾部单次调用，去除二次筛选。
- [ ] 关键 unit/regression 通过，诊断字段来源一致可追踪。
- [ ] FixedMu 双轨语义结束（或有明确、短期、可测的回退窗口）。
