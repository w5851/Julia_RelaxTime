---
title: Solver 五大痛点治理 PR-E 任务单（Entry Unification + Governance Unification）
archived: true
original: docs/dev/active/2026-04-08_solver_五大痛点治理_PR-E任务单.md
archived_date: 2026-04-08
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Solver 五大痛点治理 PR-E 任务单（Entry Unification + Governance Unification）

## 1. 目标

- [x] FixedMu 与非 FixedMu 统一走 ProblemSpec 主链，消除“特殊直通”心智分叉。
- [x] attempt 执行与 selector 选择形成单一路径，调用方不再重复筛选。

## 2. 范围

### 2.1 本期范围

- [x] E1.1 `solve_constraint(::FixedMu)` 默认改为 ProblemSpec 主链（保留短期开关以便回退）。
- [x] E1.2 统一 kwargs 解释语义（`seed_guess/seed_candidates/selector/semantic_mode/diagnostic_level`）。
- [x] E2.1 明确 selector 单点调用位置（建议只在治理尾部调用一次）。
- [x] E2.2 清理 `Solver.jl` 局部二次筛选逻辑，避免“先治理再本地再筛选”。
- [x] E2.3 统一 candidate 最小字段契约在各层传递，不再隐式补字段。

### 2.2 非范围

- [x] 不在本 PR 内拆分大文件（留给 PR-F，已在 PR-F 完成）。
- [x] 不新增新的 selector 策略类型（已满足）。

## 2.3 深层治理判据（防止保留双轨）

- [x] FixedMu 旧直通路径不长期并存；已设短期回退开关与退役期限。
- [x] selector 不在 `Solver`、`ProblemSpec`、scan 侧重复调用；已定义并测试“唯一调用点”。
- [x] kwargs 语义映射形成单一映射，避免调用层各自删 key/补 key。

## 2.4 退役计划（必须落文档）

- [x] 标注 `fixedmu_use_problem_spec` 迁移阶段：默认开启日期、旧路径移除条件、最终删除目标版本。
- [x] 标注历史兼容参数的退役窗口与 warning 文案。

### 2.4.1 `fixedmu_use_problem_spec` 退役时间线

- 默认开启日期：`2026-04-08`（PR-E 起，`solve_constraint(::FixedMu)` 默认 `fixedmu_use_problem_spec=true`）。
- 兼容窗口：`2026-04-08` ~ `2026-05-31`，仅保留显式 `fixedmu_use_problem_spec=false` 作为紧急回退。
- 旧路径移除条件（全部满足后执行）：
  1. 连续两周 nightly/regression 无 FixedMu 回归差异告警；
  2. `tests/unit/models/test_solver.jl` 与 `tests/unit/models/test_solver_diagnostic_contract.jl` 中默认主链契约稳定；
  3. 无活跃 issue 依赖旧直通路径行为。
- 目标删除版本：`v0.9`（移除 `fixedmu_use_problem_spec` 参数与旧直通分支实现）。

### 2.4.2 兼容参数 warning 策略

- 当前阶段：保留参数但在文档标注“Deprecated compatibility switch (temporary)”。
- 过渡阶段（`2026-05`）：当显式传 `fixedmu_use_problem_spec=false` 时增加 `@warn`，提示迁移到 ProblemSpec 默认主链。
- 删除阶段（`v0.9`）：参数被移除，继续传入时抛 `ArgumentError`（与 `use_problem_spec/allow_legacy_path` 的现行治理一致）。

## 3. 实施任务（可勾选）

### 3.1 入口统一

- [x] 给 `FixedMu` 增加 ProblemSpec 主链回归测试（功能等价 + 诊断字段等价）。
- [x] 将 `fixedmu_use_problem_spec` 从“可选路径”升级为“默认路径”（保留兼容退路开关，默认关闭旧路）。
- [x] 清理与旧直通路径绑定的特殊 kwargs 过滤逻辑。

### 3.2 治理统一

- [x] 审计 `solve_multi` / `ProblemSpec` / scan 侧 selector 调用次数，确保单次选择。
- [x] 删除重复的 `_select_*` 局部 helper 或将其降为 selector 适配器。
- [x] 在诊断输出中固化 `selection_reason` 来源字段，保证可追溯。

## 4. 验证清单

- [x] 单测：`tests/unit/models/test_solver.jl` 增加 FixedMu ProblemSpec 默认链路断言。
- [x] 单测：`tests/unit/models/test_problem_spec_contract.jl` 增加 selector 单次调用断言（可通过 mock selector 计数）。
- [x] 单测：增加“FixedMu 主链默认生效”断言，防止回滚为双轨默认。
- [x] 单测：增加“selector 调用次数 = 1”断言（覆盖 solve_multi 与 ProblemSpec 路径）。
- [x] 回归：
  - `models/test_solver_attempt_engine_convergence_regression.jl`
  - `models/test_solver_diagnostic_exception_regression.jl`
  - `pnjl/test_constraint_selection_regression.jl`

## 5. 风险与缓解

- 风险：FixedMu 切主链后，历史 kwargs 某些边缘行为变化。
  - 缓解：建立“兼容参数映射表”并加迁移提示 warning（短窗口）。
- 风险：selector 单点化后改变个别分支 tie-break 行为。
  - 缓解：显式固化 tie-break（residual -> pressure -> seed_index）并回归守护。

## 6. PR-E DoD

- [x] FixedMu 与非 FixedMu 在编排层走统一 ProblemSpec 主路径。
- [x] selector 只在治理尾部单次调用，去除二次筛选。
- [x] 关键 unit/regression 通过，诊断字段来源一致可追踪。
- [x] FixedMu 双轨语义结束（或有明确、短期、可测的回退窗口）。

## 7. 进展记录（2026-04-08）

- [x] `solve_constraint(::FixedMu)` 默认切换为 ProblemSpec 主链；`fixedmu_use_problem_spec` 仍保留为短期回退开关。
- [x] `diagnostic_level` 与 FixedMu 旧直通链路冲突时给出显式错误，避免语义分叉。
- [x] `solve(::non-FixedMu)` 路径改为把 seed 池一次性下发给 `solve_constraint(..., seed_candidates=...)`，移除 `Solver.jl` 内局部二次筛选。
- [x] selector 单次调用契约已加测试（ProblemSpec 直调 + Solver 包装调用）。
- [x] E2.3 candidate 最小字段契约收敛：`normalize_governance_candidate` 改为严格字段契约校验（缺字段/一致性冲突即报错），调用层补齐最小字段后再归一化。
- [x] 删除/降级 `Solver.jl` 局部 `_select_*` 逻辑：保留 `_fixedmu_multiseed_selector_adapter` 作为本地语义适配器，不再承担非统一治理职责。
- [x] 诊断输出新增 `selection_reason_source=:problem_spec_selector`（summary/full/candidates），满足来源可追溯。
- [x] 退役计划文档化已补充（默认开启日期/移除条件/目标版本 + warning 策略）。

### 本轮验证

- [x] `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_solver.jl;models/test_problem_spec_contract.jl;models/test_solver_diagnostic_contract.jl"; include("tests/unit/runtests.jl")'`
- [x] `julia --project=. -e 'include("tests/integration/models/test_solver_auto_backend_semantic_parity.jl"); include("tests/integration/models/test_phase_solver_diagnostic_adapter_smoke.jl")'`
- [x] `julia --project=. -e 'ENV["REGRESSION_FILES"]="models/test_solver_attempt_engine_convergence_regression.jl;models/test_solver_diagnostic_exception_regression.jl;pnjl/test_constraint_selection_regression.jl"; include("tests/regression/runtests.jl")'`
