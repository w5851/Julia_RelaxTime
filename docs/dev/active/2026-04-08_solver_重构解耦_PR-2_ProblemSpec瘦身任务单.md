# Solver 重构解耦 PR-2 任务单（ProblemSpec 瘦身）

日期：2026-04-08
关联总纲：`docs/dev/active/2026-04-08_solver_重构解耦总纲领与PR路线图.md`
阶段目标：拆分 `ProblemSpec.jl` 的编排与诊断职责，保留契约核心。

## 0. 协作与分支规则（本 PR 强制）

- 合并方式：`squash + delete`。
- 开工前：从最新 `main` 新建分支（建议：`refactor/solver-pr2-problemspec-slim`）。
- 主分支纪律：保持 `main` 干净，不在 `main` 直接开发。
- 工作区策略：本项目本 PR 不新建 worktree，直接在当前工作区切分支执行。
- 收尾：PR 合并后删除远端分支，并在本地删除已合并分支后回到 `main`。

## 1. 范围与非范围

### 前置条件

- 基线分支：`main` 已同步远端最新提交。
- 依赖关系：PR-1 已合并（边界硬化完成）。
- 开工检查（建议执行）：
  - `git switch main && git pull --ff-only`
  - `git switch -c refactor/solver-pr2-problemspec-slim`

### 范围

- 从 `ProblemSpec.jl` 拆出：
  - attempt plan 生成
  - governed attempt 执行编排
  - diagnostic 构建/附加
- 建议新增文件：
  - `src/models/solver/ProblemSpecOrchestrator.jl`
  - `src/models/solver/SolverDiagnostics.jl`

### 非范围

- 不改 selector 的物理语义。
- 不改 mode 级 residual 定义。
- 不在本 PR 合并治理引擎。

## 2. 任务分解（可勾选）

- [x] 标记 `ProblemSpec.jl` 中契约定义与执行编排边界（注释或 section）。
- [x] 迁移 `_execute_governed_attempt_plan` 到 orchestrator 文件。
- [x] 迁移 `_solver_diagnostic_from_candidate` 与 `_attach_solver_diagnostic` 到 diagnostics 文件。
- [x] 迁移 `_build_governed_attempt_plan` 到 orchestrator 文件并保留调用兼容。
- [x] 将 `ProblemSpec` 保持为“mode 契约 + dim + forward_solve 注册”主职责。
- [x] 更新 `src/models/Models.jl` include 顺序，确保加载依赖正确。
- [x] 补 unit/regression：覆盖 `diagnostic_level=:none/:summary/:full` 输出结构。
- [x] 更新 solver API 文档中“维护者阅读顺序”。

## 3. 验证与验收标准

- [x] `ProblemSpec.jl` 体量明显下降（目标区间约 300-450 行）。
- [x] 迁移前后诊断字段完全兼容。
- [x] 关键 integration 与 regression 测试全绿。
- [x] 无循环 include 或初始化顺序错误。

建议最小验证命令（按顺序）：

- `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_solver_diagnostic_contract.jl,models/test_solver.jl"; include("tests/unit/runtests.jl")'`
- `julia --project=. -e 'ENV["REGRESSION_FILES"]="models/test_solver_diagnostic_exception_regression.jl,models/test_solver_attempt_engine_convergence_regression.jl"; include("tests/regression/runtests.jl")'`
- `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`

## 4. 风险与回退

- 风险：函数迁移后作用域/依赖遗漏。
- 回退：先保留薄包装转调，再逐步内联删除旧实现。

## 5. PR 交付物

- 代码变更：`ProblemSpec.jl` + 新增 orchestrator/diagnostics 文件 + include 调整。
- 测试变更：diagnostic 合同回归测试。
- 文档变更：本任务单勾选状态与阅读路径说明。

## 6. 完成信号（Completion Signal）

- [x] 本任务单条目已全部完成并回写，再创建 PR。
- [x] PR 创建后等待 CI 全通过。
- [x] PR 创建后等待 GitHub Copilot 自动 review（PR 触发）并完成 review 问题处理。
- [x] PR 已通过评审并采用 `squash + delete` 合并。
- [x] 本地执行：`git switch main`，删除本地开发分支。
- [x] 本任务单状态已更新（含 `ProblemSpec.jl` 行数变化与诊断兼容证据）。

## 8. 实际结果摘要（PR-2 首轮）

- 结构拆分：
  - 新增 `src/models/solver/ProblemSpecOrchestrator.jl`（attempt plan / governed attempt 编排 / mode forward_solve 组织）
  - 新增 `src/models/solver/SolverDiagnostics.jl`（diagnostic level 解析与诊断拼装）
  - `src/models/solver/ProblemSpec.jl` 收敛为契约核心与 forward_solve 注册绑定。
- include 顺序：`src/models/Models.jl` 已在 `ProblemSpec.jl` 后、`Solver.jl` 前加入上述两个新文件。
- 文档：`docs/api/models/solver/README.md` 已补充维护者阅读顺序（ProblemSpec -> Orchestrator -> Diagnostics -> Solver）。

### 验证证据

- `ProblemSpec.jl` 行数：`113`（命令：`julia --project=. -e 'include("src/models/Models.jl"); println(length(readlines("src/models/solver/ProblemSpec.jl")))'`）
  - 说明：低于任务单中的“约 300-450 行”参考区间，按职责边界结果属于“超额瘦身”；当前更符合 PR-2 目标“ProblemSpec 仅保留契约核心，执行与诊断外置”。
- unit：
  - 命令：`julia --project=. -e 'ENV["UNIT_FILES"]="models/test_solver_diagnostic_contract.jl,models/test_solver.jl"; include("tests/unit/runtests.jl")'`
  - 结果：`Unit | Pass 76 / Total 76`
- regression：
  - 命令：`julia --project=. -e 'ENV["REGRESSION_FILES"]="models/test_solver_diagnostic_exception_regression.jl,models/test_solver_attempt_engine_convergence_regression.jl"; include("tests/regression/runtests.jl")'`
  - 结果：`Regression | Pass 37 / Total 37`
- integration smoke：
  - 命令：`julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
  - 结果：`Integration | Pass 463 / Total 463`
  - 备注：日志中出现既有 `fail-on-fallback` 提示，但测试汇总仍为全绿（仓库现存已知现象，非本次改动引入）。

## 7. 与下一 PR 的衔接检查（PR-2 -> PR-3）

- [ ] `ProblemSpec` 与 orchestrator/diagnostics 职责边界已稳定。
- [ ] mode kernel 改造所需调用点已统一，不再依赖旧内联实现细节。
- [ ] 已在 PR 描述中标注“对 PR-3 的输入契约与兼容约束”。

## 9. PR/CI/Review 实际闭环记录

- PR：`https://github.com/w5851/Julia_RelaxTime/pull/69`
- 合并方式：`squash + delete`（已合并到 `main`）
- Copilot review：3 条 inline 评论，均已修复并在线程回复：
  - FixedMu 诊断拼装复用 `SolverDiagnostics` helper，避免字段规则漂移。
  - `seed_candidates=nothing` 在 FixedRho 与 non-rho 路径均改为安全回退到空池/默认池。
- CI：全部检查通过（包括 `unit-smoke`、`integration-smoke`、`validation`、`pnjl-benchmark`）。
- 本地收尾：已回到 `main`，与 `origin/main` 同步；远端分支已删除（`git fetch --prune` 验证）。
