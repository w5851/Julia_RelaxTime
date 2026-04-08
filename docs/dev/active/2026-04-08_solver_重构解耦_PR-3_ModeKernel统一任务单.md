# Solver 重构解耦 PR-3 任务单（ModeKernel 统一）

日期：2026-04-08
关联总纲：`docs/dev/active/2026-04-08_solver_重构解耦总纲领与PR路线图.md`
阶段目标：统一 `ConstraintSolverFixed*.jl` 模式求解骨架，消除重复后处理模板。

## 0. 协作与分支规则（本 PR 强制）

- 合并方式：`squash + delete`。
- 开工前：从最新 `main` 新建分支（建议：`refactor/solver-pr3-mode-kernel-unify`）。
- 主分支纪律：保持 `main` 干净，不在 `main` 直接开发。
- 工作区策略：本项目本 PR 不新建 worktree，直接在当前工作区切分支执行。
- 收尾：PR 合并后删除远端分支，并在本地删除已合并分支后回到 `main`。

## 1. 范围与非范围

### 前置条件

- 基线分支：`main` 已同步远端最新提交。
- 依赖关系：PR-2 已合并（ProblemSpec 瘦身完成）。
- 开工检查（建议执行）：
  - `git switch main && git pull --ff-only`
  - `git switch -c refactor/solver-pr3-mode-kernel-unify`

### 范围

- 收敛各 mode 文件中的共用流程：
  - residual 构造
  - 后处理热力学计算
  - residual_norm 合成
  - 空候选/失败候选降级模板
- 目标文件：
  - `ConstraintSolverFixedMu.jl`
  - `ConstraintSolverFixedRho.jl`
  - `ConstraintSolverFixedEntropy.jl`
  - `ConstraintSolverFixedSigma.jl`
  - `ConstraintSolverFixedAsymmetricRho.jl`

### 非范围

- 不变更物理约束定义本身（只重组实现）。
- 不改 public API 签名。

## 2. 任务分解（可勾选）

- [ ] 提取共用 postprocess helper（压力、密度、熵、能量、质量、有限性检查）。
- [ ] 提取共用 residual_norm 合成策略（gap + mode-constraint 最大范数）。
- [ ] 提取失败候选构造 helper（统一字段与默认值）。
- [ ] 统一 FixedRho/Entropy/Sigma/AsymmetricRho 的 outer solve 模板。
- [ ] 统一 hard constraint 接入点，避免 mode 文件内重复拼接。
- [ ] 保留 mode 特有约束表达（例如 ud_ratio、sigma 约束）在最小局部函数。
- [ ] 补充 unit 回归：覆盖各 mode 的 converged/residual_norm 字段一致性。
- [ ] 补充 integration smoke：确保多 mode 主路径可用。

## 3. 验证与验收标准

- [ ] `ConstraintSolverFixed*.jl` 重复模板减少（目标至少 30%）。
- [ ] 各 mode 结果字段语义一致（`converged/iterations/residual_norm`）。
- [ ] 既有 solver 相关 unit/integration/regression 全绿。
- [ ] 性能不出现显著退化（必要时补充 benchmark 对比）。

建议最小验证命令（按顺序）：

- `julia --project=. -e 'ENV["UNIT_PROFILE"]="full"; include("tests/unit/runtests.jl")'`
- `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
- `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`

## 4. 风险与回退

- 风险：模板统一误伤某 mode 特殊判据。
- 回退：保留 mode 局部 override 点，必要时按 mode 单独回退。

## 5. PR 交付物

- 代码变更：`ConstraintSolverFixed*.jl` + 共用 helper（可能在 Common 或新文件）。
- 测试变更：mode 间一致性与主路径回归。
- 文档变更：本任务单状态更新。

## 6. 完成信号（Completion Signal）

- [ ] 本任务单条目已全部完成并回写，再创建 PR。
- [ ] PR 创建后等待 CI 全通过。
- [ ] PR 创建后等待 GitHub Copilot 自动 review（PR 触发）并完成 review 问题处理。
- [ ] PR 已通过评审并采用 `squash + delete` 合并。
- [ ] 本地执行：`git switch main`，删除本地开发分支。
- [ ] 本任务单状态已更新（含重复度下降依据与关键模式回归结果）。

## 7. 与下一 PR 的衔接检查（PR-3 -> PR-4）

- [ ] mode kernel 输出字段与候选结构已统一，满足治理主引擎单语义输入要求。
- [ ] 无残留“模式私有排序逻辑”阻断 PR-4 合并治理规则。
- [ ] 已在 PR 描述中标注“对 PR-4 的治理收敛准备项”。
