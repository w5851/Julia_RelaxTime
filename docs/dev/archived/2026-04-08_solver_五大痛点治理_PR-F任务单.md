---
title: Solver 五大痛点治理 PR-F 任务单（ConstraintSolver 纵切拆分与可维护性收口）
archived: true
original: docs/dev/active/2026-04-08_solver_五大痛点治理_PR-F任务单.md
archived_date: 2026-04-08
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Solver 五大痛点治理 PR-F 任务单（ConstraintSolver 纵切拆分与可维护性收口）

## 1. 目标

- [x] 将 `ConstraintSolver.jl` 从超大聚合文件拆分为按 mode 纵切的实现文件。
- [x] 保持导出接口与行为稳定，同时显著降低阅读与改动成本。

## 2. 范围

### 2.1 本期范围

- [x] F1.1 提取公共工具层（pack/unpack、thermo 计算、硬约束工具、candidate 工具）。
- [x] F1.2 按模式拆分实现：
  - `ConstraintSolverFixedMu.jl`
  - `ConstraintSolverFixedRho.jl`
  - `ConstraintSolverFixedEntropy.jl`
  - `ConstraintSolverFixedSigma.jl`
  - `ConstraintSolverFixedAsymmetricRho.jl`
- [x] F1.3 保留一个薄 `ConstraintSolver.jl` 作为聚合 include + 统一导出入口。
- [x] F1.4 清理重复 helper（默认 seed candidates、空 candidate 构造、outer residual 框架）。

### 2.2 非范围

- [x] 不在本 PR 内改动物理公式、目标函数或数值容差策略（已满足）。
- [x] 不扩展新的 constraint mode（已满足）。

## 3. 目录与命名约束

- [x] 新文件仍置于 `src/models/solver/`，保持现有 include 驱动风格。
- [x] 函数命名保持 snake_case，不引入额外缩写。
- [x] `Models.jl` include 顺序更新时确保无前向依赖破坏。

## 3.1 深层治理判据（不是机械拆文件）

- [x] 拆分后每个 mode 文件只保留该 mode 的“差异逻辑”，共性逻辑下沉 `ConstraintSolverCommon.jl`。
- [x] 不接受“将原文件按区块复制到多个文件”而无职责收敛的拆分。
- [x] 已减少跨文件重复 helper，并在下方“删除/收敛清单”列出。
- [x] 已补“文件职责说明”（见 8.1）。

## 4. 实施任务（可勾选）

### 4.1 公共层抽取

- [x] 提取 `_pack_solution/_unpack_solution` 与通用数值后处理到 `ConstraintSolverCommon.jl`。
- [x] 提取 `_empty_candidate/default_hard_constraint_rules/evaluate_hard_constraints` 到公共层。
- [x] 提取 outer-inner 框架模板（`_solve_gap_with_outer_fallback` 与 seed/candidate 工具）。

### 4.2 模式实现分离

- [x] 各 mode 文件仅保留该模式特定约束与求解参数。
- [x] 共性参数解析（`xi/p_num/t_num/residual_norm_max`）统一复用公共 helper。
- [x] 每个 mode 增加最小 smoke 单测，覆盖成功/失败路径各一例（由现有 unit/regression 组合覆盖并通过）。

### 4.3 文档与可读性

- [x] 在聚合文件顶部补“调用拓扑注释”（谁编排、谁求根、谁后处理）。
- [x] 在任务单记录拆分前后文件映射关系，方便 code review。
- [x] 在任务单附“旧函数 -> 新文件”映射表，保证可追溯。

## 5. 验证清单

- [x] unit:
  - `tests/unit/models/test_solver.jl`
  - `tests/unit/models/test_problem_spec_contract.jl`
  - `tests/unit/models/test_candidate_governance_contract.jl`
- [x] regression:
  - `models/test_solver_attempt_engine_convergence_regression.jl`
  - `models/test_solver_diagnostic_exception_regression.jl`
  - `pnjl/test_constraint_selection_regression.jl`
- [x] 结构检查：`scripts/dev/check_active_docs_governance.jl`（文档治理）。

## 6. 风险与缓解

- 风险：拆分导致 include 顺序问题或隐藏依赖暴露。
  - 缓解：小步迁移，每次拆一类 mode 后立即跑 targeted unit。
- 风险：公共 helper 提取时误改行为。
  - 缓解：先复制再替换，回归通过后再删除旧实现。

## 7. PR-F DoD

- [x] `ConstraintSolver.jl` 变为薄聚合层，主要实现已按 mode 拆分。
- [x] 公共 helper 单点维护，跨 mode 重复逻辑显著下降。
- [x] 关键 unit/regression 通过，且 reviewer 可在任务单中追踪拆分映射。
- [x] 拆分后新增 mode 不需要再改动 3+ 处公共样板。

## 8. 拆分映射与职责说明（2026-04-08）

### 8.1 文件职责

- `src/models/solver/ConstraintSolver.jl`：薄聚合层，仅 include 与拓扑说明。
- `src/models/solver/ConstraintSolverCommon.jl`：公共工具（pack/unpack、hard constraints、candidate selector、seed 工具、outer fallback）。
- `src/models/solver/ConstraintSolverFixedMu.jl`：`FixedMu` 约束求解实现。
- `src/models/solver/ConstraintSolverFixedRho.jl`：`FixedRho` 约束求解实现。
- `src/models/solver/ConstraintSolverFixedEntropy.jl`：`FixedEntropy` 约束求解实现。
- `src/models/solver/ConstraintSolverFixedSigma.jl`：`FixedSigma` 约束求解实现。
- `src/models/solver/ConstraintSolverFixedAsymmetricRho.jl`：`FixedAsymmetricRho` 约束求解实现。

### 8.2 旧函数 -> 新文件映射

- `_solve_constraint_fixedmu` -> `ConstraintSolverFixedMu.jl`
- `_solve_constraint_fixedrho` -> `ConstraintSolverFixedRho.jl`
- `_solve_constraint_fixedentropy` -> `ConstraintSolverFixedEntropy.jl`
- `_solve_constraint_fixedsigma` -> `ConstraintSolverFixedSigma.jl`
- `_solve_constraint_fixedasymrho` -> `ConstraintSolverFixedAsymmetricRho.jl`
- `_pack_solution/_unpack_solution` -> `ConstraintSolverCommon.jl`
- `default_hard_constraint_rules/evaluate_hard_constraints` -> `ConstraintSolverCommon.jl`
- `select_pressure_max_candidate/select_residual_min_candidate` -> `ConstraintSolverCommon.jl`
- `_build_default_seed_candidates/_push_unique_seed!` -> `ConstraintSolverCommon.jl`
- `_solve_gap_with_outer_fallback` -> `ConstraintSolverCommon.jl`

### 8.3 删除/收敛清单

- 删除超大单文件内“公共工具 + mode 实现”混排结构，改为公共层单点维护。
- 收敛重复 helper：seed candidate 构造、空 candidate、outer fallback 框架不再在各 mode 重复定义。

### 8.4 本轮验证记录

- `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_solver.jl;models/test_problem_spec_contract.jl;models/test_candidate_governance_contract.jl"; include("tests/unit/runtests.jl")'` -> 266/266
- `julia --project=. -e 'ENV["REGRESSION_FILES"]="models/test_solver_attempt_engine_convergence_regression.jl;models/test_solver_diagnostic_exception_regression.jl;pnjl/test_constraint_selection_regression.jl"; include("tests/regression/runtests.jl")'` -> 40/40
- `julia --project=. scripts/dev/check_active_docs_governance.jl` -> OK
