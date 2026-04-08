# Solver 五大痛点治理 PR-F 任务单（ConstraintSolver 纵切拆分与可维护性收口）

## 1. 目标

- [ ] 将 `ConstraintSolver.jl` 从超大聚合文件拆分为按 mode 纵切的实现文件。
- [ ] 保持导出接口与行为稳定，同时显著降低阅读与改动成本。

## 2. 范围

### 2.1 本期范围

- [ ] F1.1 提取公共工具层（pack/unpack、thermo 计算、硬约束工具、candidate 工具）。
- [ ] F1.2 按模式拆分实现：
  - `ConstraintSolverFixedMu.jl`
  - `ConstraintSolverFixedRho.jl`
  - `ConstraintSolverFixedEntropy.jl`
  - `ConstraintSolverFixedSigma.jl`
  - `ConstraintSolverFixedAsymmetricRho.jl`
- [ ] F1.3 保留一个薄 `ConstraintSolver.jl` 作为聚合 include + 统一导出入口。
- [ ] F1.4 清理重复 helper（默认 seed candidates、空 candidate 构造、outer residual 框架）。

### 2.2 非范围

- [ ] 不在本 PR 内改动物理公式、目标函数或数值容差策略。
- [ ] 不扩展新的 constraint mode。

## 3. 目录与命名约束

- [ ] 新文件仍置于 `src/models/solver/`，保持现有 include 驱动风格。
- [ ] 函数命名保持 snake_case，不引入额外缩写。
- [ ] `Models.jl` include 顺序更新时确保无前向依赖破坏。

## 3.1 深层治理判据（不是机械拆文件）

- [ ] 拆分后每个 mode 文件只保留该 mode 的“差异逻辑”，共性逻辑必须下沉公共层。
- [ ] 不接受“将原文件按区块复制到多个文件”而无职责收敛的拆分。
- [ ] 必须减少跨文件重复 helper 的数量，并以清单给出删除项。
- [ ] 必须补“文件职责说明”，让 reviewer 可快速定位改动影响域。

## 4. 实施任务（可勾选）

### 4.1 公共层抽取

- [ ] 提取 `_pack_solution/_unpack_solution` 与通用数值后处理到公共 helper 区。
- [ ] 提取 `_empty_candidate/default_hard_constraint_rules/evaluate_hard_constraints`（如仍在本文件维护）。
- [ ] 提取 outer-inner 框架模板，避免各 mode 复制粘贴。

### 4.2 模式实现分离

- [ ] 各 mode 文件仅保留该模式特定约束与求解参数。
- [ ] 共性参数解析（`xi/p_num/t_num/residual_norm_max`）统一由上层 helper 处理。
- [ ] 每个 mode 增加最小 smoke 单测，覆盖成功/失败路径各一例。

### 4.3 文档与可读性

- [ ] 在聚合文件顶部补“调用拓扑注释”（谁编排、谁求根、谁后处理）。
- [ ] 在任务单记录拆分前后文件映射关系，方便 code review。
- [ ] 在任务单附“旧函数 -> 新文件”映射表，保证可追溯。

## 5. 验证清单

- [ ] unit:
  - `tests/unit/models/test_solver.jl`
  - `tests/unit/models/test_problem_spec_contract.jl`
  - `tests/unit/models/test_candidate_governance_contract.jl`
- [ ] regression:
  - `models/test_solver_attempt_engine_convergence_regression.jl`
  - `models/test_solver_diagnostic_exception_regression.jl`
  - `pnjl/test_constraint_selection_regression.jl`
- [ ] 结构检查：`scripts/dev/check_active_docs_governance.jl`（文档治理）、必要时 `scripts/dev/check_docs_consistency.jl`。

## 6. 风险与缓解

- 风险：拆分导致 include 顺序问题或隐藏依赖暴露。
  - 缓解：小步迁移，每次拆一类 mode 后立即跑 targeted unit。
- 风险：公共 helper 提取时误改行为。
  - 缓解：先复制再替换，回归通过后再删除旧实现。

## 7. PR-F DoD

- [ ] `ConstraintSolver.jl` 变为薄聚合层，主要实现已按 mode 拆分。
- [ ] 公共 helper 单点维护，跨 mode 重复逻辑显著下降。
- [ ] 关键 unit/regression 通过，且 reviewer 可在任务单中追踪拆分映射。
- [ ] 拆分后新增 mode 不需要再改动 3+ 处公共样板。
