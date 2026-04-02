# Auto/Legacy 语义收敛 PR-0 微计划

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在进入 R1 功能开发前，先冻结 FixedRho 数值验收门槛、多初值候选治理规则与 auto/models/legacy 对照矩阵，避免“结构推进但数值漂移”。

**Architecture:** 仅基于 `docs/dev/active/2026-04-02_auto-legacy求解器语义同构与模式收敛任务单.md` 做前置准备，不改动主架构代码。先产出可执行验收规范，再进入下一轮功能 PR。

**Tech Stack:** Julia 1.10+, `Models`/`PNJL` 求解入口, `tests/unit|integration|regression`, `scripts/dev/export_pnjl_constraint_fixedpoint_baseline.jl`.

---

## Scope and Constraints

- 仅做前置准备与验收规范固化，不进入 R1 代码实现。
- 不重写 04-02 主任务单，只做可执行化补充。
- 不新增与 FixedRho/语义同构无关的扩展功能。

## File Responsibility Map

- Modify: `docs/dev/active/2026-04-02_auto-legacy求解器语义同构与模式收敛任务单.md`
  - 责任：补充 PR-0 硬门槛、测试矩阵、停止线。
- Create: `tests/regression/models/test_fixedrho_precision_guard_regression.jl`
  - 责任：固定关键点精度护栏（residual/rho 偏差门槛）。
- Create: `tests/integration/models/test_solver_backend_semantic_parity_guard.jl`
  - 责任：`auto/models/legacy` 同语义对照矩阵。
- Create: `tests/unit/models/test_multiseed_candidate_governance.jl`
  - 责任：多初值 selector/tie-break 行为可重复。

## Task 0: 冻结 PR-0 验收阈值（文档先行）

**Files:**
- Modify: `docs/dev/active/2026-04-02_auto-legacy求解器语义同构与模式收敛任务单.md`

- [x] **Step 1: 补充 FixedRho 精度硬门槛**
  - 关键点（`T=90/110/130 MeV`, `rho*=0.2/0.6/1.0`）要求：
    - `residual_norm <= 1e-6`（目标：`<=1e-8`）
    - `abs(rho_norm-rho_target) <= 1e-6`

- [x] **Step 2: 补充候选治理硬规则**
  - 多初值候选排序固定：
    1) `hard_constraint_ok`
    2) `residual_norm` 最小
    3) `pressure` 最大
    4) `seed_index` 最小（保证可重复）

- [x] **Step 3: 补充语义同构对照矩阵**
  - 同点同参数下比较 `solver_backend=:auto/:models/:legacy`：
    - 收敛一致性
    - 关键热力学量容差
    - 分支标签一致性（若有）

- [x] **Step 4: 定义停止线**
  - 任一触发即禁止进入 R1：
    - FixedRho 精度门槛不达标
    - auto/models/legacy 同语义不一致
    - regression smoke 退化

## Task 1: FixedRho 精度护栏测试（先失败后实现）

**Files:**
- Create: `tests/regression/models/test_fixedrho_precision_guard_regression.jl`
- Modify: `tests/regression/runtests.jl` (接线)

- [x] **Step 1: 写失败 regression 测试**
  - 断言关键点满足 PR-0 精度门槛。

- [x] **Step 2: 运行并确认失败（若当前未满足）**
  - Run: `julia --project=. -e 'ENV["REGRESSION_FILES"]="models/test_fixedrho_precision_guard_regression.jl"; include("tests/regression/runtests.jl")'`

- [x] **Step 3: 最小实现/参数校准使其转绿**
  - 仅允许修复与门槛直接相关的求解判据与候选治理。

- [x] **Step 4: 复跑验证**
  - 同 Step 2，预期 PASS。

## Task 2: 多初值候选治理护栏（先失败后实现）

**Files:**
- Create: `tests/unit/models/test_multiseed_candidate_governance.jl`
- Modify: `src/models/constraint_solver.jl` (若需)

- [x] **Step 1: 写失败 unit 测试**
  - 覆盖 selector 与 tie-break 可重复性。

- [x] **Step 2: 运行并确认失败**
  - Run: `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_multiseed_candidate_governance.jl"; include("tests/unit/runtests.jl")'`

- [x] **Step 3: 最小实现使其通过**
  - 固化排序规则与稳定 tie-break。

- [x] **Step 4: 回归 unit smoke**
  - Run: `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`

## Task 3: auto/models/legacy 语义对照护栏

**Files:**
- Create: `tests/integration/models/test_solver_backend_semantic_parity_guard.jl`

- [x] **Step 1: 写失败 integration 对照测试**
  - 对照固定点集，比较三类 backend 在同语义下结果。

- [x] **Step 2: 运行并确认失败（若存在偏差）**
  - Run: `julia --project=. -e 'include("tests/integration/models/test_solver_backend_semantic_parity_guard.jl")'`

- [x] **Step 3: 最小修复并复跑**
  - 保证 `auto` 仅路由实现，不改变语义定义。

- [x] **Step 4: integration smoke 回归**
  - Run: `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`

## Task 4: PR-0 收口与进入 R1 准入判定

**Files:**
- Modify: `docs/dev/active/2026-04-02_auto-legacy求解器语义同构与模式收敛任务单.md`

- [x] **Step 1: 运行最终验证矩阵**
  - `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
  - `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
  - `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
  - `julia --project=. scripts/dev/check_docs_consistency.jl`
  - `julia --project=. scripts/dev/check_active_docs_governance.jl`

- [x] **Step 2: 回填 PR-0 证据与结论**
  - 明确：是否达到 R1 准入。

- [x] **Step 3: 仅在准入通过后开启下一轮功能 PR**
  - 分支命名建议：`feat/models-solver-semantic-convergence-r1`。

## Done Definition

- PR-0 三类硬门槛（精度/候选治理/语义对照）全部落地并可复现。
- 文档中明确 R1 准入结果与证据。
- 未引入与 04-02 主线无关的范围膨胀。
