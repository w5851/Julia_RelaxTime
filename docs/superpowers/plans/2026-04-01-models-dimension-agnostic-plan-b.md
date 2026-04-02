# Models Dimension-Agnostic (Plan B) Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在不以性能为首要约束的前提下，将 `models` 求解/扫描主流程彻底改造为统一动态向量管线（Plan B），并在最终收口阶段移除旧架构与兼容层，形成单一新架构主线。

**Architecture:** 上层统一采用“状态 schema + 动态 flatten/unflatten + 约束组装器”模式，消除 `5/3/8` 写死切片。中层以 `ModelStateSchema` 驱动 residual 组装与输入输出校验。迁移中允许短期兼容桥接，但该桥接仅用于过渡，最终状态必须删除旧固定维度路径与兼容层。采用分波次迁移（先契约与适配层、再 solver、再 scan/workflow、最后回归文档与旧路径清除）。

**Tech Stack:** Julia 1.10+, existing `Models` entrypoints, `NLsolve` path, `ForwardDiff`, existing test layers (`tests/unit|integration|regression`).

---

## Scope and Non-Goals

- 范围内：`src/models/solver/`, `src/models/constraint_solver.jl`, `src/models/scans/`, `src/models/entrypoints.jl` 的维度解耦与契约统一。
- 范围内：新增 `ModelStateSchema` 与动态 residual 管线；修订对外结果对象的维度敏感字段。
- 范围内：覆盖当前已注册全部模型族：`NJL/NJL2/PNJL/PNJLMagnetic/RPNJL/Rotation/GasLiquid`。
- 暂不做：数值算法替换（如 trust-region/newton 策略重构）、大规模性能优化、物理模型公式变更。

## Preflight Checklist (Before Chunk 1)

- [ ] 冻结固定维度耦合点基线清单（至少覆盖 `state.jl / solver/Conditions.jl / solver/ImplicitSolver.jl / constraint_solver.jl / scans/*`）。
- [ ] 记录 7 模型最小行为快照（可复现命令 + 输出摘要），作为迁移前对照。
- [ ] 明确外部契约变更边界（`Models` 导出、扫描输出字段、API 文档更新范围）。
- [ ] 冻结回归阈值策略（残差与关键量容忍区间），避免迁移中反复调整标准。
- [ ] 确认执行在隔离工作树/分支进行，避免与并行开发流互相污染。
- [ ] 定义停止线：出现核心模型回归失败或治理检查失败时，不进入下一 Chunk。

## File Responsibility Map

- Create: `src/models/solver/StateSchema.jl`  
  责任：定义 `ModelStateSchema`、canonical fields、flatten/unflatten、维度派生 API。
- Create: `tests/unit/models/test_state_schema.jl`  
  责任：schema 构造、序列化回环、错误输入校验。
- Modify: `src/models/solver/ImplicitSolver.jl`  
  责任：`SolverResult` 从固定维度字段迁移到 schema-aware 动态表示。
- Modify: `src/models/solver/Conditions.jl`  
  责任：`gap_conditions/build_conditions` 由固定切片迁移为 schema 驱动 residual 组装。
- Modify: `src/models/constraint_solver.jl`  
  责任：候选解评估/回退路径不再写死 `SVector{5}`、`SVector{3}`。
- Modify: `src/models/scans/TmuScan.jl` and `src/models/scans/TrhoScan.jl`  
  责任：scan IO 统一使用 schema-aware result 读取。
- Modify: `src/models/entrypoints.jl`  
  责任：统一入口返回契约保持兼容（对外字段稳定，内部来源改为 schema）。
- Create: `tests/integration/models/test_dimension_agnostic_scan_smoke.jl`  
  责任：跨模型扫描 smoke，验证无固定维度假设。
- Create: `tests/regression/models/test_dimension_agnostic_solver_regression.jl`  
  责任：数值行为回归（有限基准点）。
- Modify: `docs/api/models/solver/Overview.md` and `docs/api/models/solver/ConstraintModes.md`  
  责任：更新 solver 数据契约与维度无关说明。
- Delete/Modify: legacy fixed-dimension helper paths in `src/models/solver/` and `src/models/constraint_solver.jl`  
  责任：在最终收口阶段移除旧架构代码路径与临时 compat 适配层。

## Chunk 1: Schema Groundwork

### Task 1: 建立统一状态 schema（TDD）

**Files:**
- Create: `tests/unit/models/test_state_schema.jl`
- Create: `src/models/solver/StateSchema.jl`
- Modify: `src/models/Models.jl` (or solver module include/export wiring)

- [ ] **Step 1: 写失败测试（schema 基本契约）**
  - 覆盖：
    - `schema_for_model(:PNJL)` 返回 field 列表与顺序；
    - `flatten_state(schema, named_state)` 与 `unflatten_state(schema, vec)` 回环一致；
    - 非法长度输入抛 `ArgumentError`。

- [ ] **Step 2: 运行单测确认失败**
  - Run: `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_state_schema.jl"; include("tests/unit/runtests.jl")'`
  - Expected: FAIL（StateSchema API 未实现）

- [ ] **Step 3: 实现最小 StateSchema 能力**
  - 实现最小 API：
    - `ModelStateSchema`
    - `schema_for_model(model_kind::Symbol)`
    - `flatten_state(schema, st)`
    - `unflatten_state(schema, x::AbstractVector)`
    - `state_dim(schema)`

- [ ] **Step 4: 复跑单测转绿**
  - Run: 同 Step 2
  - Expected: PASS

- [ ] **Step 5: 提交（小步）**
  - Message style suggestion: `refactor(models/solver): introduce schema-driven state flattening contract`

## Chunk 2: Solver Core Migration

### Task 2: `SolverResult` 与条件构建改为 schema-aware

**Files:**
- Modify: `src/models/solver/ImplicitSolver.jl`
- Modify: `src/models/solver/Conditions.jl`
- Modify: `tests/unit/models/test_constraint_solver.jl`
- Create: `tests/unit/models/test_solver_dimension_agnostic.jl`

- [ ] **Step 1: 写失败测试（结果对象与 residual 组装）**
  - 覆盖：
    - `SolverResult` 不要求固定 `x_state[1:5]` / `mu[1:3]` 切片；
    - `build_conditions` 能基于 schema 从 `x` 中提取 state + constraints。

- [ ] **Step 2: 运行目标单测确认失败**
  - Run: `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_solver_dimension_agnostic.jl,models/test_constraint_solver.jl"; include("tests/unit/runtests.jl")'`

- [ ] **Step 3: 实现最小通过版本**
  - 在 `ImplicitSolver` 引入 schema 输入；
  - 在 `Conditions` 用 schema 映射替代 `x[1:5]`/`x[6:8]`；
  - 保留临时兼容 accessor（仅过渡层使用）。

- [ ] **Step 4: 回归相关 unit 套件**
  - Run: `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`

- [ ] **Step 5: 提交（小步）**
  - Message style suggestion: `refactor(models/solver): remove fixed slice assumptions in residual pipeline`

## Chunk 3: Constraint + Scan/Workflow Migration

### Task 3: 约束求解与扫描路径去除固定维度常量

**Files:**
- Modify: `src/models/constraint_solver.jl`
- Modify: `src/models/scans/TmuScan.jl`
- Modify: `src/models/scans/TrhoScan.jl`
- Modify: `src/models/entrypoints.jl`
- Create: `tests/integration/models/test_dimension_agnostic_scan_smoke.jl`

- [ ] **Step 1: 写失败 integration 测试（跨模型）**
  - 覆盖 7 模型至少一条主链路（solver/scan/workflow 之一），确保流程不依赖固定 5/3。

- [ ] **Step 2: 运行测试确认失败**
  - Run: `julia --project=. -e 'include("tests/integration/models/test_dimension_agnostic_scan_smoke.jl")'`

- [ ] **Step 3: 迁移实现**
  - `constraint_solver` 中候选评估、residual 计算改 schema 驱动；
  - `TmuScan/TrhoScan` 读取结果字段时不假设固定切片；
  - `entrypoints` 保持对外字段稳定（必要时提供向后兼容映射）。

- [ ] **Step 4: 运行 integration smoke**
  - Run: `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`

- [ ] **Step 5: 提交（小步）**
  - Message style suggestion: `refactor(models/scans): migrate scan pipeline to schema-driven solver outputs`

## Chunk 4: Regression + Docs + Governance

### Task 4: 回归收口与文档同步

**Files:**
- Create: `tests/regression/models/test_dimension_agnostic_solver_regression.jl`
- Modify: `docs/api/models/solver/Overview.md`
- Modify: `docs/api/models/solver/ConstraintModes.md`
- Modify: `docs/dev/active/2026-04-01_PNJL求解器解耦Wave-E任务单.md` (或后续新任务单)

- [ ] **Step 1: 写回归测试（有限固定点）**
  - 选择少量代表点（每模型 1-2 点）验证收敛标记、关键热力学量范围与残差阈值。
  - 覆盖模型清单：`NJL/NJL2/PNJL/PNJLMagnetic/RPNJL/Rotation/GasLiquid`。

- [ ] **Step 2: 运行 regression smoke**
  - Run: `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`

- [ ] **Step 3: 更新 API 文档与任务单证据**
  - 明确：流程维度无关，但数值核仍隐式存在 `n`（自动派生，不手写常量）。

- [ ] **Step 4: 运行治理检查**
  - Run:
    - `julia --project=. scripts/dev/check_docs_consistency.jl`
    - `julia --project=. scripts/dev/check_active_docs_governance.jl`

- [ ] **Step 5: 提交（小步）**
  - Message style suggestion: `docs(models/solver): document schema-driven dimension-agnostic contracts`

## Chunk 5: Legacy/Compat Removal Finalization

### Task 5: 删除旧架构与兼容层，完成单主线收口

**Files:**
- Modify/Delete: fixed-dimension compatibility code in `src/models/solver/ImplicitSolver.jl`
- Modify/Delete: fixed-slice compatibility code in `src/models/solver/Conditions.jl`
- Modify/Delete: transitional compatibility branches in `src/models/constraint_solver.jl`
- Modify: `docs/api/models/solver/Overview.md`
- Modify: `docs/api/models/solver/ConstraintModes.md`

- [ ] **Step 1: 写失败测试（确保不存在旧路径依赖）**
  - 新增断言：主流程不通过旧 fixed-dimension helpers；兼容入口调用应触发明确迁移错误或已删除。

- [ ] **Step 2: 运行目标测试确认失败**
  - Run: `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_solver_dimension_agnostic.jl"; include("tests/unit/runtests.jl")'`

- [ ] **Step 3: 删除旧架构代码与 compat 分支**
  - 移除固定切片 helper、过渡 shim、已废弃路径。

- [ ] **Step 4: 全量 smoke + regression + governance**
  - Run:
    - `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
    - `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
    - `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
    - `julia --project=. scripts/dev/check_docs_consistency.jl`
    - `julia --project=. scripts/dev/check_active_docs_governance.jl`

- [ ] **Step 5: 提交（小步）**
  - Message style suggestion: `refactor(models/solver): remove legacy fixed-dimension compat paths`

## Verification Matrix (Final Gate)

- [ ] `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
- [ ] `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
- [ ] `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
- [ ] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [ ] `julia --project=. scripts/dev/check_active_docs_governance.jl`
- [ ] `julia --project=. -e 'include("src/constants/Constants_PNJL.jl"); include("src/models/Models.jl"); println(Models.registered_model_kinds())'`（确认覆盖矩阵与注册模型一致）

## Risk Register

- 兼容风险：外部调用仍假设 `x_state` 长度 5；需保留过渡 accessor 并给出 deprecation 提示。
- 数值风险：动态向量路径可能改变初值组装/尺度，需固定回归点监控。
- 维护风险：迁移中出现双轨（旧切片+新 schema）；每个 chunk 完成后立即删除已替代旧路径，避免长期并存。

## Done Definition

- 所有主流程不再出现固定 `5/3/8` 切片常量。
- `scan/workflow/entrypoints` 在支持模型上均通过 smoke。
- `NJL/NJL2/PNJL/PNJLMagnetic/RPNJL/Rotation/GasLiquid` 全模型覆盖测试通过（至少一条 solver/scan/workflow 主链路）。
- 文档与任务单状态同步，治理脚本通过。
- PR checks 全绿，review 线程闭环。
- 旧架构与兼容层已移除，不保留长期双轨。
